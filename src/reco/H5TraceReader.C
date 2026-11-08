#include "H5TraceReader.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>

#ifdef RK4_WITH_HDF5
#include <hdf5.h>
#endif

namespace rk4 {

namespace {

#ifdef RK4_WITH_HDF5
// 从组里按名字找 dataset；找到返回 hid_t，否则返回负数
hid_t FindDataset(hid_t group, const char* name) {
    if (H5Lexists(group, name, H5P_DEFAULT) <= 0) return -1;
    return H5Dopen2(group, name, H5P_DEFAULT);
}
#endif

} // namespace

bool H5TraceReader::Available() {
#ifdef RK4_WITH_HDF5
    return true;
#else
    return false;
#endif
}

H5TraceReader::~H5TraceReader() {
#ifdef RK4_WITH_HDF5
    if (file_ != nullptr) {
        H5Fclose(static_cast<hid_t>(reinterpret_cast<intptr_t>(file_)));
    }
#endif
}

bool H5TraceReader::Open(const std::string& path) {
#ifdef RK4_WITH_HDF5
    event_ids_.clear();
    cursor_ = 0;
    if (file_ != nullptr) {
        H5Fclose(static_cast<hid_t>(reinterpret_cast<intptr_t>(file_)));
        file_ = nullptr;
    }

    const hid_t fid = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0) {
        std::fprintf(stderr, "[H5TraceReader] 无法打开 %s\n", path.c_str());
        return false;
    }
    file_ = reinterpret_cast<void*>(static_cast<intptr_t>(fid));
    has_frib_events_ = H5Lexists(fid, "/frib", H5P_DEFAULT) > 0 &&
                       H5Lexists(fid, "/frib/evt", H5P_DEFAULT) > 0;

    // /get 组
    if (H5Lexists(fid, "/get", H5P_DEFAULT) <= 0) {
        std::fprintf(stderr, "[H5TraceReader] %s 里没有 /get 组\n", path.c_str());
        H5Fclose(fid);
        file_ = nullptr;
        return false;
    }
    const hid_t get_group = H5Gopen2(fid, "/get", H5P_DEFAULT);
    if (get_group < 0) {
        H5Fclose(fid);
        file_ = nullptr;
        return false;
    }

    // 扫描组内所有 "evt{N}_data"（H5Literate）
    struct Ctx {
        std::vector<int>* ids;
    } ctx{&event_ids_};

    auto callback = [](hid_t /*group*/, const char* name, const H5L_info2_t* /*info*/,
                       void* op_data) -> herr_t {
        Ctx* c = static_cast<Ctx*>(op_data);
        int value = -1;
        int end = 0;
        if (std::sscanf(name, "evt%d_data%n", &value, &end) == 1 &&
            end > 0 && name[end] == '\0' && value >= 0) {
            c->ids->push_back(value);
        }
        return 0;
    };

    H5Literate(get_group, H5_INDEX_NAME, H5_ITER_INC, nullptr, callback, &ctx);
    H5Gclose(get_group);

    if (event_ids_.empty()) {
        std::fprintf(stderr, "[H5TraceReader] %s 的 /get 里没有 evt*_data\n", path.c_str());
        H5Fclose(fid);
        file_ = nullptr;
        return false;
    }

    std::sort(event_ids_.begin(), event_ids_.end());
    event_ids_.erase(std::unique(event_ids_.begin(), event_ids_.end()), event_ids_.end());

    // 事件区间过滤
    if (first_event_ >= 0 || last_event_ >= 0) {
        const int lo = first_event_ >= 0 ? first_event_ : event_ids_.front();
        const int hi = last_event_ >= 0 ? last_event_ : event_ids_.back();
        std::vector<int> filtered;
        for (int v : event_ids_) {
            if (v >= lo && v <= hi) filtered.push_back(v);
        }
        event_ids_.swap(filtered);
    }

    return !event_ids_.empty();
#else
    std::fprintf(stderr,
                 "[H5TraceReader] 未启用 HDF5 支持（编译时 RK4_WITH_HDF5=OFF），"
                 "无法读取 %s\n",
                 path.c_str());
    return false;
#endif
}

bool H5TraceReader::ReadNext(TraceEvent& ev) {
#ifdef RK4_WITH_HDF5
    if (file_ == nullptr || cursor_ >= event_ids_.size()) return false;

    const hid_t fid = static_cast<hid_t>(reinterpret_cast<intptr_t>(file_));
    const int event_id = event_ids_[cursor_++];
    const auto fail = [&]() -> bool {
        throw std::runtime_error("Invalid HDF5 GET event " + std::to_string(event_id));
    };
    ev = TraceEvent();
    ev.number = event_id;
    if (has_frib_events_) {
        const std::string prefix = "/frib/evt/evt" + std::to_string(event_id);
        if (H5Lexists(fid, (prefix + "_1903").c_str(), H5P_DEFAULT) > 0 &&
            H5Lexists(fid, (prefix + "_977").c_str(), H5P_DEFAULT) > 0) {
            const hid_t ds = H5Dopen2(fid, (prefix + "_977").c_str(), H5P_DEFAULT);
            if (ds < 0) return fail();
            const hid_t sp = H5Dget_space(ds);
            const auto size = H5Sget_simple_extent_npoints(sp);
            std::vector<int32_t> bits(static_cast<size_t>(std::max<hssize_t>(size, 0)));
            const bool ok = size > 0 && H5Dread(ds, H5T_NATIVE_INT32, H5S_ALL,
                                H5S_ALL, H5P_DEFAULT, bits.data()) >= 0;
            H5Sclose(sp); H5Dclose(ds);
            if (!ok || bits[0] < 0 || bits[0] > 3) return fail();
            ev.ic_downscale = bits[0] == 2;
            const hid_t signal=H5Dopen2(fid,(prefix+"_1903").c_str(),H5P_DEFAULT);
            if(signal<0)return fail();
            const hid_t signal_space=H5Dget_space(signal);
            hsize_t shape[2]={0,0};
            const bool rank_ok=H5Sget_simple_extent_ndims(signal_space)==2;
            if(rank_ok)H5Sget_simple_extent_dims(signal_space,shape,nullptr);
            bool trace_ok=rank_ok&&shape[0]==2048&&shape[1]>0;
            if(trace_ok){
                const hsize_t start[2]={0,0},count[2]={2048,1},length=2048;
                trace_ok=H5Sselect_hyperslab(signal_space,H5S_SELECT_SET,start,nullptr,count,nullptr)>=0;
                const hid_t memory=H5Screate_simple(1,&length,nullptr);
                ev.ic_trace.resize(2048);
                trace_ok=trace_ok&&H5Dread(signal,H5T_NATIVE_FLOAT,memory,signal_space,H5P_DEFAULT,ev.ic_trace.data())>=0;
                H5Sclose(memory);
            }
            H5Sclose(signal_space);H5Dclose(signal);
            if(!trace_ok)return fail();
        }
    }

    char name[64];
    std::snprintf(name, sizeof(name), "/get/evt%d_data", event_id);

    const hid_t dset = H5Dopen2(fid, name, H5P_DEFAULT);
    if (dset < 0) return fail();

    const hid_t space = H5Dget_space(dset);
    hsize_t dims[2] = {0, 0};
    const int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims != 2) {
        H5Sclose(space);
        H5Dclose(dset);
        return fail();
    }
    H5Sget_simple_extent_dims(space, dims, nullptr);

    const hsize_t n_rows = dims[0];
    const hsize_t n_cols = dims[1];
    if (n_rows == 0 || n_cols == 0) {
        H5Sclose(space);
        H5Dclose(dset);
        return true; // An empty event does not imply EOF.
    }

    std::vector<int32_t> buffer(static_cast<size_t>(n_rows * n_cols));
    const herr_t rc = H5Dread(dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              buffer.data());
    H5Sclose(space);
    H5Dclose(dset);
    if (rc < 0) return fail();

    ev.number = event_id;
    ev.traces.clear();

    //------------------------------------------------------------------
    // 三种布局，按维度自动识别：
    //  1) 表格式（libattpc_merger 1.x，Spyral 的 MERGER_1_0）：
    //       行 = 一条迹，列 = [cobo, asad, aget, channel, pad_id, 512 个采样]
    //       维度形如 (386, 517)，也就是 cols = 5 + 512
    //  2) (n_tb, n_pad)  ：512 x 10240，第一维时间桶
    //  3) (n_pad, n_tb)  ：10240 x 512，第一维 pad
    //------------------------------------------------------------------
    const bool table_layout = n_cols == 517;
    const bool pad_major = !table_layout && n_cols == 512 && n_rows != 512;
    if (!table_layout && !pad_major && n_rows != 512) return fail();

    if (table_layout) {
        const hsize_t id_cols = 5;                  // 5 个硬件/几何 id 列
        const hsize_t n_samples = n_cols - id_cols; // 512
        ev.traces.reserve(static_cast<size_t>(n_rows));
        for (hsize_t r = 0; r < n_rows; ++r) {
            const int32_t* row = &buffer[static_cast<size_t>(r * n_cols)];
            const int pad = static_cast<int>(row[4]); // 第 5 列就是 pad_id
            if (pad == -1) continue;

            bool any = false;
            for (hsize_t i = 0; i < n_samples; ++i) {
                if (row[id_cols + i] != 0) {
                    any = true;
                    break;
                }
            }
            if (!any) continue;

            Trace t;
            t.pad = pad;
            t.samples.resize(static_cast<size_t>(n_samples));
            for (hsize_t i = 0; i < n_samples; ++i) {
                t.samples[static_cast<size_t>(i)] =
                    static_cast<float>(row[id_cols + i]);
            }
            ev.traces.push_back(std::move(t));
        }
        return true;
    }

    const hsize_t n_tb = pad_major ? n_cols : n_rows;
    const hsize_t n_pad = pad_major ? n_rows : n_cols;
    ev.traces.reserve(static_cast<size_t>(n_pad));

    // 去掉全零的 pad（真实数据里大量 pad 没有信号）
    for (hsize_t pad = 0; pad < n_pad; ++pad) {
        auto at = [&](hsize_t tb) -> int32_t {
            return pad_major ? buffer[static_cast<size_t>(pad * n_tb + tb)]
                             : buffer[static_cast<size_t>(tb * n_pad + pad)];
        };
        bool any = false;
        for (hsize_t tb = 0; tb < n_tb; ++tb) {
            if (at(tb) != 0) {
                any = true;
                break;
            }
        }
        if (!any) continue;

        Trace t;
        t.pad = static_cast<int>(pad);
        t.samples.resize(static_cast<size_t>(n_tb));
        for (hsize_t tb = 0; tb < n_tb; ++tb) {
            t.samples[static_cast<size_t>(tb)] = static_cast<float>(at(tb));
        }
        ev.traces.push_back(std::move(t));
    }

    return true;
#else
    (void)ev;
    return false;
#endif
}

} // namespace rk4
