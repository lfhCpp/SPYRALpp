//----------------------------------------------------------------------
// HDF5 原始波形读取测试
//   1) 用 HDF5 C API 造一个 /get/evt{N}_data (512, 10240) 文件
//   2) 造一个 libattpc_merger 1.x 的"表格式"文件：
//      行 = 一条迹，列 = [cobo, asad, aget, channel, pad_id, 512 个采样]
//      （真实数据 run_0151.h5 就是这种布局，517 = 5 + 512）
//   3) 用 H5TraceReader 读回两种布局，验证事件号、pad、采样值
//   root -l -b -q test/testH5Trace.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__ADD_INCLUDE_PATH(extern/hdf5/include)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/H5TraceReader.h"

#include "TSystem.h"

#include <hdf5.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

// 造一个只含指定 pad 信号的波形文件
static bool MakeTraceFile(const std::string& path) {
    const int n_tb = 512;
    const int n_pad = 10240;
    const int events[3] = {7, 8, 9};
    // 每个事件给出 (pad, 常量值)，其余 pad 全 0
    const int signal_pads[3] = {100, 200, 300};

    hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (f < 0) return false;

    hid_t g = H5Gcreate2(f, "/get", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) { H5Fclose(f); return false; }

    std::vector<int32_t> data(static_cast<size_t>(n_tb * n_pad));

    for (int e = 0; e < 3; ++e) {
        std::fill(data.begin(), data.end(), 0);

        // 让指定 pad 有一个方波：前 100 个 tb 为 500+e，其余 0
        const int pad = signal_pads[e];
        for (int tb = 0; tb < 100; ++tb) {
            data[static_cast<size_t>(tb * n_pad + pad)] = 500 + e;
        }

        hsize_t dims[2] = {static_cast<hsize_t>(n_tb), static_cast<hsize_t>(n_pad)};
        hid_t space = H5Screate_simple(2, dims, nullptr);
        char name[64];
        std::snprintf(name, sizeof(name), "evt%d_data", events[e]);
        hid_t d = H5Dcreate2(g, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT);
        if (d < 0) {
            H5Sclose(space);
            H5Gclose(g);
            H5Fclose(f);
            return false;
        }
        H5Dwrite(d, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(d);
        H5Sclose(space);
    }

    H5Gclose(g);
    H5Fclose(f);
    return true;
}

// 造一个 merger 1.x 表格式文件：/get/evt10000_data 形状 (2, 517) int16
static bool MakeTableFile(const std::string& path) {
    const int n_samples = 512;
    const int n_cols = 5 + n_samples; // 4 个硬件 id + pad + 512 采样
    const int n_rows = 2;
    const int pads[2] = {7331, 7092};

    hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (f < 0) return false;
    hid_t g = H5Gcreate2(f, "/get", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) { H5Fclose(f); return false; }

    std::vector<int16_t> data(static_cast<size_t>(n_rows * n_cols), 0);
    for (int r = 0; r < n_rows; ++r) {
        int16_t* row = &data[static_cast<size_t>(r * n_cols)];
        row[0] = 1;           // cobo
        row[1] = 2;           // asad
        row[2] = 3;           // aget
        row[3] = 4;           // channel
        row[4] = pads[r];     // pad_id
        for (int i = 0; i < n_samples; ++i) {
            row[5 + i] = static_cast<int16_t>((r + 1) * 10 + (i % 7));
        }
    }

    hsize_t dims[2] = {static_cast<hsize_t>(n_rows), static_cast<hsize_t>(n_cols)};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    hid_t d = H5Dcreate2(g, "evt10000_data", H5T_STD_I16LE, space, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    if (d < 0) {
        H5Sclose(space);
        H5Gclose(g);
        H5Fclose(f);
        return false;
    }
    H5Dwrite(d, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(d);
    H5Sclose(space);
    H5Gclose(g);
    H5Fclose(f);
    return true;
}

void testH5Trace() {
    std::printf("========== HDF5 波形读取测试 ==========\n");

    if (!H5TraceReader::Available()) {
        std::printf("========== SKIP（未启用 HDF5） ==========\n");
        gSystem->Exit(0);
    }

    const std::string path = "/tmp/rk4_h5trace_test.h5";
    gSystem->Unlink(path.c_str());
    Check(MakeTraceFile(path), "造测试文件");

    H5TraceReader reader;
    Check(reader.Open(path), "Open");
    Check(reader.NumEvents() == 3, "读到 3 个事件");

    const int expect_events[3] = {7, 8, 9};
    const int expect_pads[3] = {100, 200, 300};

    for (int i = 0; i < 3; ++i) {
        TraceEvent ev;
        Check(reader.ReadNext(ev), "ReadNext");
        Check(ev.number == expect_events[i], "事件号正确");
        Check(ev.traces.size() == 1, "只保留有信号的 pad");
        if (ev.traces.size() == 1) {
            Check(ev.traces[0].pad == expect_pads[i], "pad 号正确");
            Check(ev.traces[0].samples.size() == 512, "512 个采样");
            Check(std::fabs(ev.traces[0].samples[0] - (500 + i)) < 1e-6,
                  "采样值正确");
            Check(std::fabs(ev.traces[0].samples[200]) < 1e-6, "尾部为零");
        }
    }

    TraceEvent dummy;
    Check(!reader.ReadNext(dummy), "读完返回 false");

    gSystem->Unlink(path.c_str());

    // ---- merger 1.x 表格式（真实数据 run_0151.h5 的布局）----
    const std::string tpath = "/tmp/rk4_h5table_test.h5";
    gSystem->Unlink(tpath.c_str());
    Check(MakeTableFile(tpath), "造表格式测试文件");

    H5TraceReader treader;
    Check(treader.Open(tpath), "表格式 Open");
    Check(treader.NumEvents() == 1, "表格式 1 个事件");
    TraceEvent tev;
    Check(treader.ReadNext(tev), "表格式 ReadNext");
    Check(tev.number == 10000, "表格式事件号取自数据集名");
    Check(tev.traces.size() == 2, "表格式保留两行");
    if (tev.traces.size() == 2) {
        Check(tev.traces[0].pad == 7331, "第 5 列当作 pad_id（第一行）");
        Check(tev.traces[1].pad == 7092, "第 5 列当作 pad_id（第二行）");
        Check(tev.traces[0].samples.size() == 512, "表格式 512 个采样");
        Check(std::fabs(tev.traces[0].samples[0] - 10.0) < 1e-6, "采样从第 6 列开始");
        Check(std::fabs(tev.traces[1].samples[1] - 21.0) < 1e-6, "第二行采样正确");
    }
    gSystem->Unlink(tpath.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
