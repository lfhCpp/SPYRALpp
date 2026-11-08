#ifndef RK4_RECO_H5TRACEREADER_H
#define RK4_RECO_H5TRACEREADER_H

//----------------------------------------------------------------------
// HDF5 原始波形读取（GET 数据布局）
//
// 布局（与 Spyral 的 MergerLegacyReader 一致）：
//   /get/evt{N}_data      dataset (512, 10240) int32
//                        第一维 = 时间桶，第二维 = pad 号
//
// 这是真实 AT-TPC 数据的格式。FRIB 的 scaler/IC 数据暂不读取。
//
// 依赖 HDF5 C API（仓库内 third_party/hdf5），未启用 RK4_WITH_HDF5 时不可用。
//----------------------------------------------------------------------

#include "TraceReader.h"

#include <string>
#include <vector>

namespace rk4 {

class H5TraceReader : public TraceReader {
public:
    ~H5TraceReader() override;

    bool Open(const std::string& path) override;
    bool ReadNext(TraceEvent& ev) override;
    int NumEvents() const override { return static_cast<int>(event_ids_.size()); }

    // 可选：只读取某个事件区间
    void SetEventRange(int first, int last) {
        first_event_ = first;
        last_event_ = last;
    }

    // 是否有可用的 HDF5 支持
    static bool Available();

private:
    void* file_ = nullptr;          // hid_t（hid_t 是整数，用 void* 避免头文件泄漏）
    std::vector<int> event_ids_;    // 文件里实际存在的事件号
    size_t cursor_ = 0;
    int first_event_ = -1;
    int last_event_ = -1;
    int trace_length_ = 512;
    int num_pads_ = 10240;
    bool has_frib_events_ = false;
};

} // namespace rk4

#endif // RK4_RECO_H5TRACEREADER_H
