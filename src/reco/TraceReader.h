#ifndef RK4_RECO_TRACEREADER_H
#define RK4_RECO_TRACEREADER_H

//----------------------------------------------------------------------
// 原始波形读取
//
// 目前支持示例文本格式（data/event_waveform.txt）：
//   ========== event 0 ==========
//   1903_channel0_waveform: ...
//   pad_waveforms:
//   pads 0: <5 个元数据> <512 个采样>
//
// 后续 HDF5（GET/FRIB）只需再实现一个 TraceReader 子类。
//----------------------------------------------------------------------

#include <string>
#include <vector>
#include "../core/EventMetadata.h"

namespace rk4 {

struct Trace {
    int pad = -1;
    std::vector<float> samples;
};

struct TraceEvent {
    int number = -1;
    std::vector<Trace> traces;
    bool ic_downscale = false; // GET is still processed to preserve random draws.
    std::vector<float> ic_trace;
    EventMetadata metadata;
};

class TraceReader {
public:
    virtual ~TraceReader() = default;
    virtual bool Open(const std::string& path) = 0;
    virtual bool ReadNext(TraceEvent& ev) = 0;
    virtual int NumEvents() const = 0;
};

class TextTraceReader : public TraceReader {
public:
    bool Open(const std::string& path) override;
    bool ReadNext(TraceEvent& ev) override;
    int NumEvents() const override { return static_cast<int>(event_starts_.size()); }

private:
    std::vector<std::string> lines_;
    std::vector<size_t> event_starts_;
    size_t cursor_ = 0;
};

} // namespace rk4

#endif // RK4_RECO_TRACEREADER_H
