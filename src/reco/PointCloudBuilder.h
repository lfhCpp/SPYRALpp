#ifndef RK4_RECO_POINTCLOUDBUILDER_H
#define RK4_RECO_POINTCLOUDBUILDER_H

//----------------------------------------------------------------------
// 波形 -> 点云
//
// 每条波形先做基线扣除（FFT sinc 滤波），再寻峰；
// 峰值换成 (pad, tb) 后查 pad 平面得到 x/y，并用时间参考换算 z。
//----------------------------------------------------------------------

#include "../sim/DetectorConfig.h"
#include "../sim/PadPlane.h"
#include "../sim/Simulator.h"
#include "PulseShape.h"
#include "PeakFinder.h"
#include "TraceReader.h"
#include "WaveProcessor.h"
#include "NumpyPcg64.h"
#include "FribConfig.h"

#include <memory>
#include <random>
#include <vector>

namespace rk4 {

// 一条迹线上的一个峰（写中间过程文件用）
struct TracePeakInfo {
    size_t dither_index = 0;
    int pad = -1;
    int index = 0;
    int npeaks = 0;
    double time = 0.0;
    double amplitude = 0.0;
    double charge = 0.0;
    double chi2 = 0.0;
    int npoints = 0;
    int saturated = 0;
};

class PointCloudBuilder {
public:
    ~PointCloudBuilder();

    bool Init(const DetectorConfig& cfg, const PadPlane* pads,
              const PeakFinderParameters& peak_params,
              double baseline_window_scale = 20.0);

    bool Build(const TraceEvent& ev, PointCloud& cloud);
    void ConfigureFrib(const FribConfig& config) { frib_config_=config;frib_wave_.reset(); }
    void SetDeferredDither(bool value) { defer_dither_ = value; }
    void SetDitherSeed(uint64_t seed) { rng_.Seed(seed); }
    static void FinalizeDither(PointCloud& cloud, NumpyPcg64& rng, const DetectorConfig& cfg,
                               std::vector<TracePeakInfo>* peaks = nullptr);

    // 可选：用 GET 电子学响应做模板拟合（含饱和修复）
    void ConfigureTemplateFit(bool enable, const PulseShapeConfig& rc,
                              const TemplateFitConfig& tc);

    int SaturatedTraces() const { return saturated_traces_; }

    // 可选：把每条迹线的峰记到这个 sink 里（nullptr = 不记录）
    void SetTracePeakSink(std::vector<TracePeakInfo>* sink) { peak_sink_ = sink; }
    int AcceptedTraces() const { return accepted_traces_; }
    int TotalPeaks() const { return total_peaks_; }

private:
    DetectorConfig cfg_;
    const PadPlane* pads_ = nullptr;
    PeakFinderParameters peak_params_;
    std::unique_ptr<WaveProcessor> wave_;
    std::unique_ptr<WaveProcessor> frib_wave_;
    FribConfig frib_config_;
    NumpyPcg64 rng_;
    bool defer_dither_ = false;
    std::vector<double> input_;
    std::vector<double> filtered_;
    bool template_fit_ = false;
    PulseShape response_;
    PulseShapeConfig response_cfg_;
    TemplateFitConfig template_cfg_;
    int saturated_traces_ = 0;
    std::vector<TracePeakInfo>* peak_sink_ = nullptr;
    std::vector<double> fit_candidates_;
    std::vector<double> fit_samples_;
    int accepted_traces_ = 0;
    int total_peaks_ = 0;
    int trace_length_ = 0;
};

} // namespace rk4

#endif // RK4_RECO_POINTCLOUDBUILDER_H
