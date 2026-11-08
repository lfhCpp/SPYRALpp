#ifndef RK4_RECO_GETRESPONSE_H
#define RK4_RECO_GETRESPONSE_H

//----------------------------------------------------------------------
// Optional GET pulse templates and saturation extrapolation.
// Giovinazzo et al., NIM A 840 (2016), section 4.2 describes
// h(t) ~ exp(-3t/tau)*(t/tau)^3*sin(t/tau), but explicitly identifies it
// as a rough approximation, insufficient for precise deconvolution.
// This analytic template is a model, not a measured per-channel calibration.
// Keep template fitting optional; validate charge/time bias against calibrated
// pulser or detector waveforms before using it for quantitative reconstruction.
//
// Sections 4.3-4.4 describe deconvolution with a noise-dependent low-pass
// transfer function. A real even transfer function (forward/backward product)
// avoids phase shift; even filter order alone does not guarantee zero phase.
// The published time-resolution results apply to the tested electronics/noise
// settings and are not a universal resolution in a fixed number of time buckets.
//
// Saturation uses a 12-bit ADC ceiling of 4095 by default. Extrapolation from
// an analytic template remains model-dependent. Baseline removal is handled
// separately by WaveProcessor.
//----------------------------------------------------------------------

#include <vector>

namespace rk4 {

struct PulseShapeConfig {
    int n_samples = 512;       // 时间桶数
    double sampling_ns = 320;  // Run 0151 为 320 ns/TB (3.125 MHz)；按 run 覆盖
    double peaking_ns = 720;   // 成形/峰化时间（ns）
    double sat_adc = 4095.0;   // ADC 上限（12-bit）
};

class PulseShape {
public:
    void Configure(const PulseShapeConfig& cfg);
    bool Valid() const { return peak_index_ >= 0 && !h_.empty(); }

    const std::vector<double>& Template() const { return h_; } // 峰值归一化到 1
    int PeakIndex() const { return peak_index_; }
    double TemplateIntegral() const { return integral_; } // Σ h（离散）
    int Size() const { return static_cast<int>(h_.size()); }
    const PulseShapeConfig& Config() const { return cfg_; }

    // 模板在某（可带小数的）时间偏移处的取值：线性插值
    double At(double offset) const;

private:
    PulseShapeConfig cfg_;
    std::vector<double> h_;
    int peak_index_ = -1;
    double integral_ = 0.0;
};

//----------------------------------------------------------------------
// 饱和检测
//----------------------------------------------------------------------
struct SaturationInfo {
    bool saturated = false;
    int first = -1;   // 第一个饱和样本
    int last = -1;    // 最后一个饱和样本
    int count = 0;    // 平顶样本数
};

SaturationInfo DetectSaturation(const std::vector<double>& trace, double sat_adc,
                                int min_samples = 2);

//----------------------------------------------------------------------
// 峰参数
//----------------------------------------------------------------------
struct FittedPeak {
    double time = 0.0;      // 亚采样位置（时间桶）
    double amplitude = 0.0; // 真实峰值幅度（ADC）
    double charge = 0.0;    // A·Σh，正比于电荷
    double chi2 = 0.0;      // 拟合残差平方和（未饱和样本）
    int n_points = 0;       // 参与拟合的样本数
    bool saturation_corrected = false;
};

struct TemplateFitConfig {
    double sat_adc = 4095.0;       // 饱和电平（与 PulseShapeConfig 保持一致）
    int max_peaks = 8;             // 单条迹最多拟合几个峰
    double residual_threshold = 40.0; // 残差峰的幅度阈值（ADC）
    double separation = 3.0;       // 相邻峰最小间距（时间桶）
    bool correct_saturation = true;
    int min_fit_points = 6;        // 少于这么多未饱和样本时改用平顶宽度法
};

// 在 [begin, end) 区间内拟合峰；candidates 是候选位置（通常来自 PeakFinder）
std::vector<FittedPeak> FitPeaksWithTemplate(const std::vector<double>& trace,
                                             const PulseShape& resp,
                                             const SaturationInfo& sat,
                                             const std::vector<double>& candidates,
                                             const TemplateFitConfig& cfg);

// 平顶宽度反解幅度：找 A，使得 A·h 超过 sat_adc 的样本数等于 plateau
double AmplitudeFromPlateau(const PulseShape& resp, double sat_adc, int plateau);

} // namespace rk4

#endif // RK4_RECO_GETRESPONSE_H
