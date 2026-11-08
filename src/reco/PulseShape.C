#include "PulseShape.h"

#include <algorithm>
#include <cmath>

namespace rk4 {

//----------------------------------------------------------------------
// 解析响应（与 attpc_engine 的 get_response 同一条公式）
//   c2 = t / tau，h ∝ exp(-3 c2) c2^3 sin(c2)，负值截 0，峰值归一化
//----------------------------------------------------------------------
void PulseShape::Configure(const PulseShapeConfig& cfg) {
    cfg_ = cfg;
    h_.assign(static_cast<size_t>(std::max(8, cfg.n_samples)), 0.0);
    integral_ = 0.0;
    peak_index_ = -1;

    if (!(cfg.peaking_ns > 0.0) || !(cfg.sampling_ns > 0.0)) return;
    const double tau = cfg.peaking_ns;

    double hmax = 0.0;
    for (int i = 0; i < cfg_.n_samples; ++i) {
        const double t = i * cfg_.sampling_ns;
        const double c2 = t / tau;
        double v = std::exp(-3.0 * c2) * c2 * c2 * c2 * std::sin(c2);
        if (v < 0.0) v = 0.0;
        h_[static_cast<size_t>(i)] = v;
        hmax = std::max(hmax, v);
    }
    if (!(hmax > 0.0)) return;

    for (double& v : h_) v /= hmax; // 峰值归一化到 1
    peak_index_ = static_cast<int>(std::max_element(h_.begin(), h_.end()) - h_.begin());
    for (double v : h_) integral_ += v;
}

double PulseShape::At(double offset) const {
    if (h_.empty()) return 0.0;
    if (offset <= 0.0 || offset >= static_cast<double>(h_.size() - 1)) return 0.0;
    const int i = static_cast<int>(std::floor(offset));
    const double f = offset - i;
    return h_[static_cast<size_t>(i)] * (1.0 - f) + h_[static_cast<size_t>(i + 1)] * f;
}

SaturationInfo DetectSaturation(const std::vector<double>& trace, double sat_adc,
                                int min_samples) {
    SaturationInfo info;
    if (trace.empty() || !(sat_adc > 0.0)) return info;

    int run = 0, best_run = 0, best_start = -1, start = -1;
    for (int i = 0; i < static_cast<int>(trace.size()); ++i) {
        if (trace[static_cast<size_t>(i)] >= sat_adc) {
            if (run == 0) start = i;
            ++run;
            if (run > best_run) {
                best_run = run;
                best_start = start;
            }
        } else {
            run = 0;
        }
    }
    if (best_run >= std::max(1, min_samples)) {
        info.saturated = true;
        info.first = best_start;
        info.last = best_start + best_run - 1;
        info.count = best_run;
    }
    return info;
}

double AmplitudeFromPlateau(const PulseShape& resp, double sat_adc, int plateau) {
    if (!resp.Valid() || plateau <= 0 || !(sat_adc > 0.0)) return sat_adc;
    // h 已归一化：A·h >= sat_adc 的样本数随 A 单调增加 → 二分反解
    double lo = sat_adc;
    double hi = sat_adc * 1.0e3;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        int cnt = 0;
        for (int i = 0; i < resp.Size(); ++i) {
            if (mid * resp.Template()[static_cast<size_t>(i)] >= sat_adc) ++cnt;
        }
        if (cnt >= plateau) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    return hi;
}

namespace {

double SatLevel(const TemplateFitConfig& cfg, const SaturationInfo& sat) {
    return sat.saturated ? cfg.sat_adc : cfg.sat_adc;
}

// 单峰模板拟合：t0 整数网格 + 幅度线性最小二乘解析解。
// 饱和样本（被 ADC 削平）不参与拟合，只提供"平顶宽度"这条额外信息。
FittedPeak FitOnePeak(const std::vector<double>& trace, const PulseShape& resp,
                      const SaturationInfo& sat, double cand,
                      const TemplateFitConfig& cfg) {
    FittedPeak best;
    best.time = cand;

    const int n = static_cast<int>(trace.size());
    const int peak_idx = resp.PeakIndex();
    const int half_win = std::max(8, peak_idx * 3);
    const int lo = std::max(0, static_cast<int>(cand) - half_win);
    const int hi = std::min(n, static_cast<int>(cand) + half_win + 1);

    const int t0_center = static_cast<int>(std::floor(cand)) - peak_idx;
    double best_sse = 1.0e300;

    for (int t0 = t0_center - 3; t0 <= t0_center + 3; ++t0) {
        double num = 0.0, den = 0.0;
        int used = 0;
        for (int k = lo; k < hi; ++k) {
            if (sat.saturated && cfg.correct_saturation && k >= sat.first && k <= sat.last) {
                continue;
            }
            const double hk = resp.At(static_cast<double>(k - t0));
            const double y = trace[static_cast<size_t>(k)];
            num += y * hk;
            den += hk * hk;
            ++used;
        }
        if (used < 3 || !(den > 1.0e-12)) continue;

        const double A = num / den;
        double sse = 0.0;
        for (int k = lo; k < hi; ++k) {
            if (sat.saturated && cfg.correct_saturation && k >= sat.first && k <= sat.last) {
                continue;
            }
            const double model = A * resp.At(static_cast<double>(k - t0));
            const double r = trace[static_cast<size_t>(k)] - model;
            sse += r * r;
        }
        if (sse < best_sse) {
            best_sse = sse;
            best.time = static_cast<double>(t0 + peak_idx); // 峰的绝对位置
            best.amplitude = A;
            best.n_points = used;
        }
    }

    // 未饱和样本太少（平顶很宽）→ 用平顶宽度反解幅度
    if (sat.saturated && cfg.correct_saturation && best.n_points < cfg.min_fit_points) {
        best.amplitude = AmplitudeFromPlateau(resp, cfg.sat_adc, sat.count);
        best.time = 0.5 * (sat.first + sat.last);
        best.saturation_corrected = true;
    }

    best.charge = best.amplitude * resp.TemplateIntegral();
    best.chi2 = best_sse;
    (void)SatLevel(cfg, sat);
    return best;
}

} // namespace

std::vector<FittedPeak> FitPeaksWithTemplate(const std::vector<double>& trace,
                                             const PulseShape& resp,
                                             const SaturationInfo& sat,
                                             const std::vector<double>& candidates,
                                             const TemplateFitConfig& cfg) {
    std::vector<FittedPeak> out;
    if (!resp.Valid() || trace.empty() || candidates.empty()) return out;

    std::vector<double> cand = candidates;
    std::sort(cand.begin(), cand.end());

    // 逐候选拟合（都在原始迹线上做；不在残差上反复加峰，否则噪声会被拟成假峰）
    for (double c : cand) {
        if (static_cast<int>(out.size()) >= cfg.max_peaks) break;
        FittedPeak p = FitOnePeak(trace, resp, sat, c, cfg);
        if (p.amplitude < cfg.residual_threshold) continue;

        // 与已收下的峰太近就合并（保留幅度大的）
        bool merged = false;
        for (auto& e : out) {
            if (std::fabs(e.time - p.time) < cfg.separation) {
                if (p.amplitude > e.amplitude) e = p;
                merged = true;
                break;
            }
        }
        if (!merged) out.push_back(p);
    }

    std::sort(out.begin(), out.end(),
              [](const FittedPeak& a, const FittedPeak& b) { return a.time < b.time; });
    return out;
}

} // namespace rk4
