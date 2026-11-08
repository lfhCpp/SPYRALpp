//----------------------------------------------------------------------
// GET 电子学响应 / 饱和修复 / 模板拟合测试（不需要真实数据）
//   root -l -b -q test/testPulseShape.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/PulseShape.h"

#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

// 用已知模板造一条迹：A·h(k - start)，可叠加多个峰、可加噪、可削顶。
// 注意：resp.At() 的偏移是"模板起点"，所以**峰值时刻** = start + resp.PeakIndex()。
static std::vector<double> MakeTrace(const PulseShape& resp, const std::vector<double>& amps,
                                     const std::vector<double>& times, double noise,
                                     double sat_adc, int seed) {
    std::vector<double> trace(512, 0.0);
    for (size_t i = 0; i < amps.size(); ++i) {
        for (int k = 0; k < 512; ++k) {
            trace[static_cast<size_t>(k)] += amps[i] * resp.At(static_cast<double>(k) - times[i]);
        }
    }
    std::mt19937 rng(static_cast<unsigned>(seed));
    std::normal_distribution<double> g(0.0, noise);
    for (double& v : trace) {
        v += g(rng);
        if (sat_adc > 0.0 && v > sat_adc) v = sat_adc; // ADC 削顶
    }
    return trace;
}

void testPulseShape() {
    std::printf("========== GET 响应 / 饱和修复测试 ==========\n");

    PulseShapeConfig rc;
    rc.n_samples = 512;
    rc.sampling_ns = 160.0;
    rc.peaking_ns = 1000.0;
    rc.sat_adc = 4095.0;
    PulseShape resp;
    resp.Configure(rc);

    Check(resp.Valid(), "响应模板有效");
    const std::vector<double>& h = resp.Template();
    Check(!h.empty() && std::fabs(h[static_cast<size_t>(resp.PeakIndex())] - 1.0) < 1e-6,
          "模板峰值归一化到 1");
    Check(resp.TemplateIntegral() > 0.0, "模板积分 > 0");

    TemplateFitConfig tc;
    tc.sat_adc = rc.sat_adc;
    tc.max_peaks = 8;
    tc.residual_threshold = 40.0;

    // ---- 1) 单峰：位置/幅度 ----
    {
        const std::vector<double> amps = {1000.0};
        const std::vector<double> times = {120.0};
        const std::vector<double> trace = MakeTrace(resp, amps, times, 2.0, 0.0, 1);
        const SaturationInfo sat = DetectSaturation(trace, rc.sat_adc, 2);
        Check(!sat.saturated, "未饱和迹线不误报");

        const double true_peak = times[0] + resp.PeakIndex();
        std::vector<double> cand = {true_peak};
        const auto fitted = FitPeaksWithTemplate(trace, resp, sat, cand, tc);
        Check(fitted.size() == 1, "单峰拟合出 1 个峰");
        if (fitted.size() == 1) {
            Check(std::fabs(fitted[0].time - true_peak) <= 1.5, "峰位置误差 <= 1.5 TB");
            Check(std::fabs(fitted[0].amplitude - 1000.0) / 1000.0 < 0.05, "幅度误差 < 5%");
        }
    }

    // ---- 2) 双峰（重叠）：能分开 ----
    {
        const std::vector<double> amps = {800.0, 600.0};
        const std::vector<double> times = {150.0, 162.0};
        const std::vector<double> trace = MakeTrace(resp, amps, times, 2.0, 0.0, 2);
        const SaturationInfo sat = DetectSaturation(trace, rc.sat_adc, 2);
        std::vector<double> cand = {times[0] + resp.PeakIndex(), times[1] + resp.PeakIndex()};
        const auto fitted = FitPeaksWithTemplate(trace, resp, sat, cand, tc);
        Check(fitted.size() == 2, "重叠双峰能拟合成 2 个");
        if (fitted.size() == 2) {
            const double t0 = times[0] + resp.PeakIndex();
            const double t1 = times[1] + resp.PeakIndex();
            Check(std::fabs(fitted[0].time - t0) <= 2.0, "第一峰位置");
            Check(std::fabs(fitted[1].time - t1) <= 2.0, "第二峰位置");
        }
    }

    // ---- 3) 饱和：幅度被修复（明显高于削顶值）----
    {
        const std::vector<double> amps = {9000.0}; // 远超 4095，会被削顶
        const std::vector<double> times = {120.0};
        const std::vector<double> trace = MakeTrace(resp, amps, times, 2.0, rc.sat_adc, 3);
        const SaturationInfo sat = DetectSaturation(trace, rc.sat_adc, 2);
        Check(sat.saturated, "检测到饱和");
        Check(sat.count >= 2, "平顶样本数 >= 2");

        std::vector<double> cand = {times[0] + resp.PeakIndex()};
        const auto fitted = FitPeaksWithTemplate(trace, resp, sat, cand, tc);
        Check(fitted.size() == 1, "饱和迹线也能拟合");
        if (fitted.size() == 1) {
            // 削顶后的原始峰高只有 ~4095；修复后应明显更高
            const double peak_in_trace = 4095.0;
            Check(fitted[0].amplitude > peak_in_trace, "修复后的幅度高于削顶值");
            Check(std::fabs(fitted[0].amplitude - 9000.0) / 9000.0 < 0.15,
                  "修复后的幅度接近真值（15% 以内）");
        }
    }

    // ---- 4) 平顶宽度法单独校验 ----
    {
        const double plateau_level = 4095.0;
        const int plateau = 6;
        const double amp = AmplitudeFromPlateau(resp, plateau_level, plateau);
        int cnt = 0;
        for (int i = 0; i < resp.Size(); ++i) {
            if (amp * resp.Template()[static_cast<size_t>(i)] >= plateau_level) ++cnt;
        }
        Check(std::fabs(cnt - plateau) <= 1, "平顶宽度反解自洽");
        Check(amp > plateau_level, "反解幅度高于平顶电平");
    }

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
