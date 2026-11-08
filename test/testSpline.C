//----------------------------------------------------------------------
// 平滑样条单元测试：
//   1) 对带噪声的平滑曲线降噪（残差应下降）
//   2) 不过度失真（平滑后仍贴近真值）
//   root -l -b -q test/testSpline.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/Spline.h"

#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

void testSpline() {
    std::printf("========== 平滑样条测试 ==========\n");

    // 构造：z 从 0 到 200 mm，x = 30*sin(z/50)，加噪声
    const int n = 200;
    std::vector<double> zs, xs;
    zs.reserve(n);
    xs.reserve(n);

    unsigned int seed = 12345;
    auto rand01 = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return static_cast<double>((seed >> 16) & 0x7fff) / 32767.0;
    };

    for (int i = 0; i < n; ++i) {
        const double z = 200.0 * i / (n - 1);
        const double truth = 30.0 * std::sin(z / 50.0);
        const double noise = (rand01() - 0.5) * 6.0; // ±3 mm
        zs.push_back(z);
        xs.push_back(truth + noise);
    }

    // 噪声残差（相对真值）
    double noise_rms = 0.0;
    for (int i = 0; i < n; ++i) {
        const double truth = 30.0 * std::sin(zs[static_cast<size_t>(i)] / 50.0);
        const double d = xs[static_cast<size_t>(i)] - truth;
        noise_rms += d * d;
    }
    noise_rms = std::sqrt(noise_rms / n);

    SmoothingSpline sp;
    Check(sp.Fit(zs, xs, 1.0), "Fit 成功");

    double smooth_rms = 0.0;
    for (int i = 0; i < n; ++i) {
        const double truth = 30.0 * std::sin(zs[static_cast<size_t>(i)] / 50.0);
        const double d = sp.Eval(zs[static_cast<size_t>(i)]) - truth;
        smooth_rms += d * d;
    }
    smooth_rms = std::sqrt(smooth_rms / n);

    std::printf("    噪声残差 RMS = %.2f mm\n", noise_rms);
    std::printf("    平滑后 RMS   = %.2f mm\n", smooth_rms);

    Check(smooth_rms < noise_rms, "平滑后残差应下降");

    // 不应过度失真：平滑结果应仍大致跟随真值（RMS 明显小于噪声）
    Check(smooth_rms < noise_rms * 0.9, "平滑效果显著（降幅 > 10%）");

    // 单点求值不崩
    const double mid = sp.Eval(100.0);
    Check(std::isfinite(mid), "中点求值有限");
    std::printf("    x(100mm) = %.2f\n", mid);

    // 边界求值
    Check(std::isfinite(sp.Eval(zs.front())), "左边界求值有限");
    Check(std::isfinite(sp.Eval(zs.back())), "右边界求值有限");

    // Smooth 批量接口
    std::vector<double> smoothed = sp.Smooth(zs);
    Check(smoothed.size() == zs.size(), "Smooth 返回同样长度");
    double batch_rms = 0.0;
    for (int i = 0; i < n; ++i) {
        const double truth = 30.0 * std::sin(zs[static_cast<size_t>(i)] / 50.0);
        const double d = smoothed[static_cast<size_t>(i)] - truth;
        batch_rms += d * d;
    }
    batch_rms = std::sqrt(batch_rms / n);
    Check(std::fabs(batch_rms - smooth_rms) < 0.5, "Smooth 与 Eval 一致");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
