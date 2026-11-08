// 残差分布形状：为什么 ±15% 是 0
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

static double Q(std::vector<double> v, double f) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(f * v.size()))];
}

void pid_hist() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/spy_cloud.root", clouds)) return;
    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);

    std::vector<double> res, ratio;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        clusterer.Run(cloud, cs);
        for (const auto& c : cs) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (!(e > 0)) continue;
            ratio.push_back(r.dEdx / e);
        }
    }
    if (ratio.empty()) { std::printf("无\n"); return; }
    std::vector<double> r2 = ratio;
    std::sort(r2.begin(), r2.end());
    const double med = r2[r2.size() / 2];
    for (double v : ratio) res.push_back(std::log10(v / med));

    std::printf("簇 %zu，刻度(中位比值)=%.3g\n", ratio.size(), med);
    std::printf("log10 残差分位 5/16/50/84/95 = %.2f / %.2f / %.2f / %.2f / %.2f\n",
                Q(res, 0.05), Q(res, 0.16), Q(res, 0.50), Q(res, 0.84), Q(res, 0.95));
    std::printf("线性比值分位 16/50/84 = %.2f / %.2f / %.2f\n", Q(ratio, 0.16) / med,
                Q(ratio, 0.5) / med, Q(ratio, 0.84) / med);

    // 粗直方图
    const int nb = 20;
    const double lo = -1.0, hi = 1.0;
    int h[nb] = {0};
    for (double v : res) {
        int b = (int)((v - lo) / (hi - lo) * nb);
        if (b >= 0 && b < nb) ++h[b];
    }
    std::printf("\nlog10 残差直方图（-1.0 .. 1.0）\n");
    for (int i = 0; i < nb; ++i) {
        std::printf("%+.2f %5d %s\n", lo + (i + 0.5) * (hi - lo) / nb, h[i],
                    std::string(std::min(60, h[i] * 80 / std::max(1, *std::max_element(h, h + nb))), '#').c_str());
    }
    int in15 = 0, in25 = 0, in50 = 0;
    for (double v : res) {
        const double rr = std::pow(10.0, v);
        if (std::fabs(rr - 1) < 0.15) ++in15;
        if (std::fabs(rr - 1) < 0.25) ++in25;
        if (std::fabs(rr - 1) < 0.50) ++in50;
    }
    std::printf("\n±15%%: %.0f%%   ±25%%: %.0f%%   ±50%%: %.0f%%\n", 100.0 * in15 / res.size(),
                100.0 * in25 / res.size(), 100.0 * in50 / res.size());
}
