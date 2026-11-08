// 检查 Bρ 分布与 RMS（目标：Bρ 中位 ≈ 0.25，RMS ≤ 0.265）
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

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
static double Q(std::vector<double> v, double f) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}

void brho_check() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) return;
    Clusterer cl;
    std::vector<double> brho, dedx, res, polar;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        for (const auto& c : cs) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            brho.push_back(r.brho);
            dedx.push_back(r.dEdx);
            polar.push_back(r.polar * 180.0 / M_PI);
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (e > 0) res.push_back(std::log10(r.dEdx / e));
        }
    }
    if (brho.empty()) { std::printf("无\n"); return; }
    const double med = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - med) * (v - med);
    std::printf("N=%zu  Bρ中位=%.3f  p16/p84=%.3f/%.3f  dEdx中位=%.0f  polar中位=%.0f°  RMS=%.3f\n",
                brho.size(), Med(brho), Q(brho, 0.16), Q(brho, 0.84), Med(dedx), Med(polar),
                std::sqrt(s2 / res.size()));
}
