// 用 Estimator 内部的沿路径 catima 预测（dedx_pred）当参考 vs 顶点 Bethe
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

void pid_pred() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) return;
    Clusterer clusterer;
    std::vector<double> dedx, pred, brho, pol;
    int ntot = 0;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& c : clusters) {
            if (c.points.size() < 12) continue;
            ++ntot;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(40.0);
            est.SetTarget(1, 1, 2, 6.6168e-05);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0) || !(r.dedx_pred > 0)) continue;
            dedx.push_back(r.dEdx);
            pred.push_back(r.dedx_pred);
            brho.push_back(r.brho);
            pol.push_back(r.polar);
        }
    }
    std::printf("簇 %d -> 通过 %zu (%.0f%%)\n", ntot, dedx.size(),
                100.0 * dedx.size() / std::max(1, ntot));

    for (int ip = 0; ip < 2; ++ip) {
        const bool use_pred = (ip == 1);
        std::vector<double> ratio, res, p;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double ref =
                use_pred ? pred[i] : ParticleId::ExpectedDedx(brho[i], 1, 1) * 1e3;
            if (ref > 0) ratio.push_back(dedx[i] / ref);
        }
        const double c = 1.0 / std::max(1e-30, Med(ratio));
        int in15 = 0, in25 = 0;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double ref =
                use_pred ? pred[i] : ParticleId::ExpectedDedx(brho[i], 1, 1) * 1e3;
            if (!(ref > 0)) continue;
            const double rr = dedx[i] / (c * ref);
            res.push_back(std::log10(rr));
            p.push_back(pol[i]);
            if (std::fabs(rr - 1.0) < 0.15) ++in15;
            if (std::fabs(rr - 1.0) < 0.25) ++in25;
        }
        if (res.empty()) continue;
        const double med = Med(res);
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        double mp = 0, mr = 0;
        for (size_t i = 0; i < res.size(); ++i) { mp += p[i]; mr += res[i]; }
        mp /= res.size(); mr /= res.size();
        double sp = 0, sr = 0, spr = 0;
        for (size_t i = 0; i < res.size(); ++i) {
            const double a = p[i] - mp, b = res[i] - mr;
            sp += a * a; sr += b * b; spr += a * b;
        }
        std::printf("%-24s n=%zu RMS=%.3f dex  r(pol)=%+.2f  ±15%%=%.0f%%  ±25%%=%.0f%%\n",
                    use_pred ? "参考=路径catima平均" : "参考=Bethe(顶点)", res.size(),
                    std::sqrt(s2 / res.size()),
                    (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0,
                    100.0 * in15 / res.size(), 100.0 * in25 / res.size());
    }
}
