// 扫沿路径降能强度 pred_strength（用真实 Estimator + catima，前 300 个点云）
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

void pid_strength() {
    std::vector<PointCloud> clouds_all;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds_all)) return;
    std::vector<PointCloud> clouds(clouds_all.begin(),
                                   clouds_all.begin() + std::min<size_t>(300, clouds_all.size()));

    std::printf("%8s %8s %10s %10s %8s\n", "strength", "n", "RMS(dex)", "r(polar)", "±15%");
    for (double s : {0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0}) {
        Clusterer clusterer;
        std::vector<double> dedx, pred, brho, pol;
        for (const auto& cloud : clouds) {
            std::vector<Cluster> clusters;
            clusterer.Run(cloud, clusters);
            for (const auto& c : clusters) {
                if (c.points.size() < 12) continue;
                Estimator est;
                est.SetSmallPadRadius(152.0);
                est.SetBeamRegionRadius(40.0);
                est.SetTarget(1, 1, 2, 6.6168e-05);
                est.SetPredStrength(s);
                EstimateResult r;
                if (!est.Estimate(c, cloud, 2.85, r)) continue;
                if (!(r.brho > 0) || !(r.dEdx > 0) || !(r.dedx_pred > 0)) continue;
                dedx.push_back(r.dEdx);
                pred.push_back(r.dedx_pred);
                brho.push_back(r.brho);
                pol.push_back(r.polar);
            }
        }
        if (dedx.empty()) { std::printf("%8.1f %8d\n", s, 0); continue; }
        std::vector<double> ratio;
        for (size_t i = 0; i < dedx.size(); ++i) ratio.push_back(dedx[i] / pred[i]);
        const double c = 1.0 / std::max(1e-30, Med(ratio));
        std::vector<double> res, p;
        int in15 = 0;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double rr = dedx[i] / (c * pred[i]);
            res.push_back(std::log10(rr));
            p.push_back(pol[i]);
            if (std::fabs(rr - 1.0) < 0.15) ++in15;
        }
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
        std::printf("%8.2f %8zu %10.3f %10.2f %7.0f%%\n", s, res.size(),
                    std::sqrt(s2 / res.size()),
                    (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0,
                    100.0 * in15 / res.size());
    }
}
