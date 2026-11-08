// 对齐 Spyral 的 smoothing_factor=100 与 min_total_trajectory_points=30，看 PID 宽度和 polar 分布
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
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(f * v.size()))];
}

void pid_align_A() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/abA_cloud.root", clouds)) return;

    for (double lambda : {1.0}) {
        Clusterer clusterer;
        ClusterParams cp;
        cp.eps_mm = 20.0;
        cp.min_points = 3;
        cp.min_cluster_size = 50;
        clusterer.SetParams(cp);

        std::vector<double> dedx, brho, pol;
        int ntot = 0;
        for (const auto& cloud : clouds) {
            std::vector<Cluster> clusters;
            clusterer.Run(cloud, clusters);
            for (const auto& c : clusters) {
                if (c.points.size() < 30) continue; // min_total_trajectory_points
                ++ntot;
                Estimator est;
                est.SetSmallPadRadius(152.0);
                est.SetBeamRegionRadius(25.0);
                // 单一平滑路径：不启用 P-spline 预平滑（实测更好）
                (void)lambda;
                EstimateResult r;
                if (!est.Estimate(c, cloud, 2.85, r)) continue;
                if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
                dedx.push_back(r.dEdx);
                brho.push_back(r.brho);
                pol.push_back(r.polar * 180.0 / 3.14159265358979323846);
            }
        }
        if (dedx.empty()) { std::printf("lambda=%.0f 无\n", lambda); continue; }

        std::vector<double> ratio;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
            if (e > 0) ratio.push_back(dedx[i] / e);
        }
        const double c = 1.0 / std::max(1e-30, Q(ratio, 0.5));
        std::vector<double> res, p;
        int in15 = 0, in25 = 0;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
            if (!(e > 0)) continue;
            const double rr = dedx[i] / (c * e);
            res.push_back(std::log10(rr));
            p.push_back(pol[i]);
            if (std::fabs(rr - 1) < 0.15) ++in15;
            if (std::fabs(rr - 1) < 0.25) ++in25;
        }
        const double med = Q(res, 0.5);
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        double mp = 0, mr = 0;
        for (size_t i = 0; i < res.size(); ++i) { mp += p[i]; mr += res[i]; }
        mp /= res.size(); mr /= res.size();
        double sp = 0, sr = 0, spr = 0;
        for (size_t i = 0; i < res.size(); ++i) {
            const double x = p[i] - mp, y = res[i] - mr;
            sp += x * x; sr += y * y; spr += x * y;
        }
        std::printf("lambda=%-5.0f 通过 %zu/%d (%2.0f%%)  RMS=%.3f dex  r(pol)=%+.2f  "
                    "±15%%=%.0f%%  ±25%%=%.0f%%  polar[16/50/84]=%.0f/%.0f/%.0f deg\n",
                    lambda, res.size(), ntot, 100.0 * res.size() / std::max(1, ntot),
                    std::sqrt(s2 / res.size()),
                    (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0,
                    100.0 * in15 / res.size(), 100.0 * in25 / res.size(), Q(p, 0.16), Q(p, 0.5),
                    Q(p, 0.84));
    }
}
