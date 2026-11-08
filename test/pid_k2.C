// z 标度 k 扫描（新参数：npts>=30、beam 25、spy_cloud）
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

struct Cl {
    std::vector<double> z, x, y, q;
    double brho = 0, polar = 0;
};

static double Q(std::vector<double> v, double f) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(f * v.size()))];
}

void pid_k2() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/spy_cloud.root", clouds)) { std::printf("无点云\n"); return; }
    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20.0;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);

    std::vector<Cl> cls;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& c : clusters) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0)) continue;
            Cl cl;
            cl.brho = r.brho;
            cl.polar = r.polar;
            for (const auto& p : c.points) {
                if (std::hypot(p.x, p.y) > 152.0) continue;
                cl.z.push_back(p.z);
                cl.x.push_back(p.x);
                cl.y.push_back(p.y);
                cl.q.push_back(p.integral);
            }
            if (cl.z.size() < 20) continue;
            cls.push_back(std::move(cl));
        }
    }
    std::printf("簇 %zu\n%6s %8s %10s %10s %8s\n", cls.size(), "k", "n", "RMS(dex)", "r(polar)",
                "±15%");
    for (double k : {0.3, 0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0, 3.0}) {
        std::vector<double> ratio;
        for (const auto& cl : cls) {
            double L = 0.0;
            for (size_t i = 1; i < cl.z.size(); ++i) {
                const double dz = (cl.z[i] - cl.z[i - 1]) * k;
                const double dx = cl.x[i] - cl.x[i - 1];
                const double dy = cl.y[i] - cl.y[i - 1];
                L += std::sqrt(dz * dz + dx * dx + dy * dy);
            }
            double qs = 0.0;
            for (double v : cl.q) qs += v;
            const double e = ParticleId::ExpectedDedx(cl.brho, 1, 1);
            if (L > 1.0 && e > 0) ratio.push_back((qs / L) / e);
        }
        if (ratio.size() < 10) { std::printf("%6.2f 无\n", k); continue; }
        const double c = 1.0 / Q(ratio, 0.5);
        std::vector<double> res, pol;
        int in15 = 0;
        for (size_t i = 0; i < cls.size(); ++i) {
            // 重新算一遍（顺序与 ratio 一致）
        }
        // 用带极角的方式重算
        for (const auto& cl : cls) {
            double L = 0.0;
            for (size_t i = 1; i < cl.z.size(); ++i) {
                const double dz = (cl.z[i] - cl.z[i - 1]) * k;
                const double dx = cl.x[i] - cl.x[i - 1];
                const double dy = cl.y[i] - cl.y[i - 1];
                L += std::sqrt(dz * dz + dx * dx + dy * dy);
            }
            double qs = 0.0;
            for (double v : cl.q) qs += v;
            const double e = ParticleId::ExpectedDedx(cl.brho, 1, 1);
            if (!(L > 1.0) || !(e > 0)) continue;
            const double rr = (qs / L) / (c * e);
            res.push_back(std::log10(rr));
            pol.push_back(cl.polar * 180.0 / 3.14159265358979323846);
            if (std::fabs(rr - 1) < 0.15) ++in15;
        }
        const double med = Q(res, 0.5);
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        double mp = 0, mr = 0;
        for (size_t i = 0; i < res.size(); ++i) { mp += pol[i]; mr += res[i]; }
        mp /= res.size(); mr /= res.size();
        double sp = 0, sr = 0, spr = 0;
        for (size_t i = 0; i < res.size(); ++i) {
            const double a = pol[i] - mp, b = res[i] - mr;
            sp += a * a; sr += b * b; spr += a * b;
        }
        std::printf("%6.2f %8zu %10.3f %10.2f %7.0f%%\n", k, res.size(),
                    std::sqrt(s2 / res.size()),
                    (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0,
                    100.0 * in15 / res.size());
    }
}
