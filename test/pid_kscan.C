// z 标度（漂移速度）是不是错的？扫 k：L(k)=∫sqrt(k²dz²+dx²+dy²)，看 Bethe 残差怎么变
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

int main_impl() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) return 1;
    Clusterer clusterer;
    std::vector<Cl> cls;

    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& c : clusters) {
            if (c.points.size() < 12) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(40.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
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
            if (cl.z.size() < 8) continue;
            cls.push_back(std::move(cl));
        }
    }
    std::printf("簇 %zu\n\n%6s %8s %10s %10s %8s\n", cls.size(), "k", "N", "RMS(dex)", "corr_polar",
                "±15%");

    auto score = [&](double k) {
        std::vector<double> res, pol;
        int in15 = 0;
        for (const auto& cl : cls) {
            double L = 0.0;
            for (size_t i = 1; i < cl.z.size(); ++i) {
                const double dz = (cl.z[i] - cl.z[i - 1]) * k;
                const double dx = cl.x[i] - cl.x[i - 1];
                const double dy = cl.y[i] - cl.y[i - 1];
                L += std::sqrt(dz * dz + dx * dx + dy * dy);
            }
            double Q = 0.0;
            for (double q : cl.q) Q += q;
            if (!(L > 1.0)) continue;
            const double e = ParticleId::ExpectedDedx(cl.brho, 1, 1);
            if (!(e > 0)) continue;
            const double rr = (Q / L) / e;
            res.push_back(std::log10(rr));
            pol.push_back(cl.polar);
            if (std::fabs(rr - 1.0) < 0.15) ++in15;
        }
        if (res.size() < 10) return;
        std::vector<double> rr = res;
        std::sort(rr.begin(), rr.end());
        const double med = rr[rr.size() / 2];
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        // 残差与极角的相关系数
        double mp = 0, mr = 0;
        for (size_t i = 0; i < pol.size(); ++i) {
            mp += pol[i];
            mr += res[i];
        }
        mp /= pol.size();
        mr /= res.size();
        double sp = 0, sr = 0, spr = 0;
        for (size_t i = 0; i < pol.size(); ++i) {
            const double a = pol[i] - mp, b = res[i] - mr;
            sp += a * a;
            sr += b * b;
            spr += a * b;
        }
        const double corr = (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0;
        std::printf("%6.2f %8zu %10.3f %10.2f %7.0f%%\n", k, res.size(),
                    std::sqrt(s2 / res.size()), corr, 100.0 * in15 / res.size());
    };

    for (double k : {0.4, 0.6, 0.8, 1.0, 1.25, 1.6, 2.0, 2.5}) score(k);
    return 0;
}

void pid_kscan() { main_impl(); }
