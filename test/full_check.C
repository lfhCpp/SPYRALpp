// 全量验证：默认(template) vs P3(局部极大 prom10 thr20)，用当前默认聚类参数
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

static void One(const char* path, const char* tag) {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds(path, clouds)) { std::printf("%-16s 读不到\n", tag); return; }
    Clusterer cl;
    size_t ncl = 0, npass = 0;
    std::vector<double> sizes, res, npts;
    for (const auto& cloud : clouds) {
        if (!cloud.points.empty()) npts.push_back((double)cloud.points.size());
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        for (const auto& c : cs) {
            ++ncl;
            sizes.push_back((double)c.points.size());
            if (c.points.size() < 12) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (!(e > 0)) continue;
            ++npass;
            res.push_back(std::log10(r.dEdx / e));
        }
    }
    double rms = 0;
    if (!res.empty()) {
        const double m = Med(res);
        double s2 = 0;
        for (double v : res) s2 += (v - m) * (v - m);
        rms = std::sqrt(s2 / res.size());
    }
    std::printf("%-16s 云=%5zu 点/云=%5.1f 簇=%6zu 中位=%5.1f 估计=%5zu RMS=%.3f\n", tag,
                clouds.size(), Med(npts), ncl, Med(sizes), npass, rms);
}

void full_check() {
    One("/tmp/abB_cloud.root", "默认(template)");
    One("/tmp/full_p3.root", "P3(局部极大)");
    One("/tmp/full_p1.root", "P1(resid20)");
}
