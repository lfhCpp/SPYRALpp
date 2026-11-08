// 阶段 2：关掉"端点合并(join)"看簇是否变小、RMS 是否保持
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

static void Run(const std::vector<PointCloud>& clouds, bool join, const char* tag) {
    Clusterer cl;
    ClusterParams cp;
    cp.enable_join = join;
    cl.SetParams(cp);
    size_t ncl = 0;
    std::vector<double> sizes, res;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        for (const auto& c : cs) {
            ++ncl;
            sizes.push_back((double)c.points.size());
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (e > 0) res.push_back(std::log10(r.dEdx / e));
        }
    }
    const double m = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - m) * (v - m);
    std::printf("%-22s 簇=%5zu 中位=%5.1f 簇/云=%.2f 估计=%4zu RMS=%.3f\n", tag, ncl, Med(sizes),
                (double)ncl / std::max<size_t>(1, clouds.size()), res.size(),
                res.empty() ? 0.0 : std::sqrt(s2 / res.size()));
}

void join_test() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) return;
    Run(clouds, true, "join=on (现状)");
    Run(clouds, false, "join=off");
}
