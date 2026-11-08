// 阶段 2+3 联合扫描：聚类参数 → 簇统计 + 估计阶段 RMS（目标：簇数≈150、中位≈49、RMS≤0.286）
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

static void Run(const std::vector<PointCloud>& clouds, const ClusterParams& cp, const char* tag) {
    Clusterer cl;
    cl.SetParams(cp);
    size_t ncl = 0, ntot = 0, npass = 0;
    std::vector<double> sizes, res;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        for (const auto& c : cs) {
            ++ncl;
            sizes.push_back((double)c.points.size());
            if (c.points.size() < 12) continue;
            ++ntot;
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
        const double med = Med(res);
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        rms = std::sqrt(s2 / res.size());
    }
    std::printf("%-28s 簇=%4zu 中位=%5.1f 簇/云=%.2f  估计=%5zu/%5zu  通过=%4.0f%%  RMS=%.3f\n",
                tag, ncl, Med(sizes), (double)ncl / std::max<size_t>(1, clouds.size()), npass,
                ntot, 100.0 * npass / std::max<size_t>(1, ntot), rms);
}

void stage23() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/abB_cloud.root", clouds)) return;
    std::printf("点云 %zu | Spyral 目标：150 簇 / 中位 49 / 1.23 簇每云 | 现状 RMS=0.286\n\n",
                clouds.size());

    struct Cfg { double eps; int mp; int ms; const char* tag; };
    const Cfg cfgs[] = {
        {20, 3, 50, "eps20 mpts3 size50 (现状)"},
        {20, 5, 20, "eps20 mpts5 size20"},
        {12, 4, 25, "eps12 mpts4 size25"},
        {10, 4, 20, "eps10 mpts4 size20"},
        {10, 6, 20, "eps10 mpts6 size20"},
        {8, 4, 15, "eps8  mpts4 size15"},
    };
    for (const auto& c : cfgs) {
        ClusterParams cp;
        cp.eps_mm = c.eps;
        cp.min_points = c.mp;
        cp.min_cluster_size = c.ms;
        Run(clouds, cp, c.tag);
    }
}
