// 用【默认】ClusterParams（不显式设置）验证新的默认值
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

void defaults_check() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/abB_cloud.root", clouds)) return;
    Clusterer cl; // 默认参数
    ClusterParams cp; // 读取默认值用于打印
    size_t ncl = 0, npass = 0, ntot = 0;
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
        std::sort(res.begin(), res.end());
        const double med = res[res.size() / 2];
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        rms = std::sqrt(s2 / res.size());
    }
    std::sort(sizes.begin(), sizes.end());
    std::printf("默认参数: eps=%.0f mpts=%d size=%d → 簇=%zu 中位=%.1f 估计=%zu RMS=%.3f\n",
                cp.eps_mm, cp.min_points, cp.min_cluster_size, ncl,
                sizes.empty() ? 0.0 : sizes[sizes.size() / 2], npass, rms);
}
