// 扫 beam_region_radius：看通过率/RMS 随它的变化（决定默认值）
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

static void Score(const std::vector<PointCloud>& clouds, double radius, int& npass, double& rms,
                  double& frac15, int& ntot) {
    Clusterer clusterer;
    std::vector<double> res;
    npass = ntot = 0;
    int in15 = 0;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& cl : clusters) {
            if (cl.points.size() < 10) continue;
            ++ntot;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(radius);
            EstimateResult r;
            if (!est.Estimate(cl, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (!(e > 0)) continue;
            ++npass;
            const double rr = r.dEdx / e;
            res.push_back(std::log10(rr));
            if (std::fabs(rr - 1.0) < 0.15) ++in15;
        }
    }
    if (res.empty()) { rms = 0; frac15 = 0; return; }
    std::sort(res.begin(), res.end());
    const double med = res[res.size() / 2];
    double s2 = 0;
    for (double v : res) s2 += (v - med) * (v - med);
    rms = std::sqrt(s2 / res.size());
    frac15 = 100.0 * in15 / res.size();
}

void pid_radius_sweep() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) return;
    std::printf("%8s %8s %10s %10s %8s\n", "radius", "pass", "pass%", "RMS(dex)", "±15%");
    for (double rad : {0.0, 10.0, 20.0, 40.0, 80.0, 150.0}) {
        int npass = 0, ntot = 0;
        double rms = 0, f15 = 0;
        Score(clouds, rad, npass, rms, f15, ntot);
        std::printf("%8.0f %8d %9.1f%% %10.3f %7.0f%%\n", rad, npass,
                    100.0 * npass / std::max(1, ntot), rms, f15);
    }
}
