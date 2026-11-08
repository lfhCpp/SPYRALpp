// 阶段 1+2+3 一览：对每个点云文件打印 点数/云、z 跨度、簇数、估计数、Bethe 残差 RMS
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static void One(const char* path, const char* tag) {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds(path, clouds)) { std::printf("%-18s 读不到\n", tag); return; }
    std::vector<double> npts, zspan;
    Clusterer cl;
    size_t ncl = 0, npass = 0, ntot = 0;
    std::vector<double> sizes, res;
    for (const auto& cloud : clouds) {
        if (cloud.points.empty()) continue;
        npts.push_back((double)cloud.points.size());
        std::vector<double> z;
        for (const auto& p : cloud.points) z.push_back(p.z);
        std::sort(z.begin(), z.end());
        zspan.push_back(z.back() - z.front());
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
        const double m = Med(res);
        double s2 = 0;
        for (double v : res) s2 += (v - m) * (v - m);
        rms = std::sqrt(s2 / res.size());
    }
    std::printf("%-18s 云=%5zu 点/云=%5.1f z跨度=%6.0f 簇=%5zu 中位=%5.1f 估计=%5zu RMS=%.3f\n",
                tag, clouds.size(), Med(npts), Med(zspan), ncl, Med(sizes), npass, rms);
}

void stage1() {
    std::printf("%-18s %5s %5s %6s %5s %5s %5s %5s\n", "配置", "云", "点/云", "z跨度", "簇",
                "中位", "估计", "RMS");
    One("/tmp/adc40.root", "默认 (thr40)");
    One("/tmp/p1.root", "P1 residual20");
    One("/tmp/p2.root", "P2 激进");
    One("/tmp/p3.root", "P3 局部极大");
    std::printf("Spyral 目标：点/云 87，z跨度 1497\n");
}
