// 阶段 1（点云）+ 阶段 2（聚类）统计：我方
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
using namespace rk4;

static double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

void stage12() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_300_cloud.root", clouds)) { std::printf("read fail\n"); return; }
    std::vector<double> npts, zspan, integ, clsize;
    int nevt_with_cluster = 0;
    size_t nclusters = 0;
    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20.0;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);

    for (const auto& cloud : clouds) {
        if (cloud.points.empty()) continue;
        npts.push_back((double)cloud.points.size());
        std::vector<double> z, iv;
        for (const auto& p : cloud.points) { z.push_back(p.z); iv.push_back(p.integral); }
        std::sort(z.begin(), z.end());
        zspan.push_back(z.back() - z.front());
        integ.push_back(Median(iv));

        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        if (!clusters.empty()) ++nevt_with_cluster;
        for (const auto& c : clusters) {
            nclusters++;
            clsize.push_back((double)c.points.size());
        }
    }
    std::printf("OURS points: clouds=%zu pts/cloud(med)=%.1f zspan(med)=%.1f mm integ(med)=%.2f\n",
                npts.size(), Median(npts), Median(zspan), Median(integ));
    std::printf("OURS clusters: events_with_cluster=%d/%zu clusters=%zu size(med)=%.1f "
                "clusters/cloud=%.2f\n",
                nevt_with_cluster, npts.size(), nclusters, Median(clsize),
                (double)nclusters / std::max<size_t>(1, npts.size()));
}
