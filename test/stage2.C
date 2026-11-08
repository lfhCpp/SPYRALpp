// 阶段 2 对齐：找一组聚类参数，使簇数/簇大小接近 Spyral（150 簇、中位 49 点）
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

static void Run(const std::vector<PointCloud>& clouds, const ClusterParams& cp, const char* tag) {
    Clusterer cl;
    cl.SetParams(cp);
    size_t ncl = 0, nev = 0;
    std::vector<double> sizes;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        if (!cs.empty()) ++nev;
        for (const auto& c : cs) { ++ncl; sizes.push_back((double)c.points.size()); }
    }
    std::printf("%-34s 簇=%4zu 事件=%3zu/%zu 大小中位=%5.1f 最大=%5.0f 簇/云=%.2f\n", tag, ncl,
                nev, clouds.size(), Median(sizes), sizes.empty() ? 0.0 : *std::max_element(sizes.begin(), sizes.end()),
                (double)ncl / std::max<size_t>(1, clouds.size()));
}

void stage2() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_300_cloud.root", clouds)) return;
    std::printf("点云 %zu 个（Spyral 目标：150 簇、中位 49 点、1.23 簇/云）\n\n", clouds.size());

    ClusterParams a;
    a.eps_mm = 20; a.min_points = 3; a.min_cluster_size = 50;
    Run(clouds, a, "DBSCAN eps20 size50 (现状)");

    ClusterParams b = a;
    b.eps_mm = 10; b.min_cluster_size = 20;
    Run(clouds, b, "DBSCAN eps10 size20");

    ClusterParams c = a;
    c.use_hdbscan = true;
    c.hdbscan_min_samples = 3;
    c.min_cluster_size = 10;
    Run(clouds, c, "HDBSCAN mpts3 size10");

    ClusterParams d = c;
    d.min_cluster_size = 20;
    Run(clouds, d, "HDBSCAN mpts3 size20");

    ClusterParams e = c;
    e.min_cluster_size = 50;
    e.hdbscan_min_samples = 5;
    Run(clouds, e, "HDBSCAN mpts5 size50");

    ClusterParams f = c;
    f.hdbscan_min_samples = 2;
    f.min_cluster_size = 5;
    Run(clouds, f, "HDBSCAN mpts2 size5");
}
