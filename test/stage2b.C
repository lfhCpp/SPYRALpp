// 阶段 2 第二轮：提高密度阈值，看能不能把大点云切细（目标：簇数↑、中位大小↓）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include <algorithm>
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
    size_t ncl = 0;
    std::vector<double> sizes;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        cl.Run(cloud, cs);
        for (const auto& c : cs) { ++ncl; sizes.push_back((double)c.points.size()); }
    }
    std::printf("%-30s 簇=%4zu 中位=%5.1f 最大=%5.0f 簇/云=%.2f\n", tag, ncl, Median(sizes),
                sizes.empty() ? 0.0 : *std::max_element(sizes.begin(), sizes.end()),
                (double)ncl / std::max<size_t>(1, clouds.size()));
}

void stage2b() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_300_cloud.root", clouds)) return;
    std::printf("点云 %zu（Spyral: 150 簇 / 中位 49 / 1.23 簇每云）\n", clouds.size());
    for (int mp : {5, 8, 12, 16}) {
        ClusterParams a;
        a.eps_mm = 20; a.min_points = mp; a.min_cluster_size = 20;
        char tag[64];
        std::snprintf(tag, sizeof(tag), "DBSCAN eps20 mpts%d size20", mp);
        Run(clouds, a, tag);
    }
    for (int ms : {8, 12, 16}) {
        ClusterParams h;
        h.use_hdbscan = true; h.eps_mm = 20; h.min_points = 3;
        h.hdbscan_min_samples = ms; h.min_cluster_size = 10;
        char tag[64];
        std::snprintf(tag, sizeof(tag), "HDBSCAN ms%d size10", ms);
        Run(clouds, h, tag);
    }
    ClusterParams z;
    z.eps_mm = 12; z.min_points = 6; z.min_cluster_size = 20;
    Run(clouds, z, "DBSCAN eps12 mpts6 size20");
}
