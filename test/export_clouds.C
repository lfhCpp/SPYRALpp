// 导出我们的点云到 CSV，供 sklearn HDBSCAN 使用
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include <cstdio>
#include <vector>
using namespace rk4;

void export_clouds() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) { std::printf("no clouds\n"); return; }
    FILE* f = std::fopen("/tmp/our_clouds.csv", "w");
    std::fprintf(f, "event,x,y,z,integral\n");
    size_t n = 0;
    for (const auto& c : clouds) {
        for (const auto& p : c.points) {
            std::fprintf(f, "%lld,%.5f,%.5f,%.5f,%.5f\n", (long long)c.event, p.x, p.y, p.z,
                         p.integral);
            ++n;
        }
    }
    std::fclose(f);
    std::printf("clouds=%zu points=%zu -> /tmp/our_clouds.csv\n", clouds.size(), n);
}
