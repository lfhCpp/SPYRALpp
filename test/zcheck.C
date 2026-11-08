// 我们的点云 z 统计（前 5 个事件）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include <algorithm>
#include <cstdio>
#include <vector>
using namespace rk4;
void zcheck() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_300_cloud.root", clouds)) return;
    for (size_t i = 0; i < clouds.size() && i < 5; ++i) {
        std::vector<double> z;
        for (const auto& p : clouds[i].points) z.push_back(p.z);
        if (z.empty()) continue;
        std::sort(z.begin(), z.end());
        std::printf("ev%-5d n=%4zu z=[%7.1f, %7.1f] 中位=%7.1f\n", clouds[i].number, z.size(),
                    z.front(), z.back(), z[z.size() / 2]);
    }
}
