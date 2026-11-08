// 对比两个点云文件的电荷刻度：/tmp/spy_full_cloud.root（8993 那张图用的）vs /tmp/abB_cloud.root
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include <algorithm>
#include <cstdio>
#include <vector>
using namespace rk4;

static void Stats(const char* path, const char* tag) {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds(path, clouds)) { std::printf("%s: 读不到\n", tag); return; }
    std::vector<double> iv;
    for (const auto& c : clouds)
        for (const auto& p : c.points) iv.push_back(p.integral);
    if (iv.empty()) { std::printf("%s: 空\n", tag); return; }
    std::sort(iv.begin(), iv.end());
    std::printf("%-28s 云=%zu  点=%zu  integral 中位=%.4g  p10=%.4g  p90=%.4g\n", tag,
                clouds.size(), iv.size(), iv[iv.size() / 2], iv[iv.size() / 10],
                iv[iv.size() * 9 / 10]);
}

void chk() {
    Stats("/tmp/spy_full_cloud.root", "spy_full (8993 图用的)");
    Stats("/tmp/abB_cloud.root", "abB (现在模板拟合)");
    Stats("/tmp/abA_cloud.root", "abA (现在局部极大)");
}
