// 残差随 Bρ 的漂移：会不会"斜着离开参考线"
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static void Drift(const char* path, int ncols, const char* tag) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    std::vector<double> brho, dedx;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if ((int)v.size() < ncols) continue;
        if (ncols == 4) { brho.push_back(v[2]); dedx.push_back(v[3]); }
        else { brho.push_back(v[0]); dedx.push_back(v[1]); }
    }
    if (brho.empty()) return;
    std::vector<double> ratio;
    for (size_t i = 0; i < brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (e > 0) ratio.push_back(dedx[i] / e);
    }
    const double k = 1.0 / std::max(1e-30, Med(ratio));
    std::printf("%s 残差中位（按 Bρ 分箱，单位 dex）：\n", tag);
    const double edges[] = {0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.6};
    for (int i = 0; i + 1 < 10; ++i) {
        std::vector<double> r;
        for (size_t j = 0; j < brho.size(); ++j) {
            if (brho[j] < edges[i] || brho[j] >= edges[i + 1]) continue;
            const double e = ParticleId::ExpectedDedx(brho[j], 1, 1);
            if (e > 0) r.push_back(std::log10(dedx[j] / (k * e)));
        }
        if (r.size() >= 15)
            std::printf("   Bρ=%.2f-%.2f n=%4zu  median=%+.2f dex\n", edges[i], edges[i + 1],
                        r.size(), Med(r));
    }
}

void drift() {
    Drift("/tmp/ours_est_sklearn.csv", 4, "SpyralPP");
    Drift("/tmp/spy_est_ourclusters.csv", 3, "Spyral");
}
