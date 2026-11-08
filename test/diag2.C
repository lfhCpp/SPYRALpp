// 残差 vs 弧长/点数/极角：定位两个尾巴的成因
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

void diag2() {
    std::ifstream in("/tmp/ours_est_sklearn.csv");
    std::string line;
    std::getline(in, line);
    std::vector<double> brho, dedx, polar;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 5) continue;
        brho.push_back(v[2]);
        dedx.push_back(v[3]);
        polar.push_back(v[4]);
    }
    std::vector<double> ratio;
    for (size_t i = 0; i < brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (e > 0) ratio.push_back(dedx[i] / e);
    }
    const double k = 1.0 / std::max(1e-30, Med(ratio));
    std::printf("N=%zu\n", brho.size());

    // 按极角分箱
    const double pe[] = {0, 20, 40, 60, 70, 80, 85, 90, 180};
    std::printf("残差 vs 极角：\n");
    for (int i = 0; i + 1 < 9; ++i) {
        std::vector<double> r;
        for (size_t j = 0; j < brho.size(); ++j) {
            if (polar[j] < pe[i] || polar[j] >= pe[i + 1]) continue;
            const double e = ParticleId::ExpectedDedx(brho[j], 1, 1);
            if (e > 0) r.push_back(std::log10(dedx[j] / (k * e)));
        }
        if (r.size() >= 15)
            std::printf("   theta=%2.0f-%3.0f n=%4zu median=%+.2f dex\n", pe[i], pe[i + 1],
                        r.size(), Med(r));
    }
    // 按 Bρ 与极角的二维
    std::printf("低 Bρ 端的极角构成：\n");
    for (int i = 0; i + 1 < 8; ++i) {
        int n = 0;
        std::vector<double> r;
        for (size_t j = 0; j < brho.size(); ++j) {
            if (brho[j] >= 0.05 && brho[j] < 0.15 && polar[j] >= pe[i] && polar[j] < pe[i + 1]) {
                ++n;
                const double e = ParticleId::ExpectedDedx(brho[j], 1, 1);
                if (e > 0) r.push_back(std::log10(dedx[j] / (k * e)));
            }
        }
        if (n > 3) std::printf("   theta=%2.0f-%3.0f n=%3d median=%+.2f\n", pe[i], pe[i + 1], n,
                               Med(r));
    }
}
