// 同一极角切割下比较两软件沿参考线的宽度
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
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
static double Pct(std::vector<double> v, double f) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}

static void Score(const char* path, bool ours, double polar_cut, const char* tag) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    std::vector<double> brho, dedx, pol;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (ours) {
            if (v.size() < 5) continue;
            brho.push_back(v[2]); dedx.push_back(v[3]); pol.push_back(v[4]);
        } else {
            if (v.size() < 3) continue;
            brho.push_back(v[0]); dedx.push_back(v[1]); pol.push_back(v[2]);
        }
    }
    std::vector<double> ratio;
    for (size_t i = 0; i < brho.size(); ++i) {
        if (pol[i] < polar_cut) continue;
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (e > 0) ratio.push_back(dedx[i] / e);
    }
    if (ratio.empty()) { std::printf("%-28s 空\n", tag); return; }
    const double k = 1.0 / std::max(1e-30, Med(ratio));
    std::vector<double> res;
    int in15 = 0, in25 = 0;
    for (size_t i = 0; i < brho.size(); ++i) {
        if (pol[i] < polar_cut) continue;
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (!(e > 0)) continue;
        const double rr = dedx[i] / (k * e);
        res.push_back(std::log10(rr));
        if (std::fabs(rr - 1) < 0.15) ++in15;
        if (std::fabs(rr - 1) < 0.25) ++in25;
    }
    const double m = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - m) * (v - m);
    std::printf("%-28s N=%4zu RMS=%.3f IQR=%.3f ±15%%=%.0f%% ±25%%=%.0f%%\n", tag, res.size(),
                std::sqrt(s2 / res.size()), Pct(res, 0.75) - Pct(res, 0.25),
                100.0 * in15 / res.size(), 100.0 * in25 / res.size());
}

void cut_compare() {
    for (double th : {0.0, 20.0, 30.0}) {
        char t1[64], t2[64];
        std::snprintf(t1, sizeof(t1), "SpyralPP theta>%.0f", th);
        std::snprintf(t2, sizeof(t2), "Spyral   theta>%.0f", th);
        Score("/tmp/ours_est_sklearn.csv", true, th, t1);
        Score("/tmp/spy_est_ourclusters.csv", false, th, t2);
    }
}
