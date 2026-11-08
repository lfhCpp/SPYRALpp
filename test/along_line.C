// 正确的宽度指标：沿参考线（Bethe）的残差，而不是 Bρ 分箱内的垂直展宽
//   r = log10( dEdx / (k·Bethe(Bρ)) )，k = median(dEdx/Bethe)
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
static double Pct(std::vector<double> v, double f) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}

static void Score(const char* path, int ncols, const char* tag) {
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
    if (brho.empty()) { std::printf("%-22s 空\n", tag); return; }
    std::vector<double> ratio;
    for (size_t i = 0; i < brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (e > 0) ratio.push_back(dedx[i] / e);
    }
    const double k = 1.0 / std::max(1e-30, Med(ratio));
    std::vector<double> res;
    int in15 = 0, in25 = 0;
    for (size_t i = 0; i < brho.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (!(e > 0)) continue;
        const double rr = dedx[i] / (k * e);           // = 1 表示落在参考线上
        res.push_back(std::log10(rr));
        if (std::fabs(rr - 1.0) < 0.15) ++in15;
        if (std::fabs(rr - 1.0) < 0.25) ++in25;
    }
    const double m = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - m) * (v - m);
    const double rms = std::sqrt(s2 / res.size());
    const double iqr = Pct(res, 0.75) - Pct(res, 0.25);
    std::printf("%-22s N=%4zu  沿参考线: RMS=%.3f dex  IQR=%.3f dex  ±15%%=%.0f%%  ±25%%=%.0f%%\n",
                tag, res.size(), rms, iqr, 100.0 * in15 / res.size(), 100.0 * in25 / res.size());
}

void along_line() {
    Score("/tmp/ours_est_sklearn.csv", 4, "SpyralPP (C++)");
    Score("/tmp/spy_est_ourclusters.csv", 3, "Spyral (Python)");
}
