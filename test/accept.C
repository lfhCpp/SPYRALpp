// 在 Spyral 的簇上扫 beam_region_radius / 最小点数，找"多收簇且 RMS 不劣化"的工作点
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace rk4;

struct Cl {
    std::vector<CloudPoint> pts;
    int event = -1, index = 0;
};

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static void Score(const std::vector<Cl>& cls, double beam, int minpts, const char* tag) {
    std::vector<double> brho, res;
    for (const auto& c : cls) {
        if ((int)c.pts.size() < minpts) continue;
        Cluster cl;
        cl.direction = Direction::Forward;
        cl.points = c.pts;
        Estimator est;
        est.SetSmallPadRadius(152.0);
        est.SetBeamRegionRadius(beam);
        EstimateResult r;
        if (!est.Estimate(cl, PointCloud(), 2.85, r)) continue;
        if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
        brho.push_back(r.brho);
        const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
        if (e > 0) res.push_back(std::log10(r.dEdx / e));
    }
    const double m = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - m) * (v - m);
    std::printf("%-34s 估计=%4zu Bρ中位=%.3f RMS=%.3f\n", tag, brho.size(), Med(brho),
                res.empty() ? 0.0 : std::sqrt(s2 / res.size()));
}

void accept() {
    std::vector<Cl> cls;
    std::ifstream in("/tmp/spy_clusters.csv");
    std::string line;
    std::getline(in, line);
    Cl cur;
    int ev = -1, idx = 0;
    auto flush = [&]() {
        if (cur.pts.size() >= 10) cls.push_back(cur);
        cur.pts.clear();
    };
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 6) continue;
        const int e2 = (int)v[0], c2 = (int)v[1];
        if (e2 != ev || c2 != idx) { flush(); ev = e2; idx = c2; }
        CloudPoint p;
        p.x = v[2]; p.y = v[3]; p.z = v[4]; p.integral = v[5];
        cur.pts.push_back(p);
    }
    flush();
    std::printf("Spyral 的簇共 %zu 个（≥10 点）\n", cls.size());
    for (double b : {25.0, 40.0, 60.0, 120.0}) {
        char tag[64];
        std::snprintf(tag, sizeof(tag), "beam=%.0fmm minpts=30", b);
        Score(cls, b, 30, tag);
    }
    Score(cls, 25.0, 20, "beam=25mm minpts=20");
    Score(cls, 40.0, 20, "beam=40mm minpts=20");
}
