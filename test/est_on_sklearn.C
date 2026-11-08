// 用 sklearn HDBSCAN（Spyral 参数）切出的簇，跑我们的估计器
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/Hdbscan.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

void est_on_sklearn() {
    std::ifstream in("/tmp/our_clusters_sklearn.csv");
    std::string line;
    std::getline(in, line);
    std::map<std::pair<int, int>, std::vector<CloudPoint>> groups;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 6) continue;
        CloudPoint p;
        p.x = v[2]; p.y = v[3]; p.z = v[4]; p.integral = v[5];
        groups[{(int)v[0], (int)v[1]}].push_back(p);
    }
    std::printf("簇 %zu 个\n", groups.size());

    std::vector<double> brho, dedx, polar, res, sizes;
    FILE* out = std::fopen("/tmp/ours_est_sklearn.csv", "w");
    std::fprintf(out, "event,cluster,brho,dedx,polar_deg,arclength,npts\n");
    for (auto& kv : groups) {
        const auto& pts = kv.second;
        sizes.push_back((double)pts.size());
        if ((int)pts.size() < 30) continue;
        Cluster c;
        c.direction = Direction::Forward;
        c.points = pts;
        Estimator est;
        est.SetSmallPadRadius(152.0);
        est.SetBeamRegionRadius(25.0);
        EstimateResult r;
        if (!est.Estimate(c, PointCloud(), 2.85, r)) continue;
        if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
        if (r.arclength < 40.0 || r.arclength > 800.0) continue;  // 质量切割：去掉病态弧长

        if (r.polar * 180.0 / M_PI < 30.0) continue;  // 束流/前向本底
        if (r.polar * 180.0 / M_PI < 20.0) continue;  // 束流/前向本底
        brho.push_back(r.brho);
        dedx.push_back(r.dEdx);
        polar.push_back(r.polar * 180.0 / M_PI);
        const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
        if (e > 0) res.push_back(std::log10(r.dEdx / e));
        std::fprintf(out, "%d,%d,%.6f,%.6f,%.3f,%.3f,%d\n", kv.first.first, kv.first.second,
                     r.brho, r.dEdx, r.polar * 180.0 / M_PI, r.arclength, (int)pts.size());
    }
    std::fclose(out);

    const double m = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - m) * (v - m);
    std::printf("簇大小中位=%.1f | 估计=%zu | Bρ中位=%.3f | dEdx中位=%.0f | polar中位=%.0f° | RMS=%.3f\n",
                Med(sizes), brho.size(), Med(brho), Med(dedx), Med(polar),
                res.empty() ? 0.0 : std::sqrt(s2 / res.size()));
    std::printf("对照 Spyral: Bρ 0.255, polar 63°, dEdx 1427\n");
}
