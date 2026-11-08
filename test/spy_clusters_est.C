// 用 Spyral 的簇点集（同一 2000 事件）跑我们的 Estimator，看 Bρ 是否仍偏小
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

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

void spy_clusters_est() {
    std::ifstream in("/tmp/spy_clusters.csv");
    std::string line;
    std::getline(in, line); // header
    std::vector<double> brho, dedx, polar, res;
    std::vector<CloudPoint> cur;
    int cur_ev = -1, cur_cl = -1;
    int ncl = 0;

    auto flush = [&]() {
        if (cur.size() < 30) { cur.clear(); return; }
        ++ncl;
        Cluster c;
        c.points = cur;
        c.direction = Direction::Forward;
        Estimator est;
        est.SetSmallPadRadius(152.0);
        est.SetBeamRegionRadius(1000.0); // 这次只比较几何/电荷，不做束流轴切割
        EstimateResult r;
        if (est.Estimate(c, PointCloud(), 2.85, r)) {
            if (r.brho > 0 && r.dEdx > 0) {
                brho.push_back(r.brho);
                dedx.push_back(r.dEdx);
                polar.push_back(r.polar * 180.0 / M_PI);
                const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
                if (e > 0) res.push_back(std::log10(r.dEdx / e));
            }
        }
        cur.clear();
    };

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 6) continue;
        const int ev = (int)v[0], cl = (int)v[1];
        if (ev != cur_ev || cl != cur_cl) {
            flush();
            cur_ev = ev;
            cur_cl = cl;
        }
        CloudPoint p;
        p.x = v[2]; p.y = v[3]; p.z = v[4]; p.integral = v[5];
        cur.push_back(p);
    }
    flush();

    const double med = Med(res);
    double s2 = 0;
    for (double x : res) s2 += (x - med) * (x - med);
    std::printf("Spyral 的簇: npts>=30 共 %d 个, 我们的估计成功 %zu 个\n", ncl, brho.size());
    if (!brho.empty())
        std::printf("我们跑在 Spyral 簇上：Bρ中位=%.3f  dEdx中位=%.0f  polar中位=%.0f°  RMS=%.3f\n",
                    Med(brho), Med(dedx), Med(polar), std::sqrt(s2 / std::max<size_t>(1, res.size())));
    std::printf("对照：Spyral 自己的估计 Bρ中位=0.255、polar=63°；我们自己的簇 Bρ中位=0.154\n");
}
