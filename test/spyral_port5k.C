// 忠实移植 Spyral 的 estimate_physics_pass（逐行对照 estimator.py 232-330 行）
// 目的：先证明"能复现"，再谈改进
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "reco/Estimator.h"
#include "reco/Geometry.h"
#include "reco/Spline.h"
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

struct P { double x, y, z, q; };

static bool Port(const std::vector<P>& raw, double B, EstimateResult& out, double& L) {
    std::vector<P> d = raw;
    if (d.size() < 10) return false;
    std::sort(d.begin(), d.end(), [](const P& a, const P& b) { return a.z < b.z; });
    const size_t n = d.size();

    // 1) vertex = 首点
    double vx = d[0].x, vy = d[0].y;
    // 2) first arc = 到 rho 最大处
    auto max_rho_idx = [&](double cx0, double cy0) {
        double best = -1.0; size_t idx = 0;
        for (size_t i = 1; i < n; ++i) {
            const double r = std::hypot(d[i].x - cx0, d[i].y - cy0);
            if (r > best) { best = r; idx = i; }
        }
        return idx;
    };
    size_t maximum = max_rho_idx(vx, vy);
    // 3) 首弧圆拟合
    std::vector<double> xs, ys;
    for (size_t i = 0; i <= maximum; ++i) { xs.push_back(d[i].x); ys.push_back(d[i].y); }
    double cx = 0, cy = 0, R = 0;
    if (!FitCircle(xs, ys, cx, cy, R)) return false;
    // 4) 顶点 = 圆上离束流轴最近的点
    const double rho_c = std::hypot(cx, cy);
    if (rho_c > 1e-9) { vx = cx - R * cx / rho_c; vy = cy - R * cy / rho_c; }
    // 5) 用新顶点重选 first arc
    maximum = max_rho_idx(vx, vy);
    if (maximum < 10) return false;
    // 6) polar：rho vs z 回归，取前 test_index = max(10, maximum/2) 点
    const size_t test_index = std::max<size_t>(10, maximum / 2);
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 0; i < test_index; ++i) {
        const double zz = d[i].z;
        const double rr = std::hypot(d[i].x - vx, d[i].y - vy);
        sx += zz; sy += rr; sxx += zz * zz; sxy += zz * rr;
    }
    const double den = test_index * sxx - sx * sx;
    if (std::fabs(den) < 1e-9) return false;
    const double slope = (test_index * sxy - sx * sy) / den;
    if (std::fabs(slope) < 1e-9) return false;
    double polar = std::atan(slope);
    if (polar < 0) polar += M_PI;
    // 7) 电荷：first arc 累加，到 |(x,y)|>152 截断
    double charge = d[0].q;
    size_t cutoff = maximum;
    for (size_t i = 0; i + 1 <= maximum; ++i) {
        if (std::hypot(d[i + 1].x, d[i + 1].y) > 152.0) { cutoff = i + 1; break; }
        charge += d[i + 1].q;
    }
    // 8) 弧长：整条簇拟合 x(z),y(z) 平滑样条，在 [z0, z_cut] 上采 1000 点积分
    std::vector<double> zs, sxv, syv;
    for (const auto& p : d) { zs.push_back(p.z); sxv.push_back(p.x); syv.push_back(p.y); }
    SmoothingSpline spx, spy;
    const int nk = std::max(8, std::min(400, (int)zs.size()));  // 节点≈数据点
    const bool okx = spx.Fit(zs, sxv, 300.0, nk), oky = spy.Fit(zs, syv, 300.0, nk);
    if (!okx || !oky) return false;
    const double za = d[0].z, zb = d[cutoff].z;
    double len = 0, xp = spx.Eval(za), yp = spy.Eval(za);
    const int NS = 1000;
    for (int i = 1; i <= NS; ++i) {
        const double z = za + (zb - za) * i / NS;
        const double x = spx.Eval(z), y = spy.Eval(z);
        len += std::sqrt((x - xp) * (x - xp) + (y - yp) * (y - yp) + (z - (za + (zb - za) * (i - 1) / NS)) * (z - (za + (zb - za) * (i - 1) / NS)));
        xp = x; yp = y;
    }
    if (!(len > 1.0)) return false;
    out.brho = B * (R * 0.001) / std::max(1e-6, std::fabs(std::sin(polar)));
    out.polar = polar;
    out.dEdx = charge / len;   // Spyral 原式：总电荷 / 总弧长
    out.arclength = len;
    L = len;
    return true;
}

void spyral_port5k() {
    std::map<std::pair<int, int>, std::vector<P>> groups;
    std::ifstream in("/tmp/our_clusters5k.csv");
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string f;
        std::vector<double> v;
        while (std::getline(ss, f, ',')) v.push_back(std::atof(f.c_str()));
        if (v.size() < 6) continue;
        groups[{(int)v[0], (int)v[1]}].push_back(P{v[2], v[3], v[4], v[5]});
    }

    // Spyral 自己的估计（配对用）
    std::map<std::pair<int, int>, double> spy_brho, spy_dedx, spy_pol;
    {
        std::ifstream f2("/tmp/spy_est5k.csv");
        std::string l2;
        std::getline(f2, l2);
        while (std::getline(f2, l2)) {
            if (l2.empty()) continue;
            std::stringstream ss(l2);
            std::string g; std::vector<double> v;
            while (std::getline(ss, g, ',')) v.push_back(std::atof(g.c_str()));
            if (v.size() < 5) continue;
            const auto k = std::make_pair((int)v[0], (int)v[1]);
            spy_brho[k] = v[2]; spy_dedx[k] = v[3]; spy_pol[k] = v[4];
        }
    }

    std::vector<double> brho, dedx, pol, rb, rd;
    FILE* fc = std::fopen("/tmp/port_est5k.csv", "w");
    std::fprintf(fc, "event,cluster,brho,dedx,polar_deg\n");
    for (auto& kv : groups) {
        if (kv.second.size() < 30) continue;
        EstimateResult r;
        double L = 0;
        if (!Port(kv.second, 2.85, r, L)) continue;
        brho.push_back(r.brho); dedx.push_back(r.dEdx); pol.push_back(r.polar * 180 / M_PI);
        std::fprintf(fc, "%d,%d,%.6f,%.6f,%.3f\n", kv.first.first, kv.first.second,
                     r.brho, r.dEdx, r.polar * 180 / M_PI);
        auto it = spy_brho.find(kv.first);
        if (it != spy_brho.end()) {
            rb.push_back(r.brho / it->second);
            rd.push_back(r.dEdx / spy_dedx[kv.first]);
        }
    }
    auto med = [](std::vector<double> v) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
    };
    std::printf("移植版: N=%zu  Bρ中位=%.3f  dEdx中位=%.0f  polar中位=%.0f°\n", brho.size(),
                med(brho), med(dedx), med(pol));
    std::printf("配对 %zu 个: Bρ比中位=%.3f (p16/84 %.2f/%.2f)  dEdx比中位=%.3f\n", rb.size(),
                med(rb), med(rb), med(rb), med(rd));
    std::printf("Spyral: Bρ中位=0.255  dEdx中位=1427  polar中位=63°\n");
    std::fclose(fc);
}
