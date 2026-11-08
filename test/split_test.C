// 试验：对过大的簇按主轴做 2-means 拆分，看能否把中位簇大小降下来而不伤 RMS
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// 沿主轴 1D 2-means 拆分；返回是否拆分成功
static bool TrySplit(const std::vector<CloudPoint>& pts, std::vector<CloudPoint>& a,
                     std::vector<CloudPoint>& b) {
    const size_t n = pts.size();
    if (n < 60) return false;
    double mx = 0, my = 0, mz = 0;
    for (const auto& p : pts) { mx += p.x; my += p.y; mz += p.z; }
    mx /= n; my /= n; mz /= n;
    double cxx = 0, cxy = 0, cxz = 0, cyy = 0, cyz = 0, czz = 0;
    for (const auto& p : pts) {
        const double dx = p.x - mx, dy = p.y - my, dz = p.z - mz;
        cxx += dx * dx; cxy += dx * dy; cxz += dx * dz;
        cyy += dy * dy; cyz += dy * dz; czz += dz * dz;
    }
    // 幂迭代求主轴
    double vx = 1, vy = 1, vz = 1;
    for (int it = 0; it < 50; ++it) {
        double nx = cxx * vx + cxy * vy + cxz * vz;
        double ny = cxy * vx + cyy * vy + cyz * vz;
        double nz = cxz * vx + cyz * vy + czz * vz;
        const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (!(norm > 0)) return false;
        vx = nx / norm; vy = ny / norm; vz = nz / norm;
    }
    std::vector<double> t(n);
    for (size_t i = 0; i < n; ++i) t[i] = (pts[i].x - mx) * vx + (pts[i].y - my) * vy + (pts[i].z - mz) * vz;
    std::sort(t.begin(), t.end());
    const double tmin = t.front(), tmax = t.back();
    // 在 t 上做 1D 2-means
    double c1 = tmin, c2 = tmax;
    std::vector<int> lab(n, 0);
    for (int it = 0; it < 30; ++it) {
        double s1 = 0, s2 = 0; int n1 = 0, n2 = 0;
        for (size_t i = 0; i < n; ++i) {
            const double ti = (pts[i].x - mx) * vx + (pts[i].y - my) * vy + (pts[i].z - mz) * vz;
            const bool first = std::fabs(ti - c1) <= std::fabs(ti - c2);
            lab[i] = first ? 0 : 1;
            if (first) { s1 += ti; ++n1; } else { s2 += ti; ++n2; }
        }
        if (n1 == 0 || n2 == 0) return false;
        c1 = s1 / n1; c2 = s2 / n2;
    }
    // 判据：两簇中心距离 vs 平均簇内展宽
    double w1 = 0, w2 = 0; int n1 = 0, n2 = 0;
    for (size_t i = 0; i < n; ++i) {
        const double ti = (pts[i].x - mx) * vx + (pts[i].y - my) * vy + (pts[i].z - mz) * vz;
        if (lab[i] == 0) { w1 += (ti - c1) * (ti - c1); ++n1; }
        else { w2 += (ti - c2) * (ti - c2); ++n2; }
    }
    if (n1 < 20 || n2 < 20) return false;
    const double within = std::sqrt((w1 + w2) / (n1 + n2));
    if (within < 1e-6) return false;
    if (std::fabs(c2 - c1) / within < 2.0) return false; // 分不开就不拆
    for (size_t i = 0; i < n; ++i) (lab[i] == 0 ? a : b).push_back(pts[i]);
    return true;
}

void split_test() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) return;
    for (int split_pass = 0; split_pass < 2; ++split_pass) {
        Clusterer cl;
        size_t ncl = 0;
        std::vector<double> sizes, res;
        for (const auto& cloud : clouds) {
            std::vector<Cluster> cs;
            cl.Run(cloud, cs);
            std::vector<Cluster> out;
            for (const auto& c : cs) {
                std::vector<CloudPoint> a, b;
                if (split_pass == 1 && TrySplit(c.points, a, b)) {
                    Cluster ca;
                    ca.direction = c.direction;
                    ca.points = a;
                    out.push_back(ca);
                    Cluster cb;
                    cb.direction = c.direction;
                    cb.points = b;
                    out.push_back(cb);
                } else {
                    out.push_back(c);
                }
            }
            for (const auto& c : out) {
                ++ncl;
                sizes.push_back((double)c.points.size());
                if (c.points.size() < 30) continue;
                Estimator est;
                est.SetSmallPadRadius(152.0);
                est.SetBeamRegionRadius(25.0);
                EstimateResult r;
                if (!est.Estimate(c, cloud, 2.85, r)) continue;
                if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
                const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
                if (e > 0) res.push_back(std::log10(r.dEdx / e));
            }
        }
        const double m = Med(res);
        double s2 = 0;
        for (double v : res) s2 += (v - m) * (v - m);
        std::printf("split_pass=%d 簇=%zu 中位=%.1f 估计=%zu RMS=%.3f\n", split_pass, ncl,
                    Med(sizes), res.size(), res.empty() ? 0.0 : std::sqrt(s2 / res.size()));
    }
}
