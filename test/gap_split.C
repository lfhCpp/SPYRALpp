// 沿弧长找"真空隙"再拆分（只切真正断开的地方，不劈开连续径迹）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// 在弧长上找内部空隙（连续 bins 无点），有则拆分
static bool SplitAtGap(std::vector<CloudPoint> pts, std::vector<CloudPoint>& a,
                       std::vector<CloudPoint>& b, double bin_mm, int gap_bins, int min_side) {
    if ((int)pts.size() < 2 * min_side) return false;
    std::sort(pts.begin(), pts.end(), [](const CloudPoint& x, const CloudPoint& y) { return x.z < y.z; });
    std::vector<double> s(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) {
        const double dx = pts[i].x - pts[i - 1].x;
        const double dy = pts[i].y - pts[i - 1].y;
        const double dz = pts[i].z - pts[i - 1].z;
        s[i] = s[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    const double L = s.back();
    if (!(L > 20.0)) return false;
    const int nb = std::max(1, (int)std::ceil(L / bin_mm));
    std::vector<int> cnt(nb, 0);
    for (size_t i = 0; i < pts.size(); ++i) {
        int ib = std::min(nb - 1, (int)(s[i] / bin_mm));
        ++cnt[ib];
    }
    // 找内部连续 gap_bins 个空 bin
    int run = 0;
    int best_start = -1, best_len = 0;
    for (int i = 1; i + 1 < nb; ++i) {
        if (cnt[i] == 0) {
            if (run == 0) best_start = (best_start < 0 ? i : best_start);
            ++run;
            if (run > best_len) best_len = run;
        } else {
            if (run >= gap_bins && best_len == run) best_start = i - run;
            run = 0;
        }
    }
    if (best_len < gap_bins) return false;
    const double s_cut = (best_start + 0.5 * best_len) * bin_mm;
    for (size_t i = 0; i < pts.size(); ++i) (s[i] < s_cut ? a : b).push_back(pts[i]);
    if ((int)a.size() < min_side || (int)b.size() < min_side) { a.clear(); b.clear(); return false; }
    return true;
}

void gap_split() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) return;

    for (int mode = 0; mode < 2; ++mode) {
        Clusterer cl;
        size_t ncl = 0;
        std::vector<double> sizes, res;
        for (const auto& cloud : clouds) {
            std::vector<Cluster> cs;
            cl.Run(cloud, cs);
            std::deque<Cluster> work(cs.begin(), cs.end());
            std::vector<Cluster> out;
            if (mode == 1) {
                while (!work.empty()) {
                    Cluster c = work.front();
                    work.pop_front();
                    std::vector<CloudPoint> a, b;
                    if ((int)c.points.size() >= 60 && SplitAtGap(c.points, a, b, 10.0, 2, 15)) {
                        Cluster ca;
                        ca.direction = c.direction;
                        ca.points = a;
                        Cluster cb;
                        cb.direction = c.direction;
                        cb.points = b;
                        work.push_back(ca);
                        work.push_back(cb);
                    } else {
                        out.push_back(c);
                    }
                }
            } else {
                out = cs;
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
        std::printf("%s 簇=%zu 中位=%.1f 簇/云=%.2f 估计=%zu RMS=%.3f\n",
                    mode == 0 ? "不拆   " : "弧长空隙拆", ncl, Med(sizes),
                    (double)ncl / std::max<size_t>(1, clouds.size()), res.size(),
                    res.empty() ? 0.0 : std::sqrt(s2 / res.size()));
    }
}
