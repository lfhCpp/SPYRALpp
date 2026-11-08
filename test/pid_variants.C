// 估计阶段 dE/dx 变体对比：谁的 Bethe 残差最窄
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

static double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// 统计：Bethe 残差 RMS(dex)、±15% 内比例、±25% 内比例
static void Score(const char* tag, const std::vector<double>& dedx,
                  const std::vector<double>& brho, const std::vector<double>& len) {
    std::vector<double> r;
    int in15 = 0, in25 = 0;
    for (size_t i = 0; i < dedx.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
        if (!(e > 0) || !(dedx[i] > 0)) continue;
        const double ratio = dedx[i] / e;      // 与质子的偏差倍数
        r.push_back(std::log10(ratio));
        if (std::fabs(ratio - 1.0) < 0.15) ++in15;
        if (std::fabs(ratio - 1.0) < 0.25) ++in25;
    }
    if (r.empty()) { std::printf("%-26s 无数据\n", tag); return; }
    const double med = Median(r);
    double s2 = 0;
    for (double v : r) s2 += (v - med) * (v - med);
    const double rms = std::sqrt(s2 / r.size());
    std::printf("%-26s n=%zu  中位=%+.2f dex  RMS=%.3f dex  ±15%%:%.0f%%  ±25%%:%.0f%%\n",
                tag, r.size(), med, rms, 100.0 * in15 / r.size(), 100.0 * in25 / r.size());
    (void)len;
}

// 逐段 dE/dx（按 Δs 分箱，pad 分组：同一 pad 的命中合并成一段足迹）
static void Segments(const std::vector<CloudPoint>& pts, double ds, double trunc,
                     std::vector<double>& out) {
    if (pts.size() < 6) return;
    std::vector<double> s(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) {
        const double dx = pts[i].x - pts[i - 1].x, dy = pts[i].y - pts[i - 1].y;
        const double dz = pts[i].z - pts[i - 1].z;
        s[i] = s[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    const double L = s.back();
    if (!(L > 5.0)) return;
    const int nb = std::max(1, (int)std::ceil(L / ds));
    std::vector<double> q(nb, 0.0);
    std::vector<int> n(nb, 0);
    for (size_t i = 0; i < pts.size(); ++i) {
        int b = std::min(nb - 1, (int)(s[i] / ds));
        q[b] += pts[i].integral;
        ++n[b];
    }
    std::vector<double> v;
    for (int b = 0; b < nb; ++b) {
        if (n[b] < 1) continue;
        const double len = std::min(ds, L - b * ds);
        if (len > 1.0) v.push_back(q[b] / len);
    }
    if (v.empty()) return;
    std::sort(v.begin(), v.end());
    const int keep = std::max(1, (int)(v.size() * (1.0 - trunc)));
    double sum = 0;
    for (int i = 0; i < keep; ++i) sum += v[i];
    out.push_back(sum / keep); // 截断均值
}

void pid_variants() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) { std::printf("读不到\n"); return; }

    Clusterer clusterer;
    std::vector<double> b_brho, b_len;
    std::vector<double> v1, v2, v3, v4; // 整弧 / 段均值ds20 / 段截断ds20 / 段截断ds40
    int ncl = 0;

    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& cl : clusters) {
            if (cl.points.size() < 12) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetChargeRadius(0.0);
            EstimateResult r;
            if (!est.Estimate(cl, cloud, 2.85, r) || !(r.brho > 0) || !(r.dEdx > 0)) continue;

            std::vector<double> a, b, c;
            Segments(cl.points, 20.0, 0.0, a);   // 段均值
            Segments(cl.points, 20.0, 0.3, b);   // 段截断 30%
            Segments(cl.points, 40.0, 0.3, c);   // 段截断 30%，Δs=40

            b_brho.push_back(r.brho);
            b_len.push_back(r.arclength);
            v1.push_back(r.dEdx);
            v2.push_back(a.empty() ? r.dEdx : a[0]);
            v3.push_back(b.empty() ? r.dEdx : b[0]);
            v4.push_back(c.empty() ? r.dEdx : c[0]);
            ++ncl;
        }
    }
    std::printf("簇数 %d\n\n", ncl);
    Score("V1 整弧 Q/L(现状)", v1, b_brho, b_len);
    Score("V2 段均值 ds=20", v2, b_brho, b_len);
    Score("V3 段截断30% ds=20", v3, b_brho, b_len);
    Score("V4 段截断30% ds=40", v4, b_brho, b_len);
    std::printf("\n（±15% = 有希望分清 p/d 的量级；现状 RMS 越大越糊）\n");
}
