// 验证三条自洽性：
//  1) dE/dx 是否随 Δs 变化（10/20/40 mm，covered 归一）
//  2) dE/dx 是否随覆盖率/点数系统平移（相关系数）
//  3) 归一化到 Bethe_p 后的残差 RMS 哪个变体最小
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

struct Stat { double dedx = 0, geo = 0, cov = 0, npts = 0, brho = 0; };

static double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// 沿 z 累积弧长（折线），按 Δs 分箱，返回 covered 归一与几何归一两种 dE/dx
static void Segment(const std::vector<CloudPoint>& pts, double ds_mm, double& dedx_cov,
                    double& dedx_geo, double& cov, double& lgeo) {
    dedx_cov = dedx_geo = cov = lgeo = 0.0;
    if (pts.size() < 3) return;
    std::vector<double> s(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) {
        const double dx = pts[i].x - pts[i - 1].x;
        const double dy = pts[i].y - pts[i - 1].y;
        const double dz = pts[i].z - pts[i - 1].z;
        s[i] = s[i - 1] + std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    lgeo = s.back();
    if (!(lgeo > 1.0)) return;

    const int nb = std::max(1, (int)std::ceil(lgeo / ds_mm));
    std::vector<double> q(nb, 0.0);
    std::vector<int> n(nb, 0);
    for (size_t i = 0; i < pts.size(); ++i) {
        int b = (int)(s[i] / ds_mm);
        if (b >= nb) b = nb - 1;
        q[b] += pts[i].integral;
        ++n[b];
    }
    double qsum = 0, lcov = 0;
    for (int b = 0; b < nb; ++b) {
        if (n[b] == 0) continue;
        const double len = std::min(ds_mm, lgeo - b * ds_mm);
        qsum += q[b];
        lcov += len;
    }
    cov = lcov / lgeo;
    dedx_cov = (lcov > 1.0) ? qsum / lcov : 0.0;
    dedx_geo = qsum / lgeo;
}

void verify_dedx() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pid_B_cloud.root", clouds)) { std::printf("读不到点云\n"); return; }

    Clusterer clusterer;
    Estimator estimator;
    std::vector<Stat> st10, st20, st40;
    std::vector<double> cov_all, npts_all, brho_all, res20;

    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& cl : clusters) {
            if (cl.points.size() < 10) continue;
            EstimateResult est;
            if (!estimator.Estimate(cl, 2.85, est) || !(est.brho > 0)) continue;

            Stat a, b, c;
            a.npts = b.npts = c.npts = (double)cl.points.size();
            a.brho = b.brho = c.brho = est.brho;
            Segment(cl.points, 10.0, a.dedx, a.geo, a.cov, a.geo);
            Segment(cl.points, 20.0, b.dedx, b.geo, b.cov, b.geo);
            Segment(cl.points, 40.0, c.dedx, c.geo, c.cov, c.geo);
            if (!(b.dedx > 0)) continue;
            st10.push_back(a);
            st20.push_back(b);
            st40.push_back(c);
            cov_all.push_back(b.cov);
            npts_all.push_back((double)cl.points.size());
            brho_all.push_back(est.brho);

            // 残差：covered 归一 / 几何归一 各算一次（对 Bethe_p）
            const double exp = ParticleId::ExpectedDedx(est.brho, 1, 1);
            if (exp > 0) {
                res20.push_back(std::log10(b.dedx / exp));
            }
        }
    }
    if (st20.empty()) { std::printf("没有可用簇\n"); return; }

    auto med = [](const std::vector<Stat>& v, double Stat::*m) {
        std::vector<double> a;
        for (auto& s : v) a.push_back(s.*m);
        return Median(a);
    };
    std::printf("簇数：%zu\n", st20.size());
    std::printf("\n【1】Δs 依赖（covered 归一，中位数）：\n");
    std::printf("   Δs=10mm: %.1f   Δs=20mm: %.1f   Δs=40mm: %.1f  (a.u.)\n",
                med(st10, &Stat::dedx), med(st20, &Stat::dedx), med(st40, &Stat::dedx));
    std::printf("   相对 Δs=20：10mm %+.1f%%  40mm %+.1f%%\n",
                100.0 * (med(st10, &Stat::dedx) / med(st20, &Stat::dedx) - 1.0),
                100.0 * (med(st40, &Stat::dedx) / med(st20, &Stat::dedx) - 1.0));

    std::printf("\n【1b】几何归一 vs covered 归一（Δs=20）：\n");
    std::printf("   几何中位 %.1f   covered 中位 %.1f   覆盖率中位 %.2f\n",
                med(st20, &Stat::geo), med(st20, &Stat::dedx), Median(cov_all));

    // 【2】与覆盖率/点数的相关性（Pearson on logs）
    auto pearson = [](const std::vector<double>& x, const std::vector<double>& y) {
        const size_t n = std::min(x.size(), y.size());
        if (n < 3) return 0.0;
        double mx = 0, my = 0;
        for (size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
        mx /= n; my /= n;
        double sxy = 0, sxx = 0, syy = 0;
        for (size_t i = 0; i < n; ++i) {
            const double dx = x[i] - mx, dy = y[i] - my;
            sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
        }
        return (sxx > 0 && syy > 0) ? sxy / std::sqrt(sxx * syy) : 0.0;
    };
    std::vector<double> ldedx;
    for (auto& s : st20) ldedx.push_back(std::log(s.dedx));
    std::vector<double> lcov, lnpts;
    for (double c : cov_all) lcov.push_back(std::log(std::max(c, 1e-3)));
    for (double c : npts_all) lnpts.push_back(std::log(std::max(c, 1.0)));

    std::printf("\n【2】dE/dx 与长度的相关性（Pearson）：\n");
    std::printf("   r(log dE/dx, log coverage) = %+.2f\n", pearson(ldedx, lcov));
    std::printf("   r(log dE/dx, log npoints ) = %+.2f\n", pearson(ldedx, lnpts));

    // 【3】残差 RMS
    double s1 = 0, s2 = 0;
    for (double r : res20) { s1 += r; s2 += r * r; }
    const double mean = s1 / res20.size();
    const double rms = std::sqrt(std::max(0.0, s2 / res20.size() - mean * mean));
    std::printf("\n【3】Bethe 残差（covered 归一，Δs=20）：中心 %.2f dex，RMS %.3f dex\n",
                mean, rms);
    std::printf("   （对照：现在程序里整弧口径的 RMS = 0.116 dex）\n");
}
