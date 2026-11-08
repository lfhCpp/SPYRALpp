// 平滑窗口扫描：能不能消掉 r(log dE/dx, log L) = -0.62 并把 RMS 压下去
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

static double Q(std::vector<double> v, double f) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(f * v.size()))];
}
static double Pears(const std::vector<double>& x, const std::vector<double>& y) {
    const size_t n = std::min(x.size(), y.size());
    if (n < 3) return 0;
    double mx = 0, my = 0;
    for (size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double sxy = 0, sxx = 0, syy = 0;
    for (size_t i = 0; i < n; ++i) {
        const double a = x[i] - mx, b = y[i] - my;
        sxy += a * b; sxx += a * a; syy += b * b;
    }
    return (sxx > 0 && syy > 0) ? sxy / std::sqrt(sxx * syy) : 0;
}

// 局部二次平滑
static void Smooth(const std::vector<double>& zs, std::vector<double>& xs, int hw) {
    const int n = (int)zs.size();
    if (hw <= 0 || n < 5) return;
    std::vector<double> out(xs.size(), 0.0);
    for (int i = 0; i < n; ++i) {
        const int lo = std::max(0, i - hw), hi = std::min(n - 1, i + hw);
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, t0 = 0, t1 = 0, t2 = 0;
        for (int j = lo; j <= hi; ++j) {
            const double z = zs[j] - zs[i];
            const double y = xs[j];
            s0 += 1; s1 += z; s2 += z * z; s3 += z * z * z; s4 += z * z * z * z;
            t0 += y; t1 += y * z; t2 += y * z * z;
        }
        const double det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s3 * s2) + s2 * (s1 * s3 - s2 * s2);
        if (std::fabs(det) < 1e-12) { out[i] = xs[i]; continue; }
        const double a0 = ((s2 * s4 - s3 * s3) * t0 + (s2 * s3 - s1 * s4) * t1 + (s1 * s3 - s2 * s2) * t2) / det;
        out[i] = a0;
    }
    xs.swap(out);
}

struct Cl {
    std::vector<double> z, x, y, q;
    double brho = 0, polar = 0;
};

void pid_smoothscan() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/spy_cloud.root", clouds)) return;
    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);

    std::vector<Cl> cls;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> cs;
        clusterer.Run(cloud, cs);
        for (const auto& c : cs) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0)) continue;
            Cl cl;
            cl.brho = r.brho;
            cl.polar = r.polar;
            for (const auto& p : c.points) {
                if (std::hypot(p.x, p.y) > 152.0) continue;
                cl.z.push_back(p.z);
                cl.x.push_back(p.x);
                cl.y.push_back(p.y);
                cl.q.push_back(p.integral);
            }
            if (cl.z.size() < 20) continue;
            cls.push_back(std::move(cl));
        }
    }
    std::printf("簇 %zu\n\n%8s %10s %10s %10s %8s\n", cls.size(), "half_win", "RMS(dex)",
                "r(dEdx,L)", "r(polar)", "±15%");
    for (int hw : {0, 3, 6, 12, 20, 40}) {
        std::vector<double> dedx, L, th, brho;
        for (auto cl : cls) {
            auto xs = cl.x, ys = cl.y;
            Smooth(cl.z, xs, hw);
            Smooth(cl.z, ys, hw);
            double len = 0;
            const int N = 1000;
            const double z0 = cl.z.front(), z1 = cl.z.back();
            auto at = [&](double z, double& x, double& y) {
                size_t lo = 0;
                while (lo + 2 <= cl.z.size() && ((z1 > z0) ? (cl.z[lo + 1] < z) : (cl.z[lo + 1] > z))) ++lo;
                const size_t hi = std::min(lo + 1, cl.z.size() - 1);
                const double dz = cl.z[hi] - cl.z[lo];
                const double t = (std::fabs(dz) > 1e-9) ? (z - cl.z[lo]) / dz : 0;
                x = xs[lo] + t * (xs[hi] - xs[lo]);
                y = ys[lo] + t * (ys[hi] - ys[lo]);
            };
            double xp, yp, zp = z0;
            at(z0, xp, yp);
            for (int s = 1; s <= N; ++s) {
                const double z = z0 + (z1 - z0) * s / N;
                double x, y;
                at(z, x, y);
                len += std::sqrt((x - xp) * (x - xp) + (y - yp) * (y - yp) + (z - zp) * (z - zp));
                xp = x; yp = y; zp = z;
            }
            double qs = 0;
            for (double v : cl.q) qs += v;
            if (len > 1) { dedx.push_back(qs / len); L.push_back(len); th.push_back(cl.polar * 180 / M_PI); brho.push_back(cl.brho); }
        }
        if (dedx.empty()) continue;
        std::vector<double> ratio;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
            if (e > 0) ratio.push_back(dedx[i] / e);
        }
        const double c = 1.0 / Q(ratio, 0.5);
        std::vector<double> res, ld, lL;
        int in15 = 0;
        for (size_t i = 0; i < dedx.size(); ++i) {
            const double e = ParticleId::ExpectedDedx(brho[i], 1, 1);
            if (!(e > 0)) continue;
            const double rr = dedx[i] / (c * e);
            res.push_back(std::log10(rr));
            ld.push_back(std::log(dedx[i]));
            lL.push_back(std::log(L[i]));
            if (std::fabs(rr - 1) < 0.15) ++in15;
        }
        const double med = Q(res, 0.5);
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        std::printf("%8d %10.3f %10.2f %10.2f %7.0f%%\n", hw, std::sqrt(s2 / res.size()),
                    Pears(ld, lL), Pears(res, th), 100.0 * in15 / res.size());
    }
}
