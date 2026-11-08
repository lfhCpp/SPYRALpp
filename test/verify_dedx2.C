// 验证新 dE/dx 口径：pad 足迹并集 + 电荷收集
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
static double Pearson(const std::vector<double>& x, const std::vector<double>& y) {
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
}

struct Row {
    double dEdx = 0, arc = 0, geo = 0, npts = 0, brho = 0;
};

static std::vector<Row> Collect(const std::vector<PointCloud>& clouds, double charge_radius,
                                double pitch) {
    std::vector<Row> rows;
    Clusterer clusterer;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& cl : clusters) {
            if (cl.points.size() < 10) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetPadPitch(pitch);
            est.SetChargeRadius(charge_radius);
            EstimateResult r;
            if (!est.Estimate(cl, cloud, 2.85, r) || !(r.brho > 0) || !(r.dEdx > 0)) continue;
            Row row;
            row.dEdx = r.dEdx;
            row.arc = r.arclength;
            row.npts = (double)cl.points.size();
            row.brho = r.brho;
            // 几何总长（折线）
            double len = 0;
            for (size_t i = 1; i < cl.points.size(); ++i) {
                const double dx = cl.points[i].x - cl.points[i - 1].x;
                const double dy = cl.points[i].y - cl.points[i - 1].y;
                const double dz = cl.points[i].z - cl.points[i - 1].z;
                len += std::sqrt(dx * dx + dy * dy + dz * dz);
            }
            row.geo = len;
            rows.push_back(row);
        }
    }
    return rows;
}

static void Report(const char* tag, const std::vector<Row>& rows) {
    std::vector<double> d, arc, geo, npts, res, lcov;
    for (const auto& r : rows) {
        d.push_back(r.dEdx);
        arc.push_back(r.arc);
        geo.push_back(r.geo);
        npts.push_back(r.npts);
        lcov.push_back(std::log(std::max(r.arc / std::max(r.geo, 1.0), 1e-3)));
        const double exp = ParticleId::ExpectedDedx(r.brho, 1, 1);
        if (exp > 0) res.push_back(std::log10(r.dEdx / exp));
    }
    std::vector<double> ld;
    for (double v : d) ld.push_back(std::log(v));
    std::vector<double> ln;
    for (double v : npts) ln.push_back(std::log(v));

    double m = 0, s2 = 0;
    for (double r : res) { m += r; s2 += r * r; }
    m = res.empty() ? 0 : m / res.size();
    const double rms = res.empty() ? 0 : std::sqrt(std::max(0.0, s2 / res.size() - m * m));

    std::printf("%-18s 簇=%zu  dE/dx中位=%.1f  覆盖=%.2f  r(cov)=%+.2f  r(npts)=%+.2f  "
                "残差RMS=%.3f dex\n",
                tag, rows.size(), Median(d), Median([&] {
                    std::vector<double> v;
                    for (size_t i = 0; i < arc.size(); ++i)
                        v.push_back(arc[i] / std::max(geo[i], 1.0));
                    return v;
                }()), Pearson(ld, lcov), Pearson(ld, ln), rms);
}

void verify_dedx2() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pid_B_cloud.root", clouds)) { std::printf("读不到\n"); return; }
    std::printf("=== 旧口径（整弧几何长度，charge_radius=0） ===\n");
    Report("旧: R=0", Collect(clouds, 0.0, 20.0));
    std::printf("\n=== 新口径：足迹并集 + 电荷收集 ===\n");
    Report("新: R=15 p=10", Collect(clouds, 15.0, 10.0));
    Report("新: R=25 p=10", Collect(clouds, 25.0, 10.0));
    std::printf("\n=== pad 间距敏感性（新口径 R=15）===\n");
    Report("p=6", Collect(clouds, 15.0, 6.0));
    Report("p=10", Collect(clouds, 15.0, 10.0));
    Report("p=16", Collect(clouds, 15.0, 16.0));
}
