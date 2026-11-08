// 责任拆分：残差的 θ 相关性来自 dE/dx 还是 Bρ/Bethe？
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
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(f * v.size()))];
}
static double Pearson(const std::vector<double>& x, const std::vector<double>& y) {
    const size_t n = std::min(x.size(), y.size());
    if (n < 3) return 0.0;
    double mx = 0, my = 0;
    for (size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double sxy = 0, sxx = 0, syy = 0;
    for (size_t i = 0; i < n; ++i) {
        const double a = x[i] - mx, b = y[i] - my;
        sxy += a * b; sxx += a * a; syy += b * b;
    }
    return (sxx > 0 && syy > 0) ? sxy / std::sqrt(sxx * syy) : 0.0;
}

void pid_split() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/spy_cloud.root", clouds)) return;
    Clusterer clusterer;
    ClusterParams cp;
    cp.eps_mm = 20.0;
    cp.min_points = 3;
    cp.min_cluster_size = 50;
    clusterer.SetParams(cp);

    std::vector<double> ldedx, lbethe, lbrho, th, res, L, rad;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& c : clusters) {
            if (c.points.size() < 30) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(25.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            const double e = ParticleId::ExpectedDedx(r.brho, 1, 1);
            if (!(e > 0)) continue;
            ldedx.push_back(std::log(r.dEdx));
            lbethe.push_back(std::log(e));
            lbrho.push_back(std::log(r.brho));
            th.push_back(r.polar * 180.0 / 3.14159265358979323846);
            res.push_back(std::log10(r.dEdx / e));
            L.push_back(r.arclength);
            rad.push_back(r.radius);
        }
    }
    if (res.empty()) { std::printf("无簇\n"); return; }
    std::printf("簇 n=%zu\n\n", res.size());
    std::printf("r(log dE/dx, θ)      = %+.2f\n", Pearson(ldedx, th));
    std::printf("r(log Bethe(Bρ), θ)  = %+.2f\n", Pearson(lbethe, th));
    std::printf("r(log Bρ, θ)         = %+.2f\n", Pearson(lbrho, th));
    std::printf("r(log dE/dx, log L)  = %+.2f\n", Pearson(ldedx, [&] {
        std::vector<double> v;
        for (double x : L) v.push_back(std::log(x));
        return v;
    }()));
    std::printf("r(log R, log L)      = %+.2f\n", Pearson([&] {
        std::vector<double> v;
        for (double x : rad) v.push_back(std::log(x));
        return v;
    }(), [&] {
        std::vector<double> v;
        for (double x : L) v.push_back(std::log(x));
        return v;
    }()));

    // 残差的 θ 依赖：线性拟合 res = a + b·θ，去掉后的 RMS
    std::vector<double> rr = res;
    const double med = Q(rr, 0.5);
    double s2 = 0;
    for (double v : res) s2 += (v - med) * (v - med);
    const double rms0 = std::sqrt(s2 / res.size());

    double mt = 0, mr = 0;
    for (size_t i = 0; i < res.size(); ++i) { mt += th[i]; mr += res[i]; }
    mt /= res.size(); mr /= res.size();
    double stt = 0, str = 0;
    for (size_t i = 0; i < res.size(); ++i) {
        const double a = th[i] - mt, b = res[i] - mr;
        stt += a * a; str += a * b;
    }
    const double slope = (stt > 0) ? str / stt : 0.0;
    double s2c = 0;
    for (size_t i = 0; i < res.size(); ++i) {
        const double v = res[i] - slope * (th[i] - mt);
        s2c += (v - 0) * (v - 0);
    }
    std::printf("\n残差 RMS（原始）        = %.3f dex\n", rms0);
    std::printf("残差 RMS（去掉 θ 斜率后）= %.3f dex   （斜率 %.4f dex/deg）\n",
                std::sqrt(s2c / res.size()), slope);
}
