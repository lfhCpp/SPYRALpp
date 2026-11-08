// 用"沿路径的 Bethe 积分"当参考（而不是顶点处的 Bethe），看残差宽度和相关性能不能消掉
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

struct Cl {
    double brho = 0, polar = 0, dedx = 0, L = 0;
};

static double Median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// 沿路径积分 Bethe：从顶点 p 出发，每步按 dE/ds = c*Bethe(p) 降能（质子假设）
static double PathAveragedBethe(double brho, double L, double c_scale, double mass) {
    const int n = 40;
    const double ds = L / n;
    double p = 299.792458 * brho; // MeV/c
    if (!(p > 0) || !(L > 0)) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double e = ParticleId::ExpectedDedx(p / 299.792458, 1, 1); // Bethe at current p
        sum += e;
        // 该步损失（把 Bethe 换算成 MeV：用一个常数吸收刻度）
        const double dE = c_scale * e * ds;
        const double E = std::sqrt(p * p + mass * mass) - mass;
        double Enew = E - dE;
        if (Enew < 0.5) Enew = 0.5;
        p = std::sqrt(Enew * (Enew + 2.0 * mass));
    }
    return sum / n;
}

void pid_patheff() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidB_cloud.root", clouds)) return;
    Clusterer clusterer;
    std::vector<Cl> cls;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        for (const auto& c : clusters) {
            if (c.points.size() < 12) continue;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(40.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0) || r.arclength < 20.0) continue;
            cls.push_back(Cl{r.brho, r.polar, r.dEdx, r.arclength});
        }
    }
    const double mass = 938.272;
    std::printf("簇 %zu\n", cls.size());

    // 先定刻度 c：让"顶点 Bethe"参考的中位残差为 0（c 只影响沿程演化强度）
    std::vector<double> ratio;
    for (const auto& cl : cls) {
        const double e = ParticleId::ExpectedDedx(cl.brho, 1, 1);
        if (e > 0) ratio.push_back(cl.dedx / e);
    }
    const double c0 = 1.0 / std::max(1e-9, Median(ratio)); // measured = c0*Bethe (units)

    auto report = [&](const char* tag, bool path) {
        std::vector<double> res, pol;
        int in15 = 0;
        for (const auto& cl : cls) {
            const double ref = path ? PathAveragedBethe(cl.brho, cl.L, 0.02, mass)
                                    : ParticleId::ExpectedDedx(cl.brho, 1, 1);
            if (!(ref > 0)) continue;
            const double rr = (cl.dedx * c0) / ref;
            res.push_back(std::log10(rr));
            pol.push_back(cl.polar);
            if (std::fabs(rr - 1.0) < 0.15) ++in15;
        }
        const double med = Median(res);
        double s2 = 0;
        for (double v : res) s2 += (v - med) * (v - med);
        // 残差-极角相关
        double mp = 0, mr = 0;
        for (size_t i = 0; i < res.size(); ++i) { mp += pol[i]; mr += res[i]; }
        mp /= res.size(); mr /= res.size();
        double sp = 0, sr = 0, spr = 0;
        for (size_t i = 0; i < res.size(); ++i) {
            const double a = pol[i] - mp, b = res[i] - mr;
            sp += a * a; sr += b * b; spr += a * b;
        }
        std::printf("%-28s n=%zu RMS=%.3f dex  r(polar)=%+.2f  ±15%%=%.0f%%\n", tag, res.size(),
                    std::sqrt(s2 / res.size()), (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0,
                    100.0 * in15 / res.size());
    };

    report("参考 = Bethe(顶点)", false);
    report("参考 = 沿路径 Bethe 平均", true);
}
