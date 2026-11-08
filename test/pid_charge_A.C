// 簇纯度对比：DBSCAN vs HDBSCAN vs 单径迹子样本（都在 estimate 之前/之内）
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

struct Acc {
    std::vector<double> dedx, brho, pol;
    int ntot = 0;
};

static void Score(const char* tag, const Acc& a) {
    if (a.dedx.empty()) { std::printf("%-22s 无\n", tag); return; }
    std::vector<double> ratio;
    for (size_t i = 0; i < a.dedx.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(a.brho[i], 1, 1);
        if (e > 0) ratio.push_back(a.dedx[i] / e);
    }
    const double c = 1.0 / std::max(1e-30, Med(ratio));
    std::vector<double> res, p;
    int in15 = 0, in25 = 0;
    for (size_t i = 0; i < a.dedx.size(); ++i) {
        const double e = ParticleId::ExpectedDedx(a.brho[i], 1, 1);
        if (!(e > 0)) continue;
        const double rr = a.dedx[i] / (c * e);
        res.push_back(std::log10(rr));
        p.push_back(a.pol[i]);
        if (std::fabs(rr - 1) < 0.15) ++in15;
        if (std::fabs(rr - 1) < 0.25) ++in25;
    }
    const double med = Med(res);
    double s2 = 0;
    for (double v : res) s2 += (v - med) * (v - med);
    double mp = 0, mr = 0;
    for (size_t i = 0; i < res.size(); ++i) { mp += p[i]; mr += res[i]; }
    mp /= res.size(); mr /= res.size();
    double sp = 0, sr = 0, spr = 0;
    for (size_t i = 0; i < res.size(); ++i) {
        const double x = p[i] - mp, y = res[i] - mr;
        sp += x * x; sr += y * y; spr += x * y;
    }
    std::printf("%-22s 通过 %5zu/%5d (%4.0f%%)  RMS=%.3f dex  r(pol)=%+.2f  ±15%%=%.0f%%  ±25%%=%.0f%%\n",
                tag, res.size(), a.ntot, 100.0 * res.size() / std::max(1, a.ntot),
                std::sqrt(s2 / res.size()),
                (sp > 0 && sr > 0) ? spr / std::sqrt(sp * sr) : 0.0,
                100.0 * in15 / res.size(), 100.0 * in25 / res.size());
}

static void Run(const std::vector<PointCloud>& clouds, const ClusterParams& params,
                const char* tag, bool single_track_only) {
    Clusterer clusterer;
    clusterer.SetParams(params);
    Acc acc;
    for (const auto& cloud : clouds) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        if (clusters.empty()) continue;
        // 单径迹子样本：每事件只取点数最多的簇，且要求 ≥30 点
        size_t i0 = 0;
        if (single_track_only) {
            size_t best = 0;
            for (size_t i = 1; i < clusters.size(); ++i)
                if (clusters[i].points.size() > clusters[best].points.size()) best = i;
            if (clusters[best].points.size() < 30) continue;
            i0 = best;
        }
        for (size_t i = (single_track_only ? i0 : 0);
             i < (single_track_only ? i0 + 1 : clusters.size()); ++i) {
            const Cluster& c = clusters[i];
            if (c.points.size() < 12) continue;
            ++acc.ntot;
            Estimator est;
            est.SetSmallPadRadius(152.0);
            est.SetBeamRegionRadius(40.0);
            EstimateResult r;
            if (!est.Estimate(c, cloud, 2.85, r)) continue;
            if (!(r.brho > 0) || !(r.dEdx > 0)) continue;
            acc.dedx.push_back(r.dEdx);
            acc.brho.push_back(r.brho);
            acc.pol.push_back(r.polar);
        }
    }
    Score(tag, acc);
}

void pid_charge_A() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/pidA_cloud.root", clouds)) return;
    std::printf("点云 %zu\n\n", clouds.size());

    ClusterParams db;
    db.eps_mm = 20.0;
    db.min_points = 3;
    db.min_cluster_size = 10;
    Run(clouds, db, "DBSCAN eps=20", false);
    Run(clouds, db, "DBSCAN 单径迹子样本", true);

}
