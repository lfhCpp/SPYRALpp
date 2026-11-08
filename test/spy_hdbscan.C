// 复现 Spyral 的 HDBSCAN 聚类：z 缩放 0.584、自适应 min_size=max(10,0.05N)、min_samples=3、
// cluster_selection_epsilon=10（用"两簇最近距离<10mm 就合并"近似）
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "io/CloudIO.h"
#include "reco/Hdbscan.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <numeric>
#include <vector>
using namespace rk4;

static double Med(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

static double MinDist(const std::vector<CloudPoint>& a, const std::vector<CloudPoint>& b) {
    double best = 1e300;
    for (const auto& p : a)
        for (const auto& q : b) {
            const double d = std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) +
                                       (p.z - q.z) * (p.z - q.z));
            if (d < best) best = d;
        }
    return best;
}

void spyral_hdbscan() {
    std::vector<PointCloud> clouds;
    if (!ReadPointClouds("/tmp/our_2000.root", clouds)) return;

    size_t ncl = 0;
    std::vector<double> sizes, res;
    for (const auto& cloud : clouds) {
        if (cloud.points.empty()) continue;
        const int npts = (int)cloud.points.size();
        int min_size = (int)(0.05 * npts);
        if (min_size < 10) min_size = 10;

        PointCloud scaled = cloud;
        for (auto& p : scaled.points) p.z *= 0.584;
        HdbscanClusterer h;
        h.Run(scaled, 3, min_size, true);

        // 按 label 重建（用原始点）
        std::map<int, std::vector<CloudPoint>> groups;
        for (size_t i = 0; i < cloud.points.size(); ++i) {
            const int lab = h.point_labels_[i];
            if (lab < 0) continue;
            groups[lab].push_back(cloud.points[i]);
        }
        std::vector<std::vector<CloudPoint>> cl;
        for (auto& kv : groups) cl.push_back(kv.second);

        // cluster_selection_epsilon=10：最近距离 <10mm 的簇合并（并查集）
        std::vector<int> parent(cl.size());
        std::iota(parent.begin(), parent.end(), 0);
        std::function<int(int)> find = [&](int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); };
        for (size_t i = 0; i < cl.size(); ++i)
            for (size_t j = i + 1; j < cl.size(); ++j)
                if (MinDist(cl[i], cl[j]) < 10.0) parent[find((int)i)] = find((int)j);
        std::map<int, std::vector<CloudPoint>> merged;
        for (size_t i = 0; i < cl.size(); ++i) {
            auto& dst = merged[find((int)i)];
            dst.insert(dst.end(), cl[i].begin(), cl[i].end());
        }

        for (auto& kv : merged) {
            ++ncl;
            sizes.push_back((double)kv.second.size());
            if ((int)kv.second.size() < 30) continue;
            Cluster c;
            c.direction = Direction::Forward;
            c.points = kv.second;
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
    std::printf("Spyral式HDBSCAN: 簇=%zu 中位=%.1f 簇/云=%.2f 估计=%zu RMS=%.3f\n", ncl, Med(sizes),
                (double)ncl / std::max<size_t>(1, clouds.size()), res.size(),
                res.empty() ? 0.0 : std::sqrt(s2 / res.size()));
    std::printf("Spyral 目标(2000 事件): 1865 簇，中位 49，簇/云 0.93\n");
}
