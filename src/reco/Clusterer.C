#include "Clusterer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <utility>

#include "../core/Units.h"

#include "Hdbscan.h"

namespace rk4 {
namespace {

long long CellKey(int ix, int iy, int iz) {
    // 每个维度偏移到非负，再打包进 64 位
    const long long X = static_cast<long long>(ix) + (1 << 20);
    const long long Y = static_cast<long long>(iy) + (1 << 20);
    const long long Z = static_cast<long long>(iz) + (1 << 20);
    return (X << 42) | (Y << 21) | Z;
}

// 采样估计最近邻距离的中位数（用于自适应 eps）
double MedianNearestNeighbor(const PointCloud& cloud, int max_samples) {
    const size_t n = cloud.points.size();
    if (n < 3) return 0.0;

    const size_t step = std::max<size_t>(1, n / static_cast<size_t>(max_samples));
    std::vector<double> nn;
    nn.reserve(n / step + 1);

    for (size_t i = 0; i < n; i += step) {
        double best = 1.0e18;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            const double dx = cloud.points[i].x - cloud.points[j].x;
            const double dy = cloud.points[i].y - cloud.points[j].y;
            const double dz = cloud.points[i].z - cloud.points[j].z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best) best = d2;
        }
        nn.push_back(std::sqrt(best));
    }

    if (nn.empty()) return 0.0;
    std::sort(nn.begin(), nn.end());
    return nn[nn.size() / 2];
}

} // namespace

int Clusterer::Run(const PointCloud& cloud, std::vector<Cluster>& clusters) {
    clusters.clear();
    last_hdbscan_ = false;
    const int n = static_cast<int>(cloud.points.size());
    if (n < params_.min_points) return 0;

    if (params_.use_hdbscan) {
        last_hdbscan_=true;last_eps_=0;last_median_nn_=0;last_join_count_=0;
        if(n<params_.min_cloud_size)return 0;
        if(params_.hdbscan_max_points>0 && n>params_.hdbscan_max_points)
            throw std::runtime_error("Configured HDBSCAN point limit exceeded; refusing implicit DBSCAN fallback");
        const int min_size=std::max(params_.min_cluster_size,static_cast<int>(params_.min_size_scale_factor*n));
        HdbscanClusterer hdb;
        hdb.Run(cloud,params_.hdbscan_min_samples,min_size,params_.hdbscan_allow_single,
                params_.hdbscan_selection_epsilon,params_.hdbscan_z_scale);
        return FinalizeSpyral(hdb.clusters_,clusters);
    }

    double eps = params_.eps_mm;
    last_median_nn_ = MedianNearestNeighbor(cloud, 300);
    if (params_.adaptive_eps) {
        double adaptive = params_.eps_scale * last_median_nn_;
        if (adaptive > params_.eps_max_mm) adaptive = params_.eps_max_mm;
        if (adaptive > eps) eps = adaptive;
    }
    last_eps_ = eps;
    const double eps2 = eps * eps;

    // 空间网格：每个点进一个 cell
    std::unordered_map<long long, std::vector<int>> grid;
    grid.reserve(static_cast<size_t>(n) * 2);
    for (int i = 0; i < n; ++i) {
        const CloudPoint& p = cloud.points[i];
        const int ix = static_cast<int>(std::floor(p.x / eps));
        const int iy = static_cast<int>(std::floor(p.y / eps));
        const int iz = static_cast<int>(std::floor(p.z / eps));
        grid[CellKey(ix, iy, iz)].push_back(i);
    }

    // 邻域查询：3x3x3 个 cell 内距离 <= eps 的点
    auto region_query = [&](int idx) {
        std::vector<int> result;
        const CloudPoint& p = cloud.points[idx];
        const int ix = static_cast<int>(std::floor(p.x / eps));
        const int iy = static_cast<int>(std::floor(p.y / eps));
        const int iz = static_cast<int>(std::floor(p.z / eps));
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    const auto it = grid.find(CellKey(ix + dx, iy + dy, iz + dz));
                    if (it == grid.end()) continue;
                    for (int j : it->second) {
                        if (j == idx) continue;
                        const CloudPoint& q = cloud.points[j];
                        const double d2 = (p.x - q.x) * (p.x - q.x) +
                                          (p.y - q.y) * (p.y - q.y) +
                                          (p.z - q.z) * (p.z - q.z);
                        if (d2 <= eps2) result.push_back(j);
                    }
                }
            }
        }
        return result;
    };

    // -1 = 未访问，-2 = 噪声，>=0 = 簇编号
    std::vector<int> label(static_cast<size_t>(n), -1);
    int cluster_id = 0;

    for (int i = 0; i < n; ++i) {
        if (label[i] != -1) continue;

        std::vector<int> neighbours = region_query(i);
        if (static_cast<int>(neighbours.size()) + 1 < params_.min_points) {
            label[i] = -2;
            continue;
        }

        label[i] = cluster_id;
        std::vector<int> seeds = neighbours;

        for (size_t k = 0; k < seeds.size(); ++k) {
            const int j = seeds[k];
            if (label[j] == -2) label[j] = cluster_id; // 边界点
            if (label[j] != -1) continue;

            label[j] = cluster_id;
            std::vector<int> more = region_query(j);
            if (static_cast<int>(more.size()) + 1 >= params_.min_points) {
                seeds.insert(seeds.end(), more.begin(), more.end());
            }
        }

        ++cluster_id;
    }

    // 组装并过滤小簇
    std::vector<Cluster> built(static_cast<size_t>(cluster_id));
    for (int i = 0; i < n; ++i) {
        if (label[i] < 0) continue;
        built[static_cast<size_t>(label[i])].points.push_back(cloud.points[i]);
    }

    return Finalize(built, clusters);
}

int Clusterer::RunHdbscanSnapshot(const PointCloud& cloud,
                                  std::vector<Cluster>& pre_finalize,
                                  std::vector<int>& point_labels,
                                  ClusterSnapshotTiming* timing) {
    using Clock = std::chrono::steady_clock;
    const auto started = Clock::now();
    if (timing) *timing = {};
    pre_finalize.clear();
    point_labels.assign(cloud.points.size(), -1);
    last_hdbscan_ = false;
    last_join_count_ = 0;
    const int n = static_cast<int>(cloud.points.size());
    if (!params_.use_hdbscan)
        throw std::invalid_argument("RunHdbscanSnapshot requires HDBSCAN to be enabled");
    last_eps_ = 0.0;
    last_median_nn_ = 0.0;
    last_eps_ = 0.0;
    last_median_nn_ = 0.0;
    if (n < params_.min_points || n < params_.min_cloud_size) return 0;
    if (params_.hdbscan_max_points > 0 && n > params_.hdbscan_max_points)
        throw std::runtime_error("Configured HDBSCAN point limit exceeded; refusing implicit DBSCAN fallback");

    const int min_size = std::max(params_.min_cluster_size,
        static_cast<int>(params_.min_size_scale_factor * n));
    HdbscanClusterer hdb;
    hdb.Run(cloud, params_.hdbscan_min_samples, min_size,
            params_.hdbscan_allow_single, params_.hdbscan_selection_epsilon,
            params_.hdbscan_z_scale);
    pre_finalize = std::move(hdb.clusters_);
    point_labels = std::move(hdb.point_labels_);
    last_hdbscan_ = true;
    if (timing) timing->hdbscan_seconds = std::chrono::duration<double>(Clock::now() - started).count();
    return static_cast<int>(pre_finalize.size());
}

int Clusterer::FinalizeHdbscanSnapshot(const std::vector<Cluster>& pre_finalize,
                                       bool apply_original_merge,
                                       std::vector<Cluster>& finalized,
                                       ClusterSnapshotTiming* timing,
                                       const ExperimentalClusterMerge& experimental_merge,
                                       std::vector<OriginalMergeDecision>* original_decisions,
                                       bool apply_lof,
                                       const ExperimentalOriginalMergeVeto& original_merge_veto) const {
    if (timing) *timing = {};
    Clusterer replay = *this;
    replay.params_.enable_join = apply_original_merge;
    std::vector<Cluster> working = pre_finalize;
    finalized.clear();
    if(original_decisions)original_decisions->clear();
    return replay.FinalizeSpyral(working, finalized, timing, experimental_merge, original_decisions,
                                 apply_lof, original_merge_veto);
}

//----------------------------------------------------------------------
// 两条路径（DBSCAN / HDBSCAN）共用的收尾
//
// 小簇过滤 → 按 z 排序 → 首尾电荷判方向 → 连续编号 → continuity join。
// 放在这里是为了保证换聚类后端时下游看到的簇结构完全一致。
//----------------------------------------------------------------------
int Clusterer::Finalize(std::vector<Cluster>& built, std::vector<Cluster>& clusters) {
    for (auto& c : built) {
        if (static_cast<int>(c.points.size()) < params_.min_cluster_size) continue;

        // 按 z 排序，便于估计和方向判断
        std::sort(c.points.begin(), c.points.end(),
                  [](const CloudPoint& a, const CloudPoint& b) { return a.z < b.z; });

        // 方向：首尾各取一段，比较积分电荷
        const size_t m = c.points.size();
        size_t edge = static_cast<size_t>(std::max(1.0, m * params_.direction_fraction));
        if (edge > m / 2) edge = m / 2;

        double low_charge = 0.0, high_charge = 0.0;
        for (size_t i = 0; i < edge; ++i) low_charge += c.points[i].integral;
        for (size_t i = m - edge; i < m; ++i) high_charge += c.points[i].integral;

        c.direction = (high_charge >= low_charge) ? Direction::Forward : Direction::Backward;
        c.label = static_cast<int>(clusters.size());
        clusters.push_back(c);
    }

    if (params_.enable_join) {
        JoinClusters(clusters);
    }

    return static_cast<int>(clusters.size());
}

//----------------------------------------------------------------------
// 簇合并：把被 DBSCAN 切碎的同一径迹接起来
//
// 判据（参考 Spyral 的 continuity join）：
//   1. 两簇最近端点间距 < join_max_gap_mm
//   2. 端点切向与连接方向的夹角都小于 join_max_angle_deg
//      —— 即两条径迹"首尾相接"而不是"并排靠近"
//----------------------------------------------------------------------
namespace {

struct EdgeInfo {
    bool valid = false;
    int index = -1;
    double x = 0.0, y = 0.0, z = 0.0;
    double in_x = 0.0, in_y = 0.0, in_z = 0.0; // 单位向量，指向簇内部
};

// 找离目标点最近的端点，并估计该端点处的向内切向。
//
// 注意：簇内点按 z 排序，用"相邻 K 个点"估计切向对大 polar 径迹不可靠
// （z 变化小、横向变化大）。这里改用"端点 -> 簇质心"的方向，
// 对细长径迹等价于沿径迹指向内部，且不受排序方向影响。
EdgeInfo FindEdge(const Cluster& cluster, double tx, double ty, double tz,
                  int edge_points) {
    EdgeInfo edge;
    (void)edge_points;
    const size_t n = cluster.points.size();
    if (n < 2) return edge;

    double best = 1.0e18;
    size_t best_i = 0;
    for (size_t i = 0; i < n; ++i) {
        const double dx = cluster.points[i].x - tx;
        const double dy = cluster.points[i].y - ty;
        const double dz = cluster.points[i].z - tz;
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best) {
            best = d2;
            best_i = i;
        }
    }

    // 质心
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (const auto& p : cluster.points) {
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    cx /= static_cast<double>(n);
    cy /= static_cast<double>(n);
    cz /= static_cast<double>(n);

    const CloudPoint& p = cluster.points[best_i];
    double dx = cx - p.x;
    double dy = cy - p.y;
    double dz = cz - p.z;
    const double norm = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (norm < 1.0e-9) return edge;

    edge.valid = true;
    edge.index = static_cast<int>(best_i);
    edge.x = p.x;
    edge.y = p.y;
    edge.z = p.z;
    edge.in_x = dx / norm;
    edge.in_y = dy / norm;
    edge.in_z = dz / norm;
    return edge;
}

} // namespace

int Clusterer::JoinClusters(std::vector<Cluster>& clusters) const {
    last_join_count_ = 0;
    const size_t n = clusters.size();
    if (n < 2) return static_cast<int>(n);

    std::vector<int> parent(n);
    for (size_t i = 0; i < n; ++i) parent[i] = static_cast<int>(i);

    std::function<int(int)> find_root = [&](int i) {
        while (parent[static_cast<size_t>(i)] != i) {
            parent[static_cast<size_t>(i)] =
                parent[static_cast<size_t>(parent[static_cast<size_t>(i)])];
            i = parent[static_cast<size_t>(i)];
        }
        return i;
    };

    const double max_gap2 = params_.join_max_gap_mm * params_.join_max_gap_mm;
    const double cos_limit = std::cos(params_.join_max_angle_deg * Units::DEG2RAD);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (find_root(static_cast<int>(i)) == find_root(static_cast<int>(j))) continue;

            // 注意：不要用 direction 过滤。碎段的首尾电荷分布不可靠，
            // 同一径迹的两段可能被判成不同方向，导致永远无法合并。
            // 交叉径迹的误连风险由 gap 与角度判据共同限制。

            // i 的端点朝 j 找，j 的端点朝 i 找
            const CloudPoint& ci = clusters[i].points[clusters[i].points.size() / 2];
            const CloudPoint& cj = clusters[j].points[clusters[j].points.size() / 2];

            const EdgeInfo ei = FindEdge(clusters[i], cj.x, cj.y, cj.z, params_.join_edge_points);
            const EdgeInfo ej = FindEdge(clusters[j], ci.x, ci.y, ci.z, params_.join_edge_points);
            if (!ei.valid || !ej.valid) continue;

            double ux = ej.x - ei.x;
            double uy = ej.y - ei.y;
            double uz = ej.z - ei.z;
            const double gap2 = ux * ux + uy * uy + uz * uz;
            if (gap2 > max_gap2 || gap2 < 1.0e-12) continue;

            const double gap = std::sqrt(gap2);
            ux /= gap;
            uy /= gap;
            uz /= gap;

            // 离开 i 的方向应约为 -ei.in；进入 j 的方向应约为 ej.in
            const double leave_i = -(ei.in_x * ux + ei.in_y * uy + ei.in_z * uz);
            const double enter_j = ej.in_x * ux + ej.in_y * uy + ej.in_z * uz;

            if (leave_i < cos_limit || enter_j < cos_limit) continue;

            const int ri = find_root(static_cast<int>(i));
            const int rj = find_root(static_cast<int>(j));
            parent[static_cast<size_t>(rj)] = ri;
            ++last_join_count_;
        }
    }

    if (last_join_count_ == 0) return static_cast<int>(n);

    // 重建：按根节点分组
    std::vector<Cluster> merged;
    std::vector<int> root_to_index(n, -1);
    for (size_t i = 0; i < n; ++i) {
        const int root = find_root(static_cast<int>(i));
        if (root_to_index[static_cast<size_t>(root)] < 0) {
            root_to_index[static_cast<size_t>(root)] = static_cast<int>(merged.size());
            merged.push_back(Cluster{});
        }
        const int target = root_to_index[static_cast<size_t>(root)];
        Cluster& dst = merged[static_cast<size_t>(target)];
        dst.points.insert(dst.points.end(), clusters[i].points.begin(), clusters[i].points.end());
    }

    for (auto& c : merged) {
        std::sort(c.points.begin(), c.points.end(),
                  [](const CloudPoint& a, const CloudPoint& b) { return a.z < b.z; });
        // 方向按合并后的首尾电荷重算
        const size_t m = c.points.size();
        size_t edge = static_cast<size_t>(std::max(1.0, m * params_.direction_fraction));
        if (edge > m / 2) edge = m / 2;
        double low_charge = 0.0, high_charge = 0.0;
        for (size_t i = 0; i < edge; ++i) low_charge += c.points[i].integral;
        for (size_t i = m - edge; i < m; ++i) high_charge += c.points[i].integral;
        c.direction = (high_charge >= low_charge) ? Direction::Forward : Direction::Backward;
    }

    for (size_t i = 0; i < merged.size(); ++i) merged[i].label = static_cast<int>(i);
    clusters.swap(merged);
    return static_cast<int>(clusters.size());
}

} // namespace rk4
