#ifndef RK4_RECO_CLUSTERER_H
#define RK4_RECO_CLUSTERER_H

//----------------------------------------------------------------------
// Default: sklearn-compatible HDBSCAN, overlap joining, LOF cleanup.
// DBSCAN/continuity joining is available only when explicitly selected.
// HDBSCAN does not silently switch clustering algorithms.
//----------------------------------------------------------------------

#include "Cluster.h"

#include <chrono>
#include <functional>
#include <string>
#include <vector>

namespace rk4 {

struct ClusterSnapshotTiming {
    double hdbscan_seconds = 0.0;
    double original_merge_seconds = 0.0;
    double experimental_merge_seconds = 0.0;
    double postprocess_seconds = 0.0; // sorting, LOF and output construction
    double lof_seconds = 0.0; // LOF only; included in postprocess_seconds
};

struct OriginalMergeDecision {
    int iteration = 0;
    size_t left_index = 0;
    size_t right_index = 0;
    std::vector<int> left_raw_cluster_ids;
    std::vector<int> right_raw_cluster_ids;
    int left_points = 0;
    int right_points = 0;
    double left_circle_x = 0.0, left_circle_y = 0.0, left_circle_r = 0.0;
    double right_circle_x = 0.0, right_circle_y = 0.0, right_circle_r = 0.0;
    double overlap_fraction = 0.0;
    double threshold = 0.0;
    bool selected_for_merge = false;
    bool vetoed = false;
    std::string reason;
};

using ExperimentalClusterMerge = std::function<void(std::vector<Cluster>&)>;
using ExperimentalOriginalMergeVeto = std::function<bool(const OriginalMergeDecision&)>;

class Clusterer {
public:
    void SetParams(const ClusterParams& params) { params_ = params; }
    const ClusterParams& Params() const { return params_; }

    // 聚类；返回簇的个数。噪声点被丢弃。
    int Run(const PointCloud& cloud, std::vector<Cluster>& clusters);

    // Experimental audit interface: expose the exact HDBSCAN partition before
    // overlap merging/LOF, along with one label per input point (-1 = noise).
    // The production Run() contract and defaults are unchanged.
    int RunHdbscanSnapshot(const PointCloud& cloud,
                           std::vector<Cluster>& pre_finalize,
                           std::vector<int>& point_labels,
                           ClusterSnapshotTiming* timing = nullptr);
    // Replay the normal Spyral finalization on an immutable HDBSCAN snapshot.
    // `apply_original_merge=false` disables only OverlapStep; LOF and all other
    // finalization remain enabled. This enables paired A/D and A/B experiments.
    int FinalizeHdbscanSnapshot(const std::vector<Cluster>& pre_finalize,
                                bool apply_original_merge,
                                std::vector<Cluster>& finalized,
                                ClusterSnapshotTiming* timing = nullptr,
                                const ExperimentalClusterMerge& experimental_merge = {},
                                std::vector<OriginalMergeDecision>* original_decisions = nullptr,
                                bool apply_lof = true,
                                const ExperimentalOriginalMergeVeto& original_merge_veto = {}) const;

    // 合并端点相接、方向连续的簇；返回合并后的簇数
    int JoinClusters(std::vector<Cluster>& clusters) const;

    int LastJoinCount() const { return last_join_count_; }
    // 上次 Run 实际用的是不是 HDBSCAN（点数超限退回 DBSCAN 时为 false）
    bool LastBackendHdbscan() const { return last_hdbscan_; }
    // 本次 Run 实际使用的邻域半径（自适应开启时可能大于 eps_mm）
    double LastEps() const { return last_eps_; }
    // 本次 Run 估计的最近邻中位数
    double LastMedianNN() const { return last_median_nn_; }

private:
    // 两条路径共用的收尾：小簇过滤、按 z 排序、方向、编号、join
    int FinalizeSpyral(std::vector<Cluster>& built,std::vector<Cluster>& clusters,
                       ClusterSnapshotTiming* timing = nullptr,
                       const ExperimentalClusterMerge& experimental_merge = {},
                       std::vector<OriginalMergeDecision>* original_decisions = nullptr,
                       bool apply_lof = true,
                       const ExperimentalOriginalMergeVeto& original_merge_veto = {});
    int Finalize(std::vector<Cluster>& built, std::vector<Cluster>& clusters);

    ClusterParams params_;
    mutable bool last_hdbscan_ = false;
    mutable int last_join_count_ = 0;
    mutable double last_eps_ = 0.0;
    mutable double last_median_nn_ = 0.0;
};

} // namespace rk4

#endif // RK4_RECO_CLUSTERER_H
