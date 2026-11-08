#ifndef RK4_RECO_CLUSTER_H
#define RK4_RECO_CLUSTER_H

//----------------------------------------------------------------------
// 径迹簇
//----------------------------------------------------------------------

#include "../sim/Simulator.h"

#include <vector>

namespace rk4 {

enum class Direction {
    None = -1,
    Forward = 0,  // 沿 +z 走向 pad 面
    Backward = 1  // 沿 -z 走向靶窗
};

struct Cluster {
    int label = -1;
    Direction direction = Direction::None;
    std::vector<CloudPoint> points;
    int original_label = -1;
};

struct ClusterParams {
    // Explicit DBSCAN remains available; HDBSCAN is the production default.
    double eps_mm = 10.0;
    int min_points = 3;
    int min_cluster_size = 5;
    double direction_fraction = 0.25;
    bool use_hdbscan = true;
    int min_cloud_size = 50;
    int hdbscan_min_samples = 3;
    double min_size_scale_factor = 0.05;
    double hdbscan_z_scale = 0.584;
    double hdbscan_selection_epsilon = 10.0;
    int hdbscan_max_points = 0; // No implicit DBSCAN fallback.
    bool hdbscan_allow_single = false;
    int overlap_min_cluster_size = 15;
    double overlap_ratio = 0.25;
    double outlier_scale_factor = 0.05;

    // ---- 自适应邻域 ----
    // 真实点云的最近邻间距可能比模拟数据大一个量级（实测 9–27 mm vs 2 mm），
    // 固定 eps 会导致碎裂。开启后 eps = max(eps_mm, eps_scale * NN 中位数)。
    bool adaptive_eps = false;
    double eps_scale = 2.0;
    double eps_max_mm = 40.0; // 自适应上调的上限：超过它会把整条径迹连成一团

    // ---- 簇合并（把被切碎的同一径迹接起来）----
    bool enable_join = true;
    double join_max_gap_mm = 60.0;      // 两个端点的最大间距
    double join_max_angle_deg = 65.0;   // 端点切向与连接方向的夹角上限
    int join_edge_points = 10;          // 估计端点切向用的近端点数
};

} // namespace rk4

#endif // RK4_RECO_CLUSTER_H
