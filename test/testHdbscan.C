//----------------------------------------------------------------------
// HDBSCAN 测试：
//   1) 互达性 MST + 稳定性：两个密度差 7 倍的团
//   2) allow_single_cluster 语义：false = 不允许"整云缩成一簇"，此时该事件无输出
//   3) 默认 HDBSCAN；显式点数限制超限报错
//   4) 根簇选择的噪声点语义与 sklearn 对齐
//   root -l -b -q test/testHdbscan.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/Clusterer.h"
#include "reco/Hdbscan.h"

#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

// 造一个点云：团 A 密集（间距 2 mm），团 B 稀疏（间距 15 mm），彼此远
static PointCloud MakeTwoDensities() {
    PointCloud c;
    c.event = 1;

    // 团 A：沿 z 的一条密线
    for (int i = 0; i < 120; ++i) {
        CloudPoint p;
        p.x = 0.0;
        p.y = 0.0;
        p.z = 100.0 + 2.0 * i;
        p.pad = 1000 + i % 20;
        p.amplitude = 500.0;
        p.integral = 900.0;
        p.tb = static_cast<double>(i);
        p.scale = 1.0;
        c.points.push_back(p);
    }

    // 团 B：一条疏线（远离 A）
    for (int i = 0; i < 120; ++i) {
        CloudPoint p;
        p.x = 200.0;
        p.y = 0.0;
        p.z = 100.0 + 15.0 * i;
        p.pad = 3000 + i % 20;
        p.amplitude = 500.0;
        p.integral = 900.0;
        p.tb = static_cast<double>(i);
        p.scale = 1.0;
        c.points.push_back(p);
    }
    return c;
}

// 一条 200 点的直线径迹
static PointCloud MakeSingleTrack() {
    PointCloud single;
    single.event = 2;
    for (int i = 0; i < 200; ++i) {
        CloudPoint q;
        q.x = 0.0;
        q.y = 0.0;
        q.z = 50.0 + 2.0 * i;
        q.pad = 1000 + i % 20;
        q.amplitude = 500.0;
        q.integral = 900.0;
        q.tb = static_cast<double>(i);
        q.scale = 1.0;
        single.points.push_back(q);
    }
    return single;
}

void testHdbscan() {
    // Expected values independently generated with sklearn 1.7.2 HDBSCAN.
    g_fail = 0;
    Check(ClusterParams{}.use_hdbscan, "Default backend is HDBSCAN");
    const auto cloud = MakeTwoDensities();
    for (bool allow : {false, true}) {
        HdbscanClusterer hdb;
        Check(hdb.Run(cloud, 5, 20, allow) == 2, "Both densities retained as separate clusters");
        Check(hdb.clusters_.size() == 2 && hdb.clusters_[0].points.size() == 120 &&
              hdb.clusters_[1].points.size() == 120, "Each density has all 120 points");
    }
    ClusterParams params;
    params.hdbscan_min_samples = 5;
    params.min_cluster_size = 20;
    params.min_size_scale_factor = 0;
    params.hdbscan_selection_epsilon = 0;
    params.hdbscan_z_scale = 1;
    params.hdbscan_allow_single = true;
    params.enable_join = false;
    Clusterer via; via.SetParams(params);
    std::vector<Cluster> clusters;
    Check(via.Run(cloud, clusters) == 2, "Dispatcher preserves both clusters after LOF");
    Check(via.LastBackendHdbscan(), "Dispatcher used HDBSCAN");
    params.hdbscan_max_points = 50;
    via.SetParams(params);
    bool threw = false;
    try { via.Run(cloud, clusters); } catch (const std::runtime_error&) { threw = true; }
    Check(threw, "Explicit point guard throws instead of silently selecting DBSCAN");
    HdbscanClusterer single;
    const auto line = MakeSingleTrack();
    Check(single.Run(line, 5, 20, false) == 0, "Root-only uniform line rejected with allow_single=false");
    Check(single.Run(line, 5, 20, true) == 1, "Root-only uniform line retained with allow_single=true");
    Check(single.clusters_.size() == 1 && single.clusters_[0].points.size() == 196,
          "Python root selection retains 196 points and rejects 4 border points");
    std::printf("HDBSCAN regression failures: %d\n", g_fail);
    gSystem->Exit(g_fail ? 1 : 0);
}
