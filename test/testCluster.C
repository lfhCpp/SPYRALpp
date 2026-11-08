//----------------------------------------------------------------------
// Default HDBSCAN + spline estimator against simulated truth;
// explicit legacy DBSCAN for the continuity-join-specific fixture.
//   root -l -b -q test/testCluster.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/Geometry.h"
#include "sim/DetectorConfig.h"
#include "sim/Simulator.h"
#include "sim/TargetDensity.h"

#include "TSystem.h"

#include <algorithm>
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

//----------------------------------------------------------------------
// 1) 圆拟合
//----------------------------------------------------------------------
static void TestFitCircle() {
    std::printf("[1] 圆拟合\n");

    const double cx_true = 10.0, cy_true = -20.0, r_true = 50.0;
    std::vector<double> xs, ys;
    for (int i = 0; i < 120; ++i) {
        const double t = 2.0 * Units::PI * i / 120.0;
        xs.push_back(cx_true + r_true * std::cos(t));
        ys.push_back(cy_true + r_true * std::sin(t));
    }

    double cx = 0.0, cy = 0.0, r = 0.0;
    Check(FitCircle(xs, ys, cx, cy, r), "FitCircle 返回 true");
    Check(std::fabs(cx - cx_true) < 1e-6, "圆心 x");
    Check(std::fabs(cy - cy_true) < 1e-6, "圆心 y");
    Check(std::fabs(r - r_true) < 1e-6, "半径");
    std::printf("    拟合圆：中心 (%.4f, %.4f) 半径 %.4f\n", cx, cy, r);
}

//----------------------------------------------------------------------
// 构造粒子与事件
//----------------------------------------------------------------------
static Particle MakeParticle(int A, int Z, double ke_mev, double polar_deg,
                             double phi_deg) {
    const double mass = MassData::GetNuclearMassNuclear(A, Z) * Units::u_MeV;
    const double p = std::sqrt(ke_mev * (ke_mev + 2.0 * mass));
    const double th = polar_deg * Units::DEG2RAD;
    const double ph = phi_deg * Units::DEG2RAD;

    Particle part;
    part.A = A;
    part.Z = Z;
    part.label = 2;
    part.px = p * std::sin(th) * std::cos(ph);
    part.py = p * std::sin(th) * std::sin(ph);
    part.pz = p * std::cos(th);
    part.E = ke_mev + mass;
    return part;
}

static Event MakeEvent(std::vector<Particle> parts, double vz_m) {
    Event ev;
    ev.id = 1;
    ev.vx = 0.0;
    ev.vy = 0.0;
    ev.vz = vz_m;
    ev.particles = std::move(parts);
    return ev;
}

//----------------------------------------------------------------------
// 2) 单粒子：聚类 + 估计
//----------------------------------------------------------------------
static void TestSingleTrack(DetectorSimulator& sim, double bfield_t) {
    std::printf("[2] 单粒子聚类与估计\n");

    const double ke = 15.0;
    Event ev = MakeEvent({MakeParticle(1, 1, ke, 70.0, 30.0)}, 0.10);

    PointCloud cloud;
    Check(sim.Simulate(ev, cloud, std::vector<int>{0}), "Simulate");
    Check(cloud.points.size() > 50, "点云点数足够");
    std::printf("    点云点数 = %zu\n", cloud.points.size());

    Clusterer clusterer;
    std::vector<Cluster> clusters;
    const int n = clusterer.Run(cloud, clusters);
    Check(n == 1, "单粒子应聚成 1 个簇");
    std::printf("    簇数 = %d\n", n);
    for (const auto& c : clusters) {
        std::printf("      label=%d 点数=%zu\n", c.label, c.points.size());
    }
    if (clusters.empty()) return;

    Estimator estimator;
    EstimateResult est;
    Check(estimator.Estimate(clusters[0], bfield_t, est), "Estimate");

    const double mass = MassData::GetNuclearMassNuclear(1, 1) * Units::u_MeV;
    const double p_true = std::sqrt(ke * (ke + 2.0 * mass));
    const double brho_true = p_true / 299.8;

    std::printf("    Bρ: 估计 %.3f Tm，真值 %.3f Tm，偏差 %.1f%%\n", est.brho, brho_true,
                100.0 * std::fabs(est.brho - brho_true) / brho_true);
    std::printf("    polar=%.1f°  vertex=(%.1f, %.1f, %.1f) mm  R=%.1f mm  dEdx=%.1f\n",
                est.polar / Units::DEG2RAD, est.vertex_x, est.vertex_y, est.vertex_z,
                est.radius, est.dEdx);

    Check(est.ok, "估计成功");
    Check(std::fabs(est.brho - brho_true) / brho_true < 0.30, "Bρ 偏差 < 30%");
    Check(std::fabs(est.polar / Units::DEG2RAD - 70.0) < 10.0, "极角偏差 < 10°");
    Check(std::fabs(est.vertex_z - 100.0) < 10.0, "顶点 z 偏差 < 10 mm");
    Check(est.dEdx > 0.0, "dE/dx 为正");
}

//----------------------------------------------------------------------
// 3) 人为把一条径迹切成两段，join 应把它们接回一个簇
//----------------------------------------------------------------------
static void TestSplitAndJoin(DetectorSimulator& sim) {
    std::printf("[3] 切碎后合并\n");

    Event ev = MakeEvent({MakeParticle(1, 1, 15.0, 70.0, 30.0)}, 0.10);

    PointCloud cloud;
    Check(sim.Simulate(ev, cloud, std::vector<int>{0}), "Simulate");
    if (cloud.points.size() < 80) {
        Check(false, "点云点数不足");
        return;
    }

    const double z_mid = 0.5 * (cloud.points.front().z + cloud.points.back().z);
    PointCloud lower, upper;
    lower.event = cloud.event;
    upper.event = cloud.event;
    for (const auto& q : cloud.points) {
        if (q.z < z_mid - 8.0) lower.points.push_back(q);
        else if (q.z > z_mid + 8.0) upper.points.push_back(q);
    }
    std::printf("    原始 %zu 点 -> 下段 %zu / 上段 %zu\n", cloud.points.size(),
                lower.points.size(), upper.points.size());

    Clusterer clusterer;
    ClusterParams legacy; legacy.use_hdbscan = false;
    clusterer.SetParams(legacy);
    std::vector<Cluster> c1, c2;
    clusterer.Run(lower, c1);
    clusterer.Run(upper, c2);
    Check(c1.size() == 1 && c2.size() == 1, "两段各自成簇");
    if (c1.empty() || c2.empty()) return;

    PointCloud combined;
    combined.event = cloud.event;
    combined.points = lower.points;
    combined.points.insert(combined.points.end(), upper.points.begin(), upper.points.end());

    std::vector<Cluster> joined;
    const int n = clusterer.Run(combined, joined);
    std::printf("    合并后簇数 = %d（joins=%d，点数 %zu）\n", n, clusterer.LastJoinCount(),
                joined.empty() ? 0 : joined[0].points.size());
    Check(n == 1, "join 后应合成 1 个簇");
    Check(clusterer.LastJoinCount() >= 1, "确实发生了一次合并");
}

//----------------------------------------------------------------------
// 4) 距离很远的簇不应被 join 误合并
//----------------------------------------------------------------------
static void TestNoFalseJoin(DetectorSimulator& sim) {
    std::printf("[4] 不误合并分隔的簇\n");

    Event ev = MakeEvent({MakeParticle(1, 1, 15.0, 70.0, 30.0)}, 0.10);
    PointCloud base;
    sim.Simulate(ev, base, std::vector<int>{0});
    if (base.points.size() < 80) return;

    PointCloud shifted;
    shifted.event = base.event + 100;
    for (auto q : base.points) {
        q.x += 300.0;
        shifted.points.push_back(q);
    }

    PointCloud combined;
    combined.event = base.event;
    combined.points = base.points;
    combined.points.insert(combined.points.end(), shifted.points.begin(), shifted.points.end());

    Clusterer clusterer;
    std::vector<Cluster> clusters;
    const int n = clusterer.Run(combined, clusters);
    std::printf("    簇数 = %d（joins=%d）\n", n, clusterer.LastJoinCount());
    Check(n == 2, "相距 300 mm 的两簇不应合并");
}

//----------------------------------------------------------------------
void testCluster() {
    std::printf("========== 聚类与估计测试 ==========\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    DetectorConfig cfg;
    DetectorSimulator sim;
    Check(sim.Init(cfg, &gas), "Init simulator");

    TestFitCircle();
    TestSingleTrack(sim, cfg.bfield_t);
    TestSplitAndJoin(sim);
    TestNoFalseJoin(sim);

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
