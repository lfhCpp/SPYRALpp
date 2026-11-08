//----------------------------------------------------------------------
// M5 轨迹网格 + 插值 + 拟合测试
//   root -l -b -q test/testSolver.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/TrackFitter.h"
#include "reco/TrackMesh.h"
#include "sim/DetectorConfig.h"
#include "sim/ParticleDynamics.h"
#include "sim/Simulator.h"
#include "sim/TargetDensity.h"

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

//----------------------------------------------------------------------
void testSolver() {
    std::printf("========== 轨迹求解测试 ==========\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    // ---- 1) 生成网格 ----
    MeshConfig mc;
    mc.ke_bins = 12;
    mc.ke_min = 5.0;
    mc.ke_max = 60.0;
    mc.polar_bins = 12;
    mc.polar_min = 30.0;
    mc.polar_max = 150.0;
    mc.n_time = 250;
    mc.dt = 1.0e-10;
    mc.z_ref = 0.5;

    TrackMesh mesh;
    Check(mesh.Build(mc, 1, 1, gas, 2.85, 60000.0), "Build mesh");
    std::printf("    网格：%d ke x %d polar x %d time\n", mc.ke_bins, mc.polar_bins,
                mc.n_time);

    // ---- 2) 存盘/读回 ----
    const std::string path = "/tmp/rk4_trackmesh_test.bin";
    Check(mesh.Save(path), "Save mesh");
    TrackMesh loaded;
    Check(loaded.Load(path), "Load mesh");

    std::vector<double> a, b;
    Check(mesh.Interpolate(30.0, 70.0 * Units::DEG2RAD, a), "Interpolate(原始)");
    Check(loaded.Interpolate(30.0, 70.0 * Units::DEG2RAD, b), "Interpolate(读回)");
    double max_diff = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        max_diff = std::max(max_diff, std::fabs(a[i] - b[i]));
    }
    Check(max_diff < 1e-9, "存盘前后一致");
    gSystem->Unlink(path.c_str());

    // ---- 3) 插值精度：与直接积分比较终点 ----
    ParticleDynamics dyn;
    dyn.ConfigureTarget(1, 1, 2, gas.density());
    dyn.ConfigureFields(2.85, 60000.0);
    dyn.SetStopZ(-0.5, 1.5);
    dyn.SetStopRho(0.5);
    dyn.SetParticle(1, 1);

    const double vtx[3] = {0.0, 0.0, mc.z_ref};
    std::vector<TrackPoint> track;
    dyn.Propagate(15.0, vtx, 70.0 * Units::DEG2RAD, 0.0, track, mc.dt, mc.n_time);

    std::vector<double> interp;
    Check(mesh.Interpolate(15.0, 70.0 * Units::DEG2RAD, interp), "Interpolate(15 MeV)");

    const size_t last = std::min<size_t>(track.size(), mc.n_time) - 1;
    const double ix = interp[last * 3 + 0];
    const double iy = interp[last * 3 + 1];
    const double iz = interp[last * 3 + 2];
    const double dx = ix - track[last].x * 1000.0;
    const double dy = iy - track[last].y * 1000.0;
    const double dz = iz - track[last].z * 1000.0;
    const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    std::printf("    插值 vs 直接积分：终点偏差 %.2f mm\n", dist);
    Check(dist < 15.0, "插值终点偏差 < 15 mm");

    // ---- 4) 拟合一条模拟径迹 ----
    DetectorConfig cfg;
    DetectorSimulator sim;
    Check(sim.Init(cfg, &gas), "Init simulator");

    Event ev;
    ev.id = 1;
    ev.vx = 0.0;
    ev.vy = 0.0;
    ev.vz = 0.10;
    ev.particles = {MakeParticle(1, 1, 15.0, 70.0, 30.0)};

    PointCloud cloud;
    Check(sim.Simulate(ev, cloud, std::vector<int>{0}), "Simulate");
    std::printf("    点云点数 = %zu\n", cloud.points.size());

    Clusterer clusterer;
    std::vector<Cluster> clusters;
    const int n_clusters = clusterer.Run(cloud, clusters);
    Check(n_clusters >= 1, "聚类成功");
    if (clusters.empty()) {
        gSystem->Exit(1);
    }

    Estimator estimator;
    EstimateResult est;
    Check(estimator.Estimate(clusters[0], cfg.bfield_t, est), "Estimate");

    const double mass = MassData::GetNuclearMassNuclear(1, 1) * Units::u_MeV;
    const double p_true = std::sqrt(15.0 * (15.0 + 2.0 * mass));
    const double brho_true = p_true / 299.8;

    TrackFitter fitter;
    Check(fitter.Init(&mesh, mass, 1), "Init fitter");

    FitResult fit;
    Check(fitter.Fit(clusters[0].points, est, est.azimuthal, fit), "Fit");

    std::printf("    初值 Bρ = %.3f Tm（偏差 %.1f%%）\n", est.brho,
                100.0 * std::fabs(est.brho - brho_true) / brho_true);
    std::printf("    拟合 Bρ = %.3f ± %.3f Tm（真值 %.3f，偏差 %.1f%%）\n", fit.brho,
                fit.sigma_brho, brho_true,
                100.0 * std::fabs(fit.brho - brho_true) / brho_true);
    std::printf("    polar = %.1f°  vertex = (%.1f, %.1f, %.1f) mm  cost = %.2f mm\n",
                fit.polar / Units::DEG2RAD, fit.vx, fit.vy, fit.vz, fit.cost);

    Check(fit.ok, "拟合收敛");
    Check(std::fabs(fit.brho - brho_true) / brho_true < 0.10, "Bρ 偏差 < 10%");
    Check(std::fabs(fit.polar / Units::DEG2RAD - 70.0) < 5.0, "极角偏差 < 5°");
    Check(std::fabs(fit.vz - 100.0) < 15.0, "顶点 z 偏差 < 15 mm");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
