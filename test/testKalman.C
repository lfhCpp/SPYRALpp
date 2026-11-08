//----------------------------------------------------------------------
// Kalman（UKF + RTS）拟合测试：与 TMinuit 网格拟合、与真值对比
//   root -l -b -q test/testKalman.C
//
// 说明：Kalman 用 Estimator 的 Bρ 作动量初值，所以只测 Estimator 能给出
// 合理初值的角度（θ=50/60/70°）。θ≈90° 的径迹 z 几乎不变，当前的
// Estimator（沿 z 的圆拟合）本身给不出初值，属于已知限制，见 README。
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/KalmanFitter.h"
#include "reco/TrackFitter.h"
#include "reco/TrackMesh.h"
#include "sim/DetectorConfig.h"
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

static Particle MakeParticle(int A, int Z, double ke, double polar) {
    const double mass = MassData::GetNuclearMassNuclear(A, Z) * Units::u_MeV;
    const double p = std::sqrt(ke * (ke + 2.0 * mass));
    const double th = polar * Units::DEG2RAD;

    Particle q;
    q.A = A;
    q.Z = Z;
    q.label = 2;
    q.px = p * std::sin(th) * std::cos(0.3);
    q.py = p * std::sin(th) * std::sin(0.3);
    q.pz = p * std::cos(th);
    q.E = ke + mass;
    return q;
}

void testKalman() {
    std::printf("========== Kalman(UKF+RTS) 拟合测试 ==========\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    // 给 Minuit 用的轨迹网格
    MeshConfig mc;
    mc.ke_bins = 8;
    mc.ke_min = 5.0;
    mc.ke_max = 60.0;
    mc.polar_bins = 8;
    mc.polar_min = 30.0;
    mc.polar_max = 150.0;
    mc.n_time = 120;
    TrackMesh mesh;
    Check(mesh.Build(mc, 1, 1, gas, 2.85, 60000.0), "建网格");

    const double mass = MassData::GetNuclearMassNuclear(1, 1) * Units::u_MeV;

    TrackFitter minuit;
    Check(minuit.Init(&mesh, mass, 1), "Minuit 拟合器初始化");

    KalmanConfig kcfg; // 默认：0.5 mm 测点、30 个量测点、能损 straggling 关
    KalmanFitter kalman;
    Check(kalman.Init(gas, 1, 1, mass, 2.85, 60000.0, kcfg), "Kalman 拟合器初始化");

    DetectorConfig dcfg;
    DetectorSimulator sim;
    Check(sim.Init(dcfg, &gas), "模拟器初始化");

    Clusterer clusterer;
    Estimator estimator;

    const double kes[] = {10.0, 20.0, 40.0};
    const double pols[] = {50.0, 60.0, 70.0};

    std::printf("  %6s %6s | %8s %8s %8s | %7s %7s\n", "KE", "θ", "Bρ_truth", "Bρ_kalman",
                "Bρ_minuit", "χ²/ndf", "cost");
    std::printf("  ------------------------------------------------------------------\n");

    int n_case = 0, n_ok_k = 0, n_ok_m = 0;
    double max_rel_k = 0.0, max_rel_m = 0.0, max_dpol = 0.0, max_dbrho = 0.0;
    int case_id = 0;

    for (double ke : kes) {
        for (double pol : pols) {
            ++case_id;
            Event ev;
            ev.id = case_id;
            ev.vx = 0.0;
            ev.vy = 0.0;
            ev.vz = 0.10;
            ev.particles = {MakeParticle(1, 1, ke, pol)};

            PointCloud cloud;
            if (!sim.Simulate(ev, cloud, std::vector<int>{0}) || cloud.points.empty()) continue;
            std::vector<Cluster> clusters;
            clusterer.Run(cloud, clusters);
            if (clusters.empty()) continue;

            // 取第一个能给初值的簇（模拟事件偶尔会多出一个噪声簇）
            int ic = -1;
            EstimateResult est;
            for (size_t i = 0; i < clusters.size(); ++i) {
                if (estimator.Estimate(clusters[i], 2.85, est) && est.brho > 0.0) {
                    ic = static_cast<int>(i);
                    break;
                }
            }
            if (ic < 0) continue;
            const std::vector<CloudPoint>& pts = clusters[static_cast<size_t>(ic)].points;

            const double p_truth = std::sqrt(ke * (ke + 2.0 * mass));
            const double brho_truth = p_truth / 1000.0 / (0.299792458 * 1.0);

            KalmanFitResult kr;
            const bool ok_k = kalman.Fit(pts, est, kr);
            FitResult mr;
            const bool ok_m = minuit.Fit(pts, est, est.azimuthal, mr);

            ++n_case;
            if (ok_k) ++n_ok_k;
            if (ok_m) ++n_ok_m;

            const double rel_k = ok_k ? std::fabs(kr.brho - brho_truth) / brho_truth : 9.99;
            const double rel_m = ok_m ? std::fabs(mr.brho - brho_truth) / brho_truth : 9.99;
            if (ok_k) max_rel_k = std::max(max_rel_k, rel_k);
            if (ok_m) max_rel_m = std::max(max_rel_m, rel_m);

            if (ok_k && ok_m) {
                max_dbrho = std::max(max_dbrho, std::fabs(kr.brho - mr.brho) / brho_truth);
                max_dpol = std::max(max_dpol,
                                    std::fabs(kr.polar - mr.polar) / Units::DEG2RAD);
            }

            const double chi2ndf = (kr.ndf > 0) ? kr.chi2 / kr.ndf : -1.0;
            std::printf("  %6.1f %6.1f | %8.4f %8.4f %8.4f | %7.2f %7.2f\n", ke, pol,
                        brho_truth, ok_k ? kr.brho : -1.0, ok_m ? mr.brho : -1.0, chi2ndf,
                        ok_k ? kr.cost : -1.0);

            Check(ok_k, "Kalman 拟合成功");
            Check(rel_k < 0.10, "Kalman Bρ 与真值差 < 10%");
            Check(ok_k && std::isfinite(chi2ndf) && chi2ndf > 0.0 && chi2ndf < 20.0,
                  "χ²/ndf 在合理区间");
            Check(ok_k && kr.n_hits >= kcfg.min_clusters, "有效量测数足够");
            if (ok_m) {
                Check(std::fabs(kr.polar - mr.polar) / Units::DEG2RAD < 10.0,
                      "Kalman 与 Minuit 极角差 < 10°");
            }
        }
    }

    std::printf("\n  汇总：%d 个算例，Kalman 成功 %d，Minuit 成功 %d\n", n_case, n_ok_k, n_ok_m);
    std::printf("  最大相对偏差：Kalman %.2f%%，Minuit %.2f%%\n", max_rel_k * 100.0,
                max_rel_m * 100.0);
    std::printf("  Kalman vs Minuit：max ΔBρ/Bρ = %.2f%%，max Δθ = %.2f°\n", max_dbrho * 100.0,
                max_dpol);
    Check(n_ok_k == n_case, "所有算例 Kalman 都收敛");
    Check(max_rel_k < 0.10, "Kalman Bρ 总体精度 < 10%");
    Check(max_dpol < 10.0, "Kalman 与 Minuit 极角一致");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
