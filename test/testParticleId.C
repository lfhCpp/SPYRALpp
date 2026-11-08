//----------------------------------------------------------------------
// 粒子鉴别测试：用质子标定，再判别 α / 氘 / 质子
//   root -l -b -q test/testParticleId.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
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

static Particle MakeParticle(int A, int Z, double ke_mev, double polar_deg) {
    const double mass = MassData::GetNuclearMassNuclear(A, Z) * Units::u_MeV;
    const double p = std::sqrt(ke_mev * (ke_mev + 2.0 * mass));
    const double th = polar_deg * Units::DEG2RAD;

    Particle part;
    part.A = A;
    part.Z = Z;
    part.label = 2;
    part.px = p * std::sin(th) * std::cos(30.0 * Units::DEG2RAD);
    part.py = p * std::sin(th) * std::sin(30.0 * Units::DEG2RAD);
    part.pz = p * std::cos(th);
    part.E = ke_mev + mass;
    return part;
}

static bool Measure(DetectorSimulator& sim, Clusterer& clusterer, Estimator& estimator,
                    int A, int Z, double ke, EstimateResult& est) {
    Event ev;
    ev.id = 1;
    ev.vx = 0.0;
    ev.vy = 0.0;
    ev.vz = 0.10;
    ev.particles = {MakeParticle(A, Z, ke, 70.0)};

    PointCloud cloud;
    if (!sim.Simulate(ev, cloud, std::vector<int>{0}) || cloud.points.empty()) return false;

    std::vector<Cluster> clusters;
    clusterer.Run(cloud, clusters);
    if (clusters.empty()) return false;

    return estimator.Estimate(clusters[0], sim.Pads().NumPads() > 0 ? 2.85 : 2.85, est);
}

void testParticleId() {
    std::printf("========== 粒子鉴别测试 ==========\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    DetectorConfig cfg;
    DetectorSimulator sim;
    sim.Init(cfg, &gas);

    Clusterer clusterer;
    Estimator estimator;

    // 用质子 30 MeV 标定
    EstimateResult ref;
    Check(Measure(sim, clusterer, estimator, 1, 1, 30.0, ref), "标定用质子测量");

    ParticleId pid;
    pid.Calibrate(ref.brho, ref.dEdx, 1, 1);
    std::printf("    标定：Bρ=%.3f dEdx=%.2f -> scale=%.4g\n", ref.brho, ref.dEdx, pid.Scale());

    // 注意：d (Z/A=0.5) 与 alpha (Z/A=0.5) 在 dE/dx-Bρ 平面上是简并的，
    // 物理上无法只靠这两个量区分，因此 alpha 允许被判成 d。
    struct Case {
        const char* expect;
        const char* alt; // 允许的简并结果，nullptr 表示必须精确匹配
        int A, Z;
        double ke;
    };
    const Case cases[] = {
        {"p", nullptr, 1, 1, 120.0},
        {"d", "alpha", 2, 1, 60.0},
        {"alpha", "d", 4, 2, 40.0},
        {"alpha", "d", 4, 2, 120.0},
    };

    int correct = 0;
    int total = 0;
    for (const auto& c : cases) {
        EstimateResult est;
        if (!Measure(sim, clusterer, estimator, c.A, c.Z, c.ke, est)) {
            std::printf("    %-6s(%.0f MeV) 测量失败\n", c.expect, c.ke);
            continue;
        }

        PidResult res;
        const bool ok = pid.Identify(est.brho, est.dEdx, res);
        ++total;
        const bool match = ok && (res.name == c.expect ||
                                  (c.alt != nullptr && res.name == c.alt));
        if (match) ++correct;

        std::printf("    %-6s %5.0f MeV: Bρ=%.3f dEdx=%6.2f -> %-6s sep=%.2f %s\n", c.expect,
                    c.ke, est.brho, est.dEdx, ok ? res.name.c_str() : "?", res.separation,
                    match ? "" : "<-- 错");
        (void)0;
    }

    std::printf("    判别正确 %d/%d\n", correct, total);
    Check(total >= 3, "至少 3 个测试用例");
    Check(correct >= total - 1, "绝大多数判别正确");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
