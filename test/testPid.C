//----------------------------------------------------------------------
// 验证 dE/dx 是否可用于粒子鉴别：同一 Bρ 附近，α 的 dE/dx 应显著高于质子
//   root -l -b -q test/testPid.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "sim/DetectorConfig.h"
#include "sim/Simulator.h"
#include "sim/TargetDensity.h"

#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rk4;

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

struct Probe {
    const char* name;
    int A, Z;
    double ke;
    double polar;
};

void testPid() {
    std::printf("========== dE/dx 粒子鉴别可行性 ==========\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1; h.z = 1; h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    DetectorConfig cfg;
    DetectorSimulator sim;
    sim.Init(cfg, &gas);

    Clusterer clusterer;
    Estimator estimator;

    const Probe probes[] = {
        {"p   30MeV", 1, 1, 30.0, 70.0},
        {"p  120MeV", 1, 1, 120.0, 70.0},
        {"alpha 40MeV", 4, 2, 40.0, 70.0},
        {"alpha 120MeV", 4, 2, 120.0, 70.0},
    };

    std::printf("    %-12s %8s %8s %8s %8s\n", "粒子", "Bρ", "dEdx", "弧长", "点数");
    for (const auto& pr : probes) {
        Event ev;
        ev.id = 1;
        ev.vx = 0.0; ev.vy = 0.0; ev.vz = 0.10;
        ev.particles = {MakeParticle(pr.A, pr.Z, pr.ke, pr.polar)};

        PointCloud cloud;
        if (!sim.Simulate(ev, cloud, std::vector<int>{0}) || cloud.points.empty()) {
            std::printf("    %-12s 模拟失败\n", pr.name);
            continue;
        }

        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);
        if (clusters.empty()) {
            std::printf("    %-12s 聚类失败（%zu 点）\n", pr.name, cloud.points.size());
            continue;
        }

        EstimateResult est;
        if (!estimator.Estimate(clusters[0], cfg.bfield_t, est)) {
            std::printf("    %-12s 估计失败\n", pr.name);
            continue;
        }

        std::printf("    %-12s %8.3f %8.2f %8.1f %8zu\n", pr.name, est.brho, est.dEdx,
                    est.arclength, clusters[0].points.size());
    }

    gSystem->Exit(0);
}
