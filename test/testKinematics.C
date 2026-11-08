//----------------------------------------------------------------------
// M1 运动学生成器测试（ROOT 宏，无需单独编译）
//   root -l -b -q test/testKinematics.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/MassData.h"
#include "core/Reaction.h"
#include "core/Units.h"
#include "sim/KinematicsGenerator.h"
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

static void CheckClose(double a, double b, double tol, const char* msg) {
    if (std::fabs(a - b) > tol) {
        std::printf("  [FAIL] %s: %.6g vs %.6g (tol %.3g)\n", msg, a, b, tol);
        ++g_fail;
    }
}

//----------------------------------------------------------------------
// 1) 四动量守恒
//----------------------------------------------------------------------
static void TestConservation() {
    std::printf("[1] 四动量守恒\n");

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 12.8 * 12.0247;
    cfg.apply_target_loss = false; // 关掉能损，初态可精确写出

    KinematicsGenerator gen;
    Check(gen.Init(cfg, nullptr, 20240101), "Init");

    const int n = 500;
    for (int i = 0; i < n; ++i) {
        Event ev;
        if (!gen.Generate(ev)) {
            Check(false, "Generate 返回 false");
            return;
        }

        const double mb = gen.Reaction().BeamMass();
        const double p_beam = std::sqrt(cfg.beam_ke_mev * (cfg.beam_ke_mev + 2.0 * mb));

        // 事件里存的是 [靶, 束流, 出射, 剩余]：
        // 初态 = 0+1，末态 = 2+3，不能把四个一起加
        FourMomentum initial;
        initial.Add(ev.particles[0]);
        initial.Add(ev.particles[1]);

        FourMomentum final_state;
        final_state.Add(ev.particles[2]);
        final_state.Add(ev.particles[3]);

        CheckClose(final_state.px, initial.px, 1e-6, "px 守恒");
        CheckClose(final_state.py, initial.py, 1e-6, "py 守恒");
        CheckClose(final_state.pz, initial.pz, 1e-6 * p_beam, "pz 守恒");
        CheckClose(final_state.E, initial.E, 1e-6 * initial.E, "E 守恒");
    }
}

//----------------------------------------------------------------------
// 2) 每个粒子在壳
//----------------------------------------------------------------------
static void TestOnShell() {
    std::printf("[2] 粒子在壳\n");

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 150.0;
    cfg.ex_centroid = 2.0; // 有激发能，检查剩余核质量
    cfg.ex_fwhm = 0.0;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, nullptr, 7), "Init");

    for (int i = 0; i < 200; ++i) {
        Event ev;
        if (!gen.Generate(ev)) {
            Check(false, "Generate 返回 false");
            return;
        }
        CheckClose(ev.particles[2].Mass(), gen.Reaction().EjectileMass(), 1e-4, "出射核质量");
        CheckClose(ev.particles[3].Mass(), gen.Reaction().ResidualMass() + 2.0, 1e-4,
                   "剩余核质量（含激发能）");
    }
}

//----------------------------------------------------------------------
// 3) 顶点与角度分布
//----------------------------------------------------------------------
static void TestDistributions() {
    std::printf("[3] 顶点与角度分布\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 12.8 * 12.0247;
    cfg.rho_sigma_m = 0.007;
    cfg.z_min_m = 0.0;
    cfg.z_max_m = 1.0;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, &gas, 99), "Init(带靶)");

    const int n = 4000;
    double sum_z = 0.0, sum_rho = 0.0;
    for (int i = 0; i < n; ++i) {
        Event ev;
        if (!gen.Generate(ev)) {
            Check(false, "Generate 返回 false");
            return;
        }
        sum_z += ev.vz;
        sum_rho += std::hypot(ev.vx, ev.vy);
    }

    CheckClose(sum_z / n, 0.5, 0.02, "顶点 z 均值");
    // 半正态分布均值 = sigma * sqrt(2/pi)
    CheckClose(sum_rho / n, 0.007 * std::sqrt(2.0 / Units::PI), 0.001, "顶点 rho 均值");

    // 质心系极角分布单独测：PolarUniform 应对 cos(theta) 均匀
    core::Rng rng(4321);
    PolarUniform polar(0.0, Units::PI);
    double sum_cos_cm = 0.0;
    const int m = 20000;
    for (int i = 0; i < m; ++i) {
        sum_cos_cm += std::cos(polar.SamplePolar(rng));
    }
    CheckClose(sum_cos_cm / m, 0.0, 0.01, "cos(theta_CM) 均值（立体角均匀）");
}

//----------------------------------------------------------------------
// 4) 靶中能损
//----------------------------------------------------------------------
static void TestEnergyLoss() {
    std::printf("[4] 靶中能损\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 12.8 * 12.0247;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, &gas, 3), "Init");

    const double e0 = gen.BeamEnergyAt(0.0);
    const double e1 = gen.BeamEnergyAt(1.0);
    const double e05 = gen.BeamEnergyAt(0.5);

    CheckClose(e0, cfg.beam_ke_mev, 1e-6, "z=0 时无损失");
    Check(e1 < e0, "穿过 1 m 后能量下降");
    Check(e05 < e0 && e05 > e1, "能量随深度单调下降");
    std::printf("    E(0)=%.3f  E(0.5)=%.3f  E(1.0)=%.3f MeV\n", e0, e05, e1);
}

//----------------------------------------------------------------------
// 5) 激发能分布
//----------------------------------------------------------------------
static void TestExcitation() {
    std::printf("[5] 激发能分布\n");

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 150.0;
    cfg.ex_centroid = 2.0;
    cfg.ex_fwhm = 0.1;
    cfg.apply_target_loss = false;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, nullptr, 11), "Init");

    const int n = 3000;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; ++i) {
        Event ev;
        if (!gen.Generate(ev)) {
            Check(false, "Generate 返回 false");
            return;
        }
        const double ex = ev.particles[3].Mass() - gen.Reaction().ResidualMass();
        sum += ex;
        sum2 += ex * ex;
    }
    const double mean = sum / n;
    const double var = sum2 / n - mean * mean;
    CheckClose(mean, 2.0, 0.005, "激发能均值");
    CheckClose(std::sqrt(var), 0.1 / 2.35482, 0.005, "激发能 sigma");
}

//----------------------------------------------------------------------
// 6) 阈值以上的激发能应被拒绝
//----------------------------------------------------------------------
static void TestThreshold() {
    std::printf("[6] 能量不允许时拒绝采样\n");

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 150.0;
    cfg.ex_centroid = 500.0; // 远超运动学上限
    cfg.ex_fwhm = 1.0;
    cfg.apply_target_loss = false;
    cfg.max_trials = 50;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, nullptr, 5), "Init");

    Event ev;
    Check(!gen.Generate(ev), "不可能的激发能应返回 false");
}

//----------------------------------------------------------------------
// 7) 二体衰变
//----------------------------------------------------------------------
static void TestDecay() {
    std::printf("[7] 二体衰变 8Be -> 2alpha\n");

    TwoBodyDecay decay;
    Check(decay.Init(8, 4, 4, 2, 4, 5), "Init 8Be->4He+4He");

    // 静止的 8Be
    Particle parent;
    parent.A = 8;
    parent.Z = 4;
    parent.label = 3;
    parent.px = parent.py = parent.pz = 0.0;
    parent.E = MassData::GetNuclearMassNuclear(8, 4) * Units::u_MeV;

    Particle f1, f2;
    Check(decay.Calculate(parent, 0.0, 0.0, 0.0, f1, f2), "Calculate");

    CheckClose(f1.E + f2.E, parent.E, 1e-6 * parent.E, "衰变能量守恒");
    CheckClose(f1.px + f2.px, 0.0, 1e-6, "衰变 px 守恒");
    CheckClose(f1.pz + f2.pz, 0.0, 1e-6, "衰变 pz 守恒");
    CheckClose(f1.Mass(), MassData::GetNuclearMassNuclear(4, 2) * Units::u_MeV, 1e-4,
               "alpha 质量");
    // 8Be 的 Q 值约 92 keV，每个 alpha 约 46 keV
    CheckClose(f1.Kinetic(), 0.046, 0.01, "alpha 动能（Q/2）");

    // 飞行中的衰变：检查 boost
    Particle moving = parent;
    const double p = 200.0;
    moving.pz = p;
    moving.E = std::sqrt(p * p + parent.E * parent.E);

    Particle g1, g2;
    Check(decay.Calculate(moving, 0.0, 0.0, 0.0, g1, g2), "Calculate（飞行中）");
    CheckClose(g1.pz + g2.pz, p, 1e-6 * p, "飞行中 pz 守恒");
    CheckClose(g1.E + g2.E, moving.E, 1e-6 * moving.E, "飞行中 E 守恒");
    Check(g1.pz > 100.0, "boost 后前向 alpha 动量显著增加");
}

//----------------------------------------------------------------------
void testKinematics() {
    std::printf("========== 运动学生成器测试 ==========\n");

    if (!MassData::IsLoaded()) {
        std::printf("[FAIL] 质量表未加载\n");
        gSystem->Exit(1);
    }

    TestConservation();
    TestOnShell();
    TestDistributions();
    TestEnergyLoss();
    TestExcitation();
    TestThreshold();
    TestDecay();

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
