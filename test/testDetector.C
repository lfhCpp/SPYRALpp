//----------------------------------------------------------------------
// M2 探测器数字化测试
//   root -l -b -q test/testDetector.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "io/CloudIO.h"
#include "sim/DetectorConfig.h"
#include "sim/KinematicsGenerator.h"
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

static void CheckClose(double a, double b, double tol, const char* msg) {
    if (std::fabs(a - b) > tol) {
        std::printf("  [FAIL] %s: %.6g vs %.6g (tol %.3g)\n", msg, a, b, tol);
        ++g_fail;
    }
}

void testDetector() {
    std::printf("========== 探测器数字化测试 ==========\n");

    // H2 靶，600 Torr
    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    DetectorConfig cfg;
    DetectorSimulator sim;
    Check(sim.Init(cfg, &gas), "Init");
    Check(sim.Pads().NumPads() == 10240, "pad 数应为 10240");

    // 生成一条 12Be(p,p) 事件
    GeneratorConfig gcfg;
    gcfg.reaction.beam_A = 12;
    gcfg.reaction.beam_Z = 4;
    gcfg.reaction.target_A = 1;
    gcfg.reaction.target_Z = 1;
    gcfg.reaction.ejectile_A = 1;
    gcfg.reaction.ejectile_Z = 1;
    gcfg.beam_ke_mev = 150.0;
    gcfg.apply_target_loss = false;

    KinematicsGenerator gen;
    Check(gen.Init(gcfg, nullptr, 2025), "Kinematics Init");

    Event ev;
    Check(gen.Generate(ev), "Generate");
    ev.vx = 0.0;
    ev.vy = 0.0;
    ev.vz = 0.5; // 顶点放中间

    // 只模拟出射质子
    PointCloud pc;
    Check(sim.Simulate(ev, pc, std::vector<int>{2}), "Simulate");
    Check(!pc.points.empty(), "点云非空");

    int bad_pad = 0, bad_tb = 0, beam_hit = 0;
    double z_min = 1e9, z_max = -1e9, amp_max = 0.0;
    for (const auto& p : pc.points) {
        if (p.pad < 0 || p.pad >= sim.Pads().NumPads()) ++bad_pad;
        if (p.tb < 0 || p.tb >= cfg.num_tb) ++bad_tb;
        if (p.pad >= 0 && sim.Pads().IsBeamPad(p.pad)) ++beam_hit;
        z_min = std::min(z_min, p.z);
        z_max = std::max(z_max, p.z);
        amp_max = std::max(amp_max, p.amplitude);
    }

    Check(bad_pad == 0, "pad 编号合法");
    Check(bad_tb == 0, "时间桶合法");
    Check(beam_hit == 0, "没有命中 beam pad");
    Check(z_min >= -1e-6 && z_max <= cfg.length_m * 1000.0 + 1e-6, "z 在灵敏体积内");
    Check(amp_max <= 4095.0 + 1e-6, "ADC 不超 11 bit");
    Check(sim.TotalElectrons() > 1.0e6, "沉积电子数合理");
    std::printf("    点数=%zu  电子总数=%.3e  z=[%.1f, %.1f] mm  amp_max=%.0f\n",
                pc.points.size(), sim.TotalElectrons(), z_min, z_max, amp_max);

    // 剩余核（12Be）也应能产生点云
    PointCloud pc2;
    Check(sim.Simulate(ev, pc2, std::vector<int>{3}), "Simulate(剩余核)");
    Check(!pc2.points.empty(), "剩余核点云非空");

    // 阈值：把 ADC 阈值抬到极高，点云应为空
    DetectorConfig strict = cfg;
    strict.adc_threshold = 1.0e9;
    DetectorSimulator sim_strict;
    Check(sim_strict.Init(strict, &gas), "Init(严格阈值)");
    PointCloud pc3;
    Check(sim_strict.Simulate(ev, pc3, std::vector<int>{2}), "Simulate(严格阈值)");
    Check(pc3.points.empty(), "高于阈值时点云为空");

    // 同种子可复现
    DetectorSimulator sim_a, sim_b;
    Check(sim_a.Init(cfg, &gas), "Init(a)");
    Check(sim_b.Init(cfg, &gas), "Init(b)");
    PointCloud pa, pb;
    sim_a.Simulate(ev, pa, std::vector<int>{2});
    sim_b.Simulate(ev, pb, std::vector<int>{2});
    Check(pa.points.size() == pb.points.size(), "同种子点数一致");
    if (!pa.points.empty() && !pb.points.empty()) {
        CheckClose(pa.points[0].amplitude, pb.points[0].amplitude, 1e-9, "同种子首点一致");
    }

    // 点云写出 / 读回
    std::vector<PointCloud> clouds;
    clouds.push_back(pc);
    clouds.push_back(pc2);
    const std::string cloud_path = "/tmp/rk4_clouds_test.root";
    Check(WritePointClouds(cloud_path, clouds), "WritePointClouds");

    std::vector<PointCloud> back;
    Check(ReadPointClouds(cloud_path, back), "ReadPointClouds");
    Check(back.size() == clouds.size(), "点云事件数一致");
    if (back.size() == clouds.size()) {
        Check(back[0].points.size() == pc.points.size(), "点数一致");
        if (!back[0].points.empty()) {
            CheckClose(back[0].points[0].amplitude, pc.points[0].amplitude, 1e-9,
                       "首点幅度一致");
        }
    }
    gSystem->Unlink(cloud_path.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
