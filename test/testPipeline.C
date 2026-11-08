//----------------------------------------------------------------------
// 端到端测试：配置 -> 模拟 -> 写盘 -> 重建 -> 与真值比较
//   root -l -b -q test/testPipeline.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "core/Units.h"
#include "io/ResultIO.h"
#include "io/TreeIO.h"

#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

static void WriteConfig(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "w");
    std::fprintf(f,
                 "[run]\n"
                 "seed = 7\n"
                 "n_events = 5\n"
                 "\n"
                 "[gas]\n"
                 "A = 1\n"
                 "Z = 1\n"
                 "stoich = 2\n"
                 "pressure_torr = 600\n"
                 "\n"
                 "[reaction]\n"
                 "beam_A = 12\n"
                 "beam_Z = 4\n"
                 "target_A = 1\n"
                 "target_Z = 1\n"
                 "ejectile_A = 1\n"
                 "ejectile_Z = 1\n"
                 "\n"
                 "[beam]\n"
                 "ke_mev = 40\n"
                 "ex_centroid = 0\n"
                 "ex_fwhm = 0\n"
                 "theta_min = 50\n"
                 "theta_max = 110\n"
                 "\n"
                 "[vertex]\n"
                 "z_min = 0.05\n"
                 "z_max = 0.45\n"
                 "rho_sigma = 0.005\n"
                 "\n"
                 "[detector]\n"
                 "mode = ejectile\n"
                 "bfield_t = 2.85\n"
                 "efield_vm = 60000\n"
                 "\n"
                 "[solver]\n"
                 "A = 1\n"
                 "Z = 1\n"
                 "ke_bins = 16\n"
                 "ke_min = 2\n"
                 "ke_max = 160\n"
                 "polar_bins = 16\n"
                 "polar_min = 20\n"
                 "polar_max = 160\n"
                 "n_time = 300\n");
    std::fclose(f);
}

void testPipeline() {
    std::printf("========== 端到端管线测试 ==========\n");

    const std::string cfg_path = "/tmp/rk4_e2e.cfg";
    const std::string cloud_path = "/tmp/rk4_e2e_clouds.root";
    const std::string kin_path = "/tmp/rk4_e2e_clouds_kin.root";
    const std::string fit_path = "/tmp/rk4_e2e_fit.root";

    gSystem->Unlink(cloud_path.c_str());
    gSystem->Unlink(kin_path.c_str());
    gSystem->Unlink(fit_path.c_str());

    WriteConfig(cfg_path);

    const std::string cmd_sim =
        "build/rk4_simulate " + cfg_path + " " + cloud_path + " 5 > /tmp/rk4_sim.log 2>&1";
    const int rc_sim = gSystem->Exec(cmd_sim.c_str());
    Check(rc_sim == 0, "rk4_simulate 退出码 0");
    Check(!gSystem->AccessPathName(cloud_path.c_str()), "点云文件已生成");
    Check(!gSystem->AccessPathName(kin_path.c_str()), "真值文件已生成");

    const std::string cmd_ana =
        "build/rk4_analyze " + cloud_path + " " + fit_path + " " + cfg_path +
        " > /tmp/rk4_ana.log 2>&1";
    const int rc_ana = gSystem->Exec(cmd_ana.c_str());
    Check(rc_ana == 0, "rk4_analyze 退出码 0");
    Check(!gSystem->AccessPathName(fit_path.c_str()), "拟合结果已生成");

    std::vector<FitRecord> records;
    Check(ReadFitResults(fit_path, records), "读取拟合结果");
    Check(!records.empty(), "拟合结果非空");
    if (records.empty()) {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }

    std::vector<Event> events;
    Check(ReadKinematics(kin_path, events), "读取真值事件");
    Check(events.size() > 0, "真值事件非空");

    // 逐事件比较 Bρ
    int compared = 0, good = 0;
    double worst = 0.0;
    for (const auto& rec : records) {
        const Event* truth = nullptr;
        for (const auto& ev : events) {
            if (ev.id == rec.event) {
                truth = &ev;
                break;
            }
        }
        if (truth == nullptr || truth->particles.size() < 3) continue;

        const Particle& p = truth->particles[2]; // 出射质子
        const double p_total = p.P();
        const double brho_true = p_total / 299.8;

        const double dev = std::fabs(rec.brho - brho_true) / brho_true;
        worst = std::max(worst, dev);
        ++compared;
        if (dev < 0.15) ++good;
    }

    std::printf("    拟合记录 %zu 条，与真值比较 %d 条\n", records.size(), compared);
    std::printf("    最大 Bρ 偏差 %.1f%%\n", 100.0 * worst);

    Check(compared >= 3, "至少比较 3 个事件");
    Check(compared > 0 && good >= compared - 1, "绝大多数事件 Bρ 偏差 < 15%");

    gSystem->Unlink(cfg_path.c_str());
    gSystem->Unlink(cloud_path.c_str());
    gSystem->Unlink(kin_path.c_str());
    gSystem->Unlink(fit_path.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
