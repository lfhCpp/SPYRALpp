//----------------------------------------------------------------------
// 端到端测试（真实波形路径）：
//   文本波形 -> 点云 -> 聚类/估计/拟合
//   root -l -b -q test/testWaveformPipeline.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/CloudIO.h"
#include "io/ResultIO.h"

#include "TSystem.h"

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
                 "[gas]\nA = 1\nZ = 1\nstoich = 2\npressure_torr = 600\n"
                 "\n[detector]\n"
                 "pad_dir = data\n"
                 "bfield_t = 2.85\n"
                 "efield_vm = 60000\n"
                 "\n[peaks]\n"
                 "separation = 50\n"
                 "prominence = 20\n"
                 "max_width = 50\n"
                 "threshold = 40\n"
                 "baseline_scale = 20\n"
                 "\n[solver]\n"
                 "A = 1\nZ = 1\n"
                 "ke_bins = 12\nke_min = 5\nke_max = 60\n"
                 "polar_bins = 12\npolar_min = 30\npolar_max = 150\n"
                 "n_time = 250\n");
    std::fclose(f);
}

void testWaveformPipeline() {
    std::printf("========== 波形路径端到端测试 ==========\n");

    const std::string cfg_path = "/tmp/rk4_wf.cfg";
    const std::string cloud_path = "/tmp/rk4_wf_clouds.root";
    const std::string fit_path = "/tmp/rk4_wf_fit.root";

    gSystem->Unlink(cloud_path.c_str());
    gSystem->Unlink(fit_path.c_str());
    WriteConfig(cfg_path);

    const std::string cmd_wf = "build/rk4_waveform " + cfg_path +
                               " data/event_waveform.txt " + cloud_path +
                               " 3 > /tmp/rk4_wf.log 2>&1";
    const int rc_wf = gSystem->Exec(cmd_wf.c_str());
    Check(rc_wf == 0, "rk4_waveform 退出码 0");
    Check(!gSystem->AccessPathName(cloud_path.c_str()), "点云文件已生成");

    std::vector<PointCloud> clouds;
    Check(ReadPointClouds(cloud_path, clouds), "读取点云");
    Check(!clouds.empty(), "点云非空");

    size_t total_points = 0;
    for (const auto& c : clouds) total_points += c.points.size();
    std::printf("    %zu 个事件，共 %zu 个点\n", clouds.size(), total_points);
    Check(clouds.size() >= 1, "至少 1 个事件");
    Check(total_points > 50, "点数合理");

    // 再用重建程序跑一遍（不比对真值，只验证管线不崩）
    const std::string cmd_ana = "build/rk4_analyze " + cloud_path + " " + fit_path +
                                " " + cfg_path + " > /tmp/rk4_wf_ana.log 2>&1";
    const int rc_ana = gSystem->Exec(cmd_ana.c_str());
    Check(rc_ana == 0, "rk4_analyze 退出码 0");
    Check(!gSystem->AccessPathName(fit_path.c_str()), "拟合结果已生成");

    std::vector<FitRecord> records;
    Check(ReadFitResults(fit_path, records), "读取拟合结果");
    std::printf("    拟合记录 %zu 条\n", records.size());

    gSystem->Unlink(cfg_path.c_str());
    gSystem->Unlink(cloud_path.c_str());
    gSystem->Unlink(fit_path.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
