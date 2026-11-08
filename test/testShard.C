//----------------------------------------------------------------------
// 多进程分片测试：
//   1) 分片处理 + 合并，结果与单进程逐位一致
//   2) 分片之间不重不漏
//   root -l -b -q test/testShard.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/ResultIO.h"

#include "TSystem.h"

#include <algorithm>
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

static void WriteConfig(const std::string& path, int threads) {
    FILE* f = std::fopen(path.c_str(), "w");
    std::fprintf(f,
                 "[run]\nthreads = %d\n\n"
                 "[gas]\nA = 1\nZ = 1\nstoich = 2\npressure_torr = 600\n\n"
                 "[detector]\npad_dir = data\nbfield_t = 2.85\nefield_vm = 60000\n\n"
                 "[solver]\nA = 1\nZ = 1\n"
                 "ke_bins = 10\nke_min = 2\nke_max = 120\n"
                 "polar_bins = 10\npolar_min = 20\npolar_max = 160\n"
                 "n_time = 200\n"
                 "mesh_cache = /tmp/ts_mesh.bin\nmesh_shared = 1\n",
                 threads);
    std::fclose(f);
}

void testShard() {
    std::printf("========== 多进程分片测试 ==========\n");

    const std::string cfg = "/tmp/ts.cfg";
    const std::string clouds = "/tmp/ts_clouds.root";
    WriteConfig(cfg, 1);

    gSystem->Unlink("/tmp/ts_mesh.bin");
    gSystem->Unlink(clouds.c_str());

    // 生成 8 个事件（模拟程序会顺带写出真值文件，这里不需要）
    const int rc_sim = gSystem->Exec(
        "build/rk4_simulate /tmp/ts.cfg /tmp/ts_clouds.root 8 > /tmp/ts_sim.log 2>&1");
    Check(rc_sim == 0, "模拟生成点云");
    Check(!gSystem->AccessPathName(clouds.c_str()), "点云文件存在");

    // 单进程基线（顺带生成网格缓存）
    gSystem->Unlink("/tmp/ts_base.root");
    const int rc_base = gSystem->Exec(
        "build/rk4_analyze /tmp/ts_clouds.root /tmp/ts_base.root /tmp/ts.cfg "
        "> /tmp/ts_base.log 2>&1");
    Check(rc_base == 0, "单进程分析成功");
    Check(!gSystem->AccessPathName("/tmp/ts_mesh.bin"), "网格缓存已生成");

    // 4 进程分片（共享映射同一网格）
    for (int i = 0; i < 4; ++i) {
        gSystem->Unlink(("/tmp/ts_shard_" + std::to_string(i) + ".root").c_str());
    }
    const int rc_shard = gSystem->Exec(
        "for i in 0 1 2 3; do build/rk4_analyze /tmp/ts_clouds.root "
        "/tmp/ts_shard_$i.root /tmp/ts.cfg $i 4 > /tmp/ts_shard_$i.log 2>&1 & "
        "done; wait");
    Check(rc_shard == 0, "分片进程正常退出");

    gSystem->Unlink("/tmp/ts_merged.root");
    const int rc_merge = gSystem->Exec(
        "build/rk4_merge /tmp/ts_merged.root /tmp/ts_shard_0.root /tmp/ts_shard_1.root "
        "/tmp/ts_shard_2.root /tmp/ts_shard_3.root > /tmp/ts_merge.log 2>&1");
    Check(rc_merge == 0, "合并成功");

    std::vector<FitRecord> base, merged;
    Check(ReadFitResults("/tmp/ts_base.root", base), "读取基线");
    Check(ReadFitResults("/tmp/ts_merged.root", merged), "读取合并结果");

    std::printf("    基线 %zu 条，分片合并 %zu 条\n", base.size(), merged.size());
    Check(!base.empty(), "基线非空");
    Check(base.size() == merged.size(), "条数一致");

    if (base.size() == merged.size() && !base.empty()) {
        int mismatch = 0;
        double max_db = 0.0;
        for (size_t i = 0; i < base.size(); ++i) {
            if (base[i].event != merged[i].event || base[i].cluster != merged[i].cluster) {
                ++mismatch;
                continue;
            }
            max_db = std::max(max_db, std::fabs(base[i].brho - merged[i].brho));
            if (base[i].pid_index != merged[i].pid_index) ++mismatch;
        }
        Check(mismatch == 0, "事件/簇/粒子种类一致");
        Check(max_db < 1e-9, "Bρ 逐位一致");
        std::printf("    最大 Bρ 差异 = %.3g\n", max_db);
    }

    size_t sum = 0;
    for (int i = 0; i < 4; ++i) {
        std::vector<FitRecord> part;
        if (ReadFitResults("/tmp/ts_shard_" + std::to_string(i) + ".root", part)) {
            sum += part.size();
        }
    }
    Check(sum == merged.size(), "分片之和等于合并结果（无重复）");

    for (int i = 0; i < 4; ++i) {
        gSystem->Unlink(("/tmp/ts_shard_" + std::to_string(i) + ".root").c_str());
    }
    gSystem->Unlink("/tmp/ts_base.root");
    gSystem->Unlink("/tmp/ts_merged.root");
    gSystem->Unlink(clouds.c_str());
    gSystem->Unlink(cfg.c_str());
    gSystem->Unlink("/tmp/ts_mesh.bin");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
