//----------------------------------------------------------------------
// HDF5 输出测试（未启用 HDF5 时自动跳过）
//   root -l -b -q test/testHdf5.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/H5Writer.h"
#include "sim/KinematicsGenerator.h"

#include "TSystem.h"

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

void testHdf5() {
    std::printf("========== HDF5 输出测试 ==========\n");

    if (!Hdf5Available()) {
        std::printf("========== SKIP（编译时未启用 -DRK4_WITH_HDF5=ON） ==========\n");
        gSystem->Exit(0);
    }

    GeneratorConfig cfg;
    cfg.reaction.beam_A = 12;
    cfg.reaction.beam_Z = 4;
    cfg.reaction.target_A = 1;
    cfg.reaction.target_Z = 1;
    cfg.reaction.ejectile_A = 1;
    cfg.reaction.ejectile_Z = 1;
    cfg.beam_ke_mev = 150.0;
    cfg.apply_target_loss = false;

    KinematicsGenerator gen;
    Check(gen.Init(cfg, nullptr, 77), "Init");

    std::vector<Event> events;
    for (int i = 0; i < 10; ++i) {
        Event ev;
        Check(gen.Generate(ev), "Generate");
        events.push_back(ev);
    }

    const std::string path = "/tmp/rk4_kinematics_hdf5_test.h5";
    Check(WriteKinematicsHdf5(path, events, 6), "WriteKinematicsHdf5");

    Check(!gSystem->AccessPathName(path.c_str()), "文件已生成");

    gSystem->Unlink(path.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
