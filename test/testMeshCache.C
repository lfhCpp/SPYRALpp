//----------------------------------------------------------------------
// 轨迹网格缓存测试：存盘/读回/参数匹配
//   root -l -b -q test/testMeshCache.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Units.h"
#include "reco/TrackMesh.h"
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

void testMeshCache() {
    std::printf("========== 网格缓存测试 ==========\n");

    std::vector<Element> comp;
    Element h;
    h.a = 1;
    h.z = 1;
    h.s = 2;
    comp.push_back(h);
    GasTarget gas(comp, 600.0);

    MeshConfig mc;
    mc.ke_bins = 8;
    mc.ke_min = 5.0;
    mc.ke_max = 60.0;
    mc.polar_bins = 8;
    mc.polar_min = 30.0;
    mc.polar_max = 150.0;
    mc.n_time = 120;
    mc.z_ref = 0.5;

    TrackMesh mesh;
    Check(mesh.Build(mc, 1, 1, gas, 2.85, 60000.0), "Build mesh");
    Check(mesh.ParticleA() == 1 && mesh.ParticleZ() == 1, "记录了粒子种类");

    // 同参数：应匹配
    Check(mesh.MatchesConfig(mc, 1, 1), "同参数匹配");

    // 粒子不同：不匹配
    Check(!mesh.MatchesConfig(mc, 4, 2), "不同粒子不匹配");

    // 网格参数不同：不匹配
    MeshConfig other = mc;
    other.ke_bins = 10;
    Check(!mesh.MatchesConfig(other, 1, 1), "ke_bins 变化不匹配");

    other = mc;
    other.z_ref = 0.4;
    Check(!mesh.MatchesConfig(other, 1, 1), "z_ref 变化不匹配");

    // 存盘/读回
    const std::string path = "/tmp/rk4_mesh_cache_test.bin";
    Check(mesh.Save(path), "Save");

    TrackMesh loaded;
    Check(loaded.Load(path), "Load");
    Check(loaded.MatchesConfig(mc, 1, 1), "读回后参数匹配");
    Check(loaded.ParticleA() == 1 && loaded.ParticleZ() == 1, "读回粒子种类");

    std::vector<double> a, b;
    Check(mesh.Interpolate(20.0, 70.0 * Units::DEG2RAD, a), "Interpolate 原始");
    Check(loaded.Interpolate(20.0, 70.0 * Units::DEG2RAD, b), "Interpolate 读回");

    double max_diff = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        max_diff = std::max(max_diff, std::fabs(a[i] - b[i]));
    }
    Check(max_diff < 1e-9, "读回数据一致");

    gSystem->Unlink(path.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
