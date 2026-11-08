//----------------------------------------------------------------------
// 多线程模拟的正确性：
//   - 事件编号全局唯一（每线程生成器计数器都从 0 开始，不处理会重复）
//   - 点云与真值事件一一对应
//   root -l -b -q test/testSimParallel.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "io/CloudIO.h"
#include "io/TreeIO.h"

#include "TSystem.h"

#include <cstdio>
#include <set>
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

void testSimParallel() {
    std::printf("========== 多线程模拟测试 ==========\n");

    const std::string cfg = "/tmp/tsp.cfg";
    const std::string clouds = "/tmp/tsp_clouds.root";
    const std::string kin = "/tmp/tsp_clouds_kin.root";

    FILE* f = std::fopen(cfg.c_str(), "w");
    std::fprintf(f,
                 "[run]\nseed = 21\nn_events = 40\nthreads = 4\n\n"
                 "[gas]\nA = 1\nZ = 1\nstoich = 2\npressure_torr = 600\n\n"
                 "[reaction]\nbeam_A = 12\nbeam_Z = 4\ntarget_A = 1\ntarget_Z = 1\n"
                 "ejectile_A = 1\nejectile_Z = 1\n\n"
                 "[beam]\nke_mev = 40\ntheta_min = 50\ntheta_max = 110\n\n"
                 "[vertex]\nz_min = 0.05\nz_max = 0.45\nrho_sigma = 0.005\n\n"
                 "[detector]\nmode = ejectile\nbfield_t = 2.85\nefield_vm = 60000\n");
    std::fclose(f);

    gSystem->Unlink(clouds.c_str());
    gSystem->Unlink(kin.c_str());

    const int rc = gSystem->Exec(
        ("build/rk4_simulate " + cfg + " " + clouds + " 40 > /tmp/tsp_sim.log 2>&1").c_str());
    Check(rc == 0, "模拟成功");
    Check(!gSystem->AccessPathName(clouds.c_str()), "点云文件已生成");
    Check(!gSystem->AccessPathName(kin.c_str()), "真值文件已生成");

    std::vector<PointCloud> cloud_list;
    std::vector<Event> events;
    Check(ReadPointClouds(clouds, cloud_list), "读取点云");
    Check(ReadKinematics(kin, events), "读取真值");

    std::printf("    点云 %zu 条，真值 %zu 条\n", cloud_list.size(), events.size());
    Check(!cloud_list.empty(), "有输出");
    Check(cloud_list.size() == events.size(), "点云与真值数量一致");

    // 事件编号唯一（核心：多线程下不能重复）
    std::set<long long> cloud_ids;
    for (const auto& c : cloud_list) cloud_ids.insert(c.event);
    std::set<long long> event_ids;
    for (const auto& e : events) event_ids.insert(e.id);

    Check(cloud_ids.size() == cloud_list.size(), "点云事件编号唯一");
    Check(event_ids.size() == events.size(), "真值事件编号唯一");

    // 一一对应
    bool matched = true;
    for (const auto& c : cloud_list) {
        if (event_ids.find(c.event) == event_ids.end()) matched = false;
    }
    Check(matched, "每条点云都能找到对应真值");

    const std::string single_clouds="/tmp/tsp_single.root", single_kin="/tmp/tsp_single_kin.root";
    const int single_rc=gSystem->Exec(("SPYRALPP_THREADS=1 build/rk4_simulate "+cfg+" "+single_clouds+" 40 > /tmp/tsp_single.log 2>&1").c_str());
    Check(single_rc==0,"Single-thread simulation succeeded");
    std::vector<PointCloud> baseline_clouds;std::vector<Event> baseline_events;
    Check(ReadPointClouds(single_clouds,baseline_clouds),"Read single-thread clouds");
    Check(ReadKinematics(single_kin,baseline_events),"Read single-thread truth");
    bool exact=baseline_clouds.size()==cloud_list.size() && baseline_events.size()==events.size();
    for(size_t i=0;exact && i<cloud_list.size();++i){
        const auto&a=cloud_list[i];const auto&b=baseline_clouds[i];
        exact=a.event==b.event && a.points.size()==b.points.size();
        for(size_t j=0;exact && j<a.points.size();++j){
            const auto&x=a.points[j];const auto&y=b.points[j];
            exact=x.pad==y.pad && x.tb==y.tb && x.x==y.x && x.y==y.y && x.z==y.z &&
                  x.amplitude==y.amplitude && x.integral==y.integral && x.scale==y.scale;
        }
        const auto&e=events[i];const auto&f=baseline_events[i];
        exact=exact && e.id==f.id && e.vx==f.vx && e.vy==f.vy && e.vz==f.vz &&
              e.beam_energy==f.beam_energy && e.particles.size()==f.particles.size();
        for(size_t j=0;exact && j<e.particles.size();++j){
            const auto&x=e.particles[j];const auto&y=f.particles[j];
            exact=x.A==y.A && x.Z==y.Z && x.label==y.label && x.px==y.px && x.py==y.py && x.pz==y.pz && x.E==y.E;
        }
    }
    Check(exact,"Single/four-thread clouds and particle truth agree exactly");
    gSystem->Unlink(single_clouds.c_str());gSystem->Unlink(single_kin.c_str());

    gSystem->Unlink(cfg.c_str());
    gSystem->Unlink(clouds.c_str());
    gSystem->Unlink(kin.c_str());

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
