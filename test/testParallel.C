//----------------------------------------------------------------------
// 并行与共享映射测试：
//   1) 同一网格 Load / LoadShared 插值结果一致
//   2) 多线程并行拟合与单线程结果一致（逐位）
//   root -l -b -q test/testParallel.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "core/Event.h"
#include "io/ResultIO.h"
#include "core/MassData.h"
#include "core/Units.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/TrackFitter.h"
#include "reco/TrackMesh.h"
#include "sim/DetectorConfig.h"
#include "sim/Simulator.h"
#include "sim/TargetDensity.h"

#include "TSystem.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <thread>
#include <vector>

using namespace rk4;

static int g_fail = 0;

static void Check(bool ok, const char* msg) {
    if (!ok) {
        std::printf("  [FAIL] %s\n", msg);
        ++g_fail;
    }
}

static Particle MakeParticle(int A, int Z, double ke, double polar) {
    const double mass = MassData::GetNuclearMassNuclear(A, Z) * Units::u_MeV;
    const double p = std::sqrt(ke * (ke + 2.0 * mass));
    const double th = polar * Units::DEG2RAD;

    Particle q;
    q.A = A;
    q.Z = Z;
    q.label = 2;
    q.px = p * std::sin(th) * std::cos(0.5);
    q.py = p * std::sin(th) * std::sin(0.5);
    q.pz = p * std::cos(th);
    q.E = ke + mass;
    return q;
}

void testParallel() {
    std::printf("========== 并行与共享映射测试 ==========\n");

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

    TrackMesh mesh;
    Check(mesh.Build(mc, 1, 1, gas, 2.85, 60000.0), "Build mesh");

    const std::string path = "/tmp/rk4_parallel_mesh.bin";
    Check(mesh.Save(path), "Save");

    // ---- 1) 普通 Load 与 LoadShared 等价 ----
    TrackMesh loaded;
    Check(loaded.Load(path), "Load");
    TrackMesh shared;
    Check(shared.LoadShared(path), "LoadShared");
    Check(shared.IsShared(), "IsShared 为真");

    std::vector<double> a, b, c;
    Check(mesh.Interpolate(20.0, 70.0 * Units::DEG2RAD, a), "Interpolate 原始");
    Check(loaded.Interpolate(20.0, 70.0 * Units::DEG2RAD, b), "Interpolate Load");
    Check(shared.Interpolate(20.0, 70.0 * Units::DEG2RAD, c), "Interpolate Shared");

    double diff_load = 0.0, diff_shared = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        diff_load = std::max(diff_load, std::fabs(a[i] - b[i]));
        diff_shared = std::max(diff_shared, std::fabs(a[i] - c[i]));
    }
    Check(diff_load < 1e-9, "Load 与原始一致");
    Check(diff_shared < 1e-9, "LoadShared 与原始一致");
    gSystem->Unlink(path.c_str());

    // ---- 2) 多线程拟合与单线程结果一致 ----
    DetectorConfig cfg;
    DetectorSimulator sim;
    Check(sim.Init(cfg, &gas), "Init simulator");

    const double mass = MassData::GetNuclearMassNuclear(1, 1) * Units::u_MeV;

    // 造 8 个不同能量的质子点云
    std::vector<PointCloud> clouds;
    for (int i = 0; i < 8; ++i) {
        const double ke = 8.0 + 3.0 * i;
        Event ev;
        ev.id = i;
        ev.vx = 0.0;
        ev.vy = 0.0;
        ev.vz = 0.10;
        ev.particles = {MakeParticle(1, 1, ke, 60.0 + 2.0 * i)};

        PointCloud cloud;
        if (sim.Simulate(ev, cloud, std::vector<int>{0}) && !cloud.points.empty()) {
            clouds.push_back(cloud);
        }
    }
    Check(clouds.size() >= 4, "生成足够点云");
    if (clouds.size() < 4) {
        gSystem->Exit(1);
    }

    auto run_fits = [&](unsigned int n_threads) {
        std::vector<std::vector<FitRecord>> out(n_threads);
        std::atomic<size_t> next{0};

        auto worker = [&](unsigned int tid) {
            Clusterer clusterer;
            Estimator estimator;
            TrackFitter fitter;
            fitter.Init(&mesh, mass, 1);

            while (true) {
                const size_t idx = next.fetch_add(1);
                if (idx >= clouds.size()) break;

                std::vector<Cluster> clusters;
                clusterer.Run(clouds[idx], clusters);
                for (const auto& cl : clusters) {
                    EstimateResult est;
                    if (!estimator.Estimate(cl, cfg.bfield_t, est)) continue;
                    FitResult fit;
                    if (!fitter.Fit(cl.points, est, est.azimuthal, fit)) continue;

                    FitRecord rec;
                    rec.event = clouds[idx].event;
                    rec.cluster = cl.label;
                    rec.brho = fit.brho;
                    rec.polar = fit.polar;
                    rec.cost = fit.cost;
                    out[tid].push_back(rec);
                }
            }
        };

        if (n_threads <= 1) {
            worker(0);
        } else {
            std::vector<std::thread> pool;
            for (unsigned int t = 0; t < n_threads; ++t) pool.emplace_back(worker, t);
            for (auto& th : pool) th.join();
        }

        std::vector<FitRecord> merged;
        for (auto& v : out) merged.insert(merged.end(), v.begin(), v.end());
        std::sort(merged.begin(), merged.end(), [](const FitRecord& x, const FitRecord& y) {
            if (x.event != y.event) return x.event < y.event;
            return x.cluster < y.cluster;
        });
        return merged;
    };

    const std::vector<FitRecord> single = run_fits(1);
    const std::vector<FitRecord> multi = run_fits(4);

    std::printf("    单线程 %zu 条，4 线程 %zu 条\n", single.size(), multi.size());
    Check(single.size() == multi.size(), "结果条数一致");
    Check(!single.empty(), "有拟合结果");

    if (single.size() == multi.size()) {
        int mismatch = 0;
        double max_db = 0.0;
        for (size_t i = 0; i < single.size(); ++i) {
            if (single[i].event != multi[i].event || single[i].cluster != multi[i].cluster) {
                ++mismatch;
                continue;
            }
            max_db = std::max(max_db, std::fabs(single[i].brho - multi[i].brho));
        }
        Check(mismatch == 0, "事件/簇编号一致");
        Check(max_db < 1e-9, "Bρ 逐位一致");
        std::printf("    最大 Bρ 差异 = %.3g\n", max_db);
    }

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
