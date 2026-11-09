//----------------------------------------------------------------------
// rk4_simulate：配置 -> 事件 -> 探测器数字化 -> 点云 + 运动学真值
//
// 用法：rk4_simulate <config> [out_clouds.root] [n_events]
//----------------------------------------------------------------------

#include "core/Event.h"
#include "core/DataPath.h"
#include "io/CloudIO.h"
#include "io/H5Writer.h"
#include "io/Config.h"
#include "io/TreeIO.h"
#include "sim/DetectorConfig.h"
#include "sim/KinematicsGenerator.h"
#include "sim/Simulator.h"
#include "sim/TargetDensity.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using namespace rk4;

// 把 "xxx.root" 换成 "xxx{suffix}.h5"
static std::string DeriveH5Path(const std::string& path, const std::string& suffix) {
    std::string base = path;
    if (base.size() > 5 && base.compare(base.size() - 5, 5, ".root") == 0) {
        base.resize(base.size() - 5);
    }
    return base + suffix + ".h5";
}

static std::string DerivePath(const std::string& path, const std::string& suffix) {
    std::string base = path;
    if (base.size() > 5 && base.compare(base.size() - 5, 5, ".root") == 0) {
        base.resize(base.size() - 5);
    }
    return base + suffix + ".root";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("用法: rk4_simulate <config> [out_clouds.root] [n_events]\n");
        return 1;
    }

    Config cfg;
    if (!cfg.Load(argv[1])) {
        std::fprintf(stderr, "无法读取配置 %s\n", argv[1]);
        return 1;
    }

    const int n_events = (argc > 3) ? std::atoi(argv[3])
        : cfg.GetInt("simulation.events", cfg.GetInt("run.n_events", 100));
    if (n_events <= 0) {
        std::fprintf(stderr, "模拟事件数必须大于0\n");
        return 1;
    }

    std::string cloud_path = argc > 2 ? argv[2] : cfg.GetString("simulation.output", "");
    if (cloud_path.empty()) {
        const std::string workspace = cfg.GetString("output.workspace", "workspace");
        cloud_path = (std::filesystem::path(workspace) / "Simulation" / "simulated.root").string();
    }
    if (cloud_path.size() < 5 || cloud_path.compare(cloud_path.size() - 5, 5, ".root") != 0) {
        std::fprintf(stderr, "模拟点云输出必须使用 .root 后缀：%s\n", cloud_path.c_str());
        return 1;
    }
    const auto output_parent = std::filesystem::path(cloud_path).parent_path();
    if (!output_parent.empty()) std::filesystem::create_directories(output_parent);

    std::vector<Element> comp;
    const int component_count = cfg.GetInt("gas.components", 0);
    if (component_count > 0) {
        for (int i = 0; i < component_count; ++i) {
            const std::string key = "gas.component" + std::to_string(i);
            Element element;
            element.a = cfg.GetInt(key + ".a", 0);
            element.z = cfg.GetInt(key + ".z", 0);
            element.s = cfg.GetDouble(key + ".s", 0.0);
            comp.push_back(element);
        }
    } else {
        Element element;
        element.a = cfg.GetInt("gas.a", 1);
        element.z = cfg.GetInt("gas.z", 1);
        element.s = cfg.GetDouble("gas.stoich", 2);
        comp.push_back(element);
    }
    GasTarget gas(comp, cfg.GetDouble("gas.pressure_torr", 600.0));

    GeneratorConfig g;
    g.reaction.beam_A = cfg.GetInt("reaction.beam_a", 12);
    g.reaction.beam_Z = cfg.GetInt("reaction.beam_z", 4);
    g.reaction.target_A = cfg.GetInt("reaction.target_a", 1);
    g.reaction.target_Z = cfg.GetInt("reaction.target_z", 1);
    g.reaction.ejectile_A = cfg.GetInt("reaction.ejectile_a", 1);
    g.reaction.ejectile_Z = cfg.GetInt("reaction.ejectile_z", 1);
    g.beam_ke_mev = cfg.GetDouble("beam.ke_mev", 150.0);
    g.ex_centroid = cfg.GetDouble("beam.ex_centroid", 0.0);
    g.ex_fwhm = cfg.GetDouble("beam.ex_fwhm", 0.0);
    g.theta_min_deg = cfg.GetDouble("beam.theta_min", 0.0);
    g.theta_max_deg = cfg.GetDouble("beam.theta_max", 180.0);
    g.z_min_m = cfg.GetDouble("vertex.z_min", 0.0);
    g.z_max_m = cfg.GetDouble("vertex.z_max", 1.0);
    g.rho_sigma_m = cfg.GetDouble("vertex.rho_sigma", 0.007);

    KinematicsGenerator generator;
    const uint64_t seed = static_cast<uint64_t>(
        cfg.GetInt("simulation.seed", cfg.GetInt("run.seed", 1)));
    if (!generator.Init(g, &gas, seed)) {
        return 1;
    }

    DetectorConfig det;
    det.length_m = cfg.GetDouble("detector.length_m", det.length_m);
    det.bfield_t = cfg.GetDouble("detector.bfield_t", det.bfield_t);
    det.efield_vm = cfg.GetDouble("detector.efield_vm", det.efield_vm);
    det.w_value_ev = cfg.GetDouble("detector.w_value_ev", det.w_value_ev);
    det.fano = cfg.GetDouble("detector.fano", det.fano);
    det.diffusion = cfg.GetDouble("detector.diffusion", det.diffusion);
    det.mpgd_gain = cfg.GetInt("detector.mpgd_gain", det.mpgd_gain);
    det.adc_threshold = cfg.GetDouble("detector.adc_threshold", det.adc_threshold);
    det.dither_time_buckets = cfg.GetBool("detector.dither_time_buckets", det.dither_time_buckets);
    det.pad_dir = cfg.GetString("detector.pad_dir", rk4::DataDir());

    DetectorSimulator simulator;
    if (!simulator.Init(det, &gas)) return 1;

    std::vector<int> indices;
    const std::string mode = cfg.GetString(
        "simulation.mode", cfg.GetString("detector.mode", "ejectile"));
    if (mode == "all") {
        indices = {2, 3};
    } else {
        indices = {2};
    }

    // ---- 多线程生成 ----
    // 默认 1 线程（保证输出确定、可复现。用户显式设 threads > 1 时并行（每线程独立随机子流，此时事件集合依赖线程数。
    unsigned int n_threads = static_cast<unsigned int>(cfg.GetInt("run.threads", 1));
    // 环境变量覆盖（脚本里用 SPYRALPP_THREADS=N 即可）
    if (const char* env_threads = std::getenv("SPYRALPP_THREADS")) {
        const int v = std::atoi(env_threads);
        if (v > 0) n_threads = static_cast<unsigned int>(v);
    }
    if (n_threads == 0) {
        n_threads = std::thread::hardware_concurrency();
        if (n_threads == 0) n_threads = 1;
    }
    if (n_threads > static_cast<unsigned int>(n_events) && n_events > 0) {
        n_threads = static_cast<unsigned int>(n_events);
    }
    if (n_threads == 0) n_threads = 1;

    std::printf("并行生成：%u 线程，%d 个事件\n", n_threads, n_events);

    // 每个线程独立的生成器与模拟器：两者都持有随机数状态与工作缓冲，
    // 不能共享。实际抽样用 (run seed, event index, component) 子流，与调度无关。
    struct ThreadResult {
        std::vector<Event> events;
        std::vector<PointCloud> clouds;
    };
    std::vector<ThreadResult> results(n_threads);
    std::atomic<int> next_event{0};
    std::atomic<int> completed_events{0};
    std::mutex progress_mutex;
    int last_displayed = 0;
    const int progress_stride = std::max(1, n_events / 100);
    auto report_progress = [&]() {
        const int done = completed_events.fetch_add(1) + 1;
        if (done != n_events && done % progress_stride != 0) return;
        std::lock_guard<std::mutex> lock(progress_mutex);
        if (done <= last_displayed) return;
        last_displayed = done;
        constexpr int width = 30;
        const int filled = width * done / n_events;
        const std::string bar(static_cast<size_t>(filled), '#');
        const std::string tail(static_cast<size_t>(width - filled), '-');
        std::printf("\r模拟进度 [%s%s] %d/%d (%d%%)", bar.c_str(), tail.c_str(), done,
                    n_events, 100 * done / n_events);
        if (done == n_events) std::printf("\n");
        std::fflush(stdout);
    };

    auto worker = [&](unsigned int tid) {
        KinematicsGenerator local_gen;
        if (!local_gen.Init(g, &gas, seed + tid)) {
            return;
        }

        DetectorSimulator local_sim;
        if (!local_sim.Init(det, &gas)) return;

        while (true) {
            const int idx = next_event.fetch_add(1);
            if (idx >= n_events) break;

            Event ev;
            local_gen.SetEventSeed(seed, static_cast<uint64_t>(idx));
            if (!local_gen.Generate(ev)) {
                report_progress();
                continue;
            }

            // 每个线程的生成器有自己的 event_id_ 计数器（都从 0 开始。直接用会产生重复 ID。
            // 统一使用第 idx 个请求作为全局事件号，并使用对应的事件随机数子流。
            ev.id = idx;

            PointCloud cloud;
            local_sim.SetEventSeed(seed, static_cast<uint64_t>(idx));
            local_sim.Simulate(ev, cloud, indices);
            if (cloud.points.empty()) {
                report_progress();
                continue;
            }

            results[tid].events.push_back(ev);
            results[tid].clouds.push_back(cloud);
            report_progress();
        }
    };

    if (n_threads <= 1) {
        worker(0);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(n_threads);
        for (unsigned int t = 0; t < n_threads; ++t) pool.emplace_back(worker, t);
        for (auto& th : pool) th.join();
    }

    // 合并并按事件编号排序，保证输出顺序确定
    std::vector<Event> events;
    std::vector<PointCloud> clouds;
    for (auto& r : results) {
        events.insert(events.end(), r.events.begin(), r.events.end());
        clouds.insert(clouds.end(), r.clouds.begin(), r.clouds.end());
    }

    std::vector<size_t> order(clouds.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return clouds[a].event < clouds[b].event;
    });

    std::vector<Event> sorted_events;
    std::vector<PointCloud> sorted_clouds;
    sorted_events.reserve(order.size());
    sorted_clouds.reserve(order.size());
    for (size_t idx : order) {
        sorted_events.push_back(events[idx]);
        sorted_clouds.push_back(clouds[idx]);
    }
    events.swap(sorted_events);
    clouds.swap(sorted_clouds);

    const std::string kin_path = DerivePath(cloud_path, "_kin");
    if (!WritePointClouds(cloud_path, clouds)) return 1;
    if (!WriteKinematics(kin_path, events)) return 1;

    std::printf("模拟完成：%zu 个事件\n", clouds.size());
    std::printf("  点云（TTree）：%s\n", cloud_path.c_str());
    std::printf("  真值（TTree）：%s\n", kin_path.c_str());

    // 可选：同时写 HDF5（结构与 attpc_engine 的 SpyralWriter 一致，Python 版 Spyral 可直接读）
    if (cfg.GetBool("simulation.hdf5", cfg.GetBool("output.hdf5", false))) {
        if (!Hdf5Available()) {
            std::fprintf(stderr, "该构建未启用 HDF5，跳过 HDF5 输出\n");
        } else {
            const std::string h5_cloud = DeriveH5Path(cloud_path, "_cloud");
            const std::string h5_kin = DeriveH5Path(cloud_path, "_kin");
            if (WritePointCloudsHdf5(h5_cloud, clouds)) {
                std::printf("  点云（HDF5）  ：%s\n", h5_cloud.c_str());
            } else {
                std::fprintf(stderr, "点云 HDF5 写出失败\n");
            }
            if (WriteKinematicsHdf5(h5_kin, events)) {
                std::printf("  真值（HDF5）  ：%s\n", h5_kin.c_str());
            } else {
                std::fprintf(stderr, "真值 HDF5 写出失败\n");
            }
        }
    }
    return 0;
}
