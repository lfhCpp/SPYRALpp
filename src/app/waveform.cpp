//----------------------------------------------------------------------
// rk4_waveform：原始波形 -> 点云
//
// 用法：rk4_waveform <config> <waveform> <out_clouds.root> [max_events] [threads]
//
//   <waveform>  ：.h5/.hdf5（GET/merger 的真实数据）或文本波形
//   [max_events]：只处理前 N 个事件（默认全部）
//   [threads]   ：处理线程数；不写则取 run.threads，再退回 SPYRALPP_THREADS
//
// 线程模型：主线程顺序读 HDF5（HDF5 默认不是线程安全的），worker 线程做
// 基线扣除 + 寻峰 + 点云；结果按事件索引排序后写出，所以**多线程与单线程
// 的输出逐位一致**，只是快。
//----------------------------------------------------------------------

#include "core/DataPath.h"
#include "io/CloudIO.h"
#include "io/WorkspaceIO.h"
#include "io/StageIO.h"
#include "io/Config.h"
#include "reco/H5TraceReader.h"
#include "reco/PeakFinder.h"
#include "reco/PointCloudBuilder.h"
#include "reco/TraceReader.h"
#include "sim/DetectorConfig.h"
#include "sim/PadPlane.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <exception>
#include <stdexcept>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace rk4;

namespace {

struct Result {
    int index = -1;
    PointCloud cloud;
    std::vector<TracePeakInfo> peaks;
};

struct Item {
    int index = -1;
    TraceEvent ev;
};

int ResolveThreads(const Config& cfg, int argc, char** argv, int max_events) {
    int threads = cfg.GetInt("run.threads", 0);
    if (const char* env = std::getenv("SPYRALPP_THREADS")) {
        const int v = std::atoi(env);
        if (v > 0) threads = v;
    }
    if (argc > 5) {
        const int v = std::atoi(argv[5]);
        if (v > 0) threads = v;
    }
    if (threads <= 0) {
        threads = static_cast<int>(std::thread::hardware_concurrency());
        if (threads <= 0) threads = 1;
    }
    if (max_events > 0 && threads > max_events) threads = max_events;
    return std::max(1, threads);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::printf(
            "用法: rk4_waveform <config> <waveform(.h5|.txt)> <out_clouds.root> "
            "[max_events] [threads]\n");
        return 1;
    }

    Config cfg;
    if (!cfg.Load(argv[1])) {
        std::fprintf(stderr, "无法读取配置 %s\n", argv[1]);
        return 1;
    }

    // 按扩展名选择读取器：.h5/.hdf5 用 HDF5（真实 GET 数据），其余按文本格式
    const std::string wave_path = argv[2];
    const auto suffix = [&](const std::string& ending) {
        return wave_path.size() >= ending.size() &&
               wave_path.compare(wave_path.size()-ending.size(), ending.size(), ending) == 0;
    };
    const bool is_h5 = suffix(".h5") || suffix(".hdf") || suffix(".hdf5");

    TextTraceReader text_reader;
    H5TraceReader h5_reader;
    TraceReader* reader = nullptr;

    if (is_h5) {
        if (!H5TraceReader::Available()) {
            std::fprintf(stderr, "该构建未启用 HDF5，无法读取 %s\n", wave_path.c_str());
            return 1;
        }
        if (!h5_reader.Open(wave_path)) {
            std::fprintf(stderr, "无法打开 HDF5 波形 %s\n", wave_path.c_str());
            return 1;
        }
        reader = &h5_reader;
    } else {
        if (!text_reader.Open(wave_path)) {
            std::fprintf(stderr, "无法打开波形 %s\n", wave_path.c_str());
            return 1;
        }
        reader = &text_reader;
    }

    PadPlane pads;
    const std::string pad_dir = cfg.GetString("detector.pad_dir", rk4::DataDir());
    if (!pads.Load(pad_dir)) return 1;
    (void)pads.LoadTimeOffsets(pad_dir); // 没有时间修正文件也能继续

    DetectorConfig det;
    det.length_m = cfg.GetDouble("detector.length_m", det.length_m);
    det.micromegas_edge = cfg.GetDouble("detector.micromegas_edge", det.micromegas_edge);
    det.windows_edge = cfg.GetDouble("detector.windows_edge", det.windows_edge);
    det.adc_threshold = cfg.GetDouble("detector.adc_threshold", det.adc_threshold);

    PeakFinderParameters peaks;
    peaks.peak_separation = cfg.GetDouble("peaks.separation", peaks.peak_separation);
    peaks.peak_prominence = cfg.GetDouble("peaks.prominence", peaks.peak_prominence);
    peaks.peak_max_width = cfg.GetDouble("peaks.max_width", peaks.peak_max_width);
    peaks.peak_threshold = cfg.GetDouble("peaks.threshold", peaks.peak_threshold);

    const double baseline_scale = cfg.GetDouble("peaks.baseline_scale", 20.0);
    FribConfig frib;
    frib.baseline_scale=cfg.GetDouble("frib.baseline_scale",frib.baseline_scale);
    frib.delay=cfg.GetDouble("frib.ic_delay",frib.delay);
    frib.max_multiplicity=cfg.GetInt("frib.ic_multiplicity",frib.max_multiplicity);
    frib.peaks.peak_separation=cfg.GetDouble("frib.separation",frib.peaks.peak_separation);
    frib.peaks.peak_prominence=cfg.GetDouble("frib.prominence",frib.peaks.peak_prominence);
    frib.peaks.peak_max_width=cfg.GetDouble("frib.max_width",frib.peaks.peak_max_width);
    frib.peaks.peak_threshold=cfg.GetDouble("frib.threshold",frib.peaks.peak_threshold);

    // GET 电子学响应模板拟合（含饱和修复），默认关
    const bool template_fit = cfg.GetBool("peaks.template_fit", false);
    // 可选：写每条迹线的峰参数（<out>_peaks.root），方便查饱和修复/峰分裂
    const bool dump_peaks = cfg.GetBool("peaks.dump", false);
    std::vector<TracePeakRecord> peak_records;

    PulseShapeConfig resp_cfg;
    resp_cfg.n_samples = det.num_tb;
    resp_cfg.sampling_ns = cfg.GetDouble("peaks.sampling_ns", resp_cfg.sampling_ns);
    resp_cfg.peaking_ns = cfg.GetDouble("peaks.peaking_ns", resp_cfg.peaking_ns);
    resp_cfg.sat_adc = cfg.GetDouble("peaks.sat_adc", resp_cfg.sat_adc);
    TemplateFitConfig tfit_cfg;
    tfit_cfg.sat_adc = resp_cfg.sat_adc;
    tfit_cfg.max_peaks = cfg.GetInt("peaks.max_peaks", tfit_cfg.max_peaks);
    tfit_cfg.residual_threshold =
        cfg.GetDouble("peaks.residual_threshold", tfit_cfg.residual_threshold);

    const int max_events = (argc > 4) ? std::atoi(argv[4]) : reader->NumEvents();
    const int n_threads = ResolveThreads(cfg, argc, argv, max_events);
    const int total_event_limit = std::max(0, std::min(max_events, reader->NumEvents()));
    const int progress_interval = std::max(1, total_event_limit / 50 + (total_event_limit % 50 != 0));

    const auto t0 = std::chrono::steady_clock::now();

    // ---- 生产者 / 消费者 ----
    std::atomic<int> sat_traces{0};
    std::mutex qmtx;
    std::mutex rmtx;
    std::mutex progress_mtx;
    int completed_events = 0;
    std::atomic<uint64_t> worker_idle_ns{0};
    std::atomic<uint64_t> build_ns{0};
    std::condition_variable qcv;
    std::condition_variable qfull;
    std::deque<Item> queue;
    bool reading_done = false;
    std::atomic<bool> cancelled{false};
    std::string worker_error;
    std::vector<Result> results;
    std::vector<bool> excluded;

    // 队列上限：控制内存（一个事件约 386 条迹 × 512 float ≈ 0.8 MB）
    const size_t max_inflight = static_cast<size_t>(std::max(2, 2 * n_threads));

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(n_threads));
    for (int t = 0; t < n_threads; ++t) {
        pool.emplace_back([&]() {
            try {
            PointCloudBuilder builder; // 每线程一套 FFTW plan / 缓冲
            if (!builder.Init(det, &pads, peaks, baseline_scale))
                throw std::runtime_error("Waveform worker initialization failed");
            builder.ConfigureFrib(frib);
            builder.SetDeferredDither(true);
            if (template_fit) builder.ConfigureTemplateFit(true, resp_cfg, tfit_cfg);
            std::vector<TracePeakInfo> local_peaks;
            if (dump_peaks) builder.SetTracePeakSink(&local_peaks);
            while (true) {
                Item item;
                {
                    std::unique_lock<std::mutex> lk(qmtx);
                    if (queue.empty() && !cancelled.load() && !reading_done) {
                        const auto idle_start = std::chrono::steady_clock::now();
                        qcv.wait(lk, [&] { return cancelled.load() || reading_done || !queue.empty(); });
                        worker_idle_ns.fetch_add(static_cast<uint64_t>(
                            std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - idle_start).count()),
                            std::memory_order_relaxed);
                    }
                    if (cancelled.load()) return;
                    if (queue.empty()) {
                        sat_traces += builder.SaturatedTraces(); // 汇总饱和迹线
                        return; // reading_done 且队列空
                    }
                    item = std::move(queue.front());
                    queue.pop_front();
                    qfull.notify_one();
                }
                PointCloud cloud;
                const auto build_start = std::chrono::steady_clock::now();
                const bool ok = builder.Build(item.ev, cloud);
                build_ns.fetch_add(static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - build_start).count()),
                    std::memory_order_relaxed);
                if (!ok) throw std::runtime_error("Cannot build event " + std::to_string(item.ev.number));
                if (ok) {
                    {
                        std::lock_guard<std::mutex> lk(rmtx);
                        results.push_back(Result{item.index, std::move(cloud), std::move(local_peaks)});
                        local_peaks.clear();
                    }
                    {
                        std::lock_guard<std::mutex> lk(progress_mtx);
                        ++completed_events;
                        if (completed_events % progress_interval == 0 || completed_events == total_event_limit) {
                            const int percent = total_event_limit > 0
                                ? completed_events * 100 / total_event_limit : 100;
                            constexpr int bar_width = 24;
                            const int filled = total_event_limit > 0
                                ? completed_events * bar_width / total_event_limit : bar_width;
                            const std::string bar(static_cast<size_t>(filled), '#');
                            const std::string full_bar = bar + std::string(bar_width - filled, '.');
                            const int digits = static_cast<int>(std::to_string(total_event_limit).size());
                            std::printf("\rPointcloud [%s] %3d%% | %*d/%d events%s",
                                        full_bar.c_str(), percent, digits, completed_events,
                                        total_event_limit,
                                        completed_events == total_event_limit ? "\n" : "");
                            std::fflush(stdout);
                        }
                    }
                }
            }
            } catch (const std::exception& e) {
                {
                    std::lock_guard<std::mutex> lk(qmtx);
                    if (worker_error.empty()) worker_error = e.what();
                    cancelled.store(true);
                }
                qcv.notify_all(); qfull.notify_all();
            }
        });
    }

    std::printf("Pointcloud: processing %d events with %d threads; peak method: %s\n",
                total_event_limit, n_threads,
                template_fit ? "GET response-template fit (with saturation recovery)" : "baseline peak finder");
    std::fflush(stdout);
    int read = 0;
    uint64_t reader_ns = 0;
    uint64_t queue_wait_ns = 0;
    long long first_event=0,last_event=-1;
    TraceEvent tev;
    std::string read_error;
    try {
    while (!cancelled.load() && read < max_events) {
        const auto read_start = std::chrono::steady_clock::now();
        const bool has_event = reader->ReadNext(tev);
        reader_ns += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - read_start).count());
        if (!has_event) break;
        if(read==0)first_event=tev.number;
        last_event=tev.number;
        tev.metadata.orig_run=cfg.GetInt("run.number",-1);
        tev.metadata.orig_event=tev.number;
        excluded.push_back(tev.ic_downscale);
        {
            std::unique_lock<std::mutex> lk(qmtx);
            const auto queue_wait_start = std::chrono::steady_clock::now();
            qfull.wait(lk, [&] { return cancelled.load() || queue.size() < max_inflight; });
            queue_wait_ns += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - queue_wait_start).count());
            if (cancelled.load()) break;
            queue.push_back(Item{read, std::move(tev)});
            qcv.notify_one();
        }
        tev = TraceEvent(); // move 之后重建，下一轮 ReadNext 会重新填
        ++read;
    }
    } catch (const std::exception& e) { read_error = e.what(); }
    {
        std::lock_guard<std::mutex> lk(qmtx);
        reading_done = true;
    }
    qcv.notify_all();
    for (auto& th : pool) th.join();
    if (!worker_error.empty()) read_error = worker_error;
    if (!read_error.empty()) {
        std::fprintf(stderr, "Read failed: %s\n", read_error.c_str());
        return 1;
    }
    if (results.size() != static_cast<size_t>(read)) {
        std::fprintf(stderr, "Waveform workers did not finish all events\n");
        return 1;
    }

    const auto finalize_start = std::chrono::steady_clock::now();
    // 按事件索引排序 → 多线程与单线程输出一致
    std::sort(results.begin(), results.end(),
              [](const Result& a, const Result& b) {
                  return a.index < b.index;
              });

    std::vector<PointCloud> clouds;
    clouds.reserve(results.size());
    size_t total_points = 0;
    NumpyPcg64 dither(static_cast<uint64_t>(cfg.GetInt("peaks.dither_seed", 20261026)));
    int downscale_count = 0;
    for (auto& r : results) {
        PointCloudBuilder::FinalizeDither(r.cloud, dither, det, &r.peaks);
        if (excluded.at(r.index)) { ++downscale_count; continue; }
        for (const auto& info : r.peaks) {
            TracePeakRecord peak;
            peak.event = r.cloud.event; peak.pad = info.pad; peak.index = info.index;
            peak.npeaks = info.npeaks; peak.time = info.time; peak.amplitude = info.amplitude;
            peak.charge = info.charge; peak.chi2 = info.chi2; peak.npoints = info.npoints;
            peak.saturated = info.saturated;
            peak_records.push_back(peak);
        }
        total_points += r.cloud.points.size();
        clouds.push_back(std::move(r.cloud));
    }
    std::printf("FRIB IC downscale events excluded: %d\n", downscale_count);

    const auto finalize_end = std::chrono::steady_clock::now();
    const double finalize_seconds = std::chrono::duration<double>(finalize_end - finalize_start).count();

    const auto write_start = std::chrono::steady_clock::now();
    try {
        PointcloudData data;
        data.info.run=cfg.GetInt("run.number",-1);
        data.info.first_event=first_event;data.info.last_event=last_event;
        data.info.raw_events=read;data.info.excluded_ic=downscale_count;
        data.info.seed=static_cast<uint64_t>(cfg.GetInt("peaks.dither_seed",20261026));
        data.info.waveform_method=template_fit?"template":"baseline";
        data.events=std::move(clouds);
        WritePointcloudStage(argv[3],data);
        clouds=std::move(data.events);
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Pointcloud write failed: %s\n",e.what());return 1;
    }
    const auto write_end = std::chrono::steady_clock::now();
    const double write_seconds = std::chrono::duration<double>(write_end - write_start).count();
    const double seconds = std::chrono::duration<double>(write_end - t0).count();

    std::printf("波形转换完成：读取 %d 个事件，输出 %zu 个点云\n", read, clouds.size());
    std::printf("  线程 %d，耗时 %.2f s（%.2f 事件/s），总点数 %zu（平均 %.0f 点/点云）\n",
                n_threads, seconds, seconds > 0.0 ? read / seconds : 0.0, total_points,
                clouds.empty() ? 0.0
                               : static_cast<double>(total_points) / clouds.size());
    std::printf("  性能拆分：HDF5读取 %.2f s，队列满等待 %.2f s，工作线程空闲合计 %.2f 线程·s，事件构建合计 %.2f 线程·s，整理 %.2f s，写出 %.2f s\n",
                reader_ns / 1e9, queue_wait_ns / 1e9, worker_idle_ns.load(std::memory_order_relaxed) / 1e9,
                build_ns.load(std::memory_order_relaxed) / 1e9, finalize_seconds, write_seconds);
    if (template_fit) {
        std::printf("  饱和迹线：%d（用两侧样本/平顶宽度修复幅度）\n", sat_traces.load());
    }
    if (dump_peaks) {
        std::sort(peak_records.begin(), peak_records.end(),
                  [](const TracePeakRecord& a, const TracePeakRecord& b) {
                      if (a.event != b.event) return a.event < b.event;
                      if (a.pad != b.pad) return a.pad < b.pad;
                      return a.index < b.index;
                  });
        const std::string peaks_path = std::string(argv[3]) + "_peaks.root";
        if (WriteTracePeaks(peaks_path, peak_records)) {
            std::printf("  峰参数文件：%s（%zu 个峰）\n", peaks_path.c_str(),
                        peak_records.size());
        }
    }
    std::printf("  输出：%s\n", argv[3]);
    return 0;
}
