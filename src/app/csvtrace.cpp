//----------------------------------------------------------------------
// rk4_csvtrace：CSV 波形 -> 点云
//
// 用途：真实 GET/FRIB 数据在 HDF5 里，本项目需要 libhdf5-dev 才能直接读。
// 在没有该依赖的环境里，可以把波形导出成 CSV 后用它接入，避免整条链断掉。
//
// 每行格式：
//   <event>,<pad>,<sample0>,<sample1>,...
//
// 用法：rk4_csvtrace <config> <traces.csv> <out_clouds.root> [max_events]
//----------------------------------------------------------------------

#include "core/Event.h"
#include "core/DataPath.h"
#include "io/CloudIO.h"
#include "io/Config.h"
#include "reco/PeakFinder.h"
#include "reco/PointCloudBuilder.h"
#include "sim/DetectorConfig.h"
#include "sim/PadPlane.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace rk4;

namespace {

struct TraceRow {
    long long event = 0;
    int pad = -1;
    std::vector<float> samples;
};

bool ParseLine(const std::string& line, TraceRow& out) {
    std::stringstream ss(line);
    std::string field;

    if (!std::getline(ss, field, ',')) return false;
    out.event = std::atoll(field.c_str());

    if (!std::getline(ss, field, ',')) return false;
    out.pad = std::atoi(field.c_str());

    out.samples.clear();
    while (std::getline(ss, field, ',')) {
        out.samples.push_back(static_cast<float>(std::atof(field.c_str())));
    }
    return !out.samples.empty();
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::printf("用法: rk4_csvtrace <config> <traces.csv> <out_clouds.root> [max_events]\n");
        return 1;
    }

    Config cfg;
    if (!cfg.Load(argv[1])) {
        std::fprintf(stderr, "无法读取配置 %s\n", argv[1]);
        return 1;
    }

    std::ifstream in(argv[2]);
    if (!in.is_open()) {
        std::fprintf(stderr, "无法打开 %s\n", argv[2]);
        return 1;
    }

    PadPlane pads;
    const std::string pad_dir = cfg.GetString("detector.pad_dir", rk4::DataDir());
    if (!pads.Load(pad_dir)) {
        std::fprintf(stderr, "无法加载 pad 表\n");
        return 1;
    }
    (void)pads.LoadTimeOffsets(pad_dir);

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

    PointCloudBuilder builder;
    if (!builder.Init(det, &pads, peaks, cfg.GetDouble("peaks.baseline_scale", 20.0))) {
        std::fprintf(stderr, "点云构建器初始化失败\n");
        return 1;
    }

    const int max_events = (argc > 4) ? std::atoi(argv[4]) : 0; // 0 = 全部

    // 按事件分组
    std::vector<PointCloud> clouds;
    std::string line;
    long long current_event = -1;
    std::vector<rk4::CloudPoint> dummy; // 占位（未使用）

    std::vector<TraceRow> batch; // 当前事件的波形
    int events_done = 0;

    auto flush_event = [&]() {
        if (batch.empty()) return;
        TraceEvent tev;
        tev.number = static_cast<int>(batch.front().event);
        for (const auto& row : batch) {
            Trace t;
            t.pad = row.pad;
            t.samples = row.samples;
            tev.traces.push_back(t);
        }
        PointCloud pc;
        if (builder.Build(tev, pc) && !pc.points.empty()) {
            clouds.push_back(pc);
        }
        batch.clear();
    };

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        TraceRow row;
        if (!ParseLine(line, row)) continue;

        if (row.event != current_event && !batch.empty()) {
            flush_event();
            ++events_done;
            if (max_events > 0 && events_done >= max_events) break;
        }
        current_event = row.event;
        batch.push_back(row);
    }
    flush_event();

    if (!WritePointClouds(argv[3], clouds)) {
        std::fprintf(stderr, "写出失败\n");
        return 1;
    }

    std::printf("CSV 波形转换完成：%zu 个事件的点云\n", clouds.size());
    std::printf("  输出：%s\n", argv[3]);
    return 0;
}
