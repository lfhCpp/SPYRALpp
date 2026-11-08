//----------------------------------------------------------------------
// M3 波形 -> 点云测试
//   root -l -b -q test/testWaveform.C
//----------------------------------------------------------------------
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)

#include "reco/PointCloudBuilder.h"
#include "reco/TraceReader.h"
#include "sim/DetectorConfig.h"
#include "sim/PadPlane.h"

#include "TSystem.h"

#include <algorithm>
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

void testWaveform() {
    std::printf("========== 波形 -> 点云测试 ==========\n");

    const std::string path = "data/event_waveform.txt";

    TextTraceReader reader;
    Check(reader.Open(path), "打开波形文件");
    Check(reader.NumEvents() == 10, "事件数应为 10");
    std::printf("    事件数=%d\n", reader.NumEvents());

    TraceEvent ev;
    Check(reader.ReadNext(ev), "读取事件 0");
    Check(ev.traces.size() > 100, "波形条数合理");
    std::printf("    事件 %d：%zu 条波形，每条约 %zu 个采样\n", ev.number,
                ev.traces.size(), ev.traces.empty() ? 0 : ev.traces[0].samples.size());

    // pad 平面 + 时间修正
    PadPlane pads;
    Check(pads.Load("data"), "加载 pad 几何");
    Check(pads.LoadTimeOffsets("data"), "加载 pad 时间修正");

    DetectorConfig cfg;
    PeakFinderParameters peak_params;
    peak_params.peak_separation = 50;
    peak_params.peak_prominence = 20.0;
    peak_params.peak_max_width = 50.0;
    peak_params.peak_threshold = 40.0;

    PointCloudBuilder builder;
    Check(builder.Init(cfg, &pads, peak_params, 20.0), "Init builder");

    PointCloud cloud;
    Check(builder.Build(ev, cloud), "Build 点云");
    Check(!cloud.points.empty(), "点云非空");

    int bad_pad = 0, bad_tb = 0, beam_hit = 0, bad_z = 0;
    double amp_max = 0.0, amp_min = 1e30;
    for (const auto& p : cloud.points) {
        if (p.pad < 0 || p.pad >= pads.NumPads()) ++bad_pad;
        if (p.tb < 0.0 || p.tb >= cfg.num_tb) ++bad_tb;
        if (pads.IsBeamPad(p.pad)) ++beam_hit;
        if (p.z < -1e-6 || p.z > cfg.length_m * 1000.0 + 1e-6) ++bad_z;
        amp_max = std::max(amp_max, p.amplitude);
        amp_min = std::min(amp_min, p.amplitude);
    }

    Check(bad_pad == 0, "pad 编号合法");
    Check(bad_tb == 0, "时间桶合法");
    Check(beam_hit == 0, "没有 beam pad");
    Check(bad_z == 0, "z 在灵敏体积内");
    Check(amp_min >= peak_params.peak_threshold, "幅度已过阈值");

    std::printf("    点数=%zu  接受波形=%d  峰值=%d  amp=[%.0f, %.0f]\n",
                cloud.points.size(), builder.AcceptedTraces(), builder.TotalPeaks(),
                amp_min, amp_max);

    // 已知亮点：pads 5 和 10 的原始幅度约 4088 / 3963
    bool has_pad5 = false, has_pad10 = false;
    for (const auto& p : cloud.points) {
        if (p.pad == 5) has_pad5 = true;
        if (p.pad == 10) has_pad10 = true;
    }
    Check(has_pad5, "亮点 pad 5 有点");
    Check(has_pad10, "亮点 pad 10 有点");

    // 同参数重复构建，结果应完全一致（随机数种子固定）
    PointCloud cloud2;
    PointCloudBuilder builder2;
    builder2.Init(cfg, &pads, peak_params, 20.0);
    builder2.Build(ev, cloud2);
    Check(cloud2.points.size() == cloud.points.size(), "重复构建点数一致");

    if (g_fail == 0) {
        std::printf("========== 全部通过 ==========\n");
        gSystem->Exit(0);
    } else {
        std::printf("========== 失败 %d 项 ==========\n", g_fail);
        gSystem->Exit(1);
    }
}
