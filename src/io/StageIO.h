#ifndef RK4_IO_STAGEIO_H
#define RK4_IO_STAGEIO_H

//----------------------------------------------------------------------
// 中间过程文件（可选，用来分阶段检查/重跑）
//
// 整条链现在落盘的是：
//   clouds.root       rk4_waveform  -> 点云（每事件一个 vector）
//   *_stages.root     rk4_analyze   -> cluster / estimate 两个阶段（本文件）
//   fit.root          rk4_analyze   -> 拟合结果
//   *_peaks.root      rk4_waveform  -> 每条迹线的峰参数（含饱和标志，本文件）
//
// 打开方式：配置里 `[output] stages = 1` / `[peaks] dump = 1`。
// 关闭时整条链行为不变（不写、不占内存）。
//----------------------------------------------------------------------

#include <string>
#include <vector>

namespace rk4 {

//----------------------------------------------------------------------
// 聚类阶段：一个簇 + 它的所有点（便于离线重拟合）
//----------------------------------------------------------------------
struct ClusterRecord {
    long long event = -1;
    int cluster = -1;
    int npoints = 0;
    std::vector<double> x, y, z, amplitude, integral, tb;
    int direction = -1; // 0 = Forward, 1 = Backward, -1 = 未定
};

//----------------------------------------------------------------------
// 估计阶段：粗估参数
//----------------------------------------------------------------------
struct EstimateRecord {
    long long event = -1;
    int cluster = -1;
    int ok = 0;
    double vx = 0.0, vy = 0.0, vz = 0.0;      // mm
    double polar = 0.0, azimuthal = 0.0;      // rad
    double brho = 0.0;                        // T·m
    double radius = 0.0;                      // mm
    double dEdx = 0.0;                        // 电荷/mm
    double arclength = 0.0;                   // mm
    double charge = 0.0;
    int direction = -1;
};

//----------------------------------------------------------------------
// 波形阶段：一条迹线上的峰（模板拟合时还带残差与饱和标志）
//----------------------------------------------------------------------
struct TracePeakRecord {
    long long event = -1;
    int pad = -1;
    int npeaks = 0;    // 这条迹线上共几个峰
    int index = 0;     // 第几个峰
    double time = 0.0; // 时间桶（含 pad 时间修正）
    double amplitude = 0.0;
    double charge = 0.0;
    double chi2 = 0.0;
    int npoints = 0;   // 参与拟合的样本数
    int saturated = 0; // 该迹线是否饱和
};

bool WriteStages(const std::string& path, const std::vector<ClusterRecord>& clusters,
                 const std::vector<EstimateRecord>& estimates);

bool WriteTracePeaks(const std::string& path, const std::vector<TracePeakRecord>& peaks);

} // namespace rk4

#endif // RK4_IO_STAGEIO_H
