#ifndef RK4_SIM_DETECTORCONFIG_H
#define RK4_SIM_DETECTORCONFIG_H

//----------------------------------------------------------------------
// 探测器与电子学参数（M2）
//
// 默认值取自 attpc_engine 的 getting_started 示例与现有 ParticleDynamics 设置。
//----------------------------------------------------------------------

#include <string>

namespace rk4 {

struct DetectorConfig {
    // ---- 几何与场 ----
    double length_m = 1.0;      // 灵敏体积长度
    double bfield_t = 2.85;     // 磁场
    double efield_vm = 55000.0; // Run 0151 基线电场（与重建 fallback 一致，V/m）

    // ---- 气体与雪崩 ----
    double w_value_ev = 34.0;   // 产生一个电子-离子对所需平均能量
    double fano = 0.2;          // Fano 因子
    double diffusion = 0.277;   // 横向扩散参数（与 attpc_engine 同口径）
    int mpgd_gain = 175000;     // MPGD 总增益

    // ---- 电子学 ----
    double clock_mhz = 6.25;        // GET 采样率
    double amp_gain = 900.0;        // 放大器增益 (lsb/fC)
    double shaping_ns = 1000.0;     // 成形时间
    double micromegas_edge = 10.0;  // micromegas 时间参考（TB）
    double windows_edge = 560.0;    // 靶窗时间参考（TB）
    double adc_threshold = 10.0;    // 点云振幅阈值
    int num_tb = 512;               // 时间桶总数
    bool dither_time_buckets = true; // Represent a detected time bin by [tb,tb+1), as GET does.

    // ---- 数据 ----
    std::string pad_dir = "";       // padxy.csv / pad_scale.csv 所在目录；空 = 用 DataDir()

    // ---- 随机数 ----
    unsigned long long seed = 12345;
};

} // namespace rk4

#endif // RK4_SIM_DETECTORCONFIG_H
