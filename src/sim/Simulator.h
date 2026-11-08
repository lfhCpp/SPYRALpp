#ifndef RK4_SIM_SIMULATOR_H
#define RK4_SIM_SIMULATOR_H

//----------------------------------------------------------------------
// 探测器数字化：事件 -> 点云
//
// 流程：
//   粒子四动量 -> 轨迹（RK4 + 能损）-> 逐段能量沉积 -> 电子数（W/Fano/增益）
//   -> 漂移时间与横向扩散 -> pad 命中累加 -> GET 响应 + ADC 阈值 -> 点云
//----------------------------------------------------------------------

#include "../core/Event.h"
#include "../core/Rng.h"
#include "../core/EventMetadata.h"
#include "DetectorConfig.h"
#include "PadPlane.h"
#include "ParticleDynamics.h"
#include "Response.h"
#include "TargetDensity.h"

#include <unordered_map>
#include <vector>

namespace rk4 {

struct CloudPoint {
    int pad = -1;
    double tb = 0.0;                    // 时间桶
    double x = 0.0, y = 0.0, z = 0.0;   // mm
    double amplitude = 0.0;             // ADC
    double integral = 0.0;              // ADC·TB
    double scale = 1.0;                 // 大小 pad 标记
    size_t dither_index = 0; // Transient raw frontend ordinal; not serialized.
};

struct PointCloud {
    long long event = -1;
    std::vector<CloudPoint> points;
    size_t dither_draws = 0; // Includes rejected peaks and excluded pads.
    bool pending_dither = false; // Transient; resolved before writing clouds.
    EventMetadata metadata;
};

class DetectorSimulator {
public:
    // gas 为靶气体（提供密度与组成）；可为 nullptr 则用 ParticleDynamics 默认
    bool Init(const DetectorConfig& cfg, const Target* gas);

    // indices 为空时默认模拟出射粒子与剩余核（2,3）
    bool Simulate(const Event& ev, PointCloud& pc,
                  const std::vector<int>& indices = std::vector<int>());

    double TotalElectrons() const { return total_electrons_; }
    void SetEventSeed(uint64_t run_seed, uint64_t event) {
        cfg_.seed = run_seed;
        rng_.Seed(run_seed, 2*event+2);
    }
    const PadPlane& Pads() const { return pads_; }

private:
    void Deposit(double x_mm, double y_mm, double tb, double electrons, double sigma_mm);

    DetectorConfig cfg_;
    PadPlane pads_;
    GetResponse response_;
    ParticleDynamics dynamics_;
    core::Rng rng_;
    std::unordered_map<long long, double> charge_; // (pad, tb) -> 电子数
    double total_electrons_ = 0.0;
    bool ready_ = false;
};

} // namespace rk4

#endif // RK4_SIM_SIMULATOR_H
