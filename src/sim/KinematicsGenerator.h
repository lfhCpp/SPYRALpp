#ifndef RK4_SIM_KINEMATICSGENERATOR_H
#define RK4_SIM_KINEMATICSGENERATOR_H

//----------------------------------------------------------------------
// 运动学生成器：对应 attpc_engine 的 KinematicsPipeline
//
// 流程：
//   1. 在靶体积内抽样顶点 (rho, phi, z)
//   2. 计算束流穿到该 z 时的剩余能量（catima）
//   3. 抽样激发能、质心系极角与方位角
//   4. 做二体运动学，得到 4 个粒子
//   5. 若能量上不允许则整条事件重采（拒绝采样）
//----------------------------------------------------------------------

#include "../core/Distributions.h"
#include "../core/Event.h"
#include "../core/Reaction.h"
#include "../core/Rng.h"
#include "Stopping.h"
#include "TargetDensity.h"

#include <cstdint>
#include <memory>

namespace rk4 {

struct GeneratorConfig {
    ReactionConfig reaction;

    double beam_ke_mev = 0.0;   // 加速器能量（靶窗外），MeV

    double ex_centroid = 0.0;   // 激发能均值，MeV
    double ex_fwhm = 0.0;       // 激发能 FWHM，MeV（0 表示单能）

    double theta_min_deg = 0.0; // 质心系极角范围
    double theta_max_deg = 180.0;

    double z_min_m = 0.0;       // 顶点 z 范围，m
    double z_max_m = 1.0;
    double rho_sigma_m = 0.007; // 顶点横截面高斯 sigma，m

    bool apply_target_loss = true; // 是否计算靶中能损
    int max_trials = 1000;         // 单事件最多重采次数
};

class KinematicsGenerator {
public:
    KinematicsGenerator();
    ~KinematicsGenerator();

    // target 可为 nullptr（此时 apply_target_loss 会被忽略）
    bool Init(const GeneratorConfig& cfg, const Target* target, uint64_t seed = 1);

    // 生成一条事件；超过重采上限或未初始化时返回 false
    bool Generate(Event& ev);
    void SetEventSeed(uint64_t run_seed, uint64_t event) { rng_.Seed(run_seed, 2*event+1); }

    // 束流穿到 depth（m）时的动能，MeV
    double BeamEnergyAt(double depth_m) const;

    const TwoBodyReaction& Reaction() const { return reaction_; }
    long long Accepted() const { return accepted_; }
    long long Trials() const { return trials_; }

private:
    GeneratorConfig cfg_;
    TwoBodyReaction reaction_;
    std::unique_ptr<ExcitationDistribution> excitation_;
    std::unique_ptr<PolarDistribution> polar_;
    StoppingPower stopping_;
    const Target* target_ = nullptr;
    core::Rng rng_;
    long long event_id_ = 0;
    long long accepted_ = 0;
    long long trials_ = 0;
    bool ready_ = false;
};

} // namespace rk4

#endif // RK4_SIM_KINEMATICSGENERATOR_H
