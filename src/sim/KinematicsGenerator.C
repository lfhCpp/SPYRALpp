#include "KinematicsGenerator.h"

#include "../core/Units.h"

#include <cmath>
#include <cstdio>

namespace rk4 {

KinematicsGenerator::KinematicsGenerator() = default;
KinematicsGenerator::~KinematicsGenerator() = default;

bool KinematicsGenerator::Init(const GeneratorConfig& cfg, const Target* target,
                               uint64_t seed) {
    ready_ = false;
    event_id_ = accepted_ = trials_ = 0;
    cfg_ = cfg;
    target_ = target;
    rng_.Seed(seed);

    if (!reaction_.Init(cfg.reaction)) {
        std::fprintf(stderr, "[KinematicsGenerator] 反应初始化失败\n");
        return false;
    }

    excitation_.reset(new ExcitationGaussian(cfg.ex_centroid, cfg.ex_fwhm));
    polar_.reset(new PolarUniform(cfg.theta_min_deg * Units::DEG2RAD,
                                  cfg.theta_max_deg * Units::DEG2RAD));

    if (cfg.apply_target_loss && target_ != nullptr) {
        if (!stopping_.Init(cfg.reaction.beam_A, cfg.reaction.beam_Z, *target_)) {
            std::fprintf(stderr, "[KinematicsGenerator] 阻止本领初始化失败\n");
            return false;
        }
    }

    if (cfg.beam_ke_mev <= 0.0) {
        std::fprintf(stderr, "[KinematicsGenerator] 束流能量必须为正\n");
        return false;
    }

    ready_ = true;
    return true;
}

double KinematicsGenerator::BeamEnergyAt(double depth_m) const {
    if (cfg_.apply_target_loss && target_ != nullptr) {
        return stopping_.EnergyAfter(cfg_.beam_ke_mev, depth_m);
    }
    return cfg_.beam_ke_mev;
}

bool KinematicsGenerator::Generate(Event& ev) {
    if (!ready_) return false;

    const double two_pi = 2.0 * Units::PI;

    for (int trial = 0; trial < cfg_.max_trials; ++trial) {
        ++trials_;

        // 1) 顶点抽样
        const double z = rng_.Uniform(cfg_.z_min_m, cfg_.z_max_m);
        const double rho = std::fabs(rng_.Gauss(0.0, cfg_.rho_sigma_m));
        const double phi_v = rng_.Uniform(0.0, two_pi);

        // 2) 束流到达该 z 的剩余能量
        const double ke = BeamEnergyAt(z);

        // 3) 激发能与角度抽样
        const double ex = excitation_->Sample(rng_);
        if (!reaction_.IsAllowed(ke, ex)) continue;

        const double theta = polar_->SamplePolar(rng_);
        const double phi = rng_.Uniform(0.0, two_pi);

        // 4) 二体运动学
        if (!reaction_.Calculate(ke, theta, phi, ex, ev)) continue;

        // 5) 顶点信息
        ev.id = event_id_++;
        ev.vx = rho * std::cos(phi_v);
        ev.vy = rho * std::sin(phi_v);
        ev.vz = z;
        ev.beam_energy = ke;

        ++accepted_;
        return true;
    }

    return false;
}

} // namespace rk4
