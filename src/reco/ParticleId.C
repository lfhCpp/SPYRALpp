#include "ParticleId.h"

#include "../core/MassData.h"
#include "../core/Units.h"

#include <cmath>
#include <limits>

namespace rk4 {

namespace {

constexpr double kElectronMassMeV = 0.51099895;
constexpr double kIonizationEv = 40.0e-6; // MeV，气体近似

inline double BetheDedx(double beta, double gamma, int Z) {
    if (beta <= 1.0e-4) return 0.0;
    const double bg2 = beta * beta * gamma * gamma;
    const double arg = 2.0 * kElectronMassMeV * bg2 / kIonizationEv;
    if (arg <= 1.0) return 0.0;

    const double z2 = static_cast<double>(Z) * Z;
    const double value = z2 / (beta * beta) * (std::log(arg) - beta * beta);
    return value > 0.0 ? value : 0.0;
}

} // namespace

ParticleId::ParticleId() {
    candidates_ = {
        {"p", 1, 1},
        {"d", 2, 1},
        {"t", 3, 1},
        {"3He", 3, 2},
        {"alpha", 4, 2},
    };
}

void ParticleId::SetCandidates(const std::vector<PidCandidate>& candidates) {
    candidates_ = candidates;
}

double ParticleId::ExpectedDedx(double brho, int A, int Z) {
    if (brho <= 0.0 || A <= 0 || Z <= 0) return 0.0;

    const double mass = MassData::GetNuclearMassNuclear(A, Z) * Units::u_MeV;
    if (mass <= 0.0) return 0.0;

    // p = 299.8 * Bρ * Z (MeV/c)
    const double p = 299.8 * brho * Z;
    const double gamma = std::sqrt(p * p + mass * mass) / mass;
    if (gamma <= 1.0) return 0.0;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));

    return BetheDedx(beta, gamma, Z);
}

void ParticleId::Calibrate(double brho, double dedx, int A, int Z) {
    const double expected = ExpectedDedx(brho, A, Z);
    if (expected > 1.0e-9 && dedx > 0.0) {
        scale_ = dedx / expected;
    }
}

bool ParticleId::Identify(double brho, double dedx, PidResult& out) const {
    out = PidResult{};
    if (brho <= 0.0 || dedx <= 0.0 || candidates_.empty()) return false;

    double best = std::numeric_limits<double>::max();
    double second = std::numeric_limits<double>::max();
    int best_index = -1;
    double log_obs = std::log(dedx);

    for (size_t i = 0; i < candidates_.size(); ++i) {
        const double expected = scale_ * ExpectedDedx(brho, candidates_[i].A, candidates_[i].Z);
        if (expected <= 1.0e-9) continue;

        const double score = std::fabs(log_obs - std::log(expected));
        if (score < best) {
            second = best;
            best = score;
            best_index = static_cast<int>(i);
        } else if (score < second) {
            second = score;
        }
    }

    if (best_index < 0) return false;

    out.ok = true;
    out.index = best_index;
    out.name = candidates_[best_index].name;
    out.score = best;
    out.second_score = second;
    out.separation = second - best;
    return true;
}

} // namespace rk4
