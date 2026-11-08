#ifndef RK4_RECO_PIDSTATE_H
#define RK4_RECO_PIDSTATE_H

#include <cmath>

namespace rk4 {

enum class FitState : int { FIT_REJECTED = 0, FIT_ACCEPTED = 1, NOT_RECORDED = -1 };
enum class PidState : int { PID_CONFIRMED = 1, PID_REJECTED = 2, PID_UNKNOWN = 3, NOT_RECORDED = -1 };
enum class PidUnknownReason : int {
    NONE = 0,
    NOT_CONFIGURED = 1,
    MISSING_OBSERVATION = 2,
    INVALID_OBSERVATION = 3,
    CALIBRATION_UNAVAILABLE = 4,
    NOT_RECORDED = -1
};
enum class IcObservationState : int {
    IC_OBSERVED = 1,
    MISSING_OBSERVATION = 2,
    INVALID_OBSERVATION = 3,
    NOT_RECORDED = -1
};
enum class PidStateBasis : int { ESTIMATE = 0, FITTED_BRHO = 1, NOT_RECORDED = -1 };

struct PidStateInput {
    bool gate_configured = false;
    bool target_species_configured = false;
    bool observation_present = false;
    bool observation_valid = false;
    bool calibration_available = false;
    bool classifier_evaluated = false;
    bool target_matches = false;
    bool separation_passes = false;
};

struct PidStateResult {
    PidState state = PidState::PID_UNKNOWN;
    PidUnknownReason unknown_reason = PidUnknownReason::NOT_CONFIGURED;
};

inline IcObservationState ClassifyIcObservation(double value) {
    if (!std::isfinite(value)) return IcObservationState::INVALID_OBSERVATION;
    if (value == -1.0) return IcObservationState::MISSING_OBSERVATION;
    if (value < 0.0) return IcObservationState::INVALID_OBSERVATION;
    return IcObservationState::IC_OBSERVED;
}

inline PidStateResult ClassifyPidState(const PidStateInput& input) {
    if (!input.gate_configured || !input.target_species_configured)
        return {PidState::PID_UNKNOWN, PidUnknownReason::NOT_CONFIGURED};
    if (!input.observation_present)
        return {PidState::PID_UNKNOWN, PidUnknownReason::MISSING_OBSERVATION};
    if (!input.observation_valid)
        return {PidState::PID_UNKNOWN, PidUnknownReason::INVALID_OBSERVATION};
    if (!input.calibration_available)
        return {PidState::PID_UNKNOWN, PidUnknownReason::CALIBRATION_UNAVAILABLE};
    if (!input.classifier_evaluated)
        return {PidState::PID_UNKNOWN, PidUnknownReason::INVALID_OBSERVATION};
    if (input.target_matches && input.separation_passes)
        return {PidState::PID_CONFIRMED, PidUnknownReason::NONE};
    return {PidState::PID_REJECTED, PidUnknownReason::NONE};
}

// This helper preserves the existing pre-fit PID gate's short-circuit order.
// In particular, a disabled gate is pass-through for fitting, not PID confirmation.
inline bool PassesExistingPidGate(bool gate_enabled, bool identify_ok,
                                  double separation, double min_separation,
                                  int candidate_index, int required_index) {
    if (!gate_enabled) return true;
    if (!identify_ok) return false;
    if (separation < min_separation) return false;
    if (required_index >= 0 && candidate_index != required_index) return false;
    return true;
}

inline bool PassesExistingMinuitGate(double brho, double mean_distance_mm) {
    return std::isfinite(brho) && brho > 0.0 && mean_distance_mm < 20.0;
}

inline const char* ToString(FitState state) {
    switch (state) {
        case FitState::FIT_ACCEPTED: return "FIT_ACCEPTED";
        case FitState::FIT_REJECTED: return "FIT_REJECTED";
        default: return "NOT_RECORDED";
    }
}

inline const char* ToString(PidState state) {
    switch (state) {
        case PidState::PID_CONFIRMED: return "PID_CONFIRMED";
        case PidState::PID_REJECTED: return "PID_REJECTED";
        case PidState::PID_UNKNOWN: return "PID_UNKNOWN";
        default: return "NOT_RECORDED";
    }
}

inline const char* ToString(PidUnknownReason reason) {
    switch (reason) {
        case PidUnknownReason::NONE: return "NONE";
        case PidUnknownReason::NOT_CONFIGURED: return "NOT_CONFIGURED";
        case PidUnknownReason::MISSING_OBSERVATION: return "MISSING_OBSERVATION";
        case PidUnknownReason::INVALID_OBSERVATION: return "INVALID_OBSERVATION";
        case PidUnknownReason::CALIBRATION_UNAVAILABLE: return "CALIBRATION_UNAVAILABLE";
        default: return "NOT_RECORDED";
    }
}

inline const char* ToString(IcObservationState state) {
    switch (state) {
        case IcObservationState::IC_OBSERVED: return "IC_OBSERVED";
        case IcObservationState::MISSING_OBSERVATION: return "MISSING_OBSERVATION";
        case IcObservationState::INVALID_OBSERVATION: return "INVALID_OBSERVATION";
        default: return "NOT_RECORDED";
    }
}

} // namespace rk4

#endif // RK4_RECO_PIDSTATE_H
