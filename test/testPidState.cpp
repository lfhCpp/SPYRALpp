#include "reco/PidState.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void Check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
}

int main() {
    try {
        using namespace rk4;
        Check(ClassifyIcObservation(-1.0) == IcObservationState::MISSING_OBSERVATION,
              "IC sentinel -1 must be missing");
        Check(ClassifyIcObservation(0.0) == IcObservationState::IC_OBSERVED,
              "finite zero remains an observed IC value");
        Check(ClassifyIcObservation(-2.0) == IcObservationState::INVALID_OBSERVATION,
              "negative non-sentinel IC must be invalid");
        Check(ClassifyIcObservation(std::numeric_limits<double>::quiet_NaN()) ==
                  IcObservationState::INVALID_OBSERVATION,
              "NaN IC must be invalid");
        Check(-1.0 > -std::numeric_limits<double>::infinity() &&
                  -1.0 < std::numeric_limits<double>::infinity(),
              "legacy unbounded IC range must continue accepting the sentinel");

        Check(PassesExistingPidGate(false, false, 0.0, 1.0, -1, 0),
              "disabled PID gate must remain a fit pass-through");
        Check(!PassesExistingPidGate(true, false, 10.0, 0.0, 0, 0),
              "enabled PID gate must reject an unevaluated observation");
        Check(!PassesExistingPidGate(true, true, 0.99, 1.0, 0, 0),
              "enabled PID gate must preserve separation threshold");
        Check(!PassesExistingPidGate(true, true, 1.0, 1.0, 1, 0),
              "enabled PID gate must preserve required species check");
        Check(PassesExistingPidGate(true, true, 1.0, 1.0, 0, 0),
              "matching PID candidate must pass the existing gate");

        PidStateInput input;
        input.observation_present = true;
        input.observation_valid = true;
        input.classifier_evaluated = true;
        input.target_matches = true;
        input.separation_passes = true;
        Check(ClassifyPidState(input).unknown_reason == PidUnknownReason::NOT_CONFIGURED,
              "gate-off state must be unknown/not configured, never PID-confirmed");
        input.gate_configured = true;
        Check(ClassifyPidState(input).unknown_reason == PidUnknownReason::NOT_CONFIGURED,
              "no target species must remain not configured");
        input.target_species_configured = true;
        input.observation_present = false;
        Check(ClassifyPidState(input).unknown_reason == PidUnknownReason::MISSING_OBSERVATION,
              "missing PID observations need an explicit unknown reason");
        input.observation_present = true;
        input.observation_valid = false;
        Check(ClassifyPidState(input).unknown_reason == PidUnknownReason::INVALID_OBSERVATION,
              "invalid PID observations need an explicit unknown reason");
        input.observation_valid = true;
        input.calibration_available = false;
        Check(ClassifyPidState(input).unknown_reason == PidUnknownReason::CALIBRATION_UNAVAILABLE,
              "unavailable calibration must prevent confirmation");
        input.calibration_available = true;
        Check(ClassifyPidState(input).state == PidState::PID_CONFIRMED,
              "configured calibrated matching observation should confirm");
        input.target_matches = false;
        Check(ClassifyPidState(input).state == PidState::PID_REJECTED,
              "configured calibrated nonmatching observation should reject PID");

        Check(PassesExistingMinuitGate(0.2, 19.999), "legacy 20 mm gate must accept below threshold");
        Check(!PassesExistingMinuitGate(0.2, 20.0), "legacy 20 mm gate must reject equality");
        Check(!PassesExistingMinuitGate(0.0, 0.0), "legacy Minuit gate must reject nonpositive B-rho");
        Check(!PassesExistingMinuitGate(std::numeric_limits<double>::quiet_NaN(), 1.0),
              "legacy Minuit gate must reject nonfinite B-rho");

        std::cout << "PASS: IC, PID unknown reasons, legacy gate and independent fit/PID states\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
