#ifndef RK4_IO_PIDSTATEIO_H
#define RK4_IO_PIDSTATEIO_H

#include "../reco/PidState.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace rk4 {

struct PidCandidateAuditRow {
    long long event = -1;
    int cluster_index = -1;
    int cluster_label = -1;
    int n_points = 0;
    int estimate_ok = 0;
    double ic_amplitude = -1.0;
    IcObservationState ic_state = IcObservationState::NOT_RECORDED;
    int ic_gate_configured = 0;
    int ic_gate_pass = -1;
    int pid_gate_configured = 0;
    int pid_required_index = -1;
    double pid_min_separation = 0.0;
    int pid_calibration_available = 0;
    int pid_classifier_evaluated = 0;
    std::string pid_state_basis = "estimate_brho";
    int pid_gate_pass = -1;
    int pid_best_index = -1;
    double pid_separation = 0.0;
    int fit_attempted = 0;
    FitState fit_state = FitState::FIT_REJECTED;
    std::string fit_rejection_stage = "not_run";
    double fit_cost_mm = 0.0;
    int migrad_status = -1;
    double migrad_edm = 0.0;
    PidState pid_state = PidState::PID_UNKNOWN;
    PidUnknownReason pid_unknown_reason = PidUnknownReason::NOT_CONFIGURED;
};

inline bool WritePidCandidateAuditCsv(const std::string& path,
                                     const std::vector<PidCandidateAuditRow>& rows) {
    std::ofstream out(path);
    if (!out) return false;
    out << "event,cluster_index,cluster_label,n_points,estimate_ok,ic_amplitude,ic_state,ic_gate_configured,ic_gate_pass"
           ",pid_gate_configured,pid_required_index,pid_min_separation,pid_calibration_available"
           ",pid_classifier_evaluated,pid_gate_pass,pid_best_index,pid_separation,fit_attempted"
           ",fit_state,fit_rejection_stage,fit_cost_mm,migrad_status,migrad_edm,pid_state,pid_unknown_reason,pid_state_basis\n"
        << std::setprecision(17);
    for (const auto& row : rows) {
        out << row.event << ',' << row.cluster_index << ',' << row.cluster_label << ','
            << row.n_points << ',' << row.estimate_ok << ','
            << row.ic_amplitude << ',' << ToString(row.ic_state) << ',' << row.ic_gate_configured << ','
            << row.ic_gate_pass << ',' << row.pid_gate_configured << ',' << row.pid_required_index << ','
            << row.pid_min_separation << ',' << row.pid_calibration_available << ','
            << row.pid_classifier_evaluated << ',' << row.pid_gate_pass << ',' << row.pid_best_index << ','
            << row.pid_separation << ',' << row.fit_attempted << ',' << ToString(row.fit_state) << ','
            << row.fit_rejection_stage << ',' << row.fit_cost_mm << ',' << row.migrad_status << ','
            << row.migrad_edm << ',' << ToString(row.pid_state) << ','
            << ToString(row.pid_unknown_reason) << ',' << row.pid_state_basis << '\n';
    }
    out.close();
    return static_cast<bool>(out);
}

} // namespace rk4

#endif // RK4_IO_PIDSTATEIO_H
