#ifndef RK4_RECO_TRACKFITTER_H
#define RK4_RECO_TRACKFITTER_H

//----------------------------------------------------------------------
// 轨迹拟合（对应 Spyral 的 InterpSolverPhase）
//
// 用预生成的轨迹网格 + 双线性插值构造模型轨迹，
// 以独立的 Minuit2 实例最小化"数据点到模型轨迹的平均最近距离"。
//
// 拟合参数：Bρ、极角、顶点 (x, y, z) 与方位角。
//----------------------------------------------------------------------

#include "Cluster.h"
#include "Estimator.h"
#include "TrackMesh.h"


#include <vector>
#include <string>
#include <cstdint>

namespace rk4 {

struct FitResult {
    bool ok = false;
    double brho = 0.0;      // Tm
    double sigma_brho = 0.0;
    double polar = 0.0;     // rad
    double sigma_polar = 0.0;
    double azimuthal = 0.0; // rad
    double vx = 0.0, vy = 0.0, vz = 0.0; // mm
    double cost = 0.0;      // 平均最近距离 (mm)
    double optimizer_objective = 0.0; // Minuit search objective; may differ in experiments
    double edm = 0.0;
    int status = 0;
};

// Optional read-only diagnostics for the native Minuit2 fit. This records the
// existing two-stage/five-start search without changing its decisions.
struct CandidateSearchCounters {
    std::uint64_t objective_evaluations_with_curve = 0;
    std::uint64_t objective_evaluations_with_flat_tail = 0;
    std::uint64_t flat_tail_segments_total = 0;
    std::uint64_t flat_tail_segments_max = 0;
    std::uint64_t evaluated_hits = 0;
    std::uint64_t hits_with_tail_candidate = 0;
    std::uint64_t hits_with_multiple_tail_candidates = 0;
    std::uint64_t tail_candidate_slots = 0;
    std::uint64_t multiple_tail_hits_checked_exactly = 0;
    std::uint64_t multiple_tail_hits_missed_exact = 0;
    double multiple_tail_miss_distance_sum_mm = 0.0;
    double multiple_tail_miss_distance_max_mm = 0.0;
    std::uint64_t objective_evaluations_with_multiple_tail_miss = 0;
    double multiple_tail_miss_objective_delta_sum_mm = 0.0;
    double multiple_tail_miss_objective_delta_max_mm = 0.0;
};

struct MinuitPhaseDiagnostics {
    bool minimize_return = false;
    int status = -1;
    int covariance_status = -1;
    unsigned int function_calls = 0;
    double objective = 0.0;
    double edm = 0.0;
    double final_values[6]{}; // brho, polar, vx, vy, vz, azimuth
    double final_errors[6]{};
    CandidateSearchCounters candidate_search;
};

struct MinuitStartDiagnostics {
    int start_index = -1;
    double azimuth_offset = 0.0;
    double initial_azimuth = 0.0;
    MinuitPhaseDiagnostics fixed_vertex;
    MinuitPhaseDiagnostics free_vertex;
    bool fresh_free_vertex_ran = false;
    MinuitPhaseDiagnostics fresh_free_vertex;
    double initial_values[6]{};
    double lower_bounds[6]{};
    double upper_bounds[6]{};
    double final_values[6]{}; // brho, polar, vx, vy, vz, azimuth
    double final_errors[6]{};
    bool below_existing_cost_gate = false;
    bool selected_as_best = false;
    double mean_distance_cost = 0.0; // Exact legacy score used by the unchanged 20 mm gate
};

struct MinuitFitDiagnostics {
    bool attempted = false;
    bool accepted = false;
    std::string rejection_reason;
    size_t point_count = 0;
    double final_cost = 0.0;
    double final_optimizer_objective = 0.0;
    int best_start_index = -1;
    std::vector<MinuitStartDiagnostics> starts;
};

struct CandidateOrderObjectiveComparison {
    double insertion_order_mm = 0.0;
    double sorted_order_mm = 0.0;
};

class TrackFitter {
public:
    bool Init(const TrackMesh* mesh, double mass_mev, int charge);

    // Keep the local nearest-segment candidates in ascending distance order.
    // This is the production default; the legacy setter remains for old tools.
    void SetSortCandidates(bool enabled) { sort_candidates_on_merge_ = enabled; }
    void SetExperimentalSortCandidates(bool enabled) { SetSortCandidates(enabled); }

    // Experimental Gaussian-residual objective. The final legacy mean distance
    // is still recomputed for acceptance; disabled preserves production behavior.
    void SetExperimentalMeanSquaredDistanceObjective(bool enabled) {
        experimental_mean_squared_distance_objective_ = enabled;
    }

    // Experimental: when multiple sampled centers land in a repeated flat
    // propagation tail, rescan all segments and use the exact nearest one.
    void SetExperimentalExactTailCandidates(bool enabled) {
        exact_tail_candidates_ = enabled;
    }

    // Optional experimental circular bound on the fitted transverse vertex.
    // A value <= 0 preserves the production x/y parameterization and bounds.
    void SetExperimentalVertexRadiusLimit(double radius_mm) {
        experimental_vertex_radius_limit_mm_ = radius_mm;
    }

    // Optional tangential vertex-angle step target for the circular
    // parameterization. Zero preserves its legacy 0.02 rad Minuit step.
    void SetExperimentalVertexPhiStepScale(double tangential_step_mm) {
        experimental_vertex_phi_step_mm_ = tangential_step_mm;
    }

    // Move the bounded vertex-angle interval so its periodic seam is opposite
    // the Estimation seed. Disabled by default.
    void SetExperimentalCenteredVertexPhiBounds(bool enabled) {
        experimental_centered_vertex_phi_bounds_ = enabled;
    }

    // Center the full periodic particle-azimuth interval around each
    // multistart seed. This removes the arbitrary 0/2*pi Minuit seam.
    void SetCenteredAzimuthBounds(bool enabled) { centered_azimuth_bounds_ = enabled; }
    // Source-compatible alias for older experiment drivers.
    void SetExperimentalCenteredAzimuthBounds(bool enabled) {
        SetCenteredAzimuthBounds(enabled);
    }
    void SetExperimentalAzimuthBoundsCenterShift(double shift_rad) {
        experimental_azimuth_bounds_center_shift_rad_ = shift_rad;
    }
    void SetExperimentalUnboundedAzimuth(bool enabled) {
        experimental_unbounded_azimuth_ = enabled;
    }
    // Run a separate free-vertex Migrad from the fixed-stage solution. This
    // diagnostic never replaces or affects the ordinary fit result.
    void SetExperimentalFreshFreeVertexMinuit(bool enabled) {
        experimental_fresh_free_vertex_minuit_ = enabled;
    }

    bool Fit(const std::vector<CloudPoint>& data, const EstimateResult& guess,
             double azimuthal, FitResult& out,
             MinuitFitDiagnostics* diagnostics = nullptr);

    // Python Spyral's discrete weighted-distance objective and lmfit bounds.
    bool FitPythonEquivalent(const std::vector<CloudPoint>& data,
                             const EstimateResult& guess, double z_error_m,
                             double beam_radius_m, FitResult& out);

    // Objective function used by the per-instance Minuit2 functor
    double Evaluate(const double* par);

    // Evaluate identical fixed parameters and points with both candidate-list
    // merge policies. This diagnostic does not run Minuit or alter Fit defaults.
    CandidateOrderObjectiveComparison EvaluateBothCandidateOrdersForDiagnostics(
        const std::vector<CloudPoint>& data, const double* par);

private:
    const TrackMesh* mesh_ = nullptr;
    double mass_mev_ = 0.0;
    int charge_ = 1;
    bool sort_candidates_on_merge_ = true;
    bool exact_tail_candidates_ = false;
    double experimental_vertex_radius_limit_mm_ = 0.0;
    double experimental_vertex_phi_step_mm_ = 0.0;
    bool experimental_centered_vertex_phi_bounds_ = false;
    bool centered_azimuth_bounds_ = false;
    double experimental_azimuth_bounds_center_shift_rad_ = 0.0;
    bool experimental_unbounded_azimuth_ = false;
    bool experimental_fresh_free_vertex_minuit_ = false;
    bool experimental_mean_squared_distance_objective_ = false;

    const std::vector<CloudPoint>* data_ = nullptr;
    double azimuthal_ = 0.0;
    bool collect_candidate_search_diagnostics_ = false;
    CandidateSearchCounters candidate_search_counters_;

    void ResetCandidateSearchCounters();
};

} // namespace rk4

#endif // RK4_RECO_TRACKFITTER_H
