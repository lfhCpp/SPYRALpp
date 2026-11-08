#include "TrackFitter.h"
#include "PidState.h"

#include "../core/Units.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"
#include <algorithm>
#include <array>

#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>

namespace rk4 {

namespace {

constexpr double kBrhoToMomentum = Units::C * 1.0e-6;

// 点到线段的距离平方
inline double PointSegmentDistance2(double px, double py, double pz, double ax, double ay,
                                    double az, double bx, double by, double bz) {
    const double abx = bx - ax, aby = by - ay, abz = bz - az;
    const double apx = px - ax, apy = py - ay, apz = pz - az;
    const double ab2 = abx * abx + aby * aby + abz * abz;
    double t = 0.0;
    if (ab2 > 1e-12) t = (apx * abx + apy * aby + apz * abz) / ab2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    const double dx = apx - t * abx;
    const double dy = apy - t * aby;
    const double dz = apz - t * abz;
    return dx * dx + dy * dy + dz * dz;
}

} // namespace

bool TrackFitter::Init(const TrackMesh* mesh, double mass_mev, int charge) {
    mesh_ = nullptr;
    if (mesh == nullptr || mesh->Empty()) return false;
    if (!(mass_mev>0) || !std::isfinite(mass_mev) || charge<1) return false;
    mesh_ = mesh;
    mass_mev_ = mass_mev;
    charge_ = charge > 0 ? charge : 1;
    return true;
}

double TrackFitter::Evaluate(const double* par) {
    if (mesh_ == nullptr || data_ == nullptr) return 1.0e9;

    const double brho = par[0];
    const double polar = par[1];
    const double vx = par[2];
    const double vy = par[3];
    const double vz = par[4];
    const double azimuth = par[5];

    const double p = kBrhoToMomentum * brho * charge_; // MeV/c
    if (p <= 0.0) return 1.0e9;

    const double ke = p*p/(std::hypot(p,mass_mev_)+mass_mev_);
    if (ke < 0.0 || !mesh_->IsInRange(ke, polar)) return 1.0e9;

    std::vector<double> xyz;
    const bool python_grid = mesh_->IsPythonMesh();
    if (python_grid ? !mesh_->InterpolatePython(ke, polar, xyz)
                    : !mesh_->Interpolate(ke, polar, xyz)) return 1.0e9;
    // Python's shared trajectory grid stores positions in metres; the native
    // Minuit objective and point cloud use millimetres.
    if (python_grid) for (double& value : xyz) value *= 1000.0;

    const double cz = std::cos(azimuth);
    const double sz = std::sin(azimuth);
    // Native TrackMesh stores positions around z_ref; the Python NumPy mesh
    // stores displacement from its z=0 start and FitPythonEquivalent adds the
    // fitted vertex directly. Keep the two reference frames distinct.
    const double z_shift = python_grid ? vz : vz - mesh_->Config().z_ref * 1000.0;

    const int n_time = mesh_->TimeSteps();
    std::vector<std::array<double, 3>> model(static_cast<size_t>(n_time));
    for (int t = 0; t < n_time; ++t) {
        const size_t index = static_cast<size_t>(t) * 3;
        model[static_cast<size_t>(t)] = {
            xyz[index] * cz - xyz[index + 1] * sz + vx,
            xyz[index] * sz + xyz[index + 1] * cz + vy,
            xyz[index + 2] + z_shift};
    }

    int flat_tail_start = n_time - 1;
    int flat_tail_segments = 0;
    if (collect_candidate_search_diagnostics_ || exact_tail_candidates_) {
        constexpr double kFlatSegmentLength2Mm2 = 1.0e-12;
        while (flat_tail_start > 0) {
            const auto& a = model[static_cast<size_t>(flat_tail_start - 1)];
            const auto& b = model[static_cast<size_t>(flat_tail_start)];
            const double dx=b[0]-a[0],dy=b[1]-a[1],dz=b[2]-a[2];
            if (dx*dx+dy*dy+dz*dz > kFlatSegmentLength2Mm2) break;
            --flat_tail_start;
        }
        flat_tail_segments = (n_time - 1) - flat_tail_start;
        if (collect_candidate_search_diagnostics_)
            ++candidate_search_counters_.objective_evaluations_with_curve;
        if (flat_tail_segments > 0) {
            if (collect_candidate_search_diagnostics_) {
                ++candidate_search_counters_.objective_evaluations_with_flat_tail;
                candidate_search_counters_.flat_tail_segments_total +=
                    static_cast<std::uint64_t>(flat_tail_segments);
                candidate_search_counters_.flat_tail_segments_max = std::max(
                    candidate_search_counters_.flat_tail_segments_max,
                    static_cast<std::uint64_t>(flat_tail_segments));
            }
        }
    }

    // The reference Spyral mesh has 2000 time samples. Searching every segment
    // for every hit at every Minuit step is the dominant cost. Sample the full
    // curve at stride 4, keep the three best separated neighborhoods, then do
    // an exact segment-distance search within +/-4 samples of each candidate.
    // This retains competing minima on curling tracks while cutting most of the
    // repeated distance evaluations.
    constexpr int kSearchStride = 4;
    constexpr size_t kCandidateCount = 3;
    double sum = 0.0;
    double sum_squared = 0.0;
    double evaluation_multiple_tail_miss_sum_mm = 0.0;
    for (const auto& d : *data_) {
        std::array<int, kCandidateCount> candidate{{-1, -1, -1}};
        std::array<double, kCandidateCount> candidate_d2;
        candidate_d2.fill(std::numeric_limits<double>::infinity());
        const auto sample = [&](int t) {
            if (t < 0 || t + 1 >= n_time) return;
            const auto& a = model[static_cast<size_t>(t)];
            const auto& b = model[static_cast<size_t>(t + 1)];
            const double d2 = PointSegmentDistance2(d.x, d.y, d.z,
                a[0], a[1], a[2], b[0], b[1], b[2]);
            for (size_t i = 0; i < kCandidateCount; ++i) {
                if (candidate[i] >= 0 && std::abs(candidate[i] - t) < 2 * kSearchStride) {
                    if (d2 < candidate_d2[i]) { candidate[i] = t; candidate_d2[i] = d2; }
                    if (sort_candidates_on_merge_) {
                        for (size_t a = 0; a < kCandidateCount; ++a) {
                            for (size_t b = a + 1; b < kCandidateCount; ++b) {
                                if (candidate_d2[b] < candidate_d2[a]) {
                                    std::swap(candidate_d2[a], candidate_d2[b]);
                                    std::swap(candidate[a], candidate[b]);
                                }
                            }
                        }
                    }
                    return;
                }
            }
            for (size_t i = 0; i < kCandidateCount; ++i) {
                if (d2 < candidate_d2[i]) {
                    for (size_t j = kCandidateCount - 1; j > i; --j) {
                        candidate[j] = candidate[j - 1]; candidate_d2[j] = candidate_d2[j - 1];
                    }
                    candidate[i] = t; candidate_d2[i] = d2; return;
                }
            }
        };
        for (int t = 0; t + 1 < n_time; t += kSearchStride) sample(t);
        sample(n_time - 2);
        double best = std::numeric_limits<double>::infinity();
        for (int center : candidate) {
            if (center < 0) continue;
            const int begin = std::max(0, center - kSearchStride);
            const int end = std::min(n_time - 2, center + kSearchStride);
            for (int t = begin; t <= end; ++t) {
                const auto& a = model[static_cast<size_t>(t)];
                const auto& b = model[static_cast<size_t>(t + 1)];
                best = std::min(best, PointSegmentDistance2(d.x, d.y, d.z,
                    a[0], a[1], a[2], b[0], b[1], b[2]));
            }
        }
        if (collect_candidate_search_diagnostics_ || exact_tail_candidates_) {
            if (collect_candidate_search_diagnostics_)
                ++candidate_search_counters_.evaluated_hits;
            unsigned int tail_centers = 0;
            if (flat_tail_segments > 0) {
                for (int center : candidate) {
                    if (center >= flat_tail_start) ++tail_centers;
                }
            }
            if (collect_candidate_search_diagnostics_) {
                if (tail_centers > 0) ++candidate_search_counters_.hits_with_tail_candidate;
                if (tail_centers > 1) ++candidate_search_counters_.hits_with_multiple_tail_candidates;
                candidate_search_counters_.tail_candidate_slots += tail_centers;
            }
            if (tail_centers > 1) {
                if (collect_candidate_search_diagnostics_)
                    ++candidate_search_counters_.multiple_tail_hits_checked_exactly;
                double exact_all_segments = std::numeric_limits<double>::infinity();
                for (int t = 0; t + 1 < n_time; ++t) {
                    const auto& a = model[static_cast<size_t>(t)];
                    const auto& b = model[static_cast<size_t>(t + 1)];
                    exact_all_segments = std::min(exact_all_segments,
                        PointSegmentDistance2(d.x, d.y, d.z,
                            a[0], a[1], a[2], b[0], b[1], b[2]));
                }
                const double missed_mm = std::sqrt(best) - std::sqrt(exact_all_segments);
                if (missed_mm > 1.0e-6) {
                    if (collect_candidate_search_diagnostics_) {
                        ++candidate_search_counters_.multiple_tail_hits_missed_exact;
                        candidate_search_counters_.multiple_tail_miss_distance_sum_mm += missed_mm;
                        candidate_search_counters_.multiple_tail_miss_distance_max_mm = std::max(
                            candidate_search_counters_.multiple_tail_miss_distance_max_mm, missed_mm);
                        evaluation_multiple_tail_miss_sum_mm += missed_mm;
                    }
                    if (exact_tail_candidates_) best = exact_all_segments;
                }
            }
        }
        if (!std::isfinite(best)) return 1.0e9;
        sum += std::sqrt(best);
        if (experimental_mean_squared_distance_objective_) sum_squared += best;
    }

    if (collect_candidate_search_diagnostics_ && evaluation_multiple_tail_miss_sum_mm > 0.0) {
        ++candidate_search_counters_.objective_evaluations_with_multiple_tail_miss;
        const double objective_delta = evaluation_multiple_tail_miss_sum_mm /
            static_cast<double>(data_->size());
        candidate_search_counters_.multiple_tail_miss_objective_delta_sum_mm += objective_delta;
        candidate_search_counters_.multiple_tail_miss_objective_delta_max_mm = std::max(
            candidate_search_counters_.multiple_tail_miss_objective_delta_max_mm, objective_delta);
    }
    if (data_->empty()) return 1.0e9;
    const double count = static_cast<double>(data_->size());
    return experimental_mean_squared_distance_objective_
        ? sum_squared / count : sum / count;
}

void TrackFitter::ResetCandidateSearchCounters() {
    candidate_search_counters_ = CandidateSearchCounters{};
}

CandidateOrderObjectiveComparison TrackFitter::EvaluateBothCandidateOrdersForDiagnostics(
        const std::vector<CloudPoint>& data, const double* par) {
    CandidateOrderObjectiveComparison result;
    if (mesh_ == nullptr || par == nullptr) {
        result.insertion_order_mm = result.sorted_order_mm = 1.0e9;
        return result;
    }
    const auto* previous_data = data_;
    const bool previous_sort = sort_candidates_on_merge_;
    const bool previous_exact_tail = exact_tail_candidates_;
    const bool previous_collect = collect_candidate_search_diagnostics_;
    const CandidateSearchCounters previous_counters = candidate_search_counters_;
    data_ = &data;
    exact_tail_candidates_ = false;
    collect_candidate_search_diagnostics_ = false;
    sort_candidates_on_merge_ = false;
    ResetCandidateSearchCounters();
    result.insertion_order_mm = Evaluate(par);
    sort_candidates_on_merge_ = true;
    ResetCandidateSearchCounters();
    result.sorted_order_mm = Evaluate(par);
    data_ = previous_data;
    sort_candidates_on_merge_ = previous_sort;
    exact_tail_candidates_ = previous_exact_tail;
    collect_candidate_search_diagnostics_ = previous_collect;
    candidate_search_counters_ = previous_counters;
    return result;
}

bool TrackFitter::Fit(const std::vector<CloudPoint>& data, const EstimateResult& guess,
                      double azimuthal, FitResult& out,
                      MinuitFitDiagnostics* diagnostics) {
    out = FitResult{};
    if (diagnostics) {
        *diagnostics = MinuitFitDiagnostics{};
        diagnostics->point_count = data.size();
    }
    if (mesh_ == nullptr) {
        if (diagnostics) diagnostics->rejection_reason = "mesh_uninitialized";
        return false;
    }
    if (!std::isfinite(experimental_vertex_radius_limit_mm_)) {
        if (diagnostics) diagnostics->rejection_reason = "nonfinite_experimental_vertex_radius_limit";
        return false;
    }
    if (!std::isfinite(experimental_vertex_phi_step_mm_) ||
        experimental_vertex_phi_step_mm_ < 0.0 ||
        (experimental_vertex_phi_step_mm_ > 0.0 && experimental_vertex_radius_limit_mm_ <= 0.0) ||
        (experimental_centered_vertex_phi_bounds_ && experimental_vertex_radius_limit_mm_ <= 0.0)) {
        if (diagnostics) diagnostics->rejection_reason = "invalid_experimental_vertex_phi_step_scale";
        return false;
    }
    if (centered_azimuth_bounds_ &&
        !std::isfinite(experimental_azimuth_bounds_center_shift_rad_)) {
        if (diagnostics) diagnostics->rejection_reason = "nonfinite_experimental_azimuth_bounds_center_shift";
        return false;
    }
    if (centered_azimuth_bounds_ && experimental_unbounded_azimuth_) {
        if (diagnostics) diagnostics->rejection_reason = "conflicting_experimental_azimuth_parameterizations";
        return false;
    }
    if (data.size() < 10) {
        if (diagnostics) diagnostics->rejection_reason = "fewer_than_10_points";
        return false;
    }
    if (!guess.ok) {
        if (diagnostics) diagnostics->rejection_reason = "invalid_estimation_seed";
        return false;
    }

    data_ = &data;
    collect_candidate_search_diagnostics_ = diagnostics != nullptr;
    ResetCandidateSearchCounters();

    // 方位角起点：估计值附近多试几个，避免最近点目标函数把拟合带进局部极小
    const double two_pi = 2.0 * Units::PI;
    const double offsets[5] = {0.0, 0.6, -0.6, 1.2, -1.2};

    double best_objective = std::numeric_limits<double>::max();
    double best_mean_distance = std::numeric_limits<double>::max();
    double best_values[6] = {0.0};
    double best_errors[6] = {0.0};
    int best_status = -1;
    double best_edm = 0.0;
    int best_start_index = -1;

    if(!std::isfinite(guess.brho)||!std::isfinite(guess.polar)||
       !std::isfinite(guess.vertex_x)||!std::isfinite(guess.vertex_y)||
       !std::isfinite(guess.vertex_z)||!std::isfinite(azimuthal)) {
        data_=nullptr;
        if (diagnostics) diagnostics->rejection_reason = "nonfinite_seed_parameter";
        return false;
    }
    const double vz_min=std::max(0.0,guess.vertex_z-60.0);
    const double vz_max=std::min(1000.0,guess.vertex_z+60.0);
    if(!(vz_max>vz_min)){
        data_=nullptr;
        if (diagnostics) diagnostics->rejection_reason = "empty_vertex_z_bounds";
        return false;
    }
    const auto& config=mesh_->Config();
    const auto brho_of_ke=[&](double ke){return std::sqrt(ke*(ke+2*mass_mev_))/(kBrhoToMomentum*charge_);};
    const double brho_min=brho_of_ke(config.ke_min),brho_max=brho_of_ke(config.ke_max);
    const double polar_min=config.polar_min*Units::DEG2RAD,polar_max=config.polar_max*Units::DEG2RAD;
    const bool use_vertex_radius_limit = experimental_vertex_radius_limit_mm_ > 0.0;
    const double vertex_radius_limit = experimental_vertex_radius_limit_mm_;
    const auto legacy_gate_score = [this,use_vertex_radius_limit](const double* par) {
        const bool search_squared = experimental_mean_squared_distance_objective_;
        experimental_mean_squared_distance_objective_ = false;
        double score = 1.0e9;
        if (use_vertex_radius_limit) {
            const double physical[6] = {par[0],par[1],par[2]*std::cos(par[3]),
                                        par[2]*std::sin(par[3]),par[4],par[5]};
            score = Evaluate(physical);
        } else {
            score = Evaluate(par);
        }
        experimental_mean_squared_distance_objective_ = search_squared;
        return score;
    };
    ROOT::Math::Functor objective([this,use_vertex_radius_limit](const double* par){
        if (!use_vertex_radius_limit) return Evaluate(par);
        const double physical[6] = {par[0],par[1],par[2]*std::cos(par[3]),
                                    par[2]*std::sin(par[3]),par[4],par[5]};
        return Evaluate(physical);
    },6);
    if (diagnostics) {
        diagnostics->attempted = true;
        diagnostics->starts.reserve(5);
    }
    for(int start=0;start<5;++start){
        double az0=std::fmod(azimuthal+offsets[start],two_pi);if(az0<0)az0+=two_pi;
        MinuitStartDiagnostics start_diagnostics;
        start_diagnostics.start_index = start;
        start_diagnostics.azimuth_offset = offsets[start];
        start_diagnostics.initial_azimuth = az0;
        const double seed_vx=std::clamp(guess.vertex_x,-40.0,40.0);
        const double seed_vy=std::clamp(guess.vertex_y,-40.0,40.0);
        const double seed_vz=std::clamp(guess.vertex_z,vz_min,vz_max);
        const double seed_rho=std::min(std::hypot(seed_vx,seed_vy),vertex_radius_limit);
        double seed_vertex_phi=std::atan2(seed_vy,seed_vx);
        if(seed_vertex_phi<0.0)seed_vertex_phi+=two_pi;
        const double vertex_phi_lower=experimental_centered_vertex_phi_bounds_
            ? seed_vertex_phi-Units::PI : 0.0;
        const double vertex_phi_upper=experimental_centered_vertex_phi_bounds_
            ? seed_vertex_phi+Units::PI : two_pi;
        const double vertex_phi_step=experimental_vertex_phi_step_mm_>0.0
            ? std::min(Units::PI/2.0,experimental_vertex_phi_step_mm_/
                std::max(seed_rho,experimental_vertex_phi_step_mm_))
            : 0.02;
        const double azimuth_lower=centered_azimuth_bounds_
            ? az0+experimental_azimuth_bounds_center_shift_rad_-Units::PI : 0.0;
        const double azimuth_upper=centered_azimuth_bounds_
            ? az0+experimental_azimuth_bounds_center_shift_rad_+Units::PI : two_pi;
        const double initial_values[6] = {
            std::clamp(guess.brho,brho_min,brho_max),
            std::clamp(guess.polar,polar_min,polar_max),
            use_vertex_radius_limit?seed_rho:seed_vx,
            use_vertex_radius_limit?seed_vertex_phi:seed_vy,
            seed_vz,
            az0};
        double lower_bounds[6] = {brho_min,polar_min,
            use_vertex_radius_limit?-vertex_radius_limit:-40.0,
            use_vertex_radius_limit?-vertex_radius_limit:-40.0,vz_min,0.0};
        const double upper_bounds[6] = {brho_max,polar_max,
            use_vertex_radius_limit?vertex_radius_limit:40.0,
            use_vertex_radius_limit?vertex_radius_limit:40.0,vz_max,azimuth_upper};
        lower_bounds[5]=experimental_unbounded_azimuth_
            ? std::numeric_limits<double>::quiet_NaN() : azimuth_lower;
        double diagnostic_upper_bounds[6];
        std::copy(std::begin(upper_bounds),std::end(upper_bounds),diagnostic_upper_bounds);
        if(experimental_unbounded_azimuth_)
            diagnostic_upper_bounds[5]=std::numeric_limits<double>::quiet_NaN();
        double diagnostic_initial_values[6];
        std::copy(std::begin(initial_values),std::end(initial_values),diagnostic_initial_values);
        if(use_vertex_radius_limit){
            diagnostic_initial_values[2]=seed_rho*std::cos(seed_vertex_phi);
            diagnostic_initial_values[3]=seed_rho*std::sin(seed_vertex_phi);
        }
        for(int i=0;i<6;++i) {
            start_diagnostics.initial_values[i]=diagnostic_initial_values[i];
            start_diagnostics.lower_bounds[i]=lower_bounds[i];
            start_diagnostics.upper_bounds[i]=diagnostic_upper_bounds[i];
        }
        ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
        minimizer.SetPrintLevel(-1);minimizer.SetMaxFunctionCalls(5000);
        minimizer.SetTolerance(.1);minimizer.SetErrorDef(1);minimizer.SetFunction(objective);
        minimizer.SetLimitedVariable(0,"brho",std::clamp(guess.brho,brho_min,brho_max),.01,brho_min,brho_max);
        minimizer.SetLimitedVariable(1,"polar",std::clamp(guess.polar,polar_min,polar_max),.01,polar_min,polar_max);
        if(use_vertex_radius_limit){
            minimizer.SetLimitedVariable(2,"vertex_rho",seed_rho,1.,0.,vertex_radius_limit);
            minimizer.SetLimitedVariable(3,"vertex_phi",seed_vertex_phi,vertex_phi_step,
                                         vertex_phi_lower,vertex_phi_upper);
        }else{
            minimizer.SetLimitedVariable(2,"vx",seed_vx,1.,-40.,40.);
            minimizer.SetLimitedVariable(3,"vy",seed_vy,1.,-40.,40.);
        }
        minimizer.SetLimitedVariable(4,"vz",std::clamp(guess.vertex_z,vz_min,vz_max),1.,vz_min,vz_max);
        if(experimental_unbounded_azimuth_)
            minimizer.SetVariable(5,"azimuth",az0,.02);
        else
            minimizer.SetLimitedVariable(5,"azimuth",az0,.02,azimuth_lower,azimuth_upper);
        minimizer.FixVariable(2);minimizer.FixVariable(3);minimizer.FixVariable(4);
        ResetCandidateSearchCounters();
        const bool fixed_return = minimizer.Minimize();
        start_diagnostics.fixed_vertex.minimize_return = fixed_return;
        start_diagnostics.fixed_vertex.status = minimizer.Status();
        start_diagnostics.fixed_vertex.covariance_status = minimizer.CovMatrixStatus();
        start_diagnostics.fixed_vertex.function_calls = minimizer.NCalls();
        start_diagnostics.fixed_vertex.objective = minimizer.MinValue();
        start_diagnostics.fixed_vertex.edm = minimizer.Edm();
        const double* fixed_values=minimizer.X();
        const double* fixed_errors=minimizer.Errors();
        if(fixed_values){
            for(int i=0;i<6;++i){
                start_diagnostics.fixed_vertex.final_values[i]=fixed_values[i];
                start_diagnostics.fixed_vertex.final_errors[i]=fixed_errors?fixed_errors[i]:0.0;
            }
            if(centered_azimuth_bounds_ || experimental_unbounded_azimuth_){
                start_diagnostics.fixed_vertex.final_values[5]=
                    std::fmod(start_diagnostics.fixed_vertex.final_values[5],two_pi);
                if(start_diagnostics.fixed_vertex.final_values[5]<0.0)
                    start_diagnostics.fixed_vertex.final_values[5]+=two_pi;
            }
            if(use_vertex_radius_limit){
                const double rho=fixed_values[2],phi=fixed_values[3];
                start_diagnostics.fixed_vertex.final_values[2]=rho*std::cos(phi);
                start_diagnostics.fixed_vertex.final_values[3]=rho*std::sin(phi);
                const double sigma_rho=fixed_errors?fixed_errors[2]:0.0;
                const double sigma_phi=fixed_errors?fixed_errors[3]:0.0;
                start_diagnostics.fixed_vertex.final_errors[2]=
                    std::hypot(std::cos(phi)*sigma_rho,rho*std::sin(phi)*sigma_phi);
                start_diagnostics.fixed_vertex.final_errors[3]=
                    std::hypot(std::sin(phi)*sigma_rho,rho*std::cos(phi)*sigma_phi);
            }
        }
        start_diagnostics.fixed_vertex.candidate_search = candidate_search_counters_;
        if(experimental_fresh_free_vertex_minuit_ && diagnostics && fixed_values){
            start_diagnostics.fresh_free_vertex_ran=true;
            ROOT::Minuit2::Minuit2Minimizer fresh_minimizer(ROOT::Minuit2::kMigrad);
            fresh_minimizer.SetPrintLevel(-1);fresh_minimizer.SetMaxFunctionCalls(5000);
            fresh_minimizer.SetTolerance(.1);fresh_minimizer.SetErrorDef(1);fresh_minimizer.SetFunction(objective);
            fresh_minimizer.SetLimitedVariable(0,"brho",fixed_values[0],.01,brho_min,brho_max);
            fresh_minimizer.SetLimitedVariable(1,"polar",fixed_values[1],.01,polar_min,polar_max);
            if(use_vertex_radius_limit){
                fresh_minimizer.SetLimitedVariable(2,"vertex_rho",fixed_values[2],1.,0.,vertex_radius_limit);
                fresh_minimizer.SetLimitedVariable(3,"vertex_phi",fixed_values[3],vertex_phi_step,
                                                   vertex_phi_lower,vertex_phi_upper);
            }else{
                fresh_minimizer.SetLimitedVariable(2,"vx",fixed_values[2],1.,-40.,40.);
                fresh_minimizer.SetLimitedVariable(3,"vy",fixed_values[3],1.,-40.,40.);
            }
            fresh_minimizer.SetLimitedVariable(4,"vz",fixed_values[4],1.,vz_min,vz_max);
            if(experimental_unbounded_azimuth_)
                fresh_minimizer.SetVariable(5,"azimuth",fixed_values[5],.02);
            else
                fresh_minimizer.SetLimitedVariable(5,"azimuth",fixed_values[5],.02,azimuth_lower,azimuth_upper);
            ResetCandidateSearchCounters();
            start_diagnostics.fresh_free_vertex.minimize_return=fresh_minimizer.Minimize();
            start_diagnostics.fresh_free_vertex.status=fresh_minimizer.Status();
            start_diagnostics.fresh_free_vertex.covariance_status=fresh_minimizer.CovMatrixStatus();
            start_diagnostics.fresh_free_vertex.function_calls=fresh_minimizer.NCalls();
            start_diagnostics.fresh_free_vertex.objective=fresh_minimizer.MinValue();
            start_diagnostics.fresh_free_vertex.edm=fresh_minimizer.Edm();
            start_diagnostics.fresh_free_vertex.candidate_search=candidate_search_counters_;
            const double* fresh_values=fresh_minimizer.X();
            const double* fresh_errors=fresh_minimizer.Errors();
            if(fresh_values){
                for(int i=0;i<6;++i){
                    start_diagnostics.fresh_free_vertex.final_values[i]=fresh_values[i];
                    start_diagnostics.fresh_free_vertex.final_errors[i]=fresh_errors?fresh_errors[i]:0.0;
                }
                if(centered_azimuth_bounds_ || experimental_unbounded_azimuth_){
                    start_diagnostics.fresh_free_vertex.final_values[5]=
                        std::fmod(start_diagnostics.fresh_free_vertex.final_values[5],two_pi);
                    if(start_diagnostics.fresh_free_vertex.final_values[5]<0.0)
                        start_diagnostics.fresh_free_vertex.final_values[5]+=two_pi;
                }
                if(use_vertex_radius_limit){
                    const double rho=fresh_values[2],phi=fresh_values[3];
                    start_diagnostics.fresh_free_vertex.final_values[2]=rho*std::cos(phi);
                    start_diagnostics.fresh_free_vertex.final_values[3]=rho*std::sin(phi);
                    const double sigma_rho=fresh_errors?fresh_errors[2]:0.0;
                    const double sigma_phi=fresh_errors?fresh_errors[3]:0.0;
                    start_diagnostics.fresh_free_vertex.final_errors[2]=
                        std::hypot(std::cos(phi)*sigma_rho,rho*std::sin(phi)*sigma_phi);
                    start_diagnostics.fresh_free_vertex.final_errors[3]=
                        std::hypot(std::sin(phi)*sigma_rho,rho*std::cos(phi)*sigma_phi);
                }
            }
        }
        minimizer.ReleaseVariable(2);minimizer.ReleaseVariable(3);minimizer.ReleaseVariable(4);
        ResetCandidateSearchCounters();
        const bool free_return = minimizer.Minimize();
        start_diagnostics.free_vertex.minimize_return = free_return;
        start_diagnostics.free_vertex.status = minimizer.Status();
        start_diagnostics.free_vertex.covariance_status = minimizer.CovMatrixStatus();
        start_diagnostics.free_vertex.function_calls = minimizer.NCalls();
        start_diagnostics.free_vertex.objective = minimizer.MinValue();
        start_diagnostics.free_vertex.edm = minimizer.Edm();
        start_diagnostics.free_vertex.candidate_search = candidate_search_counters_;
        const double optimizer_objective=minimizer.MinValue();
        const double* attempt_values=minimizer.X();
        const double* attempt_errors=minimizer.Errors();
        const double mean_distance_cost = attempt_values
            ? (experimental_mean_squared_distance_objective_
                ? legacy_gate_score(attempt_values) : optimizer_objective)
            : std::numeric_limits<double>::infinity();
        start_diagnostics.mean_distance_cost = mean_distance_cost;
        if (attempt_values) {
            for(int i=0;i<6;++i) {
                start_diagnostics.final_values[i]=attempt_values[i];
                start_diagnostics.final_errors[i]=attempt_errors?attempt_errors[i]:0.0;
            }
            if(centered_azimuth_bounds_ || experimental_unbounded_azimuth_){
                start_diagnostics.final_values[5]=std::fmod(start_diagnostics.final_values[5],two_pi);
                if(start_diagnostics.final_values[5]<0.0)start_diagnostics.final_values[5]+=two_pi;
            }
            if(use_vertex_radius_limit){
                const double rho=attempt_values[2],phi=attempt_values[3];
                start_diagnostics.final_values[2]=rho*std::cos(phi);
                start_diagnostics.final_values[3]=rho*std::sin(phi);
                const double sigma_rho=attempt_errors?attempt_errors[2]:0.0;
                const double sigma_phi=attempt_errors?attempt_errors[3]:0.0;
                start_diagnostics.final_errors[2]=std::hypot(std::cos(phi)*sigma_rho,
                                                             rho*std::sin(phi)*sigma_phi);
                start_diagnostics.final_errors[3]=std::hypot(std::sin(phi)*sigma_rho,
                                                             rho*std::cos(phi)*sigma_phi);
            }
        }
        start_diagnostics.below_existing_cost_gate = std::isfinite(mean_distance_cost) && mean_distance_cost < 20.0;
        if(std::isfinite(optimizer_objective)&&optimizer_objective<best_objective){
            best_objective=optimizer_objective;
            best_mean_distance=mean_distance_cost;
            best_status=minimizer.Status();best_edm=minimizer.Edm();
            const double* values=minimizer.X();const double* errors=minimizer.Errors();
            for(int i=0;i<6;++i){best_values[i]=values[i];best_errors[i]=errors?errors[i]:0.;}
            if(centered_azimuth_bounds_ || experimental_unbounded_azimuth_){
                best_values[5]=std::fmod(best_values[5],two_pi);
                if(best_values[5]<0.0)best_values[5]+=two_pi;
            }
            best_start_index = start;
        }
        if (diagnostics) diagnostics->starts.push_back(std::move(start_diagnostics));
    }

    out.brho = best_values[0];
    out.sigma_brho = best_errors[0];
    out.polar = best_values[1];
    out.sigma_polar = best_errors[1];
    out.azimuthal = best_values[5];
    if(use_vertex_radius_limit){
        out.vx=best_values[2]*std::cos(best_values[3]);
        out.vy=best_values[2]*std::sin(best_values[3]);
    }else{
        out.vx=best_values[2];
        out.vy=best_values[3];
    }
    out.vz = best_values[4];
    out.cost = best_mean_distance;
    out.optimizer_objective = best_objective;
    out.edm = best_edm;
    out.status = best_status;
    // The legacy acceptance gate is cost-only. It does not certify Migrad's
    // numerical convergence or the physical correctness of the selected track.
    out.ok = PassesExistingMinuitGate(out.brho, best_mean_distance);

    if (diagnostics) {
        diagnostics->accepted = out.ok;
        diagnostics->final_cost = best_mean_distance;
        diagnostics->final_optimizer_objective = best_objective;
        diagnostics->best_start_index = best_start_index;
        if (best_start_index >= 0 && static_cast<size_t>(best_start_index) < diagnostics->starts.size())
            diagnostics->starts[static_cast<size_t>(best_start_index)].selected_as_best = true;
        if (!std::isfinite(best_objective)) diagnostics->rejection_reason = "no_finite_minuit_objective";
        else if (!std::isfinite(out.brho) || out.brho <= 0.0) diagnostics->rejection_reason = "invalid_final_brho";
        else if (!(best_mean_distance < 20.0)) diagnostics->rejection_reason = "mean_distance_not_below_20mm";
        else diagnostics->rejection_reason = "accepted_by_existing_gate";
    }

    data_ = nullptr;
    collect_candidate_search_diagnostics_ = false;
    ResetCandidateSearchCounters();
    return out.ok;
}

bool TrackFitter::FitPythonEquivalent(const std::vector<CloudPoint>& data,
                                     const EstimateResult& guess,double z_error_m,
                                     double beam_radius_m,FitResult& out){
    out=FitResult{};
    if(mesh_==nullptr||!mesh_->IsPythonMesh()||data.empty()||!guess.ok||
       !(z_error_m>0)||!(beam_radius_m>0))return false;
    const auto& grid=mesh_->Config();
    const double polar_min=grid.polar_min*Units::DEG2RAD;
    const double polar_max=grid.polar_max*Units::DEG2RAD;
    const double brho_min=std::max(std::sqrt(grid.ke_min*(grid.ke_min+2*mass_mev_))/
                                   (kBrhoToMomentum*charge_),guess.brho-2.0);
    const double brho_max=std::min(std::sqrt(grid.ke_max*(grid.ke_max+2*mass_mev_))/
                                   (kBrhoToMomentum*charge_),guess.brho+2.0);
    double pmin=polar_min,pmax=polar_max;
    if(guess.polar>Units::PI*.5&&pmin<Units::PI*.5)pmin=Units::PI*.5;
    else if(guess.polar<Units::PI*.5&&pmax>Units::PI*.5)pmax=Units::PI*.5;
    const double zmin=std::max(0.0,guess.vertex_z*.001-.2);
    const double zmax=std::min(grid.z_ref*2.,guess.vertex_z*.001+.2);
    if(!(brho_max>brho_min&&pmax>pmin&&zmax>zmin)||guess.brho<brho_min||guess.brho>brho_max||
       guess.polar<pmin||guess.polar>pmax)return false;
    const double gx=guess.vertex_x*.001,gy=guess.vertex_y*.001;
    const double rho=std::hypot(gx,gy);
    const double phi=std::atan2(gy,gx)<0?std::atan2(gy,gx)+2*Units::PI:std::atan2(gy,gx);
    if(rho>beam_radius_m)return false;
    std::vector<double> x,y,z,w;x.reserve(data.size());y.reserve(data.size());z.reserve(data.size());w.reserve(data.size());
    constexpr double big_pad_height_m=.01;
    for(const auto& point:data){
        const double scale=point.scale;
        const double xe=scale*big_pad_height_m*.5,ye=scale*big_pad_height_m/std::sqrt(3.0);
        const double variance=xe*xe+ye*ye+z_error_m*z_error_m;
        if(!(variance>0)||!std::isfinite(variance)||!std::isfinite(scale))return false;
        x.push_back(point.x*.001);y.push_back(point.y*.001);z.push_back(point.z*.001);w.push_back(1.0/variance);
    }
    const double sum_weights=std::accumulate(w.begin(),w.end(),0.0);
    const auto objective=[&](const double* par){
        const double brho=par[0],polar=par[1],vertex_rho=par[2],vertex_phi=par[3],vertex_z=par[4],azimuth=par[5];
        const double momentum=kBrhoToMomentum*brho*charge_;
        const double ke=std::hypot(momentum,mass_mev_)-mass_mev_;
        std::vector<double> trajectory;
        if(!(ke>=grid.ke_min&&ke<=grid.ke_max)||!mesh_->InterpolatePython(ke,polar,trajectory))return 1.0e6;
        const double vx=vertex_rho*std::cos(vertex_phi),vy=vertex_rho*std::sin(vertex_phi);
        const double ca=std::cos(azimuth),sa=std::sin(azimuth);double error=0;
        const size_t nt=static_cast<size_t>(grid.n_time);
        for(size_t i=0;i<x.size();++i){double best=std::numeric_limits<double>::infinity();
            for(size_t t=0;t<nt;++t){const double* q=&trajectory[t*3];
                const double tx=q[0]*ca-q[1]*sa+vx-x[i];
                const double ty=q[0]*sa+q[1]*ca+vy-y[i];
                const double tz=q[2]+vertex_z-z[i];
                const double d=std::sqrt(tx*tx+ty*ty+tz*tz);if(d<best)best=d;
            }
            error+=best*w[i];
        }
        return error/sum_weights;
    };
    const double two_pi=2.0*Units::PI;
    ROOT::Math::Functor functor(objective,6);
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetPrintLevel(-1);minimizer.SetMaxFunctionCalls(24000);minimizer.SetTolerance(.01);
    minimizer.SetErrorDef(1);minimizer.SetFunction(functor);
    minimizer.SetLimitedVariable(0,"brho",guess.brho,.01,brho_min,brho_max);
    minimizer.SetLimitedVariable(1,"polar",guess.polar,.01,pmin,pmax);
    minimizer.SetLimitedVariable(2,"vertex_rho",rho,.001,0.,beam_radius_m);
    minimizer.SetLimitedVariable(3,"vertex_phi",phi,.01,0.,two_pi);
    minimizer.SetLimitedVariable(4,"vertex_z",std::clamp(guess.vertex_z*.001,zmin,zmax),.01,zmin,zmax);
    minimizer.SetLimitedVariable(5,"azimuthal",std::clamp(guess.azimuthal,0.,two_pi),.01,0.,two_pi);
    minimizer.Minimize();
    const double* value=minimizer.X();
    if(!value||!std::isfinite(minimizer.MinValue())||!std::isfinite(value[0]))return false;
    out.brho=value[0];out.polar=value[1];out.vx=value[2]*std::cos(value[3])*1000.;
    out.vy=value[2]*std::sin(value[3])*1000.;out.vz=value[4]*1000.;out.azimuthal=value[5];
    out.sigma_brho=out.sigma_polar=1.0e6;out.cost=minimizer.MinValue();out.status=minimizer.Status();
    out.ok=out.brho>0&&std::isfinite(out.polar)&&std::isfinite(out.azimuthal)&&
           std::isfinite(out.vx)&&std::isfinite(out.vy)&&std::isfinite(out.vz);
    return out.ok;
}

} // namespace rk4
