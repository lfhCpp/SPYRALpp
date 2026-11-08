#include "RecoConfig.h"
//----------------------------------------------------------------------
// rk4_analyze：点云 -> 聚类 -> 估计 -> 拟合 -> 结果 TTree
//
// 用法：rk4_analyze <clouds.root> <out_fit.root> <config>
//----------------------------------------------------------------------

#include "core/MassData.h"
#include "core/Units.h"
#include "io/CloudIO.h"
#include "io/Config.h"
#include "io/ResultIO.h"
#include "io/PidStateIO.h"
#include "io/StageIO.h"
#include "reco/Clusterer.h"
#include "reco/Estimator.h"
#include "reco/ParticleId.h"
#include "reco/KalmanFitter.h"
#include "reco/TrackFitter.h"
#include "reco/TrackMesh.h"
#include "sim/DetectorConfig.h"
#include "sim/TargetDensity.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <exception>
#include <stdexcept>
#include <vector>

using namespace rk4;

int main(int argc, char** argv) {
    if (argc < 4) {
        std::printf("用法: rk4_analyze <clouds.root> <out_fit.root> <config> "
                    "[shard_index] [shard_count]\n");
        std::printf("  分片模式：只处理 index %% shard_count == shard_index 的点云，\n");
        std::printf("  多个进程可各自映射同一网格（solver.mesh_shared = 1）后并行跑不同分片。\n");
        return 1;
    }

    const int shard_index = (argc > 4) ? std::atoi(argv[4]) : 0;
    const int shard_count = (argc > 5) ? std::atoi(argv[5]) : 1;
    if (shard_count < 1 || shard_index < 0 || shard_index >= shard_count) {
        std::fprintf(stderr, "非法分片参数：index=%d count=%d\n", shard_index, shard_count);
        return 1;
    }

    Config cfg;
    if (!cfg.Load(argv[3])) {
        std::fprintf(stderr, "无法读取配置 %s\n", argv[3]);
        return 1;
    }

    const std::string estimate_method = cfg.GetString("estimate.method", "reinsch");
    const double reinsch_lambda = cfg.GetDouble("estimate.reinsch_lambda", 100.0);
    if (estimate_method != "legacy" && estimate_method != "reinsch") {
        std::fprintf(stderr, "Unknown estimate.method: %s (expected legacy or reinsch)\n", estimate_method.c_str());
        return 1;
    }
    if (estimate_method == "reinsch" && (!std::isfinite(reinsch_lambda) || reinsch_lambda < 0)) {
        std::fprintf(stderr, "estimate.reinsch_lambda must be finite and nonnegative\n");
        return 1;
    }

    std::vector<PointCloud> clouds;
    if (!ReadPointClouds(argv[1], clouds)) {
        std::fprintf(stderr, "无法读取点云 %s\n", argv[1]);
        return 1;
    }

    std::vector<Element> comp;
    Element element;
    element.a = cfg.GetInt("gas.a", 1);
    element.z = cfg.GetInt("gas.z", 1);
    element.s = cfg.GetDouble("gas.stoich", 2);
    comp.push_back(element);
    GasTarget gas(comp, cfg.GetDouble("gas.pressure_torr", 600.0));

    const int particle_A = cfg.GetInt("solver.a", 1);
    const int particle_Z = cfg.GetInt("solver.z", 1);
    const double mass = MassData::GetNuclearMassNuclear(particle_A, particle_Z) * Units::u_MeV;

    MeshConfig mc;
    mc.ke_bins = cfg.GetInt("solver.ke_bins", 256);
    mc.log_ke = cfg.GetBool("solver.log_ke", true);
    mc.ke_min = cfg.GetDouble("solver.ke_min", 0.1);
    mc.ke_max = cfg.GetDouble("solver.ke_max", 70.0);
    mc.polar_bins = cfg.GetInt("solver.polar_bins", 129);
    mc.polar_min = cfg.GetDouble("solver.polar_min", 2.0);
    mc.polar_max = cfg.GetDouble("solver.polar_max", 88.0);
    mc.n_time = cfg.GetInt("solver.n_time", 400);
    mc.dt = cfg.GetDouble("solver.dt", 3.0e-10);
    mc.z_ref = cfg.GetDouble("solver.z_ref", 0.5);

    const double bfield = cfg.GetDouble("detector.bfield_t", 2.85);
    const double efield = cfg.GetDouble("detector.efield_vm", 55000.0);

    // 可选后端：solver.fitter = minuit（默认）| kalman
    const std::string fitter_name = cfg.GetString("solver.fitter", "minuit");
    const bool use_kalman = (fitter_name == "kalman" || fitter_name == "ukf");
    const bool sort_minuit_candidates = cfg.GetBool(
        "solver.minuit_sort_candidates",
        cfg.GetBool("solver.experimental_minuit_sort_candidates", true));
    KalmanConfig kalman_cfg;
    kalman_cfg.meas_sigma_mm = cfg.GetDouble("solver.kalman_meas_sigma", kalman_cfg.meas_sigma_mm);
    kalman_cfg.per_cluster_covariance = cfg.GetBool(
        "solver.kalman_per_cluster_covariance", kalman_cfg.per_cluster_covariance);
    kalman_cfg.cluster_coef_t_cm2_us = cfg.GetDouble(
        "solver.kalman_cluster_coef_t", kalman_cfg.cluster_coef_t_cm2_us);
    kalman_cfg.cluster_coef_l_cm2_us = cfg.GetDouble(
        "solver.kalman_cluster_coef_l", kalman_cfg.cluster_coef_l_cm2_us);
    kalman_cfg.cluster_drift_velocity_cm_us = cfg.GetDouble(
        "solver.kalman_cluster_drift_velocity", kalman_cfg.cluster_drift_velocity_cm_us);
    kalman_cfg.cluster_tb_time_us = cfg.GetDouble(
        "solver.kalman_cluster_tb_time", kalman_cfg.cluster_tb_time_us);
    kalman_cfg.cluster_pad_resolution_xy_mm = cfg.GetDouble(
        "solver.kalman_cluster_pad_resolution_xy", kalman_cfg.cluster_pad_resolution_xy_mm);
    kalman_cfg.mom_sigma_frac =
        cfg.GetDouble("solver.kalman_mom_sigma_frac", kalman_cfg.mom_sigma_frac);
    kalman_cfg.pos_sigma_mm = cfg.GetDouble("solver.kalman_pos_sigma", kalman_cfg.pos_sigma_mm);
    kalman_cfg.ang_sigma_deg = cfg.GetDouble("solver.kalman_ang_sigma", kalman_cfg.ang_sigma_deg);
    kalman_cfg.min_clusters = cfg.GetInt("solver.kalman_min_clusters", kalman_cfg.min_clusters);
    kalman_cfg.max_steps = cfg.GetInt("solver.kalman_max_steps", kalman_cfg.max_steps);
    kalman_cfg.gap_max_steps = cfg.GetInt("solver.kalman_gap_max_steps", kalman_cfg.gap_max_steps);
    kalman_cfg.alpha = cfg.GetDouble("solver.kalman_alpha", kalman_cfg.alpha);
    kalman_cfg.beta = cfg.GetDouble("solver.kalman_beta", kalman_cfg.beta);
    kalman_cfg.kappa = cfg.GetDouble("solver.kalman_kappa", kalman_cfg.kappa);
    kalman_cfg.iterations = cfg.GetInt("solver.kalman_iterations", kalman_cfg.iterations);
    kalman_cfg.ref_covariance_inflation = cfg.GetDouble(
        "solver.kalman_ref_covariance_inflation", kalman_cfg.ref_covariance_inflation);
    kalman_cfg.smooth = cfg.GetBool("solver.kalman_smooth", kalman_cfg.smooth);
    kalman_cfg.straggling = cfg.GetBool("solver.kalman_straggling", false);
    kalman_cfg.mult_scattering = cfg.GetBool("solver.kalman_mult_scattering", false);
    kalman_cfg.rad_length_mm =
        cfg.GetDouble("solver.kalman_rad_length_mm", kalman_cfg.rad_length_mm);
    kalman_cfg.back_extrapolate =
        cfg.GetBool("solver.kalman_back_extrap", kalman_cfg.back_extrapolate);
    kalman_cfg.back_extrap_max_mm =
        cfg.GetDouble("solver.kalman_back_extrap_max", kalman_cfg.back_extrap_max_mm);
    kalman_cfg.chi2_max = cfg.GetDouble("solver.kalman_chi2_max", kalman_cfg.chi2_max);
    kalman_cfg.innovation_chi2_max = cfg.GetDouble("solver.kalman_innovation_chi2_max", kalman_cfg.innovation_chi2_max);
    kalman_cfg.min_measurement_fraction = cfg.GetDouble("solver.kalman_min_measurement_fraction", kalman_cfg.min_measurement_fraction);
    kalman_cfg.step_mm = cfg.GetDouble("solver.kalman_step_mm", kalman_cfg.step_mm);
    kalman_cfg.z_pad_plane_mm = cfg.GetDouble("solver.kalman_z_pad_plane_mm", kalman_cfg.z_pad_plane_mm);
    if (use_kalman) {
        std::printf("拟合后端：kalman（UKF + RTS 平滑，测点 σ=%.2f mm）\n",
                    kalman_cfg.meas_sigma_mm);
    } else {
        std::printf("拟合后端：minuit（轨迹网格插值）\n");
    }

    // 网格缓存：参数一致时直接读盘，避免每次重建（生成一次约数秒）
    // Kalman 后端逐点传播，不需要轨迹网格（也就不必花几秒钟建网格/读缓存）
    TrackMesh mesh;
    const std::string cache_path = cfg.GetString("solver.mesh_cache", "");

    if (!use_kalman) {
        bool mesh_ready = false;
        const bool use_shared = cfg.GetBool("solver.mesh_shared", false);
        if (!cache_path.empty() && use_shared && mesh.LoadShared(cache_path) &&
            mesh.MatchesConfig(mc, particle_A, particle_Z, gas, bfield, efield)) {
            std::printf("共享映射轨迹网格：%s（零拷贝）\n", cache_path.c_str());
            mesh_ready = true;
        }

        if (!mesh_ready && !cache_path.empty() && mesh.Load(cache_path) &&
            mesh.MatchesConfig(mc, particle_A, particle_Z, gas, bfield, efield)) {
            std::printf("复用轨迹网格缓存：%s\n", cache_path.c_str());
            mesh_ready = true;
        }

        if (!mesh_ready) {
            if (!mesh.Build(mc, particle_A, particle_Z, gas, bfield, efield)) {
                std::fprintf(stderr, "轨迹网格生成失败\n");
                return 1;
            }
            if (!cache_path.empty()) {
                if (mesh.Save(cache_path)) {
                    std::printf("轨迹网格已缓存：%s\n", cache_path.c_str());
                } else {
                    std::fprintf(stderr, "网格缓存写入失败（继续运行）\n");
                }
            }
        }
    }

    TrackFitter fitter;
    if (!use_kalman && !fitter.Init(&mesh, mass, particle_Z)) {
        std::fprintf(stderr, "拟合器初始化失败\n");
        return 1;
    }
    fitter.SetCenteredAzimuthBounds(cfg.GetBool(
        "solver.minuit_centered_azimuth_bounds",
        cfg.GetBool("solver.experimental_minuit_centered_azimuth_bounds", false)));
    fitter.SetSortCandidates(sort_minuit_candidates);


    const ClusterParams cluster_params=ReadClusterParams(cfg);

    const bool report = cfg.GetBool("cluster.report", false);

    // PID 预筛选：先用粗估的 (Bρ, dE/dx) 判粒子，把明显不是目标粒子的簇挡在拟合之前。
    // 粗估精度只有 ~10%，但足以挡掉大部分本底/噪声簇。
    const bool pid_gate = cfg.GetBool("pid.gate", false);
    const double pid_min_sep = cfg.GetDouble("pid.min_sep", 0.0);
    const int pid_require = cfg.GetInt("pid.require", -1); // -1 = 不限定种类
    // A scalar pid.scale is not a run-linked calibration record. Keep the
    // identity state unknown until this executable consumes verified provenance.
    const bool pid_calibration_available = false;
    const std::string pid_state_csv = cfg.GetString("output.pid_state_csv", "");

    // 可选：中间过程文件（cluster / estimate 两个阶段）
    const bool dump_stages = cfg.GetBool("output.stages", false);

    // ---- 处理一个点云的完整流程 ----
    auto process_cloud = [&](const PointCloud& cloud, Clusterer& clusterer,
                             Estimator& estimator, ParticleId& pid,
                             TrackFitter& fitter, KalmanFitter& kfitter, double pid_scale,
                             std::vector<FitRecord>& out,
                             std::vector<ClusterRecord>& clusters_out,
                             std::vector<EstimateRecord>& estimates_out,
                             std::vector<PidCandidateAuditRow>& pid_audit_out) {
        std::vector<Cluster> clusters;
        clusterer.Run(cloud, clusters);

        if (report) {
            std::fprintf(stderr,
                         "[event %lld] pts=%zu NN_median=%.1f mm eps=%.1f mm clusters=%zu "
                         "joins=%d backend=%s\n",
                         cloud.event, cloud.points.size(), clusterer.LastMedianNN(),
                         clusterer.LastEps(), clusters.size(), clusterer.LastJoinCount(),
                         clusterer.LastBackendHdbscan() ? "hdbscan" : "dbscan");
        }

        for (size_t cluster_index = 0; cluster_index < clusters.size(); ++cluster_index) {
            const auto& cluster = clusters[cluster_index];
            PidCandidateAuditRow audit;
            audit.event = cloud.event;
            audit.cluster_index = static_cast<int>(cluster_index);
            audit.cluster_label = cluster.label;
            audit.n_points = static_cast<int>(cluster.points.size());
            audit.ic_amplitude = cloud.metadata.ic_amplitude;
            audit.ic_state = ClassifyIcObservation(cloud.metadata.ic_amplitude);
            audit.pid_gate_configured = pid_gate ? 1 : 0;
            audit.pid_required_index = pid_require;
            audit.pid_min_separation = pid_min_sep;
            audit.pid_calibration_available = pid_calibration_available ? 1 : 0;
            const auto set_pid_state = [&](const PidResult& result, bool evaluated,
                                           double brho, double dedx) {
                audit.pid_classifier_evaluated = evaluated ? 1 : 0;
                audit.pid_best_index = evaluated ? result.index : -1;
                audit.pid_separation = evaluated ? result.separation : 0.0;
                PidStateInput input;
                input.gate_configured = pid_gate;
                input.target_species_configured =
                    pid_require >= 0 && static_cast<size_t>(pid_require) < pid.Candidates().size();
                input.observation_present = brho != 0.0 && dedx != 0.0;
                input.observation_valid = std::isfinite(brho) && std::isfinite(dedx) &&
                                          brho > 0.0 && dedx > 0.0;
                input.calibration_available = pid_calibration_available;
                input.classifier_evaluated = evaluated;
                input.target_matches = evaluated && result.index == pid_require;
                input.separation_passes = evaluated && result.separation >= pid_min_sep;
                const auto state = ClassifyPidState(input);
                audit.pid_state = state.state;
                audit.pid_unknown_reason = state.unknown_reason;
            };
            audit.fit_rejection_stage = "estimation";

            if (dump_stages) {
                ClusterRecord cr;
                cr.event = cloud.event;
                cr.cluster = cluster.label;
                cr.npoints = static_cast<int>(cluster.points.size());
                cr.direction = static_cast<int>(cluster.direction);
                cr.x.reserve(cluster.points.size());
                cr.y.reserve(cluster.points.size());
                cr.z.reserve(cluster.points.size());
                cr.amplitude.reserve(cluster.points.size());
                cr.integral.reserve(cluster.points.size());
                cr.tb.reserve(cluster.points.size());
                for (const auto& cpt : cluster.points) {
                    cr.x.push_back(cpt.x);
                    cr.y.push_back(cpt.y);
                    cr.z.push_back(cpt.z);
                    cr.amplitude.push_back(cpt.amplitude);
                    cr.integral.push_back(cpt.integral);
                    cr.tb.push_back(cpt.tb);
                }
                clusters_out.push_back(std::move(cr));
            }

            EstimateResult est;
            const bool est_ok = estimator.Estimate(cluster, bfield, est);
            audit.estimate_ok = est_ok ? 1 : 0;
            if (dump_stages) {
                EstimateRecord er;
                er.event = cloud.event;
                er.cluster = cluster.label;
                er.ok = est_ok ? 1 : 0;
                er.vx = est.vertex_x;
                er.vy = est.vertex_y;
                er.vz = est.vertex_z;
                er.polar = est.polar;
                er.azimuthal = est.azimuthal;
                er.brho = est.brho;
                er.radius = est.radius;
                er.dEdx = est.dEdx;
                er.arclength = est.arclength;
                er.charge = est.charge;
                er.direction = est.direction;
                estimates_out.push_back(er);
            }
            if (!est_ok) {
                PidResult empty;
                set_pid_state(empty, false, est.brho, est.dEdx);
                if (!pid_state_csv.empty()) pid_audit_out.push_back(std::move(audit));
                continue;
            }

            PidResult pre;
            const bool pre_identified = pid_gate && pid.Identify(est.brho, est.dEdx, pre);
            set_pid_state(pre, pre_identified, est.brho, est.dEdx);
            if (pid_gate) {
                const bool pass_pre_gate = PassesExistingPidGate(
                    pid_gate, pre_identified, pre.separation, pid_min_sep, pre.index, pid_require);
                audit.pid_gate_pass = pass_pre_gate ? 1 : 0;
                if (!pass_pre_gate) {
                    audit.fit_rejection_stage = "pid_gate";
                    if (!pid_state_csv.empty()) pid_audit_out.push_back(std::move(audit));
                    continue;
                }
            }

            FitResult fit;
            audit.fit_attempted = 1;
            audit.fit_rejection_stage = "fit";
            if (use_kalman) {
                // UKF + RTS：χ²/ndf 记到 edm 里（Minuit 的 edm 没有物理意义，这里复用字段）
                KalmanFitResult kres;
                const bool fit_return = kfitter.Fit(cluster.points, est, kres);
                if (!fit_return || !kres.ok) {
                    audit.fit_rejection_stage = "kalman_fit";
                    audit.fit_cost_mm = kres.cost;
                    if (!pid_state_csv.empty()) pid_audit_out.push_back(std::move(audit));
                    continue;
                }
                fit.ok = kres.ok;
                fit.brho = kres.brho;
                fit.sigma_brho = kres.sigma_brho;
                fit.polar = kres.polar;
                fit.azimuthal = kres.azimuthal;
                fit.vx = kres.vx;
                fit.vy = kres.vy;
                fit.vz = kres.vz;
                fit.cost = kres.cost;
                fit.edm = (kres.ndf > 0) ? kres.chi2 / kres.ndf : 0.0;
                fit.status = kres.n_hits;
            } else {
                const bool fit_return = fitter.Fit(cluster.points, est, est.azimuthal, fit);
                audit.fit_cost_mm = fit.cost;
                audit.migrad_status = fit.status;
                audit.migrad_edm = fit.edm;
                if (!fit_return) {
                    audit.fit_rejection_stage = "minuit_acceptance_gate";
                    if (!pid_state_csv.empty()) pid_audit_out.push_back(std::move(audit));
                    continue;
                }
            }
            audit.fit_state = FitState::FIT_ACCEPTED;
            audit.fit_rejection_stage = "accepted";
            audit.fit_cost_mm = fit.cost;
            if (use_kalman) audit.migrad_status = -1;

            PidResult pid_result;
            const bool pid_result_ok = pid.Identify(fit.brho, est.dEdx, pid_result);
            set_pid_state(pid_result, pid_result_ok, fit.brho, est.dEdx);
            audit.pid_state_basis = "fitted_brho";
            if (!pid_state_csv.empty()) pid_audit_out.push_back(audit);

            FitRecord rec;
            rec.event = cloud.event;
            rec.cluster = cluster.label;
            rec.brho = fit.brho;
            rec.sigma_brho = fit.sigma_brho;
            rec.polar = fit.polar;
            rec.azimuthal = fit.azimuthal;
            rec.vx = fit.vx;
            rec.vy = fit.vy;
            rec.vz = fit.vz;
            rec.cost = fit.cost;
            rec.chi2ndf = use_kalman ? fit.edm : 0.0;
            rec.ok = fit.ok ? 1 : 0;
            rec.dedx = est.dEdx;
            rec.fit_state = static_cast<int>(FitState::FIT_ACCEPTED);
            rec.pid_state = static_cast<int>(audit.pid_state);
            rec.pid_unknown_reason = static_cast<int>(audit.pid_unknown_reason);
            rec.ic_observation_state = static_cast<int>(audit.ic_state);
            rec.pid_state_basis = static_cast<int>(PidStateBasis::FITTED_BRHO);
            rec.ic_amplitude = audit.ic_amplitude;
            rec.migrad_status = use_kalman ? -1 : fit.status;
            rec.migrad_edm = use_kalman ? -1.0 : fit.edm;

            if (pid_result_ok) {
                rec.pid_index = pid_result.index;
                rec.pid_sep = pid_result.separation;
            }
            out.push_back(rec);
        }
    };

    // ---- 先串行标定 PID（避免各线程标定结果依赖调度顺序）----
    double pid_scale = cfg.GetDouble("pid.scale", 0.0);
    if (pid_scale <= 0.0) {
        Clusterer cal_clusterer;
        cal_clusterer.SetParams(cluster_params);
        Estimator cal_estimator;
        ConfigureEstimator(cal_estimator,cfg);
        for (const auto& cloud : clouds) {
            std::vector<Cluster> clusters;
            cal_clusterer.Run(cloud, clusters);
            bool done = false;
            for (const auto& cluster : clusters) {
                EstimateResult est;
                if (!cal_estimator.Estimate(cluster, bfield, est)) continue;
                if (est.dEdx > 0.0 && est.brho > 0.0) {
                    ParticleId cal_pid;
                    cal_pid.Calibrate(est.brho, est.dEdx, particle_A, particle_Z);
                    pid_scale = cal_pid.Scale();
                    done = true;
                    break;
                }
            }
            if (done) break;
        }
        if (pid_scale <= 0.0) pid_scale = 1.0;
    }

    // ---- 并行处理 ----
    unsigned int n_threads = static_cast<unsigned int>(cfg.GetInt("run.threads", 0));
    // 环境变量覆盖（脚本里用 SPYRALPP_THREADS=8 scripts/... 就行，不必改配置）
    if (const char* env_threads = std::getenv("SPYRALPP_THREADS")) {
        const int v = std::atoi(env_threads);
        if (v > 0) n_threads = static_cast<unsigned int>(v);
    }
    if (n_threads == 0) {
        n_threads = std::thread::hardware_concurrency();
        if (n_threads == 0) n_threads = 1;
    }
    size_t shard_events = clouds.size();
    if (shard_count > 1) {
        shard_events = 0;
        for (size_t i = 0; i < clouds.size(); ++i) {
            if (static_cast<int>(i % static_cast<size_t>(shard_count)) == shard_index) {
                ++shard_events;
            }
        }
    }
    if (n_threads > shard_events && shard_events > 0) {
        n_threads = static_cast<unsigned int>(shard_events);
    }
    if (n_threads == 0) n_threads = 1;

    std::printf("并行处理：%u 线程，%zu 个点云\n", n_threads, clouds.size());

    // 拟合前的平滑样条降噪（可选）
    const bool smooth_enabled = cfg.GetBool("cluster.smooth", true);
    const double smooth_lambda = cfg.GetDouble("cluster.smooth_lambda", 1.0);

    // 可选：中间过程文件（cluster / estimate 两个阶段）
    std::vector<std::vector<ClusterRecord>> stage_clusters(n_threads);
    std::vector<std::vector<EstimateRecord>> stage_estimates(n_threads);
    std::vector<std::vector<PidCandidateAuditRow>> pid_audit_rows(n_threads);

    std::vector<std::vector<FitRecord>> per_thread(n_threads);
    std::atomic<size_t> next_index{0};

    std::mutex error_mutex;
    std::exception_ptr worker_error;
    std::atomic<bool> cancelled{false};
    auto worker = [&](unsigned int tid) {
        try {
        // 每个线程一套独立对象：TrackFitter / Clusterer / Estimator / ParticleId
        Clusterer clusterer;
        clusterer.SetParams(cluster_params);
        Estimator estimator;
        ConfigureEstimator(estimator,cfg);
        ParticleId pid;
        pid.SetScale(pid_scale);
        TrackFitter fitter;
        fitter.SetCenteredAzimuthBounds(cfg.GetBool(
            "solver.minuit_centered_azimuth_bounds",
            cfg.GetBool("solver.experimental_minuit_centered_azimuth_bounds", false)));
        fitter.SetSortCandidates(sort_minuit_candidates);
        if (!use_kalman && !fitter.Init(&mesh, mass, particle_Z))
            throw std::runtime_error("Track fitter initialization failed");
        KalmanFitter kfitter;
        // Detector config stores positive field magnitudes; match Spyral 1.0.
        const double kalman_bfield = bfield * cfg.GetDouble("solver.kalman_bfield_sign", -1.0);
        const double kalman_efield = cfg.GetDouble("solver.kalman_efield_vm", efield) *
                                     cfg.GetDouble("solver.kalman_efield_sign", -1.0);
        if (use_kalman &&
            !kfitter.Init(gas, particle_A, particle_Z, mass, kalman_bfield, kalman_efield, kalman_cfg)) {
            throw std::runtime_error("Kalman fitter initialization failed");
        }

        while (!cancelled.load()) {
            const size_t idx = next_index.fetch_add(1);
            if (idx >= clouds.size()) break;

            // 分片：多进程分工，避免重复处理
            if (shard_count > 1 &&
                static_cast<int>(idx % static_cast<size_t>(shard_count)) != shard_index) {
                continue;
            }

            process_cloud(clouds[idx], clusterer, estimator, pid, fitter, kfitter, pid_scale,
                          per_thread[tid], stage_clusters[tid], stage_estimates[tid],
                          pid_audit_rows[tid]);
        }
        } catch (...) {
            std::lock_guard<std::mutex> lk(error_mutex);
            if (!worker_error) worker_error = std::current_exception();
            cancelled.store(true);
        }
    };

    if (n_threads <= 1) {
        worker(0);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(n_threads);
        for (unsigned int t = 0; t < n_threads; ++t) {
            pool.emplace_back(worker, t);
        }
        for (auto& th : pool) th.join();
    }

    if (worker_error) {
        try { std::rethrow_exception(worker_error); }
        catch (const std::exception& e) { std::fprintf(stderr, "Analysis failed: %s\n", e.what()); }
        return 1;
    }
    // ---- 合并并按事件排序 ----
    std::vector<FitRecord> records;
    for (const auto& part : per_thread) {
        records.insert(records.end(), part.begin(), part.end());
    }
    std::sort(records.begin(), records.end(), [](const FitRecord& a, const FitRecord& b) {
        if (a.event != b.event) return a.event < b.event;
        return a.cluster < b.cluster;
    });

    if (!pid_state_csv.empty()) {
        std::vector<PidCandidateAuditRow> all_pid_audit;
        for (auto& part : pid_audit_rows)
            all_pid_audit.insert(all_pid_audit.end(), std::make_move_iterator(part.begin()),
                                 std::make_move_iterator(part.end()));
        std::sort(all_pid_audit.begin(), all_pid_audit.end(), [](const auto& a, const auto& b) {
            if (a.event != b.event) return a.event < b.event;
            return a.cluster_index < b.cluster_index;
        });
        if (!WritePidCandidateAuditCsv(pid_state_csv, all_pid_audit)) {
            std::fprintf(stderr, "PID state audit write failed: %s\n", pid_state_csv.c_str());
            return 1;
        }
    }

    if (dump_stages) {
        std::vector<ClusterRecord> all_clusters;
        std::vector<EstimateRecord> all_estimates;
        for (auto& v : stage_clusters)
            all_clusters.insert(all_clusters.end(), std::make_move_iterator(v.begin()),
                                std::make_move_iterator(v.end()));
        for (auto& v : stage_estimates)
            all_estimates.insert(all_estimates.end(), std::make_move_iterator(v.begin()),
                                 std::make_move_iterator(v.end()));
        const std::string stage_path = std::string(argv[2]) + "_stages.root";
        if (WriteStages(stage_path, all_clusters, all_estimates)) {
            std::printf("中间过程文件：%s（%zu 个簇，%zu 条估计）\n", stage_path.c_str(),
                        all_clusters.size(), all_estimates.size());
        }
    }

    if (!WriteFitResults(argv[2], records)) {
        std::fprintf(stderr, "结果写出失败\n");
        return 1;
    }

    if (shard_count > 1) {
        std::printf("重建完成：分片 %d/%d，共 %zu 条点云（本片 %zu），%zu 条拟合结果\n",
                    shard_index, shard_count, clouds.size(), shard_events, records.size());
    } else {
        std::printf("重建完成：%zu 条点云，%zu 条拟合结果\n", clouds.size(), records.size());
    }
    std::printf("  输出：%s\n", argv[2]);
    return 0;
}
