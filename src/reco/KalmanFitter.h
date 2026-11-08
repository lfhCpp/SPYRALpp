#ifndef RK4_RECO_KALMANFITTER_H
#define RK4_RECO_KALMANFITTER_H

//----------------------------------------------------------------------
// 轨迹拟合：无迹卡尔曼滤波（UKF）+ RTS 平滑
//
// 这是一个**可选的**拟合后端，和 TrackFitter（TMinuit + 轨迹网格插值）并列，
// 由配置项 solver.fitter = kalman 选用。方法取自 ATTPCROOTv2 的
// OpenKF-Claude 分支（EventFit::AtFitterUKF + kf::TrackFitterUKF），
// 那里是"ATIMA/catima 能损 + 逐点 UKF 递推 + RTS 平滑"的组合；
// 本项目按自己的代码风格重写，去掉了 FairRoot/GEANT 依赖：
//
//   状态      x = (x, y, z, p, θ, φ)        位置 mm、动量 MeV/c、角度 rad
//   传播      catima dE/dx + 磁场（ParticleDynamics::equationMotion, RK4）
//   量测      z = (x, y, z)                 σ = meas_sigma（默认 0.5 mm）
//   过程噪声  Q_aug：dE/dx 缩放因子的方差（能损straggling，可关）
//             Q_mod：位置/动量/角度的模型噪声底，可选加多重散射
//   平滑      RTS，用 sigma 点算出的互相关 C_k
//
// 与 TrackFitter 的差别：
//   - 不需要预生成轨迹网格（TrackMesh），逐点传播；
//   - 天然给出每个点的协方差，χ²/ndf 可直接用作质量判据；
//   - 顶点由首个平滑测量状态线性外推到束流轴（最多 50 mm），并用
//     CATIMA dE/dx 近似补回能量；这里不是 RK4/CATIMA 反向传播。
//
// UKF 的 sigma 点 / 加权均值协方差 / 互相关用的是 extern/openkf 提供的
// kf::Matrix / kf::util::*（GPL-3.0，见 extern/openkf/COPYING）。
//----------------------------------------------------------------------

#include "../sim/TargetDensity.h"
#include "Cluster.h"
#include "Estimator.h"

#include <array>
#include <memory>
#include <vector>

namespace rk4 {

//----------------------------------------------------------------------
// 参数（默认值取自 a1954 12Be 分析 guide 的 UKF 一节：
// measSigma 0.5 mm、momSigmaFrac 0.1、能损 straggling 关闭）
//----------------------------------------------------------------------
struct KalmanConfig {
    // 初始化
    double meas_sigma_mm = 0.5;   // 单点位置测量 σ
    bool per_cluster_covariance = false; // ATTPCROOT-style detector-resolution covariance per smoothed cluster
    double cluster_coef_t_cm2_us = 9.0e-5;
    double cluster_coef_l_cm2_us = 9.0e-7;
    double cluster_drift_velocity_cm_us = 1.0;
    double cluster_tb_time_us = 0.320;
    double cluster_pad_resolution_xy_mm = 2.3;
    double mom_sigma_frac = 0.1;  // 动量先验的相对不确定度
    double pos_sigma_mm = 0.5;    // ATTPCROOT 初始位置协方差与量测 σ 相同
    double ang_sigma_deg = 1.0;   // ATTPCROOT 初始方向先验 σ

    // 模型噪声（加在预测协方差上的对角项）
    double q_pos = 1e-4;          // mm^2
    double q_mom = 1e-4;          // (MeV/c)^2
    double q_ang = 1e-6;          // rad^2

    // 过程噪声：能损 straggling（默认关，与 a1954 分析一致）
    bool straggling = false;
    double straggle_cap = 1.0 / 3.0;  // 缩放的方差上限，避免 sigma 点乱飞

    // 多重散射（Highland），默认关
    bool mult_scattering = false;
    double rad_length_mm = 7.6e6;     // H2 600 torr 量级

    // 拟合控制
    int measurements = 30;        // 兼容旧配置；Be12现按自适应空间聚类生成量测点
    bool verbose = false;         // 打印失败原因（调试用）
    int min_clusters = 10;        // a1954 Be12 UKF 生产配置
    int gap_max_steps = 100000;    // 无 pad 区后的扩展 z 面传播上限
    bool smooth = true;           // RTS 平滑
    bool back_extrapolate = true; // 反推到顶点
    double back_extrap_max_mm = 50.0;
    double z_pad_plane_mm = 0.0;    // >0 converts digi z to lab z: Z_lab = plane - Z_digi
    double step_mm = 0.1;         // 对齐 AtRK4Stepper 的 0.1 mm 步长
    int max_steps = 30000;        // 单次传播上限；覆盖典型无 pad 区后的螺旋接续
    double min_ke_mev = 0.01;     // ATTPCROOT AtPropagator 的停止阈值
    double chi2_max = 0.0;        // > 0 时 χ²/ndf 超过即判失败
    double innovation_chi2_max = 25.0; // Optional per-hit outlier gate; recovery pass can disable it.
    double min_measurement_fraction = 0.20; // 至少保留五分之一的有效测量，拒绝不完整拟合
    int iterations = 1;              // >1 enables ATTPCROOT-style smoothed-reference passes
    double ref_covariance_inflation = 4.0; // Keep reference-track re-pass covariance conservative

    // UKF 权重参数（alpha, beta, kappa）
    double alpha = 1.0e-3;
    double beta = 2.0;
    double kappa = 0.0;
};

//----------------------------------------------------------------------
// 拟合结果（字段与 TrackFitter::FitResult 对齐，便于写到同一棵 TTree）
//----------------------------------------------------------------------
struct KalmanFitResult {
    bool ok = false;
    // 0 success; 2/3/5 too few hits; 4 bad seed; 6 singular covariance;
    // 7 gap could not be rejoined; 8/9 too few accepted measurements; 10 χ² gate; 11 invalid final state.
    int failure_code = 0;
    double brho = 0.0;       // T*m（顶点处）
    double sigma_brho = 0.0; // T*m
    double polar = 0.0;      // rad（顶点处，从 +z 量起）
    double azimuthal = 0.0;  // rad
    double vx = 0.0, vy = 0.0, vz = 0.0; // mm
    double ke_vertex = 0.0;  // MeV
    double chi2 = 0.0;
    int ndf = 0;
    double cost = 0.0;       // 平滑轨迹到命中点的平均距离 mm
    int n_hits = 0;
};

// Experimental one-step UKF interface for covariance-guided track growing.
// Arrays use row-major covariance layout. These methods are additive: the
// production Fit() path and its defaults are unchanged.
struct KalmanTrackState {
    std::array<double, 6> mean{};       // x,y,z [mm], p [MeV/c], theta,phi [rad]
    std::array<double, 36> covariance{}; // row-major 6x6
};

struct KalmanMeasurementPrediction {
    KalmanTrackState predicted_state;
    std::array<double, 3> xyz_mean{};
    std::array<double, 9> innovation_covariance{}; // predicted XYZ covariance + R
    std::array<double, 18> state_measurement_cross_covariance{}; // row-major 6x3
    double path_mm = 0.0;
    bool valid = false;
};

struct KalmanMeasurementUpdate {
    KalmanTrackState filtered_state;
    std::array<double, 3> innovation{};
    double innovation_chi2 = 0.0;
    bool accepted = false;
};

//----------------------------------------------------------------------
// 拟合器。每个线程一个实例（内部含 ParticleDynamics 状态）。
//----------------------------------------------------------------------
class KalmanFitter {
public:
    KalmanFitter();
    ~KalmanFitter();

    KalmanFitter(const KalmanFitter&) = delete;
    KalmanFitter& operator=(const KalmanFitter&) = delete;

    // gas: 靶（决定 dE/dx 的材料与密度）；A/Z: 被拟合粒子的质量数/电荷数
    bool Init(const Target& gas, int A, int Z, double mass_mev, double bfield_t,
              double efield_vm, const KalmanConfig& cfg);

    bool Initialized() const;

    // hits: 同一个簇内的点云点（内部会按 z 排序，必要时反向）
    // guess: Estimator 的粗估（用它的 Bρ 做动量初值；方向用点序列定）
    bool Fit(const std::vector<CloudPoint>& hits, const EstimateResult& guess,
             KalmanFitResult& out);

    // Propagate the state/covariance by a nonzero signed arc length and
    // predict an XYZ measurement. Negative path integrates the deterministic
    // mean transport backward; it is unavailable when stochastic energy-loss
    // straggling is enabled. R is a row-major 3x3 measurement covariance.
    // Update() applies the predicted cross-covariance and the configured
    // Mahalanobis innovation gate. These are experimental primitives for
    // active hit association; they do not search or consume point-cloud hits.
    bool PredictPath(const KalmanTrackState& prior, double path_mm,
                     const std::array<double, 9>& measurement_covariance,
                     KalmanMeasurementPrediction& out);
    bool Update(const KalmanMeasurementPrediction& prediction,
                const std::array<double, 3>& measurement,
                KalmanMeasurementUpdate& out) const;

    const KalmanConfig& Config() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace rk4

#endif // RK4_RECO_KALMANFITTER_H
