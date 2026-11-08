#include "KalmanFitter.h"

#include "../core/MassData.h"
#include "../core/Units.h"
#include "../sim/ParticleDynamics.h"

// extern/openkf：UKF 用的矩阵类型与列操作工具（GPL-3.0）
#include "kf_util.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace rk4 {

namespace {

//----------------------------------------------------------------------
// 维度：状态 6、量测 3、过程噪声 1（dE/dx 缩放因子）
//----------------------------------------------------------------------
constexpr int NX = 6;
constexpr int NZ = 3;
constexpr int NV = 1;
constexpr int NA = NX + NV;      // 增广状态维
constexpr int NS = 2 * NA + 1;   // sigma 点数（15）

using VecX = Eigen::Matrix<double, NX, 1>;
using MatX = Eigen::Matrix<double, NX, NX>;
using VecZ = Eigen::Matrix<double, NZ, 1>;
using MatZ = Eigen::Matrix<double, NZ, NZ>;
using VecA = Eigen::Matrix<double, NA, 1>;
using MatA = Eigen::Matrix<double, NA, NA>;
using MatAxNS = Eigen::Matrix<double, NA, NS>;
using MatXNS = Eigen::Matrix<double, NX, NS>;
using MatZNS = Eigen::Matrix<double, NZ, NS>;
using MatXZ = Eigen::Matrix<double, NX, NZ>;

constexpr double kBrStiffness = 0.299792458; // Bρ[T·m] = p[GeV/c]/(0.299792458·Z)

//----------------------------------------------------------------------
// 状态换算：积分器用 (x,y,z) m + (γβx, γβy, γβz)，滤波器用 mm/MeV/c/rad
//----------------------------------------------------------------------
inline double GvFromP(double p_mev, double mass) { return p_mev / mass; }
inline double GammaOfGv(double gv) { return std::sqrt(1.0 + gv * gv); }
inline double BetaOfGv(double gv) {
    const double g = GammaOfGv(gv);
    return (g > 0.0) ? gv / g : 0.0;
}

void ToIntegratorState(const VecX& x, double mass, double s[6]) {
    s[0] = x[0] * 1.0e-3;
    s[1] = x[1] * 1.0e-3;
    s[2] = x[2] * 1.0e-3;
    const double gv = GvFromP(x[3], mass);
    const double st = std::sin(x[4]);
    s[3] = gv * st * std::cos(x[5]);
    s[4] = gv * st * std::sin(x[5]);
    s[5] = gv * std::cos(x[4]);
}

VecX FromIntegratorState(const double s[6], double mass) {
    VecX x;
    x[0] = s[0] * 1.0e3;
    x[1] = s[1] * 1.0e3;
    x[2] = s[2] * 1.0e3;
    const double gv = std::sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
    x[3] = gv * mass;
    x[4] = (gv > 0.0) ? std::acos(std::max(-1.0, std::min(1.0, s[5] / gv))) : 0.0;
    x[5] = std::atan2(s[4], s[3]);
    return x;
}

void UnitDir(const VecX& x, double out[3]) {
    const double st = std::sin(x[4]);
    out[0] = st * std::cos(x[5]);
    out[1] = st * std::sin(x[5]);
    out[2] = std::cos(x[4]);
}

// The RTS update is linear in the state vector, but polar and azimuth are
// periodic coordinates. Keep the reported direction on the physical sphere
// after each correction/smoothing update.
void NormalizeAngles(VecX& x) {
    while (x[4] < 0.0 || x[4] > M_PI) {
        if (x[4] < 0.0) x[4] = -x[4];
        else x[4] = 2.0 * M_PI - x[4];
        x[5] += M_PI;
    }
    x[5] = std::fmod(x[5], 2.0 * M_PI);
    if (x[5] < 0.0) x[5] += 2.0 * M_PI;
}

//----------------------------------------------------------------------
// Cholesky 分解（带抖动，保证 sigma 点总能算出来）
//----------------------------------------------------------------------
MatA CholeskyL(MatA P) {
    P = 0.5 * (P + P.transpose());
    Eigen::LLT<MatA> llt(P);
    if (llt.info() == Eigen::Success) return llt.matrixL();

    for (int i = 0; i < 4; ++i) {
        MatA R = 0.5 * (P + P.transpose());
        const double jitter = std::pow(10.0, -9.0 + i);
        R.diagonal().array() += jitter * (1.0 + R.diagonal().cwiseAbs().maxCoeff());
        Eigen::LLT<MatA> l2(R);
        if (l2.info() == Eigen::Success) return l2.matrixL();
        P = R;
    }
    MatA D = MatA::Zero();
    const double scale = 1.0 + P.diagonal().cwiseAbs().maxCoeff();
    D.diagonal().setConstant(std::sqrt(scale) * 1.0e-3);
    return D;
}

//----------------------------------------------------------------------
// 无迹变换：权重、sigma 点、加权均值/协方差、互相关
// （结构与 extern/openkf 的 unscented_kalman_filter.h 相同；区别只在
//   增广里只放"状态 + dE/dx 缩放"，量测噪声在校正时单独加 —— 与
//   ATTPCROOTv2 的 TrackFitterUKFBase 的做法一致）
//----------------------------------------------------------------------
struct UkfCore {
    double lambda = 0.0;
    double wm0 = 0.0, wc0 = 0.0, wi = 0.0;

    void SetWeights(double alpha, double beta, double kappa) {
        lambda = alpha * alpha * (NA + kappa) - NA;
        const double deno = lambda + NA;
        wm0 = lambda / deno;
        wc0 = wm0 + (1.0 - alpha * alpha + beta);
        wi = 0.5 / deno;
    }

    MatAxNS SigmaPoints(const VecA& xa, const MatA& Pa) const {
        const double mult = std::sqrt(NA + lambda);
        const MatA S = CholeskyL(Pa);
        MatAxNS sig;
        kf::util::copyToColumn<NA, NS>(0, sig, xa);
        for (int i = 0; i < NA; ++i) {
            const int i1 = i + 1;
            const int i2 = i + NA + 1;
            kf::util::copyToColumn<NA, NS>(i1, sig, xa);
            kf::util::copyToColumn<NA, NS>(i2, sig, xa);
            Eigen::Matrix<double, NA, 1> shift = mult * S.col(i);
            kf::util::addColumnFrom<NA, NS>(i1, sig, shift);
            kf::util::subColumnFrom<NA, NS>(i2, sig, shift);
        }
        return sig;
    }

    template <int D>
    void MeanCov(const Eigen::Matrix<double, D, NS>& sig, Eigen::Matrix<double, D, 1>& mean,
                 Eigen::Matrix<double, D, D>& cov) const {
        mean = wm0 * kf::util::getColumnAt<D, NS>(0, sig);
        for (int i = 1; i < NS; ++i) mean += wi * kf::util::getColumnAt<D, NS>(i, sig);

        Eigen::Matrix<double, D, 1> d = kf::util::getColumnAt<D, NS>(0, sig) - mean;
        cov = wc0 * d * d.transpose();
        for (int i = 1; i < NS; ++i) {
            d = kf::util::getColumnAt<D, NS>(i, sig) - mean;
            cov += wi * d * d.transpose();
        }
        cov = 0.5 * (cov + cov.transpose());
    }

    template <int D>
    Eigen::Matrix<double, NX, D> Cross(const MatXNS& sigX, const VecX& meanX,
                                       const Eigen::Matrix<double, D, NS>& sigY,
                                       const Eigen::Matrix<double, D, 1>& meanY) const {
        Eigen::Matrix<double, NX, 1> dx = kf::util::getColumnAt<NX, NS>(0, sigX) - meanX;
        Eigen::Matrix<double, D, 1> dy = kf::util::getColumnAt<D, NS>(0, sigY) - meanY;
        Eigen::Matrix<double, NX, D> P = wc0 * (dx * dy.transpose());
        for (int i = 1; i < NS; ++i) {
            dx = kf::util::getColumnAt<NX, NS>(i, sigX) - meanX;
            dy = kf::util::getColumnAt<D, NS>(i, sigY) - meanY;
            P += wi * (dx * dy.transpose());
        }
        return P;
    }
};

//----------------------------------------------------------------------
// 传播到测量点的最近通过点，复现 AtMeasurementPoint::PassedSurface：
// (point - position)·momentum 的符号在最近点前后翻转。
//----------------------------------------------------------------------
struct Propagator {
    ParticleDynamics* dyn = nullptr;
    double mass = 0.0;
    double step_mm = 5.0;
    int max_steps = 30000;
    int gap_max_steps = 100000;
    int plane_max_steps = 20000;
    double min_ke = 0.01;
    bool verbose = false;

    // Physics propagation over a signed geometric arc length. Negative path
    // integrates the deterministic equation of motion backward in time, which
    // reverses the Lorentz and CATIMA stopping terms together. It is an
    // approximate inverse of the mean transport only; stochastic straggling
    // cannot be inverted by this operation.
    bool ByPath(const VecX& in, double eloss_scale, double path_mm,
                VecX& out, double* actual_path_mm = nullptr) const {
        if (dyn == nullptr || path_mm == 0.0 || !(step_mm > 0.0) ||
            !std::isfinite(path_mm) || !std::isfinite(eloss_scale) || eloss_scale <= 0.0)
            return false;
        const double total_path_mm=std::abs(path_mm);
        const double time_sign=path_mm>0.0?1.0:-1.0;
        double s[6];
        ToIntegratorState(in, mass, s);
        dyn->SetInitialConditions();
        dyn->SetELossScale(eloss_scale);
        double t=0.0, traversed_mm=0.0;
        int steps=0;
        const int step_limit=std::min(max_steps,
            std::max(1, static_cast<int>(std::ceil(total_path_mm / step_mm)) + 4));
        while (traversed_mm + 1.0e-9 < total_path_mm && steps++ < step_limit) {
            const double gv=std::sqrt(s[3]*s[3]+s[4]*s[4]+s[5]*s[5]);
            if (!(gv>0.0)) return false;
            const double ke=mass*(GammaOfGv(gv)-1.0);
            if (!(ke>min_ke)) return false;
            const double beta=BetaOfGv(gv);
            const double v=beta*Units::C;
            if (!(v>0.0)) return false;
            const double requested_mm=std::min(step_mm,total_path_mm-traversed_mm);
            const double dt=time_sign*requested_mm*1.0e-3/v;
            rk4_solver::Integrator<6,ParticleDynamics> integ(
                *dyn,&ParticleDynamics::equationMotion,dt);
            double t_next=t, s_next[6];
            integ.step(t,s,t_next,s_next);
            for(double value:s_next) if(!std::isfinite(value)) return false;
            double actual_mm=1.0e3*std::sqrt(
                (s_next[0]-s[0])*(s_next[0]-s[0])+
                (s_next[1]-s[1])*(s_next[1]-s[1])+
                (s_next[2]-s[2])*(s_next[2]-s[2]));
            if (!(actual_mm>0.0) || !std::isfinite(actual_mm)) return false;
            // A final fractional RK4 step prevents arc-length overshoot.
            const double remain_mm=total_path_mm-traversed_mm;
            if (actual_mm>remain_mm*(1.0+1.0e-10)) {
                const double frac=std::clamp(remain_mm/actual_mm,0.0,1.0);
                rk4_solver::Integrator<6,ParticleDynamics> final_step(
                    *dyn,&ParticleDynamics::equationMotion,dt*frac);
                double t_final=t, s_final[6];
                final_step.step(t,s,t_final,s_final);
                for(double value:s_final) if(!std::isfinite(value)) return false;
                actual_mm=1.0e3*std::sqrt(
                    (s_final[0]-s[0])*(s_final[0]-s[0])+
                    (s_final[1]-s[1])*(s_final[1]-s[1])+
                    (s_final[2]-s[2])*(s_final[2]-s[2]));
                for(int i=0;i<6;++i) s[i]=s_final[i];
                t=t_final;
            } else {
                for(int i=0;i<6;++i) s[i]=s_next[i];
                t=t_next;
            }
            traversed_mm+=actual_mm;
        }
        if (traversed_mm+1.0e-6<total_path_mm) return false;
        out=FromIntegratorState(s,mass);
        if (!out.allFinite() || !(out[3]>0.0)) return false;
        if (actual_path_mm!=nullptr) *actual_path_mm=traversed_mm;
        return true;
    }

    bool ToPoint(const VecX& in, double eloss_scale, const double point_m[3],
                 VecX& out, double* path_mm = nullptr) const {
        double s[6];
        ToIntegratorState(in, mass, s);
        dyn->SetELossScale(eloss_scale);

        double t = 0.0;
        double path = 0.0;
        // Estimate the path needed to reach a reappearing hit after a pad
        // coverage gap. z advances monotonically for these forward tracks;
        // allow a modest margin for the helix, while retaining the configured
        // hard cap for points that are actually out of order/unreachable.
        const double dz_mm = std::abs(point_m[2] * 1.0e3 - in[2]);
        const double cos_theta = std::max(0.05, std::abs(std::cos(in[4])));
        const double estimated_steps = 1.25 * dz_mm / (cos_theta * step_mm) + 32.0;
        const int step_limit = std::min(max_steps,
            std::max(20000, static_cast<int>(std::min(estimated_steps, static_cast<double>(max_steps)))));
        int steps = 0;
        while (steps++ < step_limit) {
            const double gv = std::sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
            if (!(gv > 0.0)) {
                out = in;
                if (path_mm != nullptr) *path_mm = 0.0;
                return true; // ROOT leaves a stopped mean state unchanged.
            }
            const double ke = mass * (GammaOfGv(gv) - 1.0);
            if (!(ke > min_ke)) {
                out = in;
                if (path_mm != nullptr) *path_mm = 0.0;
                return true; // AtPropagator returns early below its stop tolerance.
            }

            const double v = BetaOfGv(gv) * Units::C;
            if (!(v > 0.0)) return false;
            const double dt = (step_mm * 1.0e-3) / v;

            rk4_solver::Integrator<6, ParticleDynamics> integ(*dyn, &ParticleDynamics::equationMotion, dt);
            double t_next = t;
            double s_next[6];
            integ.step(t, s, t_next, s_next);

            for (int k = 0; k < 6; ++k) {
                if (!std::isfinite(s_next[k])) return false;
            }

            const double gv_next = std::sqrt(s_next[3] * s_next[3] + s_next[4] * s_next[4] + s_next[5] * s_next[5]);
            if (!(gv_next > 0.0)) return false;
            const double dx0 = point_m[0] - s[0], dy0 = point_m[1] - s[1], dz0 = point_m[2] - s[2];
            const double dx1 = point_m[0] - s_next[0], dy1 = point_m[1] - s_next[1], dz1 = point_m[2] - s_next[2];
            const double before = (dx0 * s[3] + dy0 * s[4] + dz0 * s[5]) / gv;
            const double after = (dx1 * s_next[3] + dy1 * s_next[4] + dz1 * s_next[5]) / gv_next;
            // State positions are propagated in metres; the public path and
            // process-noise configuration are in millimetres.
            path += 1.0e3 * std::sqrt((s_next[0] - s[0]) * (s_next[0] - s[0]) +
                                      (s_next[1] - s[1]) * (s_next[1] - s[1]) +
                                      (s_next[2] - s[2]) * (s_next[2] - s[2]));

            if (before * after <= 0.0) {
                // Match AtPropagator: take a final RK4 step whose path length
                // is the distance from the previous state to the measured point.
                const double final_mm = 1.0e3 * std::sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
                const double beta = BetaOfGv(gv);
                const double final_dt = (final_mm * 1.0e-3) / (beta * Units::C);
                rk4_solver::Integrator<6, ParticleDynamics> final_step(*dyn, &ParticleDynamics::equationMotion, final_dt);
                double t_final = t;
                double s_final[6];
                final_step.step(t, s, t_final, s_final);
                for (double value : s_final) if (!std::isfinite(value)) return false;
                out = FromIntegratorState(s_final, mass);
                path += final_mm;
                if (path_mm != nullptr) *path_mm = path;
                return true;
            }

            for (int k = 0; k < 6; ++k) s[k] = s_next[k];
            t = t_next;
        }
        if (verbose) std::fprintf(stderr, "[kalman] point propagation exceeded %d steps\n", step_limit);
        return false;
    }

    // In the Spyral workspace, accepted hit sequences have a known forward or
    // backward z direction. A long pad-coverage gap can span several turns of
    // the helix, so the tangent-plane closest-approach test in ToPoint may stop
    // at an earlier turn. Propagate the mean state to the unique measured-z
    // plane instead; sigma points are still transported to the ROOT-style
    // momentum-normal plane through that measurement below.
    bool ToZPlane(const VecX& in, double eloss_scale, const double point_m[3],
                  VecX& out, double* path_mm = nullptr,
                  bool allow_local_backward = false) const {
        double s[6];
        ToIntegratorState(in, mass, s);
        dyn->SetELossScale(eloss_scale);

        const double target_z = point_m[2];
        const double z0 = s[2];
        const double dz = target_z - z0;
        const double direction = std::cos(in[4]);
        // The z-plane fallback is selected for a large wrong-turn mismatch.
        // Only that case allows a short backward local re-anchor; ordinary
        // z-order noise keeps the historic forward closest-approach behavior.
        if (std::abs(direction) < 0.03) return ToPoint(in, eloss_scale, point_m, out, path_mm);
        if (dz * direction < -1.0e-8 && !allow_local_backward)
            return ToPoint(in, eloss_scale, point_m, out, path_mm);
        if (dz * direction < -1.0e-8 && allow_local_backward) {
            // The mean state has advanced only a short distance past the
            // unique z plane; locally project it back instead of searching a
            // tangent plane on a different helix turn.
            const double gv = std::sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
            const double arc_m = dz / direction;
            if (!(gv > 0.0) || std::abs(arc_m) > 0.1)
                return ToPoint(in, eloss_scale, point_m, out, path_mm);
            if (verbose)
                std::fprintf(stderr, "[kalman] local z-plane projection dz=%.3f mm arc=%.3f mm\n",
                             dz * 1.0e3, arc_m * 1.0e3);
            s[0] += (s[3] / gv) * arc_m;
            s[1] += (s[4] / gv) * arc_m;
            s[2] = target_z;
            out = FromIntegratorState(s, mass);
            if (path_mm != nullptr) *path_mm = std::abs(arc_m) * 1.0e3;
            return true;
        }

        const double dz_mm = std::abs(dz) * 1.0e3;
        const double estimated_steps = 2.5 * dz_mm /
            (std::max(0.03, std::abs(direction)) * step_mm) + 32.0;
        const int step_limit = std::min(gap_max_steps,
            std::max(20000, static_cast<int>(std::min(estimated_steps, static_cast<double>(gap_max_steps)))));

        double t = 0.0;
        double path = 0.0;
        int steps = 0;
        bool stopped = false;
        while (steps++ < step_limit) {
            const double gv = std::sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
            if (!(gv > 0.0)) { stopped = true; break; }
            const double ke = mass * (GammaOfGv(gv) - 1.0);
            if (!(ke > min_ke)) { stopped = true; break; }
            const double v = BetaOfGv(gv) * Units::C;
            if (!(v > 0.0)) return false;
            const double dt = (step_mm * 1.0e-3) / v;

            rk4_solver::Integrator<6, ParticleDynamics> integ(*dyn, &ParticleDynamics::equationMotion, dt);
            double t_next = t;
            double s_next[6];
            integ.step(t, s, t_next, s_next);
            for (double value : s_next) if (!std::isfinite(value)) return false;

            const double dz0 = target_z - s[2];
            const double dz1 = target_z - s_next[2];
            path += 1.0e3 * std::sqrt(
                (s_next[0] - s[0]) * (s_next[0] - s[0]) +
                (s_next[1] - s[1]) * (s_next[1] - s[1]) +
                (s_next[2] - s[2]) * (s_next[2] - s[2]));
            if (std::abs(dz0) < 1.0e-12) {
                out = FromIntegratorState(s, mass);
                if (path_mm != nullptr) *path_mm = path;
                return true;
            }
            if (dz0 * dz1 <= 0.0) {
                const double frac = std::clamp(dz0 / (s[2] - s_next[2]), 0.0, 1.0);
                double t_final = t;
                double s_final[6];
                rk4_solver::Integrator<6, ParticleDynamics> final_step(
                    *dyn, &ParticleDynamics::equationMotion, dt * frac);
                final_step.step(t, s, t_final, s_final);
                for (double value : s_final) if (!std::isfinite(value)) return false;
                s_final[2] = target_z;
                out = FromIntegratorState(s_final, mass);
                if (path_mm != nullptr) *path_mm = path -
                    1.0e3 * (1.0 - frac) * std::sqrt(
                        (s_next[0] - s[0]) * (s_next[0] - s[0]) +
                        (s_next[1] - s[1]) * (s_next[1] - s[1]) +
                        (s_next[2] - s[2]) * (s_next[2] - s[2]));
                return true;
            }
            for (int k = 0; k < 6; ++k) s[k] = s_next[k];
            t = t_next;
        }
        if (verbose) {
            if (stopped)
                std::fprintf(stderr, "[kalman] particle stopped before z-plane (dz=%.1f mm, z=%.1f mm)\n",
                             dz_mm, s[2] * 1.0e3);
            else
                std::fprintf(stderr, "[kalman] z-plane propagation exceeded %d steps (dz=%.1f mm)\n",
                             step_limit, dz_mm);
        }
        return false;
    }

    // ATTPCROOT first propagates the mean to the point of closest approach,
    // then transports every UKF sigma point to the same plane through the
    // measurement, normal to that mean momentum. Propagating each sigma point
    // to its own closest-approach point changes the UKF transition function.
    bool ToPlane(const VecX& in, double eloss_scale, const double point_m[3],
                 const double normal[3], VecX& out,
                 double preferred_path_mm = -1.0) const {
        double s[6];
        ToIntegratorState(in, mass, s);
        dyn->SetELossScale(eloss_scale);

        const auto signed_distance = [&](const double state[6]) {
            return (state[0] - point_m[0]) * normal[0] +
                   (state[1] - point_m[1]) * normal[1] +
                   (state[2] - point_m[2]) * normal[2];
        };
        double t = 0.0;
        double traversed_mm = 0.0;
        double best_delta_mm = std::numeric_limits<double>::infinity();
        VecX best_crossing;
        const double path_window_mm = std::max(3.0 * step_mm,
            15.0 + 0.05 * std::max(0.0, preferred_path_mm));
        bool have_crossing = false;
        int steps = 0;
        while (steps++ < plane_max_steps) {
            const double gv = std::sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
            if (!(gv > 0.0)) { out = in; return true; }
            const double ke = mass * (GammaOfGv(gv) - 1.0);
            if (!(ke > min_ke)) { out = in; return true; }
            const double v = BetaOfGv(gv) * Units::C;
            if (!(v > 0.0)) return false;

            double t_next = t;
            double s_next[6];
            rk4_solver::Integrator<6, ParticleDynamics> integ(
                *dyn, &ParticleDynamics::equationMotion, (step_mm * 1.0e-3) / v);
            integ.step(t, s, t_next, s_next);
            for (double value : s_next) if (!std::isfinite(value)) return false;

            const double d0 = signed_distance(s);
            const double d1 = signed_distance(s_next);
            if ((d0 > 0.0) != (d1 > 0.0)) {
                // AtMeasurementPlane uses the absolute signed distance from
                // the previous RK state for its final RK step, then projects
                // the result exactly onto the plane.
                const double final_m = std::abs(d0);
                const double final_dt = final_m / v;
                rk4_solver::Integrator<6, ParticleDynamics> final_step(
                    *dyn, &ParticleDynamics::equationMotion, final_dt);
                double t_final = t;
                double s_final[6];
                final_step.step(t, s, t_final, s_final);
                for (double value : s_final) if (!std::isfinite(value)) return false;
                const double residual = signed_distance(s_final);
                for (int k = 0; k < 3; ++k) s_final[k] -= residual * normal[k];
                const double step_path_mm = 1.0e3 * std::sqrt(
                    (s_next[0] - s[0]) * (s_next[0] - s[0]) +
                    (s_next[1] - s[1]) * (s_next[1] - s[1]) +
                    (s_next[2] - s[2]) * (s_next[2] - s[2]));
                const double crossing_path_mm = traversed_mm +
                    std::min(step_path_mm, 1.0e3 * final_m);
                if (preferred_path_mm < 0.0) {
                    out = FromIntegratorState(s_final, mass);
                    return true;
                }
                const double delta = std::abs(crossing_path_mm - preferred_path_mm);
                if (delta < best_delta_mm) {
                    best_delta_mm = delta;
                    best_crossing = FromIntegratorState(s_final, mass);
                    have_crossing = true;
                }
                if (delta <= path_window_mm) {
                    out = best_crossing;
                    return true;
                }
            }
            traversed_mm += 1.0e3 * std::sqrt(
                (s_next[0] - s[0]) * (s_next[0] - s[0]) +
                (s_next[1] - s[1]) * (s_next[1] - s[1]) +
                (s_next[2] - s[2]) * (s_next[2] - s[2]));
            if (preferred_path_mm >= 0.0 && have_crossing &&
                traversed_mm > preferred_path_mm + path_window_mm)
                break;
            for (int k = 0; k < 6; ++k) s[k] = s_next[k];
            t = t_next;
        }
        if (have_crossing && best_delta_mm <= path_window_mm) {
            out = best_crossing;
            return true;
        }
        return false;
    }
};

//----------------------------------------------------------------------
// 能损 straggling：Bohr 近似（MeV²），只用来给过程噪声一个量级。
// 默认关闭（a1954 的 UKF 分析也是关的）。
//----------------------------------------------------------------------
double BohrStraggleSigma2(double beta, double charge_z, double z_target, double rho_L_gcm2) {
    if (!(beta > 0.0) || !(rho_L_gcm2 > 0.0)) return 0.0;
    const double N_A = 6.02214076e23;
    const double r_e_cm = 2.8179403262e-13;
    const double me_c2 = 0.51099895;
    const double z2 = charge_z * charge_z;
    const double s2 = 4.0 * M_PI * N_A * r_e_cm * r_e_cm * me_c2 * z2 * z_target * rho_L_gcm2 *
                      (1.0 - 0.5 * beta * beta) / (beta * beta);
    return (s2 > 0.0) ? s2 : 0.0;
}

} // namespace

//----------------------------------------------------------------------
// Impl
//----------------------------------------------------------------------
struct KalmanFitter::Impl {
    KalmanConfig cfg;
    double mass = 0.0;
    int A = 1;
    int Z = 1;
    double bfield = 2.85;
    double efield = 0.0;
    bool ready = false;

    ParticleDynamics dyn;
    Propagator prop;
    UkfCore ukf;

    bool Init(const Target& gas, int a_in, int z_in, double mass_mev, double b_t, double e_vm,
              const KalmanConfig& c) {
        ready = false;
        if (a_in<1 || z_in<1 || z_in>a_in || !(mass_mev>0) ||
            !std::isfinite(mass_mev) || !std::isfinite(b_t) || !std::isfinite(e_vm)) return false;
        cfg = c;
        if (cfg.per_cluster_covariance &&
            (!(cfg.cluster_coef_t_cm2_us >= 0.0) || !(cfg.cluster_coef_l_cm2_us >= 0.0) ||
             !(cfg.cluster_drift_velocity_cm_us > 0.0) || !(cfg.cluster_tb_time_us >= 0.0) ||
             !(cfg.cluster_pad_resolution_xy_mm > 0.0) ||
             !std::isfinite(cfg.cluster_coef_t_cm2_us) ||
             !std::isfinite(cfg.cluster_coef_l_cm2_us) ||
             !std::isfinite(cfg.cluster_drift_velocity_cm_us) ||
             !std::isfinite(cfg.cluster_tb_time_us) ||
             !std::isfinite(cfg.cluster_pad_resolution_xy_mm))) return false;
        A = a_in;
        Z = z_in;
        mass = mass_mev;
        bfield = b_t;
        efield = e_vm;

        // 触发表加载（缓存点云时这之前可能还没人用过质量表）
        (void)MassData::GetNuclearMassNuclear(1, 1);
        if (!MassData::IsLoaded()) {
            std::fprintf(stderr, "[KalmanFitter] 核素质量表未加载\n");
            return false;
        }
        const std::vector<Element> comp = gas.compound();
        if (comp.empty()) {
            std::fprintf(stderr, "[KalmanFitter] 靶材组成为空\n");
            return false;
        }

        dyn.ConfigureTarget(comp, gas.density());
        dyn.ConfigureFields(bfield, efield);
        dyn.SetParticle(A, Z);
        dyn.SetELossScale(1.0);
        dyn.SetStopZ(-10.0, 10.0);  // 边界交给滤波器自己判
        dyn.SetStopRho(10.0);

        prop.dyn = &dyn;
        prop.mass = mass;
        prop.step_mm = cfg.step_mm;
        prop.max_steps = cfg.max_steps;
        prop.gap_max_steps = cfg.gap_max_steps;
        prop.min_ke = cfg.min_ke_mev;
        prop.verbose = cfg.verbose;
        ukf.SetWeights(cfg.alpha, cfg.beta, cfg.kappa);
        ready = true;
        return true;
    }

    // 过程噪声方差：dE/dx 缩放的相对不确定度
    double ProcessNoiseVariance(const VecX& x, double path_mm) {
        path_mm=std::abs(path_mm);
        if (!cfg.straggling || path_mm == 0.0) return 1.0e-12; // 0 会让增广协方差奇异

        const double beta = BetaOfGv(GvFromP(x[3], mass));
        const double rho = dyn.GetDensity();          // g/cm^3
        const double rho_L = rho * (path_mm * 0.1);   // g/cm^2
        const double sigma2 =
            BohrStraggleSigma2(beta, std::abs(static_cast<double>(Z)), 1.0, rho_L);
        if (!(sigma2 > 0.0)) return 1.0e-12;

        const double ke = mass * (GammaOfGv(GvFromP(x[3], mass)) - 1.0);
        catima::Projectile proj(A, Z);
        // CATIMA's Projectile::T is kinetic energy per nucleon (MeV/u),
        // while the UKF state mass is in MeV/c^2.
        proj.T = ke / (mass / Units::u_MeV);
        const double dE_mev = catima::dedx(proj, dyn.GetTargetMaterial()) * rho_L;
        if (!(dE_mev > 0.0)) return 1.0e-12;

        double frac = std::sqrt(sigma2) / dE_mev;
        if (!(frac > 0.0) || !std::isfinite(frac)) return 1.0e-12;
        if (frac > cfg.straggle_cap) frac = cfg.straggle_cap;
        return frac * frac;
    }

    // 模型噪声：位置/动量/角度底噪，可选加 Highland 多重散射
    void ModelNoise(const VecX& x, double path_mm, MatX& qmod) {
        path_mm=std::abs(path_mm);
        qmod.setZero();
        qmod(0, 0) = cfg.q_pos;
        qmod(1, 1) = cfg.q_pos;
        qmod(2, 2) = cfg.q_pos;
        qmod(3, 3) = cfg.q_mom;
        qmod(4, 4) = cfg.q_ang;
        qmod(5, 5) = cfg.q_ang;

        if (!cfg.mult_scattering || cfg.rad_length_mm <= 0.0 || path_mm <= 0.0) return;

        const double p = x[3];
        const double beta = BetaOfGv(GvFromP(p, mass));
        if (!(p > 0.0) || !(beta > 0.0)) return;
        const double z = std::max(1.0, std::abs(static_cast<double>(Z)));
        const double LoverX0 = path_mm / cfg.rad_length_mm;
        if (!(LoverX0 > 0.0)) return;

        double theta0 = (13.6 / (beta * p)) * z * std::sqrt(LoverX0);
        const double logArg = LoverX0 * z * z / (beta * beta);
        if (logArg > 0.0) theta0 *= (1.0 + 0.038 * std::log(logArg));
        const double var = theta0 * theta0;
        if (!(var > 0.0) || !std::isfinite(var)) return;

        double sin2 = std::sin(x[4]);
        sin2 *= sin2;
        if (sin2 < 1.0e-4) sin2 = 1.0e-4;
        qmod(4, 4) += var;
        qmod(5, 5) += var / sin2;
    }

    bool PredictPath(const KalmanTrackState& prior, double path_mm,
                     const std::array<double,9>& measurement_covariance,
                     KalmanMeasurementPrediction& out) {
        out=KalmanMeasurementPrediction{};
        if(!ready||path_mm==0.0||!std::isfinite(path_mm)||
           (path_mm<0.0&&cfg.straggling)) return false;
        VecX x;
        MatX P;
        for(int i=0;i<NX;++i) x[i]=prior.mean[static_cast<size_t>(i)];
        for(int i=0;i<NX;++i) for(int j=0;j<NX;++j)
            P(i,j)=prior.covariance[static_cast<size_t>(i*NX+j)];
        if(!x.allFinite()||!P.allFinite()||!(x[3]>0.0)||x[4]<0.0||x[4]>M_PI)
            return false;
        P=0.5*(P+P.transpose());
        for(int i=0;i<NX;++i) if(P(i,i)<0.0) return false;

        MatZ R;
        for(int i=0;i<NZ;++i) for(int j=0;j<NZ;++j)
            R(i,j)=measurement_covariance[static_cast<size_t>(i*NZ+j)];
        if(!R.allFinite()) return false;
        R=0.5*(R+R.transpose());
        Eigen::LDLT<MatZ> rldlt(R);
        if(rldlt.info()!=Eigen::Success||!rldlt.isPositive()) return false;

        VecA xa;
        xa.head<NX>()=x;
        xa[NX]=1.0;
        MatA Pa=MatA::Zero();
        Pa.topLeftCorner<NX,NX>()=P;
        Pa(NX,NX)=ProcessNoiseVariance(x,path_mm);
        const MatAxNS sig=ukf.SigmaPoints(xa,Pa);
        const MatXNS sigPrev=sig.topRows<NX>();
        MatXNS sigProp;
        for(int i=0;i<NS;++i){
            VecX xi=kf::util::getColumnAt<NX,NS>(i,sigPrev);
            VecX xo;
            if(!prop.ByPath(xi,sig(NX,i),path_mm,xo)) return false;
            NormalizeAngles(xo);
            kf::util::copyToColumn<NX,NS>(i,sigProp,xo);
        }
        VecX meanPred;
        MatX Ppred;
        ukf.MeanCov<NX>(sigProp,meanPred,Ppred);
        MatX qmod;
        ModelNoise(x,path_mm,qmod);
        Ppred+=qmod;
        Ppred=0.5*(Ppred+Ppred.transpose());
        if(!meanPred.allFinite()||!Ppred.allFinite()||!(meanPred[3]>0.0)) return false;

        const MatZNS sigZ=sigProp.topRows<NZ>();
        VecZ zhat;
        MatZ Pzz;
        ukf.MeanCov<NZ>(sigZ,zhat,Pzz);
        Pzz+=R;
        Pzz=0.5*(Pzz+Pzz.transpose());
        Eigen::LDLT<MatZ> sldlt(Pzz);
        if(sldlt.info()!=Eigen::Success||!sldlt.isPositive()||!zhat.allFinite()) return false;
        const MatXZ Pxz=ukf.Cross<NZ>(sigProp,meanPred,sigZ,zhat);
        if(!Pxz.allFinite()) return false;

        for(int i=0;i<NX;++i){
            out.predicted_state.mean[static_cast<size_t>(i)]=meanPred[i];
            for(int j=0;j<NX;++j)
                out.predicted_state.covariance[static_cast<size_t>(i*NX+j)]=Ppred(i,j);
        }
        for(int i=0;i<NZ;++i){
            out.xyz_mean[static_cast<size_t>(i)]=zhat[i];
            for(int j=0;j<NZ;++j)
                out.innovation_covariance[static_cast<size_t>(i*NZ+j)]=Pzz(i,j);
        }
        for(int i=0;i<NX;++i) for(int j=0;j<NZ;++j)
            out.state_measurement_cross_covariance[static_cast<size_t>(i*NZ+j)]=Pxz(i,j);
        out.path_mm=path_mm;
        out.valid=true;
        return true;
    }

    bool Update(const KalmanMeasurementPrediction& prediction,
                const std::array<double,3>& measurement,
                KalmanMeasurementUpdate& out) const {
        out=KalmanMeasurementUpdate{};
        if(!ready||!prediction.valid) return false;
        VecX xp;
        MatX Pp;
        VecZ zhat,z;
        MatZ S;
        MatXZ Pxz;
        for(int i=0;i<NX;++i){
            xp[i]=prediction.predicted_state.mean[static_cast<size_t>(i)];
            for(int j=0;j<NX;++j)
                Pp(i,j)=prediction.predicted_state.covariance[static_cast<size_t>(i*NX+j)];
            for(int j=0;j<NZ;++j)
                Pxz(i,j)=prediction.state_measurement_cross_covariance[static_cast<size_t>(i*NZ+j)];
        }
        for(int i=0;i<NZ;++i){
            zhat[i]=prediction.xyz_mean[static_cast<size_t>(i)];
            z[i]=measurement[static_cast<size_t>(i)];
            for(int j=0;j<NZ;++j)
                S(i,j)=prediction.innovation_covariance[static_cast<size_t>(i*NZ+j)];
        }
        if(!xp.allFinite()||!Pp.allFinite()||!Pxz.allFinite()||!zhat.allFinite()||
           !z.allFinite()||!S.allFinite()) return false;
        Eigen::LDLT<MatZ> ldlt(S);
        if(ldlt.info()!=Eigen::Success||!ldlt.isPositive()) return false;
        const VecZ innov=z-zhat;
        const double chi2=innov.dot(ldlt.solve(innov));
        if(!std::isfinite(chi2)||chi2<0.0) return false;
        const MatXZ K=Pxz*ldlt.solve(MatZ::Identity());
        VecX corrected=xp+K*innov;
        NormalizeAngles(corrected);
        MatX Pcorrected=Pp-K*S*K.transpose();
        Pcorrected=0.5*(Pcorrected+Pcorrected.transpose());
        const bool physical=corrected.allFinite()&&Pcorrected.allFinite()&&
            corrected[3]>0.0&&corrected[4]>=0.0&&corrected[4]<=M_PI;
        if(!physical) return false;
        out.innovation_chi2=chi2;
        for(int i=0;i<NZ;++i) out.innovation[static_cast<size_t>(i)]=innov[i];
        const bool pass_gate=cfg.innovation_chi2_max<=0.0||
            !std::isfinite(cfg.innovation_chi2_max)||chi2<=cfg.innovation_chi2_max;
        const VecX& result=pass_gate?corrected:xp;
        const MatX& result_cov=pass_gate?Pcorrected:Pp;
        for(int i=0;i<NX;++i){
            out.filtered_state.mean[static_cast<size_t>(i)]=result[i];
            for(int j=0;j<NX;++j)
                out.filtered_state.covariance[static_cast<size_t>(i*NX+j)]=result_cov(i,j);
        }
        out.accepted=pass_gate;
        return true;
    }

    bool Fit(const std::vector<CloudPoint>& hits, const EstimateResult& guess,
             KalmanFitResult& out) {
        const int passes = std::clamp(cfg.iterations, 1, 4);
        std::vector<VecX> reference_states, next_reference_states;
        std::vector<MatX> reference_covariances, next_reference_covariances;
        KalmanFitResult first_result;
        bool first_ok = FitPass(hits, guess, first_result, nullptr, nullptr,
                                &reference_states, &reference_covariances);
        out = first_result;
        if (passes == 1 || reference_states.empty() || reference_covariances.empty())
            return first_ok;

        KalmanFitResult best_result = first_result;
        bool best_ok = first_ok;
        if (!first_ok && first_result.failure_code != 10) return false;
        double previous_ke = first_result.ke_vertex;
        double previous_polar = first_result.polar;
        for (int iteration = 1; iteration < passes; ++iteration) {
            KalmanFitResult candidate;
            next_reference_states.clear();
            next_reference_covariances.clear();
            const bool pass_ok = FitPass(hits, guess, candidate,
                                         &reference_states, &reference_covariances,
                                         &next_reference_states, &next_reference_covariances);
            if (!pass_ok || next_reference_states.empty() || next_reference_covariances.empty()) {
                if (!best_ok) out = candidate;
                break;
            }
            // A reference-track pass is only a refinement attempt. It must not
            // silently replace a fit with lower measurement coverage or worse
            // smoothed residual quality. ATTPCROOT can be configured to keep
            // iterating, but in this port the geometry/measurement ordering may
            // differ between inputs, so retain the best physically comparable
            // pass instead of blindly keeping the last one.
            const auto reduced_chi2 = [](const KalmanFitResult& result) {
                return result.ndf > 0 && std::isfinite(result.chi2)
                    ? result.chi2 / static_cast<double>(result.ndf)
                    : std::numeric_limits<double>::infinity();
            };
            const bool improves_best = !best_ok ||
                (candidate.n_hits >= best_result.n_hits &&
                 reduced_chi2(candidate) < reduced_chi2(best_result));
            if (improves_best) {
                best_result = candidate;
                best_ok = true;
            }
            const double dke = std::abs(candidate.ke_vertex - previous_ke) /
                std::max(1.0, std::abs(previous_ke));
            const double dpolar = std::abs(candidate.polar - previous_polar);
            if (best_ok) out = best_result;
            if (cfg.verbose)
                std::fprintf(stderr,
                             "[kalman] reference pass %d/%d KE=%.6g dKErel=%.3g dPolar=%.3g rad; %s pass (hits=%d redchi=%.4g)\n",
                             iteration + 1, passes, candidate.ke_vertex, dke, dpolar,
                             improves_best ? "keep" : "retain previous",
                             candidate.n_hits, reduced_chi2(candidate));
            reference_states.swap(next_reference_states);
            reference_covariances.swap(next_reference_covariances);
            if (dke < 1.0e-3 && dpolar < 1.0e-4) break;
            previous_ke = candidate.ke_vertex;
            previous_polar = candidate.polar;
        }
        if (best_ok) out = best_result;
        return best_ok;
    }

    bool FitPass(const std::vector<CloudPoint>& hits_in, const EstimateResult& guess,
                 KalmanFitResult& out, const std::vector<VecX>* reference_states,
                 const std::vector<MatX>* reference_covariances,
                 std::vector<VecX>* smoothed_reference_states,
                 std::vector<MatX>* smoothed_reference_covariances) {
        out = KalmanFitResult();
        if (smoothed_reference_states != nullptr) smoothed_reference_states->clear();
        if (smoothed_reference_covariances != nullptr) smoothed_reference_covariances->clear();
        if (!ready) { out.failure_code = 1; return false; }
        if (static_cast<int>(hits_in.size()) < cfg.min_clusters) {
            out.failure_code = 2;
            if (cfg.verbose) std::fprintf(stderr, "[kalman] rejected: only %zu raw hits\n", hits_in.size());
            return false;
        }

        // ---- 1) 排序：优先沿 Estimation 给出的运动方向 ----
        // Spyral Cluster 点云在 z 上有序；对近似横向轨迹仍应用 PCA 会把
        // 曲线投影折返，破坏 UKF 的测量次序。ROOT 原始点云则走后面的
        // Be12 最近邻重排，未带方向信息时才回退到 PCA。
        std::vector<CloudPoint> hits = hits_in;
        if (cfg.z_pad_plane_mm > 0.0) {
            for (auto& h : hits) h.z = cfg.z_pad_plane_mm - h.z;
        }
        auto rho_of = [](const CloudPoint& p) { return std::sqrt(p.x * p.x + p.y * p.y); };
        struct FitMeasurement {
            CloudPoint point;
            MatZ covariance = MatZ::Identity();
        };
        const double pad_resolution_xy = std::max(0.0, cfg.cluster_pad_resolution_xy_mm);
        const double pad_resolution_z = 1.5 * pad_resolution_xy;
        const double drift_velocity = cfg.cluster_drift_velocity_cm_us;
        const double tb_time = cfg.cluster_tb_time_us;
        bool direction_known = false;
        if (cfg.z_pad_plane_mm <= 0.0 && (guess.direction == 0 || guess.direction == 1)) {
            const bool forward = guess.direction == 0;
            std::stable_sort(hits.begin(), hits.end(), [&](const CloudPoint& a, const CloudPoint& b) {
                return forward ? a.z < b.z : a.z > b.z;
            });
            direction_known = true;
        } else if (cfg.z_pad_plane_mm <= 0.0) {
            double mx = 0.0, my = 0.0, mz = 0.0;
            for (const auto& p : hits) {
                mx += p.x;
                my += p.y;
                mz += p.z;
            }
            const double n = static_cast<double>(hits.size());
            mx /= n;
            my /= n;
            mz /= n;
            double cxx = 0.0, cxy = 0.0, cxz = 0.0, cyy = 0.0, cyz = 0.0, czz = 0.0;
            for (const auto& p : hits) {
                const double dx = p.x - mx, dy = p.y - my, dz = p.z - mz;
                cxx += dx * dx;
                cxy += dx * dy;
                cxz += dx * dz;
                cyy += dy * dy;
                cyz += dy * dz;
                czz += dz * dz;
            }
            double vx = cxx, vy = cxy, vz = cxz; // 从协方差第一行起步
            double vn = std::sqrt(vx * vx + vy * vy + vz * vz);
            if (vn <= 0.0) {
                vx = 0.0;
                vy = 0.0;
                vz = 1.0;
                vn = 1.0;
            }
            vx /= vn;
            vy /= vn;
            vz /= vn;
            for (int it = 0; it < 30; ++it) {
                double ax = cxx * vx + cxy * vy + cxz * vz;
                double ay = cxy * vx + cyy * vy + cyz * vz;
                double az = cxz * vx + cyz * vy + czz * vz;
                const double an = std::sqrt(ax * ax + ay * ay + az * az);
                if (!(an > 0.0)) break;
                vx = ax / an;
                vy = ay / an;
                vz = az / an;
            }
            std::sort(hits.begin(), hits.end(), [&](const CloudPoint& a, const CloudPoint& b) {
                const double ta = (a.x - mx) * vx + (a.y - my) * vy + (a.z - mz) * vz;
                const double tb = (b.x - mx) * vx + (b.y - my) * vy + (b.z - mz) * vz;
                return ta < tb;
            });
            // 未知方向时，用束轴端作为起点。
            if (rho_of(hits.front()) > rho_of(hits.back())) std::reverse(hits.begin(), hits.end());
        }

        // Cache weights and detector variances after sorting so each cached
        // entry remains aligned with its CloudPoint through the centroid pass.
        std::vector<double> hit_weights(hits.size(), 1.0);
        double total_peak_charge = 0.0;
        for (const auto& h : hits)
            if (std::isfinite(h.integral) && h.integral > 0.0) total_peak_charge += h.integral;
        if (total_peak_charge > 0.0) {
            for (size_t i = 0; i < hits.size(); ++i)
                hit_weights[i] = std::isfinite(hits[i].integral) && hits[i].integral > 0.0
                    ? hits[i].integral : 0.0;
        }
        std::vector<VecZ> per_hit_resolution_variance;
        if (cfg.per_cluster_covariance) {
            per_hit_resolution_variance.resize(hits.size());
            for (size_t hi = 0; hi < hits.size(); ++hi) {
                const double drift_z_mm = cfg.z_pad_plane_mm > 0.0
                    ? cfg.z_pad_plane_mm - hits[hi].z : hits[hi].z;
                const double drift_time_us = std::max(0.0, drift_z_mm) /
                    (10.0 * drift_velocity);
                const double var_t = 200.0 * cfg.cluster_coef_t_cm2_us * drift_time_us;
                const double var_l = 200.0 * cfg.cluster_coef_l_cm2_us * drift_time_us;
                const double tb_step_mm = 10.0 * drift_velocity * tb_time;
                const double var_tb = tb_step_mm * tb_step_mm / 12.0;
                per_hit_resolution_variance[hi] << pad_resolution_xy * pad_resolution_xy + var_t,
                    pad_resolution_xy * pad_resolution_xy + var_t,
                    pad_resolution_z * pad_resolution_z + var_tb + var_l;
            }
        }

        // ---- 1b) Be12 macro's adaptive ClusterizeSmooth3D measurements ----
        // It seeds a new cluster every 15 mm, averages all cloud points inside
        // a 20/25 mm sphere, then smooths at half radius around each centroid
        // and each adjacent midpoint. This preserves the curved local track
        // shape better than fixed equal-population bins.
        const double p_est = guess.brho * kBrStiffness * std::abs(Z) * 1000.0;
        const double ke_est = p_est > 0.0 ? std::sqrt(p_est * p_est + mass * mass) - mass : 0.0;
        const double radius = ke_est < 3.0 ? 25.0 : 20.0;
        constexpr double cluster_distance = 15.0;
        auto centroid_near = [&](const CloudPoint& center, double r, MatZ* covariance) {
            CloudPoint c = center;
            double sx = 0.0, sy = 0.0, sz = 0.0, sw = 0.0;
            double ux = 0.0, uy = 0.0, uz = 0.0;
            double var_x = 0.0, var_y = 0.0, var_z = 0.0;
            size_t count = 0;
            for (size_t hi = 0; hi < hits.size(); ++hi) {
                const auto& h = hits[hi];
                const double dx = h.x - center.x, dy = h.y - center.y, dz = h.z - center.z;
                if (dx * dx + dy * dy + dz * dz >= r * r) continue;
                const double w = hit_weights[hi];
                sx += w * h.x; sy += w * h.y; sz += h.z;
                ux += h.x; uy += h.y; uz += h.z; sw += w; ++count;
                if (covariance != nullptr && w > 0.0) {
                    // Direct port of ATTPCROOT AtTrackTransformer's per-hit
                    // pad + diffusion + time-bucket resolution, propagated to
                    // the charge-weighted centroid: sum(w_i^2 sigma_i^2)/(sum w_i)^2.
                    const double w2 = w * w;
                    var_x += w2 * per_hit_resolution_variance[hi][0];
                    var_y += w2 * per_hit_resolution_variance[hi][1];
                    var_z += w2 * per_hit_resolution_variance[hi][2];
                }
            }
            if (count == 0) return c;
            if (sw > 0.0) { c.x = sx / sw; c.y = sy / sw; c.z = sz / count; }
            else { c.x = ux / count; c.y = uy / count; c.z = uz / count; }
            c.integral = sw;
            if (covariance != nullptr) {
                covariance->setZero();
                constexpr double kMinClusterVarianceMm2 = 0.01; // ATTPCROOT's per-axis floor
                const double norm = sw * sw;
                if (norm > 0.0 && std::isfinite(norm)) {
                    (*covariance)(0, 0) = std::max(kMinClusterVarianceMm2, var_x / norm);
                    (*covariance)(1, 1) = std::max(kMinClusterVarianceMm2, var_y / norm);
                    (*covariance)(2, 2) = std::max(kMinClusterVarianceMm2, var_z / norm);
                } else {
                    covariance->diagonal().setConstant(cfg.meas_sigma_mm * cfg.meas_sigma_mm);
                }
            }
            return c;
        };
        auto has_near = [&](const CloudPoint& center, double r) {
            const double r2 = r * r;
            for (const auto& h : hits) {
                const double dx = h.x - center.x, dy = h.y - center.y, dz = h.z - center.z;
                if (dx * dx + dy * dy + dz * dz < r2) return true;
            }
            return false;
        };
        std::vector<CloudPoint> coarse;
        CloudPoint ref = hits.front();
        for (const auto& h : hits) {
            const double dx = h.x - ref.x, dy = h.y - ref.y, dz = h.z - ref.z;
            if (dx * dx + dy * dy + dz * dz < cluster_distance * cluster_distance) continue;
            if (!has_near(ref, radius)) { ref = h; continue; }
            const CloudPoint candidate = centroid_near(ref, radius, nullptr);
            bool separated = true;
            for (const auto& existing : coarse) {
                const double ex = candidate.x - existing.x;
                const double ey = candidate.y - existing.y;
                const double ez = candidate.z - existing.z;
                if (ex * ex + ey * ey + ez * ez < cluster_distance * cluster_distance) {
                    separated = false;
                    break;
                }
            }
            if (separated) coarse.push_back(candidate);
            ref = h;
        }
        if (coarse.size() < 2) coarse.push_back(centroid_near(hits.back(), radius, nullptr));

        std::vector<FitMeasurement> meas;
        if (coarse.size() > 2) {
            const double smooth_radius = radius * 0.5;
            for (size_t i = 0; i + 1 < coarse.size(); ++i) {
                if (has_near(coarse[i], smooth_radius)) {
                    FitMeasurement m;
                    m.point = centroid_near(coarse[i], smooth_radius,
                        cfg.per_cluster_covariance ? &m.covariance : nullptr);
                    meas.push_back(m);
                }
                CloudPoint mid = coarse[i];
                mid.x = 0.5 * (coarse[i].x + coarse[i + 1].x);
                mid.y = 0.5 * (coarse[i].y + coarse[i + 1].y);
                mid.z = 0.5 * (coarse[i].z + coarse[i + 1].z);
                if (has_near(mid, smooth_radius)) {
                    FitMeasurement m;
                    m.point = centroid_near(mid, smooth_radius,
                        cfg.per_cluster_covariance ? &m.covariance : nullptr);
                    meas.push_back(m);
                }
            }
            if (has_near(coarse.back(), smooth_radius)) {
                FitMeasurement m;
                m.point = centroid_near(coarse.back(), smooth_radius,
                    cfg.per_cluster_covariance ? &m.covariance : nullptr);
                meas.push_back(m);
            }
        } else {
            for (const auto& c : coarse) {
                FitMeasurement m;
                m.point = centroid_near(c, radius,
                    cfg.per_cluster_covariance ? &m.covariance : nullptr);
                meas.push_back(m);
            }
        }
        if (cfg.z_pad_plane_mm <= 0.0 && direction_known) {
            const bool forward = guess.direction == 0;
            std::stable_sort(meas.begin(), meas.end(), [&](const FitMeasurement& a, const FitMeasurement& b) {
                return forward ? a.point.z < b.point.z : a.point.z > b.point.z;
            });
        }
        if (cfg.z_pad_plane_mm > 0.0 && meas.size() > 1) {
            // ROOT's Be12 path starts at the largest digi-z cluster and then
            // walks to the nearest remaining cluster. Coordinates are already
            // reflected above, so the seed is minimum lab-z here.
            size_t seed = 0;
            for (size_t i = 1; i < meas.size(); ++i)
                if (meas[i].point.z < meas[seed].point.z) seed = i;
            std::vector<FitMeasurement> ordered;
            std::vector<uint8_t> used(meas.size(), 0);
            ordered.reserve(meas.size());
            ordered.push_back(meas[seed]);
            used[seed] = 1;
            while (ordered.size() < meas.size()) {
                const auto& prev = ordered.back().point;
                size_t best = meas.size();
                double best_d2 = std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < meas.size(); ++i) {
                    if (used[i]) continue;
                    const double dx = meas[i].point.x - prev.x;
                    const double dy = meas[i].point.y - prev.y;
                    const double dz = meas[i].point.z - prev.z;
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < best_d2) { best_d2 = d2; best = i; }
                }
                if (best == meas.size()) break;
                used[best] = 1;
                ordered.push_back(meas[best]);
            }
            meas.swap(ordered);
        }
        std::vector<CloudPoint> spaced;
        std::vector<MatZ> measurement_covariances;
        for (const auto& h : meas) {
            if (spaced.empty()) {
                spaced.push_back(h.point);
                measurement_covariances.push_back(h.covariance);
                continue;
            }
            const auto& prev = spaced.back();
            const double dx = h.point.x - prev.x, dy = h.point.y - prev.y, dz = h.point.z - prev.z;
            if (dx * dx + dy * dy + dz * dz >= 9.0) {
                spaced.push_back(h.point);
                measurement_covariances.push_back(h.covariance);
            }
        }
        const size_t n_keep = std::max(static_cast<size_t>(cfg.min_clusters),
                                       static_cast<size_t>(spaced.size() * 0.9));
        if (spaced.size() > n_keep) {
            spaced.resize(n_keep);
            measurement_covariances.resize(n_keep);
        }
        if (cfg.verbose) {
            std::fprintf(stderr, "[kalman] Be12 clusters coarse=%zu smooth=%zu spaced=%zu final=%zu\n",
                         coarse.size(), meas.size(), spaced.size(), n_keep);
            for (size_t i = 0; i < spaced.size(); ++i)
                std::fprintf(stderr, "[kalman] meas[%zu]=%.6f,%.6f,%.6f\n",
                             i, spaced[i].x, spaced[i].y, spaced[i].z);
        }
        hits.swap(spaced);
        int n = static_cast<int>(hits.size());
        if (n < cfg.min_clusters) {
            out.failure_code = 3;
            out.n_hits = n;
            if (cfg.verbose) std::fprintf(stderr, "[kalman] rejected: only %d smoothed measurements\n", n);
            return false;
        }

        // ---- 2) 初始状态：动量大小来自 Estimator 的 Bρ，方向用点序列 ----
        double p0 = 0.0;
        if (guess.brho > 0.0) p0 = guess.brho * kBrStiffness * std::abs(Z) * 1000.0; // MeV/c
        if (!(p0 > 0.0)) {
            out.failure_code = 4;
            if (cfg.verbose) std::fprintf(stderr, "[kalman] rejected: invalid momentum seed %.6g\n", p0);
            return false;
        }

        // A z-sorted cluster can begin with an isolated point from another
        // turn or a merged fragment. Seed from the estimator vertex and keep
        // the first early measurement that the physical trajectory actually
        // reaches. This avoids letting one bad endpoint define the UKF state.
        VecX transported_seed;
        bool have_transported_seed = false;
        int skipped_seed_points = 0;
        if (cfg.z_pad_plane_mm <= 0.0 && std::isfinite(guess.vertex_x) &&
            std::isfinite(guess.vertex_y) && std::isfinite(guess.vertex_z) &&
            std::isfinite(guess.polar) && std::isfinite(guess.azimuthal)) {
            VecX at_vertex;
            at_vertex << guess.vertex_x, guess.vertex_y, guess.vertex_z,
                        p0, guess.polar, guess.azimuthal;
            const double seed_gate_mm = std::max(20.0, 5.0 * cfg.meas_sigma_mm);
            const int seed_search = std::min(n, 12);
            for (int j = 0; j < seed_search; ++j) {
                const double point[3] = {hits[static_cast<size_t>(j)].x * 1.0e-3,
                                         hits[static_cast<size_t>(j)].y * 1.0e-3,
                                         hits[static_cast<size_t>(j)].z * 1.0e-3};
                VecX candidate;
                if (!prop.ToPoint(at_vertex, 1.0, point, candidate)) continue;
                const double dx = candidate[0] - hits[static_cast<size_t>(j)].x;
                const double dy = candidate[1] - hits[static_cast<size_t>(j)].y;
                const double dz = candidate[2] - hits[static_cast<size_t>(j)].z;
                const double residual_mm = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (residual_mm <= seed_gate_mm) {
                    skipped_seed_points = j;
                    transported_seed = candidate;
                    have_transported_seed = true;
                    break;
                }
            }
        }
        if (have_transported_seed && skipped_seed_points > 0) {
            hits.erase(hits.begin(), hits.begin() + skipped_seed_points);
            measurement_covariances.erase(measurement_covariances.begin(),
                                          measurement_covariances.begin() + skipped_seed_points);
            n = static_cast<int>(hits.size());
            if (cfg.verbose)
                std::fprintf(stderr, "[kalman] skipped %d unreachable leading measurements; seeded at first physical hit\n",
                             skipped_seed_points);
        }
        if (n < cfg.min_clusters) {
            out.failure_code = 5;
            out.n_hits = n;
            if (cfg.verbose) std::fprintf(stderr, "[kalman] rejected: only %d measurements after seed selection\n", n);
            return false;
        }

        const CloudPoint& h0 = hits.front();
        // Direction seed: Be12 ROOT's z-pad-plane path uses the chord between
        // its first two beam-axis-sorted clusters. Keep PCA for other inputs.
        double dir[3] = {0.0, 0.0, 1.0};
        if (cfg.z_pad_plane_mm > 0.0 && n >= 2) {
            dir[0] = hits[1].x - hits[0].x;
            dir[1] = hits[1].y - hits[0].y;
            dir[2] = hits[1].z - hits[0].z;
            const double norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
            if (norm > 0.0) for (double& component : dir) component /= norm;
        } else {
            const int K = std::min(n, 8);
            double mx = 0.0, my = 0.0, mz = 0.0;
            for (int i = 0; i < K; ++i) {
                mx += hits[static_cast<size_t>(i)].x;
                my += hits[static_cast<size_t>(i)].y;
                mz += hits[static_cast<size_t>(i)].z;
            }
            mx /= K;
            my /= K;
            mz /= K;
            double cxx = 0, cxy = 0, cxz = 0, cyy = 0, cyz = 0, czz = 0;
            for (int i = 0; i < K; ++i) {
                const double dx = hits[static_cast<size_t>(i)].x - mx;
                const double dy = hits[static_cast<size_t>(i)].y - my;
                const double dz = hits[static_cast<size_t>(i)].z - mz;
                cxx += dx * dx;
                cxy += dx * dy;
                cxz += dx * dz;
                cyy += dy * dy;
                cyz += dy * dz;
                czz += dz * dz;
            }
            double vx = cxx, vy = cxy, vz = cxz;
            double vn = std::sqrt(vx * vx + vy * vy + vz * vz);
            if (vn > 0.0) {
                vx /= vn;
                vy /= vn;
                vz /= vn;
                for (int it = 0; it < 20; ++it) {
                    const double ax = cxx * vx + cxy * vy + cxz * vz;
                    const double ay = cxy * vx + cyy * vy + cyz * vz;
                    const double az = cxz * vx + cyz * vy + czz * vz;
                    const double an = std::sqrt(ax * ax + ay * ay + az * az);
                    if (!(an > 0.0)) break;
                    vx = ax / an;
                    vy = ay / an;
                    vz = az / an;
                }
                // 符号跟随 hits 的排序方向（可能是反向径迹，θ > 90°），
                // 不能强行取 vz > 0，否则反冲径迹永远拟不上。
                const CloudPoint& last = hits[static_cast<size_t>(n - 1)];
                const double chord[3] = {last.x - h0.x, last.y - h0.y, last.z - h0.z};
                if (vx * chord[0] + vy * chord[1] + vz * chord[2] < 0.0) {
                    vx = -vx;
                    vy = -vy;
                    vz = -vz;
                }
                dir[0] = vx;
                dir[1] = vy;
                dir[2] = vz;
            }
        }

        VecX x;
        x[0] = h0.x;
        x[1] = h0.y;
        x[2] = h0.z;
        x[3] = p0;
        // ATTPCROOT's Be12 macro sets ZPadPlane=1000, which flips the
        // coordinates to lab frame and explicitly replaces the PRA angles
        // with the chord from the first two ordered clusters. Match that
        // seed here; the circle-fit estimate is useful for momentum, but its
        // angle can point to the opposite half-sphere on these tracks.
        if (n >= 2) {
            dir[0] = hits[1].x - h0.x;
            dir[1] = hits[1].y - h0.y;
            dir[2] = hits[1].z - h0.z;
            const double dn = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
            if (dn > 0.0) {
                dir[0] /= dn;
                dir[1] /= dn;
                dir[2] /= dn;
            }
        }
        x[4] = std::acos(std::max(-1.0, std::min(1.0, dir[2])));
        x[5] = std::atan2(dir[1], dir[0]);
        if (cfg.z_pad_plane_mm <= 0.0 && n >= 2 &&
            std::isfinite(guess.polar) && std::isfinite(guess.center_x) &&
            std::isfinite(guess.center_y) && !have_transported_seed) {
            // The first two z-sorted measurements can come from different
            // turns on a tightly curling track. Use the fitted circle center
            // to get the local tangent at h0; use the global polar estimate
            // because the first chord's small z increment is also noisy.
            const double rx = h0.x - guess.center_x;
            const double ry = h0.y - guess.center_y;
            const double rr = std::hypot(rx, ry);
            if (rr > 1.0) {
                const double sense = bfield < 0.0 ? 1.0 : -1.0;
                x[4] = guess.polar;
                x[5] = std::atan2(sense * rx, -sense * ry);
            } else {
                x[4] = guess.polar;
                x[5] = std::atan2(hits[1].y - h0.y, hits[1].x - h0.x);
            }
        }
        if (have_transported_seed) {
            x[3] = transported_seed[3];
            x[4] = transported_seed[4];
            x[5] = transported_seed[5];
        }
        const VecX estimator_seed = x;
        bool reference_seed_valid = reference_states != nullptr && reference_covariances != nullptr &&
            reference_states->size() == static_cast<size_t>(n) &&
            reference_covariances->size() == static_cast<size_t>(n) &&
            (*reference_states)[0].allFinite() && (*reference_states)[0][3] > 0.0;
        if (reference_seed_valid) {
            x = (*reference_states)[0];
            NormalizeAngles(x);
        }

        MatX P = MatX::Zero();
        const double ps = cfg.pos_sigma_mm;
        P(0, 0) = ps * ps;
        P(1, 1) = ps * ps;
        P(2, 2) = ps * ps;
        const double pm = cfg.mom_sigma_frac * x[3];
        P(3, 3) = pm * pm;
        const double as = cfg.ang_sigma_deg * Units::DEG2RAD;
        P(4, 4) = as * as;
        P(5, 5) = as * as;
        if (reference_seed_valid) {
            const MatX& ref_cov = (*reference_covariances)[0];
            bool covariance_ok = ref_cov.allFinite();
            for (int d = 0; d < NX && covariance_ok; ++d)
                covariance_ok = ref_cov(d, d) >= 0.0 && ref_cov(d, d) <= 1.0e10;
            if (covariance_ok) {
                P.setZero(); // Drop near-singular correlations as ATTPCROOT does.
                for (int d = 0; d < NX; ++d)
                    P(d, d) = ref_cov(d, d) * cfg.ref_covariance_inflation;
                const double pos_floor = cfg.pos_sigma_mm * cfg.pos_sigma_mm;
                const double mom_floor = std::pow(cfg.mom_sigma_frac * x[3], 2);
                const double ang_floor = std::pow(0.5 * Units::DEG2RAD, 2);
                for (int d = 0; d < 3; ++d) P(d, d) = std::max(P(d, d), pos_floor);
                P(3, 3) = std::max(P(3, 3), mom_floor);
                P(4, 4) = std::max(P(4, 4), ang_floor);
                P(5, 5) = std::max(P(5, 5), ang_floor);
            } else {
                reference_seed_valid = false;
                x = estimator_seed;
                P.setZero();
                P(0, 0) = ps * ps;
                P(1, 1) = ps * ps;
                P(2, 2) = ps * ps;
                P(3, 3) = std::pow(cfg.mom_sigma_frac * x[3], 2);
                P(4, 4) = as * as;
                P(5, 5) = as * as;
                if (cfg.verbose)
                    std::fprintf(stderr, "[kalman] reference seed covariance invalid; use estimator seed\n");
            }
        }

        const MatZ uniform_measurement_covariance =
            MatZ::Identity() * (cfg.meas_sigma_mm * cfg.meas_sigma_mm);

        // ---- 3) 前向滤波（每步：sigma 点传播 → 预测 → 校正）----
        std::vector<VecX> xFilt, xPred;
        std::vector<MatX> pFilt, pPred, cHist;
        std::vector<VecZ> measurements_used;
        std::vector<int> accepted_measurement_at_step;
        std::vector<int> state_hit_index;
        xFilt.push_back(x);
        pFilt.push_back(P);
        state_hit_index.push_back(0);

        double chi2 = 0.0;
        int n_used = 0;
        int n_steps = 0;
        int n_skip = 0;
        bool gap_pending = false;

        for (int k = 1; k < n; ++k) {
            double nrm[3];
            UnitDir(x, nrm);
            const double point[3] = {hits[k].x * 1.0e-3, hits[k].y * 1.0e-3,
                                     hits[k].z * 1.0e-3};

            // 先走一遍均值状态，拿真实路径长度（给过程噪声用）
            VecX mean_prop;
            double path_mm = 0.0;
            const double measured_gap_mm = std::abs(hits[k].z - hits[k - 1].z);
            const bool continuing_gap = gap_pending;
            bool use_z_plane = direction_known && cfg.z_pad_plane_mm <= 0.0 &&
                               (gap_pending || measured_gap_mm > 80.0);
            if (cfg.verbose && use_z_plane)
                std::fprintf(stderr, "[kalman] bridge uninstrumented z gap k=%d/%d gap=%.1f mm theta=%.2f deg\n",
                             k, n, measured_gap_mm, x[4] * Units::RAD2DEG);
            bool mean_ok = use_z_plane
                ? prop.ToZPlane(x, 1.0, point, mean_prop, &path_mm)
                : prop.ToPoint(x, 1.0, point, mean_prop, &path_mm);
            // On a curling track the tangent-plane closest-approach test can
            // select a later helix turn even when the measured z is only a
            // few millimetres ahead. For direction-known Spyral clusters,
            // z is monotonic along the observed sequence, so recover the
            // correct turn by retrying at the unique measured-z plane.
            constexpr double kWrongTurnZToleranceMm = 20.0;
            const double wrong_turn_dz = mean_ok ? mean_prop[2] - hits[k].z : 0.0;
            if (mean_ok && !use_z_plane && direction_known &&
                std::abs(wrong_turn_dz) > kWrongTurnZToleranceMm) {
                use_z_plane = true;
                mean_ok = prop.ToZPlane(x, 1.0, point, mean_prop, &path_mm,
                                        std::abs(wrong_turn_dz) > 100.0);
                if (cfg.verbose)
                    std::fprintf(stderr,
                                 "[kalman] retry wrong-turn closest approach at z plane k=%d dz=%.1f mm\n",
                                 k, wrong_turn_dz);
            }
            if (!mean_ok) {
                if (use_z_plane) gap_pending = true;
                // 真实数据里簇又碎又噪，个别量测点"落在身后"是常态：
                // TPC pad 的无读出区会造成长缺口，传播上限已提高到足以跨过
                // 多圈螺旋；若扩展搜索后仍有连续不可达点，则视为错序/停止尾段。
                if (++n_skip > 3) {
                    if (cfg.verbose)
                        std::fprintf(stderr, "[kalman] 扩展传播后仍连续不可达，截取可用轨迹段（k=%d/%d, 已用=%d）\n", k, n, n_used);
                    break;
                }
                continue;
            }
            gap_pending = false;
            n_skip = 0;

            VecX sigma_center = x;
            VecX plane_mean = mean_prop;
            double sigma_path_mm = path_mm;
            bool use_reference_step = false;
            if (reference_seed_valid && reference_states != nullptr &&
                reference_covariances != nullptr && k - 1 < static_cast<int>(reference_states->size())) {
                const VecX& ref_state = (*reference_states)[static_cast<size_t>(k - 1)];
                bool ref_valid = ref_state.allFinite() && ref_state[3] > 0.0 &&
                                 ref_state[3] < 1.0e6;
                if (ref_valid) {
                    VecX ref_prop;
                    double ref_path_mm = 0.0;
                    bool ref_ok = use_z_plane
                        ? prop.ToZPlane(ref_state, 1.0, point, ref_prop, &ref_path_mm)
                        : prop.ToPoint(ref_state, 1.0, point, ref_prop, &ref_path_mm);
                    const double ref_wrong_turn_dz = ref_ok ? ref_prop[2] - hits[k].z : 0.0;
                    if (ref_ok && !use_z_plane && direction_known &&
                        std::abs(ref_wrong_turn_dz) > 20.0) {
                        ref_ok = prop.ToZPlane(ref_state, 1.0, point, ref_prop,
                                               &ref_path_mm, std::abs(ref_wrong_turn_dz) > 100.0);
                    }
                    if (ref_ok && ref_prop.allFinite() && ref_prop[3] > 0.0) {
                        sigma_center = ref_state;
                        plane_mean = ref_prop;
                        sigma_path_mm = ref_path_mm;
                        use_reference_step = true;
                    }
                }
                if (!use_reference_step && cfg.verbose)
                    std::fprintf(stderr,
                                 "[kalman] reference state unavailable at k=%d; use current filtered state\n", k);
            }
            const double qv = ProcessNoiseVariance(sigma_center, sigma_path_mm);
            double plane_normal[3];
            UnitDir(plane_mean, plane_normal);
            VecA xa;
            xa.head<NX>() = sigma_center;
            xa[NX] = 1.0; // dE/dx 缩放因子均值 = 1
            MatA Pa = MatA::Zero();
            Pa.topLeftCorner<NX, NX>() = P;
            Pa(NX, NX) = qv;
            const MatAxNS sig = ukf.SigmaPoints(xa, Pa);
            const MatXNS sigPrev = sig.topRows<NX>();

            MatXNS sigProp;
            VecX meanPred;
            MatX Ppred;
            MatX qmod;
            ModelNoise(sigma_center, sigma_path_mm, qmod);
            bool locally_linearized = false;

            for (int i = 0; i < NS; ++i) {
                VecX xi = kf::util::getColumnAt<NX, NS>(i, sigPrev);
                const double v = sig(NX, i);
                VecX xo;
                // Match ATTPCROOT's UKF transition: the mean determines one
                // measurement plane and every sigma point is propagated to
                // that same plane. On a curled-track branch, first prefer the
                // mean's path length so sigma points stay on the same turn.
                // If that branch is unreachable or crosses a different turn
                // (detected by its z displacement), fall back to the unique z
                // plane for this z-ordered measurement instead of injecting an
                // unpropagated sigma point into the weighted transform.
                const double preferred_path_mm = use_z_plane ? sigma_path_mm : -1.0;
                bool sigma_ok = prop.ToPlane(xi, v, point, plane_normal, xo,
                                             preferred_path_mm);
                constexpr double kSigmaWrongTurnZToleranceMm = 20.0;
                const bool wrong_turn_crossing = sigma_ok && use_z_plane &&
                    std::abs(xo[2] - hits[k].z) > kSigmaWrongTurnZToleranceMm;
                if (use_z_plane && (!sigma_ok || wrong_turn_crossing))
                    sigma_ok = prop.ToZPlane(xi, v, point, xo, nullptr, true);
                if (!sigma_ok || (use_z_plane &&
                    std::abs(xo[2] - hits[k].z) > kSigmaWrongTurnZToleranceMm)) {
                    if (cfg.verbose)
                        std::fprintf(stderr,
                                     "[kalman] sigma branch ambiguous at k=%d sigma=%d z=%.1f mm; use local linear transition\n",
                                     k, i, hits[k].z);
                    locally_linearized = true;
                    break;
                }
                kf::util::copyToColumn<NX, NS>(i, sigProp, xo);
            }
            if (locally_linearized) {
                // The UT assumes one smooth transition function. If one sigma
                // point crosses a different helix turn (or cannot reach the
                // plane), mixing its unrelated state into the weighted moments
                // is invalid. Use the first-order transition around the
                // physically propagated mean for the complete sigma set.
                // This keeps the measurement in the fit; residual, coverage,
                // and final chi2 gates still assess whether it belongs there.
                for (int i = 0; i < NS; ++i) {
                    VecX xi = kf::util::getColumnAt<NX, NS>(i, sigPrev);
                    NormalizeAngles(xi);
                    VecX delta = xi - x;
                    delta[5] = std::remainder(delta[5], 2.0 * M_PI);
                    VecX xo = plane_mean + delta;
                    NormalizeAngles(xo);
                    kf::util::copyToColumn<NX, NS>(i, sigProp, xo);
                }
            }
            ukf.MeanCov<NX>(sigProp, meanPred, Ppred);
            Ppred += qmod;
            Ppred = 0.5 * (Ppred + Ppred.transpose());
            if (use_reference_step) {
                // Reference-track linearization stabilizes the covariance and
                // smoother cross-correlation. Preserve the exact nonlinear
                // propagation of the running filtered mean by recentering the
                // propagated reference sigma ensemble on mean_prop.
                const VecX shift = mean_prop - meanPred;
                for (int i = 0; i < NS; ++i) {
                    VecX yi = kf::util::getColumnAt<NX, NS>(i, sigProp) + shift;
                    kf::util::copyToColumn<NX, NS>(i, sigProp, yi);
                }
                meanPred = mean_prop;
            }

            // ATTPCROOT validates the predicted state and covariance before
            // applying a measurement. Here a non-finite or runaway sigma
            // transform must reject the candidate before it contaminates the
            // remaining hits and RTS smoother.
            bool prediction_valid = meanPred.allFinite() && Ppred.allFinite();
            for (int d = 0; d < NX && prediction_valid; ++d) {
                prediction_valid = std::abs(meanPred[d]) <= 1.0e6 &&
                                   Ppred(d, d) >= 0.0 && Ppred(d, d) <= 1.0e10;
            }
            prediction_valid = prediction_valid && meanPred[3] > 0.0 &&
                               meanPred[3] < 1000.0 && meanPred[4] >= 0.0 &&
                               meanPred[4] <= M_PI;
            if (!prediction_valid) {
                if (cfg.verbose)
                    std::fprintf(stderr,
                                 "[kalman] rejected: invalid predicted state/covariance at k=%d p=%.6g theta=%.6g\n",
                                 k, meanPred[3], meanPred[4]);
                out.failure_code = 6;
                out.n_hits = n_used + 1;
                return false;
            }

            cHist.push_back(ukf.Cross<NX>(sigPrev, sigma_center, sigProp, meanPred));

            xPred.push_back(meanPred);
            pPred.push_back(Ppred);

            // 校正
            const MatZ& measurement_covariance = cfg.per_cluster_covariance
                ? measurement_covariances[static_cast<size_t>(k)]
                : uniform_measurement_covariance;
            VecZ zhat;
            MatZ Pzz;
            const MatZNS sigZ = sigProp.topRows<NZ>();
            ukf.MeanCov<NZ>(sigZ, zhat, Pzz);
            Pzz += measurement_covariance;
            Pzz = 0.5 * (Pzz + Pzz.transpose());

            VecZ z;
            z[0] = hits[k].x;
            z[1] = hits[k].y;
            z[2] = hits[k].z;

            Eigen::LDLT<MatZ> ldlt(Pzz);
            if (ldlt.info() != Eigen::Success) {
                out.failure_code = 6;
                if (cfg.verbose) std::fprintf(stderr, "[kalman] rejected: singular innovation covariance at k=%d\n", k);
                return false;
            }
            const MatXZ Pxz = ukf.Cross<NZ>(sigProp, meanPred, sigZ, zhat);
            const MatXZ K = Pxz * ldlt.solve(MatZ::Identity());

            const VecZ innov = z - zhat;
            const double innovation_chi2 = innov.dot(ldlt.solve(innov));
            VecX corrected = meanPred + K * innov;
            NormalizeAngles(corrected);
            const bool physical_update = std::isfinite(corrected[3]) && corrected[3] > 0.0 &&
                                         std::isfinite(corrected[4]) && std::isfinite(corrected[5]);
            const bool accept_measurement = physical_update &&
                (!std::isfinite(cfg.innovation_chi2_max) || cfg.innovation_chi2_max <= 0.0 ||
                 innovation_chi2 <= cfg.innovation_chi2_max);
            if (accept_measurement) {
                x = corrected;
                P = Ppred - K * Pzz * K.transpose();
                P = 0.5 * (P + P.transpose());
                bool filtered_valid = x.allFinite() && P.allFinite();
                for (int d = 0; d < NX && filtered_valid; ++d)
                    filtered_valid = P(d, d) >= 0.0 && P(d, d) <= 1.0e10;
                filtered_valid = filtered_valid && std::abs(x[0]) <= 1.0e6 &&
                                 std::abs(x[1]) <= 1.0e6 && std::abs(x[2]) <= 1.0e6 &&
                                 x[3] > 0.0 && x[3] < 1000.0 &&
                                 x[4] >= 0.0 && x[4] <= M_PI;
                if (!filtered_valid) {
                    out.failure_code = 6;
                    out.n_hits = n_used + 1;
                    if (cfg.verbose)
                        std::fprintf(stderr,
                                     "[kalman] rejected: invalid filtered state/covariance at k=%d p=%.6g theta=%.6g\n",
                                     k, x[3], x[4]);
                    return false;
                }
                measurements_used.push_back(z);
                accepted_measurement_at_step.push_back(n_used);
                ++n_used;
            } else {
                // Keep the deterministic physical mean when the sigma-point
                // mean is nonphysical, but retain the propagated covariance.
                // Reverting to P and adding only Q here prevents uncertainty
                // from growing through a missing-pad interval, so a genuine
                // reappearing segment can never pass the innovation gate.
                // The UKF weighted mean can become nonphysical when sigma
                // trajectories land on different turns. In that case keep
                // the deterministic physical propagation while preserving
                // the sigma-point covariance transport.
                x = mean_prop;
                P = use_z_plane ? Ppred : P + qmod;
                P = 0.5 * (P + P.transpose());
                xPred.back() = mean_prop;
                pPred.back() = P;
                if (use_z_plane) cHist.back().setZero();
                accepted_measurement_at_step.push_back(-1);
                if (cfg.verbose)
                    std::fprintf(stderr,
                                 "[kalman] gated k=%d innovation_chi2=%.2f physical=%d pred=(%.1f,%.1f,%.1f) resid=(%.1f,%.1f,%.1f); continue propagation\n",
                                 k, innovation_chi2, physical_update ? 1 : 0,
                                 zhat[0], zhat[1], zhat[2], innov[0], innov[1], innov[2]);
            }

            if (cfg.verbose && accept_measurement)
                std::fprintf(stderr,
                             "[kalman] update k=%d z=%.1f pred=(%.1f,%.1f,%.1f) resid=(%.1f,%.1f,%.1f) p=%.2f->%.2f theta=%.2f\n",
                             k, hits[k].z, zhat[0], zhat[1], zhat[2], innov[0], innov[1], innov[2],
                             meanPred[3], x[3], x[4] * Units::RAD2DEG);

            ++n_steps;

            xFilt.push_back(x);
            pFilt.push_back(P);
            state_hit_index.push_back(k);
        }

        // A partial fit ending before a known dead-pad gap can have a small
        // chi-square simply because none of the reappearing segment was used.
        // Do not report that as a successful full-track Kalman fit.
        if (gap_pending) {
            out.failure_code = 7;
            out.n_hits = n_used + 1;
            if (cfg.verbose)
                std::fprintf(stderr, "[kalman] rejected incomplete track: resumed segment remained unreachable\n");
            return false;
        }
        if (n_used < 2) {
            out.failure_code = 8;
            out.n_hits = n_used + 1;
            if (cfg.verbose) std::fprintf(stderr, "[kalman] rejected: only %d accepted measurements\n", n_used);
            return false;
        }
        const int minimum_measurements = std::max(
            std::max(1, cfg.min_clusters - 1),
            static_cast<int>(std::ceil(std::max(0.0, cfg.min_measurement_fraction) * (n - 1))));
        if (n_used < minimum_measurements) {
            out.failure_code = 9;
            out.n_hits = n_used + 1;
            if (cfg.verbose)
                std::fprintf(stderr,
                             "[kalman] rejected incomplete fit: used %d/%d measurements, require %d\n",
                             n_used, n - 1, minimum_measurements);
            return false;
        }

        // ---- 4) RTS 平滑（索引：xPred/pPred/cHist 的第 j 项对应第 j+1 个量测）----
        std::vector<VecX> xS = xFilt;
        std::vector<MatX> pS = pFilt;
        if (cfg.smooth) {
            for (int k = n_steps; k >= 1; --k) {
                const MatX& PpredK = pPred[static_cast<size_t>(k - 1)];
                Eigen::LLT<MatX> llt(PpredK);
                if (llt.info() != Eigen::Success) continue;
                const MatX D = llt.solve(cHist[static_cast<size_t>(k - 1)].transpose()).transpose();
                VecX delta = xS[static_cast<size_t>(k)] - xPred[static_cast<size_t>(k - 1)];
                delta[5] = std::remainder(delta[5], 2.0 * M_PI);
                // (theta, phi) and (pi-theta, phi+pi) encode the same
                // direction. Use the nearer representation before applying
                // the linear RTS correction across angle boundaries.
                const double theta_alt = M_PI - xS[static_cast<size_t>(k)][4];
                const double phi_alt = xS[static_cast<size_t>(k)][5] + M_PI;
                const double dtheta_alt = theta_alt - xPred[static_cast<size_t>(k - 1)][4];
                const double dphi_alt = std::remainder(
                    phi_alt - xPred[static_cast<size_t>(k - 1)][5], 2.0 * M_PI);
                const double direct_angle = delta[4] * delta[4] +
                    std::pow(std::sin(xPred[static_cast<size_t>(k - 1)][4]) * delta[5], 2);
                const double alternate_angle = dtheta_alt * dtheta_alt +
                    std::pow(std::sin(xPred[static_cast<size_t>(k - 1)][4]) * dphi_alt, 2);
                if (alternate_angle < direct_angle) {
                    delta[4] = dtheta_alt;
                    delta[5] = dphi_alt;
                }
                xS[static_cast<size_t>(k - 1)] =
                    xFilt[static_cast<size_t>(k - 1)] + D * delta;
                NormalizeAngles(xS[static_cast<size_t>(k - 1)]);
                MatX Pk = pFilt[static_cast<size_t>(k - 1)] +
                          D * (pS[static_cast<size_t>(k)] - PpredK) * D.transpose();
                pS[static_cast<size_t>(k - 1)] = 0.5 * (Pk + Pk.transpose());
            }
        }

        if (smoothed_reference_states != nullptr && smoothed_reference_covariances != nullptr) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            smoothed_reference_states->assign(static_cast<size_t>(n), VecX::Constant(nan));
            smoothed_reference_covariances->assign(static_cast<size_t>(n), MatX::Constant(nan));
            for (size_t i = 0; i < state_hit_index.size() && i < xS.size() && i < pS.size(); ++i) {
                const int hit_index = state_hit_index[i];
                if (hit_index < 0 || hit_index >= n) continue;
                (*smoothed_reference_states)[static_cast<size_t>(hit_index)] = xS[i];
                (*smoothed_reference_covariances)[static_cast<size_t>(hit_index)] = pS[i];
            }
        }

        // ---- 5) 结果 ----
        const VecX x0 = xS.front();
        const MatX P0 = pS.front();

        double sum_d = 0.0;
        int n_cmp = 0;
        const double initial_dx = xS.front()[0] - h0.x;
        const double initial_dy = xS.front()[1] - h0.y;
        const double initial_dz = xS.front()[2] - h0.z;
        sum_d += std::sqrt(initial_dx * initial_dx + initial_dy * initial_dy + initial_dz * initial_dz);
        ++n_cmp;
        chi2 = 0.0;
        for (int step = 0; step < n_steps; ++step) {
            const int measurement_index = accepted_measurement_at_step[static_cast<size_t>(step)];
            if (measurement_index < 0) continue;
            const VecX& xs = xS[static_cast<size_t>(step + 1)];
            const VecZ& measurement = measurements_used[static_cast<size_t>(measurement_index)];
            const double mx = measurement[0];
            const double my = measurement[1];
            const double mz = measurement[2];
            const double dx = xs[0] - mx;
            const double dy = xs[1] - my;
            const double dz = xs[2] - mz;
            sum_d += std::sqrt(dx * dx + dy * dy + dz * dz);
            ++n_cmp;
            // ATTPCROOT optionally uses per-cluster R during the filter update,
            // but its final quality gate still uses the configured global
            // measurement sigma. Preserve that convention for ROOT comparisons.
            const double sigma2 = cfg.meas_sigma_mm * cfg.meas_sigma_mm;
            chi2 += (dx * dx + dy * dy + dz * dz) / sigma2;
        }
        out.cost = (n_cmp > 0) ? sum_d / static_cast<double>(n_cmp) : 0.0;
        // Match ATTPCROOT's fit metadata: chi2 is computed from the smoothed
        // positions, with the configured measurement sigma, not from filter
        // innovations (which include the predicted-state covariance).
        out.chi2 = chi2;
        out.ndf = 3 * n_used - 6;
        out.n_hits = n_used + 1;

        // Be12 ROOT implementation linearly back-extrapolates from the first
        // smoothed cluster to the beam axis. Its UKF does not propagate toward
        // the origin (which would incorrectly force vertex z to zero).
        VecX xv = x0;
        NormalizeAngles(xv);
        const double rho = std::sqrt(x0[0] * x0[0] + x0[1] * x0[1]);
        if (cfg.back_extrapolate && rho > 1.0 && x0[3] > 0.0 && cfg.back_extrap_max_mm > 0.0) {
            double sin_theta = std::sin(x0[4]);
            double path_mm = sin_theta > 0.1 ? rho / sin_theta : rho;
            path_mm = std::min(path_mm, cfg.back_extrap_max_mm);
            double dir[3];
            UnitDir(x0, dir);
            xv[0] -= dir[0] * path_mm;
            xv[1] -= dir[1] * path_mm;
            xv[2] -= dir[2] * path_mm;

            // Restore the energy lost along this short target-gas segment,
            // following AtFitterUKF's CATIMA dE/dx correction.
            const double ke_cluster = mass * (GammaOfGv(GvFromP(x0[3], mass)) - 1.0);
            catima::Projectile proj(A, Z);
            // Match AtELossCATIMA::GetdEdx: its projectile mass is stored in u,
            // so CATIMA expects KE divided by mass/u rather than mass in MeV.
            proj.T = ke_cluster / (mass / Units::u_MeV);
            const double dedx = catima::dedx(proj, dyn.GetTargetMaterial()) * dyn.GetDensity() / 10.0;
            const double ke_vertex = std::max(0.0, ke_cluster + dedx * path_mm);
            xv[3] = std::sqrt(ke_vertex * (ke_vertex + 2.0 * mass));
        }
        dyn.SetELossScale(1.0);

        out.vx = xv[0];
        out.vy = xv[1];
        out.vz = xv[2];
        out.polar = xv[4];
        out.azimuthal = xv[5];
        out.ke_vertex = mass * (GammaOfGv(GvFromP(xv[3], mass)) - 1.0);
        out.brho = (xv[3] / 1000.0) / (kBrStiffness * std::abs(Z));
        out.sigma_brho = (std::sqrt(std::max(0.0, P0(3, 3))) / 1000.0) /
                         (kBrStiffness * std::abs(Z));

        if (out.ndf >= 1 && cfg.chi2_max > 0.0 && chi2 / out.ndf > cfg.chi2_max) {
            out.failure_code = 10;
            out.ok = false;
            return false;
        }

        out.ok = std::isfinite(out.brho) && out.brho > 0.0 && std::isfinite(out.polar);
        if (!out.ok) out.failure_code = 11;
        if (!out.ok && cfg.verbose)
            std::fprintf(stderr, "[kalman] rejected: invalid final state brho=%.6g polar=%.6g\n", out.brho, out.polar);
        return out.ok;
    }
};

//----------------------------------------------------------------------
// 公开接口
//----------------------------------------------------------------------
KalmanFitter::KalmanFitter() : impl_(new Impl) {}
KalmanFitter::~KalmanFitter() = default;

bool KalmanFitter::Init(const Target& gas, int A, int Z, double mass_mev, double bfield_t,
                        double efield_vm, const KalmanConfig& cfg) {
    return impl_->Init(gas, A, Z, mass_mev, bfield_t, efield_vm, cfg);
}

bool KalmanFitter::Initialized() const { return impl_->ready; }

bool KalmanFitter::Fit(const std::vector<CloudPoint>& hits, const EstimateResult& guess,
                       KalmanFitResult& out) {
    return impl_->Fit(hits, guess, out);
}

bool KalmanFitter::PredictPath(const KalmanTrackState& prior, double path_mm,
                               const std::array<double,9>& measurement_covariance,
                               KalmanMeasurementPrediction& out) {
    return impl_->PredictPath(prior,path_mm,measurement_covariance,out);
}

bool KalmanFitter::Update(const KalmanMeasurementPrediction& prediction,
                          const std::array<double,3>& measurement,
                          KalmanMeasurementUpdate& out) const {
    return impl_->Update(prediction,measurement,out);
}

const KalmanConfig& KalmanFitter::Config() const { return impl_->cfg; }

} // namespace rk4
