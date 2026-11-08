#ifndef RK4_IO_RESULTIO_H
#define RK4_IO_RESULTIO_H

#include "../reco/PidState.h"

//----------------------------------------------------------------------
// 拟合结果 I/O（ROOT TTree，扁平 ntuple）
//----------------------------------------------------------------------

#include <string>
#include <vector>

namespace rk4 {

struct FitRecord {
    long long event = -1;
    int cluster = -1;
    double brho = 0.0;
    double sigma_brho = 0.0;
    double polar = 0.0;
    double azimuthal = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double cost = 0.0;
    double chi2ndf = 0.0; // Kalman 后端的 χ²/ndf；Minuit 路径为 0
    int ok = 0;
    int pid_index = -1;    // 候选粒子索引（-1 表示未判别）
    double pid_sep = 0.0;  // 与次佳候选的对数偏差差，越大越可信
    double dedx = 0.0;     // 实测比能损（用于离线检查）
    // Additive state metadata. Older ROOT trees have no such branches; readers
    // retain NOT_RECORDED defaults for those files.
    int fit_state = static_cast<int>(FitState::NOT_RECORDED);
    int pid_state = static_cast<int>(PidState::NOT_RECORDED);
    int pid_unknown_reason = static_cast<int>(PidUnknownReason::NOT_RECORDED);
    int ic_observation_state = static_cast<int>(IcObservationState::NOT_RECORDED);
    int pid_state_basis = static_cast<int>(PidStateBasis::NOT_RECORDED);
    double ic_amplitude = -1.0;
    int migrad_status = -1;
    double migrad_edm = 0.0;
};

bool WriteFitResults(const std::string& path, const std::vector<FitRecord>& records,
                     const std::string& tree_name = "fit");

bool ReadFitResults(const std::string& path, std::vector<FitRecord>& records,
                    const std::string& tree_name = "fit");

} // namespace rk4

#endif // RK4_IO_RESULTIO_H
