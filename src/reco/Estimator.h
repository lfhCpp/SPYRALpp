#ifndef RK4_RECO_ESTIMATOR_H
#define RK4_RECO_ESTIMATOR_H

//----------------------------------------------------------------------
// 径迹参数粗估（对应 Spyral 的 EstimationPhase）
//
// 输出：顶点、极角/方位角、磁刚度 Bρ、比能损 dE/dx、弧长。
// 这里只做"够给求解器当初值"的估计，精度约 10%。
//----------------------------------------------------------------------

#include "Cluster.h"

namespace rk4 {

struct EstimateResult {
    bool ok = false;
    double vertex_x = 0.0, vertex_y = 0.0, vertex_z = 0.0; // mm
    double center_x = 0.0, center_y = 0.0;                 // mm
    double polar = 0.0;                                    // rad
    double azimuthal = 0.0;                                // rad
    double brho = 0.0;                                     // T*m
    double radius = 0.0;                                   // mm
    double dEdx = 0.0;                                     // 电荷/mm
    double dedx_pred = 0.0;                                // 沿路径的 Bethe 平均（MeV/mm，质子假设）
    double arclength = 0.0;                                // mm
    double charge = 0.0;                                   // 总积分电荷
    int direction = -1;
};

class Estimator {
public:
    // bfield 单位 T
    bool Estimate(const Cluster& cluster, double bfield_t, EstimateResult& out) const;

    // 拟合前是否对簇内 x(z)/y(z) 做平滑样条降噪（默认关）
    // Spyral 的 dE/dx 只积到"小 pad 区"边界（默认 152 mm），
    // 因为大小 pad 的电荷刻度不同，混在一起会偏。
    void SetSmallPadRadius(double mm) { small_pad_radius_mm_ = mm; }


    // pad 间距（mm）：决定"一个命中点覆盖多长路径"
    void SetPadPitch(double mm) { pad_pitch_mm_ = mm; }
    // 顶点必须落在束流轴附近（Spyral 的 beam_region_radius）
    void SetBeamRegionRadius(double mm) { beam_region_radius_mm_ = mm; }
    // 靶气体（用于沿路径能量损失预测）
    // 沿路径降能强度：0 = 不降能（退化成顶点 Bethe），1 = catima 全量
    void SetPredStrength(double f) { pred_strength_ = f; }
    void SetTarget(int A, int Z, double stoich, double density_gcm3) {
        target_A_ = A;
        target_Z_ = Z;
        target_stoich_ = stoich;
        target_density_ = density_gcm3;
        target_ready_ = (A > 0 && Z > 0 && density_gcm3 > 0.0);
    }
    // 电荷收集半径（mm）：把簇外但邻近的 pad 电荷也收进来；<=0 关闭（退回旧的整弧口径）
    void SetChargeRadius(double mm) { charge_radius_mm_ = mm; }

    // 新接口：能用整条想点云做电荷收集；旧接口内部转调本函数（cloud = nullptr）
    bool Estimate(const Cluster& cluster, const PointCloud& cloud, double bfield_t,
                  EstimateResult& out) const;
    void SetMinTrajectoryPoints(int n) { min_trajectory_points_=n; }
    double SmallPadRadius() const { return small_pad_radius_mm_; }

    void SetSmoothing(bool enabled, double lambda = 1.0) {
        smooth_ = enabled;
        lambda_ = lambda;
    }
    bool Smoothing() const { return smooth_; }
    // Production default follows the Spyral 1.0 estimator; legacy remains explicit.
    void SetUseReinsch(bool enabled, double lambda = 100.0) {
        use_reinsch_ = enabled; reinsch_lambda_ = lambda;
    }

private:
    bool EstimateReinsch(const Cluster& cluster, double bfield_t, EstimateResult& out) const;
    bool use_reinsch_ = true;
    int min_trajectory_points_ = 30;
    double reinsch_lambda_ = 100.0;
    // 实测：真实波形数据上 lambda=1 使拟合残差从 13.08 mm 降到 11.23 mm（lambda=10 会过度平滑反弹到 13.04）
    // 实测：P-spline 预平滑 + 局部二次叠加会过度平滑，长度信息被抹掉
    // （r(log dE/dx, log L) 从 -0.01 回到 -0.62）。默认只保留局部二次。
    bool smooth_ = false;
    double lambda_ = 3000.0;
    double small_pad_radius_mm_ = 152.0;
    double pad_pitch_mm_ = 10.0;
    double beam_region_radius_mm_ = 25.0;
    int target_A_ = 0, target_Z_ = 0;
    double target_stoich_ = 2;
    double target_density_ = 0.0;
    bool target_ready_ = false;
    double pred_strength_ = 0.3;
    double charge_radius_mm_ = 0.0;
};

} // namespace rk4

#endif // RK4_RECO_ESTIMATOR_H
