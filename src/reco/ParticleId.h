#ifndef RK4_RECO_PARTICLEID_H
#define RK4_RECO_PARTICLEID_H

//----------------------------------------------------------------------
// 粒子鉴别（dE/dx vs Bρ）
//
// 轨迹形状只约束 Bρ，不约束粒子种类（质子与 α 在同一 Bρ 上轨迹相同）；
// 必须用比能损区分：dE/dx ∝ Z²/β² · [ln(2 m_e c² β²γ² / I) - β²]。
//
// 本实现用简化的 Bethe 公式做"候选之间的相对比较"：
// 对每个候选 (A, Z)，由 Bρ 得 p = 299.8·Bρ·Z，再由质量得 βγ，
// 计算期望 dE/dx，与实测比较，取对数偏差最小者。
//
// 注意：dE/dx 的绝对标度取决于探测器增益，用 scale 参数标定；
// 判别本身对 scale 不敏感（候选间相对比较）。
//----------------------------------------------------------------------

#include <string>
#include <vector>

namespace rk4 {

struct PidCandidate {
    std::string name;
    int A = 1;
    int Z = 1;
};

struct PidResult {
    bool ok = false;
    int index = -1;
    std::string name;
    double score = 0.0;        // 最佳候选的对数偏差
    double second_score = 0.0; // 次佳，用于判断置信度
    double separation = 0.0;   // second_score - score，越大越可信
};

class ParticleId {
public:
    ParticleId();

    // 默认候选：p, d, t, 3He, alpha
    void SetCandidates(const std::vector<PidCandidate>& candidates);
    const std::vector<PidCandidate>& Candidates() const { return candidates_; }

    // dE/dx 标定常数（把公式值映射到点云的 integral/mm）
    void SetScale(double scale) { scale_ = scale; }
    double Scale() const { return scale_; }

    // 用一条已知粒子的 (Bρ, dEdx) 标定 scale
    void Calibrate(double brho, double dedx, int A, int Z);

    bool Identify(double brho, double dedx, PidResult& out) const;

    // 理论相对 dE/dx（未乘 scale）
    static double ExpectedDedx(double brho, int A, int Z);

private:
    std::vector<PidCandidate> candidates_;
    double scale_ = 1.0;
};

} // namespace rk4

#endif // RK4_RECO_PARTICLEID_H
