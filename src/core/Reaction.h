#ifndef RK4_CORE_REACTION_H
#define RK4_CORE_REACTION_H

//----------------------------------------------------------------------
// 二体反应与二体衰变（相对论运动学）
//
// 反应书写：a(b,c)d
//   a = target（静止）
//   b = beam（弹核）
//   c = ejectile（出射粒子，抽样其质心系角度）
//   d = residual（剩余核，可带激发能）
//
// 全部使用原子核质量（u -> MeV），输入输出单位均为 MeV、rad。
//----------------------------------------------------------------------

#include "Event.h"

#include <string>

namespace rk4 {

struct ReactionConfig {
    int beam_A = 0, beam_Z = 0;
    int target_A = 0, target_Z = 0;
    int ejectile_A = 0, ejectile_Z = 0;
    double excitation_mev = 0.0; // 剩余核激发能
};

class TwoBodyReaction {
public:
    // 读取质量、推导剩余核；失败返回 false
    bool Init(const ReactionConfig& cfg);

    // 给定束流动能（MeV）与激发能，判断是否在能量上允许
    bool IsAllowed(double beam_ke_mev, double excitation_mev) const;

    // 阈能（MeV，对应给定激发能）
    double ThresholdEnergy(double excitation_mev) const;

    // 计算一条事件；theta/phi 为质心系出射粒子角度（rad）
    bool Calculate(double beam_ke_mev, double theta_cm, double phi_cm,
                   double excitation_mev, Event& ev) const;

    int ResidualA() const { return residual_a_; }
    int ResidualZ() const { return residual_z_; }
    double BeamMass() const { return m_beam_; }
    double TargetMass() const { return m_target_; }
    double EjectileMass() const { return m_ejectile_; }
    double ResidualMass() const { return m_residual_; }

    // 反应式字符串，例如 "1H(12Be,1H)12Be"
    const char* Symbol() const { return symbol_.c_str(); }

private:
    ReactionConfig cfg_;
    int residual_a_ = 0;
    int residual_z_ = 0;
    double m_beam_ = 0.0;     // MeV/c^2
    double m_target_ = 0.0;
    double m_ejectile_ = 0.0;
    double m_residual_ = 0.0;
    std::string symbol_;
};

//----------------------------------------------------------------------
// 二体衰变：parent -> frag1 + frag2
// 父核的四动量由上游给出（已包含激发能），此处只需再给定 frag2 的激发能
//----------------------------------------------------------------------
class TwoBodyDecay {
public:
    bool Init(int parent_A, int parent_Z, int frag1_A, int frag1_Z,
              int label1, int label2);

    // parent 为实验室系四动量；theta/phi 为 frag1 在父核静止系的角度
    bool Calculate(const Particle& parent, double theta, double phi,
                   double frag2_excitation_mev, Particle& frag1, Particle& frag2) const;

    int Frag1A() const { return frag1_a_; }
    int Frag1Z() const { return frag1_z_; }
    int Frag2A() const { return frag2_a_; }
    int Frag2Z() const { return frag2_z_; }
    double ParentMass() const { return m_parent_; }
    double Frag1Mass() const { return m_frag1_; }
    double Frag2Mass() const { return m_frag2_; }

private:
    int frag1_a_ = 0, frag1_z_ = 0;
    int frag2_a_ = 0, frag2_z_ = 0;
    int label1_ = -1, label2_ = -1;
    double m_parent_ = 0.0;
    double m_frag1_ = 0.0;
    double m_frag2_ = 0.0;
};

} // namespace rk4

#endif // RK4_CORE_REACTION_H
