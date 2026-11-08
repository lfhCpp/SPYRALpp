#include "Reaction.h"
#include "MassData.h"

#include <cmath>
#include <cstdio>
#include <sstream>

namespace rk4 {

namespace {

// Källén 函数 lambda(x, y, z) = x^2 + y^2 + z^2 - 2xy - 2xz - 2yz
inline double Kallen(double x, double y, double z) {
    return (x - y - z) * (x - y - z) - 4.0 * y * z;
}

inline double MomentumCM(double s, double m1, double m2) {
    const double lambda = Kallen(s, m1 * m1, m2 * m2);
    if (lambda <= 0.0) return 0.0;
    return std::sqrt(lambda) / (2.0 * std::sqrt(s));
}

} // namespace

//----------------------------------------------------------------------
// TwoBodyReaction
//----------------------------------------------------------------------
bool TwoBodyReaction::Init(const ReactionConfig& cfg) {
    cfg_ = cfg;
    m_beam_ = MassData::GetNuclearMassNuclear(cfg.beam_A, cfg.beam_Z) * Units::u_MeV;
    m_target_ = MassData::GetNuclearMassNuclear(cfg.target_A, cfg.target_Z) * Units::u_MeV;
    m_ejectile_ = MassData::GetNuclearMassNuclear(cfg.ejectile_A, cfg.ejectile_Z) * Units::u_MeV;

    if (m_beam_ <= 0.0 || m_target_ <= 0.0 || m_ejectile_ <= 0.0) {
        std::fprintf(stderr, "[TwoBodyReaction] 质量读取失败\n");
        return false;
    }

    residual_z_ = cfg.target_Z + cfg.beam_Z - cfg.ejectile_Z;
    residual_a_ = cfg.target_A + cfg.beam_A - cfg.ejectile_A;
    if (residual_z_ < 0 || residual_a_ < residual_z_) {
        std::fprintf(stderr, "[TwoBodyReaction] 非法反应：剩余核 Z=%d A=%d\n",
                     residual_z_, residual_a_);
        return false;
    }

    m_residual_ = MassData::GetNuclearMassNuclear(residual_a_, residual_z_) * Units::u_MeV;
    if (m_residual_ <= 0.0) {
        std::fprintf(stderr, "[TwoBodyReaction] 剩余核质量读取失败\n");
        return false;
    }

    std::ostringstream oss;
    oss << cfg.target_A << MassData::GetNuclearName(cfg.target_A, cfg.target_Z)
        << "(" << cfg.beam_A << MassData::GetNuclearName(cfg.beam_A, cfg.beam_Z)
        << "," << cfg.ejectile_A << MassData::GetNuclearName(cfg.ejectile_A, cfg.ejectile_Z)
        << ")" << residual_a_ << MassData::GetNuclearName(residual_a_, residual_z_);
    symbol_ = oss.str();
    return true;
}

bool TwoBodyReaction::IsAllowed(double beam_ke_mev, double excitation_mev) const {
    const double s = m_beam_ * m_beam_ + m_target_ * m_target_ +
                     2.0 * m_target_ * (m_beam_ + beam_ke_mev);
    const double m_out = m_ejectile_ + m_residual_ + excitation_mev;
    return std::sqrt(s) > m_out;
}

double TwoBodyReaction::ThresholdEnergy(double excitation_mev) const {
    const double m_out = m_ejectile_ + m_residual_ + excitation_mev;
    const double s_min = m_out * m_out;
    // s = mb^2 + mt^2 + 2 mt (mb + T)
    return (s_min - m_beam_ * m_beam_ - m_target_ * m_target_) / (2.0 * m_target_) - m_beam_;
}

bool TwoBodyReaction::Calculate(double beam_ke_mev, double theta_cm, double phi_cm,
                                double excitation_mev, Event& ev) const {
    if (!IsAllowed(beam_ke_mev, excitation_mev)) return false;

    // 初态
    const double p_beam = std::sqrt(beam_ke_mev * (beam_ke_mev + 2.0 * m_beam_));
    const double e_tot = beam_ke_mev + m_beam_ + m_target_;
    const double s = e_tot * e_tot - p_beam * p_beam;
    const double sqrt_s = std::sqrt(s);

    // 质心系出射粒子动量
    const double m_res_ex = m_residual_ + excitation_mev;
    const double p_cm = MomentumCM(s, m_ejectile_, m_res_ex);
    const double e_ej_cm = std::sqrt(p_cm * p_cm + m_ejectile_ * m_ejectile_);

    // 质心系 -> 实验室系
    const double beta = p_beam / e_tot;
    const double gamma = e_tot / sqrt_s;

    const double pz_cm = p_cm * std::cos(theta_cm);
    const double px_cm = p_cm * std::sin(theta_cm) * std::cos(phi_cm);
    const double py_cm = p_cm * std::sin(theta_cm) * std::sin(phi_cm);

    Particle target;
    target.A = cfg_.target_A;
    target.Z = cfg_.target_Z;
    target.label = 0;
    target.px = target.py = target.pz = 0.0;
    target.E = m_target_;

    Particle beam;
    beam.A = cfg_.beam_A;
    beam.Z = cfg_.beam_Z;
    beam.label = 1;
    beam.px = 0.0;
    beam.py = 0.0;
    beam.pz = p_beam;
    beam.E = beam_ke_mev + m_beam_;

    Particle ejectile;
    ejectile.A = cfg_.ejectile_A;
    ejectile.Z = cfg_.ejectile_Z;
    ejectile.label = 2;
    ejectile.px = px_cm;
    ejectile.py = py_cm;
    ejectile.pz = gamma * (pz_cm + beta * e_ej_cm);
    ejectile.E = gamma * (e_ej_cm + beta * pz_cm);

    Particle residual;
    residual.A = residual_a_;
    residual.Z = residual_z_;
    residual.label = 3;
    residual.px = -ejectile.px;
    residual.py = -ejectile.py;
    residual.pz = p_beam - ejectile.pz;
    residual.E = e_tot - ejectile.E;

    ev.particles.clear();
    ev.particles.push_back(target);
    ev.particles.push_back(beam);
    ev.particles.push_back(ejectile);
    ev.particles.push_back(residual);
    ev.beam_energy = beam_ke_mev;
    return true;
}

//----------------------------------------------------------------------
// TwoBodyDecay
//----------------------------------------------------------------------
bool TwoBodyDecay::Init(int parent_A, int parent_Z, int frag1_A, int frag1_Z,
                        int label1, int label2) {
    m_parent_ = MassData::GetNuclearMassNuclear(parent_A, parent_Z) * Units::u_MeV;
    m_frag1_ = MassData::GetNuclearMassNuclear(frag1_A, frag1_Z) * Units::u_MeV;

    frag1_a_ = frag1_A;
    frag1_z_ = frag1_Z;
    frag2_a_ = parent_A - frag1_A;
    frag2_z_ = parent_Z - frag1_Z;
    label1_ = label1;
    label2_ = label2;

    if (frag2_z_ < 0 || frag2_a_ < frag2_z_ || m_parent_ <= 0.0 || m_frag1_ <= 0.0) {
        std::fprintf(stderr, "[TwoBodyDecay] 非法衰变：父核 A=%d Z=%d -> A=%d Z=%d\n",
                     parent_A, parent_Z, frag1_A, frag1_Z);
        return false;
    }
    m_frag2_ = MassData::GetNuclearMassNuclear(frag2_a_, frag2_z_) * Units::u_MeV;
    return m_frag2_ > 0.0;
}

bool TwoBodyDecay::Calculate(const Particle& parent, double theta, double phi,
                             double frag2_excitation_mev,
                             Particle& frag1, Particle& frag2) const {
    const double mass_parent = parent.Mass();
    const double m_frag2_ex = m_frag2_ + frag2_excitation_mev;

    if (mass_parent <= m_frag1_ + m_frag2_ex) {
        return false; // 能量不够
    }

    const double p_cm = MomentumCM(mass_parent * mass_parent, m_frag1_, m_frag2_ex);
    const double e1_cm = std::sqrt(p_cm * p_cm + m_frag1_ * m_frag1_);

    const double p1x = p_cm * std::sin(theta) * std::cos(phi);
    const double p1y = p_cm * std::sin(theta) * std::sin(phi);
    const double p1z = p_cm * std::cos(theta);

    // 父核的实验室系速度
    const double p_par = parent.P();
    const double beta = (mass_parent > 0.0) ? p_par / parent.E : 0.0;
    const double gamma = (mass_parent > 0.0) ? parent.E / mass_parent : 1.0;

    if (p_par > 0.0) {
        const double ux = parent.px / p_par;
        const double uy = parent.py / p_par;
        const double uz = parent.pz / p_par;

        const double pp_par = p1x * ux + p1y * uy + p1z * uz; // 平行分量
        const double px_perp = p1x - pp_par * ux;
        const double py_perp = p1y - pp_par * uy;
        const double pz_perp = p1z - pp_par * uz;

        const double pp_par_lab = gamma * (pp_par + beta * e1_cm);
        frag1.px = px_perp + pp_par_lab * ux;
        frag1.py = py_perp + pp_par_lab * uy;
        frag1.pz = pz_perp + pp_par_lab * uz;
        frag1.E = gamma * (e1_cm + beta * pp_par);
    } else {
        // 父核静止：质心系即实验室系
        frag1.px = p1x;
        frag1.py = p1y;
        frag1.pz = p1z;
        frag1.E = e1_cm;
    }

    frag1.A = frag1_a_;
    frag1.Z = frag1_z_;
    frag1.label = label1_;

    frag2.A = frag2_a_;
    frag2.Z = frag2_z_;
    frag2.label = label2_;
    frag2.px = parent.px - frag1.px;
    frag2.py = parent.py - frag1.py;
    frag2.pz = parent.pz - frag1.pz;
    frag2.E = parent.E - frag1.E;

    return true;
}

} // namespace rk4
