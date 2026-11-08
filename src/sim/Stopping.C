#include "Stopping.h"

#include <cstdio>
#include <algorithm>

bool StoppingPower::Init(int A, int Z, const Target& target) {
    valid_ = false;
    mass_u_ = MassData::GetNuclearMassNuclear(A, Z);
    if (mass_u_ <= 0.0) {
        std::fprintf(stderr, "[StoppingPower] 核质量读取失败 A=%d Z=%d\n", A, Z);
        return false;
    }

    const std::vector<Element> comp = target.compound();
    if (comp.empty()) {
        std::fprintf(stderr, "[StoppingPower] 靶材组成为空\n");
        return false;
    }

    material_ = catima::Material();
    const double scale = Utils::stoichiometricScale(comp);
    for (const auto& e : comp) {
        // Avoid rounding the smallest count just below Catima's boundary 1.
        material_.add_element(e.a, e.z, std::max(1.0, e.s * scale));
    }

    density_ = target.density();
    if (density_ <= 0.0) {
        std::fprintf(stderr, "[StoppingPower] 靶材密度非正：%g\n", density_);
        return false;
    }
    material_.density(density_);
    material_.thickness(0.0);

    projectile_ = catima::Projectile(A, Z);
    valid_ = true;
    return true;
}

double StoppingPower::EnergyAfter(double ke_mev, double depth_m) const {
    if (!valid_) return ke_mev;
    if (depth_m <= 0.0) return ke_mev;

    // m -> cm，再乘密度得到 g/cm^2
    const double thickness = depth_m * 100.0 * density_;
    material_.thickness(thickness);

    projectile_.T = ke_mev / mass_u_; // MeV/u
    const double out = catima::energy_out(projectile_, material_) * mass_u_;
    return out > 0.0 ? out : 0.0;
}

double StoppingPower::Dedx(double ke_mev) const {
    if (!valid_) return 0.0;
    projectile_.T = ke_mev / mass_u_;
    return catima::dedx(projectile_, material_);
}
