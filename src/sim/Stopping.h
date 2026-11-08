#ifndef RK4_SIM_STOPPING_H
#define RK4_SIM_STOPPING_H

//----------------------------------------------------------------------
// 阻止本领（dE/dx）封装，内部用 catima
//
// 只暴露两个量：
//   EnergyAfter(ke, depth)  —— 穿过 depth 米靶材后剩余动能 (MeV)
//   Dedx(ke)                —— 质量阻止本领 MeV/(g/cm^2)
//
// 这样以后要换成 SRIM 表或别的实现，只需替换本文件。
//----------------------------------------------------------------------

#include "../core/MassData.h"
#include "../core/Units.h"
#include "TargetDensity.h"

#include "catima/catima.h"

class StoppingPower {
public:
    // 用靶材的化合物组成与密度初始化；失败返回 false
    bool Init(int A, int Z, const Target& target);

    // depth 单位 m，返回剩余动能 MeV（不会为负）
    double EnergyAfter(double ke_mev, double depth_m) const;

    // 质量阻止本领 MeV/(g/cm^2)
    double Dedx(double ke_mev) const;

    bool IsValid() const { return valid_; }

private:
    mutable catima::Material material_;
    mutable catima::Projectile projectile_;
    double density_ = 0.0;   // g/cm^3
    double mass_u_ = 0.0;    // 核质量 (u)
    bool valid_ = false;
};

#endif // RK4_SIM_STOPPING_H
