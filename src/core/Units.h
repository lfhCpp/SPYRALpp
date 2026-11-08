#ifndef RK4_CORE_UNITS_H
#define RK4_CORE_UNITS_H

//----------------------------------------------------------------------
// 全局单位与物理常数（唯一来源）
//
// 约定：
//   能量      MeV
//   质量      u（原子质量单位），需要 MeV/c^2 时乘 Units::u_MeV
//   长度      m（对外 I/O、pad 几何用 mm，只在边界处换算）
//   时间      s
//   角度      rad
//   磁场      T
//   电场      V/m
//   密度      g/cm^3，面密度 g/cm^2
//
// 说明：本文件不依赖 ROOT，任何模块都可以单独包含。
//----------------------------------------------------------------------

namespace Units {

// ---- 基本换算 ----
constexpr double u_MeV     = 931.49410242;      // 1 u -> MeV/c^2
constexpr double C         = 2.99792458e8;      // 光速 (m/s)
constexpr double E_CHARGE  = 1.602176634e-19;   // 基本电荷 (C)
constexpr double MEV_2_KG  = 1.782661921e-30;   // MeV/c^2 -> kg
constexpr double MEV_2_JOULE = 1.602176634e-13; // MeV -> J
constexpr double me_u      = 5.48579909065e-4;  // 电子静止质量 (u)

// ---- 角度 ----
constexpr double PI       = 3.14159265358979323846;
constexpr double DEG2RAD  = PI / 180.0;
constexpr double RAD2DEG  = 180.0 / PI;

// ---- 气体 ----
constexpr double gasConstant = 62363.7;         // cm^3·Torr/(mol·K)
constexpr double roomTemp    = 293.15;          // K (20 °C)
constexpr double avogadro    = 6.02214076e23;   // 1/mol
constexpr double boltzmann   = 1.380649e-23;    // J/K
constexpr double ug2g        = 1.0e-6;          // μg -> g

// ---- 探测器默认值（后续 DetectorConfig 会覆盖）----
constexpr double at_tpc_length      = 1.0;      // m，灵敏体积长度
constexpr double at_tpc_rho_max     = 0.28;     // m，灵敏体积半径
constexpr double at_tpc_bfield      = 2.85;     // T
constexpr double at_tpc_efield      = 45000.0;  // V/m

} // namespace Units

#endif // RK4_CORE_UNITS_H
