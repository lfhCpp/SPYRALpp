#ifndef RK4_CORE_EVENT_H
#define RK4_CORE_EVENT_H

//----------------------------------------------------------------------
// 事件与粒子
//
// 约定（与 attpc_engine 对齐）：
//   反应 a(b,c)d：0=靶 a，1=束流 b，2=出射 c，3=剩余核 d
//   后续衰变依次占 4,5 / 6,7 / ...
//
// 单位：动量 MeV/c，能量 MeV（含静止质量），位置 m
//----------------------------------------------------------------------

#include <cmath>
#include <vector>

namespace rk4 {

struct Particle {
    int A = 0;      // 质量数
    int Z = 0;      // 电荷数
    int label = -1; // 在反应链中的序号

    double px = 0.0, py = 0.0, pz = 0.0; // MeV/c
    double E = 0.0;                      // MeV

    double P2() const { return px * px + py * py + pz * pz; }
    double P() const { return std::sqrt(P2()); }

    // 不变质量；数值误差导致 m^2 < 0 时返回 0
    double Mass() const {
        const double m2 = E * E - P2();
        return m2 > 0.0 ? std::sqrt(m2) : 0.0;
    }

    double Kinetic() const { return E - Mass(); }
};

struct Event {
    long long id = 0;
    double vx = 0.0, vy = 0.0, vz = 0.0; // 顶点位置 (m)
    double beam_energy = 0.0;            // 顶点处束流动能 (MeV)
    std::vector<Particle> particles;
};

// 四动量求和，用于守恒性检查
struct FourMomentum {
    double px = 0.0, py = 0.0, pz = 0.0, E = 0.0;
    void Add(const Particle& p) {
        px += p.px;
        py += p.py;
        pz += p.pz;
        E += p.E;
    }
};

} // namespace rk4

#endif // RK4_CORE_EVENT_H
