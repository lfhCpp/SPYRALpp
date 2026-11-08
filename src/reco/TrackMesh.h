#ifndef RK4_RECO_TRACKMESH_H
#define RK4_RECO_TRACKMESH_H

//----------------------------------------------------------------------
// 轨迹网格（对应 Spyral 的 TrackMesh / 插值表）
//
// 在 (ke, polar) 二维网格上预先解好 ODE，每个格点存一条定长时间步的轨迹。
// 求解时对 (ke, polar) 做双线性插值，再旋转/平移到实际顶点，
// 避免在拟合循环里反复积分。
//
// 网格参考顶点固定在 (0, 0, z_ref)，z_ref 默认取探测器长度的一半，
// 这样向前/向后的径迹都有足够的活动空间。
//----------------------------------------------------------------------

#include "../sim/TargetDensity.h"

#include <string>
#include <vector>
#include <cstring>
#include <istream>

namespace rk4 {

struct MeshConfig {
    int ke_bins = 256;        // 动能网格；高分辨率生产默认
    int polar_bins = 129;     // 极角网格；高分辨率生产默认
    double ke_min = 0.1;      // MeV；与已验证的生产 mesh 配置一致
    double ke_max = 70.0;
    double polar_min = 2.0;   // 度
    double polar_max = 88.0;
    int n_time = 400;         // 每条轨迹的时间步数；与 0.3 ns 默认步长对应约 120 ns
    double dt = 3.0e-10;      // 时间步长 (s)，默认 0.3 ns
    double z_ref = 0.5;       // 参考顶点 z (m)
    bool log_ke = true;       // Uniform relative energy resolution, linear interpolation in KE.
};

class TrackMesh {
public:
    ~TrackMesh();
    TrackMesh() = default;
    TrackMesh(const TrackMesh&) = delete;
    TrackMesh& operator=(const TrackMesh&) = delete;

    // 用给定粒子与靶生成网格。磁场、电场参数传正的配置幅值；传播时
    // 内部使用 Spyral 1.0 约定的负 B/E 符号。
    bool Build(const MeshConfig& cfg, int A, int Z, const Target& gas,
               double bfield_t, double efield_vm);

    bool Save(const std::string& path) const;
    bool Load(const std::string& path);

    // Memory-map a Python Spyral NumPy trajectory mesh (time, xyz, KE, polar).
    bool LoadPythonNpy(const std::string& path, const MeshConfig& cfg, int A, int Z,
                       double bfield_t, double efield_vm);
    bool IsPythonMesh() const { return python_mesh_; }
    bool InterpolatePython(double ke_mev, double polar_rad, std::vector<double>& xyz) const;

    // 用 mmap(MAP_SHARED) 映射网格文件：多个进程共享同一份物理内存，
    // 不会各自拷贝几百 KB 的数据。映射区只读，生命周期由本对象管理。
    bool LoadShared(const std::string& path);

    // 是否来自共享映射
    bool IsShared() const { return shared_view_ != nullptr; }

    // 双线性插值，返回 n_time*3 个位置（mm），坐标系：参考顶点在 (0,0,z_ref)
    bool Interpolate(double ke_mev, double polar_rad, std::vector<double>& xyz) const;

    bool IsInRange(double ke_mev, double polar_rad) const;

    // Three-argument form checks grid/particle only; use the full overload to
    // decide whether a physics cache can be reused for a gas and field setting.
    bool MatchesConfig(const MeshConfig& cfg, int A, int Z) const;
    bool MatchesConfig(const MeshConfig& cfg, int A, int Z, const Target& gas,
                       double bfield_t, double efield_vm) const;

    const MeshConfig& Config() const { return cfg_; }
    int ParticleA() const { return particle_A_; }
    int ParticleZ() const { return particle_Z_; }
    int TimeSteps() const { return cfg_.n_time; }
    bool Empty() const { return Count() == 0; }

private:
    size_t Index(int t, int comp, int ike, int ipolar) const;
    void ResetStorage();
    bool ReadHeader(std::istream& input);
    double EnergyAt(int index) const;

    MeshConfig cfg_;
    int particle_A_ = 0;
    int particle_Z_ = 0;
    std::vector<double> data_; // [t][3][ke][polar]（自有数据）

    // 共享映射（LoadShared 时有效）
    const char* shared_view_ = nullptr;
    void* map_base_ = nullptr;
    size_t map_length_ = 0;
    size_t shared_count_ = 0;
    bool has_context_ = false;
    bool python_mesh_ = false;
    double bfield_ = 0, efield_ = 0, density_ = 0, particle_mass_ = 0;
    std::vector<Element> composition_;

    // 统一读取入口：优先共享视图
    double At(size_t index) const {
        if (!shared_view_) return data_[index];
        double value;
        std::memcpy(&value,shared_view_ + index*sizeof(double),sizeof(double));
        return value;
    }
    size_t Count() const {
        return shared_view_ != nullptr ? shared_count_ : data_.size();
    }
};

} // namespace rk4

#endif // RK4_RECO_TRACKMESH_H
