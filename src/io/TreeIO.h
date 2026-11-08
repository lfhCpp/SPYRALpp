#ifndef RK4_IO_TREEIO_H
#define RK4_IO_TREEIO_H

//----------------------------------------------------------------------
// 事件 I/O：ROOT TTree
//
// 内部标准格式。每个事件写成一条记录：
//   event / vx / vy / vz / beam_energy / npart
//   A[npart] Z[npart] label[npart]
//   px[npart] py[npart] pz[npart] E[npart]
//
// 头文件不直接包含 ROOT，避免把 ROOT 拖进不需要它的编译单元。
//----------------------------------------------------------------------

#include "../core/Event.h"

#include <string>
#include <vector>

namespace rk4 {

constexpr int kMaxParticles = 16;

bool WriteKinematics(const std::string& path, const std::vector<Event>& events,
                     const std::string& tree_name = "kinematics");

bool ReadKinematics(const std::string& path, std::vector<Event>& events,
                    const std::string& tree_name = "kinematics");

} // namespace rk4

#endif // RK4_IO_TREEIO_H
