#ifndef RK4_IO_H5WRITER_H
#define RK4_IO_H5WRITER_H

//----------------------------------------------------------------------
// HDF5 输出（可选，依赖 HighFive + libhdf5）
//
// 布局与 attpc_engine 的 KinematicsPipeline 输出保持一致：
//
//   /data                     attrs: n_events, proton_numbers, mass_numbers,
//                                    chunk_size, n_chunks
//     /chunk_0                attrs: min_event, max_event
//       /event_0              dataset (N,4) = [px, py, pz, E]
//                             attrs: vertex_x, vertex_y, vertex_z
//
// 这样生成的文件可以直接被现有 Python 工具链读取。
// 未启用 HDF5 时，Hdf5Available() 返回 false，写出函数返回 false。
//----------------------------------------------------------------------

#include "../core/Event.h"
#include "../sim/Simulator.h"

#include <string>
#include <vector>

namespace rk4 {

bool Hdf5Available();

bool WriteKinematicsHdf5(const std::string& path, const std::vector<Event>& events,
                         int chunk_size = 1000000);

// 点云 HDF5 输出，布局与 attpc_engine 的 SpyralWriter 一致：
//
//   /cloud                  attrs: min_event, max_event
//     /cloud_{N}            dataset (M,8) = [x, y, z, amplitude, integral, pad, tb, scale]
//                           attrs: orig_run, orig_event, ic_amplitude/integral/centroid/multiplicity = -1
//     /labels_{N}           dataset (M,) = 产生该点的粒子索引
//
// 这样 Python 版 Spyral 可以直接读。
// labels 可为空（此时写全 -1）。
bool WritePointCloudsHdf5(const std::string& path, const std::vector<PointCloud>& clouds,
                          const std::vector<std::vector<long long>>& labels = {});

} // namespace rk4

#endif // RK4_IO_H5WRITER_H
