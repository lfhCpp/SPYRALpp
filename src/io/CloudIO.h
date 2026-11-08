#ifndef RK4_IO_CLOUDIO_H
#define RK4_IO_CLOUDIO_H

//----------------------------------------------------------------------
// 点云 I/O：ROOT TTree
//
// 每个事件一条记录：
//   event / npoints
//   x[npoints] y[npoints] z[npoints]      (mm)
//   amplitude[npoints] integral[npoints]
//   pad[npoints] tb[npoints] scale[npoints]
//----------------------------------------------------------------------

#include "../sim/Simulator.h"

#include <string>
#include <vector>

namespace rk4 {

bool WritePointClouds(const std::string& path, const std::vector<PointCloud>& clouds,
                      const std::string& tree_name = "cloud");

bool ReadPointClouds(const std::string& path, std::vector<PointCloud>& clouds,
                     const std::string& tree_name = "cloud");

} // namespace rk4

#endif // RK4_IO_CLOUDIO_H
