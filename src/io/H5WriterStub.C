#include "H5Writer.h"

// 未启用 RK4_WITH_HDF5 时的占位实现：
// 保持接口可用，调用方不必写条件编译。

namespace rk4 {

bool Hdf5Available() { return false; }

bool WriteKinematicsHdf5(const std::string& path, const std::vector<Event>& events,
                         int chunk_size) {
    (void)path;
    (void)events;
    (void)chunk_size;
    return false;
}

} // namespace rk4
