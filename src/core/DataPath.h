#ifndef RK4_CORE_DATAPATH_H
#define RK4_CORE_DATAPATH_H

//----------------------------------------------------------------------
// 数据目录解析
//
// 程序运行时需要 data/ 下的核素质量表、SRIM 停止本领、pad 几何等文件。
// 为了同时支持“在源码树里跑”和“安装后跑”，统一按下面的顺序查找：
//
//   1. 环境变量 SPYRALPP_DATA_DIR
//   2. 当前目录依次尝试：data/、./、../data/、../../data/
//      （在 SpyralPP/ 或 SpyralPP/test/ 下跑宏都能命中）
//   3. 可执行文件旁边的 ../share/spyralpp/data（cmake --install 之后的布局）
//
// 判据是该目录下存在 amdc_2020。
//----------------------------------------------------------------------

#include <string>

namespace rk4 {

// 返回找到的数据目录（不含结尾 '/'）；都找不到时返回 "data"。
const std::string& DataDir();

// 返回 DataDir() + "/" + name
std::string DataFile(const std::string& name);

} // namespace rk4

#endif // RK4_CORE_DATAPATH_H
