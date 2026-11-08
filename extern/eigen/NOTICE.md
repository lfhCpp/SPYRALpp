# Eigen（vendored）

来源：<https://gitlab.com/libeigen/eigen>，版本 3.4.0，只收录了 `Eigen/` 头文件
（OpenKF 的 `types.h` 用 `#include <Eigen/Dense>`）。

许可证：MPL2（`COPYING.MPL2`），另见上游的 `COPYING.*`。

用法：`CMakeLists.txt` 把 `extern/eigen` 加进 `rk4` 的公共头文件路径，
因此 `#include <Eigen/Dense>` 直接可用，不依赖系统装没装 Eigen。
