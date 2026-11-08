# For Developers

SpyralPP 是 AT-TPC 模拟与重建流水线的 C++ 实现。欢迎参与改进！如果你希望加入协作，请联系维护者，或直接提交 Pull Request。下面是在本项目中工作的约定。

## Required Development Tools

### 编译器与构建系统

需要支持 C++17 的编译器（GCC ≥ 9 / Clang ≥ 10）与 CMake ≥ 3.16。

### ROOT

使用 ROOT 6.30 或更高版本。本项目依赖 ROOT 的 `TTree`、`TF1`、`TGraph`、`TMinuit`。
注意：`TMinuit` 不在 `root-config --libs` 的默认列表里，CMake 会单独查找 `libMinuit`。

### 内置依赖（extern/）

除了 ROOT，本项目**不要求任何预装的第三方库**。所有依赖都以源码（或预编译静态库）形式放在 `extern/` 下，
由 `CMakeLists.txt` 直接 `add_subdirectory()` 编译，产物静态链接进 `librk4.so`：

| 依赖 | 位置 | 说明 |
|---|---|---|
| catima 1.8.0 | `extern/catima/` | 阻止本领 / 射程，上游源码随仓库分发，编成静态库（GPL-3.0） |
| FFTW 3.3.10 | `extern/fftw/` | 双精度 FFT，保留其自带 CMakeLists，编成静态库（GPL-2.0-or-later） |
| rk4_solver | `extern/rk4_solver/` | header-only 的 RK4 积分器（MIT） |
| HighFive | `extern/highfive/` | header-only 的 HDF5 C++ 封装（BSL-1.0） |
| HDF5 | `extern/hdf5/` | 从源码编译的 1.12.0 静态库与头文件 |
| OpenKF | `extern/openkf/` | UKF 的 sigma 点/互相关工具（GPL-3.0，见 `NOTICE.md`） |
| Eigen 3.4.0 | `extern/eigen/` | header-only 矩阵库（MPL2，OpenKF 依赖） |

注意：即使系统里装了 `/usr/local/lib/libcatima.so` 或 `libfftw3`，构建也**不会**用到它们；
升级这两个库就是替换 `extern/` 下对应的目录。

### 推荐的构建命令

```bash
export ROOTSYS=/path/to/root
export PATH=$ROOTSYS/bin:$PATH

cmake -S . -B build -DCMAKE_PREFIX_PATH=$ROOTSYS
cmake --build build -j

bash test/run_tests.sh
```

## Pull Requests

### Statement of Intent

请在 Pull Request 中说明这次修改的目标是什么，以及为此改动了什么。没有意图说明的 PR 不会被接受。

### Statement of Dependencies

请说明是否新增了依赖。本项目的基本原则是**尽量不引入新依赖**：能自己写的小工具就自己写，必须引入时优先考虑 header-only 或可 vendor 的库，并把它放进 `extern/`，同时在本文档的「内置依赖」表格里登记。

新增依赖必须与 GPL-3.0 兼容。

### Review

任何 Pull Request 至少需要一位 reviewer。

## Code Requirements

下面是对加进 SpyralPP 的代码的要求。

### 头文件自洽

每个 `.h` 必须自己 include 它需要的头文件，不能依赖调用方的包含顺序。
这条是硬性要求，因为同一份代码要同时满足两种使用方式：

- 被 ROOT 宏 `#include` 后由 Cling 解释执行
- 被 CMake 编译成共享库

### 命名与风格

- 类名 `PascalCase`，方法 `GetXxx` / `SetXxx` / `CalcXxx`，成员变量 `lower_snake_case`
- `.h` 与 `.C` 成对出现，实现放在 `.C`
- 单位：内部统一 `m / s / MeV / rad / T / V/m`；mm 只出现在 pad 几何、点云与拟合边界
- 数学与物理量的含义写在注释里，不要靠变量名猜

### 格式化

以仓库根的 [`.clang-format`](.clang-format) 为准（4 空格、100 列、指针靠类型）：

```bash
clang-format -i src/reco/KalmanFitter.C
```

只格式化你改动的文件，不要提交"整仓重排"这种 diff——那会让 review 无从下手。

### 注释

请在关键算法处写中文注释，说明**为什么**这样做，而不是复述代码在做什么。
尤其是那些"看起来可以简化、但简化后会出错"的地方（例如 `Clusterer` 里的方向判定、`TrackMesh` 的头解析）。

### 文件大小

请尽量不要写出上千行的巨型文件。一个文件应当承载一个相对独立的单元。
如果某个类开始承担多种职责，考虑用编译开关或拆分子模块（`ParticleDynamics` 的遗留分析代码就用
`SPYRALPP_WITH_LEGACY_ANALYSIS` 隔离）。

### 类型与数值

- 浮点比较用容差，不要用 `==`
- 涉及坐标变换时，明确坐标系（实验室系 / 质心系 / pad 平面）
- 新增物理量时，先在 `core/Units.h` 里确认单位约定

### Unit Tests

请在 `test/` 下为新代码补充测试。测试是 ROOT 宏，**不需要单独编译**：

```bash
root -l -b -q test/testXxx.C
```

宏的固定写法：

```cpp
R__ADD_INCLUDE_PATH(src)
R__LOAD_LIBRARY(build/librk4.so)
#include "core/Event.h"
```

测试应打印 `[FAIL]` 并最终以非零退出码结束（用 `gSystem->Exit(1)`），
`test/run_tests.sh` 会据此汇总结果。

新增测试后，记得把它加入 `test/run_tests.sh` 的 `TESTS` 数组。

### Documentation

新增主要功能时，请补充文档：

- 用户向的用法与配置项 → `README.md`
- 模块职责、数据流、扩展点 → `docs/ARCHITECTURE.md`

配置项要写清楚**名称、单位、默认值**以及它对结果的影响。

## Final Thoughts

下面是一个 Pull Request 描述的示例：

```text
Intent:
  修复 TrackMesh 在共享映射模式下读取头部时的偏移错误。
  当二进制数据首字节恰为 0x0A 时会多跳一个字节，导致插值轨迹整体错位。

Dependencies Added:
  无
```

你可以在提交前对照本文件自查一遍。如果某项约定在本项目中确实不适用，请在 PR 中说明理由。

如果你需要对某个具体实验做深度定制，建议 fork 本仓库，这样可以保留版本控制与协作能力，
而不必让改动符合本项目「实验无关」的定位。
