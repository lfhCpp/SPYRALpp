# OpenKF（vendored）

来源：<https://github.com/hippo91/OpenKF> 的 UKF/KF 头文件，经 ATTPCROOTv2
`OpenKF-Claude` 分支（<https://github.com/Yassid/ATTPCROOTv2>，commit
`fbd3695a70cd8280d4b4ab8a38758dafc118e4f7`）转引，路径
`AtReconstruction/AtFitter/OpenKF/`。

收录的文件（未做修改）：

| 文件 | 用途 |
|---|---|
| `types.h` | `kf::Matrix` / `kf::Vector` 的 Eigen 别名 |
| `kf_util.h` | sigma 点矩阵的列操作工具（`copyToColumn` / `getColumnAt` / …） |
| `kalman_filter/kalman_filter.h` | 线性 KF 基类（状态、协方差、`predictLKF`/`correctLKF`） |
| `kalman_filter/unscented_kalman_filter.h` | UKF（sigma 点、加权均值/协方差、互相关） |
| `kalman_filter/unscented_transform.h` | 无迹变换的独立实现 |

许可证：GPL-3.0（见 `COPYING`，每个头文件里也有声明）。作者
Mohanad Youssef (Al-khwarizmi)。本项目的 `src/reco/KalmanFitter.C` 按同一套
sigma 点/权重/互相关写法实现了轨迹拟合所需的 UKF + RTS 平滑，改动之处：
增广状态只含"状态 + dE/dx 缩放"（量测噪声在校正时单独加），并加了模型噪声
与平滑器。
