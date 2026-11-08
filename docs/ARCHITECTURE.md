# SpyralPP 架构说明

## 分层

```
core/   物理内核（不依赖 ROOT、不依赖 catima）
  ├── Units.h          单位与常数（唯一来源）
  ├── Element.h        化合物元素项
  ├── Rng.h            PCG32 随机数（可复现）
  ├── Event.h          Particle / Event / FourMomentum
  ├── Distributions.h  激发能、极角分布
  ├── DataPath         数据目录解析（环境变量 / 源码树 / 安装目录）
  ├── MassData         核素质量表（AMDC 2020，缓存）
  └── Reaction         二体反应、二体衰变

sim/    粒子与探测器（catima 静态链入）
  ├── TargetDensity    气体/固体/混合靶
  ├── Stopping         catima 封装：剩余能量、dE/dx
  ├── KinematicsGenerator  顶点抽样 + 靶中能损 + 拒绝采样
  ├── ParticleDynamics 运动方程、RK4 积分、停止条件
  ├── DetectorConfig   探测器与电子学参数
  ├── PadPlane         pad 几何、时间修正、桶网格查找
  ├── BeamPads         束流区 pad（由 attpc_engine 生成）
  ├── Response         GET 响应函数
  └── Simulator        事件 -> 点云

reco/   重建
  ├── WaveProcessor    FFT sinc 基线扣除
  ├── PulseShape       GET 响应模板（attpc_engine 同款公式）+ 饱和修复 + 模板拟合（FFTW 静态链入）
  ├── PeakFinder       局部极大 + 突出度 + 宽度
  ├── TraceReader      文本 / CSV 波形读取
  ├── H5TraceReader    真实 GET 数据（HDF5 `/get/evt{N}_data`）
  ├── PointCloudBuilder 波形 -> 点云
  ├── Cluster / Clusterer  默认 HDBSCAN → 圆重叠合并 → LOF；显式旧 DBSCAN 可选
  ├── Hdbscan          HDBSCAN（默认后端；显式点数限制超限报错）
  ├── Spline           P-spline 平滑
  ├── Geometry         最小二乘圆拟合
  ├── Estimator        顶点 / 极角 / Bρ / dE/dx
  ├── ParticleId       Bethe 近似 + 候选判别
  ├── TrackMesh        (ke, polar) 网格 + 双线性插值 + 缓存（Minuit 后端用）
  ├── KalmanFitter     UKF + RTS 平滑（可选后端，catima 能损直接进过程模型）
  └── TrackFitter      TMinuit 两阶段拟合

io/     输入输出
  ├── Config           key=value 配置
  ├── TreeIO           运动学事件 TTree
  ├── CloudIO          点云 TTree
  ├── ResultIO         拟合结果 TTree
  ├── StageIO          中间过程（cluster / estimate / 峰参数）
  └── H5Writer         HDF5（内置 extern/hdf5 + extern/highfive）

scripts/ 运维脚本（bash，无额外依赖）
  ├── spyralpp.sh       单个 run：simulate / analyze（分片并行）/ merge / status
  └── run_all.sh        读 configs/runs.txt 连续跑，跳过已完成，失败不中断，结尾汇总

app/    命令行程序
  ├── simulate         config -> 事件 -> 点云 + 真值
  ├── waveform         波形 -> 点云
  ├── csvtrace         CSV 波形 -> 点云
  ├── analyze          点云 -> 聚类 -> 估计 -> 拟合
  └── merge            分片结果合并
```

## 第三方库（extern/）

第三方依赖全部以源码形式随仓库分发，`CMakeLists.txt` 用 `add_subdirectory()` 编成静态库，
再静态链接进 `librk4.so`。因此运行时除了 ROOT 没有其他动态依赖，也不会误用系统里的同名库。

| 目录 | 版本 | 用法 |
|---|---|---|
| `extern/catima` | 1.8.0 | 保留上游 `CMakeLists.txt`，编成 `libcatima.a` |
| `extern/fftw` | 3.3.10 | 保留自带 CMakeLists（双精度、SSE2/AVX/AVX2），编成 `libfftw3.a` |
| `extern/hdf5` | 1.12.0 | 预编译静态库 + 头文件 |
| `extern/highfive` | header-only | HDF5 的 C++ 封装 |
| `extern/rk4_solver` | header-only | RK4 积分器 |
| `extern/openkf` | GPL-3.0 | UKF 的 sigma 点 / 加权均值协方差 / 互相关工具 |
| `extern/eigen` | 3.4.0，MPL2 | 矩阵库（header-only，OpenKF 依赖） |

`extern/openkf` 与 `extern/eigen` 的出处、收录文件、许可证见各自的 `NOTICE.md`。

两个子项目都用 `EXCLUDE_FROM_ALL` 加入，构建时照常参与，但它们自己的 install 规则不生效，
所以 `cmake --install` 只装本项目的东西。

## 数据流

```
config ──► rk4_simulate ──► clouds.root + clouds_kin.root
                                  │
waveform / CSV / HDF5 ──► rk4_waveform / rk4_csvtrace ──► clouds.root
                                  │
                                  ▼
                            rk4_analyze ──► fit.root
                                            (brho/polar/vertex/pid/cost)
```

## 关键约定

| 项 | 约定 |
|---|---|
| 单位 | 内部 m / s / MeV / rad / T / V/m；mm 只在 pad 几何、点云、拟合边界出现 |
| 质量 | `GetNuclearMassNuclear()` = 原子质量 − Z·mₑ |
| 坐标 | 束流从 z=0 入射，pad 面在 z=length |
| 时间桶 | `tb = (length − z)/dv + micromegas_edge`，`dv = length/(windows−micromegas)` |
| 网格参考点 | `(0, 0, z_ref)`，默认 z_ref = 0.5 m，使用时平移到实际顶点 |
| 事件索引 | `[靶, 束流, 出射, 剩余]`，后续衰变 4,5 / 6,7 … |
| 数据目录 | `DataDir()`：`SPYRALPP_DATA_DIR` → `data/`、`./`、`../data/` → 安装后的 `share/spyralpp/data` |

## ParticleDynamics 的边界

`sim/ParticleDynamics` 只负责：运动方程、RK4 积分、停止条件、靶/核素配置、二体反应与衰变。

文件里另外约 1500 行是遗留的实验分析代码（DWBA/截面/绘图），用 `RK4_WITH_LEGACY_ANALYSIS`
编译开关隔离，默认不参与构建。默认配置下该编译单元的目标文件从 420 KB 降到 49 KB。

## 并行模型

- **线程**：`rk4_analyze` 用 `std::thread` 池按事件并行；每线程一套独立处理对象，
  共享只读 `TrackMesh`。TMinuit 的回调指针是 `thread_local`。
- **共享内存**：`TrackMesh::LoadShared` 用 `mmap(MAP_SHARED)` 映射网格文件，
  多进程共享同一份物理页；`Count()`/`At()` 屏蔽了"自有 vector"与"映射视图"两种后端。
- **多进程**：`rk4_analyze ... <shard_index> <shard_count>` 只处理
  `idx % shard_count == shard_index` 的点云；各进程用 `LoadShared` 映射同一网格，
  结果由 `rk4_merge` 排序去重合并。
- **波形并行**：`rk4_waveform` 是生产者/消费者——主线程读 HDF5（HDF5 默认非线程安全），
  N 个 worker 各自一套 `PointCloudBuilder`；结果按事件号排序，输出与单线程逐位一致。
- **线程安全约定**（踩过坑后的硬规矩）：
  - `MassData::Load()` 用函数内静态锁串行化（首次加载会写静态 map）；
  - `ParticleDynamics` 的构造 / `ConfigureTarget` / `SetParticle` 串行化（catima 全局存储 +
    ROOT 对象注册表都不是线程安全的），`equationMotion` 保持并行；
  - FFTW 只串行化"规划"（`fftw_plan_*`、`fftw_destroy_plan`），`fftw_execute` 并行；
  - 寻峰的亚采样抖动按事件号定种子，避免线程调度改变物理结果；
  - `TMinuit` 通过全局 `gMinuit` 桥接 Fortran MINUIT，`TrackFitter::Fit` 用全局锁串行化
    （Minuit 后端因此不随线程数加速；并行请用 Kalman 后端）。
- **模拟并行**：`rk4_simulate` 每线程独立的生成器+模拟器；事件编号用全局请求索引
  （避免每线程计数器重复）。
- **结果确定性**：
  - `rk4_analyze` 并行只改变处理顺序，不改变数值（测试：单线程 vs 4 线程、
    单进程 vs 4 进程分片，两组逐位比对）；
  - `rk4_simulate` 并行会改变生成的事件集合（每线程随机子流不同），因此默认 1 线程；
    开启后结果依赖 `threads` 值，但编号唯一、真值一一对应（`testSimParallel` 保护）。

## I/O 与格式

| 数据 | 格式 | 位置 |
|---|---|---|
| 运动学事件 | TTree / HDF5 | `io/TreeIO`、`io/H5Writer` |
| 点云 | TTree / HDF5 | `io/CloudIO`、`io/H5Writer` |
| 拟合结果 | TTree | `io/ResultIO` |
| 原始波形（真实数据） | HDF5 `/get/evt{N}_data`（merger 1.x 表格式 / 两种老布局自动识别） | `reco/H5TraceReader` |
| 原始波形（示例/导出） | 文本、CSV | `reco/TraceReader`、`app/csvtrace` |

HDF5 的两套布局都与 Python 生态对齐：

- 运动学：`/data/chunk_N/event_M`（对应 `attpc_engine` 的 KinematicsPipeline）
- 点云：`/cloud/cloud_N` + `/labels_N`（对应 `attpc_engine` 的 SpyralWriter，Python 版 Spyral 可直接读）

## 验收

干净重建 → 全量测试 → 安装 → 端到端实跑，均通过：

```
configure rc=0   build rc=0   warnings=0        （首次构建约 20 s，含编译外部库）
测试 21/21 通过
安装产物：bin/（5 个程序）、lib/librk4.so、include/spyralpp/、share/spyralpp/
ldd librk4.so 只剩 ROOT：无 libcatima / libfftw3 / libhdf5 动态依赖
安装后的 rk4_simulate 在任意目录可直接运行（自动定位 share/spyralpp/data）
模拟 8 事件（2 线程 + HDF5 双输出）→ 重建 8 条拟合 → 第二次复用共享网格缓存
```

## 扩展点

- **换积分器**：`ParticleDynamics::Propagate` 是唯一入口，内部用 `rk4_solver`。
- **换聚类算法**：`Clusterer::Run` 输出 `std::vector<Cluster>`，替换实现不影响下游。
- **换拟合器**：`solver.fitter` 在 `minuit`（`TrackFitter`，轨迹网格 + TMinuit）与
  `kalman`（`KalmanFitter`，UKF + RTS）之间切换，输出字段一致；两条路径都在
  `app/analyze.cpp` 里分发。
- **换聚类后端**：`ClusterParams::use_hdbscan` 在 DBSCAN 与 HDBSCAN 之间切换，
  HDBSCAN 使用 `FinalizeSpyral`（圆重叠合并 / z 排序 / LOF / 编号）；
  显式 DBSCAN 使用旧 `Finalize`。不能将两种清理流程混为同一个实验。
- **加粒子种类**：`ParticleId::SetCandidates`。
- **加 I/O 格式**：`io/` 下新增读写器；`H5Writer` 已留出 `RK4_WITH_HDF5` 开关。

## 测试

`test/` 下是 ROOT 宏（无需单独编译），`bash test/run_tests.sh` 统一运行：

| 宏 | 覆盖 |
|---|---|
| testKinematics | 四动量守恒、在壳、分布、能损、激发、阈值、衰变 |
| testIO | 运动学 TTree 往返 |
| testDataPath | 数据目录解析、质量表加载 |
| testHdf5 | HDF5 输出（未启用时跳过） |
| testDetector | 电子数、扩散、pad、响应、ADC、点云 I/O |
| testWaveform | 波形 -> 点云（含 CR 行尾 CSV） |
| testCluster | 圆拟合、单径迹、切碎后合并、不误合并 |
| testSolver | 网格插值精度、两阶段拟合精度 |
| testKalman | UKF+RTS 拟合：与真值、与 Minuit 对比（Bρ、θ、χ²/ndf） |
| testParticleId | dE/dx 判别与 d–α 简并 |
| testMeshCache | 网格缓存参数匹配与失效 |
| testPipeline | 模拟 -> 点云 -> 重建，与真值比对 |
| testWaveformPipeline | 真实波形 -> 点云 -> 重建 |
| testRK4 | 单条轨迹积分（原始用例） |
