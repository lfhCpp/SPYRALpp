# Workspace 阶段接口

公共入口为项目根目录 `./spyralpp`。`scripts/spyralpp` 也转发到该入口，旧配置优先的命令仍由内部 legacy 脚本支持。历史结果保留，不自动搬动。

| 阶段 | 输入 | Python 输出 | ROOT 输出 / TTree |
|---|---|---|---|
| Pointcloud | 原始 GET HDF5 或项目文本波形 | `Pointcloud/run_0151.h5` | `Pointcloud/run_0151.root` / `cloud` |
| Cluster | Pointcloud HDF5 / ROOT | `Cluster/run_0151.h5` | `Cluster/run_0151.root` / `cluster` + `events` |
| Estimation | Cluster HDF5 / ROOT | `Estimation/run_0151.parquet` | `Estimation/run_0151.root` / `estimation` |
| InterpSolver | Cluster HDF5 / ROOT + Estimation Parquet / ROOT | `InterpSolver/run_0151_1H.parquet` | `InterpSolver/run_0151_1H.root` / `fit` |

原始波形 ROOT 的 schema 没有统一约定；当前 ROOT 接口从 Pointcloud 阶段产物开始，不能把任意 ROOT 波形文件当成 GET HDF5 读入。

## 独立接续

`--from` 是开始**计算**的阶段，输入是该阶段所需的上游文件。省略 `--input` 时自动从同一 workspace 读取，优先所选输出格式，另一种格式也可直接读取。

```bash
# Python 点云 -> C++ 簇 + 估计（结果独立保存）
./spyralpp 151 --from Cluster --input /old/workspace/Pointcloud/run_0151.h5 --workspace cpp_workspace

# ROOT 簇 -> Parquet 估计
./spyralpp 151 --from Estimation --input clusters.root --workspace cpp_workspace

# Python 估计与簇 -> C++ Kalman；不调用 Clusterer 或 Estimator
./spyralpp 151 --from InterpSolver --to InterpSolver \
  --input /old/workspace/Estimation/run_0151.parquet \
  --clusters /old/workspace/Cluster/run_0151.h5 \
  --config /path/to/main.py --fitter kalman --workspace cpp_workspace

# 单独导出；反向转换只需交换输入、输出格式
./spyralpp export Pointcloud clouds.h5 clouds.root
./spyralpp export Cluster clusters.root clusters.h5
./spyralpp export Estimation estimates.root estimates.parquet
./spyralpp export InterpSolver fitted.parquet fitted.root
```

入口不会执行 main.py 的导入、副作用、Pipeline 或批量 run 范围，只读取常量参数声明。动态配置表达式、Garfield 校正、continuity_join 等尚不支持的导入配置会明确报错。自定义 pad CSV 按原文件复制到临时目录使用。当前 target.json 导入支持 `compound` 气体格式（顺序为 Z,A,计数）。

## 数据契约

Pointcloud HDF5 为 `/cloud/cloud_<event>`，每行 8 列：x,y,z,amplitude,integral,pad,time,scale。Cluster HDF5 为 `/cluster/event_<event>/cluster_<index>/cloud`，每行 5 列：x,y,z,integral,scale。原事件与 IC 属性沿用 Python；簇组的 `label` 是原 HDBSCAN label，不同于连续的 cluster_index。保留空点云和零簇事件。

ROOT Pointcloud tree 使用 event、npoints、x/y/z、amplitude、integral、pad、tb、scale；向量一行一个事件。Cluster tree 一行一个簇，包含 event、cluster_index、cluster_label、npoints、direction、x/y/z/integral/scale；`events` tree 保存零簇事件。两者均携带 orig_run、orig_event、ic_amplitude、ic_centroid、ic_integral、ic_multiplicity。Pointcloud 与 Cluster 坐标为 **mm**。

Estimation 的 23 列严格沿用原 Python：event、cluster_index、cluster_label、orig_run、orig_event、四个 IC 字段、vertex_x/y/z、center_x/y/z、polar、azimuthal、brho、dEdx、sqrt_dEdx、dE、arclength、direction。ID 和 direction 是 Int64，其余 Float64；仅写有效估计，空表也保留 schema。位置 **mm**，角度 **rad**。

InterpSolver 使用 Python SolverResult 的 20 列：五个 ID、vertex_x/sigma_vx、vertex_y/sigma_vy、vertex_z/sigma_vz、brho/sigma_brho、ke/sigma_ke、polar/sigma_polar、azimuthal/sigma_azimuthal、redchisq。位置为 **m**，能量 MeV，角度 rad，Bρ 为 T·m。附加 solver_method、objective_kind 两列。未知协方差写 NaN/Parquet null，不能当作精确测量误差。

`redchisq` 沿用 Python 兼容字段名：Kalman 记录 χ²/ndf；原生 Minuit2 记录平均最近距离 mm；`python-objective` 记录 Python 加权最近点距离 m。由 `objective_kind` 标记定义，三种统计量不可直接混比。`python-objective` 使用 Python 官方 `.npy` 网格、Python 权重和参数范围，但由 Minuit2 优化；它不等同于 Python 的 SciPy L-BFGS-B。Python 参考目标对每个点独立求整条轨迹的最近样本，不保留击中顺序，也没有基于该目标的质量拒绝门；`redchisq` 是 lmfit 对标量距离目标计算的 `redchi`，不能按常规 χ²/ndf 解读。命令行通过 `--fitter python-objective --python-mesh <mesh.npy>` 启用。当前三个 Minuit2 路径都要求 fit_vertex_rho/phi/azimuthal 全为 true；Python `fit_method` 不决定 C++ 后端。

## C++ 调用

头文件为 [`src/io/WorkspaceIO.h`](../src/io/WorkspaceIO.h)，命名空间 rk4：

```cpp
auto clouds = ReadPointcloudStage("Pointcloud/run_0151.h5");
WritePointcloudStage("Pointcloud/run_0151.root", clouds);
auto clusters = ReadClusterStage("Cluster/run_0151.root");
WriteClusterStage("Cluster/run_0151.h5", clusters);
auto estimates = ReadScalarStage("estimates.root", "Estimation");
WriteScalarStage("estimates.csv", "Estimation", estimates);
```

编译后独立引擎：

```bash
build/rk4_workspace cluster clouds.root clusters.root configs/local.cfg
build/rk4_workspace estimate clusters.h5 estimates.root configs/local.cfg
build/rk4_workspace solve clusters.root estimates.root fitted.root workspace/.spyralpp/run_0151_config.cfg
build/rk4_workspace convert Cluster clusters.h5 clusters.root
```

原生引擎负责 ROOT/HDF5/内部 CSV；用户入口负责 Parquet 边界以及 target/PID JSON 和 main.py 导入，使用这些配置时优先通过入口。native solve 的 cfg 应为入口生成的 `.spyralpp/run_0151_config.cfg`，其中已展开物理材料及门参数。

## 复用与写出

`.spyralpp/run_0151.json` 记录各阶段输入身份、参数、二进制 SHA256、输出大小和修改时间。原始大 H5 用路径、大小和纳秒修改时间标识，未对整份原始数据做内容散列。阶段专属参数变化会重算该阶段及其依赖；指定 `--from` 时以所提供的已有文件为起点。显式切换输出格式复用已知等价的上游文件，仅做转换。

输出先写同一 workspace 的隐藏临时目录，完成后原子替换；失败不产生半成品。每个 run 加文件锁，避免同一 workspace 同时写同一 run。计时包含文件读写与进程启动；复用记录为 0 秒，不作为处理性能。
