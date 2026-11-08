# 怎么跑

这份是"操作手册"：从零到跑完一批 run 的全部命令，以及跑挂/想重跑时怎么办。
（方法细节见 [ARCHITECTURE.md](ARCHITECTURE.md)，参数逐项解释见 [README](../README.md)。）

## 0. 一次性准备（只做一次）

```bash
# ROOT：每次新开 shell 都要 source
export ROOTSYS=/path/to/root
source $ROOTSYS/bin/thisroot.sh

# 构建（约 20 s，含从源码编译 catima/FFTW）
cd SpyralPP
cmake -S . -B build -DCMAKE_PREFIX_PATH=$ROOTSYS
cmake --build build -j

# 可选：安装到别处，之后用 /opt/spyralpp/bin/rk4_* 也一样
cmake --install build --prefix /opt/spyralpp
```

自检（跑一遍，21 个宏测试，约 2 分钟）：

```bash
bash test/run_tests.sh
```

## 1. 先手动跑一次（3 条命令）

```bash
./build/rk4_simulate configs/example.cfg /tmp/demo/clouds.root 50
./build/rk4_analyze  /tmp/demo/clouds.root /tmp/demo/fit.root configs/example.cfg
```

想分片并行（4 进程处理同一份点云，结果用 `rk4_merge` 合并）：

```bash
for i in 0 1 2 3; do
  ./build/rk4_analyze /tmp/demo/clouds.root /tmp/demo/fit_shard_$i.root \
      configs/example.cfg $i 4 &
done
wait
./build/rk4_merge /tmp/demo/fit.root /tmp/demo/fit_shard_*.root
```

## 2. 一条命令跑完一个 run（推荐日常用法）

`scripts/spyralpp.sh` 把上面这些包起来，带日志、跳过已完成产物：

```bash
scripts/spyralpp.sh all demo configs/example.cfg 200
```

它做的事：模拟 → 4 分片并行重建 → 合并 → 产物落在 `output/demo/`：

```
output/demo/clouds.root        点云（TTree）
output/demo/clouds_kin.root    运动学真值
output/demo/fit.root           拟合结果
output/logs/demo.simulate.log  / demo.analyze.0..3.log / demo.merge.log
```

其他子命令：

```bash
scripts/spyralpp.sh simulate demo configs/example.cfg 200   # 只模拟
scripts/spyralpp.sh analyze  demo configs/example.cfg       # 只重建（自动分片+合并）
scripts/spyralpp.sh merge    demo                           # 只合并
scripts/spyralpp.sh status                                  # 所有 run 的进度表
scripts/spyralpp.sh status   demo                           # 单个 run
```

## 3. 一只脚本跑到尾（像 Python 那样配一次就长期跑）

把所有要分析的 run 写进清单 `configs/runs.txt`：

```text
# run_id        config                  n_events
demo            configs/example.cfg     200
h2_600torr      configs/example.cfg     2000
```

然后：

```bash
scripts/run_all.sh                      # 读 configs/runs.txt，一条条跑
scripts/run_all.sh my_runs.txt          # 换清单
scripts/run_all.sh --force              # 全部重跑
```

它的行为（`run_all.sh` 结尾会打印汇总，有失败则退出码非 0）：

- **断点续跑**：`clouds.root` / `fit.root` 已存在的 run 直接跳过，所以中途 Ctrl-C、断网、
  机器重启后重新执行同一条命令即可，不会重复算；
- **失败不中断**：某个 run 挂了会记下来继续跑下一个，最后列出失败的 run_id；
- **日志分离**：每个 run、每个阶段一个日志，失败时脚本自动打印日志尾部 5 行。

想后台挂着跑：

```bash
nohup scripts/run_all.sh configs/runs.txt > output/run_all.log 2>&1 &
```

装在 `<prefix>` 的话，脚本和数据一起被装到 `<prefix>/share/spyralpp/`，直接执行即可：

```bash
SPYRALPP_BIN=/opt/spyralpp/bin \
/opt/spyralpp/share/spyralpp/scripts/run_all.sh \
    /opt/spyralpp/share/spyralpp/configs/runs.txt
```

## 4. 常用开关（环境变量）

| 变量 | 作用 | 默认 |
|---|---|---|
| `SPYRALPP_THREADS` | 每个进程内部用几个线程（覆盖配置里的 `run.threads`） | 0 = 自动 |
| `SPYRALPP_SHARDS` | 重建分几个进程分片 | 1 |
| `SPYRALPP_OUT` | 输出根目录 | `<repo>/output` |
| `SPYRALPP_BIN` | 可执行文件目录（装到 /opt 后指过去即可） | `<repo>/build` |
| `SPYRALPP_FORCE` | 1 = 忽略已有产物重跑（等价 `--force`） | 0 |
| `SPYRALPP_CLEAN_SHARDS` | 1 = 合并后删掉 `fit_shard_*.root` | 0 |
| `SPYRALPP_DATA_DIR` | 数据表目录（默认自动找 `data/` 或安装后的 `share/`） | 自动 |

例：24 核机器上跑 4 分片 × 6 线程：

```bash
SPYRALPP_SHARDS=4 SPYRALPP_THREADS=6 scripts/run_all.sh configs/runs.txt
```

分片时脚本会**先串行跑 0 号分片**（若配置了 `solver.mesh_cache`，这步把网格缓存建好），
其余分片再并行加载；建议分片时配置里加 `solver.mesh_shared = 1`（mmap 零拷贝共享网格）。

## 5. 换拟合后端 / 换参数

改配置，不改命令。`configs/example.cfg` 里关键的几项：

```text
[solver]
fitter = minuit          # 或 kalman（UKF + RTS，不需要轨迹网格）
kalman_meas_sigma = 0.5  # kalman 后端的测点 σ，mm
mesh_cache = /tmp/proton_h2.mesh
mesh_shared = 1
[cluster]
use_hdbscan = true       # 默认；false 显式选择旧 DBSCAN
[run]
threads = 1              # 也可被 SPYRALPP_THREADS 覆盖
```

公共 workspace 入口还支持 `--fitter python-objective --python-mesh <mesh.npy>`：它读取官方 Spyral 网格并复用 Python 目标函数，但以 C++ Minuit2 优化，不是 SciPy L-BFGS-B。

## 6. 输入不是模拟数据时

| 输入 | 命令 |
|---|---|
| 真实 GET 的 HDF5（`/get/evt{N}_data`，merger 1.x 表格式或老布局） | `rk4_waveform <cfg> <file.h5> <out/clouds.root> [max_events] [threads]`，之后照常 `rk4_analyze` |
| 导出的 CSV 波形 | `rk4_csvtrace <cfg> <traces.csv> <out/clouds.root> [max_events]` |
| 已有的点云（比如 Python 版 Spyral 写的） | 直接 `rk4_analyze`，只要 TTree/HDF5 布局与 README 一致 |

把 `clouds.root` 放到 `output/<run_id>/` 下，就能直接 `scripts/spyralpp.sh analyze <run_id> <cfg>`。

真实数据跑批（`run_0151.h5` 实测）：

```bash
# 波形 -> 点云：500 个事件，16 线程；新流程耗时见 FIDELITY_AUDIT.md
rk4_waveform configs/example.cfg /path/run_0151.h5 output/run0151/clouds.root 500 16
# 重建：8 线程
SPYRALPP_THREADS=8 rk4_analyze output/run0151/clouds.root output/run0151/fit.root configs/example.cfg
```

想把中间过程也落盘（查饱和、查峰分裂、离线重拟合）：

```ini
[peaks]
dump = 1        # -> <out>_peaks.root：每条迹线的峰参数（含 saturated 标志）

[output]
stages = 1      # -> <fit>_stages.root：cluster 树（每簇的点）+ estimate 树（粗估参数）
```

想用模板拟合（含饱和修复）和 PID 预筛，在配置里加：

模板拟合仍是可选的诊断功能：Run 0151 尚无逐通道 GET 响应标定，生产默认保持关闭。

```ini
[peaks]
template_fit = 1
sampling_ns = 320     # Run 0151: 3.125 MHz
peaking_ns = 1000
sat_adc = 4095

[pid]
gate = 1
min_sep = 0.3
```

两个必查项：**B 场符号**（这份文件在我们约定下要 `bfield_t = 2.85`，写 −2.85 会一条都拟不上）
和 **Kalman 的 `kalman_meas_sigma`**（真实簇要按分辨率放大，实测 5 mm 比 0.5 mm 好一个量级）。

## 7. 出问题怎么办

| 现象 | 原因 / 处理 |
|---|---|
| `找不到 ROOT` | 新 shell 忘了 `source thisroot.sh`；脚本会直接报错退出 |
| `找不到 data/amdc_2020` | 在别的目录跑且没设 `SPYRALPP_DATA_DIR`；构建树/安装树都会自动找，手工拷贝的场景要显式指 |
| 重建结果 0 条 | 看 `cluster.report = 1` 的诊断行（点数/簇数/后端）；簇太少或初值失败 |
| 某个 run 想重跑 | `scripts/spyralpp.sh all <id> <cfg> <n>` 前先删 `output/<id>/`，或 `SPYRALPP_FORCE=1` |
| 想看某阶段在干嘛 | `tail -f output/logs/<id>.analyze.0.log` |
| 网格缓存坏了 | 删掉配置里 `solver.mesh_cache` 指向的文件（带 `mesh_shared=1` 时尤其注意别半途 Ctrl-C） |

## 8. 产物与下游

- `clouds.root`：`event / npoints / x / y / z / amplitude / integral / pad / tb / scale`
- `fit.root`：`event / cluster / brho / sigma_brho / polar / azimuthal / vx / vy / vz / cost /
  chi2ndf / ok / pid_index / pid_sep / dedx`
- 配置里开 `[output] hdf5 = 1` 时会额外产出与 `attpc_engine` 兼容的 HDF5
  （Python 版 Spyral 可直接读）。

用 ROOT 看一眼：

```bash
root -l output/demo/fit.root
root [0] fit->Draw("brho")
```
