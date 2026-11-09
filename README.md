# SpyralPP

C++ reconstruction pipeline for AT-TPC data, organized into four stages:

`Pointcloud → Cluster → Estimation → InterpSolver`

The reconstruction algorithms and ROOT/HDF5 I/O run in C++. The Python entry point handles configuration, stage scheduling, caching, and Parquet conversion.

## Build

Requirements: C++17, CMake, ROOT with Minuit and Minuit2, HDF5 development files, and Python 3. Install the Python Parquet dependency with `python3 -m pip install -r requirements.txt`. CATIMA, FFTW, and other source dependencies are included under `extern/`.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

If ROOT is not discoverable, initialize its environment or set `ROOTSYS` / `CMAKE_PREFIX_PATH` before configuring.

## Configure and run

Create a local configuration from the portable template, then set the paths for your input data and optional PID file:

```bash
cp configs/local.cfg.example configs/local.cfg
# Edit configs/local.cfg
./spyral
```

Alternatively, pass the configuration explicitly:

```bash
./spyral configs/local.cfg
```

The same configuration can generate simulated events:

```bash
./spyral simulate configs/local.cfg
```

Simulation settings live in `[simulation]`, `[reaction]`, `[beam]`, and
`[vertex]`. The default output is `workspace/Simulation/simulated.root` plus
`simulated_kin.root`; `[simulation] output = ...` can select another ROOT file.
Use `--events N` to override the configured event count and `--force` to
replace existing outputs. The output directory is created automatically.
Simulation writes ROOT trees by default; set `simulation.hdf5 = true` for
additional HDF5 copies.

`configs/local.cfg` is machine-specific and ignored by Git. The example configuration is safe to commit.

### Main settings

```ini
[pipeline]
pointcloud = true
cluster = true
estimation = true
interp_solver = false

[run]
start = 151
end = 160
skip_missing = true
threads = 8

[input]
trace_dir = /path/to/run_hdf5_files

[output]
workspace = workspace
format = python

[pointcloud]
waveform_fit = false
```

The run interval includes both endpoints; missing `run_NNNN.h5` files are skipped when `skip_missing = true`. Each pipeline stage can be enabled independently. Disabled stages reuse existing workspace outputs when downstream stages need them.

`waveform_fit = true` enables GET response-template waveform fitting, including saturation recovery. `false` uses the faster baseline peak finder.

`format = python` writes Pointcloud/Cluster as HDF5 and Estimation/InterpSolver as Parquet. `format = root` writes all stages as ROOT. Missing workspace and stage directories are created automatically. The command-line options `--threads`, `--waveform-fit`, `--no-waveform-fit`, and `--format` override the configuration for one invocation.

Progress is printed to the terminal and saved under `workspace/.spyralpp/`. Example:

```bash
./spyral --threads 8 --events 100
```

## Stage files

```text
workspace/
├── Pointcloud/     run_0151.h5
├── Cluster/        run_0151.h5
├── Estimation/     run_0151.parquet
└── InterpSolver/   run_0151_1H.parquet
```

Each stage can read/write ROOT files too. See [workspace formats](docs/WORKSPACE_FORMATS.md), [pipeline architecture](docs/ARCHITECTURE.md), and [running guide](docs/RUNNING.md).
