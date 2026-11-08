"""SpyralPP workspace entry. Physics runs in native binaries; no main.py execution."""
from __future__ import annotations
import argparse
import ast
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import shlex
import signal
import subprocess
import sys
import tempfile
import time

PROJECT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("SPYRALPP_BIN", PROJECT / "build"))
ENGINE = BIN / "rk4_workspace"
STAGES = ("Pointcloud", "Cluster", "Estimation", "InterpSolver")
PARAMS = {"PadParameters", "GetParameters", "FribParameters", "DetectorParameters", "ClusterParameters", "OverlapJoinParameters", "EstimateParameters", "SolverParameters"}

def say(text):
    print(text, flush=True)

def run_logged(command, log):
    """Keep a child process log while showing its progress in the terminal."""
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace", bufsize=1,
    )
    try:
        with open(log, "a", encoding="utf-8") as handle:
            for line in process.stdout:
                handle.write(line)
                handle.flush()
                sys.stdout.write(line)
                sys.stdout.flush()
        return process.wait()
    except KeyboardInterrupt:
        if process.poll() is None:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
        raise
    finally:
        if process.stdout is not None:
            process.stdout.close()

def native(*args, log=None):
    command = [str(ENGINE), *map(str, args)]
    if log:
        status = run_logged(command, log)
        if status:
            raise RuntimeError(f"阶段失败，日志：{log}\n{Path(log).read_text()[-3500:]}")
        return None
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    return result.stdout

def literal(node, values):
    """Restricted AST interpreter: literal parameter declarations, never eval/import."""
    if isinstance(node, ast.Constant): return node.value
    if isinstance(node, ast.Name):
        if node.id == "DEFAULT_MAP": return "DEFAULT_MAP"
        if node.id in values: return values[node.id]
    if isinstance(node, (ast.List, ast.Tuple)): return [literal(x, values) for x in node.elts]
    if isinstance(node, ast.Dict): return {literal(k, values): literal(v, values) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        v = literal(node.operand, values)
        return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.BinOp):
        a, b = literal(node.left, values), literal(node.right, values)
        if isinstance(node.op, ast.Div) and isinstance(a, str) and isinstance(b, str): return str(Path(a) / b)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Div): return a / b
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        name = node.func.id
        if name == "Path" and len(node.args) == 1 and not node.keywords: return str(literal(node.args[0], values))
        if name in PARAMS and not node.args:
            if any(k.arg is None for k in node.keywords): raise ValueError("不支持参数 ** 展开")
            return {"__type__": name, **{k.arg: literal(k.value, values) for k in node.keywords}}
    raise ValueError("不支持动态配置表达式；请使用常量参数或 .cfg 配置")

def python_config(path):
    values = {}
    tree = ast.parse(path.read_text(), filename=str(path))
    required = {"get_params", "frib_params", "det_params", "cluster_params", "estimate_params", "solver_params", "pad_params", "trace_path"}
    for statement in tree.body:
        if isinstance(statement, ast.Assign) and len(statement.targets) == 1 and isinstance(statement.targets[0], ast.Name):
            name = statement.targets[0].id
            try: values[name] = literal(statement.value, values)
            except (ValueError, TypeError, ZeroDivisionError) as exc:
                if name in required: raise ValueError(f"{path}:{statement.lineno} {name}: {exc}") from exc
    found = {v.get("__type__") for v in values.values() if isinstance(v, dict)}
    missing = {"GetParameters", "DetectorParameters", "ClusterParameters", "EstimateParameters"} - found
    if missing: raise ValueError(f"main.py 缺少可静态读取的参数声明：{', '.join(sorted(missing))}")
    cfg = {}
    section_names = {"PadParameters": "pad", "GetParameters": "get", "FribParameters": "frib", "DetectorParameters": "detector", "ClusterParameters": "cluster", "EstimateParameters": "estimate", "SolverParameters": "solver"}
    for value in values.values():
        if isinstance(value, dict) and value.get("__type__") in section_names:
            prefix = section_names[value["__type__"]]
            for k, v in value.items():
                if k == "__type__": continue
                if isinstance(v, dict):
                    for sub, item in v.items():
                        if sub != "__type__": cfg[f"{prefix}.{k}.{sub}"] = item
                elif v is not None: cfg[f"{prefix}.{k}"] = v
    for name in ("trace_path",):
        if name in values: cfg["input.trace_dir"] = values[name]
    for python_name, config_name in (
        ("run_number", "run.number"),
        ("run_start", "run.start"),
        ("run_end", "run.end"),
        ("run_skip_missing", "run.skip_missing"),
        ("pipeline_from", "pipeline.from"),
        ("pipeline_to", "pipeline.to"),
        ("pipeline_pointcloud", "pipeline.pointcloud"),
        ("pipeline_cluster", "pipeline.cluster"),
        ("pipeline_estimation", "pipeline.estimation"),
        ("pipeline_interp_solver", "pipeline.interp_solver"),
    ):
        if python_name in values: cfg[config_name] = values[python_name]
    if truth(cfg.get("detector.do_garfield_correction", False)):
        raise ValueError("当前原生接口尚不支持 Garfield 校正；不能静默忽略此设置")
    if "cluster_params" in values:
        cfg["cluster.enable_join"] = values["cluster_params"].get("overlap_join") is not None
    if values.get("cluster_params", {}).get("continuity_join") is not None:
        raise ValueError("main.py 的 continuity_join 暂不支持；当前默认使用 overlap_join")
    return cfg

def truth(value):
    return str(value).lower() in ("true", "1", "yes", "on")

def write_config(path, cfg):
    lines = []
    for k, v in sorted(cfg.items()):
        if isinstance(v, bool): v = str(v).lower()
        if isinstance(v, (float, int)) and not math.isfinite(v): raise ValueError(f"非有限配置值 {k}")
        v = str(v)
        if "\n" in v or "\r" in v or '"' in v: raise ValueError(f"配置中含不支持的字符：{k}")
        lines.append(f'{k} = "{v}"')
    path.write_text("\n".join(lines) + "\n")

def load_config(path, scratch):
    if path.suffix == ".py":
        cfg = python_config(path)
        temporary = scratch / "import.cfg"
        write_config(temporary, cfg)
        cfg = json.loads(native("config", temporary))
    else: cfg = json.loads(native("config", path))
    aliases = {
        "frib.baseline_window_scale": "frib.baseline_scale", "frib.peak_separation": "frib.separation",
        "frib.peak_prominence": "frib.prominence", "frib.peak_max_width": "frib.max_width", "frib.peak_threshold": "frib.threshold",
        "frib.ic_delay_time_bucket": "frib.ic_delay", "cluster.min_size_lower_cutoff": "cluster.min_cluster_size",
        "cluster.overlap_join.min_cluster_size_join": "cluster.overlap_min_cluster_size",
        "cluster.overlap_join.circle_overlap_ratio": "cluster.overlap_ratio",
    }
    for key, alias in aliases.items():
        if key in cfg: cfg[alias] = cfg[key]
    for key in ("solver.gas_data_path", "solver.particle_id_filename", "input.trace_dir", "detector.pad_dir"):
        if key in cfg:
            p = Path(cfg[key]).expanduser()
            cfg[key] = str((path.parent / p).resolve() if not p.is_absolute() else p)
    if "solver.gas_data_path" in cfg:
        data = json.loads(Path(cfg["solver.gas_data_path"]).read_text())
        compound = data.get("compound")
        if compound is None: raise ValueError("此入口的 target.json 目前需要 compound 气体格式")
        cfg["gas.components"] = len(compound)
        for i, (z, a, s) in enumerate(compound):
            for k, v in (("a", a), ("z", z), ("s", s)): cfg[f"gas.component{i}.{k}"] = v
        cfg["gas.pressure_torr"] = data["pressure(Torr)"]
    if "solver.particle_id_filename" in cfg:
        pid = json.loads(Path(cfg["solver.particle_id_filename"]).read_text())
        cfg["solver.A"], cfg["solver.Z"] = pid["A"], pid["Z"]
        cfg["pid.xaxis"], cfg["pid.yaxis"] = pid["xaxis"], pid["yaxis"]
        cfg["pid.vertices"] = len(pid["vertices"])
        for i, (x, y) in enumerate(pid["vertices"]): cfg[f"pid.vertex{i}.x"], cfg[f"pid.vertex{i}.y"] = x, y
    # Preserve custom pad maps by staging the same CSVs, without editing originals.
    pad_names = {"pad.pad_geometry_path": "padxy.csv", "pad.pad_time_path": "pad_time_correction.csv", "pad.pad_scale_path": "pad_scale.csv"}
    if any(cfg.get(k, "DEFAULT_MAP") != "DEFAULT_MAP" for k in pad_names):
        pad_dir = scratch / "pads"; pad_dir.mkdir()
        for key, filename in pad_names.items():
            source = cfg.get(key, "DEFAULT_MAP")
            source = PROJECT / "data" / filename if source == "DEFAULT_MAP" else path.parent / Path(source)
            if not source.is_file(): raise ValueError(f"找不到 pad map：{source}")
            shutil.copyfile(source, pad_dir / filename)
        cfg["detector.pad_dir"] = str(pad_dir)
    return cfg

def read_config_values(config_path):
    if config_path.suffix == ".py":
        return python_config(config_path)
    return json.loads(native("config", config_path))

def infer_run_number(config_path, cfg=None):
    """Read run.number, or infer it only when the config points to one HDF file."""
    if cfg is None: cfg = read_config_values(config_path)
    configured = cfg.get("run.number")
    if configured is not None:
        try: run = int(configured)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"配置里的 run.number 不是整数：{configured!r}") from exc
        if run < 0: raise ValueError("run.number 不能为负")
        return run

    trace_dir_value = cfg.get("input.trace_dir")
    if not trace_dir_value:
        raise ValueError("配置没有 run.number 或 input.trace_dir；请加 run.number = 151，或使用 --run 151")
    trace_dir = Path(trace_dir_value).expanduser()
    if not trace_dir.is_absolute(): trace_dir = config_path.parent / trace_dir
    if not trace_dir.exists() and str(trace_dir).startswith("/home/DATA/"):
        mapped = Path(str(trace_dir)[5:])
        if mapped.exists(): trace_dir = mapped
    if trace_dir.is_file():
        candidates = [trace_dir]
    elif trace_dir.is_dir():
        candidates = sorted(trace_dir.glob("run_*.h5"))
    else:
        candidates = []
    runs = []
    for candidate in candidates:
        match = re.fullmatch(r"run_(\d+)\.h5", candidate.name)
        if match: runs.append(int(match.group(1)))
    if len(runs) == 1: return runs[0]
    if len(runs) > 1:
        raise ValueError(
            f"{trace_dir} 下有 {len(runs)} 个 run 文件，无法从配置唯一判断要跑哪一个；"
            "请在配置中加 run.number = 151，或使用 --run 151"
        )
    raise ValueError(
        f"在 {trace_dir} 没找到 run_*.h5；请在配置中加 run.number = 151，或使用 --run 151"
    )

def schema(stage): return json.loads(native("schema", stage))

def parquet_to_csv(path, output, stage):
    import polars as pl
    table = pl.read_parquet(path)
    cols = schema(stage)
    # Original Python SolverResult lacks optional native diagnostic text columns.
    for name, kind in cols:
        if name not in table.columns:
            if stage == "InterpSolver" and name in ("solver_method", "objective_kind"):
                table = table.with_columns(pl.lit("python" if name == "solver_method" else "python_objective").alias(name))
            elif stage == "InterpSolver" and name == "particle_hypothesis":
                table = table.with_columns(pl.lit("unknown").alias(name))
            else: raise ValueError(f"{path} 缺少列 {name}")
    table.select([name for name, _ in cols]).write_csv(output, float_scientific=True, float_precision=17, null_value="nan")

def csv_to_parquet(path, output, stage):
    import polars as pl
    types = {"int64": pl.Int64, "float64": pl.Float64, "string": pl.String}
    table = pl.read_csv(path, schema_overrides={name: types[kind] for name, kind in schema(stage)}, null_values=["nan", "NaN"])
    table.write_parquet(output, compression="zstd")

def require_polars():
    try:
        import polars  # noqa: F401
    except ImportError as exc:
        install = (
            f"{shlex.quote(sys.executable)} -m pip install --user -r "
            f"{shlex.quote(str(PROJECT / 'requirements.txt'))}"
        )
        raise RuntimeError(
            f"Parquet 阶段需要 Polars，但当前 Python ({sys.executable}) 未安装。\n"
            f"请先运行：{install}"
        ) from exc

def to_native(path, stage, scratch):
    if path.suffix != ".parquet": return path
    csv_path = scratch / f"input_{stage}.csv"
    parquet_to_csv(path, csv_path, stage)
    return csv_path

def convert(stage, source, target, scratch, log):
    if source.resolve() == target.resolve(): raise ValueError("转换输入与输出必须不同")
    native_input = to_native(source, stage, scratch)
    native_output = scratch / "converted.csv" if target.suffix == ".parquet" else target
    native("convert", stage, native_input, native_output, log=log)
    if target.suffix == ".parquet": csv_to_parquet(native_output, target, stage)

def file_id(path):
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}

def digest(value): return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()

def input_id(path, manifest):
    actual = file_id(path)
    for stage, record in manifest["stages"].items():
        variants = record.get("outputs", [record.get("output", {})])
        if actual in variants:
            return {"stage": stage, "signature": record["signature"]}
    return actual

def binary_id():
    result = {}
    for p in (ENGINE, BIN / "rk4_waveform", BIN / "librk4.so", Path(__file__)):
        result[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    return result

def phase_config(cfg, stage):
    prefixes = {"Pointcloud": ("detector.", "peaks.", "frib.", "pad."), "Cluster": ("cluster.",), "Estimation": ("estimate.", "detector.", "cluster.smooth"), "InterpSolver": ("solver.", "gas.", "pid.", "detector.")}[stage]
    return {k: v for k, v in cfg.items() if k.startswith(prefixes) and not k.endswith(("mesh_cache", "pad_dir"))}

def suffix(stage, fmt): return ".root" if fmt == "root" else ".h5" if stage in STAGES[:2] else ".parquet"

def base_name(run, stage, cfg):
    base = f"run_{run:04d}"
    if stage == "InterpSolver":
        z, a = int(cfg.get("solver.Z", 1)), int(cfg.get("solver.A", 1))
        symbols = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne"}
        if z not in symbols: raise ValueError("请扩展 isotope symbol 映射以生成此核素的文件名")
        base += f"_{a}{symbols[z]}"
    return base

def stage_path(workspace, run, stage, cfg, fmt): return workspace / stage / (base_name(run, stage, cfg) + suffix(stage, fmt))

def existing(workspace, run, stage, cfg, fmt, manifest):
    candidates = [stage_path(workspace, run, stage, cfg, fmt), stage_path(workspace, run, stage, cfg, "python" if fmt == "root" else "root")]
    record = manifest["stages"].get(stage)
    if record:
        variants = record.get("outputs", [record.get("output", {})])
        for path in candidates:
            if path.is_file() and file_id(path) in variants: return path
        raise RuntimeError(f"{stage} 的已有文件与当前来源记录不一致；指定 --input 采用该文件，或从上游重算")
    for path in candidates:
        if path.is_file(): return path
    raise FileNotFoundError(f"缺少 {stage} 输入：{candidates[0]}（也可为 ROOT）")

def atomic_json(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)

def run_pipeline(args, scratch):
    config_path = Path(args.config).expanduser().resolve()
    cfg = load_config(config_path, scratch)
    workspace = Path(args.workspace).expanduser().resolve()
    for stage in STAGES: (workspace / stage).mkdir(parents=True, exist_ok=True)
    hidden = workspace / ".spyralpp"; hidden.mkdir(exist_ok=True)
    lock = open(hidden / f"run_{args.run:04d}.lock", "a")
    try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError: raise RuntimeError("此 workspace 的同一 run 正在运行")
    cfg["run.number"] = args.run
    cfg["run.threads"] = args.threads
    cfg.setdefault("cluster.use_hdbscan", "true")
    cfg.setdefault("peaks.dither_seed", 20261026)
    if args.waveform_fit is not None:
        cfg["peaks.template_fit"] = args.waveform_fit
    elif "pointcloud.waveform_fit" in cfg:
        cfg["peaks.template_fit"] = truth(cfg["pointcloud.waveform_fit"])
    if args.fitter: cfg["solver.fitter"] = args.fitter
    if args.python_mesh: cfg["solver.python_mesh"] = str(Path(args.python_mesh).expanduser().resolve())
    for item in args.set:
        if "=" not in item: raise ValueError("--set 需要 key=value")
        key, value = item.split("=", 1); cfg[key.lower()] = value
    cfg["solver.mesh_cache"] = str(hidden / "track_mesh.bin")
    if args.threads < 1 or (args.events is not None and args.events < 1): raise ValueError("threads/events 必须为正整数")
    selected_stages = args.enabled_stages
    if not selected_stages: raise ValueError("pipeline 至少需要启用一个阶段")
    first_index = min(STAGES.index(stage) for stage in selected_stages)
    resolved_cfg = scratch / "resolved.cfg"; write_config(resolved_cfg, cfg)
    stamp_path = hidden / f"run_{args.run:04d}.json"
    manifest = json.loads(stamp_path.read_text()) if stamp_path.exists() else {"run": args.run, "stages": {}}
    binaries = binary_id()
    invocation = {"started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "config": str(config_path), "from": selected_stages[0], "to": selected_stages[-1], "threads": args.threads, "events_limit": args.events, "format": args.format, "stages": {}, "binaries": binaries}
    invocation_start = time.perf_counter()
    for stage in selected_stages:
        i = STAGES.index(stage)
        inputs = []
        if stage == "Pointcloud":
            raw = raw_path_for_run(args.run, args, cfg, Path(args.config))
            if not raw.is_file(): raise FileNotFoundError(f"找不到原始数据 {raw}；使用 --trace-dir 或 --input")
            inputs.append(raw)
        else:
            inputs.append(Path(args.input).expanduser().resolve() if i == first_index and args.input else existing(workspace, args.run, STAGES[i-1], cfg, args.format, manifest))
            if stage == "InterpSolver": inputs.insert(0, Path(args.clusters).expanduser().resolve() if args.clusters else existing(workspace, args.run, "Cluster", cfg, args.format, manifest))
        output = stage_path(workspace, args.run, stage, cfg, args.format)
        phase_settings = phase_config(cfg, stage)
        # Custom map contents matter; temporary staging directory names do not.
        if stage == "Pointcloud":
            pad_dir = Path(cfg.get("detector.pad_dir", PROJECT / "data"))
            phase_settings["pad_checksums"] = {name: hashlib.sha256((pad_dir / name).read_bytes()).hexdigest() for name in ("padxy.csv", "pad_scale.csv", "pad_time_correction.csv")}
        signature = digest({"schema": 1, "settings": phase_settings, "inputs": [input_id(p, manifest) for p in inputs], "events_limit": args.events if i == 0 else None, "binaries": binaries})
        previous = manifest["stages"].get(stage, {})
        if not args.force and previous.get("signature") == signature and previous.get("output") and Path(previous["output"]["path"]).exists() and file_id(Path(previous["output"]["path"])) == previous["output"]:
            prior = Path(previous["output"]["path"])
            if prior == output:
                say(f"{stage}：复用 {output}")
                invocation["stages"][stage] = {"reused": True, "wall_seconds": 0., "output": str(output)}
                continue
            say(f"{stage}：已有结果转为 {args.format}，无需重算")
            do_convert = prior
        else: do_convert = None
        if output.exists() and stage not in manifest["stages"] and not args.force:
            raise RuntimeError(f"{output} 是未记录来源的已有文件。接续请用 --from 下一阶段；确定重算覆盖可加 --force")
        log = hidden / f"run_{args.run:04d}_{stage}.log"
        log.write_text(f"inputs: {[str(p) for p in inputs]}\n")
        stage_start = time.perf_counter()
        action = "转换已有结果" if do_convert else "开始处理"
        say(f"Run {args.run:04d} · {stage}：{action}，输出 {output}")
        with tempfile.TemporaryDirectory(prefix=f"{stage}_", dir=hidden) as task:
            task = Path(task)
            temporary = task / output.name
            if do_convert:
                convert(stage, do_convert, temporary, task, log)
            elif stage == "Pointcloud":
                command = [str(BIN / "rk4_waveform"), str(resolved_cfg), str(inputs[0]), str(temporary), str(args.events if args.events else 2147483647), str(args.threads)]
                status = run_logged(command, log)
                if status: raise RuntimeError(f"Pointcloud 失败：{log}\n{log.read_text()[-3500:]}")
            else:
                native_inputs = [to_native(path, "Estimation" if stage == "InterpSolver" and j == 1 else STAGES[i-1], task) for j, path in enumerate(inputs)]
                native_output = task / "output.csv" if temporary.suffix == ".parquet" else temporary
                command = {"Cluster": "cluster", "Estimation": "estimate", "InterpSolver": "solve"}[stage]
                native(command, *native_inputs, native_output, resolved_cfg, log=log)
                if temporary.suffix == ".parquet": csv_to_parquet(native_output, temporary, stage)
            if stage == "Pointcloud":
                peaks = Path(str(temporary) + "_peaks.root")
                if peaks.exists(): peaks.replace(hidden / f"run_{args.run:04d}_peaks.root")
            temporary.replace(output)
        elapsed = time.perf_counter() - stage_start
        variants = [entry for entry in previous.get("outputs", [previous.get("output", {})]) if entry and Path(entry["path"]).exists() and file_id(Path(entry["path"])) == entry] if previous.get("signature") == signature else []
        output_id = file_id(output)
        variants = [entry for entry in variants if entry["path"] != output_id["path"]] + [output_id]
        record = {"signature": signature, "output": output_id, "outputs": variants, "inputs": [file_id(p) for p in inputs], "settings": phase_settings, "wall_seconds": elapsed, "converted": bool(do_convert)}
        manifest["stages"][stage] = record
        atomic_json(stamp_path, manifest)
        invocation["stages"][stage] = {"reused": False, "wall_seconds": elapsed, "output": str(output), "converted": bool(do_convert)}
        say(f"{stage}：{elapsed:.3f} s → {output}")
    invocation["wall_seconds"] = time.perf_counter() - invocation_start
    atomic_json(hidden / f"run_{args.run:04d}_last_timing.json", invocation)
    # Keep effective configuration reviewable after temporary data are cleaned.
    write_config(hidden / f"run_{args.run:04d}_config.cfg", {k: v for k, v in cfg.items() if k != "detector.pad_dir" or not str(v).startswith(str(scratch))})
    say(f"完成阶段 {', '.join(selected_stages)}，本次 {invocation['wall_seconds']:.3f} s；时间记录：{hidden / f'run_{args.run:04d}_last_timing.json'}")

def parser():
    p = argparse.ArgumentParser(
        prog=os.environ.get("SPYRALPP_PROG", "spyralpp"),
        description="SpyralPP：一个 workspace，四个可独立接续的阶段。",
        epilog="示例：./spyral configs/local.cfg；./spyral 151 --config /path/to/main.py",
    )
    p.add_argument("target", nargs="?", help="run 号或配置文件；配置可包含 run.number 或 run.start/run.end")
    p.add_argument("--run", dest="run_override", type=int, help="强制只处理一个 run，覆盖配置区间")
    p.add_argument("--config", help="配置文件（原 main.py 或 .cfg/.yaml）；旧式 run 参数下仍可使用")
    p.add_argument("--workspace", help="覆盖配置里的 output.workspace（默认 ./workspace）")
    p.add_argument("--from", dest="start", choices=STAGES, help="选择连续阶段的起点，覆盖配置")
    p.add_argument("--to", choices=STAGES, help="选择连续阶段的终点，覆盖配置")
    p.add_argument("--format", choices=("python", "root"), help="覆盖 output.format：Python HDF5/Parquet 或 ROOT")
    p.add_argument("--input", help="所选起始阶段的输入文件；自动识别格式")
    p.add_argument("--clusters", help="InterpSolver 的 Cluster 文件（不在当前 workspace 时指定）")
    p.add_argument("--trace-dir", help="原始 run H5 所在目录")
    p.add_argument("--threads", type=int, help="覆盖 run.threads（默认 1）")
    p.add_argument("--events", type=int, help="只处理前 N 个原始事件；省略则全 run")
    waveform = p.add_mutually_exclusive_group()
    waveform.add_argument("--waveform-fit", dest="waveform_fit", action="store_true", help="Pointcloud 使用 GET 响应模板拟合")
    waveform.add_argument("--no-waveform-fit", dest="waveform_fit", action="store_false", help="Pointcloud 使用快速基线峰值法")
    p.set_defaults(waveform_fit=None)
    p.add_argument("--fitter", choices=("minuit", "python-objective", "kalman"), help="InterpSolver 后端；python-objective 复用 Python 目标函数但用 Minuit2 优化")
    p.add_argument("--python-mesh", help="python-objective 使用的官方 Spyral .npy 轨迹网格")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="覆盖个别配置项，可重复")
    p.add_argument("--force", action="store_true", help="重算所选阶段，可覆盖其旧文件")
    return p

def raw_path_for_run(run, args, cfg, config_path):
    if args.input: return Path(args.input).expanduser().resolve()
    if args.trace_dir:
        trace_dir = Path(args.trace_dir).expanduser()
    else:
        trace_dir_value = cfg.get("input.trace_dir")
        if not trace_dir_value:
            raise ValueError("请在配置的 [input] trace_dir 中设置原始 run 数据目录，或使用 --trace-dir")
        trace_dir = Path(trace_dir_value).expanduser()
        if not trace_dir.is_absolute(): trace_dir = config_path.resolve().parent / trace_dir
    if not trace_dir.exists() and str(trace_dir).startswith("/home/DATA/"):
        mapped = Path(str(trace_dir)[5:])
        if mapped.exists(): trace_dir = mapped
    return trace_dir / f"run_{run:04d}.h5"

def main():
    if not ENGINE.is_file(): raise RuntimeError(f"缺少 {ENGINE}，先运行 cmake --build build -j8")
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        p = argparse.ArgumentParser(description="阶段文件格式互换，无需重算")
        p.add_argument("stage", choices=STAGES); p.add_argument("input", type=Path); p.add_argument("output", type=Path)
        p.add_argument("--force", action="store_true")
        args = p.parse_args(sys.argv[2:]); args.output.parent.mkdir(parents=True, exist_ok=True)
        allowed = (".root", ".h5", ".hdf5") if args.stage in STAGES[:2] else (".root", ".parquet")
        if args.input.suffix not in allowed or args.output.suffix not in allowed: raise ValueError(f"{args.stage} 支持 {allowed}")
        if args.input.suffix == ".parquet" or args.output.suffix == ".parquet": require_polars()
        if args.output.exists() and not args.force: raise ValueError("输出已存在；需要覆盖时加 --force")
        if args.input.resolve() == args.output.resolve(): raise ValueError("转换输入与输出必须不同")
        with tempfile.TemporaryDirectory(prefix=".spyralpp-export-", dir=args.output.parent) as temp:
            scratch = Path(temp); target = scratch / args.output.name
            convert(args.stage, args.input, target, scratch, scratch / "export.log")
            target.replace(args.output)
        say(f"{args.stage} 已转换：{args.output}"); return
    cli = parser()
    args = cli.parse_args()
    default_config = PROJECT / "configs/local.cfg"
    explicit_number = args.target is not None and re.fullmatch(r"\d+", args.target)
    if explicit_number:
        config_path = Path(args.config).expanduser() if args.config else default_config
        if args.run_override is not None and args.run_override != int(args.target):
            cli.error("位置 run 号与 --run 不一致")
        runs = [int(args.target)]
    else:
        if args.target is not None and args.config:
            cli.error("配置文件已作为位置参数给出，不要再重复指定 --config")
        config_path = Path(args.target if args.target is not None else args.config or default_config).expanduser()
        config_path = config_path.resolve()
        config_values = read_config_values(config_path)
        if args.run_override is not None:
            runs = [args.run_override]
        else:
            has_start = "run.start" in config_values
            has_end = "run.end" in config_values
            if has_start != has_end:
                cli.error("run.start 和 run.end 必须同时设置")
            if has_start:
                try:
                    first, last = int(config_values["run.start"]), int(config_values["run.end"])
                except (TypeError, ValueError):
                    cli.error("run.start 和 run.end 必须是整数")
                if first < 0 or last < first:
                    cli.error("run 区间必须满足 0 <= start <= end")
                runs = list(range(first, last + 1))
            else:
                try: runs = [infer_run_number(config_path, config_values)]
                except ValueError as exc: cli.error(str(exc))
        if len(runs) > 1 and args.input:
            cli.error("run 区间不能共用单个 --input 文件；请让每个 run 从配置的数据目录读取")
    if not explicit_number:
        config_values = read_config_values(config_path)
    else:
        config_path = config_path.resolve()
        config_values = read_config_values(config_path)
    if args.format is None:
        args.format = str(config_values.get("output.format", "python")).strip().lower()
        if args.format not in ("python", "root"):
            cli.error("[output] format 必须是 python 或 root")
    if args.threads is None:
        try:
            args.threads = int(config_values.get("run.threads", 1))
        except (TypeError, ValueError):
            cli.error("[run] threads 必须是正整数")
    if args.threads < 1:
        cli.error("[run] threads 必须是正整数")
    stage_names = {stage.lower(): stage for stage in STAGES}
    stage_flags = {
        "Pointcloud": "pipeline.pointcloud",
        "Cluster": "pipeline.cluster",
        "Estimation": "pipeline.estimation",
        "InterpSolver": "pipeline.interp_solver",
    }
    has_stage_flags = any(key in config_values for key in stage_flags.values())
    if args.start is not None or args.to is not None:
        first = args.start or "Pointcloud"
        last = args.to or "Estimation"
        if STAGES.index(first) > STAGES.index(last):
            cli.error("--from 不能晚于 --to")
        args.enabled_stages = STAGES[STAGES.index(first):STAGES.index(last) + 1]
    elif has_stage_flags:
        args.enabled_stages = [
            stage for stage, key in stage_flags.items()
            if truth(config_values.get(key, False))
        ]
        if not args.enabled_stages:
            cli.error("[pipeline] 至少要把一个阶段设为 true")
    else:
        # Backward compatibility for older configs with pipeline.from / pipeline.to.
        raw_first = str(config_values.get("pipeline.from", "Pointcloud")).strip().lower()
        raw_last = str(config_values.get("pipeline.to", "Estimation")).strip().lower()
        first, last = stage_names.get(raw_first), stage_names.get(raw_last)
        if first is None: cli.error(f"pipeline.from 阶段无效：{config_values.get('pipeline.from')}")
        if last is None: cli.error(f"pipeline.to 阶段无效：{config_values.get('pipeline.to')}")
        if STAGES.index(first) > STAGES.index(last):
            cli.error("pipeline.from 不能晚于 pipeline.to")
        args.enabled_stages = STAGES[STAGES.index(first):STAGES.index(last) + 1]
    args.start, args.to = args.enabled_stages[0], args.enabled_stages[-1]
    if any(run < 0 for run in runs): cli.error("run 号不能为负")
    parquet_input = any(Path(value).suffix == ".parquet" for value in (args.input, args.clusters) if value)
    parquet_output = args.format == "python" and any(
        stage in args.enabled_stages for stage in ("Estimation", "InterpSolver")
    )
    if parquet_input or parquet_output:
        require_polars()
    args.config = str(config_path)
    if args.workspace is None:
        args.workspace = str(config_values.get("output.workspace", "workspace"))
    skip_missing = truth(config_values.get("run.skip_missing", True))
    workspace = Path(args.workspace).expanduser().resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    for stage in STAGES:
        (workspace / stage).mkdir(parents=True, exist_ok=True)
    processed = 0
    skipped = 0
    for run in runs:
        args.run = run
        if "Pointcloud" in args.enabled_stages and not args.input:
            raw = raw_path_for_run(run, args, config_values, config_path)
            if not raw.is_file():
                if skip_missing:
                    say(f"跳过 Run {run:04d}：原始文件不存在 ({raw})")
                    skipped += 1
                    continue
                raise FileNotFoundError(f"找不到原始数据 {raw}；使用 --trace-dir 或 --input")
        hidden = workspace / ".spyralpp"
        hidden.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f"config_{run:04d}_", dir=hidden) as temp:
            run_pipeline(args, Path(temp))
        processed += 1
    if len(runs) > 1 or skipped:
        say(f"Run 区间完成：处理 {processed} 个，跳过 {skipped} 个不存在的 run。")

if __name__ == "__main__":
    try: main()
    except KeyboardInterrupt:
        print("\nSpyralPP：收到中断，当前阶段已停止。", file=sys.stderr)
        sys.exit(130)
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print(f"SpyralPP：{exc}", file=sys.stderr); sys.exit(1)
