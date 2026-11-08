#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SpyralPP 网页控制器（只用 Python 标准库，不用装任何东西）

    python3 scripts/webui.py            # 打开 http://127.0.0.1:8765

功能：
  · 教程标签页：从"怎么跑"到"参数含义"的速查
  · 参数表单：按分组列出全部可调参数，带默认值与说明，保存成本项目认的 .cfg
  · 运行按钮：直接调 build/rk4_waveform / rk4_analyze / rk4_merge，日志实时滚动
  · 产物列表：点云/拟合/中间文件的大小与时间

只监听 127.0.0.1，不做鉴权；要给别人用得自己加。
"""

import json
import os
import subprocess
import sys
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.environ.get("SPYRALPP_BIN", os.path.join(ROOT, "build"))
OUT_DIR = os.path.join(ROOT, "web_output")
CFG_PATH = os.path.join(OUT_DIR, "web.cfg")
HTML_PATH = os.path.join(HERE, "webui.html")
PORT = int(os.environ.get("SPYRALPP_WEB_PORT", "8765"))

# ---------------------------------------------------------------- 参数表
# type: int / float / bool / text / select
SCHEMA = [
    ("输入输出", [
        ("run", "input", "text", "", "原始波形 .h5（GET/merger）或文本波形；留空表示只用已有想点云"),
        ("run", "clouds", "text", "", "已有想点云 .root；给了就跳过波形转换"),
        ("run", "tag", "text", "web", "输出文件名前缀"),
        ("run", "n_events", "int", 200, "只处理前 N 个事件（0 = 全部）"),
        ("run", "threads", "int", 8, "线程数（环境变量 SPYRALPP_THREADS）"),
        ("run", "shards", "int", 1, "重建分片数（多进程）"),
    ]),
    ("探测器", [
        ("detector", "length_m", "float", 1.0, "灵敏体积长度（m）"),
        ("detector", "bfield_t", "float", 2.85, "磁场（T）。真实数据手性不对时改 -2.85"),
        ("detector", "efield_vm", "float", 55000.0, "电场（V/m）"),
        ("detector", "micromegas_edge", "float", 10.0, "micromegas 时间参考（TB）"),
        ("detector", "windows_edge", "float", 560.0, "靶窗时间参考（TB）"),
        ("detector", "adc_threshold", "float", 10.0, "点云振幅阈值"),
        ("detector", "num_tb", "int", 512, "时间桶数"),
        ("detector", "pad_dir", "text", "", "pad 表目录（留空 = 自动找 data/）"),
    ]),
    ("波形 / 峰参数", [
        ("peaks", "template_fit", "bool", False, "1 = 用 GET 电子学响应模板拟合（含饱和修复）；0 = 局部极大"),
        ("peaks", "sampling_ns", "float", 320.0, "Run 0151 采样周期（ns）；3.125 MHz = 320，按 run 覆盖"),
        ("peaks", "peaking_ns", "float", 1000.0, "成形/峰化时间（ns）"),
        ("peaks", "sat_adc", "float", 4095.0, "ADC 削顶电平（12-bit = 4095）"),
        ("peaks", "max_peaks", "int", 12, "单条迹最多拟几个峰"),
        ("peaks", "residual_threshold", "float", 40.0, "残差峰幅度阈值（ADC）"),
        ("peaks", "separation", "int", 20, "峰最小间距（TB）"),
        ("peaks", "prominence", "float", 50.0, "峰最小突出度"),
        ("peaks", "max_width", "float", 100.0, "峰最大宽度（TB）"),
        ("peaks", "threshold", "float", 40.0, "峰振幅阈值"),
        ("peaks", "baseline_scale", "float", 20.0, "FFT sinc 基线滤波尺度"),
        ("peaks", "dump", "bool", False, "1 = 写每条迹线的峰参数 <out>_peaks.root"),
    ]),
    ("聚类", [
        ("cluster", "eps_mm", "float", 20.0, "DBSCAN 邻域半径（mm）"),
        ("cluster", "min_points", "int", 3, "核心点所需邻居数"),
        ("cluster", "min_cluster_size", "int", 10, "最小簇点数"),
        ("cluster", "direction_fraction", "float", 0.25, "用首尾各多少比例的点判方向"),
        ("cluster", "enable_join", "bool", True, "把切碎的同一径迹接起来"),
        ("cluster", "join_gap_mm", "float", 60.0, "端点最大间距"),
        ("cluster", "join_angle_deg", "float", 65.0, "端点夹角上限"),
        ("cluster", "join_edge_points", "int", 10, "估计端点切向用的点数"),
        ("cluster", "adaptive_eps", "bool", False, "eps 随最近邻中位数自适应"),
        ("cluster", "eps_scale", "float", 2.0, "自适应倍率"),
        ("cluster", "eps_max_mm", "float", 40.0, "自适应上限"),
        ("cluster", "smooth", "bool", True, "拟合前平滑样条降噪"),
        ("cluster", "smooth_lambda", "float", 1.0, "平滑强度"),
        ("cluster", "report", "bool", False, "打印每事件的 NN/eps/簇数/后端"),
        ("cluster", "use_hdbscan", "bool", True, "换 HDBSCAN 后端"),
        ("cluster", "hdbscan_min_samples", "int", 3, "HDBSCAN 核心距离 k"),
        ("cluster", "hdbscan_max_points", "int", 0, "HDBSCAN 点数限制；0 不限制，超过时显式失败"),
        ("cluster", "hdbscan_allow_single", "bool", False, "允许整云单簇"),
    ]),
    ("求解 / 拟合", [
        ("solver", "A", "int", 1, "被拟合粒子质量数"),
        ("solver", "Z", "int", 1, "被拟合粒子电荷数"),
        ("solver", "ke_bins", "int", 256, "网格动能格点数（高分辨率默认）"),
        ("solver", "ke_min", "float", 0.1, "网格最小动能（MeV）"),
        ("solver", "ke_max", "float", 70.0, "网格最大动能（MeV）"),
        ("solver", "polar_bins", "int", 129, "网格极角格点数（高分辨率默认）"),
        ("solver", "polar_min", "float", 2.0, "网格极角下限（度）"),
        ("solver", "polar_max", "float", 88.0, "网格极角上限（度）"),
        ("solver", "n_time", "int", 400, "网格时间点数"),
        ("solver", "dt", "float", 3e-10, "轨迹网格时间步长（秒；0.3 ns 默认，400 点约覆盖 120 ns）"),
        ("solver", "mesh_cache", "text", "/tmp/proton_h2.mesh", "网格缓存路径"),
        ("solver", "mesh_shared", "bool", True, "mmap 共享网格"),
        ("solver", "fitter", "select", "minuit", "拟合后端", ["minuit", "kalman"]),
        ("solver", "minuit_sort_candidates", "bool", True, "Minuit 最近轨迹段候选按距离排序（推荐默认）"),
        ("solver", "minuit_centered_azimuth_bounds", "bool", True, "每个 Minuit 起点使用居中的完整周期方位角区间"),
        ("solver", "kalman_meas_sigma", "float", 0.5, "Kalman 测点 σ（mm）"),
        ("solver", "kalman_mom_sigma_frac", "float", 0.1, "动量先验相对不确定度"),
        ("solver", "kalman_min_clusters", "int", 5, "Kalman 最少量测点数"),
        ("solver", "kalman_smooth", "bool", True, "RTS 平滑"),
        ("solver", "kalman_straggling", "bool", False, "能损 straggling 过程噪声"),
        ("solver", "kalman_mult_scattering", "bool", False, "多重散射过程噪声"),
        ("solver", "kalman_chi2_max", "float", 0.0, "χ²/ndf 上限（0 = 不限）"),
    ]),
    ("PID 预筛", [
        ("pid", "gate", "bool", False, "拟合前先用粗估 (Bρ, dE/dx) 判粒子"),
        ("pid", "min_sep", "float", 0.3, "对数似然间隔阈值，越大越严"),
        ("pid", "require", "int", -1, "只保留某候选（-1 不限，0=p 1=d 2=t）"),
    ]),
    ("输出", [
        ("output", "hdf5", "bool", False, "同时导出 HDF5"),
        ("output", "stages", "bool", False, "写中间过程文件 <fit>_stages.root"),
    ]),
]

JOBS = {}


def default_values():
    vals = {}
    for _, items in SCHEMA:
        for item in items:
            sec, key, kind, default = item[0], item[1], item[2], item[3]
            vals[sec + "." + key] = default
    return vals


def schema_json():
    groups = []
    for name, items in SCHEMA:
        fields = []
        for item in items:
            sec, key, kind, default = item[0], item[1], item[2], item[3]
            entry = {
                "key": sec + "." + key,
                "kind": kind,
                "default": default,
                "help": item[4] if len(item) > 4 else "",
            }
            if kind == "select":
                entry["options"] = item[5]
            fields.append(entry)
        groups.append({"group": name, "fields": fields})
    return groups


def write_cfg(values):
    os.makedirs(OUT_DIR, exist_ok=True)
    by_section = {}
    for k, v in values.items():
        sec, _, key = k.partition(".")
        if not key:
            continue
        by_section.setdefault(sec, []).append((key, v))
    with open(CFG_PATH, "w", encoding="utf-8") as f:
        for sec in by_section:
            f.write("[%s]\n" % sec)
            for key, v in by_section[sec]:
                if isinstance(v, bool):
                    v = "1" if v else "0"
                f.write("%s = %s\n" % (key, v))
            f.write("\n")
    return CFG_PATH


def out_path(name):
    return os.path.join(OUT_DIR, name)


def spawn(args, env_extra=None):
    env = os.environ.copy()
    env.update(env_extra or {})
    return subprocess.Popen(args, cwd=ROOT, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1, env=env)


def run_job(job_id, steps):
    job = JOBS[job_id]
    for args, env_extra in steps:
        job["log"].append("$ " + " ".join(args) + "\n")
        proc = spawn(args, env_extra)
        for line in proc.stdout:
            job["log"].append(line)
        proc.wait()
        job["log"].append("[退出码 %d]\n" % proc.returncode)
        if proc.returncode != 0:
            job["status"] = "failed"
            return
    job["status"] = "done"


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        data = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        url = urlparse(self.path)
        if url.path in ("/", "/index.html"):
            with open(HTML_PATH, "rb") as f:
                self._send(200, f.read(), "text/html; charset=utf-8")
            return
        if url.path == "/api/schema":
            self._send(200, json.dumps({"groups": schema_json(),
                                        "defaults": default_values(),
                                        "cfg_path": CFG_PATH,
                                        "out_dir": OUT_DIR,
                                        "bin": BIN}, ensure_ascii=False))
            return
        if url.path == "/api/log":
            q = parse_qs(url.query)
            jid = q.get("id", [""])[0]
            job = JOBS.get(jid)
            if not job:
                self._send(404, json.dumps({"error": "no such job"}))
                return
            self._send(200, json.dumps({"status": job["status"],
                                        "log": "".join(job["log"])}, ensure_ascii=False))
            return
        if url.path == "/api/files":
            items = []
            if os.path.isdir(OUT_DIR):
                for name in sorted(os.listdir(OUT_DIR)):
                    p = os.path.join(OUT_DIR, name)
                    if os.path.isfile(p):
                        items.append({"name": name, "size": os.path.getsize(p),
                                      "mtime": time.strftime("%H:%M:%S",
                                                             time.localtime(os.path.getmtime(p)))})
            self._send(200, json.dumps({"files": items}, ensure_ascii=False))
            return
        self._send(404, json.dumps({"error": "unknown"}))

    def do_POST(self):
        url = urlparse(self.path)
        length = int(self.headers.get("Content-Length", "0"))
        body = json.loads(self.rfile.read(length) or b"{}")

        if url.path == "/api/config":
            path = write_cfg(body.get("values", {}))
            self._send(200, json.dumps({"ok": True, "path": path}))
            return

        if url.path == "/api/run":
            action = body.get("action")
            vals = body.get("values", {})
            write_cfg(vals)
            threads = str(int(vals.get("run.threads", 8) or 8))
            n_events = int(vals.get("run.n_events", 0) or 0)
            tag = vals.get("run.tag", "web") or "web"
            input_path = vals.get("run.input", "") or ""
            clouds_path = vals.get("run.clouds", "") or ""
            cloud_out = out_path(tag + "_cloud.root")
            fit_out = out_path(tag + "_fit.root")
            os.makedirs(OUT_DIR, exist_ok=True)

            steps = []
            if action in ("clouds", "all"):
                if not input_path:
                    self._send(400, json.dumps({"error": "请先填 --input（原始波形文件）"}))
                    return
                args = [os.path.join(BIN, "rk4_waveform"), CFG_PATH, input_path, cloud_out]
                if n_events > 0:
                    args.append(str(n_events))
                args.append(threads)
                steps.append((args, {"SPYRALPP_THREADS": threads}))
                clouds_path = cloud_out

            if action in ("fit", "all"):
                if not clouds_path:
                    self._send(400, json.dumps({"error": "请先填点云路径，或先用『生成点云』"}))
                    return
                steps.append(([os.path.join(BIN, "rk4_analyze"), clouds_path, fit_out, CFG_PATH],
                              {"SPYRALPP_THREADS": threads}))

            if not steps:
                self._send(400, json.dumps({"error": "未知动作"}))
                return

            jid = uuid.uuid4().hex[:12]
            JOBS[jid] = {"status": "running", "log": []}
            threading.Thread(target=run_job, args=(jid, steps), daemon=True).start()
            self._send(200, json.dumps({"job": jid}))
            return

        if url.path == "/api/stop":
            jid = body.get("id", "")
            job = JOBS.get(jid)
            if job:
                job["status"] = "stopped"
            self._send(200, json.dumps({"ok": True}))
            return

        self._send(404, json.dumps({"error": "unknown"}))


def main():
    if not os.path.isfile(HTML_PATH):
        print("找不到 %s" % HTML_PATH)
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    server = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    print("SpyralPP 网页控制器：http://127.0.0.1:%d" % PORT)
    print("  可执行目录：%s" % BIN)
    print("  输出目录  ：%s" % OUT_DIR)
    print("  Ctrl-C 退出")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n再见")
    return 0


if __name__ == "__main__":
    sys.exit(main())
