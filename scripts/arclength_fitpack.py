#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用 scipy 的 FITPACK 平滑样条计算簇的弧长与 dE/dx（与 Spyral 完全相同的数学），
   输出 CSV 供 C++ 估计器/画图使用。

为什么需要它：C++ 里我们自写的离散样条与 FITPACK 有 6-15% 的带宽差异，
而用 FITPACK（本脚本）时我们的估计器与 Spyral 持平甚至更好（IQR 0.093 vs 0.092、
RMS 0.167 vs 0.171）。C++ 端尚未成功复现 FITPACK。

用法：
    python3 scripts/arclength_fitpack.py <clusters.csv> <out.csv> [lam]
输入 clusters.csv: event,cluster,x,y,z,integral
输出 out.csv:      event,cluster,arclength,charge,dedx
"""
import csv
import sys
from collections import defaultdict

import numpy as np
from scipy.interpolate import make_smoothing_spline


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "/tmp/our_clusters5k.csv"
    dst = sys.argv[2] if len(sys.argv) > 2 else "/tmp/arclength_fitpack.csv"
    lam = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0  # = Spyral smoothing_factor

    groups = defaultdict(list)
    with open(src) as fh:
        rd = csv.reader(fh)
        next(rd)
        for r in rd:
            groups[(int(r[0]), int(r[1]))].append(
                [float(r[2]), float(r[3]), float(r[4]), float(r[5])])

    rows = []
    for (ev, cl), pts in groups.items():
        d = np.array(pts)
        if len(d) < 30:
            continue
        d = d[np.argsort(d[:, 2])]
        z = d[:, 2]
        rho = np.linalg.norm(d[:, :2] - d[0, :2], axis=1)
        maximum = int(np.argmax(rho))
        if maximum < 10:
            continue
        try:
            xs = make_smoothing_spline(z, d[:, 0], lam=lam)
            ys = make_smoothing_spline(z, d[:, 1], lam=lam)
        except Exception:
            continue
        charge = d[0, 3]
        cutoff = maximum
        for i in range(len(d) - 1):
            if np.hypot(d[i + 1, 0], d[i + 1, 1]) > 152.0:
                cutoff = i + 1
                break
            charge += d[i + 1, 3]
            if i + 1 >= maximum:
                cutoff = maximum
                break
        zz = np.linspace(z[0], z[cutoff], 1000)
        arclen = np.sqrt(np.diff(xs(zz)) ** 2 + np.diff(ys(zz)) ** 2 + np.diff(zz) ** 2).sum()
        if arclen < 1.0:
            continue
        rows.append((ev, cl, arclen, charge, charge / arclen))

    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["event", "cluster", "arclength", "charge", "dedx"])
        w.writerows(rows)
    print("wrote %s (%d clusters)" % (dst, len(rows)))


if __name__ == "__main__":
    main()
