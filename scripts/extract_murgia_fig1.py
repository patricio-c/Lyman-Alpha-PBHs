#!/usr/bin/env python3
"""
Recover the six curves of Murgia+2019 Figure 1, exactly.

The arXiv e-print (https://arxiv.org/e-print/1903.10509) ships the figure as a
vector PDF, `relative_spectra.pdf`.  `pdftocairo -svg` exposes the polyline
vertices matplotlib wrote, so this is not digitisation by eye: it is the same
numbers the authors plotted.  The calibration constants below are the
tick-mark paths of that same file.

Two checks the output has to pass, both printed at the end:

  1. the vertices land on exact multiples of k_fund = 2 pi / 20 = 0.3142 h/Mpc,
     the fundamental mode of their box, so these are their P1D bins and not a
     smoothed trace;
  2. the LINEAR M3/M2 ratio is 10.00 at every k, which is exactly
     P_iso proportional to M_PBH f_PBH (their Eqs. 1-2).

If either drifts, the calibration below no longer matches the file and nothing
downstream should be trusted.

Usage
-----
    python scripts/extract_murgia_fig1.py [--workdir DIR]

Writes `murgia2019_fig1.json` in the working directory: one entry per curve,
`k_hMpc` in comoving h/Mpc and `pct` as the percentage excess over LCDM.
Curve names are {linear,nonlinear,flux}_{M2,M3} with M2 = 10^2 M_sun and
M3 = 10^3 M_sun, both at f_PBH = 1, z = 5.

Needs `pdftocairo` (poppler-utils) on PATH.  Downloads ~0.9 MB from arXiv the
first time; Clementina has no outbound internet, so run this on the laptop.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import urllib.request

import numpy as np

EPRINT = "https://arxiv.org/e-print/1903.10509"

# Calibration, read off the black tick paths inside relative_spectra.pdf.
# Do not adjust these by hand; if they stop working, re-read the tick paths.
X0, K0 = 147.684162, 1.0        # x of the 10^0 h/Mpc major tick
X1, K1 = 331.247154, 10.0       # x of the 10^1 h/Mpc major tick
Y0, V0 = 58.725157, 0.0         # y of the 0 per cent tick
Y1, V1 = 248.985527, 140.0      # y of the 140 per cent tick

# Line style encodes the field, colour encodes the PBH mass.
STYLE = {"5.55,2.4": "linear", "1.5,2.475": "nonlinear", "-": "flux"}
COLOR = {"0%,0%,100%": "M2", "100%,0%,0%": "M3"}

K_FUND = 2.0 * np.pi / 20.0     # their box is 20 Mpc/h


def fetch(workdir):
    """Download and unpack the e-print, then rasterise Figure 1 to SVG."""
    pdf = os.path.join(workdir, "relative_spectra.pdf")
    if not os.path.exists(pdf):
        tgz = os.path.join(workdir, "src.tar.gz")
        if not os.path.exists(tgz):
            print(f"downloading {EPRINT}")
            urllib.request.urlretrieve(EPRINT, tgz)
        subprocess.run(["tar", "xzf", tgz], cwd=workdir, check=True)
    svg = os.path.join(workdir, "fig1.svg")
    if not os.path.exists(svg):
        subprocess.run(["pdftocairo", "-svg", "relative_spectra.pdf", "fig1.svg"],
                       cwd=workdir, check=True)
    return svg


def read_curves(svg):
    text = open(svg).read()

    def to_k(x):
        f = (np.log10(K1) - np.log10(K0)) / (X1 - X0)
        return 10.0 ** (np.log10(K0) + (np.asarray(x) - X0) * f)

    def to_v(y):
        return V0 + (np.asarray(y) - Y0) * (V1 - V0) / (Y1 - Y0)

    curves = {}
    for path in re.finditer(r"<path([^>]*)/>", text, re.S):
        attrs = path.group(1)
        colour = re.search(r"stroke:rgb\(([^)]*)\)", attrs)
        data = re.search(r'\sd="([^"]*)"', attrs)
        if not (colour and data) or colour.group(1) not in COLOR:
            continue
        dash = re.search(r"stroke-dasharray:([^;\"]*)", attrs)
        style = STYLE[dash.group(1).strip() if dash else "-"]
        pts = np.array([[float(u), float(v)] for u, v in
                        re.findall(r"([-\d.]+)\s+([-\d.]+)", data.group(1))])
        if len(pts) < 3:
            continue                      # the swatches in the legend box
        k, v = to_k(pts[:, 0]), to_v(pts[:, 1])
        order = np.argsort(k)
        curves[f"{style}_{COLOR[colour.group(1)]}"] = (k[order], v[order])
    return curves


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default=".")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    curves = read_curves(fetch(args.workdir))
    missing = {f"{s}_{m}" for s in ("linear", "nonlinear", "flux")
               for m in ("M2", "M3")} - set(curves)
    if missing:
        raise SystemExit(f"did not recover {sorted(missing)}; the calibration "
                         "or the file layout changed, do not use the output")

    out = os.path.join(args.workdir, "murgia2019_fig1.json")
    json.dump({n: {"k_hMpc": list(map(float, k)), "pct": list(map(float, v))}
               for n, (k, v) in curves.items()}, open(out, "w"), indent=1)

    print(f"{'curve':16s} {'npts':>5s} {'k range [h/Mpc]':>24s} {'excess [%]':>22s}")
    for name in sorted(curves):
        k, v = curves[name]
        print(f"{name:16s} {len(k):5d}   {k.min():8.4f} - {k.max():8.4f}   "
              f"{v.min():9.2f} - {v.max():9.2f}")

    k = curves["flux_M3"][0]
    resid = np.abs(k / K_FUND - np.round(k / K_FUND)).max()
    print(f"\ncheck 1  vertices on multiples of k_fund={K_FUND:.4f}: "
          f"worst residual {resid:.4f}  ({'ok' if resid < 0.02 else 'FAIL'})")

    def at(name, q):
        kk, vv = curves[name]
        return float(np.interp(np.log10(q), np.log10(kk), vv))

    ratios = [at("linear_M3", q) / at("linear_M2", q) for q in (5, 8.22, 10, 15, 19)]
    worst = max(abs(r - 10.0) for r in ratios)
    print("check 2  linear M3/M2 (must be 10.00): "
          + " ".join(f"{r:.2f}" for r in ratios)
          + f"  ({'ok' if worst < 0.05 else 'FAIL'})")
    print(f"\nwritten -> {out}")


if __name__ == "__main__":
    main()
