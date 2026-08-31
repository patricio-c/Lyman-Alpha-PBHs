#!/usr/bin/env python3
"""
Stage 04 - the 1D flux power spectrum and the ratio between runs.

This is the figure the whole analysis exists to produce.  Any number of
caches can be given; the first one is the denominator of the ratio panel.

The normalisation is a command-line switch, not a hard-coded choice,
because the first question anyone asks is whether the normalisation made
the result:

    --norm taueff   rescale every run to a common tau_eff (default)
    --norm none     A = 1 everywhere, use the field as the simulation made it
    --norm own      each run rescaled to its own Turner+24 value at its own z
                    (only differs from `none` for runs at different z)

Usage
-----
    python stages/04_p1d.py cache/cache_cdm.npz cache/cache_fct.npz \\
        --labels CDM FCT --out figures/p1d

    python stages/04_p1d.py cache/*.npz --norm none --out figures/p1d_raw

Options
-------
    caches...          one or more .npz; the first is the reference
    --labels ...       display names, same order (default: file stems)
    --norm MODE        taueff | none | own
    --tau-eff X        override the common target
    --kmax-frac F      cut at F x Nyquist (default 0.5)
    --nbins N          log-k bins (default 28)
    --sherwood PATH    overlay a Sherwood cache as an external reference
    --ylim-ratio A B   y range of the ratio panel (default 0.6 1.4)
    --out PREFIX       writes PREFIX.png and PREFIX.txt
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import cache as cachelib, units  # noqa: E402
from common.p1d import (logbin, nyquist_cut, p1d_from_tau,  # noqa: E402
                        tau_eff)

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd",
           "#ff7f0e", "#8c564b", "#17becf"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("caches", nargs="+")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--norm", default="taueff",
                    choices=["taueff", "none", "own"])
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--kmax-frac", type=float, default=0.5)
    ap.add_argument("--nbins", type=int, default=28)
    ap.add_argument("--sherwood", default=None)
    ap.add_argument("--ylim-ratio", type=float, nargs=2, default=[0.6, 1.4])
    ap.add_argument("--out", default="figures/p1d")
    args = ap.parse_args()

    labels = args.labels or [os.path.basename(p).replace("cache_", "")
                             .replace(".npz", "") for p in args.caches]
    if len(labels) != len(args.caches):
        raise SystemExit("--labels must have one entry per cache")

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    runs_ = []
    for p, lab in zip(args.caches, labels):
        runs_.append((lab, cachelib.load(p)))

    z0 = runs_[0][1].z
    dv0 = runs_[0][1].dv
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(z0)

    say("=" * 70)
    say(f"Stage 04 - P1D    norm = {args.norm}")
    say("=" * 70)
    say(f"reference run: {labels[0]}    z = {z0:.4f}   dv = {dv0:.5f} km/s")
    if args.norm == "taueff":
        say(f"common tau_eff target = {target:.5f}")
    say()
    say(f"{'run':14s} {'nlos':>6s} {'z':>7s} {'tau_eff raw':>12s} "
        f"{'A':>9s} {'tau_eff out':>12s}")

    curves = {}
    for lab, c in runs_:
        if args.norm == "taueff":
            k, P, e, A = p1d_from_tau(c.tau, c.dv, target=target)
        elif args.norm == "own":
            k, P, e, A = p1d_from_tau(
                c.tau, c.dv, target=units.tau_eff_turner24(c.z))
        else:
            k, P, e, A = p1d_from_tau(c.tau, c.dv, A=1.0)
        say(f"{lab:14s} {c.tau.shape[0]:6d} {c.z:7.3f} "
            f"{tau_eff(c.tau, 1.0):12.5f} {A:9.5f} "
            f"{tau_eff(c.tau, A):12.5f}")
        curves[lab] = (k, P, e, c.dv, c.z)

    if args.sherwood:
        cs = cachelib.load(args.sherwood)
        k, P, e, A = p1d_from_tau(cs.tau, cs.dv, target=target)
        curves["Sherwood"] = (k, P, e, cs.dv, cs.z)
        labels = labels + ["Sherwood"]
        say(f"{'Sherwood':14s} {cs.tau.shape[0]:6d} {cs.z:7.3f} "
            f"{tau_eff(cs.tau, 1.0):12.5f} {A:9.5f} "
            f"{tau_eff(cs.tau, A):12.5f}")
    say()

    kmax = args.kmax_frac * np.pi / dv0
    binned = {}
    for lab in labels:
        k, P, e, dv, z = curves[lab]
        m = nyquist_cut(k, dv, args.kmax_frac)
        kb, pb, eb = logbin(k[m], P[m], e[m], nbins=args.nbins)
        binned[lab] = (kb, pb, eb)

    ref = labels[0]
    kref = binned[ref][0]

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 8.2), sharex=True,
                                 gridspec_kw={"height_ratios": [2.1, 1],
                                              "hspace": 0.06})
    for i, lab in enumerate(labels):
        kb, pb, eb = binned[lab]
        col = PALETTE[i % len(PALETTE)]
        ax.plot(kb, kb * pb / np.pi, color=col, lw=2.0, label=lab)
        ax.fill_between(kb, kb * (pb - eb) / np.pi, kb * (pb + eb) / np.pi,
                        color=col, alpha=0.18, lw=0)

    ax.axvspan(*units.desi_window(z0), color="green", alpha=0.07, zorder=0)
    ax.axvline(2 * np.pi / (dv0 * runs_[0][1].tau.shape[1]),
               color="0.4", ls=":", lw=1)
    ax.set(xscale="log", yscale="log", ylabel=r"$k\,P_{\rm 1D}(k)\,/\,\pi$",
           title=(rf"$z = {z0:.1f}$,  normalisation: {args.norm}"
                  + (rf",  $\tau_{{\rm eff}} = {target:.4f}$"
                     if args.norm == "taueff" else "")))
    ax.legend(frameon=False, loc="lower left", fontsize=9)
    ax.grid(alpha=0.2, which="both")

    say(f"{'k [s/km]':>10s} " + " ".join(f"{l:>10s}" for l in labels[1:]))
    anchors = [0.003, 0.005, 0.01, 0.02, 0.03, 0.06]
    rows = {lab: [] for lab in labels[1:]}
    for i, lab in enumerate(labels[1:], start=1):
        kb, pb, eb = binned[lab]
        pr = np.interp(np.log(kref), np.log(kb), np.log(np.maximum(pb, 1e-300)))
        p0 = np.log(np.maximum(binned[ref][1], 1e-300))
        r = np.exp(pr - p0)
        col = PALETTE[i % len(PALETTE)]
        bx.plot(kref, r, color=col, lw=2.0, label=f"{lab} / {ref}")
        for ka in anchors:
            if kref.min() <= ka <= kref.max():
                rows[lab].append(float(np.interp(ka, kref, r)))
            else:
                rows[lab].append(np.nan)
    for j, ka in enumerate(anchors):
        vals = " ".join(f"{rows[l][j]:10.4f}" for l in labels[1:])
        say(f"{ka:10.4f} {vals}")

    bx.axhline(1.0, color="k", ls="--", lw=1)
    bx.axvspan(*units.desi_window(z0), color="green", alpha=0.07, zorder=0)
    bx.set(xscale="log", xlabel=r"$k$  [s km$^{-1}$]",
           ylabel=f"ratio to {ref}", ylim=tuple(args.ylim_ratio))
    bx.legend(frameon=False, fontsize=9)
    bx.grid(alpha=0.2, which="both")

    y0, y1 = ax.get_ylim()
    lo, hi = units.desi_window(z0)
    ax.text(np.sqrt(lo * hi), y1 * 0.45, "DESI window", color="green",
            ha="center", fontsize=9)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    say()
    say(f"kmax used = {kmax:.4f} s/km  ({args.kmax_frac} x Nyquist)")
    print(f"written -> {args.out}.png")
    print(f"written -> {args.out}.txt")


if __name__ == "__main__":
    main()
