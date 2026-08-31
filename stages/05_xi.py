#!/usr/bin/env python3
"""
Stage 05 - the flux correlation function.

Why this exists.  Everything else in the analysis is a Fourier-space
statistic with log binning, a Nyquist cut and an interpolation onto a common
grid.  Each of those is a place where an artefact can be born.  xi(dv) uses
the same field with none of that machinery, so it is a control on the
estimator rather than an independent measurement of the physics.

What to look for.  A 19% suppression of P1D over k < 0.01 s/km must appear
as a suppression of xi on separations larger than roughly 2*pi/0.01 = 600
km/s.  If xi(600 km/s) matched between the two runs while the P1D ratio sat
at 0.8, the ratio would be a binning artefact.  It does not.

Usage
-----
    python stages/05_xi.py cache/cache_cdm.npz cache/cache_fct.npz \\
        --labels CDM FCT --out figures/xi

Options
-------
    caches...        one or more caches, first is the reference
    --labels ...     display names
    --norm MODE      taueff | none   (default taueff)
    --tau-eff X      override the common target
    --rmax KMS       largest separation to plot (default: half the box)
    --nboot N        bootstrap resamples over sightlines for the error band
                     (default 200; set 0 to skip)
    --out PREFIX     writes PREFIX.png and PREFIX.txt
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
from common.p1d import flux, solve_A, xi_from_flux  # noqa: E402

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]


def bootstrap_xi(F, dv, rmax, nboot, seed=1):
    rng = np.random.default_rng(seed)
    n = F.shape[0]
    out = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        _, xi = xi_from_flux(F[idx], dv, rmax_kms=rmax)
        out.append(xi)
    return np.std(np.array(out), axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("caches", nargs="+")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--norm", default="taueff", choices=["taueff", "none"])
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--rmax", type=float, default=None)
    ap.add_argument("--nboot", type=int, default=200)
    ap.add_argument("--out", default="figures/xi")
    args = ap.parse_args()

    labels = args.labels or [os.path.basename(p).replace("cache_", "")
                             .replace(".npz", "") for p in args.caches]
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    data = [(lab, cachelib.load(p)) for lab, p in zip(labels, args.caches)]
    z0, dv0 = data[0][1].z, data[0][1].dv
    npix = data[0][1].tau.shape[1]
    rmax = args.rmax if args.rmax else 0.5 * npix * dv0
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(z0)

    say("=" * 70)
    say(f"Stage 05 - flux correlation function    norm = {args.norm}")
    say("=" * 70)
    say(f"z = {z0:.3f}   dv = {dv0:.5f} km/s   rmax = {rmax:.1f} km/s")
    say()

    curves = {}
    for lab, c in data:
        A = solve_A(c.tau, target) if args.norm == "taueff" else 1.0
        F = flux(c.tau, A)
        r, xi = xi_from_flux(F, c.dv, rmax_kms=rmax)
        err = (bootstrap_xi(F, c.dv, rmax, args.nboot)
               if args.nboot else np.zeros_like(xi))
        curves[lab] = (r, xi, err)
        say(f"{lab:10s} A = {A:.5f}   xi(0) = {xi[0]:.5f}   "
            f"xi(100) = {np.interp(100.0, r, xi):.5f}   "
            f"xi(600) = {np.interp(600.0, r, xi):.5f}")
    say()

    ref = labels[0]
    r0 = curves[ref][0]
    say(f"{'dv [km/s]':>10s} " + " ".join(f"{l + '/' + ref:>14s}"
                                          for l in labels[1:]))
    for rr in (50, 100, 200, 400, 600, 1000, 1500):
        if rr > r0.max():
            continue
        row = []
        for lab in labels[1:]:
            a = np.interp(rr, *curves[lab][:2])
            b = np.interp(rr, r0, curves[ref][1])
            row.append(f"{a / b:14.4f}" if abs(b) > 1e-12 else f"{'--':>14s}")
        say(f"{rr:10.0f} " + " ".join(row))
    say()
    say("Read this against the P1D: a suppression over k < 0.01 s/km has to")
    say("show up on separations above ~600 km/s. If it does, the Fourier")
    say("result is not a product of the binning.")

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 8.0), sharex=True,
                                 gridspec_kw={"height_ratios": [2, 1],
                                              "hspace": 0.06})
    for i, lab in enumerate(labels):
        r, xi, err = curves[lab]
        col = PALETTE[i % len(PALETTE)]
        ax.plot(r[1:], xi[1:], color=col, lw=2.0, label=lab)
        if err.any():
            ax.fill_between(r[1:], (xi - err)[1:], (xi + err)[1:],
                            color=col, alpha=0.2, lw=0)
    ax.set(xscale="log", yscale="log", ylabel=r"$\xi_F(\Delta v)$",
           title=f"flux correlation function, z = {z0:.1f}, "
                 f"norm = {args.norm}")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2, which="both")

    for i, lab in enumerate(labels[1:], start=1):
        r, xi, _ = curves[lab]
        base = np.interp(r, r0, curves[ref][1])
        good = np.abs(base) > 1e-10
        bx.plot(r[good][1:], (xi[good] / base[good])[1:],
                color=PALETTE[i % len(PALETTE)], lw=2.0,
                label=f"{lab} / {ref}")
    bx.axhline(1.0, color="k", ls="--", lw=1)
    bx.set(xscale="log", xlabel=r"$\Delta v$ [km/s]",
           ylabel=f"ratio to {ref}", ylim=(0.5, 1.5))
    bx.legend(frameon=False)
    bx.grid(alpha=0.2, which="both")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"written -> {args.out}.png")
    print(f"written -> {args.out}.txt")


if __name__ == "__main__":
    main()
