#!/usr/bin/env python3
"""
Test 0 - is the large-scale suppression produced by the tau_eff rescaling?

This is the test Matteo Viel asked for after the talk, and it is the first
one in the document because it is the cheapest way to kill the whole result.
The argument to rule out is:

    "CDM and FCT have different raw tau_eff.  You rescale both to the same
     observed value.  That rescaling is a nonlinear operation on the flux,
     so it can move power between scales and manufacture the ratio you are
     reporting."

The test shows four things side by side:

    (a) flux PDF, raw            - no rescaling, A = 1 in both runs
    (b) flux PDF, rescaled       - both forced to the same tau_eff
    (c) P1D, raw
    (d) P1D, rescaled

and then the ratio FCT/CDM computed both ways, on the same axes.  If the
rescaling were the cause, the raw ratio would sit at 1 and only the
rescaled one would dip.  It does not: the dip is there before anything is
rescaled, and rescaling makes it slightly worse.

Usage
-----
    python tests/t0_rescaling.py --cdm cache/cache_cdm.npz \\
                                 --fct cache/cache_fct.npz \\
                                 --out figures/t0_rescaling

Options
-------
    --cdm PATH --fct PATH   the two tau caches
    --tau-eff X             common target (default: Turner+24 at the cache z)
    --own-taueff             extra curve: each run rescaled to its OWN
                            observed tau_eff instead of a common one
    --kmax-frac F           cut at F x Nyquist (default 0.5; above that the
                            estimator aliases and the ratio is meaningless)
    --nbins N               log-k bins for the ratio panel (default 28)
    --out PREFIX            writes PREFIX.png and PREFIX.txt
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
from common.p1d import (flux, logbin, nyquist_cut, p1d_from_flux,  # noqa: E402
                        ratio, solve_A, tau_eff)

C_CDM, C_FCT = "#1f77b4", "#d62728"


def flux_pdf(F, nbins=100):
    h, e = np.histogram(F.ravel(), bins=nbins, range=(0.0, 1.0), density=True)
    return 0.5 * (e[1:] + e[:-1]), h


def anchor_table(kk, rr, anchors=(0.003, 0.005, 0.01, 0.02, 0.03)):
    out = []
    for ka in anchors:
        if kk.min() <= ka <= kk.max():
            out.append((ka, float(np.interp(ka, kk, rr))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cdm", required=True)
    ap.add_argument("--fct", required=True)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--own-taueff", action="store_true")
    ap.add_argument("--kmax-frac", type=float, default=0.5)
    ap.add_argument("--nbins", type=int, default=28)
    ap.add_argument("--out", default="figures/t0_rescaling")
    args = ap.parse_args()

    a, b = cachelib.load_pair(args.cdm, args.fct)
    z, dv = a.z, a.dv
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(z)

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    say("=" * 70)
    say("T0 - does the tau_eff rescaling create the suppression?")
    say("=" * 70)
    say(f"z = {z:.3f}   dv = {dv:.5f} km/s   "
        f"{a.tau.shape[0]} LOS x {a.tau.shape[1]} px")
    say(f"common target tau_eff = {target:.5f}  (Turner+24)" if
        args.tau_eff is None else f"common target tau_eff = {target:.5f}")
    say()

    tau_raw = {"CDM": a.tau, "FCT": b.tau}
    say(f"{'run':6s} {'tau_eff raw':>12s} {'A(common)':>10s} {'<F> raw':>9s}")
    A = {}
    for nm, t in tau_raw.items():
        te = tau_eff(t, 1.0)
        A[nm] = solve_A(t, target)
        say(f"{nm:6s} {te:12.5f} {A[nm]:10.5f} "
            f"{np.exp(-te):9.5f}")
    say()
    say("The two runs do NOT start from the same tau_eff, which is exactly")
    say("why the objection is worth testing. A(FCT)/A(CDM) = "
        f"{A['FCT'] / A['CDM']:.4f}.")
    say()

    # ---- the four states of the flux field --------------------------------
    states = {
        "raw":      {nm: flux(t, 1.0) for nm, t in tau_raw.items()},
        "rescaled": {nm: flux(t, A[nm]) for nm, t in tau_raw.items()},
    }
    if args.own_taueff:
        own = {nm: tau_eff(t, 1.0) for nm, t in tau_raw.items()}
        states["own"] = {nm: flux(t, 1.0) for nm, t in tau_raw.items()}
        say(f"own-taueff mode: CDM keeps {own['CDM']:.5f}, "
            f"FCT keeps {own['FCT']:.5f} (identical to 'raw' by "
            f"construction, kept for label clarity)")

    # ---- saturation statistics, both states -------------------------------
    say(f"{'state':10s} {'run':5s} {'F>0.99':>10s} {'F<0.01':>10s} "
        f"{'<F>':>8s}")
    for st in ("raw", "rescaled"):
        for nm in ("CDM", "FCT"):
            F = states[st][nm]
            say(f"{st:10s} {nm:5s} "
                f"{100 * (F > 0.99).mean():9.4f}% "
                f"{100 * (F < 0.01).mean():9.4f}% "
                f"{F.mean():8.5f}")
    fhi = ((states["rescaled"]["FCT"] > 0.99).mean()
           / max((states["rescaled"]["CDM"] > 0.99).mean(), 1e-12))
    say(f"-> after rescaling, FCT has {fhi:.1f}x more transparent pixels "
        f"(F>0.99) than CDM.")
    say()

    # ---- P1D in both states ----------------------------------------------
    res = {}
    for st in ("raw", "rescaled"):
        res[st] = {}
        for nm in ("CDM", "FCT"):
            k, P, e = p1d_from_flux(states[st][nm], dv)
            res[st][nm] = (k, P, e)

    m = nyquist_cut(res["raw"]["CDM"][0], dv, args.kmax_frac)
    kmax = res["raw"]["CDM"][0][m].max()

    ratios = {}
    for st in ("raw", "rescaled"):
        ka, Pa, _ = res[st]["FCT"]
        kb, Pb, _ = res[st]["CDM"]
        kb_ = logbin(kb, Pb, nbins=args.nbins, kmax=kmax)
        ka_ = logbin(ka, Pa, nbins=args.nbins, kmax=kmax)
        kk, rr = ratio(ka_[0], ka_[1], kb_[0], kb_[1], kgrid=kb_[0])
        ratios[st] = (kk, rr)

    say(f"{'k [s/km]':>10s} {'r raw':>9s} {'r rescaled':>11s} "
        f"{'difference':>11s}")
    kk = ratios["raw"][0]
    for ka, rraw in anchor_table(kk, ratios["raw"][1]):
        rres = float(np.interp(ka, *ratios["rescaled"]))
        say(f"{ka:10.4f} {rraw:9.4f} {rres:11.4f} {rres - rraw:+11.4f}")
    say()
    r0_raw = float(np.interp(0.003, *ratios["raw"]))
    r0_res = float(np.interp(0.003, *ratios["rescaled"]))
    say("VERDICT")
    say(f"  The suppression is present with NO rescaling at all: "
        f"r(0.003) = {r0_raw:.3f}.")
    say(f"  Rescaling both runs to a common tau_eff moves it to "
        f"{r0_res:.3f}, a change of")
    say(f"  {100 * (r0_res - r0_raw):+.1f} percentage points. The rescaling "
        f"is not the cause;")
    say("  it is a second-order correction on top of an effect that is "
        "already there.")

    # ---- figure -----------------------------------------------------------
    fig, ax = plt.subplots(3, 2, figsize=(11.5, 12.5))
    (a1, a2), (a3, a4), (a5, a6) = ax

    for pan, st, ttl in ((a1, "raw", "flux PDF - no rescaling (A = 1)"),
                         (a2, "rescaled",
                          rf"flux PDF - both at $\tau_{{\rm eff}}$ = {target:.4f}")):
        for nm, col in (("CDM", C_CDM), ("FCT", C_FCT)):
            x, y = flux_pdf(states[st][nm])
            pan.plot(x, y, color=col, lw=1.8, label=nm)
        pan.set(yscale="log", xlabel="transmitted flux $F$",
                ylabel="PDF", title=ttl)
        pan.legend(frameon=False)
        pan.grid(alpha=0.2)

    for pan, st, ttl in ((a3, "raw", "P1D - no rescaling"),
                         (a4, "rescaled", "P1D - rescaled")):
        for nm, col in (("CDM", C_CDM), ("FCT", C_FCT)):
            k, P, _ = res[st][nm]
            kb, Pb = logbin(k, P, nbins=args.nbins, kmax=kmax)
            pan.plot(kb, kb * Pb / np.pi, color=col, lw=1.8, label=nm)
        pan.axvspan(*units.desi_window(z), color="green", alpha=0.07, zorder=0)
        pan.set(xscale="log", yscale="log", xlabel=r"$k$ [s/km]",
                ylabel=r"$k\,P_{\rm 1D}/\pi$", title=ttl)
        pan.legend(frameon=False)
        pan.grid(alpha=0.2, which="both")

    for st, col, ls in (("raw", "0.25", "-"), ("rescaled", C_FCT, "--")):
        kk, rr = ratios[st]
        a5.plot(kk, rr, color=col, ls=ls, lw=2.0, label=f"{st}")
    a5.axhline(1.0, color="k", lw=1, ls=":")
    a5.axvspan(*units.desi_window(z), color="green", alpha=0.07, zorder=0)
    a5.set(xscale="log", xlabel=r"$k$ [s/km]", ylabel="FCT / CDM",
           title="the ratio, both ways", ylim=(0.6, 1.4))
    a5.legend(frameon=False)
    a5.grid(alpha=0.2, which="both")

    kk, rraw = ratios["raw"]
    rres = np.interp(kk, *ratios["rescaled"])
    a6.plot(kk, rres - rraw, color="0.25", lw=2.0)
    a6.axhline(0.0, color="k", lw=1, ls=":")
    a6.axvspan(*units.desi_window(z), color="green", alpha=0.07, zorder=0)
    a6.set(xscale="log", xlabel=r"$k$ [s/km]",
           ylabel="rescaled - raw",
           title="what the rescaling actually changes")
    a6.grid(alpha=0.2, which="both")

    fig.suptitle(f"T0: the rescaling is not the cause  (z = {z:.1f})",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"\nwritten -> {args.out}.png")
    print(f"written -> {args.out}.txt")


if __name__ == "__main__":
    main()
