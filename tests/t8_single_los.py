#!/usr/bin/env python3
"""
Test 8 - what does a single line of sight look like in both runs?

The two boxes start from IDENTICAL initial conditions and the sightlines are
drawn at the same transverse positions, so line i in CDM and line i in FCT
pass through the same structures.  That makes a pixel-by-pixel comparison
meaningful, which it would not be for two independent realisations.

The point of this test is to see the mechanism with your eyes before
trusting the ensemble average.  If the large-scale suppression comes from
gas being removed at the density peaks, then on a single sightline you
should see the FCT absorption troughs sitting where the CDM ones are but
shallower or missing, and the single-line P1D should lose power at low k
while roughly tracking CDM at high k.

Three outputs:

    <out>_profiles.png   tau and F along each selected sightline, both runs
    <out>_p1d.png        P1D of each single sightline, and its ratio, with
                         the full-sample ratio drawn behind for reference
    <out>.txt            the peak census: for every strong CDM absorber,
                         what FCT has at the same pixel

Usage
-----
    python tests/t8_single_los.py --cdm cache/cache_cdm.npz \\
                                  --fct cache/cache_fct.npz \\
                                  --los 12 340 900 --out figures/t8

    python tests/t8_single_los.py --cdm ... --fct ... --pick extreme

Options
-------
    --los I J K       explicit sightline indices (default: use --pick)
    --pick MODE       median | extreme | random | strongest  (default median)
                        median    lines whose individual ratio is closest to
                                  the ensemble ratio - the typical case
                        extreme   lines with the most suppressed ratio
                        strongest lines with the deepest CDM absorption
                        random    reproducible random draw (--seed)
    --n N             how many sightlines when using --pick (default 3)
    --tau-eff X       common rescaling target (default Turner+24 at cache z)
    --no-rescale      leave A = 1 in both runs
    --peak-tau X      what counts as a strong CDM absorber (default 2.0)
    --smooth-kms X    boxcar smoothing of the profiles for display only
    --seed N          for --pick random
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
                        solve_A)

C_CDM, C_FCT = "#1f77b4", "#d62728"


def per_los_ratio(Fa, Fb, dv, kmax, nbins=18):
    """Individual-sightline P1D ratio FCT/CDM, log-binned."""
    ka, Pa, _ = p1d_from_flux(Fa[None, :], dv)
    kb, Pb, _ = p1d_from_flux(Fb[None, :], dv)
    kba, pba = logbin(ka, Pa, nbins=nbins, kmax=kmax)
    kbb, pbb = logbin(kb, Pb, nbins=nbins, kmax=kmax)
    return kbb, pba / np.maximum(pbb, 1e-300), (kba, pba), (kbb, pbb)


def boxcar(x, dv, width_kms):
    if not width_kms or width_kms <= dv:
        return x
    n = max(int(round(width_kms / dv)), 1)
    ker = np.ones(n) / n
    return np.convolve(x, ker, mode="same")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cdm", required=True)
    ap.add_argument("--fct", required=True)
    ap.add_argument("--los", type=int, nargs="+", default=None)
    ap.add_argument("--pick", default="median",
                    choices=["median", "extreme", "random", "strongest"])
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--no-rescale", action="store_true")
    ap.add_argument("--peak-tau", type=float, default=2.0)
    ap.add_argument("--smooth-kms", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="figures/t8")
    args = ap.parse_args()

    a, b = cachelib.load_pair(args.cdm, args.fct)
    dv, z = a.dv, a.z
    nlos, npix = a.tau.shape
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(z)

    if args.no_rescale:
        A_cdm = A_fct = 1.0
    else:
        A_cdm, A_fct = solve_A(a.tau, target), solve_A(b.tau, target)

    F_cdm = flux(a.tau, A_cdm)
    F_fct = flux(b.tau, A_fct)

    kfull, Pc, _ = p1d_from_flux(F_cdm, dv)
    _, Pf, _ = p1d_from_flux(F_fct, dv)
    kmax = kfull[nyquist_cut(kfull, dv, 0.5)].max()
    kb_ens, pc_ens = logbin(kfull, Pc, nbins=24, kmax=kmax)
    _, pf_ens = logbin(kfull, Pf, nbins=24, kmax=kmax)
    r_ens = pf_ens / pc_ens

    # ---- pick the sightlines ---------------------------------------------
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    say("=" * 70)
    say("T8 - single sightlines, matched between the two runs")
    say("=" * 70)
    say(f"z = {z:.3f}   dv = {dv:.5f} km/s   {nlos} LOS x {npix} px")
    say(f"rescaling: A(CDM) = {A_cdm:.5f}   A(FCT) = {A_fct:.5f}"
        + ("   (none)" if args.no_rescale else f"   target {target:.5f}"))
    say()

    if args.los:
        sel = list(args.los)
    else:
        klow = kb_ens < 0.01
        rr = np.empty(nlos)
        for i in range(nlos):
            _, r_i, _, _ = per_los_ratio(F_fct[i], F_cdm[i], dv, kmax)
            rr[i] = np.nanmean(r_i[klow])
        if args.pick == "median":
            tgt = np.nanmean(r_ens[klow])
            sel = list(np.argsort(np.abs(rr - tgt))[:args.n])
        elif args.pick == "extreme":
            sel = list(np.argsort(rr)[:args.n])
        elif args.pick == "strongest":
            sel = list(np.argsort(-a.tau.max(axis=1))[:args.n])
        else:
            sel = list(np.random.default_rng(args.seed)
                       .choice(nlos, size=args.n, replace=False))
        say(f"pick = {args.pick}: sightlines {sel}")
        say(f"low-k ratio of the full sample = "
            f"{np.nanmean(r_ens[klow]):.4f}, "
            f"spread across sightlines = {np.nanstd(rr):.4f}")
        say()

    v = np.arange(npix) * dv

    # ---- profiles figure --------------------------------------------------
    fig1, axes = plt.subplots(2 * len(sel), 1,
                              figsize=(11, 2.6 * 2 * len(sel)), sharex=True)
    axes = np.atleast_1d(axes)
    for j, i in enumerate(sel):
        ax_t, ax_f = axes[2 * j], axes[2 * j + 1]
        tc = boxcar(a.tau[i], dv, args.smooth_kms)
        tf = boxcar(b.tau[i], dv, args.smooth_kms)
        ax_t.plot(v, tc, color=C_CDM, lw=1.0, label="CDM")
        ax_t.plot(v, tf, color=C_FCT, lw=1.0, label="FCT")
        ax_t.set(yscale="log", ylabel=r"$\tau$")
        ax_t.set_title(f"line of sight {i}", fontsize=10, loc="left")
        ax_t.legend(frameon=False, fontsize=8, ncol=2)
        ax_t.grid(alpha=0.2)

        peaks = a.tau[i] > args.peak_tau
        ax_f.plot(v, F_cdm[i], color=C_CDM, lw=1.0)
        ax_f.plot(v, F_fct[i], color=C_FCT, lw=1.0)
        if peaks.any():
            ax_f.fill_between(v, 0, 1, where=peaks, color="0.85", zorder=0,
                              label=rf"CDM $\tau>${args.peak_tau:g}")
            ax_f.legend(frameon=False, fontsize=8, loc="lower right")
        ax_f.set(ylim=(-0.05, 1.05), ylabel="$F$")
        ax_f.grid(alpha=0.2)
    axes[-1].set_xlabel(r"$v$ [km/s]")
    fig1.suptitle("T8: the same structures seen by both runs "
                  "(shaded = strong CDM absorbers)", fontsize=12)
    fig1.tight_layout(rect=(0, 0, 1, 0.985))

    # ---- P1D figure -------------------------------------------------------
    fig2, (bx, cx) = plt.subplots(2, 1, figsize=(8.2, 8.6), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1],
                                               "hspace": 0.06})
    bx.plot(kb_ens, kb_ens * pc_ens / np.pi, color=C_CDM, lw=2.6,
            label=f"CDM, all {nlos} LOS")
    bx.plot(kb_ens, kb_ens * pf_ens / np.pi, color=C_FCT, lw=2.6,
            label=f"FCT, all {nlos} LOS")
    styles = ["-", "--", ":", "-."]
    for j, i in enumerate(sel):
        kk, r_i, (ka, pa), (kbb, pb) = per_los_ratio(F_fct[i], F_cdm[i],
                                                     dv, kmax)
        ls = styles[j % len(styles)]
        bx.plot(kbb, kbb * pb / np.pi, color=C_CDM, lw=1.0, ls=ls, alpha=0.75)
        bx.plot(ka, ka * pa / np.pi, color=C_FCT, lw=1.0, ls=ls, alpha=0.75)
        cx.plot(kk, r_i, color="0.35", lw=1.3, ls=ls, label=f"LOS {i}")
    bx.axvspan(*units.desi_window(z), color="green", alpha=0.07, zorder=0)
    bx.set(xscale="log", yscale="log", ylabel=r"$k\,P_{\rm 1D}/\pi$",
           title="T8: single sightlines (thin) against the ensemble (thick)")
    bx.legend(frameon=False, fontsize=9)
    bx.grid(alpha=0.2, which="both")

    cx.plot(kb_ens, r_ens, color=C_FCT, lw=2.6, label="ensemble")
    cx.axhline(1.0, color="k", ls=":", lw=1)
    cx.axvspan(*units.desi_window(z), color="green", alpha=0.07, zorder=0)
    cx.set(xscale="log", xlabel=r"$k$ [s/km]", ylabel="FCT / CDM",
           ylim=(0.0, 2.2))
    cx.legend(frameon=False, fontsize=9, ncol=2)
    cx.grid(alpha=0.2, which="both")

    # ---- peak census ------------------------------------------------------
    say("PEAK CENSUS")
    say(f"For every pixel with tau_CDM > {args.peak_tau:g}, what does FCT "
        f"have at the same pixel?")
    say()
    say(f"{'LOS':>6s} {'N peaks':>8s} {'<tau_CDM>':>10s} {'<tau_FCT>':>10s} "
        f"{'ratio':>7s} {'FCT<0.1':>8s}")
    for i in list(sel) + ["ALL"]:
        if i == "ALL":
            mc, mf = a.tau, b.tau
            label = "ALL"
        else:
            mc, mf = a.tau[i], b.tau[i]
            label = str(i)
        pk = mc > args.peak_tau
        n = int(pk.sum())
        if n == 0:
            say(f"{label:>6s} {0:8d}")
            continue
        tc, tf = mc[pk].mean(), mf[pk].mean()
        weak = 100.0 * (mf[pk] < 0.1).mean()
        say(f"{label:>6s} {n:8d} {tc:10.3f} {tf:10.3f} "
            f"{tf / tc:7.3f} {weak:7.2f}%")
    say()
    say("A ratio well below 1 in that table is the mechanism in one number:")
    say("the absorber is at the same place in both runs, but in FCT the gas")
    say("that produced it has been converted, so the trough is shallower.")
    say("The 'FCT<0.1' column counts absorbers that have essentially")
    say("disappeared.")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig1.savefig(args.out + "_profiles.png", dpi=140, bbox_inches="tight")
    fig2.savefig(args.out + "_p1d.png", dpi=150, bbox_inches="tight")
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    for s in ("_profiles.png", "_p1d.png", ".txt"):
        print(f"written -> {args.out}{s}")


if __name__ == "__main__":
    main()
