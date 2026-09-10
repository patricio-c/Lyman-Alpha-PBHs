#!/usr/bin/env python3
"""
Stage 11 - is the low-k deficit the gas drain, propagated correctly?

Stage 10 fitted DP1D = -f P1D_CDM + c and found f = 0.267 +- 0.008 with a
chi^2/dof of 13.3 and residuals whose signs run - - - + + + + + -.  That is
not scatter.  The constant-f assumption is the part that is wrong: removing
half the baryons changes the SHAPE of the gas distribution, not only its
amplitude, so the suppression should be scale dependent.

This stage stops assuming the shape and measures it.  Stage 09 gives the
3D gas power of both runs on the same grid, so the suppression

    g(k) = 1 - P_gas_FCT(k) / P_gas_CDM(k)

is a measurement, not a parameter.  The question is whether that measured
suppression, pushed through the P1D integral, reproduces the low-k deficit.

The chain, and why each step is what it is
------------------------------------------
1.  P1D and P3D are related by

        P1D(k_par) = (1/2pi) int_{k_par}^{inf} k P3D(k) dk

    which inverts exactly:

        P3D(k) = -(2 pi / k) dP1D/dk

    So P3D of the FLUX is recoverable from the measured P1D_CDM.  This
    matters: the drain must be applied to the flux power, not to the gas
    power, and nothing in the pipeline measures 3D flux power directly.
    Differentiating noisy data amplifies noise, so P1D_CDM is first fitted
    by a smoothing spline in log-log and the derivative taken from the fit.
    --spline-s controls the smoothing and --diagnose plots the round trip.

2.  Suppress it by the measured gas suppression and integrate back:

        DP1D_drain(k_par) = -(1/2pi) int_{k_par}^{inf} k P3D_flux(k) g(k) dk

    The one physical assumption, stated plainly: the fractional suppression
    of flux power at a 3D scale k equals the fractional suppression of gas
    power at the same k.  That is an assumption about how a density change
    maps into an absorption change, and it is exactly what this stage
    tests.  It is not guaranteed - saturation makes the flux respond
    sublinearly in dense gas - so a failure here is informative rather than
    fatal.

3.  Fit

        M3   DP1D = alpha * DP1D_drain(k) + c

    against the same bootstrap covariance stage 10 used.  Two parameters,
    the same count as M2, but the drain now has a MEASURED shape and only
    its amplitude floats.

What the answer means
---------------------
    alpha ~ 1 and chi^2 improves    the deficit IS the gas drain, correctly
                                    propagated. The two-regime picture is
                                    confirmed with the drain measured
                                    rather than assumed, and c is a clean
                                    estimate of the primordial boost.
    alpha far from 1                the shape is right but the amplitude is
                                    not, so the flux does not respond to
                                    gas removal in proportion. Report the
                                    number; it is the saturation correction.
    chi^2 no better than M2         the low-k deficit is not the gas drain.
                                    The mechanism is then still open and
                                    the census result does not explain the
                                    forest result.

Truncation
----------
The integral runs to infinity; the data stop at the Nyquist frequency of
the stage 09 grid.  This is not bookkeeping.  For the measured 3D slope of
-2.55 and a grid Nyquist of 27 Mpc^-1, the part of the integral above the
grid is 5% of the total at the bottom of the DESI window and 24% at the
top, so simply dropping it tilts the drain shape across exactly the band
where the shape is being tested.  The tail is therefore continued as a
power law fitted to the top decade, which is exact when the power law
holds, and refused outright when the fitted slope is >= -2 and the
integral does not converge.  --no-extrapolate reverts to truncation so
that the size of the bias can be measured; it is not a way to produce a
result.

Usage
-----
    python stages/11_drain_shape.py \\
        figures/dp1d_40_boot.npz figures/pk_z3_CDM.txt figures/pk_z3_FCT.txt \\
        --out figures/drain_40

Options
-------
    boot_npz               PREFIX_boot.npz written by stage 10
    pk_ref pk_test         stage 09 tables, reference first
    --spline-s S           smoothing of the log-log P1D fit (default 0.002)
    --no-extrapolate       truncate the integral at the last measured k
                           instead of continuing the tail as a power law.
                           Biases the drain shape; use it to measure that
                           bias, not to produce a result.
    --diagnose             extra panel showing the P1D -> P3D -> P1D round trip
    --out PREFIX           writes PREFIX.png and PREFIX.txt
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.interpolate import UnivariateSpline

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import units  # noqa: E402


def read_pk(path):
    """
    Stage 09 table -> dict of columns.

    Header is
        # k  Nmodes  P_gas  P_conv  P_dm  P_matter  P_baryon  P_gc
    with k in inverse internal length units, i.e. Mpc^-1 for these runs.
    The 8th column is absent when the run was given --no-cross, so the
    reader takes what is there rather than insisting on the full width.
    """
    d = np.loadtxt(path)
    names = ["k", "nmodes", "gas", "conv", "dm", "matter", "baryon", "gc"]
    if d.shape[1] < 7:
        raise SystemExit(f"{path}: expected at least 7 columns, "
                         f"got {d.shape[1]}")
    out = {n: d[:, i] for i, n in enumerate(names[:d.shape[1]])}
    meta = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "label" in line:
                t = line.split()
                meta["label"] = t[t.index("label") + 1]
                meta["z"] = float(t[t.index("z") + 1])
                meta["box"] = float(t[t.index("box") + 1])
    out["_meta"] = meta
    return out


def p3d_from_p1d(k, P1D, s):
    """
    Invert P1D(k) = (1/2pi) int_k^inf k' P3D(k') dk' by differentiation.

    Differentiating measured data is the step that goes wrong quietly, so
    it is done on a smoothing spline of log P1D against log k and the
    derivative is taken analytically:

        dP1D/dk = (P1D / k) * dlnP1D/dlnk
        P3D(k)  = -(2 pi / k) dP1D/dk
                = -(2 pi / k^2) P1D dlnP1D/dlnk

    P1D falls with k, so the log slope is negative and P3D comes out
    positive.  A positive slope anywhere means the spline is following
    noise and the caller is told.
    """
    x, y = np.log(k), np.log(P1D)
    sp = UnivariateSpline(x, y, s=s, k=3)
    slope = sp.derivative()(x)
    P3D = -(2.0 * np.pi / k ** 2) * P1D * slope
    return P3D, slope, sp


def p1d_from_p3d(k, P3D, tail=0.0):
    """
    Forward integral, cumulative from the top of the band downwards.

    Trapezoid in log k, because the integrand k^2 P3D is what is smooth
    there.  Accurate to ~1e-5 against the analytic result for power-law
    P3D with slopes -2.6 to -3.5.

    `tail` is the part of the integral above the last k, which is a
    CONSTANT added at every k below it.  It is not optional bookkeeping:
    for the measured 3D slope of -2.55 and a grid Nyquist of 27 Mpc^-1,
    the missing tail is 5% of the integral at the bottom of the DESI
    window and 24% at the top, so dropping it tilts the drain shape in
    exactly the band where the shape is being tested.
    """
    integrand = k ** 2 * P3D          # k P3D dk = k^2 P3D dlnk
    lk = np.log(k)
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(lk)
    cum = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
    return cum / (2.0 * np.pi) + tail


def powerlaw_tail(k, P3D, ndec=1.0):
    """
    (tail integral, fitted slope) for int_{kmax}^inf k P3D dk / 2pi.

    P3D is continued above the last measured k as a power law fitted to
    the top `ndec` decades.  For P3D = P_max (k/kmax)^s with s < -2,

        int_{kmax}^inf k P3D dk  =  kmax^2 P_max / |s + 2|

    The slope decides whether this is meaningful at all: at s >= -2 the
    integral diverges, which means the band has not reached the regime
    where the power is falling fast enough to be truncated safely.  In
    that case NaN comes back and the caller must say so rather than
    quietly using a number.
    """
    m = k >= k[-1] / 10.0 ** ndec
    if int(m.sum()) < 3 or (P3D[m] <= 0).any():
        return 0.0, np.nan
    s = float(np.polyfit(np.log(k[m]), np.log(P3D[m]), 1)[0])
    if s >= -2.0:
        return np.nan, s
    return float(k[-1] ** 2 * P3D[-1] / (2.0 * np.pi * abs(s + 2.0))), s


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("boot_npz")
    ap.add_argument("pk_ref")
    ap.add_argument("pk_test")
    ap.add_argument("--spline-s", type=float, default=0.002)
    ap.add_argument("--diagnose", action="store_true")
    ap.add_argument("--no-extrapolate", action="store_true")
    ap.add_argument("--out", default="figures/drain")
    args = ap.parse_args()

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    B = np.load(args.boot_npz)
    kb = B["k"]
    Pa, dP, sd = B["P_ref"], B["dP"], B["err_dP"]
    fm = B["fit_mask"]
    z = float(B["z"])
    boot = B["boot_dP"].astype(np.float64)
    exact = bool(B["exact_bootstrap"]) if "exact_bootstrap" in B else False

    A = read_pk(args.pk_ref)
    C = read_pk(args.pk_test)

    say("=" * 74)
    say("Stage 11 - the drain shape, measured rather than assumed")
    say("=" * 74)
    say(f"P1D   : {args.boot_npz}   z = {z:.4f}   "
        f"{'exact' if exact else 'FROZEN-A (provisional)'} bootstrap")
    say(f"P3D   : {args.pk_ref}  [{A['_meta'].get('label','?')}]")
    say(f"        {args.pk_test}  [{C['_meta'].get('label','?')}]")
    if not exact:
        say("  WARNING: the stage 10 npz was written without --exact-boot, "
            "so every chi2 below inherits errors that are known to be "
            "underestimated.")
    if not np.allclose(A["k"], C["k"]):
        raise SystemExit("the two stage 09 tables are on different k grids")
    say()

    # --- the measured suppression -----------------------------------------
    k3 = A["k"]
    ratio_gas = C["gas"] / A["gas"]
    g = 1.0 - ratio_gas
    say("measured gas-power suppression  g(k) = 1 - P_gas_test / P_gas_ref")
    say(f"{'k [1/Mpc]':>10s} {'P_gas ref':>12s} {'P_gas test':>12s} "
        f"{'ratio':>8s} {'g':>8s}")
    for j in np.linspace(0, k3.size - 1, 10).astype(int):
        say(f"{k3[j]:10.4f} {A['gas'][j]:12.5e} {C['gas'][j]:12.5e} "
            f"{ratio_gas[j]:8.4f} {g[j]:8.4f}")
    say(f"g at the largest scale measured (k = {k3[0]:.4f} Mpc^-1): "
        f"{g[0]:.4f}")
    say()

    # --- P1D -> P3D --------------------------------------------------------
    kfac = units.DEFAULT.k_kms_to_mpc(z)
    k1_mpc = kb * kfac
    P3D, slope, sp = p3d_from_p1d(k1_mpc, Pa, args.spline_s)
    bad = slope >= 0
    say(f"P1D -> P3D by differentiating a smoothing spline (s = "
        f"{args.spline_s})")
    say(f"  log-slope range: {slope.min():.3f} to {slope.max():.3f}")
    if bad.any():
        say(f"  WARNING: the fitted slope is >= 0 in {int(bad.sum())} bins, "
            f"so P3D is negative there. Raise --spline-s.")
    t_ref, s_ref = powerlaw_tail(k1_mpc, P3D)
    if args.no_extrapolate:
        t_ref = 0.0
    round_trip = p1d_from_p3d(k1_mpc, P3D, tail=0.0 if np.isnan(t_ref)
                              else t_ref)
    resid = round_trip / np.maximum(Pa, 1e-300)
    say(f"  round trip P1D -> P3D -> P1D recovers "
        f"{float(np.median(resid[:-1])):.4f} of P1D (median)")
    say(f"  tail slope above k = {k1_mpc[-1]:.2f} Mpc^-1: s = {s_ref:.3f}")
    say()

    # --- the drain, propagated --------------------------------------------
    gi = np.interp(np.log(k1_mpc), np.log(k3), g,
                   left=g[0], right=g[-1])
    n_extrap = int((k1_mpc < k3[0]).sum() + (k1_mpc > k3[-1]).sum())
    if n_extrap:
        say(f"  note: {n_extrap} of {kb.size} P1D bins lie outside the "
            f"measured 3D band [{k3[0]:.3f}, {k3[-1]:.3f}] Mpc^-1; "
            f"g is held flat there.")
    t_dr, s_dr = powerlaw_tail(k1_mpc, P3D * gi)
    if args.no_extrapolate:
        say("  --no-extrapolate: the integral is truncated at the last "
            "measured k. The drain shape below is biased; this mode is for "
            "measuring how much, not for producing a result.")
        t_dr = 0.0
    elif np.isnan(t_dr):
        say(f"  WARNING: the drain integrand's tail slope is s = {s_dr:.3f} "
            f">= -2, so the integral above k = {k1_mpc[-1]:.2f} Mpc^-1 does "
            f"not converge.")
        say("           The band has not reached the falling regime and the "
            "tail cannot be continued. Truncating instead; treat every "
            "number below as a lower bound.")
        t_dr = 0.0
    else:
        say(f"  tail continued as a power law, slope s = {s_dr:.3f}, "
            f"contributing {t_dr:.5e} km/s at every k")

    cum = p1d_from_p3d(k1_mpc, P3D * gi, tail=0.0)
    drain = -(cum + t_dr)
    say("  fraction of the drain integral carried by that tail:")
    for jj in np.flatnonzero(fm):
        tot = cum[jj] + t_dr
        say(f"    k = {kb[jj]:8.5f} s/km = {k1_mpc[jj]:7.3f} Mpc^-1   "
            f"{t_dr / max(tot, 1e-300) * 100:5.1f}%")
    say()

    # --- fit ---------------------------------------------------------------
    nf = int(fm.sum())
    Cov = np.cov(boot[:, fm], rowvar=False)
    nbt = boot.shape[0]
    hl = (nbt - nf - 2) / (nbt - 1.0)
    if nbt < 3 * nf or hl <= 0:
        Cinv = np.diag(1.0 / sd[fm] ** 2)
        say(f"covariance: diagonal ({nbt} resamples is too few for {nf} bins)")
    else:
        Cinv = hl * np.linalg.inv(Cov)
        say(f"covariance: full, Hartlap {hl:.4f}, {nbt} resamples")

    y = dP[fm]
    X3 = np.column_stack([drain[fm], np.ones(nf)])
    cov3 = np.linalg.inv(X3.T @ Cinv @ X3)
    beta3 = cov3 @ (X3.T @ Cinv @ y)
    r3 = y - X3 @ beta3
    chi3 = float(r3 @ Cinv @ r3)
    alpha, c3 = beta3

    # M2 refitted here on the same covariance, so the comparison is fair
    X2 = np.column_stack([-Pa[fm], np.ones(nf)])
    cov2 = np.linalg.inv(X2.T @ Cinv @ X2)
    beta2 = cov2 @ (X2.T @ Cinv @ y)
    r2 = y - X2 @ beta2
    chi2 = float(r2 @ Cinv @ r2)

    say()
    say("M2  DP1D = -f P1D_ref + c        (constant drain, refitted here)")
    say(f"    f     = {beta2[0]:.5f} +- {np.sqrt(cov2[0,0]):.5f}")
    say(f"    c     = {beta2[1]:.5e} +- {np.sqrt(cov2[1,1]):.3e} km/s")
    say(f"    chi2/dof = {chi2:.2f} / {nf-2} = {chi2/(nf-2):.3f}")
    say()
    say("M3  DP1D = alpha * drain(k) + c  (MEASURED drain shape)")
    say(f"    alpha = {alpha:.5f} +- {np.sqrt(cov3[0,0]):.5f}")
    say(f"    c     = {c3:.5e} +- {np.sqrt(cov3[1,1]):.3e} km/s")
    say(f"    chi2/dof = {chi3:.2f} / {nf-2} = {chi3/(nf-2):.3f}")
    say()

    say(f"M2 vs M3: delta chi2 = {chi2 - chi3:.2f} at equal parameter count")
    if chi3 < chi2 - 4:
        say("  -> the measured drain shape beats a constant fraction. The "
            "low-k deficit tracks the gas suppression.")
    elif chi3 > chi2 + 4:
        say("  -> the measured drain shape is WORSE than a constant "
            "fraction. The deficit does not follow the gas power.")
    else:
        say("  -> the two are not distinguishable with these errors.")

    say()
    say(f"alpha = {alpha:.3f}")
    if abs(alpha - 1.0) < 2 * np.sqrt(cov3[0, 0]):
        say("  consistent with 1. The flux power responds to gas removal in "
            "proportion, so the propagation needs no correction factor.")
    elif alpha < 1.0:
        say("  below 1. The flux is LESS suppressed than the gas, which is "
            "what saturation in dense absorbers would do.")
    else:
        say("  above 1. The flux is MORE suppressed than the gas; the "
            "removal is biased towards the gas the forest actually sees.")
    say()

    say(f"{'k [s/km]':>10s} {'DP1D':>12s} {'err':>10s} {'M2':>12s} "
        f"{'M3':>12s} {'(d-M3)/e':>9s}")
    for j, jj in enumerate(np.flatnonzero(fm)):
        say(f"{kb[jj]:10.5f} {dP[jj]:12.5e} {sd[jj]:10.3e} "
            f"{(X2 @ beta2)[j]:12.5e} {(X3 @ beta3)[j]:12.5e} "
            f"{r3[j]/sd[jj]:9.2f}")
    say()
    say("residual signs M3: " + "".join("+" if v > 0 else "-" for v in r3))
    say("residual signs M2: " + "".join("+" if v > 0 else "-" for v in r2))
    say("  a run of like signs means the model is missing a shape, not "
        "scatter, whatever the chi2 says.")

    # --- figure ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        npan = 3 if args.diagnose else 2
        fig, axes = plt.subplots(npan, 1, figsize=(7.6, 3.2 * npan + 1.4),
                                 sharex=False)
        ax = axes[0]
        ax.axhline(0, color="k", lw=1, ls="--")
        ax.errorbar(kb, dP, yerr=sd, fmt="o", ms=4, color="#1f77b4", lw=1.2,
                    capsize=2, label=r"$\Delta P_{\rm 1D}$")
        kk = kb[fm]
        ax.plot(kk, X2 @ beta2, color="#d62728", lw=1.8,
                label=rf"M2 const $f$, $\chi^2/\nu$={chi2/(nf-2):.2f}")
        ax.plot(kk, X3 @ beta3, color="#2ca02c", lw=2.0,
                label=rf"M3 measured drain, $\chi^2/\nu$={chi3/(nf-2):.2f}")
        ax.set(xscale="log", ylabel=r"$\Delta P_{\rm 1D}$ [km s$^{-1}$]",
               title=rf"$z={z:.1f}$,  $\alpha={alpha:.3f}$")
        ax.legend(frameon=False, fontsize=9)
        ax.grid(alpha=0.2, which="both")

        bx = axes[1]
        bx.axhline(1, color="k", lw=1, ls="--")
        bx.plot(k3, ratio_gas, color="#9467bd", lw=2,
                label=r"$P_{\rm gas}$ test / ref  (stage 09)")
        bx.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]",
               ylabel="3D gas power ratio")
        bx.legend(frameon=False, fontsize=9)
        bx.grid(alpha=0.2, which="both")

        if args.diagnose:
            cx = axes[2]
            cx.plot(k1_mpc, Pa, "o", ms=3, color="#1f77b4", label="P1D data")
            cx.plot(k1_mpc, round_trip, "-", color="#ff7f0e",
                    label="P1D from P3D round trip")
            cx.set(xscale="log", yscale="log", xlabel=r"$k$ [Mpc$^{-1}$]",
                   ylabel=r"$P_{\rm 1D}$")
            cx.legend(frameon=False, fontsize=9)
            cx.grid(alpha=0.2, which="both")

        fig.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                    exist_ok=True)
        fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
        print(f"written -> {args.out}.png")
    except Exception as exc:
        say(f"figure skipped: {exc.__class__.__name__}: {exc}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"written -> {args.out}.txt")


if __name__ == "__main__":
    main()
