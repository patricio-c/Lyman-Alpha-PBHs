#!/usr/bin/env python3
"""
Stage 11 - is the P1D difference just the 3D gas difference, propagated?

Stage 10 fitted DP1D = -f P1D_CDM + c, assuming the drain is a constant
fraction of the flux power.  It got chi^2/dof = 11.5 with residual signs
running - - - + + + + + -, which is a missing shape rather than scatter.

This stage stops assuming the shape and takes it from the data.  Stage 09
measures the 3D gas power of both runs, so

    g(k) = 1 - P_gas_test(k) / P_gas_ref(k)

is a measurement, not a parameter.  On the 40 Mpc/h pair g is POSITIVE at
large scales (the test run has less gas power) and NEGATIVE at small
scales (it has more).  That sign change matters: a single template forces
one amplitude to describe both effects, so the models below are built in
increasing order of how much they let the data speak.

    M2   DP1D = -f P1D_ref + c                   constant fraction
    M3   DP1D = alpha T[g] + c                   one measured template
    M4   DP1D = a_sup T[g+] + a_exc T[g-] + c    split, with a constant
    M5   DP1D = a_sup T[g+] + a_exc T[g-]        split, no constant

where T[.] is the propagation described below.

M4 and M5 exist as a pair for a numerical reason worth knowing before
reading their output.  T[g-] turns out to be nearly constant across the
fit window, so in M4 it is almost degenerate with the free constant, the
normal matrix is close to singular, and the fitted errors come back as
NaN.  Because g = g+ + g- exactly, span{T[g], 1} is a SUBSPACE of
span{T[g+], T[g-], 1}, so chi^2(M4) <= chi^2(M3) has to hold as algebra;
the stage checks that inequality explicitly and says so when it fails,
which is the signature of the degeneracy rather than of any physics.  M5
drops the redundant constant and is what the stage reports when M4 is not
conditioned.

The fit-free bound, which is the strongest thing here
-----------------------------------------------------
P1D is a positive-kernel integral of P3D,

    P1D(kp) = (1/2pi) int_kp^inf k P3D(k) dk

so the RATIO of two P1D curves is a weighted average of the ratio of their
P3D curves, with weights w(k) = k P3D_ref(k) > 0:

    r_1D(kp) = int w r_3D / int w

A weighted average with positive weights lies between the extremes of what
it averages:

    min_{k>kp} r_3D  <=  r_1D(kp)  <=  max_{k>kp} r_3D

Two things follow with no fitting and no free parameter.

*   The integral cannot manufacture a ratio below 1 out of 3D ratios that
    are all at or above 1.  A deficit in P1D therefore REQUIRES a deficit
    in 3D.  It is not an artefact of the integral weighting small scales.

*   If the measured r_1D falls OUTSIDE that band, then flux power and gas
    power are not changed by the same fraction, and no fit can argue that
    away.

The second point only counts when the violation is larger than the error
on r_1D, so the distance outside the band is reported in sigma and
nothing below 3 sigma is called a demonstration.  Note also that the
lowest P1D bin - the one carrying the deepest deficit, and so the most
informative - cannot be tested at all when the 3D band starts above it.
That is a box-size limitation, not an analysis choice, and a 3D
measurement reaching lower k would sharpen this test considerably.

The propagation
---------------
1.  Invert P1D_ref to 3D FLUX power.  The relation above inverts exactly,

        P3D(k) = -(2 pi / k) dP1D/dk = -(2 pi / k^2) P1D dlnP1D/dlnk

    which matters because the suppression has to act on flux power, and
    nothing in the pipeline measures 3D flux power directly.
    Differentiating data amplifies noise, so P1D is fitted by a smoothing
    spline in log-log and the derivative taken from the fit.  Bins whose
    fitted slope comes out >= 0 would give a NEGATIVE P3D, which is a
    corruption inside the integral rather than a small error; they are
    named and take the log-slope of their nearest well-behaved neighbour,
    since the bin that needs it is typically the endpoint, where a spline
    overshoots and where zeroing would throw away real signal.

2.  Weight by g(k) and integrate back.  The one physical assumption,
    stated plainly: the fractional change of flux power at a 3D scale k
    equals the fractional change of gas power at the same k.  Everything
    this stage measures is a test of that assumption.

--diagnose plots the closure of step 1 as the RATIO of the round trip to
the data, on a linear axis.  Plotting the two curves against each other
over eight decades of log axis, as an earlier version did, looks perfect
no matter what: the whole result lives at the percent level and a log
axis that wide cannot show a percent.

Bands and tails
---------------
The stage 09 tables run out to sqrt(3) k_Nyquist, the CORNERS of the
Fourier cube.  Those are real modes but poorly sampled, and the MAS
correction there is unreliable; on the 40 Mpc/h pair they produced a last
point that broke an otherwise monotonic trend.  The grid is read from the
table header and the band is cut at k_Nyquist.  --kcut-3d overrides.

The integral itself runs to infinity.  For a 3D slope of -2.55 and a grid
Nyquist of 27 Mpc^-1 the part above the grid is 5% of the total at the
bottom of the DESI window and 24% at the top, so truncating tilts the
result across exactly the band being tested.  The tail is continued as a
power law fitted to the top decade, and each template says which of three
things happened to it: continued, truncated because the integrand changes
sign, or exactly zero because the weight had already switched it off.

Usage
-----
    python stages/11_drain_shape.py \\
        figures/dp1d_40_boot.npz figures/pk_z3_CDM.txt figures/pk_z3_FCT.txt \\
        --diagnose --out figures/drain_40

Options
-------
    boot_npz               PREFIX_boot.npz written by stage 10
    pk_ref pk_test         stage 09 tables, reference first
    --spline-s S           smoothing of the log-log P1D fit (default 0.002)
    --kcut-3d X            cut the 3D band here [Mpc^-1] instead of at the
                           grid Nyquist read from the table header
    --no-extrapolate       truncate the integral at the last measured k
                           instead of continuing the tail. Biases the
                           result; use it to measure that bias.
    --diagnose             extra panel with the inversion closure
    --out PREFIX           writes PREFIX.png and PREFIX.txt
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.interpolate import UnivariateSpline

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import prov  # noqa: E402
from common import units  # noqa: E402


def read_pk(path):
    """
    Stage 09 table -> dict of columns plus a _meta dict.

    Header is
        # stage 09   label CDM   z 3.000000   box 58.737151
        # grid 512  MAS CIC
        ...
        # k  Nmodes  P_gas  P_conv  P_dm  P_matter  P_baryon  P_gc
    with k in inverse internal length units, i.e. Mpc^-1 for these runs.
    The 8th column is absent when the run was given --no-cross, so the
    reader takes what is there rather than insisting on the full width.

    The grid is parsed because k_Nyquist = pi * grid / box is the only way
    to know where the trustworthy part of the table ends: Pylians returns
    bins out to sqrt(3) k_Nyq, and the last 40% of the k range is cube
    corners.
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
            t = line.split()
            if "label" in t:
                meta["label"] = t[t.index("label") + 1]
                meta["z"] = float(t[t.index("z") + 1])
                meta["box"] = float(t[t.index("box") + 1])
            if "grid" in t:
                meta["grid"] = int(t[t.index("grid") + 1])
    if "grid" in meta and "box" in meta:
        meta["k_nyq"] = np.pi * meta["grid"] / meta["box"]
    out["_meta"] = meta
    return out


def p3d_from_p1d(k, P1D, s):
    """
    Invert P1D(k) = (1/2pi) int_k^inf k' P3D(k') dk' by differentiation.

    Differentiating measured data is the step that goes wrong quietly, so
    it is done on a smoothing spline of log P1D against log k and the
    derivative is taken analytically:

        dP1D/dk = (P1D / k) * dlnP1D/dlnk
        P3D(k)  = -(2 pi / k^2) P1D dlnP1D/dlnk

    P1D falls with k, so the log slope is negative and P3D comes out
    positive.  A slope >= 0 anywhere means the spline is following noise
    there, and the caller is given the offending k values.
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

    `tail` is the part of the integral above the last k, a CONSTANT added
    at every k below it.  Not optional bookkeeping: for a 3D slope of
    -2.55 and a grid Nyquist of 27 Mpc^-1 it is 5% of the integral at the
    bottom of the DESI window and 24% at the top.
    """
    integrand = k ** 2 * P3D          # k P3D dk = k^2 P3D dlnk
    lk = np.log(k)
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(lk)
    cum = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
    return cum / (2.0 * np.pi) + tail


def powerlaw_tail(k, P3D, ndec=1.0):
    """
    (tail integral, fitted slope, status) for int_{kmax}^inf k P3D dk/2pi.

    P3D is continued above the last measured k as a power law fitted to
    the top `ndec` decades.  For P3D = P_max (k/kmax)^s with s < -2,

        int_{kmax}^inf k P3D dk  =  kmax^2 P_max / |s + 2|

    Four outcomes, and they are genuinely different, which an earlier
    version got wrong by collapsing several of them into a bare NaN that
    the caller then misread:

        "ok"          fitted and convergent
        "zero"        the weight switches the integrand off before the top
                      of the band, so the tail is exactly zero. This is
                      the normal case for the suppression template, and it
                      is not a failure - nothing is being dropped.
        "divergent"   s >= -2, the integral does not converge; the band
                      has not reached the falling regime
        "unfittable"  the integrand changes sign over the top decade, so
                      no power law describes it. This is the normal case
                      for a SIGNED template.

    In the last two cases the tail is zero and the caller must say which
    happened, because they call for different things: more resolution in
    one, a split template in the other.
    """
    m = k >= k[-1] / 10.0 ** ndec
    if int(m.sum()) < 3:
        return 0.0, np.nan, "unfittable"
    scale = float(np.max(np.abs(P3D))) or 1.0
    if np.all(np.abs(P3D[m]) < 1e-12 * scale):
        return 0.0, np.nan, "zero"
    if (P3D[m] <= 0).any():
        return 0.0, np.nan, "unfittable"
    s = float(np.polyfit(np.log(k[m]), np.log(P3D[m]), 1)[0])
    if s >= -2.0:
        return 0.0, s, "divergent"
    return (float(k[-1] ** 2 * P3D[-1] / (2.0 * np.pi * abs(s + 2.0))),
            s, "ok")


def bounds_check(k3, r3, k1, r1):
    """
    The fit-free test described in the module docstring.

        min_{k>kp} r_3D  <=  r_1D(kp)  <=  max_{k>kp} r_3D

    Returns (lo, hi, inside) per k1, with NaN where the 3D band does not
    reach down to that k1 and the test cannot be applied.
    """
    lo = np.full(k1.size, np.nan)
    hi = np.full(k1.size, np.nan)
    for i, kp in enumerate(k1):
        m = k3 >= kp
        if not m.any() or kp < k3[0]:
            continue
        lo[i], hi[i] = float(r3[m].min()), float(r3[m].max())
    inside = (r1 >= lo) & (r1 <= hi)
    return lo, hi, inside


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
    ap.add_argument("--kcut-3d", type=float, default=None)
    ap.add_argument("--out", default="figures/drain")
    args = ap.parse_args()

    log = list(prov.header(inputs=[args.boot_npz, args.pk_ref, args.pk_test]))

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
    say("Stage 11 - the 3D gas difference, propagated into P1D")
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
    kfac = units.DEFAULT.k_kms_to_mpc(z)
    k1_mpc = kb * kfac
    say(f"k conversion at this z: k[s/km] x {kfac:.4f} = k[Mpc^-1]")
    say()

    # --- cut the 3D tables where the modes stop being trustworthy ---------
    k3_all = A["k"]
    knyq = A["_meta"].get("k_nyq")
    if args.kcut_3d:
        kcut, why = args.kcut_3d, "--kcut-3d"
    elif knyq:
        kcut = knyq
        why = (f"grid Nyquist, grid {A['_meta']['grid']}, "
               f"box {A['_meta']['box']:.4f}")
    else:
        kcut, why = k3_all[-1], "no grid in header, nothing cut"
    keep = k3_all <= kcut
    say(f"3D band cut at k = {kcut:.4f} Mpc^-1 [{why}]")
    say(f"  keeping {int(keep.sum())} of {k3_all.size} bins; "
        f"the table runs to {k3_all[-1]:.4f} Mpc^-1 "
        f"= {k3_all[-1]/max(kcut,1e-30):.3f} x the cut")
    if k3_all[-1] > 1.6 * kcut:
        say("  the bins above the cut are the CORNERS of the Fourier cube "
            "(out to sqrt(3) k_Nyq). They are real modes but poorly "
            "sampled and the MAS correction is unreliable there, so they "
            "are not used.")
    say()

    k3 = k3_all[keep]
    ratio_gas = (C["gas"] / A["gas"])[keep]
    g = 1.0 - ratio_gas
    say("measured gas-power suppression  g(k) = 1 - P_gas_test / P_gas_ref")
    say(f"{'k [1/Mpc]':>10s} {'P_gas ref':>12s} {'P_gas test':>12s} "
        f"{'ratio':>8s} {'g':>8s}")
    for jj in np.linspace(0, k3.size - 1, 12).astype(int):
        say(f"{k3[jj]:10.4f} {A['gas'][keep][jj]:12.5e} "
            f"{C['gas'][keep][jj]:12.5e} {ratio_gas[jj]:8.4f} {g[jj]:8.4f}")
    say(f"  ratio range over the kept band: {ratio_gas.min():.4f} "
        f"(at k = {k3[ratio_gas.argmin()]:.4f}) to {ratio_gas.max():.4f} "
        f"(at k = {k3[ratio_gas.argmax()]:.4f})")
    nsign = int((g > 0).sum())
    say(f"  g > 0 in {nsign} bins, g < 0 in {int((g < 0).sum())} bins")
    if nsign and nsign < k3.size:
        kx = k3[np.argmax(g < 0)] if (g < 0).any() else np.nan
        say(f"  g CHANGES SIGN near k = {kx:.3f} Mpc^-1: the test run has "
            f"less gas power at large scales and more at small scales, so "
            f"a single template mixes the two effects.")
    say()

    # --- the fit-free bound -----------------------------------------------
    R1 = B["ratio"]
    sR = B["err_ratio"] if "err_ratio" in B else np.full(R1.size, np.nan)
    lo, hi, inside = bounds_check(k3, ratio_gas, k1_mpc, R1)
    say("FIT-FREE BOUND.  P1D is a positive-kernel integral of P3D, so the")
    say("P1D ratio is a weighted average of the 3D ratio above that k and")
    say("must lie between its extremes there. If the flux were suppressed")
    say("by the same fraction as the gas, every row below would be INSIDE.")
    say("A violation only counts if it is larger than the error on r_1D,")
    say("so the distance outside the band is reported in sigma.")
    say(f"{'k [s/km]':>10s} {'k [1/Mpc]':>10s} {'r_1D':>8s} {'err':>7s} "
        f"{'min r_3D':>9s} {'max r_3D':>9s} {'n_sigma':>8s} {'verdict':>9s}")
    nsig = np.zeros(R1.size)
    n_out = 0
    for jj in np.flatnonzero(fm):
        if np.isnan(lo[jj]):
            v, ns = "no 3D", np.nan
        else:
            d = max(lo[jj] - R1[jj], R1[jj] - hi[jj], 0.0)
            ns = d / sR[jj] if sR[jj] > 0 else np.nan
            nsig[jj] = 0.0 if np.isnan(ns) else ns
            v = "inside" if d <= 0 else "OUTSIDE"
            if d > 0:
                n_out += 1
        say(f"{kb[jj]:10.5f} {k1_mpc[jj]:10.4f} {R1[jj]:8.4f} "
            f"{sR[jj]:7.4f} {lo[jj]:9.4f} {hi[jj]:9.4f} {ns:8.2f} "
            f"{v:>9s}")
    worst = float(np.nanmax(nsig[fm])) if fm.any() else 0.0
    if n_out == 0:
        say("  -> every testable bin is inside. The data are consistent "
            "with flux and gas changing by the same fraction; any fitted "
            "amplitude away from 1 is doing something else.")
    elif worst < 3.0:
        say(f"  -> {n_out} bin(s) outside, worst at {worst:.1f} sigma. "
            f"SUGGESTIVE, NOT PROOF. Below 3 sigma this is a hint that "
            f"flux and gas respond differently, not a demonstration.")
        say("     The fitted amplitude is the stronger statement here, "
            "because it uses every bin at once.")
    else:
        say(f"  -> {n_out} bin(s) outside, worst at {worst:.1f} sigma. "
            f"With no fitting and no free parameter, flux power is NOT "
            f"suppressed by the same fraction as gas power.")
    if np.isnan(lo[np.flatnonzero(fm)[0]]):
        j0 = np.flatnonzero(fm)[0]
        say(f"     NOTE: the lowest P1D bin (k = {k1_mpc[j0]:.4f} Mpc^-1) "
            f"carries the deepest deficit and CANNOT be tested, "
            f"because the 3D band starts at {k3[0]:.4f} Mpc^-1.")
        say("     A 3D measurement reaching lower k would make this test "
            "far sharper; that is a box-size question, not an analysis one.")
    say()

    # --- P1D -> P3D --------------------------------------------------------
    P3D, slope, sp = p3d_from_p1d(k1_mpc, Pa, args.spline_s)
    say(f"P1D -> P3D by differentiating a smoothing spline "
        f"(s = {args.spline_s})")
    say(f"  log-slope range: {slope.min():.3f} to {slope.max():.3f}")
    bad = slope >= 0
    if bad.any():
        say(f"  {int(bad.sum())} bins have slope >= 0, so P3D would be "
            f"negative there. They are at k = "
            + ", ".join(f"{v:.4f}" for v in k1_mpc[bad]) + " Mpc^-1.")
        good = np.flatnonzero(~bad)
        if good.size:
            # A negative P3D inside the integral is a corruption, not a
            # small error. Zeroing the bin is also wrong when the bin is an
            # endpoint carrying real signal, which is where a spline
            # overshoots. Borrowing the nearest fitted slope keeps the
            # magnitude sane and is stated rather than hidden.
            for j in np.flatnonzero(bad):
                jn = good[np.argmin(np.abs(good - j))]
                P3D[j] = -(2.0 * np.pi / k1_mpc[j] ** 2) * Pa[j] * slope[jn]
            say(f"    Those bins take the log-slope of their nearest "
                f"well-behaved neighbour instead. This is a patch on an "
                f"endpoint overshoot, not a measurement; if more than one "
                f"or two bins need it, raise --spline-s.")
        else:
            raise SystemExit("every bin has a non-negative slope; the "
                             "spline is not usable, raise --spline-s")
    t_ref, s_ref, st_ref = powerlaw_tail(k1_mpc, P3D)
    if args.no_extrapolate:
        t_ref = 0.0
    round_trip = p1d_from_p3d(k1_mpc, P3D, tail=t_ref)
    say(f"  round trip P1D -> P3D -> P1D recovers "
        f"{float(np.median((round_trip / np.maximum(Pa, 1e-300))[:-1])):.4f}"
        f" of P1D (median); tail status {st_ref}, slope {s_ref:.3f}")
    say()

    # --- the templates -----------------------------------------------------
    gi = np.interp(np.log(k1_mpc), np.log(k3), g, left=g[0], right=g[-1])
    n_ex = int((k1_mpc < k3[0]).sum() + (k1_mpc > k3[-1]).sum())
    if n_ex:
        say(f"  {n_ex} of {kb.size} P1D bins lie outside the kept 3D band "
            f"[{k3[0]:.3f}, {k3[-1]:.3f}] Mpc^-1; g is held flat there.")

    def template(weight, name):
        t, sl, st = powerlaw_tail(k1_mpc, P3D * weight)
        if args.no_extrapolate or st != "ok":
            if st == "divergent":
                say(f"  {name}: tail slope {sl:.3f} >= -2, the integral "
                    f"does not converge. Truncated; treat as a bound.")
            elif st == "zero":
                say(f"  {name}: the weight is zero over the top decade, so "
                    f"the tail is exactly zero. Nothing is being dropped.")
            elif st == "unfittable":
                say(f"  {name}: the integrand changes sign over the top "
                    f"decade, so no power-law tail exists. Truncated; the "
                    f"missing piece is bounded by the excess template.")
            t = 0.0
        else:
            say(f"  {name}: tail continued, slope {sl:.3f}, "
                f"contributing {t:.4e} km/s")
        return -(p1d_from_p3d(k1_mpc, P3D * weight, tail=0.0) + t)

    drain = template(gi, "combined g")
    t_sup = template(np.maximum(gi, 0.0), "suppression part (g > 0)")
    t_exc = template(np.minimum(gi, 0.0), "excess part (g < 0)")
    say()

    # --- fits ---------------------------------------------------------------
    nf = int(fm.sum())
    Cov = np.cov(boot[:, fm], rowvar=False)
    nbt = boot.shape[0]
    hl = (nbt - nf - 2) / (nbt - 1.0)
    if nbt < 3 * nf or hl <= 0:
        Cinv = np.diag(1.0 / sd[fm] ** 2)
        say(f"covariance: diagonal ({nbt} resamples is too few for {nf} "
            f"bins)")
    else:
        Cinv = hl * np.linalg.inv(Cov)
        say(f"covariance: full, Hartlap {hl:.4f}, {nbt} resamples")
    say()

    y = dP[fm]

    def fit(X, name, pnames):
        """
        Generalised least squares, guarded against a singular design.

        The split template turned out to be very nearly degenerate: the
        excess piece is almost constant across the fit window, so it and
        the free constant span nearly the same direction.  `inv` on that
        normal matrix returns negative variances and NaN errors, and the
        chi^2 it produces is meaningless.  `pinv` plus a reported
        condition number makes the degeneracy visible instead of letting
        it through as numbers.
        """
        N = X.T @ Cinv @ X
        cond = float(np.linalg.cond(N))
        cv = np.linalg.pinv(N, rcond=1e-12)
        bt = cv @ (X.T @ Cinv @ y)
        r = y - X @ bt
        c2 = float(r @ Cinv @ r)
        dof = nf - X.shape[1]
        say(name)
        var = np.diag(cv)
        for nm, v, e in zip(pnames, bt, var):
            err = np.sqrt(e) if e > 0 else np.nan
            say(f"    {nm:<8s} = {v:12.5e} +- {err:.3e}")
        say(f"    chi2/dof = {c2:.2f} / {dof} = {c2/dof:.3f}"
            f"    cond(N) = {cond:.2e}")
        if cond > 1e8 or (var <= 0).any():
            say("    WARNING: the design is degenerate at this precision. "
                "The parameters are not separately determined and the "
                "chi2 is not reliable. Read the next model instead.")
        say()
        return bt, cv, c2, r, dof, cond

    X2 = np.column_stack([-Pa[fm], np.ones(nf)])
    X3 = np.column_stack([drain[fm], np.ones(nf)])
    X4 = np.column_stack([t_sup[fm], t_exc[fm], np.ones(nf)])
    X5 = np.column_stack([t_sup[fm], t_exc[fm]])

    b2, c2v, chi2, r2, d2, q2 = fit(
        X2, "M2  DP1D = -f P1D_ref + c            (constant fraction)",
        ["f", "c"])
    b3, c3v, chi3, r3, d3, q3 = fit(
        X3, "M3  DP1D = alpha T[g] + c            (one combined template)",
        ["alpha", "c"])
    b4, c4v, chi4, r4, d4, q4 = fit(
        X4, "M4  DP1D = a_sup T[g+] + a_exc T[g-] + c   (split + constant)",
        ["a_sup", "a_exc", "c"])
    b5, c5v, chi5, r5, d5, q5 = fit(
        X5, "M5  DP1D = a_sup T[g+] + a_exc T[g-]       (split, NO constant)",
        ["a_sup", "a_exc"])

    # M3's model space is a strict subspace of M4's, because
    # T[g] = T[g+] + T[g-] exactly. chi2 therefore CANNOT rise from M3 to
    # M4. If it does, the fit is numerically broken, not physically
    # interesting, and saying so is the whole point of checking.
    say("nesting check: span{T[g], 1} is contained in "
        "span{T[g+], T[g-], 1},")
    say("so chi2(M4) <= chi2(M3) must hold as algebra.")
    if chi4 > chi3 + 1e-6:
        say(f"  VIOLATED: chi2(M4) = {chi4:.2f} > chi2(M3) = {chi3:.2f}. "
            f"The M4 fit is numerically broken; ignore its numbers and "
            f"use M5, which drops the redundant constant.")
    else:
        say(f"  holds: chi2(M4) = {chi4:.2f} <= chi2(M3) = {chi3:.2f}")
    say()

    say(f"M2 -> M3: delta chi2 = {chi2-chi3:8.2f} at equal parameters")
    say(f"M3 -> M4: delta chi2 = {chi3-chi4:8.2f} for 1 extra parameter")
    say(f"M3 -> M5: delta chi2 = {chi3-chi5:8.2f} at equal parameters")
    say()
    # Which split fit is usable is a property of the design, not a
    # preference: pick the conditioned one and say which was picked.
    ok4 = q4 < 1e8 and (np.diag(c4v) > 0).all() and chi4 <= chi3 + 1e-6
    bs, cs, chis, ds, rs, tag = ((b4, c4v, chi4, d4, r4, "M4")
                                 if ok4 else
                                 (b5, c5v, chi5, d5, r5, "M5"))
    Xs = X4 if ok4 else X5
    say(f"split-template result taken from {tag} "
        f"({'well conditioned' if ok4 else 'M4 was degenerate'})")
    say(f"  a_sup = {bs[0]:.4f} +- {np.sqrt(max(cs[0,0],0)):.4f}   "
        f"(flux response to the large-scale gas deficit)")
    say(f"  a_exc = {bs[1]:.4f} +- {np.sqrt(max(cs[1,1],0)):.4f}   "
        f"(flux response to the small-scale gas excess)")
    say("  These are separately meaningful only because the template was "
        "split: with one template the two responses are forced equal and "
        "the single amplitude absorbs both.")
    say()

    if ok4:
        csig = abs(bs[2]) / np.sqrt(cs[2, 2])
        say(f"M4's constant is {csig:.1f} sigma from zero.")
        if csig < 2:
            say("  -> consistent with zero. The two MEASURED 3D components "
                "account for the whole P1D difference; nothing else is "
                "needed.")
        else:
            say("  -> not zero. Something in DP1D is not carried by the 3D "
                "gas power difference, and its size is that constant.")
    else:
        say("The 'is a constant needed' question is answered instead by "
            "M5 against M3, since M5 has no constant at all:")
        say(f"  chi2(M5) = {chi5:.2f} / {d5}   vs   "
            f"chi2(M3) = {chi3:.2f} / {d3}")
        if chi5 <= chi3 + 4:
            say("  -> dropping the constant costs nothing. The two measured "
                "3D components account for the P1D difference on their "
                "own, with no free additive term.")
        else:
            say("  -> dropping the constant hurts. Something in DP1D is not "
                "carried by the 3D gas power difference.")
    say()

    say(f"{'k [s/km]':>10s} {'DP1D':>12s} {'err':>10s} {'M2':>12s} "
        f"{'M3':>12s} {tag:>12s} {'(d-' + tag + ')/e':>10s}")
    for jx, jj in enumerate(np.flatnonzero(fm)):
        say(f"{kb[jj]:10.5f} {dP[jj]:12.5e} {sd[jj]:10.3e} "
            f"{(X2@b2)[jx]:12.5e} {(X3@b3)[jx]:12.5e} {(Xs@bs)[jx]:12.5e} "
            f"{rs[jx]/sd[jj]:10.2f}")
    say()
    for nm, r in (("M2", r2), ("M3", r3), (tag, rs)):
        say(f"residual signs {nm}: "
            + "".join("+" if v > 0 else "-" for v in r))
    say("  a run of like signs means a missing shape, not scatter, whatever")
    say("  the chi2 says. Compare the runs, not just the numbers.")
    say()

    # --- figure ------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        npan = 3 if args.diagnose else 2
        fig, axes = plt.subplots(npan, 1, figsize=(7.8, 3.3 * npan + 1.4))
        ax = axes[0]
        ax.axhline(0, color="k", lw=1, ls="--")
        ax.errorbar(kb, dP, yerr=sd, fmt="o", ms=4, color="#1f77b4", lw=1.2,
                    capsize=2, label=r"$\Delta P_{\rm 1D}$")
        kk = kb[fm]
        ax.plot(kk, X2 @ b2, color="#d62728", lw=1.5,
                label=rf"M2 const $f$, $\chi^2/\nu$={chi2/d2:.2f}")
        ax.plot(kk, X3 @ b3, color="#ff7f0e", lw=1.5,
                label=rf"M3 combined, $\chi^2/\nu$={chi3/d3:.2f}")
        ax.plot(kk, Xs @ bs, color="#2ca02c", lw=2.2,
                label=rf"{tag} split, $\chi^2/\nu$={chis/ds:.2f}")
        ax.set(xscale="log", ylabel=r"$\Delta P_{\rm 1D}$ [km s$^{-1}$]",
               title=(rf"$z={z:.1f}$,  $a_{{\rm sup}}={bs[0]:.2f}$,  "
                      rf"$a_{{\rm exc}}={bs[1]:.2f}$"))
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.2, which="both")

        bx = axes[1]
        bx.axhline(1, color="k", lw=1, ls="--")
        bx.plot(k3, ratio_gas, color="#9467bd", lw=2,
                label=r"$P_{\rm gas}$ test/ref (3D, stage 09)")
        bx.plot(k1_mpc, R1, color="#1f77b4", lw=2,
                label=r"$P_{\rm 1D}$ test/ref")
        okb = ~np.isnan(lo)
        bx.fill_between(k1_mpc[okb], lo[okb], hi[okb], color="#9467bd",
                        alpha=0.15, lw=0,
                        label=r"allowed band for $r_{\rm 1D}$")
        if knyq:
            bx.axvline(knyq, color="0.5", ls=":", lw=1.2)
            bx.text(knyq, bx.get_ylim()[1], r" $k_{\rm Nyq}$", fontsize=8,
                    color="0.4", va="top")
        bx.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]", ylabel="ratio")
        bx.legend(frameon=False, fontsize=8)
        bx.grid(alpha=0.2, which="both")

        if args.diagnose:
            # The round trip plotted as P1D against P1D over eight decades
            # of log axis looks perfect no matter what, because the whole
            # result lives at the percent level and a log axis that wide
            # cannot show a percent. The ratio on a linear axis can.
            cx = axes[2]
            cx.axhline(1.0, color="k", lw=1, ls="--")
            cx.axhspan(0.95, 1.05, color="0.6", alpha=0.15, lw=0)
            cx.plot(k1_mpc, round_trip / np.maximum(Pa, 1e-300), "o-",
                    ms=4, lw=1.4, color="#ff7f0e")
            cx.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]",
                   ylabel="round trip / data",
                   title="closure of the P1D -> P3D -> P1D inversion "
                         "(shaded: $\\pm5\\%$)")
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
