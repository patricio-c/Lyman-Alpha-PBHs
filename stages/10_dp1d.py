#!/usr/bin/env python3
"""
Stage 10 - the additive difference DP1D, its error bars, and the fits that
decide the shape of the result.

Everything upstream of this stage reports ratios.  A ratio is the natural
thing to plot and the wrong thing to fit, because the model that produces
the signal is additive:

    P1D_FCT(k)  =  P1D_CDM(k)  +  DP1D(k)

and for a discrete population the Poisson isocurvature term is white noise
in 3D, so to O((k/k_max)^2) the prediction is that DP1D is a CONSTANT:

    DP1D(k) = (A_P / 4 pi) (k_max^2 - k^2)  ~  const

That is a falsifiable statement about the additive difference, and it is
invisible in the ratio panel, where the same constant divided by a steeply
falling P1D_CDM(k) looks like a tilt.

Three models are fitted and compared with a proper chi^2:

    M0   DP1D = c                      the white-noise prediction
    M1   DP1D = a (k_max^2 - k^2)      the same, with the k^2 correction
                                       (needs --kmax-mpc; skipped otherwise)
    M2   DP1D = -f P1D_CDM + c         drain + boost, two regimes

M2 is the one worth explaining.  If the star-formation scheme removes a
fraction of the absorbing gas in FCT that it does not remove in CDM, the
effect is multiplicative on the flux power, so it enters as -f P1D_CDM;
the primordial small-scale power enters additively as +c.  Because
P1D_CDM falls steeply with k, the first term dominates at low k and the
second at high k, and the difference CHANGES SIGN at P1D_CDM = c/f.  A
single constant c cannot change sign, so M0 and M2 are distinguishable by
exactly the feature already seen in the ratio.  This stage measures which
one the data prefer, with a number rather than an impression.

Error bars
----------
Bootstrap over sightlines, never analytic.  Two things make it honest:

*   The resampling is PAIRED.  Both runs were extracted on the same
    sightline positions from boxes that share Panphasia phases, so ray i
    of FCT and ray i of CDM sample the same comoving structure.  Drawing
    the same index set from both runs keeps that correlation, and the
    error on the DIFFERENCE is far smaller than the errors on the two
    curves added in quadrature.  `--unpaired` recomputes without it, and
    the ratio of the two is printed, because if pairing buys nothing then
    the sightlines are not actually matched and something upstream is
    wrong.

*   The chi^2 uses the full bootstrap COVARIANCE across k bins, with the
    Hartlap correction to the inverse.  Neighbouring log-k bins are
    correlated; a diagonal chi^2 would overstate the significance of every
    fit.  If there are too few resamples to estimate the covariance
    (n_boot < 3 n_bins) the code falls back to diagonal and says so.

Approximation, and its MEASURED FAILURE.  By default the flux rescaling A
and the global mean flux are held at their full-sample values inside the
bootstrap rather than re-solved per resample.  This file originally
claimed A would move by ~0.2%, far below the scatter of P1D itself, so
that the neglected term was invisible.  `--resample-A` was written to
check that claim and the claim did not survive it: on the 40 Mpc/h pair
A moves by 1.15% (CDM) and 0.87% (FCT) against a P1D bootstrap scatter of
1.27%.  Comparable, not negligible.  The frozen-A errors are therefore
UNDERESTIMATED by a factor that has to be measured, not assumed.

`--exact-boot N` does the honest thing: it re-solves A and recomputes the
per-sightline spectra inside each of N resamples, reports the inflation
factor against the frozen-A errors, and uses the exact errors for
everything downstream -- the table, the covariance, the fits and the
figure.  It costs about a second per resample, so a few hundred is the
practical range; that is enough to measure the inflation even when it is
too few for a full covariance, in which case the code falls back to
diagonal and says so.

**No chi^2 from this stage should be quoted outside this repository unless
it came from an --exact-boot run.**  The stage prints a warning at the end
when it did not.

Precondition
------------
Both runs MUST be renormalised to the same tau_eff before the difference
means anything, otherwise DP1D measures the difference in mean flux.  This
is enforced: `--tau-eff` or the default target is printed and written to
the log together with the A actually used for each run.

Usage
-----
    python stages/10_dp1d.py cache/cache_cdm40.npz cache/cache_fct40.npz \\
        --labels CDM FCT --tau-eff 0.42461 --nboot 2000 \\
        --out figures/dp1d_40

Options
-------
    cache_ref cache_test   reference first; DP1D = test - ref
    --labels A B           display names
    --tau-eff X            common target (default: Turner+24 at the ref z)
    --nboot N              bootstrap resamples (default 2000)
    --nbins N              log-k bins (default 20)
    --kmax-frac F          cut at F x Nyquist (default 0.5)
    --window desi|all      k range used for the FITS (default desi)
    --kmax-mpc X           filtering-scale cutoff [Mpc^-1], enables M1
    --unpaired             also report unpaired errors
    --resample-A N         re-solve A on N resamples as a check only
    --exact-boot N         re-solve A INSIDE N resamples and use those
                           errors for the table, covariance, fits and
                           figure. Required for any quotable chi^2.
    --seed S               bootstrap seed (default 20260908)
    --out PREFIX           writes PREFIX.png, PREFIX.txt and PREFIX_boot.npz
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import cache as cachelib, units  # noqa: E402
from common.p1d import flux, solve_A, tau_eff  # noqa: E402


# --- estimator, kept per sightline so the bootstrap is cheap ---------------

def per_los_pk(F, dv):
    """
    P1D of every sightline separately, using the GLOBAL mean flux.

    The mean over rows of what this returns is exactly common.p1d's
    estimator; keeping the rows lets the bootstrap be a re-average instead
    of a re-FFT, which is what makes 2000 resamples take seconds.
    """
    F = np.atleast_2d(np.asarray(F, dtype=np.float64))
    n_los, n_pix = F.shape
    d = F / F.mean() - 1.0
    fk = np.fft.rfft(d, axis=-1) * dv
    Pk = (fk.real ** 2 + fk.imag ** 2) / (n_pix * dv)
    return 2.0 * np.pi * np.fft.rfftfreq(n_pix, d=dv), Pk


def bin_rows(k, Pmat, nbins, kmin, kmax):
    """Log-bin along k, preserving the sightline axis."""
    m = (k > 0) & (k >= kmin) & (k <= kmax)
    if m.sum() < nbins:
        raise SystemExit(f"only {m.sum()} modes in [{kmin:.4g}, {kmax:.4g}] "
                         f"s/km, cannot make {nbins} bins")
    kk, Pm = k[m], Pmat[:, m]
    edges = np.geomspace(kk.min() * 0.999, kk.max() * 1.001, nbins + 1)
    idx = np.clip(np.digitize(kk, edges) - 1, 0, nbins - 1)
    cnt = np.bincount(idx, minlength=nbins)
    ok = cnt > 0
    kb = np.bincount(idx, weights=kk, minlength=nbins)[ok] / cnt[ok]
    B = np.empty((Pm.shape[0], int(ok.sum())))
    for j, jj in enumerate(np.flatnonzero(ok)):
        B[:, j] = Pm[:, idx == jj].mean(axis=1)
    return kb, B, cnt[ok]


# --- generalised least squares with a bootstrap covariance -----------------

def gls(X, y, Cinv):
    """
    beta, cov(beta), chi2 for a linear model X beta = y.

    Written out rather than pulled from a fitting package because the
    covariance is the whole point: with Cinv diagonal these numbers are
    optimistic by a factor of a few, and the reader has to be able to see
    that Cinv is what it claims to be.
    """
    XtC = X.T @ Cinv
    cov = np.linalg.inv(XtC @ X)
    beta = cov @ (XtC @ y)
    r = y - X @ beta
    return beta, cov, float(r @ Cinv @ r)


def hartlap(nboot, nbins):
    """
    Unbiasing factor for the inverse of a covariance estimated from nboot
    realisations (Hartlap, Simon & Schneider 2007).  Without it the inverse
    is biased high and every chi^2 comes out too good.
    """
    num = nboot - nbins - 2
    if num <= 0:
        return None
    return num / (nboot - 1.0)


def solve_A_near(tau, target, A0):
    """
    solve_A bracketed tightly around a known solution.

    Inside a bootstrap the resample's A sits within a percent or so of the
    full-sample value, so a tight bracket cuts brentq from ~40 evaluations
    to ~15 and each evaluation is an exp() over three million pixels.  The
    fallback exists because a pathological resample can push the solution
    outside the bracket, and a slow correct answer beats a fast crash.
    """
    try:
        return solve_A(tau, target, lo=0.6 * A0, hi=1.7 * A0)
    except ValueError:
        return solve_A(tau, target)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cache_ref")
    ap.add_argument("cache_test")
    ap.add_argument("--labels", nargs=2, default=None)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--kmax-frac", type=float, default=0.5)
    ap.add_argument("--window", default="desi", choices=["desi", "all"])
    ap.add_argument("--kmax-mpc", type=float, default=None)
    ap.add_argument("--unpaired", action="store_true")
    ap.add_argument("--resample-A", type=int, default=0)
    ap.add_argument("--exact-boot", type=int, default=0)
    ap.add_argument("--seed", type=int, default=20260908)
    ap.add_argument("--out", default="figures/dp1d")
    args = ap.parse_args()

    lab_a, lab_b = args.labels or [
        os.path.basename(p).replace("cache_", "").replace(".npz", "")
        for p in (args.cache_ref, args.cache_test)]

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    # load_pair refuses mismatched z, dv or shape, which is exactly the
    # precondition an additive difference needs.
    a, b = cachelib.load_pair(args.cache_ref, args.cache_test)
    n_los, n_pix = a.tau.shape
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(a.z)

    say("=" * 74)
    say("Stage 10 - additive difference, bootstrap errors, and model fits")
    say("=" * 74)
    say(f"reference : {lab_a}   {args.cache_ref}")
    say(f"test      : {lab_b}   {args.cache_test}")
    say(f"z = {a.z:.4f}   dv = {a.dv:.5f} km/s   "
        f"n_los = {n_los}   n_pix = {n_pix}")
    say(f"common tau_eff target = {target:.5f}")

    A_a, A_b = solve_A(a.tau, target), solve_A(b.tau, target)
    say(f"{lab_a:>10s}  A = {A_a:.6f}   tau_eff out = "
        f"{tau_eff(a.tau, A_a):.6f}")
    say(f"{lab_b:>10s}  A = {A_b:.6f}   tau_eff out = "
        f"{tau_eff(b.tau, A_b):.6f}")
    say()

    k, Pa_rows = per_los_pk(flux(a.tau, A_a), a.dv)
    _, Pb_rows = per_los_pk(flux(b.tau, A_b), b.dv)

    k_nyq = args.kmax_frac * np.pi / a.dv
    kb, Ba, cnt = bin_rows(k, Pa_rows, args.nbins, k[1], k_nyq)
    _, Bb, _ = bin_rows(k, Pb_rows, args.nbins, k[1], k_nyq)
    nb = kb.size

    # --- bootstrap ---------------------------------------------------------
    rng = np.random.default_rng(args.seed)
    boot_d = np.empty((args.nboot, nb))
    boot_r = np.empty((args.nboot, nb))
    boot_a = np.empty((args.nboot, nb))
    for i in range(args.nboot):
        idx = rng.integers(0, n_los, n_los)
        pa = Ba[idx].mean(axis=0)
        pb = Bb[idx].mean(axis=0)
        boot_a[i] = pa
        boot_d[i] = pb - pa
        boot_r[i] = pb / pa

    Pa = Ba.mean(axis=0)
    Pb = Bb.mean(axis=0)
    dP = Pb - Pa
    R = Pb / Pa
    sd = boot_d.std(axis=0, ddof=1)
    sr = boot_r.std(axis=0, ddof=1)
    sa = boot_a.std(axis=0, ddof=1)

    say(f"bootstrap: {args.nboot} paired resamples of {n_los} sightlines, "
        f"seed {args.seed}")

    if args.unpaired:
        bu = np.empty((args.nboot, nb))
        for i in range(args.nboot):
            ia = rng.integers(0, n_los, n_los)
            ib = rng.integers(0, n_los, n_los)
            bu[i] = Bb[ib].mean(axis=0) - Ba[ia].mean(axis=0)
        su = bu.std(axis=0, ddof=1)
        gain = float(np.median(su / np.maximum(sd, 1e-300)))
        say(f"pairing gain: unpaired errors are {gain:.2f}x larger "
            f"(median over bins)")
        if gain < 1.05:
            say("  WARNING: pairing buys nothing. Either the two extractions "
                "do NOT share sightline positions, or the phases are not "
                "matched. Check before trusting any error bar below.")
    say()

    # --- optional check that freezing A is harmless ------------------------
    if args.resample_A:
        say(f"check: re-solving A on {args.resample_A} resamples")
        As = []
        for _ in range(args.resample_A):
            idx = rng.integers(0, n_los, n_los)
            As.append((solve_A(a.tau[idx], target),
                       solve_A(b.tau[idx], target)))
        As = np.asarray(As)
        say(f"  {lab_a}: A = {A_a:.6f} full sample, "
            f"resample scatter {As[:, 0].std(ddof=1) / A_a * 100:.3f}%")
        say(f"  {lab_b}: A = {A_b:.6f} full sample, "
            f"resample scatter {As[:, 1].std(ddof=1) / A_b * 100:.3f}%")
        say(f"  P1D bootstrap scatter for comparison: "
            f"{float(np.median(sa / Pa)) * 100:.3f}%")
        say("  the first two must be well below the third for the frozen-A "
            "approximation to hold")
        say()

    # --- the honest bootstrap, if asked for ---------------------------------
    exact = False
    if args.exact_boot:
        say(f"exact bootstrap: {args.exact_boot} resamples with A re-solved "
            f"and the spectra recomputed")
        be = np.empty((args.exact_boot, nb))
        br = np.empty((args.exact_boot, nb))
        for i in range(args.exact_boot):
            idx = rng.integers(0, n_los, n_los)
            ta, tb = a.tau[idx], b.tau[idx]
            _, ra = per_los_pk(flux(ta, solve_A_near(ta, target, A_a)), a.dv)
            _, rb = per_los_pk(flux(tb, solve_A_near(tb, target, A_b)), b.dv)
            _, ea, _ = bin_rows(k, ra, args.nbins, k[1], k_nyq)
            _, eb, _ = bin_rows(k, rb, args.nbins, k[1], k_nyq)
            pa, pb = ea.mean(axis=0), eb.mean(axis=0)
            be[i] = pb - pa
            br[i] = pb / pa
        se = be.std(axis=0, ddof=1)
        infl = se / np.maximum(sd, 1e-300)
        say(f"  error inflation over the frozen-A bootstrap: "
            f"median {float(np.median(infl)):.2f}x, "
            f"range {float(infl.min()):.2f} - {float(infl.max()):.2f}")
        say("  everything below uses these errors, not the frozen-A ones.")
        boot_d, sd, sr = be, se, br.std(axis=0, ddof=1)
        exact = True
        say()

    # --- the table ---------------------------------------------------------
    say(f"{'k [s/km]':>10s} {'k [1/Mpc]':>10s} {'P_' + lab_a:>12s} "
        f"{'DP1D':>12s} {'err':>10s} {'DP/err':>8s} "
        f"{'ratio':>8s} {'err':>8s}")
    kfac = units.DEFAULT.k_kms_to_mpc(a.z)
    for j in range(nb):
        say(f"{kb[j]:10.5f} {kb[j]*kfac:10.4f} {Pa[j]:12.5e} "
            f"{dP[j]:12.5e} {sd[j]:10.3e} {dP[j]/sd[j]:8.2f} "
            f"{R[j]:8.4f} {sr[j]:8.4f}")
    say()

    # --- where the ratio crosses one --------------------------------------
    cross = None
    s = np.sign(R - 1.0)
    w = np.flatnonzero(s[:-1] * s[1:] < 0)
    if w.size:
        j = w[0]
        t = (1.0 - R[j]) / (R[j + 1] - R[j])
        cross = float(np.exp(np.log(kb[j]) + t * (np.log(kb[j + 1])
                                                  - np.log(kb[j]))))
        say(f"ratio crosses 1 at k = {cross:.5f} s/km "
            f"= {cross*kfac:.3f} Mpc^-1")
    else:
        say("ratio does not cross 1 inside the measured band")
    say()

    # --- fits --------------------------------------------------------------
    lo, hi = units.desi_window(a.z) if args.window == "desi" \
        else (kb.min() * 0.999, kb.max() * 1.001)
    fm = (kb >= lo) & (kb <= hi)
    nf = int(fm.sum())
    say(f"fit window ({args.window}): {lo:.5f} - {hi:.5f} s/km, "
        f"{nf} of {nb} bins")
    if nf < 4:
        raise SystemExit("fewer than 4 bins in the fit window; widen it")

    nb_used = boot_d.shape[0]
    C = np.cov(boot_d[:, fm], rowvar=False)
    h = hartlap(nb_used, nf)
    if h is None or nb_used < 3 * nf:
        say(f"  n_boot = {nb_used} is too few for a {nf}x{nf} covariance; "
            f"falling back to DIAGONAL errors. chi2 below is optimistic.")
        Cinv = np.diag(1.0 / sd[fm] ** 2)
        used_cov = "diagonal"
    else:
        Cinv = h * np.linalg.inv(C)
        used_cov = f"full, Hartlap factor {h:.4f}"
    say(f"  covariance: {used_cov}   "
        f"errors: {'exact' if exact else 'frozen A (provisional)'}")
    dc = np.diag(C)
    off = C / np.sqrt(np.outer(dc, dc))
    say(f"  median |correlation| between neighbouring bins: "
        f"{float(np.median(np.abs(np.diag(off, 1)))):.3f}")
    say()

    y = dP[fm]
    results = {}

    beta, cov, chi2 = gls(np.ones((nf, 1)), y, Cinv)
    results["M0"] = (chi2, nf - 1)
    say(f"M0  DP1D = c")
    say(f"    c        = {beta[0]:.5e} +- {np.sqrt(cov[0,0]):.3e} km/s")
    say(f"    chi2/dof = {chi2:.2f} / {nf-1} = {chi2/(nf-1):.3f}")
    say()

    if args.kmax_mpc:
        kmax_kms = args.kmax_mpc / kfac
        X = ((kmax_kms ** 2 - kb[fm] ** 2)).reshape(-1, 1)
        beta1, cov1, chi21 = gls(X, y, Cinv)
        results["M1"] = (chi21, nf - 1)
        say(f"M1  DP1D = a (k_max^2 - k^2),  k_max = {args.kmax_mpc:.3f} "
            f"Mpc^-1 = {kmax_kms:.4f} s/km")
        say(f"    a        = {beta1[0]:.5e} +- {np.sqrt(cov1[0,0]):.3e}")
        say(f"    chi2/dof = {chi21:.2f} / {nf-1} = {chi21/(nf-1):.3f}")
        say(f"    (k/k_max)^2 at the top of the window = "
            f"{(kb[fm].max()/kmax_kms)**2:.2e}")
        say()
    else:
        say("M1  skipped: needs --kmax-mpc from the pressure filtering "
            "scale.")
        say("    k_max is physics, not a knob; it is not invented here.")
        say("    Note that for k << k_max, M1 is M0 with a different name.")
        say()

    X2 = np.column_stack([-Pa[fm], np.ones(nf)])
    beta2, cov2, chi22 = gls(X2, y, Cinv)
    results["M2"] = (chi22, nf - 2)
    f_hat, c_hat = beta2
    say(f"M2  DP1D = -f P1D_{lab_a} + c")
    say(f"    f        = {f_hat:.5f} +- {np.sqrt(cov2[0,0]):.5f}   "
        f"(drain fraction)")
    say(f"    c        = {c_hat:.5e} +- {np.sqrt(cov2[1,1]):.3e} km/s   "
        f"(boost)")
    say(f"    chi2/dof = {chi22:.2f} / {nf-2} = {chi22/(nf-2):.3f}")
    if f_hat > 0 and c_hat > 0:
        pcross = c_hat / f_hat
        say(f"    predicted sign change where P1D_{lab_a} = {pcross:.5e}")
        if Pa[fm].min() <= pcross <= Pa[fm].max():
            kpred = float(np.exp(np.interp(np.log(pcross),
                                           np.log(Pa[fm][::-1]),
                                           np.log(kb[fm][::-1]))))
            say(f"    i.e. at k = {kpred:.5f} s/km")
            if cross:
                verdict = ("consistent" if abs(np.log(kpred / cross)) < 0.15
                           else "INCONSISTENT")
                say(f"    measured crossing      k = {cross:.5f} s/km"
                    f"   -> {verdict}")
        else:
            say("    that P1D value is outside the fitted band")
    say()

    d01 = results["M0"][0] - results["M2"][0]
    say(f"M0 vs M2: delta chi2 = {d01:.2f} for 1 extra parameter")
    if d01 > 9:
        say("  -> the two-regime model is preferred at >3 sigma. The "
            "difference is NOT a single constant, so the white-noise")
        say("     prediction alone does not describe this pair.")
    elif d01 < 1:
        say("  -> no preference. DP1D is consistent with a constant, "
            "which is the white-noise prediction.")
    else:
        say("  -> mild preference for M2. Not decisive; more sightlines "
            "or more redshifts are needed.")
    say()

    # --- effective slope of the reference -----------------------------------
    p = np.polyfit(np.log(kb[fm]), np.log(Pa[fm]), 1)
    say(f"effective slope of P1D_{lab_a} in the window: "
        f"n_eff = {p[0]:.4f}")
    say(f"  implied 3D slope n_3D = n_eff - 2 = {p[0]-2:.4f}")
    if p[0] - 2 >= -2:
        say("  WARNING: n_3D >= -2, the P1D integral does not converge as a "
            "power law here; the logarithmic case applies.")
    else:
        say(f"  the ratio denominator therefore behaves as k^{p[0]:.2f}")
    say()

    if not exact:
        say("WARNING: these errors froze A at its full-sample value, an "
            "approximation this stage has already measured to be invalid "
            "(see --resample-A).")
        say("         Every chi2 above is therefore provisional. Re-run with "
            "--exact-boot 300 before quoting any of them.")
        say()

    # --- figure -------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax, bx) = plt.subplots(
            2, 1, figsize=(7.6, 8.4), sharex=True,
            gridspec_kw={"height_ratios": [1.5, 1], "hspace": 0.07})

        ax.axhline(0, color="k", lw=1, ls="--")
        ax.errorbar(kb, dP, yerr=sd, fmt="o", ms=4, color="#1f77b4",
                    lw=1.2, capsize=2,
                    label=f"$\\Delta P_{{\\rm 1D}}$ ({lab_b} $-$ {lab_a})")
        kk = kb[fm]
        ax.plot(kk, np.full(nf, beta[0]), color="#d62728", lw=1.8,
                label=f"M0 const, $\\chi^2/\\nu$={chi2/(nf-1):.2f}")
        ax.plot(kk, X2 @ beta2, color="#2ca02c", lw=1.8,
                label=f"M2 drain+boost, $\\chi^2/\\nu$={chi22/(nf-2):.2f}")
        ax.axvspan(*units.desi_window(a.z), color="green", alpha=0.07,
                   zorder=0)
        ax.set(xscale="log", ylabel=r"$\Delta P_{\rm 1D}$  [km s$^{-1}$]",
               title=(f"$z={a.z:.1f}$, "
                      rf"$\tau_{{\rm eff}}={target:.5f}$, "
                      f"{boot_d.shape[0]} "
                      f"{'exact' if exact else 'frozen-A'} bootstraps"))
        ax.legend(frameon=False, fontsize=9, loc="best")
        ax.grid(alpha=0.2, which="both")

        bx.axhline(1, color="k", lw=1, ls="--")
        bx.fill_between(kb, R - sr, R + sr, color="#1f77b4", alpha=0.25, lw=0)
        bx.plot(kb, R, color="#1f77b4", lw=2)
        if cross:
            bx.axvline(cross, color="0.4", ls=":", lw=1.2)
            bx.text(cross, bx.get_ylim()[1], f" {cross:.4f}", fontsize=8,
                    color="0.3", va="top")
        bx.axvspan(*units.desi_window(a.z), color="green", alpha=0.07,
                   zorder=0)
        bx.set(xscale="log", xlabel=r"$k$  [s km$^{-1}$]",
               ylabel=f"{lab_b} / {lab_a}")
        bx.grid(alpha=0.2, which="both")

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

    np.savez_compressed(
        args.out + "_boot.npz", k=kb, P_ref=Pa, P_test=Pb, dP=dP, ratio=R,
        err_dP=sd, err_ratio=sr, boot_dP=boot_d.astype(np.float32),
        fit_mask=fm, z=a.z, tau_eff_target=target, A_ref=A_a, A_test=A_b,
        exact_bootstrap=exact, n_boot=boot_d.shape[0])
    print(f"written -> {args.out}_boot.npz")


if __name__ == "__main__":
    main()
