#!/usr/bin/env python3
"""
Test 9 - the significance of the P1D suppression, without a pairing.

Stage 02 answered open question #3 on 2026-08-31: the sightlines of the two
runs are NOT the same lines.  0% match by index, no valid permutation, a
median transverse offset of about a third of the box.  The two runs drew
independent random positions.  The paired bootstrap that produced the 20
sigma quoted at COSMO-26 therefore paired unrelated lines, and the number
has to be re-derived from something that does not assume a pairing.

What this script is actually testing
------------------------------------
The expectation written into HANDOFF section 0 is that the unpaired number
comes out weaker.  That does not follow on its own.  A bootstrap *measures*
the covariance between the two samples, it does not assume one: resampling
a common index set and taking the difference gives a variance of
var_A + var_B - 2cov, and when the lines are unrelated the cov in the data
is already zero.  A paired bootstrap over a fake pairing should therefore
degrade to the unpaired answer by itself, and the 20 sigma may be close to
unchanged.  That is a claim about the data, so this script measures it
instead of arguing it: it runs both and prints them side by side.

Four error bars, in increasing order of how much a referee should believe
them:

    analytic          std over sightlines / sqrt(n_los), what
                      common/p1d.py returns.  Assumes independent lines.
    bootstrap, common index
                      resample one index set, apply it to both runs.  This
                      is what the 20 sigma was, as far as the pairing is
                      concerned.
    bootstrap, unpaired
                      resample each run independently.  The honest version
                      given stage 02.
    block jackknife   delete one spatial tile of the transverse plane from
                      BOTH runs at once.  Quote this one.

Why the block jackknife is the one to quote
-------------------------------------------
6144 sightlines through a 40 Mpc/h box sit about 0.5 Mpc/h apart, far below
the correlation length of the field.  They are not 6144 independent
measurements, and every estimator above that resamples individual lines
inherits that assumption.  The block jackknife does not: it resamples
regions.  The script reports the effective number of independent lines
implied by the ratio of the two, which is the number worth putting in the
paper.

The tiling also recovers part of what the broken pairing cost.  Line i of
CDM is not line i of FCT, but tile b of CDM and tile b of FCT are the same
region of the same initial conditions.  The large-scale variance that a
pairing was supposed to cancel is a property of the region, not of the
line, so deleting the same tile from both runs cancels it without any
sightline ever being matched.

Usage
-----
    python tests/t9_unpaired_significance.py \\
        --cdm cache/cache_cdm.npz --fct cache/cache_fct.npz \\
        --los-a /path/cdm/los_0010.hdf5 --los-b /path/fct/los_0010.hdf5 \\
        --out figures/t9

Without --los-a/--los-b the bootstraps still run and the block jackknife is
skipped with a warning.  It is the part worth having, so pass them.

Options
-------
    --cdm / --fct PATH    the two caches
    --los-a / --los-b     the SWIFT LOS files, for the sightline positions
    --pos-a / --pos-b     [n_los, 2] .npy of transverse positions instead,
                          with --box-int
    --box-int X           box size in the same internal units as --pos-*
    --band LO HI          the k band the quoted number refers to, s/km
                          (default 0 to 0.01, the low-k definition used by
                          t8).  If the original 20 sigma used a different
                          statistic, match it here before comparing.
    --nside N             tiles per side, N^2 blocks total (default 8)
    --n-boot N            bootstrap draws (default 4000)
    --nbins N             log bins for the curve (default 20)
    --tau-eff X           common rescaling target (default Turner+24 at z)
    --no-rescale          leave A = 1 in both runs
    --seed N              rng seed (default 1)
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2 as chi2_dist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from common import boot, cache as cachelib, units  # noqa: E402
from common.p1d import flux, nyquist_cut, solve_A, tau_eff  # noqa: E402

C_CDM, C_FCT = "#1f77b4", "#d62728"


def load_stage02():
    """
    Import stages/02_check_los_match.py by path.

    The module name starts with a digit so it cannot be imported normally.
    Loading it by path rather than copying `ray_positions` here is
    deliberate: the ray position has to be read with exactly the logic the
    extractor uses, and there must be one copy of that logic.
    """
    p = os.path.join(ROOT, "stages", "02_check_los_match.py")
    if not os.path.exists(p):
        raise SystemExit(f"cannot find {p}")
    spec = importlib.util.spec_from_file_location("stage02", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def positions_from_los(path_a, path_b):
    """Transverse ray positions of both runs, in internal units."""
    sys.path.insert(0, os.path.join(ROOT, "legacy"))
    import swift_extract as sw
    st = load_stage02()
    ra, ma, _, _ = st.ray_positions(path_a, sw)
    rb, mb, _, _ = st.ray_positions(path_b, sw)
    if abs(ma.boxsize_int - mb.boxsize_int) / ma.boxsize_int > 1e-9:
        raise SystemExit(f"different box sizes: {ma.boxsize_int} vs "
                         f"{mb.boxsize_int}. The tiling would not be the "
                         f"same physical region in the two runs.")
    return ra, rb, float(ma.boxsize_int)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cdm", required=True)
    ap.add_argument("--fct", required=True)
    ap.add_argument("--los-a", default=None)
    ap.add_argument("--los-b", default=None)
    ap.add_argument("--pos-a", default=None)
    ap.add_argument("--pos-b", default=None)
    ap.add_argument("--box-int", type=float, default=None)
    ap.add_argument("--band", type=float, nargs=2, default=(0.0, 0.01))
    ap.add_argument("--nside", type=int, default=8)
    ap.add_argument("--n-boot", type=int, default=4000)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--no-rescale", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="figures/t9")
    args = ap.parse_args()

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    # ---- load -------------------------------------------------------------
    # load() twice, not load_pair(): load_pair enforces equal shapes because
    # everything it was written for compares LOS for LOS. This script exists
    # precisely because that comparison is invalid, and it does not need the
    # two runs to have the same number of sightlines. z and dv still have to
    # agree, so that is checked here instead.
    a = cachelib.load(args.cdm)
    b = cachelib.load(args.fct)
    if abs(a.z - b.z) > 0.02:
        raise SystemExit(f"redshift mismatch: {a.z:.4f} vs {b.z:.4f}")
    if abs(a.dv - b.dv) / a.dv > 1e-6:
        raise SystemExit(f"dv mismatch: {a.dv:.6f} vs {b.dv:.6f} km/s")
    dv, z = a.dv, a.z

    say("=" * 72)
    say("T9 - significance of the P1D suppression with no sightline pairing")
    say("=" * 72)
    say(f"z = {z:.4f}   dv = {dv:.5f} km/s")
    say(f"CDM  {a.tau.shape[0]} LOS x {a.tau.shape[1]} px   {args.cdm}")
    say(f"FCT  {b.tau.shape[0]} LOS x {b.tau.shape[1]} px   {args.fct}")
    if a.tau.shape[0] != b.tau.shape[0]:
        say("[i] different numbers of sightlines. Fine for everything below "
            "except the common-index bootstrap, which is skipped.")
    say()

    # ---- rescaling and the mean flux -------------------------------------
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(z)
    if args.no_rescale:
        A_cdm = A_fct = 1.0
    else:
        A_cdm, A_fct = solve_A(a.tau, target), solve_A(b.tau, target)
    F_cdm, F_fct = flux(a.tau, A_cdm), flux(b.tau, A_fct)
    mF_c, mF_f = float(F_cdm.mean()), float(F_fct.mean())

    say("--- mean flux " + "-" * 58)
    say(f"  target tau_eff  {target:.6f}"
        + ("   (A = 1, no rescaling)" if args.no_rescale else ""))
    say(f"  CDM   A = {A_cdm:.6f}   tau_eff = {tau_eff(a.tau, A_cdm):.6f}"
        f"   Fbar = {mF_c:.6f}")
    say(f"  FCT   A = {A_fct:.6f}   tau_eff = {tau_eff(b.tau, A_fct):.6f}"
        f"   Fbar = {mF_f:.6f}")
    say(f"  Fbar mismatch   {abs(mF_c - mF_f) / mF_c:.3e}")
    say()

    # ---- periodograms, once ----------------------------------------------
    k, Qc, mfc = boot.periodogram(F_cdm, dv)
    _, Qf, mff = boot.periodogram(F_fct, dv)
    dc = boot.check_against_p1d(k, Qc, mfc, dv, F_cdm)
    df = boot.check_against_p1d(k, Qf, mff, dv, F_fct)
    say(f"  estimator cross-check against common/p1d.py: "
        f"max deviation {max(dc, df):.2e}")
    say()

    kmax = k[nyquist_cut(k, dv, 0.5)].max()
    lo = args.band[0] if args.band[0] > 0 else None
    hi = args.band[1]
    kb, sel_c, Wc = boot.bin_operator(k, nbins=args.nbins, kmax=kmax)
    k_band, sel_b, Wb, n_modes = boot.band_operator(k, lo, hi)

    # one bootstrap covers curve and band: stack them side by side
    Qc_b = np.hstack([Qc[:, sel_c] @ Wc, Qc[:, sel_b] @ Wb])
    Qf_b = np.hstack([Qf[:, sel_c] @ Wc, Qf[:, sel_b] @ Wb])
    nb = kb.size
    del Qc, Qf

    Pc_full = boot.full_sample(Qc_b, mfc)
    Pf_full = boot.full_sample(Qf_b, mff)
    R_curve = Pf_full[:nb] / Pc_full[:nb]
    R_band = float(Pf_full[nb] / Pc_full[nb])

    say("--- the statistic " + "-" * 54)
    say(f"  band            {lo if lo else 0:.5g} < k < {hi:.5g} s/km"
        f"   ({n_modes} Fourier modes)")
    say(f"  k[s/km] -> k[h/Mpc]  x {units.DEFAULT.k_kms_to_hmpc(z):.4f}")
    say(f"  R_band = P1D(FCT)/P1D(CDM) = {R_band:.6f}")
    say(f"  suppression = {100 * (1 - R_band):.3f}%")
    say()

    rng = np.random.default_rng(args.seed)
    n_c, n_f = Qc_b.shape[0], Qf_b.shape[0]

    def band_ratio(Pf, Pc):
        return Pf[:, nb] / Pc[:, nb]

    results = {}

    # ---- 1. analytic ------------------------------------------------------
    # var of the mean of a ratio, propagated from the two band variances.
    # Independent lines assumed, which is the whole problem with it.
    qc, qf = Qc_b[:, nb] / mfc.mean() ** 2, Qf_b[:, nb] / mff.mean() ** 2
    s_c = qc.std(ddof=1) / np.sqrt(n_c) / Pc_full[nb]
    s_f = qf.std(ddof=1) / np.sqrt(n_f) / Pf_full[nb]
    results["analytic (std/sqrt N)"] = R_band * np.hypot(s_c, s_f)

    # ---- 2. bootstrap, common index (what the 20 sigma was) --------------
    if n_c == n_f:
        cnt = boot.draw_counts(n_c, args.n_boot, rng)
        Pc_p = boot.apply_counts(cnt, Qc_b, mfc)
        Pf_p = boot.apply_counts(cnt, Qf_b, mff)
        r_paired = band_ratio(Pf_p, Pc_p)
        results["bootstrap, common index"] = float(r_paired.std(ddof=1))
        del cnt
    else:
        r_paired = None
        say("[i] common-index bootstrap skipped: different n_los.")

    # ---- 3. bootstrap, unpaired ------------------------------------------
    Pc_u = boot.apply_counts(boot.draw_counts(n_c, args.n_boot, rng),
                             Qc_b, mfc)
    Pf_u = boot.apply_counts(boot.draw_counts(n_f, args.n_boot, rng),
                             Qf_b, mff)
    r_unpaired = band_ratio(Pf_u, Pc_u)
    results["bootstrap, unpaired"] = float(r_unpaired.std(ddof=1))

    # curve covariance, for the whole-band chi2 against R = 1
    R_draws = Pf_u[:, :nb] / Pc_u[:, :nb]
    C, Cinv, hart = boot.covariance(R_draws, hartlap=True)
    resid = R_curve - 1.0
    chi2_R = float(resid @ Cinv @ resid)
    p_R = float(chi2_dist.sf(chi2_R, nb))

    # ---- 4. block jackknife ----------------------------------------------
    sig_jk = None
    pos_a = pos_b = None
    box = args.box_int
    if args.los_a and args.los_b:
        pos_a, pos_b, box = positions_from_los(args.los_a, args.los_b)
    elif args.pos_a and args.pos_b:
        if box is None:
            raise SystemExit("--pos-a/--pos-b need --box-int")
        pos_a, pos_b = np.load(args.pos_a), np.load(args.pos_b)

    if pos_a is None:
        say("[!] BLOCK JACKKNIFE SKIPPED - no sightline positions given.")
        say("    Everything above resamples individual lines, which assumes")
        say("    they are independent. They are not. Pass --los-a/--los-b.")
        say()
    else:
        if len(pos_a) != n_c or len(pos_b) != n_f:
            raise SystemExit(
                f"positions do not match the caches: {len(pos_a)} vs {n_c} "
                f"(CDM), {len(pos_b)} vs {n_f} (FCT). The cache was probably "
                f"written from a different LOS file or with max_los set.")
        la = boot.tile_labels(pos_a, box, args.nside)
        lb = boot.tile_labels(pos_b, box, args.nside)
        tiles = np.arange(args.nside ** 2)
        occ_a = np.bincount(la, minlength=tiles.size)
        occ_b = np.bincount(lb, minlength=tiles.size)
        live = tiles[(occ_a > 0) & (occ_b > 0)]
        theta = []
        for t in live:
            rc = boot.subsample(Qc_b, mfc, la != t)
            rf = boot.subsample(Qf_b, mff, lb != t)
            theta.append(rf[nb] / rc[nb])
        sig_jk = float(np.sqrt(boot.jackknife_var(theta, len(live))))
        results[f"block jackknife ({len(live)} tiles)"] = sig_jk
        say("--- spatial blocks " + "-" * 53)
        say(f"  {args.nside} x {args.nside} tiles over the transverse plane, "
            f"{len(live)} with lines in both runs")
        say(f"  occupancy CDM  min {occ_a[live].min()}  "
            f"median {int(np.median(occ_a[live]))}  max {occ_a[live].max()}")
        say(f"  occupancy FCT  min {occ_b[live].min()}  "
            f"median {int(np.median(occ_b[live]))}  max {occ_b[live].max()}")
        say("  the same tile is deleted from both runs: region-level "
            "pairing, which survives stage 02")
        say()

    # ---- report -----------------------------------------------------------
    say("--- significance of R_band != 1 " + "-" * 41)
    say(f"{'method':<32s} {'sigma(R)':>10s} {'significance':>13s}")
    for name, s in results.items():
        say(f"{name:<32s} {s:>10.6f} {abs(1 - R_band) / s:>12.1f}s")
    say()

    if sig_jk is not None:
        s_line = results["bootstrap, unpaired"]
        n_eff = n_c * (s_line / sig_jk) ** 2
        say(f"  effective independent sightlines: {n_eff:.0f} of {n_c} "
            f"({100 * n_eff / n_c:.1f}%)")
        say("  that is the ratio of the two variances. Resampling lines")
        say("  overstates the information by roughly that factor, whatever")
        say("  the pairing does.")
        say()

    if r_paired is not None:
        rp, ru = results["bootstrap, common index"], \
            results["bootstrap, unpaired"]
        say("--- did the broken pairing matter? " + "-" * 38)
        say(f"  sigma(common index) / sigma(unpaired) = {rp / ru:.4f}")
        say(f"  correlation between the two runs' band powers under a")
        say(f"  common index draw: "
            f"{np.corrcoef(Pf_p[:, nb], Pc_p[:, nb])[0, 1]:+.4f}")
        if abs(rp / ru - 1.0) < 0.05:
            say("  VERDICT: the two agree to better than 5%. The bootstrap")
            say("  measured a covariance of zero because there was one, so")
            say("  the paired number was already the unpaired number. The")
            say("  quoted significance does not change materially, and the")
            say("  correction to make is to the *description* of the test,")
            say("  not to the number.")
        else:
            say("  VERDICT: they differ. The common-index draw is not")
            say("  equivalent to the unpaired one here, so the previously")
            say("  quoted number does change. Use the unpaired or, better,")
            say("  the block jackknife.")
        say()

    say("--- whole-band test against R(k) = 1 " + "-" * 36)
    say(f"  chi2 = {chi2_R:.1f} for {nb} bins   p = {p_R:.3e}")
    say(f"  bootstrap covariance, {args.n_boot} draws, "
        f"Hartlap factor {hart:.4f}")
    say("  this is the curve, not the band average, and it uses errors from")
    if sig_jk is not None:
        infl = (sig_jk / results["bootstrap, unpaired"]) ** 2
        say(f"  resampled lines, so divide it by about {infl:.1f} for the")
        say(f"  block-jackknife errors: chi2 ~ {chi2_R / infl:.1f} for {nb} bins.")
    else:
        say("  resampled lines, so it is an upper bound on the evidence.")
    say()

    # ---- figure -----------------------------------------------------------
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.2, 8.8),
                                 gridspec_kw={"height_ratios": [3, 2],
                                              "hspace": 0.28})
    e_boot = R_draws.std(axis=0, ddof=1)
    ax.fill_between(kb, R_curve - e_boot, R_curve + e_boot, color=C_FCT,
                    alpha=0.25, lw=0, label="bootstrap, unpaired")
    if sig_jk is not None:
        scale = sig_jk / results["bootstrap, unpaired"]
        ax.fill_between(kb, R_curve - scale * e_boot,
                        R_curve + scale * e_boot, color="0.4", alpha=0.22,
                        lw=0, label=f"block jackknife (x{scale:.1f})")
    ax.plot(kb, R_curve, color=C_FCT, lw=2.2)
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.axvspan(*units.desi_window(z), color="green", alpha=0.07, zorder=0)
    if lo or hi:
        ax.axvspan(max(lo or kb.min(), kb.min()), hi, color=C_CDM,
                   alpha=0.10, zorder=0, label="quoted band")
    ax.set(xscale="log", xlabel=r"$k$ [s/km]", ylabel="FCT / CDM",
           title=f"T9: P1D ratio at z = {z:.2f}, errors with no pairing")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2, which="both")

    bx.hist(r_unpaired, bins=60, color=C_FCT, alpha=0.55, density=True,
            label=f"unpaired  $\\sigma$={results['bootstrap, unpaired']:.5f}")
    if r_paired is not None:
        bx.hist(r_paired, bins=60, color=C_CDM, alpha=0.45, density=True,
                label=f"common index  $\\sigma$="
                      f"{results['bootstrap, common index']:.5f}")
    if sig_jk is not None:
        xs = np.linspace(r_unpaired.min(), r_unpaired.max(), 400)
        bx.plot(xs, np.exp(-0.5 * ((xs - R_band) / sig_jk) ** 2)
                / (sig_jk * np.sqrt(2 * np.pi)), color="0.2", lw=1.8,
                label=f"block jackknife  $\\sigma$={sig_jk:.5f}")
    bx.axvline(R_band, color="k", ls="--", lw=1)
    bx.set(xlabel=r"$R_{\rm band}$", ylabel="density",
           title="band statistic under resampling")
    bx.legend(frameon=False, fontsize=9)
    bx.grid(alpha=0.2)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    np.savez(args.out + "_draws.npz", kb=kb, R_curve=R_curve, C=C,
             r_unpaired=r_unpaired,
             r_paired=r_paired if r_paired is not None else np.array([]),
             R_band=R_band, z=z)
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    for s in (".png", ".txt", "_draws.npz"):
        print(f"written -> {args.out}{s}")


if __name__ == "__main__":
    main()
