#!/usr/bin/env python3
"""
Test 10 - is the FCT signal in P1D shape information, or only amplitude?

The question decides the strategy of the paper, so it is worth stating
precisely.  At linear order an excess of 3D power concentrated above k_t
enters the 1D spectrum as

    Delta P1D(k_par) = (1 / 2 pi) * int_{k_t}^{inf} k P3D_excess(k) dk

which is a constant for k_par << k_t.  The linear signature is an additive
offset that is FLAT in k, not a ratio with structure in it.  The ratio
R(k) = P1D_FCT / P1D_CDM looks like it has shape, but most of that is
P1D_CDM falling with k underneath it: a constant offset divided by a
falling denominator produces a rising curve out of nothing.

So the diagnostic is: subtract, do not divide, and look for departures from
flatness.  Any real departure is non-linearity and thermal state, which is
the part no amplitude rescaling can imitate.

What "flat" is and is not
-------------------------
Flatness is the null hypothesis of this test, not a theorem.  The relation
above is linear theory for a field with scale-independent bias; the flux is
a non-linear function of the density and its bias is neither constant nor
local.  A flat Delta P1D says the data give no reason to claim shape
information.  It does not say linear theory has been verified.  Write it
that way in the paper.

Expected sign
-------------
There is a known amplitude deficit of about 0.83 in FCT, from the QLA
scheme converting ~50% of the FCT baryons plus pressure filtering.  That is
opposite in sign to the positive offset the linear term predicts, so
Delta P1D will most likely come out NEGATIVE.  That is not a bug and it is
not to be corrected away: it is the evidence that the removed-gas effect
dominates over the linear term.  Report it.

The condition the whole test rests on
-------------------------------------
The mean flux must be renormalised to the same value in both runs before
P1D is computed.  If it is not, Delta P1D is measuring the difference in
Fbar and the diagnostic is worthless: P1D of the contrast scales as
1/Fbar^2, so a few percent difference in Fbar is a few percent offset in
Delta P1D all by itself, flat in k, exactly mimicking the signal being
looked for.  This script rescales both runs to a common tau_eff by default,
prints the A and the Fbar it ended up with for each, and refuses to run if
they do not match.  --no-rescale exists for the cross-check against t0, and
prints a warning saying the result is not a shape diagnostic.

Errors
------
Bootstrap over sightlines, resampled independently in the two runs, because
stage 02 showed on 2026-08-31 that the sightlines are not paired.  The
covariance across k bins comes from the same draws, and the constant fit
uses it: neighbouring log bins of a P1D are correlated, and a chi2 built
from diagonal errors alone will overstate the evidence for structure.  Both
are reported so the difference is visible.

The caveat from t9 carries over: sightlines through one 40 Mpc/h box are
not independent, so these errors are optimistic.  Run t9 with the same
caches to get the block-jackknife inflation factor and pass it here with
--inflate, or read the chi2 as an upper bound on the evidence for shape.

Usage
-----
    python tests/t10_delta_p1d.py \\
        --pair cache/cache_cdm.npz cache/cache_fct.npz \\
        --out figures/t10

    # several redshifts, overlaid
    python tests/t10_delta_p1d.py \\
        --pair cache/cache_cdm_z3.0.npz cache/cache_fct_z3.0.npz \\
        --pair cache/cache_cdm_z2.6.npz cache/cache_fct_z2.6.npz \\
        --out figures/t10

Options
-------
    --pair CDM FCT     a cache pair, repeatable, one per redshift
    --nbins N          log bins in k (default 20)
    --n-boot N         bootstrap draws (default 4000)
    --kt-mpc X         k_t of the broken spectrum, Mpc^-1 (default 10.0).
                       Converted to s/km per redshift and drawn on the plot;
                       flatness is only predicted well below it.
    --fit-kmax-frac F  fit the constant below F * k_t (default 0.5)
    --inflate X        multiply every error bar by X, e.g. the ratio
                       sigma(block jackknife)/sigma(bootstrap) from t9
    --tau-eff X        common rescaling target (default Turner+24 at each z)
    --no-rescale       leave A = 1. Not a shape diagnostic, see above.
    --seed N           rng seed (default 1)
"""

from __future__ import annotations

import argparse
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


def analyse(path_cdm, path_fct, args, say):
    """Everything for one redshift.  Returns a dict for the plots."""
    a = cachelib.load(path_cdm)
    b = cachelib.load(path_fct)
    if abs(a.z - b.z) > 0.02:
        raise SystemExit(f"redshift mismatch: {a.z:.4f} vs {b.z:.4f}")
    if abs(a.dv - b.dv) / a.dv > 1e-6:
        raise SystemExit(f"dv mismatch: {a.dv:.6f} vs {b.dv:.6f} km/s")
    dv, z = a.dv, a.z

    say("=" * 72)
    say(f"z = {z:.4f}   dv = {dv:.5f} km/s")
    say(f"  CDM  {a.tau.shape[0]} LOS x {a.tau.shape[1]} px   {path_cdm}")
    say(f"  FCT  {b.tau.shape[0]} LOS x {b.tau.shape[1]} px   {path_fct}")

    # ---- the critical condition ------------------------------------------
    target = args.tau_eff if args.tau_eff is not None \
        else units.tau_eff_turner24(z)
    if args.no_rescale:
        A_cdm = A_fct = 1.0
    else:
        A_cdm, A_fct = solve_A(a.tau, target), solve_A(b.tau, target)
    F_cdm, F_fct = flux(a.tau, A_cdm), flux(b.tau, A_fct)
    mF_c, mF_f = float(F_cdm.mean()), float(F_fct.mean())
    rel = abs(mF_c - mF_f) / mF_c

    say("  --- mean flux, the condition this test rests on ---")
    say(f"    target tau_eff  {target:.6f}")
    say(f"    CDM   A = {A_cdm:.6f}   tau_eff = "
        f"{tau_eff(a.tau, A_cdm):.6f}   Fbar = {mF_c:.8f}")
    say(f"    FCT   A = {A_fct:.6f}   tau_eff = "
        f"{tau_eff(b.tau, A_fct):.6f}   Fbar = {mF_f:.8f}")
    say(f"    fractional Fbar mismatch  {rel:.3e}")
    if args.no_rescale:
        say("    [!] --no-rescale: Fbar differs between the runs, so "
            "Delta P1D")
        say("        contains the mean-flux difference and is NOT a shape")
        say("        diagnostic. Cross-check against t0 only.")
    elif rel > 1e-6:
        raise SystemExit(
            f"the two runs were rescaled to a common tau_eff but their mean "
            f"fluxes still differ by {rel:.2e}. Delta P1D would be measuring "
            f"that difference. Stop and find out why before going on.")
    else:
        say("    both runs sit at the same Fbar. Delta P1D is a shape "
            "diagnostic.")

    # ---- power ------------------------------------------------------------
    k, Qc, mfc = boot.periodogram(F_cdm, dv)
    _, Qf, mff = boot.periodogram(F_fct, dv)
    worst = max(boot.check_against_p1d(k, Qc, mfc, dv, F_cdm),
                boot.check_against_p1d(k, Qf, mff, dv, F_fct))
    say(f"  estimator cross-check against common/p1d.py: {worst:.2e}")

    kmax = k[nyquist_cut(k, dv, 0.5)].max()
    kb, sel, W = boot.bin_operator(k, nbins=args.nbins, kmax=kmax)
    Qc_b, Qf_b = Qc[:, sel] @ W, Qf[:, sel] @ W
    del Qc, Qf

    Pc = boot.full_sample(Qc_b, mfc)
    Pf = boot.full_sample(Qf_b, mff)
    dP = Pf - Pc
    R = Pf / Pc

    # ---- bootstrap --------------------------------------------------------
    rng = np.random.default_rng(args.seed)
    Pc_d = boot.apply_counts(boot.draw_counts(Qc_b.shape[0], args.n_boot, rng),
                             Qc_b, mfc)
    Pf_d = boot.apply_counts(boot.draw_counts(Qf_b.shape[0], args.n_boot, rng),
                             Qf_b, mff)
    dP_d = (Pf_d - Pc_d) * args.inflate
    R_d = Pf_d / Pc_d
    e_dP = dP_d.std(axis=0, ddof=1)
    e_R = R_d.std(axis=0, ddof=1) * args.inflate
    if args.inflate != 1.0:
        say(f"  error bars inflated by x{args.inflate:g} "
            f"(from a block jackknife, presumably t9's)")

    # ---- flatness ---------------------------------------------------------
    kt_kms = args.kt_mpc / units.DEFAULT.k_kms_to_mpc(z)
    fit = kb <= args.fit_kmax_frac * kt_kms
    say(f"  k_t = {args.kt_mpc:g} Mpc^-1 = {kt_kms:.5f} s/km at this z")
    say(f"  constant fitted over k < {args.fit_kmax_frac:g} k_t = "
        f"{args.fit_kmax_frac * kt_kms:.5f} s/km  ({int(fit.sum())} bins)")
    if fit.sum() < 3:
        raise SystemExit(
            f"only {int(fit.sum())} bins below {args.fit_kmax_frac:g} k_t. "
            f"Raise --nbins or --fit-kmax-frac; a constant fit to two bins "
            f"is not a measurement.")

    C, Cinv, hart = boot.covariance(dP_d[:, fit], hartlap=True)
    c_f, sc_f, chi2_f, dof_f = boot.fit_constant(dP[fit], Cinv=Cinv)
    c_d, sc_d, chi2_d, dof_d = boot.fit_constant(dP[fit], err=e_dP[fit])

    say("  --- flatness of Delta P1D ---")
    say(f"    {'':<22s} {'const [km/s]':>14s} {'sigma':>11s} "
        f"{'chi2/dof':>10s} {'p':>10s}")
    say(f"    {'full covariance':<22s} {c_f:>14.5e} {sc_f:>11.3e} "
        f"{chi2_f / dof_f:>10.3f} "
        f"{chi2_dist.sf(chi2_f, dof_f):>10.3e}")
    say(f"    {'diagonal errors':<22s} {c_d:>14.5e} {sc_d:>11.3e} "
        f"{chi2_d / dof_d:>10.3f} "
        f"{chi2_dist.sf(chi2_d, dof_d):>10.3e}")
    say(f"    dof = {dof_f}, bootstrap covariance from {args.n_boot} draws, "
        f"Hartlap {hart:.4f}")
    say(f"    offset significance (covariance): "
        f"{abs(c_f) / sc_f:.1f} sigma, sign "
        f"{'NEGATIVE' if c_f < 0 else 'POSITIVE'}")

    p_flat = float(chi2_dist.sf(chi2_f, dof_f))
    say("  --- reading ---")
    if p_flat > 0.05:
        say("    Delta P1D is flat within the errors. On this evidence the")
        say("    signal is pure amplitude, degenerate with Fbar bin by bin,")
        say("    and the z evolution is the only handle. That is the thesis")
        say("    of the paper, so this outcome supports the current plan.")
    else:
        say("    Delta P1D is NOT flat. There is shape information, and it")
        say("    is the part no amplitude rescaling can imitate. This")
        say("    changes the introduction: check first that it is not the")
        say("    top of the k range doing it (aliasing above half-Nyquist)")
        say("    by rerunning with a tighter --fit-kmax-frac.")
    if c_f < 0:
        say("    The offset is negative, opposite to the linear-term")
        say("    prediction. Expected: the QLA gas removal dominates over")
        say("    the linear excess. Report it, do not correct it.")
    else:
        say("    The offset is positive, the sign of the linear term. That")
        say("    is NOT what the known 0.83 amplitude deficit predicts.")
        say("    Check the rescaling before believing it.")
    say()

    # The bootstrap draws go out with the result, not just their spread:
    # t11 needs the covariance between k bins to fit anything to this curve,
    # and recomputing them there would mean a second definition of what
    # Delta P1D is. There is one, and it is here.
    return dict(z=z, kb=kb, dP=dP, e_dP=e_dP, R=R, e_R=e_R, fit=fit,
                c=c_f, sc=sc_f, chi2=chi2_f, dof=dof_f, kt=kt_kms,
                Fbar_cdm=mF_c, Fbar_fct=mF_f, A_cdm=A_cdm, A_fct=A_fct,
                dv=dv, nlos=(Qc_b.shape[0], Qf_b.shape[0]),
                Pc=Pc, Pf=Pf, dP_draws=dP_d, Pc_draws=Pc_d, R_draws=R_d)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", nargs=2, action="append", metavar=("CDM", "FCT"),
                    required=True)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--n-boot", type=int, default=4000)
    ap.add_argument("--kt-mpc", type=float, default=10.0)
    ap.add_argument("--fit-kmax-frac", type=float, default=0.5)
    ap.add_argument("--inflate", type=float, default=1.0)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--no-rescale", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="figures/t10")
    args = ap.parse_args()

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    say("T10 - additive Delta P1D: amplitude or shape?")
    say()
    out = [analyse(c, f, args, say) for c, f in args.pair]

    # ---- figure -----------------------------------------------------------
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8.4, 9.2),
                                 gridspec_kw={"hspace": 0.26})
    cmap = plt.get_cmap("viridis")
    zs = [o["z"] for o in out]
    norm = (lambda z: 0.5) if len(zs) < 2 else \
        (lambda z: (z - min(zs)) / (max(zs) - min(zs)))

    for o in out:
        col = cmap(0.15 + 0.7 * norm(o["z"]))
        lbl = f"z = {o['z']:.2f}"
        ax.errorbar(o["kb"], o["dP"], yerr=o["e_dP"], color=col, lw=1.6,
                    marker="o", ms=3.5, capsize=2, label=lbl)
        kf = o["kb"][o["fit"]]
        ax.plot([kf.min(), kf.max()], [o["c"], o["c"]], color=col, ls="--",
                lw=1.4)
        ax.fill_between([kf.min(), kf.max()], o["c"] - o["sc"],
                        o["c"] + o["sc"], color=col, alpha=0.18, lw=0)
        ax.axvline(o["kt"], color=col, ls=":", lw=1.0, alpha=0.7)
        bx.errorbar(o["kb"], o["R"], yerr=o["e_R"], color=col, lw=1.6,
                    marker="o", ms=3.5, capsize=2, label=lbl)

    ax.axhline(0.0, color="k", lw=1)
    ax.set(xscale="log", xlabel=r"$k$ [s/km]",
           ylabel=r"$\Delta P_{\rm 1D} = P^{\rm FCT}_{\rm 1D} - "
                  r"P^{\rm CDM}_{\rm 1D}$  [km/s]",
           title=r"T10: additive difference. Dashed = best-fit constant, "
                 r"dotted = $k_t$")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2, which="both")

    bx.axhline(1.0, color="k", ls=":", lw=1)
    bx.set(xscale="log", xlabel=r"$k$ [s/km]", ylabel="FCT / CDM",
           title="the same data as a ratio, for contrast")
    bx.legend(frameon=False, fontsize=9)
    bx.grid(alpha=0.2, which="both")

    say("=" * 72)
    say("SUMMARY")
    say(f"{'z':>7s} {'const [km/s]':>14s} {'sigma':>8s} {'chi2/dof':>10s} "
        f"{'p(flat)':>10s} {'Fbar':>10s}")
    for o in out:
        say(f"{o['z']:>7.3f} {o['c']:>14.5e} {abs(o['c']) / o['sc']:>7.1f}s "
            f"{o['chi2'] / o['dof']:>10.3f} "
            f"{chi2_dist.sf(o['chi2'], o['dof']):>10.3e} "
            f"{o['Fbar_cdm']:>10.6f}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    np.savez(args.out + "_data.npz",
             **{f"z{i}_{key}": np.asarray(o[key])
                for i, o in enumerate(out)
                for key in ("z", "kb", "dP", "e_dP", "R", "e_R", "c", "sc")})
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    for s in (".png", ".txt", "_data.npz"):
        print(f"written -> {args.out}{s}")


if __name__ == "__main__":
    main()
