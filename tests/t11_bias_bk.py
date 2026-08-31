#!/usr/bin/env python3
"""
Test 11 - the bias b(k,z), and everything the linear calculation omits.

This is GOAL.md Stage A past the model-free part. t10 asks whether
Delta P1D has shape. t11 asks what it is worth against the analytic
prediction, which means committing to A_P, A_b and k_max and dividing.

    b(k,z)  =  Delta P1D_flux(k,z) |_sim  /  Delta P1D_linear(k,z) |_theory

b is not a fudge factor. It is the quantity that absorbs, and therefore
measures, the non-linear density-to-flux mapping, thermal and pressure
smoothing, mode coupling in the Poisson term, and the QLA scheme's removal
of ~50% of the FCT baryons.  Whether it is constant in k decides the shape
of Paper I, and either answer is publishable.  A negative b is a legitimate
result: it means baryon removal and pressure filtering beat the linear
boost.  It is reported, not corrected.

The theory
----------
A discrete PBH population contributes white noise in 3D, P3D = A_P, and the
broken spectrum contributes A_b k^n_b above k_t.  Putting both through
P1D(k) = (1/2pi) int_k^inf k' P3D(k') dk', truncated at k_max:

    Delta P1D(k) = (A_P / 4pi) (k_max^2 - k^2)
                 + (A_b / (2pi (n_b+2))) (k_max^(n_b+2) - k^(n_b+2))

both multiplied by D^2(z)/D^2(z_ref), the only part of the redshift
evolution that is calculable rather than calibrated.

With n_b = 2 the broken term scales as k_max^4.  It is the most fragile
number in the whole programme, and --kmax-scan exists because of it: the
script reports d ln b / d ln k_max so the sensitivity is on the page rather
than implied.  Nothing from that term is quotable until resolution
convergence is shown (GOAL Stage E).

What this script computes, in the order GOAL asks for it
--------------------------------------------------------
    1. fit Delta P1D to (k_max^2 - k^2) with free amplitude, and to a pure
       constant, at equal dof.  Failure to distinguish them confirms
       k << k_max, which is the assumption the whole compression rests on.
    2. b(k,z), from A_P and A_b as GIVEN (from the initial conditions,
       known, not fitted) and k_max from the filtering scale.
    3. whether b is constant in k: the weighted mean, its error, the raw
       dispersion, and chi2/dof against a constant.
    4. the effective slope of P1D_CDM in the window, d ln P1D / d ln k,
       which says whether the ratio's denominator goes as k^-1 (n ~ -3) or
       logarithmically (n ~ -2).  Measured, not assumed.
    5. all of the above at every redshift given, with b(z) split into the
       D^2(z) part and the residual.

Units, because this is where it goes wrong
------------------------------------------
The k inside the P1D integral is comoving, in Mpc^-1.  The k we measure is
observed k_parallel, in s/km.  They differ by H(z)/(1+z), which is a
function of redshift, so the conversion cannot be done once and forgotten.
This script converts the measurement into comoving units before it touches
the theory, prints both, and refuses to guess the h convention of A_P and
A_b: pass --units to declare it.

    P1D[km/s] = P1D[Mpc] * H(z)/(1+z)          k[Mpc^-1] = k[s/km] * H/(1+z)

The window
----------
Fixed across all redshifts, per GOAL section 3.  desi_window() returns a
k_max that grows with z, and integrating over it bin by bin would put an
instrumental component into the measured redshift evolution.  The default
here is units.common_window() over exactly the redshifts given.

Usage
-----
    python tests/t11_bias_bk.py \\
        --pair cache/cache_cdm.npz cache/cache_fct.npz \\
        --ap 1.23e-3 --kmax-mpc 30 --units mpc \\
        --out figures/t11

Options
-------
    --pair CDM FCT     a cache pair, repeatable, one per redshift
    --ap A_P           white-noise amplitude of the Poisson term. REQUIRED.
                       From the initial conditions, not fitted.
    --ab A_b           amplitude of the broken term (default 0, Poisson
                       only). With --nb, its slope.
    --nb N             broken-spectrum slope (default 2.0)
    --kmax-mpc X       truncation scale, Mpc^-1. REQUIRED. Physically the
                       pressure filtering scale k_F(z).
    --kmax-scan LO HI N   also report b for N values of k_max in [LO, HI],
                       and d ln b / d ln k_max
    --units mpc|hmpc   the h convention of --ap, --ab and --kmax-mpc.
                       REQUIRED, no default: guessing this wrong is a
                       silent factor of h^3.
    --z-ref Z          redshift at which A_P and A_b are quoted (default 0)
    --window K1 K2     override the fixed window, s/km
    --nbins / --n-boot / --seed / --tau-eff / --no-rescale / --inflate
                       passed through to the t10 measurement
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
sys.path.insert(0, os.path.join(ROOT, "tests"))
from common import boot, units  # noqa: E402
import t10_delta_p1d as t10  # noqa: E402


def theory_dp1d(k_mpc, a_p, a_b, n_b, kmax, grow):
    """
    Delta P1D from linear theory, in Mpc (i.e. with k in Mpc^-1).

    Returns (total, poisson, broken) so the decomposition can be printed:
    with n_b = 2 the broken term carries k_max^4 and the Poisson term
    k_max^2, and knowing which one dominates is the difference between a
    robust number and one that moves with the resolution.
    """
    pois = (a_p / (4.0 * np.pi)) * (kmax ** 2 - k_mpc ** 2)
    if a_b:
        brok = (a_b / (2.0 * np.pi * (n_b + 2.0))) \
            * (kmax ** (n_b + 2.0) - k_mpc ** (n_b + 2.0))
    else:
        brok = np.zeros_like(k_mpc)
    return grow * (pois + brok), grow * pois, grow * brok


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pair", nargs=2, action="append", metavar=("CDM", "FCT"),
                    required=True)
    ap.add_argument("--ap", type=float, required=True)
    ap.add_argument("--ab", type=float, default=0.0)
    ap.add_argument("--nb", type=float, default=2.0)
    ap.add_argument("--kmax-mpc", type=float, required=True)
    ap.add_argument("--kmax-scan", type=float, nargs=3, default=None,
                    metavar=("LO", "HI", "N"))
    ap.add_argument("--units", choices=["mpc", "hmpc"], required=True)
    ap.add_argument("--z-ref", type=float, default=0.0)
    ap.add_argument("--window", type=float, nargs=2, default=None)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--n-boot", type=int, default=4000)
    ap.add_argument("--kt-mpc", type=float, default=10.0)
    ap.add_argument("--fit-kmax-frac", type=float, default=0.5)
    ap.add_argument("--inflate", type=float, default=1.0)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--no-rescale", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="figures/t11")
    args = ap.parse_args()

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    cos = units.DEFAULT
    h = cos.h

    # ---- units, declared and converted, once ------------------------------
    a_p, a_b, kmax = args.ap, args.ab, args.kmax_mpc
    say("T11 - the bias b(k,z) against linear theory")
    say("=" * 72)
    say("--- inputs and their units " + "-" * 45)
    say(f"  declared convention: {args.units}   h = {h:.4f}")
    if args.units == "hmpc":
        # A_P is a 3D power, [length^3]; A_b multiplies k^n_b, so it carries
        # [length^(3+n_b)]; k_max is an inverse length.
        a_p = a_p / h ** 3
        a_b = a_b / h ** (3.0 + args.nb)
        kmax = kmax * h
        say(f"  converted from h units:")
        say(f"    A_P  {args.ap:.6e} (Mpc/h)^3        -> {a_p:.6e} Mpc^3")
        if args.ab:
            say(f"    A_b  {args.ab:.6e} (Mpc/h)^(3+n_b) -> "
                f"{a_b:.6e} Mpc^(3+n_b)")
        say(f"    k_max {args.kmax_mpc:.5g} h/Mpc         -> "
            f"{kmax:.5g} Mpc^-1")
    else:
        say(f"    A_P   = {a_p:.6e} Mpc^3")
        if args.ab:
            say(f"    A_b   = {a_b:.6e} Mpc^(3+n_b),  n_b = {args.nb:g}")
        say(f"    k_max = {kmax:.5g} Mpc^-1")
    say(f"  A_P and A_b are taken as GIVEN, from the initial conditions.")
    say(f"  If they were fitted, b would be meaningless.")
    say(f"  reference redshift for the amplitudes: z = {args.z_ref:g}")
    say()

    # ---- the fixed window -------------------------------------------------
    out = []
    zs_probe = []
    for c, f in args.pair:
        import numpy as _np
        with _np.load(c) as d:
            zs_probe.append(float(d["z"]))
    if args.window:
        k1, k2 = args.window
        src = "given on the command line"
    else:
        k1, k2 = units.common_window(zs_probe)
        src = f"common DESI intersection over z = " \
              f"{', '.join(f'{z:.2f}' for z in sorted(zs_probe))}"
    say("--- the window, fixed for every redshift " + "-" * 31)
    say(f"  k1 = {k1:.6g} s/km   k2 = {k2:.6g} s/km")
    say(f"  source: {src}")
    say("  fixed on purpose: desi_window()'s upper edge grows with z, and")
    say("  letting it float would put an instrumental component into the")
    say("  measured redshift evolution (GOAL section 3).")
    say()

    # ---- per redshift -----------------------------------------------------
    for c, f in args.pair:
        o = t10.analyse(c, f, args, say)
        z = o["z"]
        conv = cos.k_kms_to_mpc(z)          # H(z)/(1+z), km/s/Mpc
        grow = (units.growth_factor(z) / units.growth_factor(args.z_ref)) ** 2

        w = (o["kb"] >= k1) & (o["kb"] <= k2)
        if w.sum() < 4:
            raise SystemExit(
                f"z = {z:.2f}: only {int(w.sum())} bins inside the fixed "
                f"window. Raise --nbins, or the window is wrong.")

        k_mpc = o["kb"][w] * conv
        dP_mpc = o["dP"][w] / conv
        dPd_mpc = o["dP_draws"][:, w] / conv

        say(f"  --- comoving frame at z = {z:.3f} ---")
        say(f"    H(z)/(1+z) = {conv:.4f} km/s/Mpc")
        say(f"    window     = {k_mpc.min():.5g} to {k_mpc.max():.5g} Mpc^-1"
            f"   ({int(w.sum())} bins)")
        say(f"    D^2(z)/D^2({args.z_ref:g}) = {grow:.6e}")
        say(f"    (k_max/k)^-2 at the top of the window: "
            f"{(k_mpc.max() / kmax) ** 2:.3e}")

        C, Cinv, hart = boot.covariance(dPd_mpc, hartlap=True)

        # --- 1. constant vs (k_max^2 - k^2), equal dof --------------------
        tmpl = kmax ** 2 - k_mpc ** 2
        c_c, sc_c, x2_c, dof_c = boot.fit_constant(dP_mpc, Cinv=Cinv)
        c_t, sc_t, x2_t, dof_t = boot.fit_template(dP_mpc, tmpl, Cinv=Cinv)
        say("    --- 1. shape of Delta P1D, at equal dof ---")
        say(f"      constant           {c_c:+.6e} +- {sc_c:.3e} Mpc"
            f"   chi2/dof = {x2_c / dof_c:.4f}")
        say(f"      A (k_max^2 - k^2)  {c_t:+.6e} +- {sc_t:.3e}"
            f"        chi2/dof = {x2_t / dof_t:.4f}")
        say(f"      delta chi2 (const - template) = {x2_c - x2_t:+.6f} "
            f"for {dof_c} dof")
        if abs(x2_c - x2_t) < 1.0:
            say("      the two are indistinguishable, which is the")
            say("      measurement that k << k_max. The compression is safe.")
        else:
            say("      the two are distinguishable. Either k is not far")
            say("      below k_max, or something else has k-dependence.")

        # --- 4. effective slope of P1D_CDM --------------------------------
        Pc, Pcd = o["Pc"][w], o["Pc_draws"][:, w]
        lp, lpd = np.log(Pc), np.log(Pcd)
        Cl, Clinv, _ = boot.covariance(lpd, hartlap=True)
        M = np.column_stack([np.ones(k_mpc.size), np.log(k_mpc)])
        th, cov_th, x2_s, dof_s = boot.fit_linear(lp, M, Cinv=Clinv)
        slope, e_slope = th[1], float(np.sqrt(cov_th[1, 1]))
        say("    --- 4. effective slope of P1D_CDM in the window ---")
        say(f"      dlnP1D/dlnk = {slope:+.4f} +- {e_slope:.4f}"
            f"   chi2/dof = {x2_s / dof_s:.3f}")
        say(f"      implied 3D slope n = dlnP1D/dlnk - 2 = "
            f"{slope - 2:+.4f}")
        if slope < -0.7:
            say("      the ratio's denominator falls roughly as k^-1 or")
            say("      faster: the apparent shape of R(k) is mostly this.")
        elif slope > -0.3:
            say("      the denominator is nearly flat, so the integral is")
            say("      closer to logarithmic. R(k) shape is NOT explained")
            say("      by the denominator alone.")
        else:
            say("      intermediate. Quote the number, not a regime name.")

        # --- 2 and 3. b(k) -------------------------------------------------
        th_tot, th_p, th_b = theory_dp1d(k_mpc, a_p, a_b, args.nb, kmax, grow)
        if np.any(th_tot == 0):
            raise SystemExit("theoretical Delta P1D is zero somewhere; "
                             "check A_P, A_b and k_max")
        bk = dP_mpc / th_tot
        bd = dPd_mpc / th_tot[None, :]
        e_bk = bd.std(axis=0, ddof=1)
        Cb, Cbinv, _ = boot.covariance(bd, hartlap=True)
        b_const, e_b, x2_b, dof_b = boot.fit_constant(bk, Cinv=Cbinv)
        p_b = float(chi2_dist.sf(x2_b, dof_b))

        say("    --- 2/3. the bias b(k) ---")
        say(f"      theory at the window centre: total {np.median(th_tot):+.5e}"
            f" Mpc")
        say(f"        Poisson term {100 * np.median(th_p) / np.median(th_tot):6.2f}%"
            f"   broken term {100 * np.median(th_b) / np.median(th_tot):6.2f}%")
        if a_b and abs(np.median(th_b)) > abs(np.median(th_tot)) * 0.2:
            say("      [!] the broken term carries more than 20% of the")
            say("          prediction and it scales as k_max^(n_b+2). Until")
            say("          resolution convergence is shown (GOAL stage E),")
            say("          treat b as provisional.")
        say(f"      b weighted mean = {b_const:+.5f} +- {e_b:.5f}")
        say(f"      raw dispersion across bins = {bk.std(ddof=1):.5f}"
            f"   min {bk.min():+.5f}  max {bk.max():+.5f}")
        say(f"      constant in k? chi2/dof = {x2_b / dof_b:.4f} "
            f"({dof_b} dof), p = {p_b:.3e}")
        if p_b > 0.05:
            say("      b is constant in k within the errors: a predictive")
            say("      one-parameter-per-redshift model. Strong result.")
        else:
            say("      b varies with k. That variation IS the shape")
            say("      information no amplitude rescaling can mimic.")
            say("      Also strong, but it changes the paper.")
        if b_const < 0:
            say("      b is NEGATIVE: baryon removal and pressure filtering")
            say("      beat the linear boost. Legitimate and publishable.")
            say("      Report it; do not correct it away.")

        # --- k_max sensitivity --------------------------------------------
        scan = None
        if args.kmax_scan:
            slo, shi, sn = args.kmax_scan
            grid = np.geomspace(slo, shi, int(sn))
            if args.units == "hmpc":
                grid = grid * h
            bs = []
            for km in grid:
                tt, _, _ = theory_dp1d(k_mpc, a_p, a_b, args.nb, km, grow)
                bs.append(float(np.mean(dP_mpc / tt)))
            bs = np.array(bs)
            dlnb = np.gradient(np.log(np.abs(bs)), np.log(grid))
            scan = (grid, bs, dlnb)
            say("    --- k_max sensitivity ---")
            say(f"      dln|b|/dln k_max over the scan: "
                f"{dlnb.min():+.3f} to {dlnb.max():+.3f}")
            say(f"      (pure Poisson would give -2, pure broken with "
                f"n_b={args.nb:g} would give {-(args.nb + 2):+.0f})")
            say(f"      a 10% error on k_max moves b by about "
                f"{100 * abs(np.median(dlnb)) * 0.1:.1f}%")
        say()

        o.update(k_mpc=k_mpc, dP_mpc=dP_mpc, e_dP_mpc=dPd_mpc.std(axis=0,
                                                                  ddof=1),
                 w=w, conv=conv, grow=grow,
                 bk=bk, e_bk=e_bk, b_const=b_const, e_b=e_b,
                 x2_b=x2_b, dof_b=dof_b, p_b=p_b, slope=slope,
                 e_slope=e_slope, th_tot=th_tot, scan=scan,
                 x2_c=x2_c, dof_c=dof_c, x2_t=x2_t, c_c=c_c, c_t=c_t)
        out.append(o)

    # ---- 5. b(z), and what D^2 already explains ---------------------------
    out.sort(key=lambda d: d["z"])
    say("=" * 72)
    say("--- 5. redshift dependence " + "-" * 45)
    say(f"{'z':>7s} {'b':>11s} {'+-':>10s} {'chi2/dof':>9s} "
        f"{'dlnP/dlnk':>10s} {'D^2 ratio':>11s} {'dP_sim ratio':>13s} "
        f"{'residual':>10s}")
    ref = out[0]
    for o in out:
        d2 = o["grow"] / ref["grow"]
        sim = float(np.mean(o["dP_mpc"]) / np.mean(ref["dP_mpc"]))
        say(f"{o['z']:>7.3f} {o['b_const']:>+11.5f} {o['e_b']:>10.5f} "
            f"{o['x2_b'] / o['dof_b']:>9.3f} {o['slope']:>+10.4f} "
            f"{d2:>11.4f} {sim:>13.4f} {sim / d2:>10.4f}")
    say(f"  ratios are relative to z = {ref['z']:.3f}")
    say("  'D^2 ratio' is the part of the evolution that is calculable.")
    say("  'residual' is what is left after dividing it out: filtering")
    say("  scale, thermal state, and the QLA removal. That residual is the")
    say("  thing Paper I has to model, and it is not the growth factor.")
    if len(out) < 2:
        say()
        say("  [!] ONE redshift. b(z) is the whole idea and it cannot be")
        say("      measured from one point. Extract more caches before")
        say("      reading anything into the table above.")
    say()

    # ---- figures ----------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(12.4, 9.4))
    (ax, bx), (cx, dx) = axs
    cmap = plt.get_cmap("viridis")
    zs = [o["z"] for o in out]
    nrm = (lambda z: 0.5) if len(zs) < 2 else \
        (lambda z: (z - min(zs)) / (max(zs) - min(zs)))

    for o in out:
        col = cmap(0.15 + 0.7 * nrm(o["z"]))
        lbl = f"z = {o['z']:.2f}"
        ax.errorbar(o["k_mpc"], o["bk"], yerr=o["e_bk"], color=col, lw=1.7,
                    marker="o", ms=4, capsize=2, label=lbl)
        ax.axhline(o["b_const"], color=col, ls="--", lw=1.0, alpha=0.6)
        bx.errorbar(o["k_mpc"], o["dP_mpc"], yerr=o["e_dP_mpc"],
                    color=col, lw=1.6, marker="o", ms=3.5, capsize=2,
                    label=lbl + " sim")
        bx.plot(o["k_mpc"], o["th_tot"], color=col, ls=":", lw=2.0,
                label=lbl + " linear")
        if o["scan"] is not None:
            g, bs, _ = o["scan"]
            cx.plot(g, bs, color=col, lw=1.8, marker=".", label=lbl)

    ax.axhline(0.0, color="k", lw=1)
    ax.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]", ylabel=r"$b(k)$",
           title="T11: the central figure. Dashed = weighted mean per $z$")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2, which="both")

    bx.axhline(0.0, color="k", lw=1)
    bx.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]",
           ylabel=r"$\Delta P_{\rm 1D}$ [Mpc]",
           title="measurement against linear theory")
    bx.legend(frameon=False, fontsize=8, ncol=2)
    bx.grid(alpha=0.2, which="both")

    if any(o["scan"] is not None for o in out):
        cx.set(xscale="log", xlabel=r"$k_{\rm max}$ [Mpc$^{-1}$]",
               ylabel=r"$\bar b$",
               title=r"sensitivity to $k_{\rm max}$ (the fragile input)")
        cx.legend(frameon=False, fontsize=9)
        cx.grid(alpha=0.2, which="both")
    else:
        cx.text(0.5, 0.5, "no --kmax-scan given", ha="center", va="center",
                transform=cx.transAxes, color="0.5")
        cx.set_axis_off()

    zz = np.array([o["z"] for o in out])
    bb = np.array([o["b_const"] for o in out])
    ee = np.array([o["e_b"] for o in out])
    dx.errorbar(zz, bb, yerr=ee, color="#d62728", lw=1.8, marker="o", ms=5,
                capsize=3, label=r"$b(z)$, $D^2(z)$ already divided out")
    dx.axhline(0.0, color="k", lw=1)
    dx.set(xlabel="$z$", ylabel=r"$\bar b$",
           title=r"what the growth factor does NOT explain")
    dx.legend(frameon=False, fontsize=9)
    dx.grid(alpha=0.2)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    np.savez(args.out + "_data.npz",
             **{f"z{i}_{key}": np.asarray(o[key])
                for i, o in enumerate(out)
                for key in ("z", "k_mpc", "dP_mpc", "bk", "e_bk", "b_const",
                            "e_b", "slope", "th_tot", "grow")})
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    for s in (".png", ".txt", "_data.npz"):
        print(f"written -> {args.out}{s}")


if __name__ == "__main__":
    main()
