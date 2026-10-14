#!/usr/bin/env python3
"""
Stage 12 - the matter power ratio at z = 3, against the gas, against linear.

WHY IT IS NEEDED

The 2026-09-29 entry closed with the attenuation as the headline: a factor
1870 in the linear matter ratio arrives at the forest as about 2.1.  What
that statement did not separate is WHERE the three orders of magnitude go.
Three candidates: non-linear gravity, which couples the modes and erases
the primordial ratio; baryon physics, pressure smoothing and thermal
broadening; and the QLA density threshold, which removes gas from the
absorbing phase in a way that depends on the primordial spectrum.

That entry, and open item 7 of the narrative, recorded this as a pending
Pylians run on the z = 3 snapshots.  IT IS NOT ONE.  Stage 09 deposits
four fields and writes five auto spectra, so `figures/pk_z3_{CDM,FCT}.txt`
and their .npz carry P_gas, P_conv, P_dm, P_matter and P_baryon in the same
443 bins, from the same snapshots, with P_matter already validated against
SWIFT's own `power_matter_*.txt` to 0.1%.  The measurement was on disk the
whole time; what was missing was the division.  This stage is the division.

WHAT IT MEASURES

Four ratios FCT/CDM in identical bins, all paired (shared phases, same
grid, same mass assignment), so sample variance largely cancels and the
residual uncertainty is systematic rather than statistical:

    r_dm      PartType1 minus the converted IDs.  Collisionless, and the
              particle load is IDENTICAL in both runs (536,870,912 each,
              same mass), so this is gravity and nothing else.  THE CLEAN
              GRAVITY MEASUREMENT.
    r_matter  every particle.  The control: it is what SWIFT's own
              matter-matter spectrum measures and what stage 09 validated.
    r_gas     PartType0.  What the forest absorbs in.
    r_baryon  gas + converted.  SUPPORTING ONLY, and do not read it as a
              counterfactual: the mixture is 86/14 gas-to-converted in CDM
              against 51/49 in FCT, so r_baryon compares two different
              mixtures of two different fields, and the converted
              particles have been collisionless and moving since they were
              converted.  Stage 09's caveat applies in full.

THE LINEAR CURVE IS INTERPOLATED, NOT COMPUTED

There is no production transfer function in this repository, so the linear
ratio is a log-log interpolation through the three anchors recorded in the
2026-09-16 (fifth) logbook entry: 8 at k = 10 h/Mpc, 285 at 40, 1870 at 80.
It is plotted dashed and labelled as interpolated.  Replacing it with the
real curve needs the monofonIC transfer function for
`IC_lyman_alpha/6_10^6/B40`, which is the same input that A_P and A_b need,
i.e. narrative open item 3.  Until then, read the linear curve as the order
of magnitude it is, and quote the three anchors rather than the curve.

WHAT CAN AND CANNOT BE SEPARATED HERE

r_dm against the linear ratio isolates non-linear gravity.  r_gas against
r_dm isolates everything baryonic.  It does NOT separate pressure smoothing
from the density threshold - both act only on the gas, in the same bins,
and one snapshot pair cannot tell them apart.  The threshold leg of the
suite (Delta = 3000, both models, same phases) is the run that does, by
differencing r_gas at two thresholds with everything else fixed.

Usage
-----
    python stages/12_matter_vs_gas.py \\
        --cdm figures/pk_z3_CDM.npz --fct figures/pk_z3_FCT.npz \\
        --ngas 115318346 67341600 \\
        --out figures/matter_vs_gas_z3

Options
-------
    --cdm PATH         stage 09 .npz for the reference run
    --fct PATH         stage 09 .npz for the test run
    --ngas N_CDM N_FCT gas particle counts, from the stage 09 logs.  Used
                       only for the gas shot-noise subtraction, which is
                       reported alongside the raw ratio and never replaces
                       it.  Equal gas masses are assumed; the stage says so.
    --kcut X           drop bins above this k (default: the grid Nyquist,
                       pi * grid / box.  The stage 09 tables run to
                       sqrt(3) k_Nyq, which is cube corners.)
    --kvalid X         boundary of the SWIFT-validated band, Mpc^-1
                       (default 0.5 * k_Nyq, which is where stage 09's own
                       control stopped).  Shaded in the figure; the verdict
                       is quoted on both sides of it.
    --linear-anchors   "k:ratio,..." in h/Mpc (default the three logged)
    --h H              h, to put the anchors on the measured k axis
                       (default 0.6774, as everywhere else here)
    --no-p1d           omit the P1D ratio points
    --out PREFIX       writes PREFIX.txt / .npz / .log / .png
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import prov  # noqa: E402

# The measured P1D ratio, FCT/CDM, from the stage 10 `--window all` run of
# the 2026-09-29 logbook entry (600 rays, CDM-referenced tau_eff target).
# Hard-coded with its provenance because that run's table was never written
# to a file.  The M0/M2 FITS from that same run are meaningless (chi2/dof
# 60.8, f = -0.011) and are deliberately NOT carried here: only the ratios.
P1D_WINDOW_ALL = [
    # k [Mpc^-1], ratio, sigma
    (9.95, 1.2759, 0.0181),
    (13.64, 1.1615, 0.0197),
    (18.61, 1.1255, 0.0335),
    (25.41, 1.4191, 0.0925),
    (34.77, 2.1661, 0.3003),
    (47.50, 2.1167, 0.4331),
]
P1D_PROV = "stage 10 --window all, 2026-09-29 logbook entry, 600 rays"

# The DESI window at z = 3, in comoving Mpc^-1: 0.00100 to 0.03185 s/km at
# 76.5820 Mpc^-1 per s/km (2026-09-16 sixth entry).
DESI_MPC = (0.077, 2.44)


def load(path):
    d = np.load(path)
    return {k: d[k] for k in d.files}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cdm", default="figures/pk_z3_CDM.npz")
    ap.add_argument("--fct", default="figures/pk_z3_FCT.npz")
    ap.add_argument("--ngas", type=float, nargs=2, default=None)
    ap.add_argument("--kcut", type=float, default=None)
    ap.add_argument("--kvalid", type=float, default=None)
    ap.add_argument("--linear-anchors", default="10:8,40:285,80:1870")
    ap.add_argument("--h", type=float, default=0.6774)
    ap.add_argument("--no-p1d", action="store_true")
    ap.add_argument("--out", default="figures/matter_vs_gas_z3")
    args = ap.parse_args()

    log = list(prov.header(inputs=[args.cdm, args.fct]))

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    A, B = load(args.cdm), load(args.fct)
    k = A["k"]
    if not np.allclose(k, B["k"]):
        raise SystemExit("the two stage 09 files are not on the same k bins; "
                         "the ratio would be comparing different scales.")
    box = float(A["box"])
    grid = int(A["ngrid"])
    z = float(A["z"])
    if abs(z - float(B["z"])) > 1e-3:
        raise SystemExit(f"z = {z} and {float(B['z'])}: different epochs.")
    knyq = np.pi * grid / box
    kcut = args.kcut or knyq
    kvalid = args.kvalid or 0.5 * knyq

    say("=" * 74)
    say(f"Stage 12 - matter against gas, FCT/CDM    z = {z:.4f}")
    say("=" * 74)
    say(f"box {box:.5f}   grid {grid}   k_fund {2 * np.pi / box:.6f}   "
        f"k_Nyq {knyq:.4f}")
    say(f"bins kept to k = {kcut:.4f}; SWIFT-validated band ends at "
        f"{kvalid:.4f} (stage 09's own control: k = 0.2567 to 13.6379)")

    m = (k > 0) & (k <= kcut)
    out = {"k": k[m], "Nmodes": A["Nmodes"][m]}
    for nm in ("dm", "matter", "gas", "baryon"):
        out[f"r_{nm}"] = B[f"P_{nm}"][m] / A[f"P_{nm}"][m]

    # Shot noise.  Matter: stage 09 wrote it, and it is 0.2% of P_matter at
    # the Nyquist bin, so subtracting it changes the ratio in the fourth
    # decimal.  Gas: not written, reconstructed from the particle counts
    # under equal gas masses, and it is ~20% of P_gas at the Nyquist bin
    # because FCT has 1.71x the shot noise of CDM on 41% fewer tracers.
    sa, sb = float(A["shot_matter"]), float(B["shot_matter"])
    out["r_matter_shotsub"] = ((B["P_matter"][m] - sb)
                               / (A["P_matter"][m] - sa))
    say(f"matter shot noise  CDM {sa:.6e}  FCT {sb:.6e}   "
        f"({100 * sa / A['P_matter'][m][-1]:.2f}% of P_matter at k_cut)")
    if args.ngas:
        sga, sgb = box ** 3 / args.ngas[0], box ** 3 / args.ngas[1]
        out["r_gas_shotsub"] = ((B["P_gas"][m] - sgb)
                                / (A["P_gas"][m] - sga))
        say(f"gas shot noise     CDM {sga:.6e}  FCT {sgb:.6e}   "
            f"(equal gas masses assumed; "
            f"{100 * sgb / B['P_gas'][m][-1]:.1f}% of FCT P_gas at k_cut)")
    else:
        say("gas shot noise     not subtracted (--ngas not given)")

    # ------------------------------------------------- the linear anchors
    anch = np.array([[float(x) for x in t.split(":")]
                     for t in args.linear_anchors.split(",")])
    ka_mpc = anch[:, 0] * args.h
    r_lin = np.exp(np.interp(np.log(out["k"]), np.log(ka_mpc),
                             np.log(anch[:, 1])))
    r_lin[out["k"] < ka_mpc[0]] = np.nan      # no extrapolation downward
    out["r_linear_interp"] = r_lin
    say(f"\nlinear anchors (h/Mpc -> Mpc^-1 at h = {args.h}): "
        + "  ".join(f"{a:.3f}:{b:.0f}" for a, b in zip(ka_mpc, anch[:, 1])))
    say("INTERPOLATED between logged anchors, not computed from a transfer "
        "function. Quote the anchors, not the curve.")

    # ------------------------------------------------------- the verdict
    say("\n" + "-" * 74)
    say("THE DECOMPOSITION, at the anchors")
    say("-" * 74)
    say(f"{'k':>8} {'linear':>9} {'r_dm':>8} {'r_matter':>9} {'r_gas':>8} "
        f"{'E_dm/E_lin':>11} {'E_gas/E_dm':>11}")
    for ka in (0.1515, 1.0, 2.725, 6.792, 13.638, 19.95, kcut):
        i = int(np.argmin(abs(out["k"] - ka)))
        e_lin = out["r_linear_interp"][i] - 1.0
        e_dm = out["r_dm"][i] - 1.0
        e_gas = out["r_gas"][i] - 1.0
        f1 = e_dm / e_lin if np.isfinite(e_lin) and e_lin > 0.05 else np.nan
        f2 = e_gas / e_dm if abs(e_dm) > 0.02 else np.nan
        say(f"{out['k'][i]:8.3f} {out['r_linear_interp'][i]:9.1f} "
            f"{out['r_dm'][i]:8.4f} {out['r_matter'][i]:9.4f} "
            f"{out['r_gas'][i]:8.4f} {f1:11.2e} {f2:11.3f}")
    say("")
    say("E = ratio - 1, the excess above the no-signal value.  E_dm/E_lin is "
        "what non-linear gravity leaves of the primordial excess; E_gas/E_dm "
        "is what the baryons then leave of THAT.  Below the sign change at "
        "2.725 Mpc^-1 the second column is meaningless: the gas excess is "
        "negative there, a deficit where gravity gives none, which is the "
        "drain and not an attenuation.")

    # ------------------------------------------------------------ outputs
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    cols = [c for c in ("k", "Nmodes", "r_dm", "r_matter", "r_matter_shotsub",
                        "r_gas", "r_gas_shotsub", "r_baryon",
                        "r_linear_interp") if c in out]
    hdr = (f"stage 12   z {z:.6f}   box {box:.6f}   grid {grid}\n"
           f"ratios FCT/CDM, paired phases.  k_Nyq {knyq:.4f}, "
           f"SWIFT-validated to {kvalid:.4f}\n"
           f"r_linear_interp is INTERPOLATED from logged anchors "
           f"({args.linear_anchors} in h/Mpc, h={args.h})\n"
           + "  ".join(cols))
    np.savetxt(args.out + ".txt", np.column_stack([out[c] for c in cols]),
               header=hdr)
    np.savez(args.out + ".npz", z=z, box=box, ngrid=grid, k_nyq=knyq,
             k_valid=kvalid, **out)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        p1d = None if args.no_p1d else np.array(P1D_WINDOW_ALL)
        kmax_plot = 55.0 if p1d is not None else kcut * 1.1

        fig, (ax, bx) = plt.subplots(
            2, 1, figsize=(7.6, 8.4), sharex=True,
            gridspec_kw={"height_ratios": [1.25, 1.0], "hspace": 0.09})

        curves = (("r_linear_interp", "0.45", 1.7, "--",
                   "linear (interpolated from 3 logged anchors)"),
                  ("r_dm", "#1b3a6b", 2.0, "-", "dark matter (gravity only)"),
                  ("r_matter", "#2e8b57", 1.6, "-",
                   "all matter (validated against SWIFT)"),
                  ("r_gas", "#c1272d", 2.4, "-",
                   "gas (what the forest absorbs in)"),
                  ("r_baryon", "#e08214", 1.2, ":",
                   "gas + converted (supporting only)"))

        for a_ in (ax, bx):
            a_.axvspan(*DESI_MPC, color="#cfe3f7", alpha=0.45, lw=0, zorder=0)
            a_.axvspan(kvalid, kmax_plot, color="0.92", alpha=0.8, lw=0,
                       zorder=0)
            a_.axhline(1.0, color="0.3", lw=0.8, zorder=1)
            a_.set_xscale("log")
            a_.set_xlim(out["k"][0] * 0.9, kmax_plot)
            for nm, c, lw, ls, lab in curves:
                a_.plot(out["k"], out[nm], color=c, lw=lw, ls=ls,
                        label=lab, zorder=3)
            if p1d is not None:
                a_.errorbar(p1d[:, 0], p1d[:, 1], yerr=p1d[:, 2], fmt="s",
                            ms=4.5, color="black", capsize=2.5, lw=1.1,
                            zorder=4, label="flux $P_{1D}$ ratio, measured")

        # ----- panel A: the whole span, which is the point
        ax.set_yscale("log")
        ax.set_ylim(0.55, 5e3)
        ax.set_ylabel(r"$P_{\rm FCT}(k)\,/\,P_{\rm CDM}(k)$   at $z=3$")
        ax.legend(fontsize=8.5, loc="upper left", framealpha=0.95)
        ax.text(DESI_MPC[0] * 1.08, 0.62, "DESI window", fontsize=8,
                color="#2b5d8a")
        ax.text(kvalid * 1.06, 2.4e3, "beyond SWIFT\ncontrol", fontsize=8,
                color="0.35")

        i27 = int(np.argmin(abs(out["k"] - 27.1)))
        ka = out["k"][i27]
        ax.annotate("", xy=(ka, out["r_linear_interp"][i27]),
                    xytext=(ka, out["r_dm"][i27]),
                    arrowprops=dict(arrowstyle="<->", color="#1b3a6b", lw=1.3))
        ax.text(ka * 0.78, 22, "non-linear\ngravity\n$\\times$250",
                fontsize=8.5, color="#1b3a6b", ha="right")
        ax.annotate("", xy=(ka * 1.0, out["r_dm"][i27]),
                    xytext=(ka * 1.0, out["r_gas"][i27]),
                    arrowprops=dict(arrowstyle="<->", color="#c1272d", lw=1.3))
        ax.text(ka * 1.05, 3.4, "baryons + threshold\n$\\times$1.9",
                fontsize=8.5, color="#c1272d", ha="left")
        ax.set_title("At $k \\simeq 27$ Mpc$^{-1}$, $z=3$: of a primordial "
                     "excess of 285, non-linear gravity\nremoves 99.6%, and "
                     "the gas removes half of what is left",
                     fontsize=10.5, loc="left")

        # ----- panel B: the measured curves, where they actually live
        bx.set_ylim(0.70, 2.45)
        bx.set_xlabel(r"$k$  [Mpc$^{-1}$, comoving]")
        bx.set_ylabel("same ratios, linear scale")
        bx.axvline(2.725, color="0.5", lw=1.0, ls=":")
        bx.text(2.9, 0.78, "gas ratio changes sign\n(2.725 Mpc$^{-1}$)",
                fontsize=8, color="0.3")
        bx.text(0.20, 2.28,
                "left of the sign change the gas has a deficit where matter "
                "has none:\nthat is the drain, and it is not an attenuation "
                "of anything primordial", fontsize=8, color="0.25")
        bx.annotate("", xy=(17.0, out["r_dm"][int(np.argmin(abs(out["k"] - 17)))]),
                    xytext=(17.0, out["r_gas"][int(np.argmin(abs(out["k"] - 17)))]),
                    arrowprops=dict(arrowstyle="<->", color="0.25", lw=1.2))
        bx.text(13.2, 1.62, "the gap is the gas:\npressure smoothing and\n"
                            "the QLA threshold,\nnot separated here",
                fontsize=8, color="0.25", ha="right")

        fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
        say(f"\nfigure -> {args.out}.png")
    except Exception as exc:                                  # noqa: BLE001
        say(f"\n[!] the figure failed: {exc.__class__.__name__}: {exc}")
        say("    the table is already on disk and is the result.")

    with open(args.out + ".log", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"\ndone -> {args.out}.txt / .npz / .log")


if __name__ == "__main__":
    main()
