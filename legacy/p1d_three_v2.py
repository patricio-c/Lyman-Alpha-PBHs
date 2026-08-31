#!/usr/bin/env python
"""
p1d_three_v2.py - P1D of Sherwood + CDM + FCT, talk-ready version.

Changes vs v1:
  * ratio panel binned in log k (the unbinned version is single-mode noise
    at low k and looks like structure that is not there)
  * x range truncated at --kmax (default 0.3 s/km): above that the pixel
    deconvolution and the particle shot noise dominate and the curves are
    not physical
  * errors combine both runs in quadrature, not just one
  * k_t of the FCT model marked

Usage:
    python p1d_three_v2.py --cdm cache_cdm.npz --fct cache_fct.npz \\
        --sherwood tauH1_lya_z3.0.dat --kt 10 --out p1d_three_v2.png
"""
import argparse

import numpy as np

from sherwood_los import (flux_power_1d, read_sherwood_spectra,
                          rescale_tau, tau_eff)

LAMBDA_LYA = 1215.6701
C_KMS = 2.99792458e5


def pk_per_los(tau, dv, A):
    tau = np.asarray(tau, dtype=np.float64)
    F = np.exp(-A * tau)
    npix = F.shape[1]
    d = F / F.mean() - 1.0
    k = 2.0 * np.pi * np.fft.rfftfreq(npix, d=dv)
    pk = (npix * dv / npix ** 2) * np.abs(np.fft.rfft(d, axis=1)) ** 2
    return k, pk / np.sinc(k * dv / (2.0 * np.pi)) ** 2


def binned(k, pk, edges):
    """Mean P and error on the mean, per log-k bin, over LOS and modes."""
    idx = np.digitize(k, edges) - 1
    kb, pb, eb = [], [], []
    nlos = pk.shape[0]
    for b in range(len(edges) - 1):
        m = idx == b
        if not m.any():
            continue
        per_los = pk[:, m].mean(axis=1)
        kb.append(k[m].mean())
        pb.append(per_los.mean())
        eb.append(per_los.std() / np.sqrt(nlos - 1))
    return np.array(kb), np.array(pb), np.array(eb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdm", required=True)
    ap.add_argument("--fct", required=True)
    ap.add_argument("--sherwood", required=True)
    ap.add_argument("--tau-eff", type=float, default=None)
    ap.add_argument("--kmax", type=float, default=0.3,
                    help="cut the plot here [s/km]; above this the "
                         "deconvolution and shot noise dominate")
    ap.add_argument("--per-decade", type=int, default=20)
    ap.add_argument("--kt", type=float, default=None,
                    help="FCT break scale in Mpc^-1 (e.g. 10)")
    ap.add_argument("--out", default="p1d_three_v2.png")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sh = read_sherwood_spectra(args.sherwood)
    tau_sh = np.asarray(sh.tau, dtype=np.float64)
    target = args.tau_eff if args.tau_eff is not None else tau_eff(tau_sh)

    A_sh = rescale_tau(tau_sh, target)
    k_sh, pk_sh = pk_per_los(tau_sh, sh.dv, A_sh)
    print(f"common tau_eff = {target:.5f}\nA(Sherwood) = {A_sh:.4f}")

    runs = {}
    for nm, path, col in [("CDM", args.cdm, "C0"), ("FCT", args.fct, "C3")]:
        d = np.load(path)
        if abs(float(d["z"]) - sh.z) > 0.02:
            raise SystemExit(f"{nm}: z mismatch with Sherwood")
        A = rescale_tau(np.asarray(d["tau"], np.float64), target)
        k, pk = pk_per_los(d["tau"], float(d["dv"]), A)
        runs[nm] = (k, pk, col, d["tau"].shape[0], A)
        print(f"A({nm}) = {A:.4f}   ({d['tau'].shape[0]} LOS)")

    kmin = min(k_sh[1], min(r[0][1] for r in runs.values()))
    nb = max(4, int(args.per_decade * np.log10(args.kmax / kmin)))
    edges = np.logspace(np.log10(kmin * 0.99), np.log10(args.kmax), nb + 1)

    kb_sh, pb_sh, eb_sh = binned(k_sh, pk_sh, edges)
    R_z = C_KMS * 0.8 / ((1.0 + sh.z) * LAMBDA_LYA)
    desi = (1.0e-3, 0.5 * np.pi / R_z)

    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(7.4, 7.6), sharex=True,
        gridspec_kw={"height_ratios": [1.7, 1], "hspace": 0.06})

    ax.plot(kb_sh, kb_sh * pb_sh / np.pi, color="0.35", lw=2.4,
            label=f"Sherwood ({sh.nlos} LOS)")
    for nm, (k, pk, col, nlos, A) in runs.items():
        kb, pb, eb = binned(k, pk, edges)
        ax.plot(kb, kb * pb / np.pi, color=col, lw=1.8,
                label=f"{nm} ({nlos} LOS)")
    ax.axvspan(*desi, color="green", alpha=0.07, zorder=0)
    ax.set(xscale="log", yscale="log",
           ylabel=r"$k\,P_{\rm 1D}(k)\,/\,\pi$",
           title=rf"z = {sh.z:.1f},  all rescaled to "
                 rf"$\tau_{{\rm eff}} = {target:.3f}$")
    ax.legend(frameon=False, loc="lower left", fontsize=9)
    ax.grid(alpha=0.2, which="both")

    bx.axhline(1.0, color="k", ls="--", lw=1)
    bx.axvspan(*desi, color="green", alpha=0.07, zorder=0)
    for nm, (k, pk, col, nlos, A) in runs.items():
        kb, pb, eb = binned(k, pk, edges)
        r = pb / pb_sh
        rel = np.sqrt((eb / pb) ** 2 + (eb_sh / pb_sh) ** 2)
        bx.fill_between(kb, r * (1 - rel), r * (1 + rel), color=col, alpha=0.20)
        bx.plot(kb, r, color=col, lw=1.8, label=nm)
    if args.kt:
        kt_kms = args.kt / (sh.hz / (1.0 + sh.z))
        for a in (ax, bx):
            a.axvline(kt_kms, color="0.4", ls=":", lw=1.2)
        bx.text(kt_kms, 1.42, r"  $k_{\rm t}$", color="0.35", fontsize=9)
    bx.set(xscale="log", xlim=(kmin * 0.9, args.kmax),
           ylim=(0.5, 1.5), xlabel=r"$k$  [s km$^{-1}$]",
           ylabel="ratio to Sherwood")
    bx.grid(alpha=0.2, which="both")
    bx.text(desi[0] * 1.15, 0.545, "DESI window", color="green", fontsize=9)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"written -> {args.out}")

    print(f"\n{'k':>9}{'CDM/Sh':>9}{'FCT/Sh':>9}{'FCT/CDM':>9}")
    kb_c, pb_c, _ = binned(*runs["CDM"][:2], edges)
    kb_f, pb_f, _ = binned(*runs["FCT"][:2], edges)
    for i in range(0, len(kb_sh), max(1, len(kb_sh) // 14)):
        print(f"{kb_sh[i]:9.4f}{pb_c[i]/pb_sh[i]:9.3f}"
              f"{pb_f[i]/pb_sh[i]:9.3f}{pb_f[i]/pb_c[i]:9.3f}")


if __name__ == "__main__":
    main()
