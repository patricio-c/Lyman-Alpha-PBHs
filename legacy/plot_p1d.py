#!/usr/bin/env python
"""
plot_p1d.py - Grafica tu P1D contra el de Sherwood, con panel de cociente.

Lleva ambos al MISMO tau_eff antes de comparar (sin eso la comparacion no
significa nada) y les aplica el MISMO tratamiento. Marca las escalas donde la
comparacion deja de ser confiable: el modo fundamental de cada caja, la
frecuencia de Nyquist de cada grilla y la ventana en k que recomienda DESI.

Uso:
    python plot_p1d.py prueba.hdf5 tauH1_lya_z3.0.dat \
        --treecool TREECOOL_HM12_G+Q --out p1d_vs_sherwood.png

    # con el post-proceso instrumental de Sherwood (LSF 7 km/s + rebin 3 km/s)
    python plot_p1d.py ... --instrumental
"""

import argparse

import numpy as np

from sherwood_los import (flux_power_1d, gamma_hi, read_sherwood_spectra,
                          rescale_tau, tau_eff)
from swift_extract import extract_all, open_los_file
from matplotlib.ticker import MultipleLocator

LAMBDA_LYA = 1215.6701
C_KMS = 2.99792458e5


def prep(tau, dv, target, instrumental):
    """Reescala a target tau_eff y, si se pide, aplica LSF + rebin."""
    A = rescale_tau(tau, target)
    F = np.exp(-A * np.asarray(tau, dtype=np.float64))
    if instrumental:
        from sherwood_postprocess import convolve_lsf, rebin_velocity
        F = convolve_lsf(F, dv, 7.0)
        F, dv = rebin_velocity(F, dv, 3.0)
    return F, dv, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("swift_los")
    ap.add_argument("sherwood")
    ap.add_argument("--treecool", required=True)
    ap.add_argument("--npix", type=int, default=2048)
    ap.add_argument("--max-los", type=int, default=None)
    ap.add_argument("--tau-eff", type=float, default=None,
                    help="tau_eff comun; por defecto el nominal de Sherwood")
    ap.add_argument("--instrumental", action="store_true",
                    help="aplicar LSF 7 km/s y rebin 3 km/s a AMBOS")
    ap.add_argument("--out", default="p1d_vs_sherwood.png")
    ap.add_argument("--label", default="Esta corrida")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- Sherwood ---------------------------------------------------------
    sh = read_sherwood_spectra(args.sherwood)
    tau_sh = np.asarray(sh.tau, dtype=np.float64)
    target = args.tau_eff if args.tau_eff is not None else tau_eff(tau_sh)

    # --- la corrida propia ------------------------------------------------
    meta = open_los_file(args.swift_los)
    if abs(meta.z - sh.z) > 0.02:
        raise SystemExit(f"Redshifts distintos: {meta.z} vs {sh.z}")
    G = gamma_hi(args.treecool, meta.z)
    print(f"z = {meta.z:.3f}, Gamma_HI = {G:.4e} s^-1")
    print(f"tau_eff comun = {target:.5f}")
    tau_me, dv_me, _ = extract_all(args.swift_los, args.npix, G,
                                   max_los=args.max_los)
    nlos_me = tau_me.shape[0]
    print(f"tau_eff nominal tuyo = {tau_eff(tau_me):.5f}  "
          f"(Sherwood {tau_eff(tau_sh):.5f})")

    F_me, dvo_me, A_me = prep(tau_me, dv_me, target, args.instrumental)
    F_sh, dvo_sh, A_sh = prep(tau_sh, sh.dv, target, args.instrumental)
    print(f"factor A: tuyo {A_me:.4f}, Sherwood {A_sh:.4f}")

    k_me, p_me, e_me = flux_power_1d(F_me, dvo_me)
    k_sh, p_sh, e_sh = flux_power_1d(F_sh, dvo_sh)

    # --- escalas de referencia -------------------------------------------
    k_fund_me = 2.0 * np.pi / meta.box_kms
    k_fund_sh = 2.0 * np.pi / sh.box_kms
    R_z = C_KMS * 0.8 / ((1.0 + meta.z) * LAMBDA_LYA)
    k_desi_lo, k_desi_hi = 1.0e-3, 0.5 * np.pi / R_z

    # --- figura -----------------------------------------------------------
    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(7.2, 8.0), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.06})

    m_me, m_sh = k_me > 0, k_sh > 0
    ax.errorbar(k_sh[m_sh], k_sh[m_sh] * p_sh[m_sh] / np.pi,
                yerr=k_sh[m_sh] * e_sh[m_sh] / np.pi,
                color="C1", lw=1.4, label=f"Sherwood ({sh.nlos} LOS)")
    ax.errorbar(k_me[m_me], k_me[m_me] * p_me[m_me] / np.pi,
                yerr=k_me[m_me] * e_me[m_me] / np.pi,
                color="C0", lw=1.4, label=f"{args.label} ({nlos_me} LOS)")

    ax.axvspan(k_desi_lo, k_desi_hi, color="green", alpha=0.07, zorder=0)
    ax.axvline(k_fund_me, color="C0", ls=":", lw=1)
    ax.axvline(k_fund_sh, color="C1", ls=":", lw=1)
    ax.set(xscale="log", yscale="log",
           ylabel=r"$k\,P_{\rm 1D}(k)/\pi$",
           title=(f"z = {meta.z:.1f},  "
                  + (r"$\tau_{\rm eff}$" + f" = {target:.3f}")
                  + ("  (con LSF + rebin)" if args.instrumental else "")))
    ax.legend(frameon=False, loc="lower left")
    ax.grid(alpha=0.2, which="both")
    ax.text(np.sqrt(k_desi_lo * k_desi_hi), ax.get_ylim()[1] * 0.4,
            "ventana DESI", color="green", ha="center", fontsize=9)

    lo = max(k_me[1], k_sh[1])
    hi = min(k_me[-1], k_sh[-1])
    kk = k_sh[(k_sh >= lo) & (k_sh <= hi)]
    r = np.interp(kk, k_me, p_me) / np.interp(kk, k_sh, p_sh)
    # error relativo propagado, dominado por la corrida chica
    rel = np.interp(kk, k_me, e_me) / np.maximum(np.interp(kk, k_me, p_me), 1e-300)

    bx.fill_between(kk, r * (1 - rel), r * (1 + rel),
                    color="C0", alpha=0.2)
    bx.plot(kk, r, color="C0", lw=1.4)

    bx.axhline(1.0, color="k", ls="--", lw=1)
    bx.axvspan(k_desi_lo, k_desi_hi, color="green", alpha=0.07, zorder=0)

    # Valores centrales
    bx.axvline(k_fund_me, color="C0", ls=":", lw=1)
    bx.axvline(k_fund_sh, color="C1", ls=":", lw=1)

    # ±20% alrededor de cada valor
    for k, color in [(k_fund_me, "C0"), (k_fund_sh, "C1")]:
        bx.axvline(0.8 * k, color=color, ls="--", lw=0.9, alpha=0.6)
        bx.axvline(1.2 * k, color=color, ls="--", lw=0.9, alpha=0.6)

    bx.set(
        xscale="log",
        ylim=(0, 2.0),
        xlabel=r"$k$ [s/km]",
        ylabel="cociente / Sherwood"
    )

    # 20 intervalos entre 0 y 2 → tick cada 0.1
    bx.yaxis.set_major_locator(MultipleLocator(0.1))

    bx.grid(alpha=0.2, which="both")


    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nescrito -> {args.out}")

    # --- resumen numerico -------------------------------------------------
    print(f"\nmodo fundamental: tuyo {k_fund_me:.5f}, Sherwood {k_fund_sh:.5f} s/km")
    print(f"ventana DESI a este z: {k_desi_lo:.5f} < k < {k_desi_hi:.5f} s/km")
    for lo_, hi_, nm in [(k_fund_me, 0.01, "escalas grandes"),
                         (0.01, 0.05, "intermedias"),
                         (0.05, 0.2, "chicas")]:
        s = (kk > lo_) & (kk < hi_)
        if s.any():
            print(f"  cociente medio en {nm:<16} "
                  f"({lo_:.4f} < k < {hi_:.3f}): {r[s].mean():.3f}")
    if nlos_me < 100:
        print(f"\n  OJO: con solo {nlos_me} LOS el error a k chico es enorme.")
        print("  La banda sombreada del panel inferior lo muestra. No leas")
        print("  nada de las escalas grandes hasta tener ~1000 LOS.")


if __name__ == "__main__":
    main()
