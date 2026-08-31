#!/usr/bin/env python
"""
ratio_fct_cdm.py - Cociente P1D(FCT)/P1D(CDM) desde dos archivos de LOS.

Por que el cociente y no el P1D absoluto:

  - La desviacion por resolucion de masa (Bolton et al. 2017, figura A4, curva
    40-512-cc) es identica en las dos corridas y se cancela.
  - Lo mismo el tamano de caja, el estimador, la ventana de pixel y, si los
    rayos comparten semilla, buena parte de la varianza de muestreo.
  - Lo que NO se cancela es el reescalado de tau, porque F = exp(-A*tau) es no
    lineal. Por eso ambas se llevan al MISMO tau_eff y se reportan los dos A:
    la diferencia entre ellos ya es informacion fisica (cuanta absorcion extra
    produce el modelo).

Uso:
    python ratio_fct_cdm.py cdm.hdf5 fct.hdf5 --treecool TREECOOL_HM12_G+Q \\
        --tau-eff 0.3719 --out ratio_z3.txt --plot ratio_z3.png

    # con los datos de DESI encima (texto con columnas k, P1D, error)
    python ratio_fct_cdm.py cdm.hdf5 fct.hdf5 --treecool ... \\
        --tau-eff 0.3719 --desi desi_dr1_z3.0.txt --plot ratio_z3.png
"""

import argparse

import numpy as np

from sherwood_los import flux_power_1d, gamma_hi, rescale_tau, tau_eff
from swift_extract import extract_all, open_los_file

LAMBDA_LYA = 1215.6701
C_KMS = 2.99792458e5


def p1d_of(path, npix, treecool, target, max_los):
    meta = open_los_file(path)
    G = gamma_hi(treecool, meta.z)
    tau, dv, _ = extract_all(path, npix, G, max_los=max_los)
    tau = np.asarray(tau, dtype=np.float64)
    nom = tau_eff(tau)
    A = rescale_tau(tau, target)
    k, p, e = flux_power_1d(np.exp(-A * tau), dv, deconvolve_pixel=True)
    return meta, k, p, e, A, nom, tau.shape[0], dv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdm")
    ap.add_argument("fct")
    ap.add_argument("--treecool", required=True)
    ap.add_argument("--npix", type=int, default=2048)
    ap.add_argument("--max-los", type=int, default=None)
    ap.add_argument("--tau-eff", type=float, required=True,
                    help="tau_eff comun. Para comparar con DESI usa el "
                         "observado: 0.3719 a z=3 (Turner et al. 2024).")
    ap.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    ap.add_argument("--desi", default=None,
                    help="texto con k, P1D, error de DESI al mismo z, solo "
                         "para marcar la ventana util y la precision")
    ap.add_argument("--out", default=None)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    la, lb = args.labels
    res = {}
    for path, lab in [(args.cdm, la), (args.fct, lb)]:
        print(f"\n--- {lab}: {path}")
        res[lab] = p1d_of(path, args.npix, args.treecool, args.tau_eff,
                          args.max_los)
        meta, k, p, e, A, nom, nlos, dv = res[lab]
        print(f"    {nlos} LOS, dv = {dv:.4f} km/s")
        print(f"    tau_eff nominal = {nom:.5f}  ->  A = {A:.5f}")

    ma = res[la][0]
    z = ma.z
    A_a, A_b = res[la][4], res[lb][4]
    print(f"\nAmbas llevadas a tau_eff = {args.tau_eff:.5f}")
    print(f"  A({la}) = {A_a:.5f}   A({lb}) = {A_b:.5f}   "
          f"cociente = {A_b/A_a:.5f}")
    print(f"  Interpretacion: {lb} necesita un factor {A_b/A_a:.3f} veces el de")
    print(f"  {la} para dar el mismo flujo medio. Eso ya es una medida de")
    print(f"  cuanta absorcion extra (o menos) produce el modelo.")

    # cociente en la grilla comun
    ka, pa, ea = res[la][1], res[la][2], res[la][3]
    kb, pb, eb = res[lb][1], res[lb][2], res[lb][3]
    lo = max(ka[1], kb[1])
    hi = min(ka[-1], kb[-1]) * 0.5      # media Nyquist, convencion de DESI
    m = (ka >= lo) & (ka <= hi)
    kk = ka[m]
    pb_i = np.interp(kk, kb, pb)
    eb_i = np.interp(kk, kb, eb)
    r = pb_i / pa[m]
    # error relativo combinado
    rel = np.sqrt((eb_i / np.maximum(pb_i, 1e-300)) ** 2
                  + (ea[m] / np.maximum(pa[m], 1e-300)) ** 2)

    k_fund = 2.0 * np.pi / ma.box_kms
    R_z = C_KMS * 0.8 / ((1.0 + z) * LAMBDA_LYA)
    k_desi = (1.0e-3, 0.5 * np.pi / R_z)

    print(f"\n  {'k [s/km]':>10}{'cociente':>11}{'error':>10}")
    for i in np.unique(np.geomspace(1, kk.size - 1, 14).astype(int)):
        print(f"  {kk[i]:10.5f}{r[i]:11.4f}{r[i]*rel[i]:10.4f}")

    for a_, b_, nm in [(k_fund, 0.01, "escalas grandes"),
                       (0.01, 0.02, "intermedias"),
                       (0.02, k_desi[1], "hasta k_max DESI")]:
        s = (kk > a_) & (kk < b_)
        if s.any():
            w = 1.0 / np.maximum(rel[s], 1e-6) ** 2
            print(f"  cociente medio en {nm:<18} ({a_:.4f} < k < {b_:.4f}): "
                  f"{np.sum(r[s]*w)/np.sum(w):.4f}")

    if args.out:
        np.savetxt(args.out, np.column_stack([kk, pa[m], pb_i, r, r*rel]),
                   header=(f"z={z:.4f} tau_eff={args.tau_eff:.6f} "
                           f"A_{la}={A_a:.6f} A_{lb}={A_b:.6f}\n"
                           f"k[s/km]  P1D_{la}  P1D_{lb}  cociente  error"))
        print(f"\nescrito -> {args.out}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.axhline(1.0, color="k", ls="--", lw=1)
        ax.axvspan(*k_desi, color="green", alpha=0.07, zorder=0)
        ax.axvline(k_fund, color="grey", ls=":", lw=1)
        ax.fill_between(kk, r*(1-rel), r*(1+rel), color="C3", alpha=0.25)
        ax.plot(kk, r, color="C3", lw=1.8)
        if args.desi:
            d = np.loadtxt(args.desi)
            kd, pd = d[:, 0], d[:, 1]
            ed = d[:, 2] if d.shape[1] > 2 else np.zeros_like(pd)
            s = (kd >= lo) & (kd <= hi)
            ax.errorbar(kd[s], np.ones(s.sum()), yerr=ed[s]/pd[s],
                        fmt="o", ms=3, color="k", alpha=0.5, capsize=2,
                        label="precision DESI DR1")
            ax.legend(frameon=False, loc="upper left")
        ax.set(xscale="log", xlabel=r"$k$ [s/km]",
               ylabel=f"$P_{{\\rm 1D}}$({lb}) / $P_{{\\rm 1D}}$({la})",
               title=(f"z = {z:.1f}, mismo " + r"$\tau_{\rm eff}$"
                      + f" = {args.tau_eff:.4f}"))
        ax.grid(alpha=0.2, which="both")
        ax.text(np.sqrt(k_desi[0]*k_desi[1]), ax.get_ylim()[1]*0.97,
                "ventana DESI", color="green", ha="center", va="top",
                fontsize=9)
        fig.savefig(args.plot, dpi=150, bbox_inches="tight")
        print(f"escrito -> {args.plot}")


if __name__ == "__main__":
    main()
