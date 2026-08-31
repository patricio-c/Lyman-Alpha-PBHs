#!/usr/bin/env python
"""
PIPELINE A. Validacion: reproducir el P1D de Sherwood desde tus propias LOS.

Objetivo: NO es hacer ciencia, es comprobar que tu cadena esta sana. Si tu
CDM no reproduce Sherwood dentro de lo esperable por resolucion y varianza
cosmica, cualquier cosa que digas despues sobre materia oscura primordial es
ruido de pipeline.

Que compara y que no:

  - Ambos se llevan al MISMO tau_eff. Sin esto la comparacion no significa
    nada: el P1D es dominado por el flujo medio.
  - Ambos pasan por el MISMO post-proceso instrumental de Bolton y
    colaboradores 2017: perfil instrumental de 7 km/s, rebin a 3 km/s.
  - Ambos usan el MISMO estimador y la misma deconvolucion de ventana de pixel.
  - Se comparan solo en el rango de k comun a las dos cajas.

Diferencias que vas a ver y que NO son errores:
  - Tu caja da 5080 km/s y la de Sherwood 5102 km/s, por Omega_m distinto.
    Los k no caen en los mismos valores; se interpola.
  - Tu resolucion de masa es mucho peor que la de Sherwood. Esperá deficit de
    potencia a k grande. Eso es fisica de resolucion, no un bug.

Uso:
    python pipeline_sherwood.py los_0004.hdf5 tauH1_lya_z4_2.dat \
        --treecool TREECOOL_HM12_G_Q --npix 2048 --max-los 200
"""

import argparse

import numpy as np

from sherwood_los import (flux_power_1d, gamma_hi, read_sherwood_spectra,
                          rescale_tau, tau_eff)
from sherwood_postprocess import convolve_lsf, rebin_velocity
from swift_extract import extract_all


def process(tau, dv, target_tau_eff, fwhm=7.0, dv_out=3.0):
    """Reescalado + post-proceso instrumental de Sherwood. Devuelve (F, dv)."""
    A = rescale_tau(tau, target_tau_eff)
    F = np.exp(-A * np.asarray(tau, dtype=np.float64))
    F = convolve_lsf(F, dv, fwhm)
    F, dv2 = rebin_velocity(F, dv, dv_out)
    return F, dv2, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("swift_los", help="archivo HDF5 de LOS de SWIFT")
    ap.add_argument("sherwood", help="binario de espectros de Sherwood al mismo z")
    ap.add_argument("--treecool", required=True)
    ap.add_argument("--npix", type=int, default=2048)
    ap.add_argument("--max-los", type=int, default=None,
                    help="limitar el numero de LOS (util para una prueba rapida)")
    ap.add_argument("--tau-eff", type=float, default=None,
                    help="tau_eff comun. Por defecto usa el nominal de Sherwood.")
    ap.add_argument("--out", default="validacion_sherwood.txt")
    args = ap.parse_args()

    # --- Sherwood ---------------------------------------------------------
    sh = read_sherwood_spectra(args.sherwood)
    tau_sh = np.asarray(sh.tau, dtype=np.float64)
    target = args.tau_eff if args.tau_eff is not None else tau_eff(tau_sh)

    print(f"\nSherwood: z={sh.z:.3f}  nlos={sh.nlos}  npix={sh.npix}")
    print(f"  caja = {sh.box_kms:.2f} km/s, dv = {sh.dv:.4f} km/s")
    print(f"  tau_eff objetivo (comun a ambos) = {target:.5f}")

    # --- SWIFT ------------------------------------------------------------
    print(f"\nExtrayendo tau de {args.swift_los} ...")
    from swift_extract import open_los_file
    meta = open_los_file(args.swift_los)
    G = gamma_hi(args.treecool, meta.z)
    print(f"  z={meta.z:.3f}  h={meta.h:.4f}  Omega_m={meta.omega_m:.5f}")
    print(f"  caja = {meta.box_comoving_mpc:.4f} Mpc com = "
          f"{meta.box_comoving_mpc * meta.h:.3f} Mpc/h")
    print(f"  H(z) = {meta.hz:.3f} km/s/Mpc  ->  caja = {meta.box_kms:.2f} km/s")
    print(f"  dv = {meta.dv(args.npix):.4f} km/s")
    print(f"  Gamma_HI = {G:.4e} s^-1")
    if abs(meta.z - sh.z) > 0.02:
        raise SystemExit(f"Redshifts distintos: SWIFT z={meta.z}, Sherwood z={sh.z}")

    tau_sw, dv_sw, _ = extract_all(args.swift_los, args.npix, G,
                                   max_los=args.max_los)
    print(f"  tau_eff nominal de tu corrida = {tau_eff(tau_sw):.5f}")

    # --- mismo post-proceso para los dos ----------------------------------
    F_sh, dvo_sh, A_sh = process(tau_sh, sh.dv, target)
    F_sw, dvo_sw, A_sw = process(tau_sw, dv_sw, target)
    print(f"\n  factor de reescalado Sherwood A = {A_sh:.5f}")
    print(f"  factor de reescalado tuyo     A = {A_sw:.5f}")

    k_sh, p_sh, e_sh = flux_power_1d(F_sh, dvo_sh)
    k_sw, p_sw, e_sw = flux_power_1d(F_sw, dvo_sw)

    # --- rango comun e interpolacion --------------------------------------
    kmin = max(k_sh[1], k_sw[1])
    kmax = min(k_sh[-1], k_sw[-1])
    m = (k_sh >= kmin) & (k_sh <= kmax)
    kk = k_sh[m]
    p_sw_i = np.interp(kk, k_sw, p_sw)
    ratio = p_sw_i / p_sh[m]

    print(f"\nCociente tuyo/Sherwood, ambos a tau_eff={target:.4f}:")
    print(f"  {'k [s/km]':>10}  {'kP/pi tuyo':>12}  {'kP/pi Sherwood':>15}  {'cociente':>9}")
    for i in np.unique(np.geomspace(1, kk.size - 1, 14).astype(int)):
        print(f"  {kk[i]:10.5f}  {kk[i]*p_sw_i[i]/np.pi:12.4e}  "
              f"{kk[i]*p_sh[m][i]/np.pi:15.4e}  {ratio[i]:9.4f}")

    lo = (kk > 0.002) & (kk < 0.01)
    hi = (kk > 0.05) & (kk < 0.2)
    if lo.any():
        print(f"\n  cociente medio en 0.002 < k < 0.01  : {ratio[lo].mean():.4f}")
    if hi.any():
        print(f"  cociente medio en 0.05  < k < 0.2   : {ratio[hi].mean():.4f}")
    print("\n  Interpretacion: a k chico el cociente deberia dar cerca de 1 salvo")
    print("  varianza cosmica. Un deficit creciente a k grande es tu resolucion")
    print("  de masa, no un error. Un exceso a k grande si es sospechoso.")

    np.savetxt(args.out, np.column_stack([kk, p_sw_i, p_sh[m], ratio]),
               header=(f"z={sh.z:.4f} tau_eff={target:.6f} "
                       f"dv_tuyo={dvo_sw:.5f} dv_sherwood={dvo_sh:.5f}\n"
                       "k[s/km]  P1D_tuyo  P1D_sherwood  cociente"))
    print(f"\n  escrito -> {args.out}")


if __name__ == "__main__":
    main()
