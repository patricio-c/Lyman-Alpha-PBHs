#!/usr/bin/env python
"""
Diagnostico de la extraccion. Corre esto ANTES de volver a la pipeline de
validacion. Responde cuatro preguntas concretas:

  A. La geometria cierra? (b <= H para todas las particulas, todas las LOS)
  B. El muestreo SPH es sano? (suma de pesos ~ 1 en todos los pixeles)
  C. La DISTRIBUCION de tau se parece a la de Sherwood, no solo su media?
     El sintoma del bug era una PDF bimodal; la media (tau_eff) era casi
     correcta y por eso engañaba.
  D. Cuanto vale el factor A para llegar al tau_eff de Sherwood?
     Sano: A ~ 0.8-1.3. Enfermo: A >> 1 o << 1.

Uso:
    python diagnose_extraction.py los_0004.hdf5 tauH1_lya_z4_2.dat \
        --treecool TREECOOL_HM12_G_Q --max-los 30
"""

import argparse

import numpy as np

from sherwood_los import gamma_hi, read_sherwood_spectra, rescale_tau, tau_eff
from swift_extract import extract_all, open_los_file


def pdf_percentiles(t, label):
    q = np.percentile(t, [1, 5, 25, 50, 75, 95, 99])
    print(f"  {label:<10s} p1={q[0]:.3e} p5={q[1]:.3e} p25={q[2]:.3e} "
          f"p50={q[3]:.3e} p75={q[4]:.3e} p95={q[5]:.3e} p99={q[6]:.3e}")
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("swift_los")
    ap.add_argument("sherwood")
    ap.add_argument("--treecool", required=True)
    ap.add_argument("--npix", type=int, default=2048)
    ap.add_argument("--max-los", type=int, default=30)
    args = ap.parse_args()

    meta = open_los_file(args.swift_los)
    G = gamma_hi(args.treecool, meta.z)
    sh = read_sherwood_spectra(args.sherwood)
    tau_sh = np.asarray(sh.tau, dtype=np.float64)

    print(f"\nExtrayendo {args.max_los} LOS con diagnosticos ...")
    tau, dv, _, diags = extract_all(args.swift_los, args.npix, G,
                                    max_los=args.max_los, collect_diag=True)
    tau = np.asarray(tau, dtype=np.float64)

    # --- A: geometria ------------------------------------------------------
    mb = max(d["max_b_over_H"] for d in diags)
    fb = max(d["frac_b_gt_H"] for d in diags)
    srcs = {d["ray_source"] for d in diags}
    print(f"\nA) Geometria: rayo via {srcs}")
    print(f"   max(b/H) sobre todas las LOS = {mb:.4f}  (debe ser <= 1)")
    print(f"   peor fraccion con b > H     = {100 * fb:.2f}%  (debe ser ~0)")

    # --- B: muestreo -------------------------------------------------------
    wmin = min(d["wsum_min"] for d in diags)
    wmed = float(np.median([d["wsum_med"] for d in diags]))
    wmax = max(d["wsum_max"] for d in diags)
    print(f"\nB) Suma de pesos SPH por pixel (ideal ~ 1):")
    print(f"   min = {wmin:.3f}   mediana = {wmed:.3f}   max = {wmax:.3f}")
    print(f"   min << 1 => pixeles submuestreados (voids sin cobertura);")
    print(f"   sistematicamente > 1 => doble conteo o error de unidades.")

    # --- C: distribucion de tau -------------------------------------------
    print(f"\nC) Percentiles de tau (comparar fila a fila):")
    q_tu = pdf_percentiles(tau.ravel(), "tuyo")
    q_sh = pdf_percentiles(tau_sh.ravel(), "Sherwood")
    print(f"   cociente p50 tuyo/Sherwood = {q_tu[3] / q_sh[3]:.3f}")
    print(f"   fraccion de pixeles con tau < 0.01: "
          f"tuyo {100 * (tau < 0.01).mean():.2f}%  "
          f"Sherwood {100 * (tau_sh < 0.01).mean():.2f}%")

    # --- D: curva tau_eff(A) ----------------------------------------------
    target = tau_eff(tau_sh)
    print(f"\nD) tau_eff nominal tuyo = {tau_eff(tau):.4f}  "
          f"(Sherwood: {target:.4f})")
    print(f"   tau_eff(A) - si es plana, el campo es binario:")
    for A in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]:
        print(f"     A = {A:5.2f}  ->  tau_eff = {tau_eff(tau, A):.4f}")
    try:
        A = rescale_tau(tau, target)
        verdict = "SANO" if 0.5 < A < 2.0 else "REVISAR: sigue habiendo un problema"
        print(f"   A para igualar Sherwood = {A:.4f}   [{verdict}]")
    except ValueError as e:
        print(f"   rescale fallo: {e}")


if __name__ == "__main__":
    main()
