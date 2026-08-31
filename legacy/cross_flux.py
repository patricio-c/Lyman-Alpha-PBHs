#!/usr/bin/env python
"""
cross_flux.py - Espectro cruzado entre los campos de flujo de las dos corridas.

El problema
-----------
El cociente P1D da r(k=0.005) = 0.84 mientras el P(k) 3D de materia da 1.000
exacto a esa escala. Ya se descarto cosmologia, reescalado de tau_eff,
muestreo SPH y sistemas saturados. Queda la pregunta de fondo: ese 16% que
falta, es una DIFERENCIA DE AMPLITUD sobre los mismos modos, o es que las dos
simulaciones estan trazando estructura distinta?

Son dos afirmaciones muy distintas para defender en una charla:

  "el bosque de la FCT tiene un bias 8% menor porque QLA le saco la mitad
   del gas"                                                    <- resultado

  "a escalas grandes mis dos corridas no se parecen"            <- bug

El test
-------
Como las LOS son las mismas (mismo seed, mismas posiciones, mismo eje), los
dos campos de flujo son dos mediciones del MISMO volumen. Entonces se puede
calcular el espectro cruzado y descomponer:

    r_c(k) = P_ab / sqrt(P_aa P_bb)        coeficiente de correlacion
    b(k)   = P_ab / P_aa                    bias relativo (regresion de b en a)
    P_bb   = b^2 P_aa  +  P_est             con P_est = P_bb (1 - r_c^2)

Lectura:

  r_c -> 1 a k chico  =>  los mismos modos, misma fase. La diferencia es
                          PURA AMPLITUD. b(k) es el bias relativo y es un
                          numero fisico y citable.
  r_c < 1 a k chico   =>  hay estructura decorrelacionada. P_est te dice
                          cuanta. Si P_est es grande, algo esta mal.

Ademas hay un chequeo de consistencia automatico: si la diferencia es pura
amplitud, debe cumplirse b(k)^2 = P_bb/P_aa. La brecha entre esas dos curvas
ES la estocasticidad, y verla es mas informativo que el cociente solo.

Normalizacion
-------------
Los P se calculan aca adentro, no con flux_power_1d, para que las tres
combinaciones (aa, ab, bb) compartan exactamente la misma convencion. La
normalizacion absoluta puede diferir de flux_power_1d; los cocientes no. Y
la ventana de pixel se cancela en todos los cocientes, asi que no se
deconvoluciona.

Ejemplo
-------
    python cross_flux.py --a cache_cdm.npz --b cache_fct.npz \\
        --tau-eff 0.3719 --nboot 200 --kt 10 --out cross.png
"""

import argparse
import sys

import numpy as np

from sherwood_los import rescale_tau

LAMBDA_LYA = 1215.6701
C_KMS = 2.99792458e5


def load_cache(path):
    d = np.load(path, allow_pickle=True)
    return dict(tau=np.asarray(d["tau"], dtype=np.float64),
                dv=float(d["dv"]), z=float(d["z"]),
                box_kms=float(d["box_kms"]), h=float(d["h"]),
                hz=float(d["hz"]))


def delta_fft(F, dv):
    """rfft de delta_F = F/<F> - 1, con <F> global. Devuelve (k, fk)."""
    d = F / F.mean() - 1.0
    fk = np.fft.rfft(d, axis=1) * dv
    k = 2.0 * np.pi * np.fft.rfftfreq(F.shape[1], d=dv)
    return k, fk


def spectra(fa, fb, npix, dv):
    """Auto y cruzado POR LINEA DE VISION, sin promediar todavia.

    Se devuelve sin promediar para poder bootstrapear sobre las LOS despues,
    que es la unica fuente de error que tiene sentido aca: los modos dentro
    de una misma LOS no son independientes entre corridas, pero las LOS si
    lo son entre si (aproximadamente, ignorando que comparten la caja).
    """
    n = npix * dv
    paa = (fa.real ** 2 + fa.imag ** 2) / n
    pbb = (fb.real ** 2 + fb.imag ** 2) / n
    pab = (fa.real * fb.real + fa.imag * fb.imag) / n     # Re(fa conj(fb))
    return paa, pbb, pab


def logbin(k, arrs, nbin):
    """Rebinea en log k. arrs es una lista de arrays con la forma de k."""
    m = k > 0
    edges = np.logspace(np.log10(k[m][0] * 0.999),
                        np.log10(k[m][-1] * 1.001), nbin + 1)
    idx = np.digitize(k, edges) - 1
    kb, out = [], [[] for _ in arrs]
    for b in range(nbin):
        s = (idx == b) & m
        if s.sum() == 0:
            continue
        kb.append(k[s].mean())
        for j, a in enumerate(arrs):
            out[j].append(a[s].mean())
    return np.asarray(kb), [np.asarray(o) for o in out]


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", required=True)
    p.add_argument("--b", required=True)
    p.add_argument("--tau-eff", type=float, required=True)
    p.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    p.add_argument("--nbin", type=int, default=25)
    p.add_argument("--nboot", type=int, default=200,
                   help="remuestreos sobre las LOS para las barras de error")
    p.add_argument("--seed", type=int, default=1234)
    p.add_argument("--kt", type=float, default=None,
                   help="escala de ruptura en Mpc^-1 comovil")
    p.add_argument("--xlim", nargs=2, type=float, default=None)
    p.add_argument("--table", default=None)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    A, B = load_cache(args.a), load_cache(args.b)
    la, lb = args.labels
    if A["tau"].shape != B["tau"].shape:
        sys.exit(f"Los caches tienen distinta forma: {A['tau'].shape} vs "
                 f"{B['tau'].shape}. El cruzado necesita las MISMAS LOS.")
    if abs(A["dv"] / B["dv"] - 1.0) > 1e-3:
        sys.exit(f"dv distinto: {A['dv']:.5f} vs {B['dv']:.5f} km/s. Los "
                 "pixeles no corresponden al mismo volumen y el cruzado no "
                 "significa nada. Revisa que las dos corridas tengan la misma "
                 "caja y el mismo H(z).")

    nlos, npix = A["tau"].shape
    dv, z = A["dv"], A["z"]
    to_mpc = A["hz"] / (1.0 + z)
    to_hmpc = to_mpc / A["h"]

    aA = rescale_tau(A["tau"], args.tau_eff)
    aB = rescale_tau(B["tau"], args.tau_eff)
    print(f"z = {z:.2f}   {nlos} LOS x {npix} px   dv = {dv:.4f} km/s")
    print(f"A({la}) = {aA:.5f}   A({lb}) = {aB:.5f}")

    k, fa = delta_fft(np.exp(-aA * A["tau"]), dv)
    _, fb = delta_fft(np.exp(-aB * B["tau"]), dv)
    paa, pbb, pab = spectra(fa, fb, npix, dv)

    def derived(sel):
        Paa, Pbb, Pab = paa[sel].mean(0), pbb[sel].mean(0), pab[sel].mean(0)
        kb, (Paa, Pbb, Pab) = logbin(k, [Paa, Pbb, Pab], args.nbin)
        with np.errstate(divide="ignore", invalid="ignore"):
            rc = Pab / np.sqrt(Paa * Pbb)
            bias = Pab / Paa
            amp = np.sqrt(Pbb / Paa)
            stoch = Pbb * (1.0 - rc ** 2) / Paa
        return kb, rc, bias, amp, stoch

    kb, rc, bias, amp, stoch = derived(np.arange(nlos))

    # bootstrap sobre LOS
    rng = np.random.default_rng(args.seed)
    boot = np.empty((args.nboot, 4, kb.size))
    for i in range(args.nboot):
        sel = rng.integers(0, nlos, nlos)
        _, boot[i, 0], boot[i, 1], boot[i, 2], boot[i, 3] = derived(sel)
    lo, hi = np.percentile(boot, [16, 84], axis=0)

    # ---- figura ----------------------------------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), sharex=True,
                             gridspec_kw={"hspace": 0.08})
    R_z = C_KMS * 0.8 / ((1.0 + z) * LAMBDA_LYA)
    for ax in axes:
        ax.axvspan(1.0e-3, 0.5 * np.pi / R_z, color="green", alpha=0.07,
                   zorder=0)
        ax.axvline(2.0 * np.pi / A["box_kms"], color="grey", ls=":", lw=1.2)
        if args.kt:
            ax.axvline(args.kt / to_mpc, color="C3", ls=":", lw=1.2)
        ax.grid(alpha=0.2, which="both")
        ax.set_xscale("log")

    ax = axes[0]
    ax.fill_between(kb, lo[0], hi[0], color="C2", alpha=0.3)
    ax.plot(kb, rc, color="C2", lw=2.0)
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.set(ylabel=r"$r_c(k) = P_{ab}/\sqrt{P_{aa}P_{bb}}$")
    ax.set_title(rf"cross-spectrum {la} $\times$ {lb},  $z$ = {z:.1f},  "
                 rf"$\tau_{{\rm eff}}$ = {args.tau_eff:.4f}", fontsize=10)

    ax = axes[1]
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.fill_between(kb, lo[1], hi[1], color="C0", alpha=0.3)
    ax.plot(kb, bias, color="C0", lw=2.0, label=r"$b = P_{ab}/P_{aa}$")
    ax.fill_between(kb, lo[2], hi[2], color="C3", alpha=0.25)
    ax.plot(kb, amp, color="C3", lw=1.6, ls="--",
            label=r"$\sqrt{P_{bb}/P_{aa}}$")
    ax.set(ylabel="relative amplitude")
    ax.legend(frameon=False, fontsize=9, loc="upper left")

    ax = axes[2]
    ax.fill_between(kb, lo[3], hi[3], color="C4", alpha=0.3)
    ax.plot(kb, stoch, color="C4", lw=2.0)
    ax.axhline(0.0, color="k", ls="--", lw=1)
    ax.set(yscale="log", ylabel=r"$P_{\rm stoch}/P_{aa}$",
           xlabel=r"$k$  [s km$^{-1}$]")
    if args.xlim:
        ax.set_xlim(*args.xlim)

    tx = axes[0].secondary_xaxis("top", functions=(lambda x: x * to_hmpc,
                                                   lambda x: x / to_hmpc))
    tx.set_xlabel(r"$k$  [$h$ Mpc$^{-1}$]")

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"escrito -> {args.out}")

    # ---- numeros ---------------------------------------------------------
    print(f"\n{'k[s/km]':>10}{'k[1/Mpc]':>10}{'r_c':>9}{'b':>9}"
          f"{'sqrt(Pbb/Paa)':>15}{'b^2':>9}{'Pbb/Paa':>10}{'stoch':>10}")
    for kk in [0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
        if kk < kb[0] or kk > kb[-1]:
            continue
        f = lambda a: np.interp(kk, kb, a)
        print(f"{kk:10.4f}{kk*to_mpc:10.3f}{f(rc):9.4f}{f(bias):9.4f}"
              f"{f(amp):15.4f}{f(bias)**2:9.4f}{f(amp)**2:10.4f}"
              f"{f(stoch):10.2e}")

    i0 = np.argmin(np.abs(kb - 0.005))
    print(f"\nEn k = {kb[i0]:.4f} s/km ({kb[i0]*to_mpc:.3f} Mpc^-1):")
    print(f"  r_c = {rc[i0]:.4f}  [{lo[0][i0]:.4f}, {hi[0][i0]:.4f}]")
    print(f"  bias relativo = {bias[i0]:.4f} "
          f"[{lo[1][i0]:.4f}, {hi[1][i0]:.4f}]")
    print(f"  estocasticidad = {100*stoch[i0]:.2f}% de P_aa")
    if rc[i0] > 0.98:
        print(f"\n  => r_c ~ 1: mismos modos, misma fase. La diferencia es")
        print(f"     PURA AMPLITUD. El bosque de {lb} tiene un bias")
        print(f"     {100*(1-bias[i0]):.1f}% menor que el de {la} a esa escala.")
        print(f"     Eso es un resultado fisico, no un problema de pipeline.")
    else:
        print(f"\n  => r_c = {rc[i0]:.3f} < 1: hay estructura decorrelacionada.")
        print(f"     No es solo amplitud. Revisar antes de citar el cociente.")

    if args.table:
        np.savetxt(args.table,
                   np.column_stack([kb, rc, bias, amp, stoch]),
                   header=f"z={z:.4f} tau_eff={args.tau_eff:.6f} "
                          f"A_{la}={aA:.6f} A_{lb}={aB:.6f}\n"
                          "k[s/km]  r_c  bias  sqrt(Pbb/Paa)  stoch/Paa")
        print(f"escrito -> {args.table}")


if __name__ == "__main__":
    main()
