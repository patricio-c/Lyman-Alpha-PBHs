#!/usr/bin/env python
"""
deficit_from_caches.py - Diagnostico del deficit a k bajo, SOLO desde los caches.
Corre en ~1 min, sin re-extraer nada. Uso:
    python deficit_from_caches.py cache_cdm.npz cache_fct.npz --out deficit_diag.png

Tests (todos sobre los tau reales de los caches):
  1. cociente con reescalado a tau_eff comun (lo de la figura)
  2. cociente con A = 1 en ambas  -> aisla el reescalado
  3. cociente con el MISMO A en ambas
  4. coherencia r(k) entre pares de LOS (misma semilla)
  5. bootstrap PAREADO del cociente (error real, no el de LOS independientes)
  6. fracciones de pixeles ultra-transparentes (huella de la conversion QLA)
"""
import argparse
import numpy as np

from sherwood_los import flux_power_1d, rescale_tau, tau_eff


def pk_per_los(F, dv):
    npix = F.shape[1]
    d = F / F.mean() - 1.0
    k = 2.0 * np.pi * np.fft.rfftfreq(npix, d=dv)
    pk = (npix * dv / npix ** 2) * np.abs(np.fft.rfft(d, axis=1)) ** 2
    w2 = np.sinc(k * dv / (2.0 * np.pi)) ** 2
    return k, pk / w2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a"); ap.add_argument("b")
    ap.add_argument("--tau-eff", type=float, default=0.3719)
    ap.add_argument("--out", default="deficit_diag.png")
    args = ap.parse_args()

    A, B = np.load(args.a), np.load(args.b)
    ta = A["tau"].astype(np.float64); tb = B["tau"].astype(np.float64)
    dv = float(A["dv"])
    assert ta.shape == tb.shape, "caches con distinta forma: no hay pares"

    Aa = rescale_tau(ta, args.tau_eff); Ab = rescale_tau(tb, args.tau_eff)
    Am = np.sqrt(Aa * Ab)
    print(f"A(CDM) = {Aa:.5f}   A(FCT) = {Ab:.5f}   cociente = {Ab/Aa:.4f}")

    cases = {
        rf"rescaled, $\tau_{{\rm eff}}$ = {args.tau_eff:.3f}": (Aa, Ab),
        "A = 1 (sin reescalar)": (1.0, 1.0),
        f"mismo A = {Am:.3f}": (Am, Am),
    }
    curves = {}
    for lab, (sa, sb) in cases.items():
        k, pka = pk_per_los(np.exp(-sa * ta), dv)
        _, pkb = pk_per_los(np.exp(-sb * tb), dv)
        curves[lab] = (k, pkb.mean(0) / pka.mean(0), pka, pkb)

    # bootstrap pareado sobre el caso reescalado
    lab0 = list(cases)[0]
    k, r0, pka, pkb = curves[lab0]
    rng = np.random.default_rng(1)
    nlos = ta.shape[0]
    iks = [np.argmin(abs(k - x)) for x in (0.0028, 0.0042, 0.0056, 0.0084)]
    boots = np.empty((400, len(iks)))
    for i in range(400):
        s = rng.integers(0, nlos, nlos)
        boots[i] = pkb[s][:, iks].mean(0) / pka[s][:, iks].mean(0)
    print("\ncociente reescalado, bootstrap PAREADO (misma LOS en num y den):")
    for j, ik in enumerate(iks):
        sig = (1.0 - r0[ik]) / boots[:, j].std()
        print(f"  k = {k[ik]:.4f}: {r0[ik]:.4f} +- {boots[:, j].std():.4f} "
              f"({sig:.1f} sigma debajo de 1)")

    # coherencia (misma semilla => la materia de gran escala es identica)
    Fa = np.exp(-Aa * ta); Fb = np.exp(-Ab * tb)
    fa = np.fft.rfft(Fa / Fa.mean() - 1, axis=1)
    fb = np.fft.rfft(Fb / Fb.mean() - 1, axis=1)
    coh = (fa * np.conj(fb)).real.mean(0) / np.sqrt(
        (np.abs(fa) ** 2).mean(0) * (np.abs(fb) ** 2).mean(0))
    print("\ncoherencia r(k) del flujo entre corridas (pares de LOS):")
    for x in (0.0028, 0.0056, 0.02, 0.05, 0.13):
        ik = np.argmin(abs(k - x))
        print(f"  k = {k[ik]:.4f}: r = {coh[ik]:.3f}   "
              f"sqrt(cociente) = {np.sqrt(r0[ik]):.3f}")

    # huella de la conversion: pixeles ultra-transparentes
    print("\nfraccion de pixeles transparentes (flujo reescalado):")
    for thr in (0.95, 0.98, 0.99):
        fa_, fb_ = (Fa > thr).mean(), (Fb > thr).mean()
        print(f"  F > {thr}: CDM {fa_:.4f}   FCT {fb_:.4f}   "
              f"cociente {fb_/fa_:.2f}")

    # figura
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 7.6), sharex=True,
                                 gridspec_kw={"height_ratios": [2, 1],
                                              "hspace": 0.06})
    m = (k > 0) & (k < 0.36)
    for (lab, (kk, r, _, _)), c in zip(curves.items(), ["C3", "C0", "C2"]):
        ax.plot(kk[m], r[m], color=c, lw=1.6, label=lab)
    ax.axhline(1, color="k", ls="--", lw=1)
    ax.axvspan(1e-3, 0.0318, color="green", alpha=0.07)
    ax.set(xscale="log", ylabel="P(FCT) / P(CDM)", ylim=(0.65, 1.6))
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2, which="both")
    bx.plot(k[m], coh[m], color="C4", lw=1.6)
    bx.axhline(1, color="k", ls="--", lw=1)
    bx.set(xscale="log", xlabel=r"$k$ [s km$^{-1}$]",
           ylabel="coherencia r(k)", ylim=(0, 1.05))
    bx.grid(alpha=0.2, which="both")
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nescrito -> {args.out}")


if __name__ == "__main__":
    main()
