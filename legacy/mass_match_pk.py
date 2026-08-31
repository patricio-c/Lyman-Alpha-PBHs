#!/usr/bin/env python
"""
mass_match_pk.py - Control: cuanto del deficit se explica SOLO por la
fraccion de masa removida, sin nada de fisica de FCT.

Le saca a CDM el gas mas denso hasta dejar la misma masa de gas que
sobrevivio en FCT, y mide el P(k) 3D de ese campo recortado contra el
campo de gas completo de CDM.

    cociente ~ 0.81  -> el deficit se explica solo por la fraccion removida
    cociente ~ 0.92  -> la mitad es la fraccion, la mitad es DONDE remueve FCT
    cociente ~ 1.00  -> la fraccion no importa; es la geometria de la remocion

Uso (la masa objetivo es la que imprimio gas_vs_baryon_pk.py para FCT):

    python mass_match_pk.py \\
        --snap ../lyman/cdm-box-40-1024/cdm-40-m6-lyman_0003/cdm-40-m6-lyman_0003.hdf5 \\
        --keep-mass 6.5094e4 --ngrid 256 --out mass_match.png
"""
import argparse

import h5py
import numpy as np

from gas_vs_baryon_pk import CHUNK, cic, pk_from_grid


def density_threshold(path, keep_mass, nbins=4000):
    """Umbral de densidad tal que la masa con rho <= umbral es keep_mass.

    Histograma en log(rho) para no cargar 10^8 densidades en memoria.
    """
    lo, hi, total = np.inf, -np.inf, 0.0
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        n = g["Masses"].shape[0]
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            r = g["Densities"][s:e].astype(np.float64)
            r = r[r > 0]
            if r.size:
                lo = min(lo, float(np.log10(r.min())))
                hi = max(hi, float(np.log10(r.max())))
            total += float(g["Masses"][s:e].sum())
    edges = np.linspace(lo, hi * 1.000001, nbins + 1)
    hist = np.zeros(nbins)
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        n = g["Masses"].shape[0]
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            r = g["Densities"][s:e].astype(np.float64)
            m = g["Masses"][s:e].astype(np.float64)
            k = r > 0
            idx = np.clip(np.digitize(np.log10(r[k]), edges) - 1, 0, nbins - 1)
            np.add.at(hist, idx, m[k])
    cum = np.cumsum(hist)
    if keep_mass >= total:
        raise SystemExit(f"keep_mass={keep_mass:.4e} >= masa total "
                         f"{total:.4e}: no hay nada que recortar.")
    j = int(np.searchsorted(cum, keep_mass))
    j = min(j, nbins - 1)
    thr = 10.0 ** edges[j + 1]
    print(f"  masa total de gas   = {total:.4e}")
    print(f"  masa objetivo       = {keep_mass:.4e}  "
          f"({100 * keep_mass / total:.1f}% de la de CDM)")
    print(f"  umbral rho          = {thr:.4e} (unidades internas)")
    return thr, total


def gas_iter(path, rho_max=None):
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        n = g["Masses"].shape[0]
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            m = g["Masses"][s:e].astype(np.float64)
            p = g["Coordinates"][s:e].astype(np.float64)
            if rho_max is not None:
                k = g["Densities"][s:e].astype(np.float64) <= rho_max
                m, p = m[k], p[k]
            if m.size:
                yield p, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snap", required=True, help="snapshot de CDM")
    ap.add_argument("--keep-mass", type=float, required=True,
                    help="masa de gas objetivo, en las mismas unidades "
                         "internas que imprime gas_vs_baryon_pk.py")
    ap.add_argument("--ngrid", type=int, default=256)
    ap.add_argument("--out", default="mass_match.png")
    args = ap.parse_args()

    with h5py.File(args.snap, "r") as f:
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])

    print("buscando el umbral de densidad:")
    thr, total = density_threshold(args.snap, args.keep_mass)

    print("\ndepositando:")
    gF, tF, t2F = cic(gas_iter(args.snap), args.ngrid, box)
    print(f"  gas completo  : M = {tF:.4e}")
    gT, tT, t2T = cic(gas_iter(args.snap, thr), args.ngrid, box)
    print(f"  gas recortado : M = {tT:.4e}  "
          f"(objetivo {args.keep_mass:.4e}, error "
          f"{100 * (tT / args.keep_mass - 1):+.2f}%)")

    kF, pF, sF = pk_from_grid(gF, tF, t2F, box, args.ngrid)
    kT, pT, sT = pk_from_grid(gT, tT, t2T, box, args.ngrid)

    print(f"\n{'k [1/Mpc]':>10}{'recortado/completo':>20}")
    for i in range(0, len(kF), max(1, len(kF) // 16)):
        print(f"{kF[i]:10.4f}{pT[i] / pF[i]:20.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(kF, pT / pF, color="C0", lw=1.8,
            label="CDM gas trimmed / CDM gas full")
    ax.axhline(1, color="0.5", ls="--", lw=1)
    ax.axhline(0.808, color="C2", ls=":", lw=1.5,
               label="FCT gas / CDM gas (measured)")
    ax.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]", ylabel="P ratio",
           ylim=(0.6, 1.15),
           title="Control: mass removal alone, z = 3")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2, which="both")
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nescrito -> {args.out}")


if __name__ == "__main__":
    main()
