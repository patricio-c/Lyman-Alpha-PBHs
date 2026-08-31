#!/usr/bin/env python
"""
removal_curve_pk.py - Curva de "masa removida vs potencia perdida" para CDM,
para ubicar el punto de FCT sobre ella.

Tres modos de remocion, todos aplicados al gas de CDM:

  threshold    corte duro: se va TODO el gas con Delta > valor.
               Es lo que hace QLA en un instante dado.

  random-above se va una fraccion aleatoria `valor` de las particulas con
               Delta > --delta-min. Se parece mas a QLA acumulado: las
               particulas cruzaron el umbral en momentos distintos y el
               campo de hoy no es un corte nitido.

  random-all   se va una fraccion aleatoria en TODAS partes.
               Control nulo: no debe cambiar P(k), solo el shot noise.
               Si este modo mueve la curva, hay un bug.

Un solo pase por el archivo: deposita todas las variantes a la vez.

    python removal_curve_pk.py --snap <cdm_snapshot> \\
        --mode random-above --delta-min 10 \\
        --values 0.2 0.4 0.6 0.8 0.95 --ngrid 256 --out curve_random.png
"""
import argparse

import h5py
import numpy as np

from gas_vs_baryon_pk import CHUNK, pk_from_grid

RHO_CRIT_INTERNAL = 27.754      # 2.7754e11 Msun/Mpc^3 en unidades de 1e10 Msun


def cic_add(grid, pos, m, ngrid, box):
    f = (pos / box) * ngrid
    f -= np.floor(f / ngrid) * ngrid
    i0 = np.floor(f).astype(np.int64)
    d = f - i0
    i0 %= ngrid
    i1 = (i0 + 1) % ngrid
    flat = grid.reshape(-1)
    n2 = ngrid * ngrid
    for a in (0, 1):
        wx = d[:, 0] if a else 1.0 - d[:, 0]
        ix = i1[:, 0] if a else i0[:, 0]
        for b in (0, 1):
            wy = d[:, 1] if b else 1.0 - d[:, 1]
            iy = i1[:, 1] if b else i0[:, 1]
            wxy = wx * wy
            base = ix * n2 + iy * ngrid
            for c in (0, 1):
                wz = d[:, 2] if c else 1.0 - d[:, 2]
                iz = i1[:, 2] if c else i0[:, 2]
                flat += np.bincount(base + iz, weights=m * wxy * wz,
                                    minlength=flat.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snap", required=True)
    ap.add_argument("--mode", required=True,
                    choices=["threshold", "random-above", "random-all"])
    ap.add_argument("--values", type=float, nargs="+", required=True,
                    help="Delta de corte (threshold) o fraccion removida "
                         "(random-*)")
    ap.add_argument("--delta-min", type=float, default=10.0,
                    help="solo para random-above: densidad minima de las "
                         "candidatas a remover")
    ap.add_argument("--ngrid", type=int, default=256)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--ref-ratio", type=float, default=0.808,
                    help="cociente de P(k) medido de FCT/CDM a k bajo")
    ap.add_argument("--ref-massfrac", type=float, default=0.404,
                    help="fraccion del gas de CDM que le falta a FCT")
    ap.add_argument("--out", default="removal_curve.png")
    args = ap.parse_args()

    with h5py.File(args.snap, "r") as f:
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        c = f["Cosmology"].attrs
        h = float(np.atleast_1d(c["h"])[0])
        ob = float(np.atleast_1d(c["Omega_b"])[0])
        npart = f["PartType0"]["Masses"].shape[0]
    rho_mean = ob * RHO_CRIT_INTERNAL * h ** 2
    print(f"densidad bariónica media comovil = {rho_mean:.4f} (internas)")
    print(f"{npart} particulas de gas, caja {box:.4f} Mpc\n")

    nv = len(args.values)
    grids = [np.zeros((args.ngrid,) * 3) for _ in range(nv + 1)]
    tot = np.zeros(nv + 1)
    tot2 = np.zeros(nv + 1)
    rng = np.random.default_rng(args.seed)

    with h5py.File(args.snap, "r") as f:
        g = f["PartType0"]
        for s in range(0, npart, CHUNK):
            e = min(s + CHUNK, npart)
            m = g["Masses"][s:e].astype(np.float64)
            p = g["Coordinates"][s:e].astype(np.float64)
            delta = g["Densities"][s:e].astype(np.float64) / rho_mean
            u = rng.random(m.size)
            cic_add(grids[0], p, m, args.ngrid, box)
            tot[0] += m.sum(); tot2[0] += (m ** 2).sum()
            for j, v in enumerate(args.values):
                if args.mode == "threshold":
                    keep = delta <= v
                elif args.mode == "random-above":
                    keep = ~((delta > args.delta_min) & (u < v))
                else:
                    keep = u >= v
                if keep.any():
                    cic_add(grids[j + 1], p[keep], m[keep], args.ngrid, box)
                    tot[j + 1] += m[keep].sum()
                    tot2[j + 1] += (m[keep] ** 2).sum()
            if s % (CHUNK * 8) == 0:
                print(f"  {e}/{npart}")

    res = [pk_from_grid(grids[i], tot[i], tot2[i], box, args.ngrid)
           for i in range(nv + 1)]
    k, p_full = res[0][0], res[0][1]
    lowk = k < 0.7                       # zona donde el cociente es plano

    label = "Delta_cut" if args.mode == "threshold" else "frac_rand"
    print(f"\n{label:>11}{'M removida':>12}{'P/P_full (k<0.7)':>18}")
    curve = []
    for j, v in enumerate(args.values):
        removed = 1.0 - tot[j + 1] / tot[0]
        r = float(np.mean(res[j + 1][1][lowk] / p_full[lowk]))
        curve.append((removed, r))
        print(f"{v:11.4g}{removed:12.4f}{r:18.4f}")
    print(f"\npunto de FCT: masa removida {args.ref_massfrac:.3f}, "
          f"P/P_full {args.ref_ratio:.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    xs = [c[0] for c in curve]; ys = [c[1] for c in curve]
    ax.plot(xs, ys, "o-", color="C0", lw=1.8,
            label=f"CDM, {args.mode}")
    ax.plot(args.ref_massfrac, args.ref_ratio, "*", color="C3", ms=18,
            label="FCT (measured)")
    ax.axhline(1, color="0.5", ls="--", lw=1)
    ax.set(xlabel="fraction of CDM gas mass removed",
           ylabel=r"$P_{\rm gas}$ / $P_{\rm gas,full}$  ($k<0.7$ Mpc$^{-1}$)",
           title="Cost in large-scale power of removing gas, z = 3")
    ax.legend(frameon=False); ax.grid(alpha=0.2)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nescrito -> {args.out}")
    print("\nSi FCT cae SOBRE la curva -> presupuesto de masa.")
    print("Si cae POR ENCIMA -> FCT remueve de lugares mas baratos.")


if __name__ == "__main__":
    main()
