#!/usr/bin/env python
"""
Donde estan los agujeros? Test de falsacion de la hipotesis del yml.

La hipotesis (range_when_shooting_down_* = [0,40] en caja de 58.7372 Mpc
internos) predice DOS cosas exactas:

  P1. En cada LOS, la coordenada axial de las particulas termina en
      ~40.0 + gamma*h (derrame de kernels de borde), nunca cerca de 58.7.
  P2. El intervalo de pixeles vacios es UNO solo por LOS, ubicado en
      [~40.3, 58.74] Mpc internos = [~27.5, 40.0] Mpc/h, IDENTICO en todas
      las LOS y en los tres ejes.

La hipotesis alternativa (resolucion insuficiente / voids sin cobertura)
predice lo contrario: varios huecos por LOS, en posiciones DISTINTAS en cada
linea, sin relacion con la coordenada 40.

Uso:
    python where_holes.py los_0004.hdf5 --n-los 12
"""

import argparse

import h5py
import numpy as np

from swift_extract import (_cgs_factor, _ray_position, detect_los_axis,
                           kernel_m4, open_los_file)


def holes_of(path, los_name, meta, npix=2048):
    with h5py.File(path, "r") as f:
        g = f[los_name]
        coords = g["Coordinates"][:].astype(np.float64)
        hsml = g["SmoothingLengths"][:].astype(np.float64)
        rho = g["Densities"][:].astype(np.float64)
        mass = g["Masses"][:].astype(np.float64)
        axis = detect_los_axis(coords, meta.boxsize_int)
        tr = [i for i in range(3) if i != axis]
        ray, _ = _ray_position(g, coords[:, tr],
                               hsml * meta.kernel_gamma, meta.boxsize_int)

    box = meta.boxsize_int
    Hsup = hsml * meta.kernel_gamma
    x = coords[:, axis]
    dp = coords[:, tr] - ray[None, :]
    dp -= box * np.round(dp / box)
    b = np.hypot(dp[:, 0], dp[:, 1])

    dR = box / npix
    xg = (np.arange(npix) + 0.5) * dR
    dx = xg[None, :] - x[:, None]
    dx -= box * np.round(dx / box)
    w = kernel_m4(np.sqrt(b[:, None] ** 2 + dx ** 2), Hsup[:, None])
    wsum = w.T @ (mass / rho)

    empty = wsum <= 0.0
    ivals = []
    if empty.any() and not empty.all():
        i0 = int(np.argmin(empty))
        e = np.roll(empty, -i0)
        edges = np.flatnonzero(np.diff(np.concatenate(
            [[0], e.view(np.int8), [0]])))
        for s, t in zip(edges[::2], edges[1::2]):
            lo = ((s + i0) % npix) * dR
            hi = lo + (t - s) * dR
            ivals.append((lo, hi if hi <= box else hi - box))
    return axis, float(x.min()), float(x.max()), ivals, float(Hsup.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("swift_los")
    ap.add_argument("--n-los", type=int, default=12)
    ap.add_argument("--npix", type=int, default=2048)
    args = ap.parse_args()

    meta = open_los_file(args.swift_los)
    h = meta.h
    box = meta.boxsize_int
    print(f"\nCaja = {box:.4f} Mpc internos = {box * h:.2f} Mpc/h")
    print(f"Prediccion del yml: particulas hasta ~40.0+H, "
          f"agujero unico en [~40.3, {box:.2f}] Mpc "
          f"= [{40.3 * h:.2f}, {box * h:.2f}] Mpc/h\n")

    hdr = (f"{'LOS':<10}{'eje':>4}{'x_min':>9}{'x_max':>9}"
           f"{'agujeros [Mpc internos] (y en Mpc/h)':>45}")
    print(hdr)
    starts, ends, nholes = [], [], []
    for nm in meta.los_names[: args.n_los]:
        axis, xmin, xmax, ivals, Hmax = holes_of(
            args.swift_los, nm, meta, args.npix)
        txt = "  ".join(f"[{a:.2f},{b:.2f}] ([{a*h:.1f},{b*h:.1f}]h)"
                        for a, b in ivals) or "ninguno"
        print(f"{nm:<10}{'xyz'[axis]:>4}{xmin:>9.2f}{xmax:>9.2f}   {txt}")
        nholes.append(len(ivals))
        for a, b in ivals:
            starts.append(a)
            ends.append(b)

    print("\nResumen:")
    print(f"  agujeros por LOS: {sorted(set(nholes))}  "
          f"(P2 exige exactamente [1])")
    if starts:
        print(f"  inicio del agujero: media = {np.mean(starts):.3f}, "
              f"desvio = {np.std(starts):.3f} Mpc  "
              f"(P2 exige ~40.3 con desvio ~0.1, el derrame de H)")
        print(f"  fin del agujero:    media = {np.mean(ends):.3f}, "
              f"desvio = {np.std(ends):.3f} Mpc  "
              f"(P2 exige ~{box:.2f} con desvio ~0)")
    print(f"  x_max de particulas por LOS deberia rondar 40.0-40.5 (P1).")
    print("\n  Si inicio~40.3 y fin~58.7 con desvios chicos en TODAS las")
    print("  LOS y los tres ejes: hipotesis del yml CONFIRMADA.")
    print("  Si los agujeros caen en posiciones distintas por LOS: refutada,")
    print("  y volvemos a sospechar de resolucion o del escritor de LOS.")


if __name__ == "__main__":
    main()
