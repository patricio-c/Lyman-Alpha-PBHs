#!/usr/bin/env python
"""
debug_tau.py - Localiza el origen de un tau negativo.

Rehace extract_tau paso a paso para UNA linea de visión, reportando min, max,
NaN e Inf de cada arreglo intermedio, y despues descompone los pixeles con tau
negativo en las contribuciones que los formaron.

Todos los factores del bucle de Voigt son positivos por construccion, asi que
un tau negativo significa que alguno NO lo es: o hay NaN/Inf propagandose, o
un intermedio se vuelve negativo por una via que no previmos. Este script dice
cual.

Uso:
    python debug_tau.py prueba.hdf5 LOS_0002 --treecool TREECOOL_HM12_G_Q
    python debug_tau.py prueba.hdf5 --scan --treecool TREECOOL_HM12_G_Q
"""

import argparse

import h5py
import numpy as np

from ionization import neutral_fraction
from sherwood_los import gamma_hi
from swift_extract import (_cgs_factor, _ray_position, GAMMA_LYA,
                           LAMBDA_LYA_CM, SIGMA_LYA_DV, detect_los_axis,
                           doppler_b, kernel_m4, open_los_file,
                           voigt_hjerting)


def report(name, arr, fatal_if_neg=False):
    a = np.asarray(arr, dtype=np.float64)
    nan = int(np.isnan(a).sum())
    inf = int(np.isinf(a).sum())
    fin = a[np.isfinite(a)]
    mn = fin.min() if fin.size else float("nan")
    mx = fin.max() if fin.size else float("nan")
    neg = int((fin < 0).sum())
    flag = ""
    if nan or inf:
        flag = "  <<< NaN/Inf"
    elif neg and fatal_if_neg:
        flag = "  <<< NEGATIVOS"
    print(f"  {name:<22} min={mn:>12.4e}  max={mx:>12.4e}  "
          f"neg={neg:>7}  nan={nan}  inf={inf}{flag}")
    return nan or inf or (neg and fatal_if_neg)


def debug_one(path, los_name, gamma_HI, npix=2048, X_H=0.76,
              he_state="HeII", n_sigma=8.0, verbose=True):
    meta = open_los_file(path)

    with h5py.File(path, "r") as f:
        g = f[los_name]
        coords_int = g["Coordinates"][:].astype(np.float64)
        rho_c = _cgs_factor(g["Densities"], physical=False)
        rho = g["Densities"][:].astype(np.float64) * rho_c * (1.0 + meta.z) ** 3
        T = g["Temperatures"][:].astype(np.float64) * _cgs_factor(g["Temperatures"])
        vel = g["Velocities"][:].astype(np.float64) * _cgs_factor(g["Velocities"]) / 1.0e5
        hsml_int = g["SmoothingLengths"][:].astype(np.float64)
        mass = g["Masses"][:].astype(np.float64) * _cgs_factor(g["Masses"])
        c_len = _cgs_factor(g["Coordinates"], physical=False)
        axis = detect_los_axis(coords_int, meta.boxsize_int)
        tr = [i for i in range(3) if i != axis]
        ray_int, ray_src = _ray_position(g, coords_int[:, tr],
                                         hsml_int * meta.kernel_gamma,
                                         meta.boxsize_int)

    if verbose:
        print(f"\n=== {los_name}  (eje {'xyz'[axis]}, rayo via {ray_src}, "
              f"{len(rho)} particulas) ===\n")
        print("--- entradas por particula ---")
    bad = False
    if verbose:
        bad |= report("Densities [g/cm3]", rho, True)
        bad |= report("Temperatures [K]", T, True)
        bad |= report("Masses [g]", mass, True)
        bad |= report("SmoothingLengths", hsml_int, True)
        bad |= report("Velocities_par", vel[:, axis])

    to_cm = c_len * meta.a
    x_par = coords_int[:, axis] * to_cm
    box_cm = meta.boxsize_int * to_cm
    Hsup = hsml_int * meta.kernel_gamma * to_cm
    dperp = coords_int[:, tr] - ray_int[None, :]
    dperp -= meta.boxsize_int * np.round(dperp / meta.boxsize_int)
    b_perp = np.hypot(dperp[:, 0], dperp[:, 1]) * to_cm

    n_frac = neutral_fraction(rho * X_H / 1.67262192e-24, T, gamma_HI,
                              X_H=X_H, he_state=he_state)
    mHI = mass * X_H * n_frac / 1.67262192e-24

    if verbose:
        print("\n--- derivados por particula ---")
        bad |= report("b_perp / Hsup", b_perp / Hsup)
        bad |= report("fraccion neutra", n_frac, True)
        bad |= report("atomos HI", mHI, True)

    dR = box_cm / npix
    x_grid = (np.arange(npix) + 0.5) * dR
    dx = x_grid[None, :] - x_par[:, None]
    dx -= box_cm * np.round(dx / box_cm)
    r3 = np.sqrt(b_perp[:, None] ** 2 + dx ** 2)
    w = kernel_m4(r3, Hsup[:, None])

    n_HI = w.T @ mHI
    norm = np.where(n_HI > 0, n_HI, 1.0)
    T_g = np.where(n_HI > 0, (w.T @ (mHI * T)) / norm, 1.0e4)
    v_g = np.where(n_HI > 0, (w.T @ (mHI * vel[:, axis])) / norm, 0.0)

    if verbose:
        print("\n--- campos depositados en la grilla ---")
        bad |= report("kernel w", w, True)
        bad |= report("n_HI [cm-3]", n_HI, True)
        bad |= report("T_g [K]", T_g, True)
        bad |= report("v_g [km/s]", v_g)
        nz = int((n_HI <= 0).sum())
        print(f"  pixeles con n_HI <= 0: {nz} de {npix} ({100*nz/npix:.1f}%)")
        tiny = int(((n_HI > 0) & (n_HI < 1e-300)).sum())
        if tiny:
            print(f"  pixeles con n_HI subnormal (<1e-300): {tiny}  "
                  "<<< puede dar T_g = 0 y b = 0")

    dv = meta.box_kms / npix
    v_hub = (np.arange(npix) + 0.5) * dv
    u = v_hub + v_g
    bth = doppler_b(T_g) / 1.0e5
    a_damp = GAMMA_LYA * LAMBDA_LYA_CM / (4.0 * np.pi * bth * 1.0e5)

    if verbose:
        print("\n--- espacio de velocidades ---")
        bad |= report("b termico [km/s]", bth, True)
        bad |= report("a amortiguamiento", a_damp, True)
        z0 = int((bth <= 0).sum())
        if z0:
            print(f"  b termico <= 0 en {z0} pixeles  <<< division por cero")

    reach = (np.abs(v_g).max() + n_sigma * bth.max()) / dv
    hw = int(min(np.ceil(reach), npix // 2))

    tau = np.zeros(npix, dtype=np.float64)
    worst = np.zeros(npix, dtype=np.float64)
    pref = SIGMA_LYA_DV * dR / np.sqrt(np.pi)
    for off in range(-hw, hw + 1):
        j = (np.arange(npix) - off) % npix
        dvel = v_hub - u[j]
        dvel -= meta.box_kms * np.round(dvel / meta.box_kms)
        H = voigt_hjerting(a_damp[j], dvel / bth[j])
        term = pref * n_HI[j] / (bth[j] * 1.0e5) * H
        tau += term
        worst = np.minimum(worst, term)

    if verbose:
        print(f"\n--- integral de Voigt (banda +/-{hw} pixeles) ---")
        bad |= report("tau", tau, True)
        bad |= report("termino mas negativo", worst, True)

        neg = np.flatnonzero(tau < 0)
        if neg.size:
            print(f"\n  {neg.size} pixeles con tau < 0. Los 5 peores:")
            for i in neg[np.argsort(tau[neg])][:5]:
                print(f"    pixel {i:5d}: tau={tau[i]:+.4e}  n_HI={n_HI[i]:.3e}  "
                      f"T_g={T_g[i]:.3e}  b={bth[i]:.3e}  a={a_damp[i]:.3e}  "
                      f"v_g={v_g[i]:+.2f}")
            print("\n  Si n_HI, T_g y b son positivos y aun asi tau < 0, el")
            print("  problema esta en la aproximacion del perfil para ESOS")
            print("  valores de (a, x). Pasame esta tabla.")
        elif not bad:
            print("\n  Todo positivo en esta LOS.")
    return tau, bad, np.flatnonzero(tau < 0).size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("swift_los")
    ap.add_argument("los_name", nargs="?", default=None)
    ap.add_argument("--treecool", required=True)
    ap.add_argument("--npix", type=int, default=2048)
    ap.add_argument("--scan", action="store_true",
                    help="recorrer todas las LOS y listar cuales fallan")
    args = ap.parse_args()

    meta = open_los_file(args.swift_los)
    G = gamma_hi(args.treecool, meta.z)
    print(f"z = {meta.z:.4f}, Gamma_HI = {G:.4e} s^-1")

    if args.scan:
        print("\nRecorriendo todas las LOS:")
        malas = []
        for nm in meta.los_names:
            try:
                _, bad, nneg = debug_one(args.swift_los, nm, G, args.npix,
                                         verbose=False)
                if nneg:
                    malas.append((nm, nneg))
                    print(f"  {nm}: {nneg} pixeles con tau < 0")
            except Exception as ex:
                print(f"  {nm}: excepcion {ex}")
        print(f"\n{len(malas)} de {len(meta.los_names)} LOS con tau negativo")
        if malas:
            print(f"Ahora mira una en detalle:")
            print(f"  python debug_tau.py {args.swift_los} {malas[0][0]} "
                  f"--treecool {args.treecool}")
        return

    nm = args.los_name or meta.los_names[0]
    debug_one(args.swift_los, nm, G, args.npix)


if __name__ == "__main__":
    main()
