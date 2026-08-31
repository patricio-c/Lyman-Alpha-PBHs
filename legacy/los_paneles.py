#!/usr/bin/env python
"""
los_panels.py - La figura clasica del bosque: una linea de vision, tres paneles.

Delta(v), T(v) y F(v) sobre el mismo eje, CDM en azul y FCT en rojo. Como las
dos corridas comparten seed y posiciones de rayo, la LOS numero i es
LITERALMENTE la misma en las dos: se ven los mismos picos y se ve cual le
falta a cual. Es la figura que cuenta el argumento de Bolton et al. 2017 sin
tener que explicarlo.

tau NO se recalcula: sale del cache, junto con el factor de reescalado A, asi
que el flujo del panel de abajo es exactamente el que entra en el P1D. Lo
unico que hay que depositar de nuevo es Delta, T y v_par, que el cache no
guarda. Eso cuesta un par de segundos por LOS.

Subcomandos
-----------
  pick   Recorre los caches y rankea las LOS por cuanto se diferencian, con
         varios criterios. El util para la charla es --criterion lost-sat:
         busca la LOS donde la CDM tiene un sistema saturado que la FCT
         perdio. Esa es la figura.

  plot   Dibuja los paneles para un indice dado.

Ejemplos
--------
    python los_panels.py pick --a cache_cdm.npz --b cache_fct.npz \\
        --tau-eff 0.3719 --criterion lost-sat --top 10

    python los_panels.py plot --a cache_cdm.npz --b cache_fct.npz \\
        --los-a regen/cdm40.hdf5 --los-b regen/fct40.hdf5 \\
        --tau-eff 0.3719 --index 1234 --out los_1234.png
"""

import argparse
import sys

import h5py
import numpy as np

from sherwood_los import rescale_tau
from swift_extract import (MPC_CM, _cgs_factor, _ray_position, detect_los_axis,
                           kernel_m4, open_los_file)

RHO_CRIT_H2 = 1.87847e-29        # [g/cm^3] rho_crit,0 / h^2
M_H = 1.67262192e-24
K_B = 1.380649e-16


# ---------------------------------------------------------------------------

def load_cache(path):
    d = np.load(path, allow_pickle=True)
    return dict(tau=np.asarray(d["tau"], dtype=np.float64),
                dv=float(d["dv"]), z=float(d["z"]),
                box_kms=float(d["box_kms"]), h=float(d["h"]),
                hz=float(d["hz"]), source=str(d["source"]))


# ---------------------------------------------------------------------------
# deposicion de Delta, T y v_par sobre la grilla
# ---------------------------------------------------------------------------

def deposit_fields(path, los_name, npix, meta=None, w_floor=0.05,
                   normalize=True):
    """(Delta, T, v_par) por pixel para una LOS.

    Misma geometria que swift_extract.extract_tau: kernel evaluado en la
    distancia 3D r = sqrt(b_perp^2 + dx^2), y misma correccion tipo Shepard.
    Se replica aqui en vez de importarse porque extract_tau solo devuelve tau
    y no vale la pena tocar ese archivo, que ya esta validado.

    Delta se pondera por masa (es una densidad, se deposita y se divide por la
    media). T y v se ponderan por masa tambien, no por HI: el objetivo del
    panel es mostrar el estado del gas, no lo que ve la linea. Si se quisiera
    lo segundo habria que pesar por m*X_H*f_HI, como hace extract_tau.
    """
    if meta is None:
        meta = open_los_file(path)

    with h5py.File(path, "r") as f:
        g = f[los_name]
        coords_int = g["Coordinates"][:].astype(np.float64)
        T = g["Temperatures"][:].astype(np.float64) * _cgs_factor(g["Temperatures"])
        vel = g["Velocities"][:].astype(np.float64) * _cgs_factor(g["Velocities"]) / 1.0e5
        hsml_int = g["SmoothingLengths"][:].astype(np.float64)
        mass = g["Masses"][:].astype(np.float64) * _cgs_factor(g["Masses"])
        rho_c = _cgs_factor(g["Densities"], physical=False)
        rho = g["Densities"][:].astype(np.float64) * rho_c * (1.0 + meta.z) ** 3
        c_len = _cgs_factor(g["Coordinates"], physical=False)

        axis = detect_los_axis(coords_int, meta.boxsize_int)
        tr = [i for i in range(3) if i != axis]
        ray_int, _ = _ray_position(g, coords_int[:, tr],
                                   hsml_int * meta.kernel_gamma,
                                   meta.boxsize_int)

    to_cm = c_len * meta.a
    x_par = coords_int[:, axis] * to_cm
    box_cm = meta.boxsize_int * to_cm
    Hsup = hsml_int * meta.kernel_gamma * to_cm

    dperp = coords_int[:, tr] - ray_int[None, :]
    dperp -= meta.boxsize_int * np.round(dperp / meta.boxsize_int)
    b_perp = np.hypot(dperp[:, 0], dperp[:, 1]) * to_cm

    dR = box_cm / npix
    x_grid = (np.arange(npix) + 0.5) * dR
    dx = x_grid[None, :] - x_par[:, None]
    dx -= box_cm * np.round(dx / box_cm)

    w = kernel_m4(np.sqrt(b_perp[:, None] ** 2 + dx ** 2), Hsup[:, None])

    rho_g = w.T @ mass                       # [g/cm^3] fisico
    wsum = w.T @ (mass / rho)
    if normalize:
        rho_g = rho_g / np.maximum(wsum, w_floor)

    wm = w.T @ mass
    norm = np.where(wm > 0, wm, 1.0)
    T_g = np.where(wm > 0, (w.T @ (mass * T)) / norm, np.nan)
    v_g = np.where(wm > 0, (w.T @ (mass * vel[:, axis])) / norm, np.nan)

    rho_bar = (RHO_CRIT_H2 * meta.h ** 2 * meta.omega_b
               * (1.0 + meta.z) ** 3)        # fisico, mismo (1+z)^3 que rho
    return rho_g / rho_bar, T_g, v_g


# ---------------------------------------------------------------------------
# pick
# ---------------------------------------------------------------------------

def cmd_pick(args):
    A, B = load_cache(args.a), load_cache(args.b)
    if A["tau"].shape != B["tau"].shape:
        sys.exit(f"Caches de distinto tamaño: {A['tau'].shape} vs "
                 f"{B['tau'].shape}. Las LOS tienen que ser las mismas.")
    aA = rescale_tau(A["tau"], args.tau_eff)
    aB = rescale_tau(B["tau"], args.tau_eff)
    tA, tB = aA * A["tau"], aB * B["tau"]
    th = args.tau_thresh

    if args.criterion == "lost-sat":
        # pixeles saturados en CDM que la FCT no tiene: el efecto de Bolton
        score = ((tA > th) & (tB < th)).sum(axis=1)
        desc = f"pixeles con tau_{args.labels[0]} > {th} y tau_{args.labels[1]} < {th}"
    elif args.criterion == "gained-sat":
        score = ((tB > th) & (tA < th)).sum(axis=1)
        desc = f"pixeles con tau_{args.labels[1]} > {th} y tau_{args.labels[0]} < {th}"
    elif args.criterion == "flux-diff":
        score = np.abs(np.exp(-tA) - np.exp(-tB)).sum(axis=1)
        desc = "suma de |F_a - F_b|"
    else:
        score = np.maximum(tA.max(axis=1), tB.max(axis=1))
        desc = "tau maximo"

    order = np.argsort(score)[::-1][:args.top]
    print(f"A({args.labels[0]}) = {aA:.5f}   A({args.labels[1]}) = {aB:.5f}")
    print(f"criterio: {desc}\n")
    print(f"  {'idx':>7}{'score':>12}{'tauMax_a':>12}{'tauMax_b':>12}"
          f"{'<F>_a':>9}{'<F>_b':>9}")
    for i in order:
        print(f"  {i:7d}{score[i]:12.3f}{tA[i].max():12.3f}{tB[i].max():12.3f}"
              f"{np.exp(-tA[i]).mean():9.4f}{np.exp(-tB[i]).mean():9.4f}")
    print("\nElegi uno con score alto pero tauMax_a no absurdo (un DLA gigante")
    print("satura el panel de Delta y no se ve nada mas). Despues:")
    print(f"  python los_panels.py plot --index {order[0]} ...")


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------

def cmd_plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    A, B = load_cache(args.a), load_cache(args.b)
    npix = A["tau"].shape[1]
    i = args.index
    la, lb = args.labels
    ca, cb = args.colors

    aA = rescale_tau(A["tau"], args.tau_eff)
    aB = rescale_tau(B["tau"], args.tau_eff)
    FA = np.exp(-aA * A["tau"][i])
    FB = np.exp(-aB * B["tau"][i])
    
    lost = (aA * A["tau"][i] > args.tau_thresh) & \
       (aB * B["tau"][i] < args.tau_lo)

    metaA, metaB = open_los_file(args.los_a), open_los_file(args.los_b)
    if len(metaA.los_names) != A["tau"].shape[0]:
        print(f"AVISO: el archivo tiene {len(metaA.los_names)} LOS y el cache "
              f"{A['tau'].shape[0]}. Si el cache se hizo con --max-los el "
              f"indice sigue siendo valido; si no, revisá.")
    nmA, nmB = metaA.los_names[i], metaB.los_names[i]
    print(f"depositando {nmA} y {nmB} ...")
    DA, TA, VA = deposit_fields(args.los_a, nmA, npix, meta=metaA)
    DB, TB, VB = deposit_fields(args.los_b, nmB, npix, meta=metaB)

    z = A["z"]
    v = (np.arange(npix) + 0.5) * A["dv"]
    to_hmpc = metaA.h / (metaA.hz / (1.0 + z))     # km/s -> cMpc/h

    npan = 4 if args.with_vel else 3
    fig, axes = plt.subplots(npan, 1, figsize=(9.0, 2.1 * npan + 1.2),
                             sharex=True,
                             gridspec_kw={"hspace": 0.08})
    ax_d, ax_t = axes[0], axes[1]
    ax_f = axes[-1]

    # sombrear donde CDM satura y FCT no: es el punto de la figura
    if args.mark_lost:
        lost = (aA * A["tau"][i] > args.tau_thresh) & \
               (aB * B["tau"][i] < args.tau_thresh)
        if lost.any():
            for ax in axes:
                ax.fill_between(v, 0, 1, where=lost, transform=ax.get_xaxis_transform(),
                                color="gold", alpha=0.25, lw=0, zorder=0)
            print(f"{lost.sum()} pixeles saturados en {la} y no en {lb}")

    ax_d.plot(v, DA, color=ca, lw=1.4, label=la)
    ax_d.plot(v, DB, color=cb, lw=1.2, label=lb)
    ax_d.set(yscale="log", ylabel=r"$\Delta$")
    ax_d.legend(frameon=False, ncol=2, fontsize=9, loc="upper right")

    ax_t.plot(v, TA, color=ca, lw=1.4)
    ax_t.plot(v, TB, color=cb, lw=1.2)
    ax_t.set(yscale="log", ylabel=r"$T$  [K]")

    if args.with_vel:
        axes[2].plot(v, VA, color=ca, lw=1.4)
        axes[2].plot(v, VB, color=cb, lw=1.2)
        axes[2].axhline(0, color="k", lw=0.6, ls=":")
        axes[2].set(ylabel=r"$v_{\parallel}$  [km s$^{-1}$]")

    ax_f.plot(v, FA, color=ca, lw=1.4)
    ax_f.plot(v, FB, color=cb, lw=1.2)
    ax_f.axhline(1.0, color="k", lw=0.6, ls=":")
    ax_f.set(ylim=(-0.05, 1.08), ylabel=r"$F$",
             xlabel=r"$v_{\rm H}$  [km s$^{-1}$]", xlim=(v[0], v[-1]))

    for ax in axes:
        ax.grid(alpha=0.15, which="both")
    if args.xlim:
        for ax in axes:
            ax.set_xlim(*args.xlim)

    axes[0].set_title(
        rf"LOS {i},  $z$ = {z:.1f},  both at $\tau_{{\rm eff}}$ = "
        rf"{args.tau_eff:.4f}   ($A_{{\rm {la}}}$ = {aA:.3f}, "
        rf"$A_{{\rm {lb}}}$ = {aB:.3f})", fontsize=10)
    tx = axes[0].secondary_xaxis("top", functions=(lambda x: x * to_hmpc,
                                                   lambda x: x / to_hmpc))
    tx.set_xlabel(r"comoving  [$h^{-1}$ Mpc]", fontsize=9)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"escrito -> {args.out}")
    print(f"\nDelta max: {la} {DA.max():.1f}   {lb} {DB.max():.1f}")
    print(f"T max:     {la} {np.nanmax(TA):.3e}   {lb} {np.nanmax(TB):.3e}")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    k = sub.add_parser("pick")
    k.add_argument("--a", required=True)
    k.add_argument("--b", required=True)
    k.add_argument("--tau-eff", type=float, required=True)
    k.add_argument("--tau-thresh", type=float, default=10.0)
    k.add_argument("--criterion",
                   choices=["lost-sat", "gained-sat", "flux-diff", "tau-max"],
                   default="lost-sat")
    k.add_argument("--top", type=int, default=10)
    k.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    k.set_defaults(func=cmd_pick)

    q = sub.add_parser("plot")
    q.add_argument("--a", required=True)
    q.add_argument("--b", required=True)
    q.add_argument("--los-a", required=True)
    q.add_argument("--los-b", required=True)
    q.add_argument("--tau-eff", type=float, required=True)
    q.add_argument("--index", type=int, required=True)
    q.add_argument("--tau-thresh", type=float, default=10.0)
    q.add_argument("--mark-lost", action="store_true", default=True)
    q.add_argument("--no-mark-lost", dest="mark_lost", action="store_false")
    q.add_argument("--with-vel", action="store_true")
    q.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    q.add_argument("--colors", nargs=2, default=["C0", "C3"])
    q.add_argument("--xlim", nargs=2, type=float, default=None)
    q.add_argument("--tau-lo", type=float, default=2.0)
    q.add_argument("--out", required=True)
    q.set_defaults(func=cmd_plot)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
