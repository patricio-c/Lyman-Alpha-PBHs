#!/usr/bin/env python
"""
plot_3d_1d.py - Figura combinada: cociente FCT/CDM en 3D (materia) y en 1D
(flujo del forest), con el eje k compartido en Mpc^-1 comovil y un eje
secundario en s/km.

Por que juntarlos: el P1D es una integral sobre modos transversales,

    P_1D(k_par) = int dk_perp k_perp P_3D(k) / (2 pi)

asi que el P1D en un dado k recibe contribucion de TODOS los modos 3D con
k >= k_par. La señal aparece entonces a k menores en 1D que en 3D. Alinear
los dos paneles hace visible ese corrimiento, que es fisica y no un desajuste.

Entradas:
  --ratio3d   texto con k [Mpc^-1], P_CDM, P_FCT  (o k y cociente con --ratio3d-is-ratio)
  --ratio1d   la tabla que escribe ratio_fct_cdm.py: k[s/km], P_CDM, P_FCT, cociente, error

Uso:
    python plot_3d_1d.py --ratio3d pk3d.txt --ratio1d ratio_z3.txt \\
        --z 3.0 --hz 306.31 --kt 10 --box-mpc 58.7372 --ngas 512 \\
        --out fct_cdm_3d_1d.png
"""

import argparse

import numpy as np

LAMBDA_LYA = 1215.6701
C_KMS = 2.99792458e5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio3d", required=True,
                    help="P(k) 3D de la primera corrida (CDM). Si trae los "
                         "dos espectros en columnas, usar solo este.")
    ap.add_argument("--ratio3d-fct", default=None,
                    help="P(k) 3D de la segunda corrida (FCT), archivo aparte")
    ap.add_argument("--cols3d", nargs=2, type=int, default=None,
                    help="columnas (k, P) del archivo 3D; por defecto se "
                         "autodetectan")
    ap.add_argument("--ratio3d-is-ratio", action="store_true",
                    help="el archivo 3D ya trae k y cociente")
    ap.add_argument("--ratio1d", required=True)
    ap.add_argument("--z", type=float, required=True)
    ap.add_argument("--hz", type=float, required=True,
                    help="H(z) [km/s/Mpc]; a z=3 con tu cosmologia, 306.31")
    ap.add_argument("--kt", type=float, default=None,
                    help="escala de ruptura del modelo, en Mpc^-1 comovil")
    ap.add_argument("--box-mpc", type=float, default=None,
                    help="caja comovil [Mpc], para marcar el modo fundamental")
    ap.add_argument("--ngas", type=int, default=None,
                    help="lado efectivo de la grilla de gas, para marcar el "
                         "Nyquist de particula")
    ap.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    ap.add_argument("--single", action="store_true",
                    help="un solo panel con las dos curvas superpuestas")
    ap.add_argument("--xscale", choices=["log", "linear"], default="log")
    ap.add_argument("--yscale", choices=["log", "linear"], default="linear")
    ap.add_argument("--xlim", nargs=2, type=float, default=None,
                    help="rango en Mpc^-1; con --xscale linear conviene "
                         "acotarlo, por ejemplo 0 30")
    ap.add_argument("--ylim", nargs=2, type=float, default=None)
    ap.add_argument("--out", default="fct_cdm_3d_1d.png")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # k[Mpc^-1 comovil] = k[s/km] * H(z)/(1+z)
    conv = args.hz / (1.0 + args.z)
    la, lb = args.labels

    # --- 3D ----------------------------------------------------------------
    def load_pk(path, cols):
        head = []
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    head.append(line.rstrip())
                elif line.strip():
                    break
        d = np.loadtxt(path)
        if cols is not None:
            return d[:, cols[0]], d[:, cols[1]], head
        # Autodeteccion. SWIFT escribe (0) Redshift, (1) k, (2) potencia,
        # asi que NO se puede asumir que la columna 0 sea k: la de redshift es
        # constante. k es la primera columna positiva y estrictamente
        # creciente; P es la siguiente positiva que decrece en promedio.
        kcol = None
        for j in range(d.shape[1]):
            v = d[:, j]
            if np.all(v > 0) and np.all(np.diff(v) > 0):
                kcol = j
                break
        if kcol is None:
            raise SystemExit(
                f"{path}: no encontre una columna de k (positiva y creciente). "
                "Pasala a mano con --cols3d K P. Cabecera:\n  "
                + "\n  ".join(head[:15]))
        pcol = None
        for j in range(kcol + 1, d.shape[1]):
            v = d[:, j]
            if np.all(v > 0) and v[0] > v[-1]:
                pcol = j
                break
        if pcol is None:
            raise SystemExit(f"{path}: no encontre la columna de P(k). "
                             "Usa --cols3d.")
        print(f"  autodeteccion: k = columna {kcol}, P = columna {pcol}")
        return d[:, kcol], d[:, pcol], head

    if args.ratio3d_fct:
        k3, p_a, h_a = load_pk(args.ratio3d, args.cols3d)
        k3b, p_b, _ = load_pk(args.ratio3d_fct, args.cols3d)
        print("Cabecera del archivo 3D:")
        for l in h_a[:12]:
            print("  " + l)
        print(f"  -> usando {len(k3)} puntos, k de {k3.min():.4g} a "
              f"{k3.max():.4g}")
        if not np.allclose(k3, k3b, rtol=1e-6):
            p_b = np.interp(k3, k3b, p_b)
            print("  (grillas de k distintas: se interpolo la segunda)")
        r3 = p_b / p_a
    else:
        d3 = np.loadtxt(args.ratio3d)
        if args.ratio3d_is_ratio:
            k3, r3 = d3[:, 0], d3[:, 1]
        else:
            k3, r3 = d3[:, 0], d3[:, 2] / d3[:, 1]

    # --- 1D ----------------------------------------------------------------
    d1 = np.loadtxt(args.ratio1d)
    kv, r1 = d1[:, 0], d1[:, 3]
    e1 = d1[:, 4] if d1.shape[1] > 4 else np.zeros_like(r1)
    k1 = kv * conv

    # --- escalas de referencia --------------------------------------------
    marks = []
    if args.kt:
        marks.append((args.kt, f"$k_t$ = {args.kt:g} Mpc$^{{-1}}$", "C3"))
    if args.box_mpc:
        marks.append((2*np.pi/args.box_mpc, "modo fundamental", "grey"))
    if args.ngas:
        marks.append((np.pi*args.ngas/args.box_mpc,
                      "Nyquist de partícula", "grey"))
    R_z = C_KMS * 0.8 / ((1.0 + args.z) * LAMBDA_LYA)
    desi = (1.0e-3 * conv, 0.5 * np.pi / R_z * conv)

    lo = max(k3.min(), k1.min()) * 0.7
    hi = min(k3.max(), k1.max() * 1.6)

    if args.single:
        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.axhline(1.0, color="k", ls="--", lw=1)
        ax.axvspan(*desi, color="green", alpha=0.07, zorder=0)
        for kk, txt, col in marks:
            if lo < kk < hi:
                ax.axvline(kk, color=col, ls=":", lw=1.2)
        ax.plot(k3, r3, color="k", lw=1.8,
                label=f"materia 3D  ({lb}/{la})")
        rel1 = e1 / np.maximum(r1, 1e-30)
        ax.fill_between(k1, r1 * (1 - rel1), r1 * (1 + rel1),
                        color="C3", alpha=0.25)
        ax.plot(k1, r1, color="C3", lw=1.9,
                label=f"flujo del bosque 1D  ({lb}/{la})")
        ax.set(xscale=args.xscale, yscale=args.yscale,
               xlabel=r"$k$ [Mpc$^{-1}$ comóvil]",
               ylabel=f"$P$({lb}) / $P$({la})",
               title=f"z = {args.z:.1f}")
        ax.set_xlim(*(args.xlim if args.xlim else (lo, hi)))
        if args.ylim:
            ax.set_ylim(*args.ylim)
        ax.grid(alpha=0.2, which="both")
        ax.legend(frameon=False, loc="upper left")

        y0, y1 = ax.get_ylim()
        for kk, txt, col in marks:
            if ax.get_xlim()[0] < kk < ax.get_xlim()[1]:
                ax.text(kk, y1 - 0.02 * (y1 - y0), " " + txt, color=col,
                        rotation=90, va="top", ha="left", fontsize=8)
        d0, d1 = desi
        ax.text(np.sqrt(d0 * d1) if args.xscale == "log" else (d0 + d1) / 2,
                y0 + 0.03 * (y1 - y0), "ventana DESI", color="green",
                ha="center", va="bottom", fontsize=9)

        # segundo eje x ABAJO, desplazado, en s/km
        sx = ax.secondary_xaxis(-0.18, functions=(lambda x: x / conv,
                                                  lambda x: x * conv))
        sx.set_xlabel(r"$k$ [s/km]")

        fig.savefig(args.out, dpi=150, bbox_inches="tight")
        print(f"escrito -> {args.out}")
        print(f"\nconversion usada: k[Mpc^-1] = k[s/km] x {conv:.2f}")
        for kk, txt, _ in marks:
            print(f"  {txt}: {kk:.3f} Mpc^-1 = {kk/conv:.4f} s/km")
        print(f"\n  {'k [Mpc^-1]':>12}{'k [s/km]':>11}{'3D':>9}{'1D':>9}")
        for kt_ in [0.5, 1, 2, 5, 10, 20, 40]:
            if k1.min() <= kt_ <= k1.max() and k3.min() <= kt_ <= k3.max():
                print(f"  {kt_:12.2f}{kt_/conv:11.4f}"
                      f"{np.interp(kt_, k3, r3):9.3f}"
                      f"{np.interp(kt_, k1, r1):9.3f}")
        return

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 7.4), sharex=True,
                                 gridspec_kw={"hspace": 0.08})

    for a_ in (ax, bx):
        a_.axhline(1.0, color="k", ls="--", lw=1)
        a_.axvspan(*desi, color="green", alpha=0.07, zorder=0)
        for kk, _, col in marks:
            if lo < kk < hi:
                a_.axvline(kk, color=col, ls=":", lw=1.2)
        a_.set(xscale="log", xlim=(lo, hi))
        a_.grid(alpha=0.2, which="both")

    ax.plot(k3, r3, color="k", lw=1.6)
    ax.set(yscale=args.yscale,
           ylabel=f"$P_{{\\rm 3D}}$({lb}) / $P_{{\\rm 3D}}$({la})")
    ax.set_title(f"z = {args.z:.1f}   ·   materia (3D) arriba, "
                 f"flujo del bosque (1D) abajo", fontsize=11)

    bx.fill_between(k1, r1*(1-e1/np.maximum(r1, 1e-30)),
                    r1*(1+e1/np.maximum(r1, 1e-30)), color="C3", alpha=0.25)
    bx.plot(k1, r1, color="C3", lw=1.7)
    bx.set(ylabel=f"$P_{{\\rm 1D}}$({lb}) / $P_{{\\rm 1D}}$({la})",
           xlabel=r"$k$ [Mpc$^{-1}$ comóvil]")

    # etiquetas de las escalas marcadas
    y0, y1 = ax.get_ylim()
    for kk, txt, col in marks:
        if lo < kk < hi:
            ax.text(kk, y1*0.6, " " + txt, color=col, rotation=90,
                    va="top", ha="left", fontsize=8)
    yb0, yb1 = bx.get_ylim()
    bx.text(np.sqrt(desi[0]*desi[1]), yb1 - 0.03*(yb1-yb0), "ventana DESI",
            color="green", ha="center", va="top", fontsize=9)

    # eje secundario en s/km
    tx = ax.secondary_xaxis("top", functions=(lambda x: x/conv,
                                              lambda x: x*conv))
    tx.set_xlabel(r"$k$ [s/km]")

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"escrito -> {args.out}")
    print(f"\nconversion usada: k[Mpc^-1] = k[s/km] x {conv:.2f}")
    print(f"ventana DESI: {desi[0]:.4f} - {desi[1]:.3f} Mpc^-1")
    for kk, txt, _ in marks:
        print(f"  {txt}: {kk:.3f} Mpc^-1 = {kk/conv:.4f} s/km")

    # comparacion cuantitativa en escalas comunes
    print(f"\n  {'k [Mpc^-1]':>12}{'k [s/km]':>11}{'3D':>9}{'1D':>9}")
    for kt_ in [0.5, 1, 2, 5, 10, 20, 40]:
        if k1.min() <= kt_ <= k1.max() and k3.min() <= kt_ <= k3.max():
            print(f"  {kt_:12.2f}{kt_/conv:11.4f}"
                  f"{np.interp(kt_, k3, r3):9.3f}"
                  f"{np.interp(kt_, k1, r1):9.3f}")


if __name__ == "__main__":
    main()

