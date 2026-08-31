#!/usr/bin/env python
"""
fix_los_attrs.py - Agrega a los archivos de LOS regenerados los atributos de
grupo que SWIFT escribe y que spectWizard necesita.

Motivo: relos.py, al generar rayos nuevos con --uniform, escribia solo Xpos e
Ypos. SWIFT ademas escribe NumParts, Xaxis, Yaxis y Zaxis, y spectWizard los
lee (KeyError: 'Xaxis' en reading_simulations.py). Toda esa informacion ya
esta en el archivo, asi que se puede completar in situ: NO hace falta
regenerar nada.

Convencion de SWIFT (line_of_sight.c): Zaxis es el eje de INTEGRACION,
Xaxis e Yaxis son los dos transversales, y Xpos/Ypos son las coordenadas del
rayo a lo largo de Xaxis e Yaxis respectivamente.

Uso:
    # 1. ver que atributos tiene un archivo ORIGINAL de SWIFT, como referencia
    python fix_los_attrs.py --inspect ../lyman/cdm-box-40-1024/los_0010.hdf5

    # 2. ver que le falta a los regenerados, sin tocar nada
    python fix_los_attrs.py regen/cdm40_z3.0_uni512_seed12345.hdf5

    # 3. aplicarlo
    python fix_los_attrs.py regen/*.hdf5 --apply \
        --reference ../lyman/cdm-box-40-1024/los_0010.hdf5
"""

import argparse
import glob
import sys

import h5py
import numpy as np

WANTED = ["NumParts", "Xaxis", "Yaxis", "Zaxis", "Xpos", "Ypos"]


def inspect(path):
    with h5py.File(path, "r") as f:
        names = sorted(k for k in f if k.startswith("LOS_"))
        if not names:
            sys.exit(f"{path}: no hay grupos LOS_XXXX")
        print(f"\n{path}   ({len(names)} lineas de visión)")
        for nm in names[:3]:
            print(f"  {nm}:")
            for k, v in f[nm].attrs.items():
                a = np.atleast_1d(v)
                print(f"    {k:<12} = {a[0]!r:<16} dtype={a.dtype}")
        print(f"  datasets: {sorted(f[names[0]].keys())[:6]} ...")


def reference_attrs(path):
    """Lee de un archivo ORIGINAL los nombres y dtypes exactos que usa SWIFT."""
    with h5py.File(path, "r") as f:
        nm = sorted(k for k in f if k.startswith("LOS_"))[0]
        return {k: np.atleast_1d(v).dtype for k, v in f[nm].attrs.items()}


def detect_axis(coords, box):
    """Eje de integracion por concentracion circular (robusta al borde)."""
    ang = 2.0 * np.pi * np.asarray(coords, dtype=np.float64) / box
    R = np.hypot(np.cos(ang).mean(axis=0), np.sin(ang).mean(axis=0))
    axis = int(np.argmin(R))
    if R[axis] > 0.6 or np.delete(R, axis).min() < 0.9:
        raise ValueError(f"eje ambiguo, R = {R}")
    return axis


def fix(path, apply_it, dtypes):
    with h5py.File(path, "r+" if apply_it else "r") as f:
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        names = sorted(k for k in f if k.startswith("LOS_"))
        missing = [k for k in WANTED if k not in f[names[0]].attrs]
        if not missing:
            print(f"  {path}: ya tiene los {len(WANTED)} atributos")
            return True
        print(f"  {path}: faltan {missing}  ({len(names)} LOS)")
        if not apply_it:
            return True

        n_ax = np.zeros(3, dtype=int)
        for nm in names:
            g = f[nm]
            c = g["Coordinates"][:].astype(np.float64)
            zax = detect_axis(c, box)
            tr = [i for i in range(3) if i != zax]
            n_ax[zax] += 1

            vals = {
                "NumParts": c.shape[0],
                "Xaxis": tr[0],
                "Yaxis": tr[1],
                "Zaxis": zax,
            }
            # Xpos/Ypos: si ya estan se respetan; si no, se estiman con la
            # mediana circular de las particulas de menor radio de suavizado
            if "Xpos" not in g.attrs or "Ypos" not in g.attrs:
                hs = g["SmoothingLengths"][:]
                idx = np.argsort(hs)[: max(10, len(hs) // 10)]
                ang = 2.0 * np.pi * c[idx][:, tr] / box
                ma = np.arctan2(np.sin(ang).mean(axis=0),
                                np.cos(ang).mean(axis=0))
                pos = (ma / (2.0 * np.pi) * box) % box
                vals["Xpos"], vals["Ypos"] = float(pos[0]), float(pos[1])

            for k, v in vals.items():
                if k in g.attrs:
                    continue
                dt = dtypes.get(k)
                if dt is not None:
                    g.attrs[k] = np.array([v], dtype=dt)
                elif k in ("NumParts", "Xaxis", "Yaxis", "Zaxis"):
                    g.attrs[k] = np.array([v], dtype=np.int32)
                else:
                    g.attrs[k] = np.array([v], dtype=np.float64)

        print(f"    aplicado. Lineas por eje: x={n_ax[0]} y={n_ax[1]} "
              f"z={n_ax[2]}")
        return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--inspect", default=None,
                    help="mostrar los atributos de un archivo y salir")
    ap.add_argument("--reference", default=None,
                    help="archivo ORIGINAL de SWIFT del que copiar los dtypes "
                         "exactos de los atributos")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    if args.inspect:
        inspect(args.inspect)
        return

    dtypes = {}
    if args.reference:
        dtypes = reference_attrs(args.reference)
        print(f"dtypes de referencia ({args.reference}):")
        for k, v in sorted(dtypes.items()):
            print(f"  {k:<12} {v}")
    else:
        print("Sin --reference: se usan int32 para los ejes y float64 para las "
              "posiciones.\nSi spectWizard es quisquilloso con los tipos, "
              "pasá un archivo original con --reference.")

    files = []
    for p in args.paths:
        files.extend(sorted(glob.glob(p)) or [p])
    if not files:
        sys.exit("No pasaste archivos")

    print(f"\n{'APLICANDO' if args.apply else 'SIMULACION (usa --apply)'}\n")
    for p in files:
        fix(p, args.apply, dtypes)
    if not args.apply:
        print("\nNada se modifico. Volve a correr con --apply.")
    else:
        print("\nListo. Verifica con --inspect sobre uno de los archivos.")


if __name__ == "__main__":
    main()
