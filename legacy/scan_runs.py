#!/usr/bin/env python
"""
scan_runs.py - Inventario de corridas: casa cada archivo de LOS con su
snapshot por redshift, y detecta el bug de range_when_shooting_down.

Uso:
    python scan_runs.py ~/lyman
    python scan_runs.py ~/lyman --only cdm-box-40-1024
"""

import argparse
import glob
import os

import h5py
import numpy as np


def z_of(path):
    """Redshift de un archivo HDF5 de SWIFT (snapshot o LOS)."""
    try:
        with h5py.File(path, "r") as f:
            for grp, key in [("Cosmology", "Redshift"), ("Header", "Redshift")]:
                if grp in f and key in f[grp].attrs:
                    return float(np.atleast_1d(f[grp].attrs[key])[0])
    except Exception:
        pass
    return None


def snapshot_files(entry):
    """Devuelve la lista de archivos de un snapshot (distribuido o unico)."""
    if os.path.isdir(entry):
        return sorted(glob.glob(os.path.join(entry, "*.hdf5")))
    return [entry] if entry.endswith(".hdf5") else []


def los_info(path):
    """(z, box_interno, n_los, x_max_axial, truncado?)"""
    with h5py.File(path, "r") as f:
        z = float(np.atleast_1d(f["Cosmology"].attrs["Redshift"])[0])
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        names = sorted(k for k in f.keys() if k.startswith("LOS_"))
        xmax = 0.0
        for nm in names[:20]:
            c = f[nm]["Coordinates"][:]
            ang = 2.0 * np.pi * c / box
            R = np.hypot(np.cos(ang).mean(axis=0), np.sin(ang).mean(axis=0))
            xmax = max(xmax, float(c[:, int(np.argmin(R))].max()))
    return z, box, len(names), xmax, xmax < 0.95 * box


def scan_run(run_dir):
    name = os.path.basename(run_dir.rstrip("/"))
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")

    ymls = [p for p in glob.glob(os.path.join(run_dir, "*.yml"))
            if "used_parameters" not in p and "unused" not in p]
    box = None

    los_files = sorted(glob.glob(os.path.join(run_dir, "los_*.hdf5")))
    if los_files:
        print(f"\n  Archivos de LOS ({len(los_files)}):")
        print(f"    {'archivo':<18}{'z':>8}{'nLOS':>8}{'x_max':>10}{'caja':>10}  estado")
        for p in los_files:
            try:
                z, box, n, xmax, trunc = los_info(p)
                flag = (f"TRUNCADO ({100 * xmax / box:.0f}% de la caja)"
                        if trunc else "completo")
                print(f"    {os.path.basename(p):<18}{z:>8.3f}{n:>8}"
                      f"{xmax:>10.2f}{box:>10.2f}  {flag}")
            except Exception as ex:
                print(f"    {os.path.basename(p):<18}  error: {ex}")

    # snapshots: directorios o archivos sueltos
    cands = [p for p in glob.glob(os.path.join(run_dir, "*"))
             if (os.path.isdir(p) and glob.glob(os.path.join(p, "*.hdf5")))
             or (p.endswith(".hdf5") and not os.path.basename(p).startswith("los_"))]
    snaps = []
    for c in sorted(cands):
        fl = snapshot_files(c)
        if not fl:
            continue
        z = z_of(fl[0])
        if z is not None:
            snaps.append((c, z, len(fl)))

    if snaps:
        print(f"\n  Snapshots ({len(snaps)}):")
        print(f"    {'entrada':<32}{'z':>8}{'archivos':>10}")
        for c, z, nf in snaps:
            print(f"    {os.path.basename(c):<32}{z:>8.3f}{nf:>10}")
    else:
        print("\n  Snapshots: NINGUNO encontrado. Sin snapshot no se puede "
              "regenerar; habria que relanzar la corrida.")

    if ymls:
        print(f"\n  Bloque LineOfSight de {os.path.basename(ymls[0])}:")
        inside, shown = False, 0
        for line in open(ymls[0]):
            if line.strip().startswith("LineOfSight:"):
                inside = True
            elif inside and line.strip() and not line.startswith((" ", "\t")):
                inside = False
            if inside and shown < 14:
                s = line.rstrip()
                bad = ""
                if "range" in s and box:
                    nums = [float(x) for x in
                            s.replace("[", " ").replace("]", " ")
                             .replace(",", " ").split() if _isnum(x)]
                    if nums and abs(nums[-1] - box) > 0.01 * box:
                        bad = f"   <-- deberia ser {box:.4f}"
                print(f"    {s}{bad}")
                shown += 1

    # emparejamiento
    if los_files and snaps:
        print(f"\n  Emparejamiento LOS <-> snapshot:")
        for p in los_files:
            try:
                z, box, n, xmax, trunc = los_info(p)
            except Exception:
                continue
            best = min(snaps, key=lambda s: abs(s[1] - z))
            ok = abs(best[1] - z) < 0.02
            print(f"    {os.path.basename(p):<18} z={z:6.3f}  ->  "
                  f"{os.path.basename(best[0]) if ok else 'SIN SNAPSHOT A ESE z'}"
                  f"{'' if ok else f' (mas cercano z={best[1]:.3f})'}")


def _isnum(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="directorio con las corridas, ej ~/lyman")
    ap.add_argument("--only", default=None, help="una sola corrida")
    args = ap.parse_args()

    root = os.path.expanduser(args.root)
    runs = ([os.path.join(root, args.only)] if args.only else
            sorted(d for d in glob.glob(os.path.join(root, "*"))
                   if os.path.isdir(d) and glob.glob(os.path.join(d, "los_*.hdf5"))))
    if not runs:
        raise SystemExit(f"No encontre corridas con los_*.hdf5 en {root}")
    for r in runs:
        scan_run(r)
    print(f"\n{'=' * 78}")
    print("Para cada corrida usa el snapshot emparejado al z que te interese.")


if __name__ == "__main__":
    main()
