#!/usr/bin/env python
"""relosz.py - resuelve LOS y snapshot por redshift y llama a relos.py.
Uso: python relosz.py ../lyman/cdm-box-40-1024 3.0 --out o.hdf5 --uniform 4 --seed 12345
"""
import glob, os, subprocess, sys
import h5py, numpy as np

def zof(p):
    try:
        with h5py.File(p, "r") as f:
            for g in ("Cosmology", "Header"):
                if g in f and "Redshift" in f[g].attrs:
                    return float(np.atleast_1d(f[g].attrs["Redshift"])[0])
    except Exception:
        return None

run, ztgt, rest = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
run = os.path.expanduser(run)

los = [(p, zof(p)) for p in sorted(glob.glob(os.path.join(run, "los_*.hdf5")))]
los = [(p, z) for p, z in los if z is not None]

snaps = []
for d in sorted(glob.glob(os.path.join(run, "*"))):
    if os.path.isdir(d):
        pcs = sorted(glob.glob(os.path.join(d, "*.hdf5")))
        if pcs and zof(pcs[0]) is not None:
            snaps.append((d, zof(pcs[0])))

hit = [d for d, z in snaps if abs(z - ztgt) < 0.02]
if not hit:
    print(f"No hay snapshot a z={ztgt} en {run}")
    print("  con snapshot:", ", ".join(f"{z:g}" for _, z in sorted(snaps, key=lambda t: -t[1])))
    sys.exit(1)

lo = [p for p, z in los if abs(z - ztgt) < 0.02]
if not lo:
    lo = [los[0][0]]
    print(f"  (sin LOS a z={ztgt}; uso {os.path.basename(lo[0])} solo como plantilla)")

base = os.path.basename(hit[0])
pat = os.path.join(hit[0], f"{base}.[0-9]*.hdf5")
if not glob.glob(pat):
    pat = os.path.join(hit[0], "*.hdf5")

print(f"z={ztgt}:\n  plantilla: {lo[0]}\n  snapshot : {pat}\n")
sys.exit(subprocess.call([sys.executable,
                          os.path.join(os.path.dirname(os.path.abspath(__file__)), "relos.py"),
                          "--old-los", lo[0], "--snapshot", pat] + rest))
