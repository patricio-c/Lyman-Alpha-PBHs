"""
Where the simulations live, and how to find the right snapshot inside them.

The two production boxes have snapshot indices that do NOT line up
(los_0010 is z=3 in one, <base>_0003 is z=3 in the other, because they come
from different output_list files).  Rather than maintain a table of indices
that goes stale the moment a new run appears, we glob the snapshots and read
Header/Redshift.  That is what makes `--run NB_1 --z 3.0` work for runs that
did not exist when this file was written.
"""

from __future__ import annotations

import glob
import os
import re

import h5py
import numpy as np

BASE = os.environ.get("LYA_BASE", "/data/contrib/pad_140/pcolazo/lyman")

# name -> path relative to BASE.  Add new ones here or pass --run-dir.
REGISTRY = {
    "cdm40":  "cdm-box-40-1024",
    "fct40":  "2-fct-box-40-1024",
    "cdm80":  "cdm-box-80-1024",
    "fct80":  "fct-box-80-1024",
}
# the more_power batch: NB_1..NB_6 and poisson_1..poisson_6
for _i in range(1, 7):
    REGISTRY[f"NB_{_i}"] = f"more_power/NB_{_i}"
    REGISTRY[f"poisson_{_i}"] = f"more_power/poisson_{_i}"


def run_dir(name):
    if os.path.isdir(name):
        return os.path.abspath(name)
    if name in REGISTRY:
        return os.path.join(BASE, REGISTRY[name])
    raise SystemExit(f"unknown run '{name}'. Known: {sorted(REGISTRY)}, "
                     f"or pass a directory path.")


def _virtual_first(paths):
    """
    With SWIFT's distributed:1 output a snapshot is a virtual .hdf5 plus
    pieces named .0.hdf5, .1.hdf5, ...  Reading both the virtual file and its
    pieces counts every particle twice.  Keep the virtual file and drop the
    pieces belonging to it.
    """
    virtual, pieces = [], []
    for p in paths:
        (pieces if re.search(r"\.\d+\.hdf5$", p) else virtual).append(p)
    if not virtual:
        return sorted(pieces)
    stems = {p[:-5] for p in virtual}
    keep = list(virtual)
    for p in pieces:
        stem = re.sub(r"\.\d+\.hdf5$", "", p)
        if stem not in stems:
            keep.append(p)
    return sorted(set(keep))


def list_hdf5(directory, pattern="**/*.hdf5"):
    """Every hdf5 under the directory, with distributed pieces collapsed."""
    return _virtual_first(glob.glob(os.path.join(directory, pattern),
                                    recursive=True))


def file_redshift(path):
    """Redshift of a snapshot OR a SWIFT LOS file (they differ in layout)."""
    try:
        with h5py.File(path, "r") as f:
            for grp in ("Cosmology", "Header"):
                if grp in f:
                    at = f[grp].attrs
                    for key in ("Redshift", "redshift"):
                        if key in at:
                            return float(np.ravel(at[key])[0])
                    for key in ("Scale-factor", "Time"):
                        if key in at:
                            a = float(np.ravel(at[key])[0])
                            if 0 < a <= 1.0:
                                return 1.0 / a - 1.0
    except Exception:
        return None
    return None


def is_los_file(path):
    """True if the file holds LOS_XXXX groups rather than PartTypeN."""
    try:
        with h5py.File(path, "r") as f:
            return any(k.startswith("LOS_") for k in f)
    except Exception:
        return False


def resolve_los_file(run, z=None, path=None, tol=0.05, verbose=True):
    """
    Find the SWIFT LOS file of a run at a given redshift.

    Snapshot indices are not comparable between runs, so nothing here uses
    one: every hdf5 in the directory is opened, its redshift read, the LOS
    files kept, and the closest to `z` returned.
    """
    if path and os.path.isfile(path):
        return os.path.abspath(path)
    d = run_dir(run)
    cands = [p for p in list_hdf5(d) if is_los_file(p)]
    if not cands:
        raise SystemExit(
            f"no SWIFT LOS files (groups LOS_XXXX) found under {d}.\n"
            f"hdf5 files present: {[os.path.basename(x) for x in list_hdf5(d)][:20]}")
    zs = [(p, file_redshift(p)) for p in cands]
    zs = [(p, zz) for p, zz in zs if zz is not None]
    if z is None:
        if verbose:
            for p, zz in sorted(zs, key=lambda t: t[1]):
                print(f"[runs] {os.path.basename(p):40s} z = {zz:.4f}")
        raise SystemExit("give --z (the list above shows what is available)")
    zs.sort(key=lambda t: abs(t[1] - z))
    best, zbest = zs[0]
    if verbose:
        print(f"[runs] {run}: LOS files at z = "
              f"{sorted({round(t[1], 3) for t in zs})}")
        print(f"[runs] picked {os.path.basename(best)}  z = {zbest:.4f}")
    if abs(zbest - z) > tol:
        raise SystemExit(f"closest LOS file is z={zbest:.4f}, you asked for "
                         f"z={z:.4f}. Pass --los-file explicitly.")
    return best


def list_snapshots(directory, pattern="**/*.hdf5"):
    paths = glob.glob(os.path.join(directory, pattern), recursive=True)
    paths = [p for p in paths if "snap" in os.path.basename(p).lower()
             or "output" in os.path.basename(p).lower()]
    return _virtual_first(paths)


def snapshot_redshift(path):
    try:
        with h5py.File(path, "r") as f:
            h = f["Header"].attrs
            for key in ("Redshift", "redshift"):
                if key in h:
                    return float(np.ravel(h[key])[0])
            if "Scale-factor" in h:
                a = float(np.ravel(h["Scale-factor"])[0])
                return 1.0 / a - 1.0
    except Exception:
        return None
    return None


def resolve_snapshot(run, z=None, snap=None, tol=0.05, verbose=True):
    """
    Return the path of the snapshot to use.

    --snap wins if given (an index or a full path).  Otherwise every
    snapshot in the run directory is opened, its redshift read, and the
    closest one to `z` returned.  Fails loudly if the closest is farther
    than `tol`, because silently analysing z=2 while thinking it is z=3 is
    the kind of thing that costs a week.
    """
    d = run_dir(run)
    if snap is not None and os.path.isfile(str(snap)):
        return os.path.abspath(str(snap))

    cands = list_snapshots(d)
    if not cands:
        raise SystemExit(f"no snapshots found under {d}")

    if snap is not None:
        tag = f"{int(snap):04d}"
        hit = [p for p in cands if tag in os.path.basename(p)]
        if not hit:
            raise SystemExit(f"no snapshot matching index {tag} in {d}")
        return hit[0]

    if z is None:
        raise SystemExit("give either --snap or --z")

    zs = [(p, snapshot_redshift(p)) for p in cands]
    zs = [(p, zz) for p, zz in zs if zz is not None]
    if not zs:
        raise SystemExit(f"could not read Header/Redshift from any file in {d}")
    zs.sort(key=lambda t: abs(t[1] - z))
    best, zbest = zs[0]
    if verbose:
        print(f"[runs] {run}: {len(zs)} snapshots, "
              f"z available = {sorted({round(t[1], 3) for t in zs})}")
        print(f"[runs] picked {os.path.basename(best)}  z = {zbest:.4f}")
    if abs(zbest - z) > tol:
        raise SystemExit(f"closest snapshot in {d} is z={zbest:.4f}, "
                         f"you asked for z={z:.4f} (tol {tol}). "
                         f"Pass --snap explicitly if this is what you want.")
    return best
