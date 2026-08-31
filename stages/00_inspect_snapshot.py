#!/usr/bin/env python3
"""
Stage 00 - what is actually inside this snapshot?

Run this before anything else on a run you have not touched before, and
always before t7_stars_back.py.  It answers one question:

    is the gas that the star-formation scheme removed still in the file,
    and if so, how do I tell it apart from the dark matter?

In the current 40 Mpc/h pair the answer is "yes, by mass": converted gas
sits in PartType1 carrying m_gas = 9.4445e6 Msun while the dark matter
carries 4.975e-3 internal units, so a mass cut separates them cleanly.  In
the runs coming out of more_power the two species have the SAME particle
mass, so that trick dies and we need either a PartType4 or a flag field.
This script tells you which case you are in instead of you guessing.

Usage
-----
    python stages/00_inspect_snapshot.py --run fct40 --z 3.0
    python stages/00_inspect_snapshot.py --run NB_1 --z 3.0 --deep
    python stages/00_inspect_snapshot.py --snap /path/to/snap_0003.hdf5

Options
-------
    --run NAME      registry name or a directory (see common/runs.py)
    --z Z           pick the snapshot closest to this redshift
    --snap PATH     explicit file or index, overrides --z
    --deep          also histogram the masses (reads the Masses arrays)
    --max-datasets  how many dataset names to print per PartType (default 40)
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import runs  # noqa: E402

PARTTYPE_MEANING = {
    0: "gas",
    1: "dark matter (and, in QLA runs, the converted gas)",
    2: "background / boundary",
    3: "sinks",
    4: "stars",
    5: "black holes",
    6: "neutrinos",
}


def show_attrs(g, indent="    ", limit=60):
    for k in list(g.attrs)[:limit]:
        v = g.attrs[k]
        if isinstance(v, (bytes, np.bytes_)):
            v = v.decode(errors="replace")
        elif isinstance(v, np.ndarray) and v.size > 8:
            v = f"array{v.shape} [{v.ravel()[0]} ... {v.ravel()[-1]}]"
        print(f"{indent}{k:38s} = {v}")


def mass_populations(m, rtol=1e-3, max_groups=12):
    """
    Cluster a mass array into distinct populations.

    Particle masses in these runs are exactly equal within a species, so
    rounding to a relative tolerance and counting is enough and is far
    cheaper than np.unique on 10^8 float64.
    """
    m = np.asarray(m, dtype=np.float64)
    m = m[np.isfinite(m) & (m > 0)]
    if m.size == 0:
        return []
    key = np.round(np.log10(m) / rtol).astype(np.int64)
    vals, counts = np.unique(key, return_counts=True)
    order = np.argsort(-counts)[:max_groups]
    out = []
    for i in order:
        sel = key == vals[i]
        out.append((float(m[sel].mean()), int(counts[i])))
    return sorted(out, key=lambda t: -t[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=None)
    ap.add_argument("--z", type=float, default=None)
    ap.add_argument("--snap", default=None)
    ap.add_argument("--deep", action="store_true",
                    help="read Masses and histogram the populations")
    ap.add_argument("--max-datasets", type=int, default=40)
    args = ap.parse_args()

    if args.run:
        path = runs.resolve_snapshot(args.run, z=args.z, snap=args.snap)
    elif args.snap:
        path = args.snap
    else:
        ap.error("give --run (with --z or --snap) or --snap PATH")

    print("=" * 74)
    print(f"FILE   {path}")
    print(f"SIZE   {os.path.getsize(path) / 1e9:.2f} GB")
    print("=" * 74)

    with h5py.File(path, "r") as f:
        print("\n--- Header " + "-" * 62)
        show_attrs(f["Header"])

        for grp in ("Cosmology", "Units", "HydroScheme", "SubgridScheme",
                    "Policy", "Code", "StarsScheme", "GravityScheme"):
            if grp in f:
                print(f"\n--- {grp} " + "-" * (70 - len(grp)))
                show_attrs(f[grp])

        virtual = "Virtual" in f["Header"].attrs
        if virtual:
            print("\n[!] this is a VIRTUAL file (distributed output). "
                  "Read THIS file, never glob the .N.hdf5 pieces as well.")

        print("\n--- particle types " + "-" * 55)
        present = sorted(int(k.replace("PartType", ""))
                         for k in f if k.startswith("PartType"))
        npart = f["Header"].attrs.get("NumPart_Total", None)

        for pt in present:
            g = f[f"PartType{pt}"]
            first = next(iter(g))
            n = g[first].shape[0]
            print(f"\nPartType{pt}  ({PARTTYPE_MEANING.get(pt, '?')})"
                  f"   N = {n:,}")
            names = list(g)
            print(f"  datasets ({len(names)}): "
                  f"{', '.join(names[:args.max_datasets])}"
                  f"{' ...' if len(names) > args.max_datasets else ''}")

            flagged = [x for x in names
                       if any(t in x.lower() for t in
                              ("birth", "formation", "progenitor", "converted",
                               "star", "sfr", "flag", "type", "lastsf"))]
            if flagged:
                print(f"  [!] possible provenance fields: {', '.join(flagged)}")

            if args.deep and "Masses" in g:
                pops = mass_populations(g["Masses"][:])
                print("  mass populations (internal units):")
                for mval, cnt in pops:
                    print(f"      m = {mval:.6e}   N = {cnt:,}"
                          f"   ({100.0 * cnt / n:.2f}%)")

        print("\n" + "=" * 74)
        print("VERDICT for t7 (putting the converted gas back in the forest)")
        print("=" * 74)

        has_stars = 4 in present and (
            f["PartType4"][next(iter(f["PartType4"]))].shape[0] > 0)
        if has_stars:
            n4 = f["PartType4"][next(iter(f["PartType4"]))].shape[0]
            print(f"  PartType4 exists with {n4:,} particles.")
            print("  -> converted gas is a separate species. t7 can read it "
                  "directly:")
            print("       --conv-parttype 4")
        elif 1 in present and args.deep:
            pops = mass_populations(f["PartType1"]["Masses"][:])
            if len(pops) >= 2:
                lo = min(p[0] for p in pops)
                hi = max(p[0] for p in pops)
                print(f"  PartType1 has {len(pops)} mass populations, "
                      f"ratio {hi / lo:.3f}.")
                print("  -> converted gas is inside PartType1 and separable "
                      "by mass. Use:")
                print(f"       --conv-parttype 1 --conv-mass-max "
                      f"{0.5 * (lo + hi):.6e}")
            else:
                print("  PartType1 has a SINGLE mass population.")
                print("  -> converted gas is NOT separable by mass in this "
                      "snapshot.")
                print("     This is the case we expect for the new runs where "
                      "DM and baryons")
                print("     share a particle mass. Options, in order of "
                      "preference:")
                print("       1. re-run with a StarsScheme that keeps "
                      "PartType4")
                print("       2. write the QLA conversion flag to the "
                      "snapshot")
                print("       3. match particle IDs against the ICs: any ID "
                      "that was gas")
                print("          in the ICs and is now in PartType1 was "
                      "converted")
                print("     Option 3 needs no re-run and is what t7 falls "
                      "back to:")
                print("       --conv-from-ids <ICs.hdf5>")
        else:
            print("  run again with --deep to decide (needs the Masses "
                  "arrays).")

        print()
        if npart is not None:
            print(f"  NumPart_Total from Header: {np.ravel(npart)}")


if __name__ == "__main__":
    main()
