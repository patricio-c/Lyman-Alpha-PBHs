#!/usr/bin/env python3
"""
Cut a small, self-contained LOS file out of a big one.

Why: the SpectWizard cross-check (HANDOFF open question 5) needs Maria to
run her extractor on EXACTLY the same sightlines as ours.  The production
LOS files are 2-4 GB, which is not something to move around or ask a
colleague to wade through.  A hundred sightlines is a few tens of MB, holds
the same groups with the same attributes and the same unit-conversion
metadata, and is enough to compare tau pixel by pixel.

The output is a valid SWIFT LOS file: every root group is copied verbatim
(Header, Cosmology, Units, HydroScheme, InternalCodeUnits, ...), the LOS
groups keep their original names and attributes, and the particle counts in
the header are corrected to match what was actually written.  Anything that
reads the original reads this.

Selection is strided by default rather than the first N in a row.  The
sightlines are written in the order SWIFT drew them, which for a truncated
file is a corner of the box; taking the first hundred would hand over a
hundred neighbours instead of a hundred samples.

Usage
-----
    python scripts/make_los_subset.py --in los_0010.hdf5 \\
        --out los_0010_sample100.hdf5 --n 100

    # the first N in file order, if you specifically want that
    python scripts/make_los_subset.py --in ... --out ... --n 100 --first

Options
-------
    --in FILE     the full LOS file
    --out FILE    where to write
    --n N         how many sightlines (default 100)
    --first       take the first N instead of striding through the file
    --seed S      with --random, the draw (default 1)
    --random      pick N at random instead of striding
"""

from __future__ import annotations

import argparse
import os

import h5py
import numpy as np

ROOT_GROUPS = ["Header", "Cosmology", "Units", "HydroScheme",
               "InternalCodeUnits", "GravityScheme", "ICs_parameters",
               "Code", "Parameters", "PolicyFlags", "RuntimePars",
               "SubgridScheme", "StarsScheme", "Policy"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--first", action="store_true")
    ap.add_argument("--random", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    with h5py.File(args.inp, "r") as src:
        names = sorted(k for k in src.keys() if k.startswith("LOS_"))
        if not names:
            raise SystemExit(f"{args.inp}: no LOS_XXXX groups")
        n_all = len(names)
        if args.n >= n_all:
            pick = names
        elif args.first:
            pick = names[:args.n]
        elif args.random:
            rng = np.random.default_rng(args.seed)
            pick = [names[i] for i in
                    sorted(rng.choice(n_all, size=args.n, replace=False))]
        else:
            step = n_all / float(args.n)
            pick = [names[int(i * step)] for i in range(args.n)]

        print(f"source      {args.inp}")
        print(f"sightlines  {n_all} available, taking {len(pick)}"
              + ("  (first)" if args.first else
                 "  (random)" if args.random else "  (strided)"))

        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                    exist_ok=True)
        with h5py.File(args.out, "w") as out:
            copied = []
            for grp in src.keys():
                if grp.startswith("LOS_"):
                    continue
                if not isinstance(src[grp], h5py.Group):
                    continue
                src.copy(grp, out)
                copied.append(grp)
            missing = [g for g in ROOT_GROUPS if g in src and g not in copied]
            if missing:
                print(f"[!] root groups not copied: {missing}")

            total = 0
            counts = []
            for nm in pick:
                src.copy(nm, out)
                c = out[nm]["Coordinates"].shape[0]
                counts.append(c)
                total += c

            if "Header" in out:
                npt = np.zeros(7, dtype=np.int64)
                npt[0] = total
                out["Header"].attrs["NumPart_ThisFile"] = npt
                out["Header"].attrs["NumPart_Total"] = npt.astype(np.uint32)
                if "TotalNumberOfParticles" in out["Header"].attrs:
                    out["Header"].attrs["TotalNumberOfParticles"] = npt

    counts = np.array(counts)
    size_mb = os.path.getsize(args.out) / 1024 ** 2
    print(f"root groups copied: {copied}")
    print(f"particles per LOS   min {counts.min()}  "
          f"median {int(np.median(counts))}  max {counts.max()}")
    print(f"total particles     {total}")
    print(f"written -> {args.out}   ({size_mb:.1f} MB)")
    print()
    print("This file is self-contained: it carries the same Header,")
    print("Cosmology and Units groups as the original, so any code that")
    print("reads the full file reads this one unchanged.")


if __name__ == "__main__":
    main()
