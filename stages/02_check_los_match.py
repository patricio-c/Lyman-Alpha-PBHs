#!/usr/bin/env python3
"""
Stage 02 - are the two runs sampling the same sightlines?

Everything paired in this analysis rests on an assumption that was never
checked: that LOS_0007 in CDM and LOS_0007 in FCT are the same line through
the same structure.  The sightline positions come from SWIFT's LOS output
configuration, not from anything we control at analysis time, so the
assumption can fail in three ways - different positions, the same positions
in a different order, or a different number of lines.

This matters for two results:

    t8_single_los.py    compares line i against line i pixel by pixel. If
                        the lines are not matched the comparison is noise.
    the paired bootstrap  quoted 20 sigma. Pairing removes cosmic variance
                        only if the pairs are real.

The check reads the ray position of every line from both files using the
same logic as the extractor (the group attributes when present, the median
transverse position of the smallest-kernel particles otherwise), and then:

    1. compares them index by index, with periodic wrapping;
    2. if that fails, looks for a permutation that does match, so you can
       reorder rather than re-extract.

Usage
-----
    python stages/02_check_los_match.py --a <cdm LOS.hdf5> --b <fct LOS.hdf5>
    python stages/02_check_los_match.py --run-a cdm40 --run-b fct40 --z 3.0

Options
-------
    --run-a / --run-b NAME   registry names, resolved by redshift
    --a / --b PATH           explicit LOS files, override --run-*
    --z Z                    redshift, for --run-*
    --tol-frac F             a pair counts as matched when the transverse
                             offset is below F times the box (default 1e-4)
    --save-permutation PATH  if a permutation is found, write it as .npy so
                             the caches can be reordered instead of redone
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "legacy"))

from common import runs  # noqa: E402


def ray_positions(path, sw):
    """Transverse ray position of every LOS group, in internal units."""
    meta = sw.open_los_file(path)
    out, srcs, axes = [], set(), set()
    with h5py.File(path, "r") as f:
        for nm in meta.los_names:
            g = f[nm]
            coords = g["Coordinates"][:].astype(np.float64)
            hsml = g["SmoothingLengths"][:].astype(np.float64)
            axis = sw.detect_los_axis(coords, meta.boxsize_int)
            tr = [i for i in range(3) if i != axis]
            pos, src = sw._ray_position(g, coords[:, tr],
                                        hsml * meta.kernel_gamma,
                                        meta.boxsize_int)
            out.append(pos)
            srcs.add(src)
            axes.add(axis)
    return np.array(out), meta, srcs, axes


def periodic_sep(pa, pb, box):
    d = pa - pb
    d -= box * np.round(d / box)
    return np.hypot(d[..., 0], d[..., 1])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-a", default=None)
    ap.add_argument("--run-b", default=None)
    ap.add_argument("--a", default=None)
    ap.add_argument("--b", default=None)
    ap.add_argument("--z", type=float, default=None)
    ap.add_argument("--tol-frac", type=float, default=1e-4)
    ap.add_argument("--save-permutation", default=None)
    args = ap.parse_args()

    pa = args.a or (runs.resolve_los_file(args.run_a, z=args.z)
                    if args.run_a else None)
    pb = args.b or (runs.resolve_los_file(args.run_b, z=args.z)
                    if args.run_b else None)
    if not (pa and pb):
        ap.error("give --a/--b or --run-a/--run-b with --z")

    sys.path.insert(0, os.path.join(ROOT, "legacy"))
    import swift_extract as sw

    print("=" * 72)
    print(f"A  {pa}")
    print(f"B  {pb}")
    print("=" * 72)

    ra, ma, sa, xa = ray_positions(pa, sw)
    rb, mb, sb, xb = ray_positions(pb, sw)
    box = ma.boxsize_int
    tol = args.tol_frac * box

    print(f"A: {len(ra)} lines, ray source {sa}, axis {xa}, "
          f"box {ma.boxsize_int:.5f}, z {ma.z:.4f}")
    print(f"B: {len(rb)} lines, ray source {sb}, axis {xb}, "
          f"box {mb.boxsize_int:.5f}, z {mb.z:.4f}")
    if abs(ma.boxsize_int - mb.boxsize_int) / box > 1e-9:
        print("[!] different box sizes. Nothing below is meaningful.")
    if xa != xb:
        print("[!] different integration axes. The lines cannot be matched.")
    print(f"tolerance = {tol:.6g} internal units "
          f"({args.tol_frac:g} of the box)")
    print()

    if len(ra) != len(rb):
        print(f"[!] different number of sightlines: {len(ra)} vs {len(rb)}.")
        print("    Truncate to the common length only if you can show the "
              "first N are the same lines - the check below does that.")
    n = min(len(ra), len(rb))

    d = periodic_sep(ra[:n], rb[:n], box)
    frac_ok = float((d < tol).mean())
    print("--- index by index " + "-" * 52)
    print(f"  matched fraction   {100 * frac_ok:.2f}%")
    print(f"  offset  median     {np.median(d):.6g}")
    print(f"          90th pct   {np.percentile(d, 90):.6g}")
    print(f"          max        {d.max():.6g}   "
          f"({100 * d.max() / box:.4f}% of the box)")

    if frac_ok > 0.999:
        print("\nVERDICT: the sightlines are matched line by line.")
        print("  t8_single_los.py is valid as written, and the paired "
              "bootstrap is a genuine pairing.")
        return

    print("\n[!] the lines are NOT matched by index. Looking for a "
          "permutation.")
    perm = np.empty(n, dtype=np.int64)
    dist = np.empty(n)
    for i in range(n):
        dd = periodic_sep(ra[i][None, :], rb[:n], box)
        j = int(np.argmin(dd))
        perm[i] = j
        dist[i] = dd[j]
    unique = len(np.unique(perm)) == n
    good = float((dist < tol).mean())
    print(f"  nearest-neighbour match: {100 * good:.2f}% within tolerance, "
          f"one-to-one = {unique}")

    if unique and good > 0.999:
        print("\nVERDICT: same sightlines, different order.")
        print("  Reorder B rather than re-extracting: tau_B[perm] lines up "
              "with tau_A.")
        if args.save_permutation:
            np.save(args.save_permutation, perm)
            print(f"  written -> {args.save_permutation}")
        else:
            print("  Pass --save-permutation to write it out.")
        return

    print("\nVERDICT: the two runs do not sample the same sightlines.")
    print("  Consequences, in order of severity:")
    print("   - t8_single_los.py compares unrelated lines. Do not use it "
          "until this is fixed.")
    print("   - the paired bootstrap is not paired; the 20 sigma has to be "
          "requoted from an unpaired test, which will be weaker.")
    print("   - the ENSEMBLE results (P1D, its ratio, the 3D power) are "
          "unaffected: they never assumed a pairing.")
    print("  Fix by re-running the SWIFT LOS output for both boxes with the "
          "same seed / same position list.")


if __name__ == "__main__":
    main()
