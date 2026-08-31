#!/usr/bin/env python3
"""
Test 12 - does our LOS regeneration reproduce SWIFT's own output?

Everything downstream starts from a file of LOS_XXXX groups holding the gas
particles near each sightline.  Two things can produce that file: SWIFT
itself, on the fly, or `legacy/relos.py`, which walks a snapshot and selects
the particles with impact parameter b < gamma*h around each ray.

t7 needs the second one - it shoots sightlines through an augmented
snapshot, which SWIFT never saw.  So the question "do we have full control
of the sightlines once we hold the snapshot?" has to be answered before any
t7 number is worth anything.  This test answers it by taking a file SWIFT
wrote itself, regenerating it from the snapshot at the same redshift and the
same ray positions, and requiring the two to hold the SAME PARTICLES.

Read this before choosing the input files
-----------------------------------------
The comparison is only meaningful against a SWIFT LOS file that is itself
correct.  The 40 Mpc/h pair's original files were written with
`range_when_shooting_down_* = [0, 40]` in a 58.7372 box, so each sightline is
missing about 32% of its particles, and the current files were regenerated
with relos.py.  Comparing relos.py against relos.py output proves nothing,
and comparing it against the truncated originals fails by construction.
Use a run whose LOS output was written correctly on the fly.

What it reports
---------------
For each file on its own:

    number of sightlines, box, redshift, kernel gamma, fields present,
    particles per sightline, and the RANGE of the ray positions.

That last one is a diagnostic in its own right.  If the ray positions of a
run span [0, 40] while the box is 58.74, that run sampled 46% of the
transverse face, and two runs where one is truncated and the other is not do
not sample the same volume - which is a candidate explanation for stage 02's
verdict that the sightlines of the two production boxes are unrelated.

For the pair:

    do the group names match; do the ray positions match; do the particle
    COUNTS match; and, the real test, are the sets of ParticleIDs per
    sightline identical.  For the particles present in both, the field
    values are compared as well, so a selection that is right but a units
    conversion that is wrong still shows up.

If the LOS files carry no ParticleIDs the comparison falls back to sorted
coordinates, which is weaker: it cannot tell two particles at the same place
apart.  Say so in the write-up if that is the path taken.

Usage
-----
    # describe one file, including the ray-position range
    python tests/t12_relos_roundtrip.py --a los_0010.hdf5

    # the round trip
    python legacy/relos.py --old-los NATIVE.hdf5 \\
        --snapshot 'snap_0003*.hdf5' --out /tmp/regen.hdf5 --max-los 100
    python tests/t12_relos_roundtrip.py --a NATIVE.hdf5 --b /tmp/regen.hdf5

Options
-------
    --a FILE          the reference file (SWIFT's own output)
    --b FILE          the regenerated file. Omit to describe --a only.
    --max-los N       compare only the first N sightlines
    --tol-pos F       ray positions count as equal below F * box (1e-6)
    --tol-val F       field values count as equal below this relative
                      difference (1e-6; they are copies, so exact is the
                      expectation and anything else is a conversion bug)
"""

from __future__ import annotations

import argparse
import os

import h5py
import numpy as np


def describe(path, max_los=None):
    """Everything about one LOS file that does not need a second one."""
    out = {"path": path}
    with h5py.File(path, "r") as f:
        names = sorted(k for k in f.keys() if k.startswith("LOS_"))
        if max_los:
            names = names[:max_los]
        if not names:
            raise SystemExit(f"{path}: no LOS_XXXX groups")
        out["names"] = names
        out["box"] = float(np.ravel(f["Header"].attrs["BoxSize"])[0])
        out["z"] = (float(np.ravel(f["Cosmology"].attrs["Redshift"])[0])
                    if "Cosmology" in f else float("nan"))
        try:
            out["gamma"] = float(np.ravel(
                f["HydroScheme"].attrs["Kernel gamma"])[0])
        except (KeyError, ValueError):
            out["gamma"] = float("nan")
        out["fields"] = sorted(f[names[0]].keys())
        out["has_ids"] = "ParticleIDs" in out["fields"]

        counts, pos, axes = [], [], []
        for nm in names:
            g = f[nm]
            c = g["Coordinates"][:].astype(np.float64)
            counts.append(c.shape[0])
            # the integration axis is the one the particles run along
            ang = 2.0 * np.pi * c / out["box"]
            R = np.hypot(np.cos(ang).mean(axis=0), np.sin(ang).mean(axis=0))
            ax = int(np.argmin(R))
            axes.append(ax)
            tr = [i for i in range(3) if i != ax]
            p = None
            for kx, ky in (("Xpos", "Ypos"), ("Xproj", "Yproj"),
                           ("x_proj", "y_proj")):
                if kx in g.attrs and ky in g.attrs:
                    p = np.array([float(np.ravel(g.attrs[kx])[0]),
                                  float(np.ravel(g.attrs[ky])[0])])
                    break
            if p is None:
                h = g["SmoothingLengths"][:]
                idx = np.argsort(h)[:max(10, len(h) // 10)]
                p = np.median(c[idx][:, tr], axis=0)
            pos.append(p)
        out["counts"] = np.array(counts)
        out["pos"] = np.array(pos)
        out["axes"] = np.array(axes)
    return out


def report_one(d, say):
    say(f"  file        {d['path']}")
    say(f"  sightlines  {len(d['names'])}")
    say(f"  z           {d['z']:.5f}   box = {d['box']:.5f} internal   "
        f"kernel gamma = {d['gamma']:.6f}")
    say(f"  axes        {sorted(set(d['axes'].tolist()))}")
    say(f"  fields      {d['fields']}")
    ids = ("present" if d["has_ids"] else
           "ABSENT - the comparison falls back to coordinates, weaker")
    say(f"  ParticleIDs {ids}")
    c = d["counts"]
    say(f"  particles per LOS   min {c.min()}  median {int(np.median(c))}  "
        f"max {c.max()}  total {c.sum()}")
    lo, hi = d["pos"].min(axis=0), d["pos"].max(axis=0)
    say(f"  ray positions       axis0 [{lo[0]:.4f}, {hi[0]:.4f}]   "
        f"axis1 [{lo[1]:.4f}, {hi[1]:.4f}]   of a {d['box']:.4f} box")
    frac = float(np.prod((hi - lo) / d["box"]))
    say(f"  transverse coverage {100 * frac:.1f}% of the box face")
    if hi.max() < 0.85 * d["box"]:
        say(f"  [!] the rays do not reach the far side of the box. This is")
        say(f"      the signature of range_when_shooting_down_* being set")
        say(f"      narrower than the box. Such a file samples a sub-volume,")
        say(f"      and two runs where only one is truncated do not sample")
        say(f"      the same region even if both are fair samples.")


def compare(A, B, tol_pos, tol_val, say):
    """The round trip itself.  Returns True if it is a clean reproduction."""
    ok = True
    box = A["box"]

    say("--- structure " + "-" * 58)
    if abs(A["box"] - B["box"]) / box > 1e-9:
        say(f"  [FAIL] different box: {A['box']} vs {B['box']}")
        return False
    if abs(A["z"] - B["z"]) > 0.02:
        say(f"  [FAIL] different redshift: {A['z']:.5f} vs {B['z']:.5f}")
        return False
    n = min(len(A["names"]), len(B["names"]))
    if A["names"][:n] != B["names"][:n]:
        say("  [FAIL] the LOS group names differ")
        return False
    say(f"  OK   {n} sightlines, same names, same box, same z")

    d = A["pos"][:n] - B["pos"][:n]
    d -= box * np.round(d / box)
    sep = np.hypot(d[:, 0], d[:, 1])
    if sep.max() > tol_pos * box:
        say(f"  [FAIL] ray positions differ, max {sep.max():.6g} internal "
            f"({100 * sep.max() / box:.4f}% of the box)")
        say("         relos.py was probably run with --uniform, which draws")
        say("         new rays instead of reusing the ones in --old-los.")
        ok = False
    else:
        say(f"  OK   ray positions identical, max offset {sep.max():.3e}")

    say()
    say("--- particle counts " + "-" * 52)
    ca, cb = A["counts"][:n], B["counts"][:n]
    same = int((ca == cb).sum())
    say(f"  identical on {same}/{n} sightlines ({100 * same / n:.2f}%)")
    if same != n:
        diff = (cb.astype(np.int64) - ca.astype(np.int64))
        say(f"  count difference B-A   min {diff.min()}  "
            f"median {int(np.median(diff))}  max {diff.max()}")
        say(f"  total A {ca.sum()}   total B {cb.sum()}   "
            f"ratio {cb.sum() / max(ca.sum(), 1):.6f}")
        ok = False

    say()
    say("--- particle identity, the actual test " + "-" * 33)
    if not (A["has_ids"] and B["has_ids"]):
        say("  ParticleIDs absent in at least one file. Falling back to")
        say("  sorted coordinates: this cannot distinguish two particles at")
        say("  the same position, so a pass here is weaker evidence.")
    exact, only_a, only_b, worst_val = 0, 0, 0, 0.0
    worst_field = ""
    with h5py.File(A["path"], "r") as fa, h5py.File(B["path"], "r") as fb:
        shared = [x for x in A["fields"] if x in B["fields"]
                  and x != "ParticleIDs"]
        for nm in A["names"][:n]:
            ga, gb = fa[nm], fb[nm]
            if A["has_ids"] and B["has_ids"]:
                ia = np.asarray(ga["ParticleIDs"][:]).ravel()
                ib = np.asarray(gb["ParticleIDs"][:]).ravel()
                sa, sb = set(ia.tolist()), set(ib.tolist())
                only_a += len(sa - sb)
                only_b += len(sb - sa)
                if sa == sb:
                    exact += 1
                # field values on the particles both files agree about
                common = np.array(sorted(sa & sb), dtype=ia.dtype)
                if common.size:
                    oa, ob = np.argsort(ia), np.argsort(ib)
                    pa = oa[np.searchsorted(ia[oa], common)]
                    pb = ob[np.searchsorted(ib[ob], common)]
                    for fl in shared:
                        va = np.asarray(ga[fl][:])[pa].astype(np.float64)
                        vb = np.asarray(gb[fl][:])[pb].astype(np.float64)
                        den = np.maximum(np.abs(va), 1e-300)
                        r = float(np.max(np.abs(va - vb) / den))
                        if r > worst_val:
                            worst_val, worst_field = r, fl
            else:
                xa = np.sort(np.asarray(ga["Coordinates"][:]).astype(
                    np.float64), axis=0)
                xb = np.sort(np.asarray(gb["Coordinates"][:]).astype(
                    np.float64), axis=0)
                if xa.shape == xb.shape and np.allclose(xa, xb, rtol=0,
                                                        atol=1e-8 * box):
                    exact += 1

    say(f"  sightlines whose particle set matches exactly: "
        f"{exact}/{n} ({100 * exact / n:.2f}%)")
    if A["has_ids"] and B["has_ids"]:
        say(f"  particles only in A (SWIFT, missed by us): {only_a}")
        say(f"  particles only in B (ours, SWIFT did not have): {only_b}")
        if shared:
            say(f"  worst relative field difference on shared particles: "
                f"{worst_val:.3e}  ({worst_field})")
            if worst_val > tol_val:
                say("  [FAIL] the two files disagree on the VALUES of "
                    "particles they both contain.")
                say("         The selection may be right and a units or")
                say("         conversion step wrong. Look there, not at the")
                say("         geometry.")
                ok = False
    if exact != n:
        ok = False

    say()
    say("--- verdict " + "-" * 60)
    if ok:
        say("  PASS. The regeneration reproduces SWIFT's own sightlines")
        say("  particle for particle. We have full control of the LOS once")
        say("  we hold the snapshot, so t7 can shoot through an augmented")
        say("  snapshot and the result means what it says.")
    else:
        say("  FAIL. Do not use relos.py output as if it were SWIFT's until")
        say("  this is understood. The numbers above say which of the three")
        say("  it is: wrong rays, wrong selection, or wrong values.")
    return ok


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", default=None)
    ap.add_argument("--max-los", type=int, default=None)
    ap.add_argument("--tol-pos", type=float, default=1e-6)
    ap.add_argument("--tol-val", type=float, default=1e-6)
    ap.add_argument("--out", default=None, help="also write the log here")
    args = ap.parse_args()

    log = []

    def say(s=""):
        print(s)
        log.append(s)

    say("=" * 72)
    say("T12 - does our LOS regeneration reproduce SWIFT's own output?")
    say("=" * 72)

    A = describe(args.a, args.max_los)
    say("--- file A (the reference) " + "-" * 45)
    report_one(A, say)
    say()

    ok = None
    if args.b:
        B = describe(args.b, args.max_los)
        say("--- file B (the regeneration) " + "-" * 41)
        report_one(B, say)
        say()
        ok = compare(A, B, args.tol_pos, args.tol_val, say)
    else:
        say("  No --b given: description only, nothing compared.")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                    exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("\n".join(log) + "\n")
        print(f"written -> {args.out}")
    raise SystemExit(0 if ok in (None, True) else 1)


if __name__ == "__main__":
    main()
