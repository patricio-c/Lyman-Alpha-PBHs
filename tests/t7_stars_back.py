#!/usr/bin/env python3
"""
Test 7 - put the converted gas back into the forest.

This is the test that closes the argument.  Everything upstream says the
suppression lives in the gas field and not in the matter field: the 3D
matter power ratio is 1.000 while the gas ratio is 0.808.  The difference
between those two numbers is entirely the baryons that the quick-Lyman-alpha
scheme turned into collisionless particles - 49.8% of them in FCT against
14.1% in CDM.

So: take those particles, call them gas again, and re-extract the forest.
If the ratio moves back towards 1, the signal we measured is a property of
the star-formation prescription, not of the primordial spectrum.  QLA
converts half the baryons while the observed stellar fraction at z=3 is a
couple of percent, so this is not a fringe systematic.

WHAT THIS SCRIPT DOES
    It writes an augmented snapshot in which PartType0 contains the
    surviving gas PLUS the converted particles, with a smoothing length and
    a temperature assigned to the latter.  It does NOT compute optical
    depths.  You then point the normal extraction at the new file:

        python tests/t7_stars_back.py --run fct40 --z 3.0 --out aug_fct.hdf5
        python stages/01_extract_los.py --snap aug_fct.hdf5 ...

    Keeping it this way means t7 is validated by exactly the same extractor
    as everything else, with no parallel code path to trust.

BEFORE YOU RUN IT
    python stages/00_inspect_snapshot.py --run fct40 --z 3.0 --deep

    That tells you which of the three identification modes applies.  In the
    current 40 Mpc/h pair the converted gas is inside PartType1 and
    separable by mass.  In the runs where dark matter and baryons share a
    particle mass it is not, and you need --conv-from-ids.

Usage
-----
    # dry run first - reports what it would do, writes nothing
    python tests/t7_stars_back.py --run fct40 --z 3.0 --dry-run --deep

    # mass-separated case (current runs)
    python tests/t7_stars_back.py --run fct40 --z 3.0 \\
        --conv-parttype 1 --conv-mass-max 2.0e-3 --out aug_fct.hdf5

    # PartType4 case
    python tests/t7_stars_back.py --run NB_1 --z 3.0 \\
        --conv-parttype 4 --out aug_NB1.hdf5

    # equal-mass case: recover the converted set from the ICs
    python tests/t7_stars_back.py --run NB_1 --z 3.0 \\
        --conv-from-ids /path/to/ics.hdf5 --out aug_NB1.hdf5

Options
-------
    --conv-parttype N     where the converted particles live (1 or 4)
    --conv-mass-max X     in internal units; particles in --conv-parttype
                          lighter than this are treated as converted gas
    --conv-from-ids PATH  ICs file; anything whose ParticleID was PartType0
                          in the ICs and is not PartType0 now was converted
    --h-mode grid|knn     how to assign a smoothing length (default grid)
    --ngrid N             grid size for --h-mode grid (default 512)
    --nngb N              neighbours for --h-mode knn (default 48)
    --t-mode trho|fixed   temperature of the reinjected gas (default trho)
    --t0 K                T at mean density for --t-mode trho (default 1.0e4)
    --gamma G             slope of the T-rho relation (default 1.55)
    --t-fixed K           temperature for --t-mode fixed
    --t0-sweep A B C      write one file per factor on T0, to show the
                          result does not depend on this choice
    --frac F              reinject only a random fraction F of the converted
                          particles, to build a curve instead of a point
    --dry-run             report and exit
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import runs  # noqa: E402

# SPH kernel normalisation: h = ETA * (m/rho)^(1/3) for a cubic spline with
# the SWIFT default resolution eta = 1.2348.
ETA = 1.2348


def first_dataset(g):
    return next(iter(g))


def count(g):
    return g[first_dataset(g)].shape[0]


# --- identifying the converted particles -----------------------------------

def converted_mask(f, args):
    """Return (parttype, boolean mask) selecting the converted gas."""
    if args.conv_from_ids:
        with h5py.File(args.conv_from_ids, "r") as ic:
            if "PartType0" not in ic:
                raise SystemExit("the ICs have no PartType0")
            gas_ids = np.asarray(ic["PartType0"]["ParticleIDs"][:])
        gas_ids.sort()
        pt = args.conv_parttype or 1
        ids = np.asarray(f[f"PartType{pt}"]["ParticleIDs"][:])
        idx = np.searchsorted(gas_ids, ids)
        idx[idx >= gas_ids.size] = 0
        mask = gas_ids[idx] == ids
        print(f"[t7] ID match against the ICs: {mask.sum():,} of {ids.size:,} "
              f"PartType{pt} particles were gas in the ICs")
        return pt, mask

    pt = args.conv_parttype
    if pt is None:
        raise SystemExit("give --conv-parttype, or --conv-from-ids")
    g = f[f"PartType{pt}"]
    n = count(g)
    if pt == 4:
        print(f"[t7] taking all {n:,} PartType4 particles as converted gas")
        return pt, np.ones(n, dtype=bool)
    if args.conv_mass_max is None:
        raise SystemExit("--conv-parttype 1 needs --conv-mass-max "
                         "(run stage 00 --deep to get the value)")
    m = np.asarray(g["Masses"][:])
    mask = m < args.conv_mass_max
    print(f"[t7] mass cut m < {args.conv_mass_max:.4e}: "
          f"{mask.sum():,} of {n:,} PartType{pt} particles selected "
          f"({100 * mask.mean():.2f}%)")
    if mask.all() or not mask.any():
        raise SystemExit("the mass cut selected everything or nothing - "
                         "check the value against stage 00 --deep")
    return pt, mask


# --- density and smoothing length for the reinjected particles -------------

def density_on_grid(pos_all, mass_all, box, ngrid, pos_query):
    """
    CIC density at the query positions from the full baryon distribution.

    Two passes over the particles, no tree, O(N).  Accurate enough for a
    smoothing length: h goes as rho^(-1/3), so a 20% density error is a 6%
    error on h, which is far below the effect we are testing.
    """
    dx = box / ngrid
    grid = np.zeros((ngrid, ngrid, ngrid), dtype=np.float64)

    chunk = 8_000_000
    for s in range(0, pos_all.shape[0], chunk):
        p = pos_all[s:s + chunk] / dx
        m = mass_all[s:s + chunk]
        i0 = np.floor(p).astype(np.int64)
        w1 = p - i0
        w0 = 1.0 - w1
        i0 %= ngrid
        i1 = (i0 + 1) % ngrid
        for a in (0, 1):
            wa = w0[:, 0] if a == 0 else w1[:, 0]
            ia = i0[:, 0] if a == 0 else i1[:, 0]
            for bb in (0, 1):
                wb = w0[:, 1] if bb == 0 else w1[:, 1]
                ib = i0[:, 1] if bb == 0 else i1[:, 1]
                for c in (0, 1):
                    wc = w0[:, 2] if c == 0 else w1[:, 2]
                    ic = i0[:, 2] if c == 0 else i1[:, 2]
                    np.add.at(grid,
                              (ia, ib, ic),
                              m * wa * wb * wc)
        print(f"    CIC {min(s + chunk, pos_all.shape[0]):,}/"
              f"{pos_all.shape[0]:,}", flush=True)

    grid /= dx ** 3
    q = pos_query / dx
    i0 = np.floor(q).astype(np.int64)
    w1 = q - i0
    w0 = 1.0 - w1
    i0 %= ngrid
    i1 = (i0 + 1) % ngrid
    out = np.zeros(q.shape[0])
    for a in (0, 1):
        wa = w0[:, 0] if a == 0 else w1[:, 0]
        ia = i0[:, 0] if a == 0 else i1[:, 0]
        for bb in (0, 1):
            wb = w0[:, 1] if bb == 0 else w1[:, 1]
            ib = i0[:, 1] if bb == 0 else i1[:, 1]
            for c in (0, 1):
                wc = w0[:, 2] if c == 0 else w1[:, 2]
                ic = i0[:, 2] if c == 0 else i1[:, 2]
                out += grid[ia, ib, ic] * wa * wb * wc
    return out


def density_knn(pos_all, mass_all, pos_query, nngb, box):
    from scipy.spatial import cKDTree
    print(f"    building tree on {pos_all.shape[0]:,} points", flush=True)
    tree = cKDTree(pos_all, boxsize=box)
    d, _ = tree.query(pos_query, k=nngb, workers=-1)
    hsml = d[:, -1]
    return mass_all.mean() * nngb / (4.0 / 3.0 * np.pi * hsml ** 3), hsml


# --- writing ---------------------------------------------------------------

def build_augmented(path_in, path_out, args, t0_factor=1.0):
    with h5py.File(path_in, "r") as f:
        box = float(np.ravel(f["Header"].attrs["BoxSize"])[0])
        pt, mask = converted_mask(f, args)
        gc = f[f"PartType{pt}"]
        g0 = f["PartType0"]
        n_gas = count(g0)
        n_conv = int(mask.sum())

        if args.frac < 1.0:
            rng = np.random.default_rng(args.seed)
            keep = rng.random(n_conv) < args.frac
            idx = np.flatnonzero(mask)[keep]
            n_conv = idx.size
            print(f"[t7] reinjecting a random {100 * args.frac:.0f}% "
                  f"-> {n_conv:,} particles")
        else:
            idx = np.flatnonzero(mask)

        print(f"[t7] gas now      : {n_gas:,}")
        print(f"[t7] to reinject  : {n_conv:,}")
        print(f"[t7] baryons total: {n_gas + n_conv:,}")

        if args.dry_run:
            print("[t7] --dry-run, nothing written")
            return None

        pos_conv = np.asarray(gc["Coordinates"][:])[idx]
        mass_conv = np.asarray(gc["Masses"][:])[idx] \
            if "Masses" in gc else np.full(n_conv,
                                           float(np.asarray(g0["Masses"][:1])[0]))
        vel_conv = np.asarray(gc["Velocities"][:])[idx] \
            if "Velocities" in gc else np.zeros((n_conv, 3))

        pos_gas = np.asarray(g0["Coordinates"][:])
        mass_gas = np.asarray(g0["Masses"][:])
        pos_all = np.concatenate([pos_gas, pos_conv])
        mass_all = np.concatenate([mass_gas, mass_conv])

        print(f"[t7] estimating density at the reinjected positions "
              f"({args.h_mode})", flush=True)
        if args.h_mode == "knn":
            rho_conv, hsml_conv = density_knn(pos_all, mass_all, pos_conv,
                                              args.nngb, box)
        else:
            rho_conv = density_on_grid(pos_all, mass_all, box, args.ngrid,
                                       pos_conv)
            rho_conv = np.maximum(rho_conv, 1e-12)
            hsml_conv = ETA * (mass_conv / rho_conv) ** (1.0 / 3.0)

        rho_mean = mass_all.sum() / box ** 3
        delta = rho_conv / rho_mean
        print(f"[t7] reinjected overdensity: median {np.median(delta):.2f}, "
              f"90th pct {np.percentile(delta, 90):.2f}")

        if args.t_mode == "fixed":
            T_conv = np.full(n_conv, args.t_fixed)
        else:
            T_conv = (args.t0 * t0_factor) * delta ** (args.gamma - 1.0)
        print(f"[t7] T of the reinjected gas: median "
              f"{np.median(T_conv):.3e} K   (T0 x {t0_factor:g})")

        # --- write ---------------------------------------------------------
        print(f"[t7] writing {path_out}", flush=True)
        with h5py.File(path_out, "w") as o:
            for grp in f:
                if grp.startswith("PartType"):
                    continue
                f.copy(grp, o)
            og = o.create_group("PartType0")

            for name in g0:
                d = g0[name]
                shape = (n_gas + n_conv,) + d.shape[1:]
                out = og.create_dataset(name, shape=shape, dtype=d.dtype,
                                        compression="gzip" if args.compress
                                        else None, shuffle=args.compress)
                out[:n_gas] = d[:]
                for k, v in d.attrs.items():
                    out.attrs[k] = v

                low = name.lower()
                if low == "coordinates":
                    out[n_gas:] = pos_conv
                elif low == "velocities":
                    out[n_gas:] = vel_conv
                elif low == "masses":
                    out[n_gas:] = mass_conv
                elif low in ("smoothinglengths", "smoothinglength"):
                    out[n_gas:] = hsml_conv
                elif low in ("densities", "density"):
                    out[n_gas:] = rho_conv
                elif low in ("temperatures", "temperature"):
                    out[n_gas:] = T_conv
                elif low in ("internalenergies", "internalenergy"):
                    u_gas = np.asarray(d[:n_gas])
                    T_gas_proxy = np.median(u_gas)
                    out[n_gas:] = T_gas_proxy * (T_conv / np.median(T_conv))
                    print("    [!] InternalEnergies filled by scaling the "
                          "median gas value; if the extractor uses u rather "
                          "than T, check the mu assumption")
                elif low == "particleids":
                    out[n_gas:] = np.asarray(gc["ParticleIDs"][:])[idx]
                else:
                    med = np.median(np.asarray(d[:min(n_gas, 100000)]), axis=0)
                    out[n_gas:] = med
                    print(f"    [.] {name}: filled with the gas median")

            hdr = o["Header"].attrs
            npart = np.array(np.ravel(hdr.get("NumPart_Total",
                                              np.zeros(7, dtype=np.int64))),
                             dtype=np.int64).copy()
            if npart.size > pt:
                npart[0] = n_gas + n_conv
                npart[pt] = max(int(npart[pt]) - n_conv, 0)
                hdr["NumPart_Total"] = npart
                if "NumPart_ThisFile" in hdr:
                    hdr["NumPart_ThisFile"] = npart
            hdr["T7_ReinjectedParticles"] = n_conv
            hdr["T7_SourcePartType"] = pt
            hdr["T7_TemperatureMode"] = np.bytes_(args.t_mode)
            hdr["T7_T0"] = args.t0 * t0_factor
            hdr["T7_Gamma"] = args.gamma
            hdr["T7_Fraction"] = args.frac
            hdr.pop("Virtual", None)

    print(f"[t7] done -> {path_out}")
    return path_out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default=None)
    ap.add_argument("--z", type=float, default=None)
    ap.add_argument("--snap", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--conv-parttype", type=int, default=None)
    ap.add_argument("--conv-mass-max", type=float, default=None)
    ap.add_argument("--conv-from-ids", default=None)
    ap.add_argument("--h-mode", default="grid", choices=["grid", "knn"])
    ap.add_argument("--ngrid", type=int, default=512)
    ap.add_argument("--nngb", type=int, default=48)
    ap.add_argument("--t-mode", default="trho", choices=["trho", "fixed"])
    ap.add_argument("--t0", type=float, default=1.0e4)
    ap.add_argument("--gamma", type=float, default=1.55)
    ap.add_argument("--t-fixed", type=float, default=1.0e4)
    ap.add_argument("--t0-sweep", type=float, nargs="*", default=None)
    ap.add_argument("--frac", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--compress", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--deep", action="store_true",
                    help="ignored, accepted so the stage-00 flag can be "
                         "pasted through")
    args = ap.parse_args()

    if args.run:
        path_in = runs.resolve_snapshot(args.run, z=args.z, snap=args.snap)
    elif args.snap:
        path_in = args.snap
    else:
        ap.error("give --run or --snap")

    if not args.dry_run and not args.out:
        ap.error("--out is required unless --dry-run")

    factors = args.t0_sweep if args.t0_sweep else [1.0]
    for fac in factors:
        out = args.out
        if out and len(factors) > 1:
            stem, ext = os.path.splitext(out)
            out = f"{stem}_T0x{fac:g}{ext}"
        print("=" * 70)
        print(f"T0 factor {fac:g}  ->  {out}")
        print("=" * 70)
        build_augmented(path_in, out, args, t0_factor=fac)


if __name__ == "__main__":
    main()
