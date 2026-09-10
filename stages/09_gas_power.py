#!/usr/bin/env python3
"""
Stage 09 - the 3D power spectrum of the gas, which SWIFT did not write.

The one measurement that closes the tilt argument with no free parameters.

WHY IT IS NEEDED

The 40 Mpc/h pair was run with `requested_spectra: ["matter-matter",
"cdm-cdm"]`, so there is no gas spectrum on disk.  It cannot be recovered
from the two that exist either: with

    P_mm = f_g^2 P_gg + 2 f_g f_c P_gc + f_c^2 P_cc

two measured spectra leave two unknowns.  But the snapshots hold every
particle, so the field can simply be built.

Do NOT use SWIFT's `cdm-cdm` to compare runs.  In a QLA run `PartType1` is
dark matter PLUS the converted baryons - 18.9M particles in CDM and 66.9M
in FCT - so it is a different physical quantity in each run.  Only
`matter-matter` divides cleanly, and this stage reproduces it as a check.

WHAT IT MEASURES

Four fields, deposited separately from the same snapshot:

    gas        PartType0                          -> P_gg   THE measurement
    converted  PartType1, IDs that were gas in the ICs -> P_cc
    dm         PartType1, everything else
    matter     all of the above                   -> P_mm   the control

and from them the cross spectrum P_gc and the all-baryon field
P_bb = (gas + converted), which is what the baryons would look like had
the conversion sink never removed anything from the gas phase.

    P_gg(FCT)/P_gg(CDM)   should be suppressed at low k if the large-scale
                          P1D deficit is the baryon drain
    P_bb(FCT)/P_bb(CDM)   should sit near the matter ratio, ~1 at low k

The difference between those two is the effect of the conversion on the
absorbing field, scale by scale, with no tau_eff, no assumed temperature,
no smoothing length and no sightlines involved.

CAVEAT ON P_cc AND P_bb

The converted particles are collisionless and have moved since they were
converted, so P_cc is not "what that gas would have done" and P_bb is an
approximation, not a counterfactual.  **P_gg is the exact one** - it is
the real gas, measured. Quote that; treat the rest as supporting.

THE CONTROL THAT COMES FIRST

P_mm is compared against SWIFT's own `power_matter_*.txt` for the same
snapshot.  If they do not agree, this pipeline is wrong and nothing below
it should be believed.  Pass --swift-pk to enable it; the run says loudly
whether the check passed.

ORDER OF OPERATIONS

The deposits and the five auto spectra are the expensive part.  They are
written to disk, and the control is run, BEFORE the cross spectrum is
attempted - a first version did the cross spectrum first and lost two
completed runs to it.  P_gc feeds only the all-baryon decomposition, which
is supporting; losing it costs nothing that matters.

Usage
-----
    python stages/09_gas_power.py \\
        --snap .../cdm-40-m6-lyman_0003/cdm-40-m6-lyman_0003.hdf5 \\
        --ics  .../cdm-40-m6-lyman_0000/cdm-40-m6-lyman_0000.hdf5 \\
        --label CDM --ngrid 512 \\
        --swift-pk .../cdm-box-40-1024/power_spectra/power_matter_0010.txt \\
        --out figures/pk_gas_cdm_z3

Options
-------
    --snap PATH       the snapshot to measure
    --ics PATH        an earlier snapshot whose PartType0 IDs define the
                      baryons; the ICs give 100% coverage
    --label NAME      display name
    --ngrid N         grid side (default 512).  4 bytes per cell per field,
                      and four fields are held at once
    --mas SCHEME      NGP | CIC | TSC | PCS (default CIC)
    --threads N       FFT threads (default 8)
    --chunk N         particles read at a time (default 8e6)
    --swift-pk PATH   SWIFT's power_matter_*.txt for the same snapshot, to
                      validate P_mm against
    --no-cross        skip the gas x converted cross spectrum
    --out PREFIX      writes PREFIX.txt (and PREFIX.npz with every array)
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import prov  # noqa: E402

try:
    import MAS_library as MASL
    import Pk_library as PKL
except ImportError as exc:                                   # noqa: BLE001
    raise SystemExit(
        f"Pylians is required for this stage and did not import: {exc}\n"
        f"conda activate astro, then check "
        f"`python -c 'import MAS_library, Pk_library'`.")


def sorted_ids(path, group, chunk):
    """Every ParticleID of `group`, sorted."""
    with h5py.File(path, "r") as f:
        d = f[group]["ParticleIDs"]
        out = np.empty(d.shape[0], dtype=np.uint64)
        for i0 in range(0, d.shape[0], chunk):
            i1 = min(i0 + chunk, d.shape[0])
            out[i0:i1] = d[i0:i1]
    out.sort()
    return out


def in_sorted(ids, ref):
    """Boolean mask: which of `ids` appear in the sorted array `ref`."""
    if len(ref) == 0:
        return np.zeros(len(ids), dtype=bool)
    pos = np.searchsorted(ref, ids)
    np.clip(pos, 0, len(ref) - 1, out=pos)
    return ref[pos] == ids


def deposit(path, group, grid, box, mas, chunk, mask_ids=None, want=True):
    """
    CIC-deposit particle MASS onto `grid`.  If `mask_ids` is given, only
    particles whose ID is in it (want=True) or not in it (want=False) are
    deposited.  Returns the mass actually deposited and sum of m^2.
    """
    m_tot = 0.0
    m2_tot = 0.0
    n_used = 0
    with h5py.File(path, "r") as f:
        g = f[group]
        n = g["Coordinates"].shape[0]
        for i0 in range(0, n, chunk):
            i1 = min(i0 + chunk, n)
            p = np.asarray(g["Coordinates"][i0:i1], dtype=np.float32)
            m = np.asarray(g["Masses"][i0:i1], dtype=np.float32)
            if mask_ids is not None:
                ids = np.asarray(g["ParticleIDs"][i0:i1], dtype=np.uint64)
                sel = in_sorted(ids, mask_ids)
                if not want:
                    sel = ~sel
                p, m = p[sel], m[sel]
                del ids, sel
            if len(p):
                MASL.MA(p, grid, box, mas, W=m)
                m_tot += float(np.float64(m).sum())
                m2_tot += float((np.float64(m) ** 2).sum())
                n_used += len(p)
            del p, m
    return m_tot, m2_tot, n_used


def to_delta(mass_grid):
    """Mass grid -> overdensity, in place. Fails loudly on an empty field."""
    mean = mass_grid.mean(dtype=np.float64)
    if not np.isfinite(mean) or mean <= 0:
        raise SystemExit("a field is empty or non-finite; nothing deposited. "
                         "Refusing to divide by it.")
    mass_grid /= np.float32(mean)
    mass_grid -= np.float32(1.0)
    return mass_grid


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snap", required=True)
    ap.add_argument("--ics", required=True)
    ap.add_argument("--label", default="run")
    ap.add_argument("--ngrid", type=int, default=512)
    ap.add_argument("--mas", default="CIC",
                    choices=["NGP", "CIC", "TSC", "PCS"])
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--swift-pk", default=None)
    ap.add_argument("--no-cross", action="store_true",
                    help="skip the gas x converted cross spectrum")
    ap.add_argument("--out", default="figures/pk_gas")
    args = ap.parse_args()

    log = list(prov.header())

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    for p in (args.snap, args.ics):
        if not os.path.isfile(p):
            raise SystemExit(f"{p} is not a file")

    with h5py.File(args.snap, "r") as f:
        box = float(np.ravel(f["Header"].attrs["BoxSize"])[0])
        z = float(np.ravel(f["Header"].attrs["Redshift"])[0])
        npt = np.atleast_1d(f["Header"].attrs["NumPart_Total"])
    with h5py.File(args.ics, "r") as f:
        z_ic = float(np.ravel(f["Header"].attrs["Redshift"])[0])
        npi = np.atleast_1d(f["Header"].attrs["NumPart_Total"])

    if z_ic <= z:
        raise SystemExit(f"--ics is at z={z_ic} and --snap at z={z}; the ICs "
                         f"file has to be the earlier one.")
    if int(npi[0]) + int(npi[1]) != int(npt[0]) + int(npt[1]):
        raise SystemExit(
            f"gas + DM is {int(npi[0]) + int(npi[1]):,} in the ICs and "
            f"{int(npt[0]) + int(npt[1]):,} now. QLA conversion moves "
            f"particles between types and never creates or destroys them, "
            f"so the converted set would not mean what it claims.")

    say("=" * 74)
    say(f"Stage 09 - 3D power spectra    [{args.label}]   z = {z:.4f}")
    say("=" * 74)
    say(f"box {box:.5f} (internal length units; k comes out in their inverse)")
    say(f"grid {args.ngrid}^3   MAS {args.mas}   ICs at z = {z_ic:.4f}")
    say(f"k_fund = {2 * np.pi / box:.6f}   "
        f"k_Nyquist = {np.pi * args.ngrid / box:.4f}")

    say("\nidentifying the converted particles by ID ...")
    ic_gas = sorted_ids(args.ics, "PartType0", args.chunk)
    now_gas = sorted_ids(args.snap, "PartType0", args.chunk)
    doomed = ic_gas[~in_sorted(ic_gas, now_gas)]
    del ic_gas, now_gas
    doomed.sort()
    say(f"  converted since the ICs: {len(doomed):,}"
        f"   ({100 * len(doomed) / max(int(npi[0]), 1):.2f}% of the IC gas)")

    shape = (args.ngrid, args.ngrid, args.ngrid)
    say(f"\ndepositing (each field is "
        f"{4 * args.ngrid ** 3 / 1024 ** 3:.2f} GB) ...")

    fields = {}
    stats = {}

    g_gas = np.zeros(shape, dtype=np.float32)
    stats["gas"] = deposit(args.snap, "PartType0", g_gas, box, args.mas,
                           args.chunk)
    say(f"  gas       {stats['gas'][2]:>12,} particles")

    g_conv = np.zeros(shape, dtype=np.float32)
    stats["conv"] = deposit(args.snap, "PartType1", g_conv, box, args.mas,
                            args.chunk, mask_ids=doomed, want=True)
    say(f"  converted {stats['conv'][2]:>12,} particles")
    if stats["conv"][2] != len(doomed):
        say(f"  [!] only {stats['conv'][2]:,} of {len(doomed):,} converted "
            f"IDs were found in PartType1. Something other than QLA "
            f"conversion removed gas; the split below is incomplete.")

    g_dm = np.zeros(shape, dtype=np.float32)
    stats["dm"] = deposit(args.snap, "PartType1", g_dm, box, args.mas,
                          args.chunk, mask_ids=doomed, want=False)
    say(f"  dark      {stats['dm'][2]:>12,} particles")
    del doomed

    m_gas, m_conv, m_dm = (stats[k][0] for k in ("gas", "conv", "dm"))
    m_all = m_gas + m_conv + m_dm
    m2_all = sum(stats[k][1] for k in ("gas", "conv", "dm"))
    say(f"\nmass fractions   gas {m_gas / m_all:.6f}   "
        f"converted {m_conv / m_all:.6f}   dark {m_dm / m_all:.6f}")
    say(f"baryon fraction  {(m_gas + m_conv) / m_all:.6f}"
        f"   (Omega_b/Omega_m, a sanity check on the split)")

    # matter and all-baryons are sums of the mass grids, before normalising
    g_mat = g_gas + g_conv + g_dm
    g_bar = g_gas + g_conv

    p_shot_mat = box ** 3 * m2_all / (m_all ** 2)
    say(f"matter shot noise = {p_shot_mat:.6e} (same length units cubed)")

    for nm, arr in (("gas", g_gas), ("conv", g_conv), ("dm", g_dm),
                    ("matter", g_mat), ("baryon", g_bar)):
        fields[nm] = to_delta(arr)

    say("\ncomputing spectra ...")
    out = {}

    def auto(nm):
        pk = PKL.Pk(fields[nm], box, 0, args.mas, args.threads, False)
        out["k"] = np.asarray(pk.k3D, dtype=np.float64)
        out["Nmodes"] = np.asarray(pk.Nmodes3D, dtype=np.float64)
        out[f"P_{nm}"] = np.asarray(pk.Pk[:, 0], dtype=np.float64)
        say(f"  P_{nm} done")

    for nm in ("gas", "conv", "dm", "matter", "baryon"):
        auto(nm)

    # Everything above this line is the expensive part - the deposits and
    # five FFTs over 6.7e8 particles. Free what the cross spectrum does not
    # need, run the control, and get the results onto disk BEFORE touching
    # anything else. A first version of this stage did the cross spectrum
    # first and lost a completed run to it.
    for nm in ("dm", "matter", "baryon"):
        fields.pop(nm, None)

    # ------------------------------------------------------- the control
    if args.swift_pk:
        say("\n" + "-" * 74)
        say("CONTROL: our P_matter against SWIFT's own")
        say("-" * 74)
        s = np.loadtxt(args.swift_pk)
        zs = float(s[0, 0])
        if abs(zs - z) > 1e-3:
            raise SystemExit(f"--swift-pk is at z={zs} but the snapshot is at "
                             f"z={z}. Comparing different epochs would make "
                             f"the control meaningless.")
        ks, ps = s[:, 1], s[:, 2]              # column 2 is shot-subtracted
        k = out["k"]
        band = (k > 2.0 * 2 * np.pi / box) & (k < 0.5 * np.pi * args.ngrid / box)
        ours = out["P_matter"][band] - p_shot_mat
        theirs = np.interp(k[band], ks, ps)
        r = ours / theirs
        say(f"  compared over {band.sum()} bins, "
            f"k = {k[band][0]:.4f} to {k[band][-1]:.4f}")
        say(f"  ratio ours/SWIFT   median {np.median(r):.4f}   "
            f"5th {np.percentile(r, 5):.4f}   95th {np.percentile(r, 95):.4f}")
        if abs(np.median(r) - 1.0) < 0.02:
            say("  PASS - the deposit, the MAS correction and the shot-noise "
                "subtraction all reproduce SWIFT.")
        else:
            say("  [!] FAIL - do not use any spectrum from this run until "
                "this is understood. Check MAS, the shot-noise term, and "
                "whether SWIFT's column 2 is really shot-subtracted.")

    def write_out():
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                    exist_ok=True)
        cols = [c for c in ("k", "Nmodes", "P_gas", "P_conv", "P_dm",
                            "P_matter", "P_baryon", "P_gc") if c in out]
        arr = np.column_stack([out[c] for c in cols])
        hdr = (f"stage 09   label {args.label}   z {z:.6f}   box {box:.6f}\n"
               f"grid {args.ngrid}  MAS {args.mas}\n"
               f"mass fractions gas {m_gas / m_all:.8f} "
               f"conv {m_conv / m_all:.8f} dm {m_dm / m_all:.8f}\n"
               f"matter shot noise {p_shot_mat:.8e}\n"
               + "  ".join(cols))
        np.savetxt(args.out + ".txt", arr, header=hdr)
        np.savez(args.out + ".npz", z=z, box=box, ngrid=args.ngrid,
                 f_gas=m_gas / m_all, f_conv=m_conv / m_all,
                 f_dm=m_dm / m_all, shot_matter=p_shot_mat, **out)
        with open(args.out + ".log", "w") as fh:
            fh.write("\n".join(log) + "\n")
        print(f"written -> {args.out}.txt / .npz / .log  "
              f"({len(cols)} columns)", flush=True)

    write_out()

    # --------------------------------------------- the cross spectrum, last
    # P_gc only feeds the all-baryon decomposition, which is the supporting
    # number and not the measurement. It goes after the results are safe.
    if args.no_cross:
        say("\n--no-cross: skipping the gas x converted cross spectrum.")
    else:
        say("\ncross spectrum (gas x converted), the last and least "
            "important step ...")
        try:
            xpk = PKL.XPk([fields["gas"], fields["conv"]], box, 0,
                          MAS=[args.mas, args.mas], threads=args.threads)
            xk = np.asarray(xpk.XPk)
            out["P_gc"] = np.asarray(
                xk[:, 0, 0] if xk.ndim == 3 else xk[:, 0],
                dtype=np.float64)
            say("  P_gc done")
            write_out()
        except Exception as exc:                            # noqa: BLE001
            say(f"  [!] the cross spectrum failed: "
                f"{exc.__class__.__name__}: {exc}")
            say("      Everything else is already on disk. P_gc is only "
                "needed for the all-baryon decomposition; P_gg, the "
                "measurement, does not depend on it. Re-run with "
                "--no-cross, or with a larger --mem, to silence this.")
    print(f"\ndone -> {args.out}.txt")


if __name__ == "__main__":
    main()
