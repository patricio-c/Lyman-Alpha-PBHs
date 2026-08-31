#!/usr/bin/env python3
"""
Stage 01 - SWIFT LOS file to a tau cache.

This is the only file in the repository that touches the extraction code.
It reimplements nothing: it imports `open_los_file` and `extract_tau` from
legacy/swift_extract.py, which is the code every published number came
from, and drives them.

WHAT IT ADDS ON TOP

  1. File resolution by redshift, so a run that did not exist when this was
     written works without editing code (--run NB_1 --z 3.0). Snapshot and
     LOS indices are not comparable between runs and nothing here uses one.

  2. The UV background as an input (--treecool / --gamma-hi) instead of a
     constant buried in a script.

  3. Two convention switches that were previously implicit:

     --hubble-mode file|fixed|ref
         H(z) from Cosmology/H in the LOS file, from a cosmology stated on
         the command line, or from a reference run so two boxes share a k
         axis by construction. H sets box_kms and therefore dv, so it
         rescales the whole k axis: if two runs were extracted with
         different H, part of their ratio is a stretch and not physics.

     --vel-convention swift|gadget|none
         The extractor converts Velocities with the file's own "physical
         CGS including cosmological corrections" attribute. That is the
         SWIFT convention and it is the correct one. `gadget` multiplies by
         sqrt(a) - a factor 2 at z=3 - and exists so the size of that error
         can be measured instead of argued about. `none` drops peculiar
         velocities altogether: not physical, but it isolates how much of
         the ratio redshift-space distortions carry.

     Both work by overriding the metadata object and wrapping the unit
     helper. Nothing in legacy/ is edited.

  4. Provenance, and the extractor's own diagnostics aggregated over all
     sightlines instead of scrolling past.

Usage
-----
    # which LOS files does this run have, and at what redshift?
    python stages/01_extract_los.py --run fct40

    # the geometry, without reading a particle
    python stages/01_extract_los.py --run fct40 --z 3.0 --geometry-only

    # production
    python stages/01_extract_los.py --run fct40 --z 3.0 --npix 2048 \\
        --treecool data/TREECOOL_HM12_G+Q --out cache/cache_fct.npz

    # the convention tests
    python stages/01_extract_los.py --run fct40 --z 3.0 \\
        --vel-convention gadget --out cache/cache_fct_gadget.npz
    python stages/01_extract_los.py --run fct40 --z 3.0 \\
        --hubble-mode ref --reference-run cdm40 \\
        --out cache/cache_fct_Hcdm.npz

    # the physics cuts, passed straight through to extract_tau
    python stages/01_extract_los.py --run fct40 --z 3.0 --delta-max 100 \\
        --out cache/cache_fct_d100.npz
    python stages/01_extract_los.py --run fct40 --z 3.0 \\
        --impose-trho 1.0e4 1.40 --out cache/cache_fct_iso140.npz
    python stages/01_extract_los.py --run fct40 --z 3.0 --no-normalize \\
        --out cache/cache_fct_nonorm.npz
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "legacy"))

from common import cache as cachelib, runs, units, uvb  # noqa: E402


def load_extractor():
    try:
        import swift_extract
    except ImportError as exc:
        raise SystemExit(f"could not import legacy/swift_extract.py: {exc}\n"
                         "Run scripts/migrate.sh first.")
    for fn in ("open_los_file", "extract_tau", "_cgs_factor"):
        if not hasattr(swift_extract, fn):
            raise SystemExit(f"legacy/swift_extract.py has no {fn}(); this "
                             "stage was written against the v2 API.")
    return swift_extract


def patch_velocity_convention(mod, factor):
    """
    Scale the Velocities unit conversion by `factor`.

    extract_tau reads velocities through _cgs_factor(g["Velocities"]), so
    wrapping that one function changes the peculiar velocity field and
    nothing else. That is why this switch needs no fork of the extractor.
    """
    if factor == 1.0:
        return
    original = mod._cgs_factor

    def wrapped(dset, physical=True):
        f = original(dset, physical=physical)
        if dset.name.rsplit("/", 1)[-1] == "Velocities":
            f *= factor
        return f

    mod._cgs_factor = wrapped


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    g = ap.add_argument_group("input")
    g.add_argument("--run", default=None, help="registry name or directory")
    g.add_argument("--los-file", default=None, help="explicit LOS hdf5")
    g.add_argument("--z", type=float, default=None)
    g.add_argument("--npix", type=int, default=2048)
    g.add_argument("--max-los", type=int, default=None,
                   help="stop after N sightlines, for a quick check")

    g = ap.add_argument_group("UV background")
    g.add_argument("--treecool", default=None)
    g.add_argument("--gamma-hi", type=float, default=None,
                   help="explicit Gamma_HI in s^-1, overrides --treecool")

    g = ap.add_argument_group("conventions")
    g.add_argument("--hubble-mode", default="file",
                   choices=["file", "fixed", "ref"])
    g.add_argument("--reference-run", default="cdm40")
    g.add_argument("--h0", type=float, default=68.1)
    g.add_argument("--omega-m", type=float, default=0.3053)
    g.add_argument("--vel-convention", default="swift",
                   choices=["swift", "gadget", "none"])

    g = ap.add_argument_group("physics, passed through to extract_tau")
    g.add_argument("--xh", type=float, default=0.76)
    g.add_argument("--he-state", default="HeII")
    g.add_argument("--n-sigma", type=float, default=8.0)
    g.add_argument("--exact-voigt", action="store_true")
    g.add_argument("--no-normalize", action="store_true",
                   help="skip the Shepard correction to the SPH deposition")
    g.add_argument("--w-floor", type=float, default=0.05)
    g.add_argument("--delta-max", type=float, default=None)
    g.add_argument("--impose-trho", type=float, nargs=2, default=None,
                   metavar=("T0", "GAMMA"))
    g.add_argument("--max-b-violation", type=float, default=0.02)

    g = ap.add_argument_group("output")
    g.add_argument("--geometry-only", action="store_true")
    g.add_argument("--progress-every", type=int, default=100)
    g.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.los_file:
        path = os.path.abspath(args.los_file)
    elif args.run:
        path = runs.resolve_los_file(args.run, z=args.z)
    else:
        ap.error("give --run or --los-file")

    sw = load_extractor()
    meta = sw.open_los_file(path)

    hz_file = meta.hz
    if args.hubble_mode == "fixed":
        cosmo = units.Cosmology(h=args.h0 / 100.0, Om=args.omega_m)
        meta.hz = cosmo.H(meta.z)
        hsrc = f"fixed h={cosmo.h:.4f} Om={cosmo.Om:.4f}"
    elif args.hubble_mode == "ref":
        rpath = runs.resolve_los_file(args.reference_run, z=meta.z,
                                      verbose=False)
        meta.hz = sw.open_los_file(rpath).hz
        hsrc = f"reference run {args.reference_run}"
    else:
        hsrc = "Cosmology/H in the LOS file"

    a = meta.a
    vfac = {"swift": 1.0, "gadget": float(np.sqrt(a)),
            "none": 0.0}[args.vel_convention]
    gHI, gsrc = uvb.gamma_hi(meta.z, table=args.treecool, value=args.gamma_hi)
    dv_pred = meta.dv(args.npix)

    print("=" * 72)
    print(f"file        {path}")
    print(f"z           {meta.z:.5f}   a = {a:.6f}   h = {meta.h:.4f}")
    print(f"Omega_m     {meta.omega_m:.5f}   Omega_b = {meta.omega_b:.5f}")
    print(f"box         {meta.box_comoving_mpc:.5f} Mpc comoving = "
          f"{meta.box_comoving_mpc * meta.h:.4f} Mpc/h")
    print(f"H(z)        {meta.hz:.4f} km/s/Mpc   source: {hsrc}")
    if abs(meta.hz - hz_file) / hz_file > 1e-9:
        print(f"            file says {hz_file:.4f}; ratio "
              f"{meta.hz / hz_file:.6f}, and the k axis stretches by that")
    print(f"v_box       {meta.box_kms:.4f} km/s")
    print(f"dv          {dv_pred:.6f} km/s  ({args.npix} px)")
    print(f"k_fund      {2 * np.pi / meta.box_kms:.6f} s/km")
    print(f"k_nyq       {np.pi / dv_pred:.6f} s/km  "
          f"(analyse below {0.5 * np.pi / dv_pred:.6f})")
    print(f"k[s/km] -> Mpc^-1  x {meta.hz / (1 + meta.z):.4f}     "
          f"-> h/Mpc  x {meta.hz / (1 + meta.z) / meta.h:.4f}")
    print(f"DESI window {units.desi_window(meta.z)[0]:.5f} - "
          f"{units.desi_window(meta.z)[1]:.5f} s/km")
    print(f"Gamma_HI    {gHI:.5e} s^-1   source: {gsrc}")
    print(f"velocities  {args.vel_convention}  (factor {vfac:.6f})")
    print(f"sightlines  {len(meta.los_names)}"
          + (f", using the first {args.max_los}" if args.max_los else ""))
    print(f"normalize   {not args.no_normalize}"
          + (f"   delta_max={args.delta_max}" if args.delta_max else "")
          + (f"   impose_trho={tuple(args.impose_trho)}"
             if args.impose_trho else ""))
    print("=" * 72)

    if args.geometry_only:
        return

    patch_velocity_convention(sw, vfac)

    names = meta.los_names[:args.max_los] if args.max_los else meta.los_names
    if not names:
        raise SystemExit(f"{path}: no LOS_XXXX groups")

    kw = dict(gamma_HI=gHI, meta=meta, X_H=args.xh, he_state=args.he_state,
              n_sigma=args.n_sigma, exact_voigt=args.exact_voigt,
              max_b_violation=args.max_b_violation,
              normalize=not args.no_normalize, w_floor=args.w_floor,
              delta_max=args.delta_max,
              impose_trho=tuple(args.impose_trho) if args.impose_trho
              else None)

    tau = np.empty((len(names), args.npix), dtype=np.float32)
    diags = []
    dv = dv_pred
    t0 = time.time()
    for i, nm in enumerate(names):
        t, dv, d = sw.extract_tau(path, nm, args.npix, **kw)
        tau[i] = t
        diags.append(d)
        if args.progress_every and (i % args.progress_every == 0
                                    or i == len(names) - 1):
            el = time.time() - t0
            eta = el / (i + 1) * (len(names) - i - 1)
            print(f"  {i + 1}/{len(names)}   {el / 60:.1f} min elapsed, "
                  f"{eta / 60:.1f} min left", flush=True)

    def agg(key):
        v = np.array([d[key] for d in diags], dtype=np.float64)
        return v.min(), np.median(v), v.max()

    print("\n--- extractor diagnostics over all sightlines " + "-" * 25)
    for key, why in (
            ("frac_b_gt_H", "particles outside their own kernel support, "
                            "must be ~0"),
            ("max_b_over_H", "worst impact parameter over support"),
            ("wsum_raw_med", "SPH partition of unity before Shepard; 1.0 is "
                             "perfect sampling"),
            ("frac_below_floor", "pixels where the Shepard floor bit"),
            ("frac_delta_cut", "fraction removed by --delta-max")):
        lo, md, hi = agg(key)
        print(f"  {key:20s} min {lo:10.4g}  med {md:10.4g}  max {hi:10.4g}"
              f"   {why}")
    print(f"  ray position source: {sorted({d['ray_source'] for d in diags})}")
    print(f"  integration axis   : {sorted({d['axis'] for d in diags})}")

    te = -np.log(np.mean(np.exp(-tau.astype(np.float64))))
    print(f"\ntau_eff (raw, A=1) = {te:.5f}")
    print(f"target Turner+24   = {units.tau_eff_turner24(meta.z):.5f}")

    out = args.out or f"cache/cache_{args.run or 'run'}.npz"
    cachelib.save(
        out, tau, dv, meta.z, inputs=[path],
        box_mpc=meta.box_comoving_mpc, run=str(args.run),
        engine="swift_extract_v2", uvb=gsrc, gamma_hi=gHI,
        hubble_mode=args.hubble_mode, hubble=meta.hz, hubble_file=hz_file,
        vel_convention=args.vel_convention,
        normalize=not args.no_normalize,
        delta_max=-1.0 if args.delta_max is None else args.delta_max,
        impose_trho=np.array(args.impose_trho if args.impose_trho
                             else [-1.0, -1.0]),
        los_names=np.array(names, dtype="S16"), tau_eff_raw=te)
    print(f"written -> {out}")


if __name__ == "__main__":
    main()
