#!/usr/bin/env python3
"""
Stage 08 - where are the converted particles now?

Stage 07 asks where the doomed gas was before it was converted.  This asks
the other half of Joop's question: where did those same particles END UP.
Are they locked inside collapsed haloes, or spread through the field?

THE PROBLEM, AND THE WAY ROUND IT

Once converted, a particle lives in `PartType1`, which carries no
`Densities` and no `Temperatures` - only positions, masses, IDs,
potentials, softenings and velocities.  So "where are they" cannot be read
off the snapshot; it has to be measured.  This stage builds a
cloud-in-cell density field from ALL the matter at the late snapshot and
samples it at the converted particles' positions.

Everything is expressed as Delta = rho / rho_mean of the same field, so no
unit conversion enters anywhere: the mean is just the total deposited mass
over the number of cells.  A bug in the units cannot bias the answer
because there are no units.

THE COMPARISON THAT MATTERS

The converted Delta distribution on its own says little - matter is
clustered, so any tracer sits at Delta > 1 on average.  What matters is
the converted distribution AGAINST a random sample of PartType1 measured
the same way, on the same grid.  If the converted particles sit at much
higher Delta than ordinary matter, they are in haloes.  If the two
distributions lie on top of each other, conversion has left them
dynamically indistinguishable from the rest and "which haloes" is the
wrong question.

Sampling at particle positions is mass-weighted and carries shot noise
from the finite number of particles per cell.  That affects the converted
population and the reference identically, which is why the RATIO of the
two medians is the statistic to quote and the raw medians are not.

RESOLUTION IS A FLOOR, NOT A DETAIL

Delta measured on a grid saturates at the cell scale: a halo smaller than
a cell is smeared out and reads low.  The cell size is printed, and the
verdict should never be read past it.  Raise --ngrid to push the floor
down, at 8 bytes per cell.

Usage
-----
    python stages/08_where_they_went.py \\
        --early .../murgia-cdm-lyman_0000/murgia-cdm-lyman_0000.hdf5 \\
                .../murgia-M3-lyman_0000/murgia-M3-lyman_0000.hdf5 \\
        --late  .../murgia-cdm-lyman_0002/murgia-cdm-lyman_0002.hdf5 \\
                .../murgia-M3-lyman_0002/murgia-M3-lyman_0002.hdf5 \\
        --labels cdm M3 --ngrid 512 --out figures/went_murgia_z5

Options
-------
    --early ...       one snapshot per run; defines the doomed set together
                      with --late, exactly as stage 07 does
    --late ...        one per run, same order; the epoch everything is
                      measured at
    --labels ...      display names, same order
    --ngrid N         cells per side of the density grid (default 512)
    --nsample N       PartType1 particles sampled for the reference
                      distribution (default 5e6)
    --slab F          thickness of the projection slab, as a fraction of
                      the box (default 0.02)
    --chunk N         particles read at a time (default 8e6)
    --out PREFIX      writes PREFIX.png and PREFIX.txt
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

D_LO, D_HI, NB = 1e-2, 1e6, 160


def h5ids(path, group, chunk):
    """Every ParticleID of `group`, sorted."""
    with h5py.File(path, "r") as f:
        d = f[group]["ParticleIDs"]
        out = np.empty(d.shape[0], dtype=np.uint64)
        for i0 in range(0, d.shape[0], chunk):
            i1 = min(i0 + chunk, d.shape[0])
            out[i0:i1] = d[i0:i1]
    out.sort()
    return out


def missing_from(early, late_sorted, chunk):
    """Gas IDs present at `early` and absent from the sorted `late_sorted`."""
    got = []
    with h5py.File(early, "r") as f:
        d = f["PartType0"]["ParticleIDs"]
        for i0 in range(0, d.shape[0], chunk):
            i1 = min(i0 + chunk, d.shape[0])
            ids = np.asarray(d[i0:i1], dtype=np.uint64)
            pos = np.searchsorted(late_sorted, ids)
            np.clip(pos, 0, len(late_sorted) - 1, out=pos)
            got.append(ids[late_sorted[pos] != ids])
    out = np.concatenate(got) if got else np.empty(0, dtype=np.uint64)
    out.sort()
    return out


def cic_deposit(grid, pos, box, weight=1.0):
    """Add particles at `pos` (N,3) onto `grid` (n,n,n) with CIC weights."""
    n = grid.shape[0]
    x = (pos / box) * n
    np.mod(x, n, out=x)
    i0 = np.floor(x).astype(np.int64)
    d = x - i0
    i0 %= n
    i1 = (i0 + 1) % n
    flat = grid.reshape(-1)
    for cx in (0, 1):
        wx = d[:, 0] if cx else 1.0 - d[:, 0]
        ix = i1[:, 0] if cx else i0[:, 0]
        for cy in (0, 1):
            wy = d[:, 1] if cy else 1.0 - d[:, 1]
            iy = i1[:, 1] if cy else i0[:, 1]
            wxy = wx * wy
            for cz in (0, 1):
                wz = d[:, 2] if cz else 1.0 - d[:, 2]
                iz = i1[:, 2] if cz else i0[:, 2]
                idx = (ix * n + iy) * n + iz
                flat += np.bincount(idx, weights=wxy * wz * weight,
                                    minlength=flat.size)


def cic_sample(grid, pos, box):
    """CIC-interpolate `grid` at `pos` (N,3). Same kernel as the deposit."""
    n = grid.shape[0]
    x = (pos / box) * n
    np.mod(x, n, out=x)
    i0 = np.floor(x).astype(np.int64)
    d = x - i0
    i0 %= n
    i1 = (i0 + 1) % n
    out = np.zeros(len(pos))
    for cx in (0, 1):
        wx = d[:, 0] if cx else 1.0 - d[:, 0]
        ix = i1[:, 0] if cx else i0[:, 0]
        for cy in (0, 1):
            wy = d[:, 1] if cy else 1.0 - d[:, 1]
            iy = i1[:, 1] if cy else i0[:, 1]
            wxy = wx * wy
            for cz in (0, 1):
                wz = d[:, 2] if cz else 1.0 - d[:, 2]
                iz = i1[:, 2] if cz else i0[:, 2]
                out += wxy * wz * grid[ix, iy, iz]
    return out


def build_field(path, ngrid, chunk):
    """CIC density grid of all the matter, in units of the mean."""
    with h5py.File(path, "r") as f:
        box = float(np.ravel(f["Header"].attrs["BoxSize"])[0])
        grid = np.zeros((ngrid, ngrid, ngrid), dtype=np.float64)
        for pt in ("PartType0", "PartType1"):
            if pt not in f:
                continue
            g = f[pt]
            n = g["Coordinates"].shape[0]
            for i0 in range(0, n, chunk):
                i1 = min(i0 + chunk, n)
                p = np.asarray(g["Coordinates"][i0:i1], dtype=np.float64)
                m = np.asarray(g["Masses"][i0:i1], dtype=np.float64)
                cic_deposit(grid, p, box, m)
                del p, m
    mean = grid.mean()
    if mean <= 0:
        raise SystemExit(f"{path}: the density grid is empty. Nothing was "
                         f"deposited, so every Delta below would be garbage.")
    grid /= mean
    return grid, box


def coords_of(path, group, want_sorted, chunk):
    """Coordinates of the particles in `group` whose IDs are in `want`."""
    out = []
    with h5py.File(path, "r") as f:
        g = f[group]
        n = g["ParticleIDs"].shape[0]
        for i0 in range(0, n, chunk):
            i1 = min(i0 + chunk, n)
            ids = np.asarray(g["ParticleIDs"][i0:i1], dtype=np.uint64)
            pos = np.searchsorted(want_sorted, ids)
            np.clip(pos, 0, len(want_sorted) - 1, out=pos)
            hit = want_sorted[pos] == ids
            if hit.any():
                out.append(np.asarray(g["Coordinates"][i0:i1],
                                      dtype=np.float64)[hit])
            del ids, pos, hit
    return (np.concatenate(out) if out
            else np.empty((0, 3), dtype=np.float64))


def sample_ref(path, grid, box, nsample, chunk, rng):
    """Delta at the position of a random subset of PartType1."""
    with h5py.File(path, "r") as f:
        n = f["PartType1"]["Coordinates"].shape[0]
        take = min(nsample, n)
        # sampling WITH replacement on purpose: choice(replace=False) builds
        # a permutation of n, which is 4 GB when PartType1 has 5.6e8 entries.
        # For estimating a distribution the duplicates are harmless.
        sel = np.sort(rng.integers(0, n, size=take))
        out = []
        for i0 in range(0, n, chunk):
            i1 = min(i0 + chunk, n)
            k = sel[(sel >= i0) & (sel < i1)]
            if not len(k):
                continue
            p = np.asarray(f["PartType1"]["Coordinates"][i0:i1],
                           dtype=np.float64)[k - i0]
            out.append(cic_sample(grid, p, box))
            del p
    return np.concatenate(out) if out else np.empty(0)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--early", nargs="+", required=True)
    ap.add_argument("--late", nargs="+", required=True)
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--ngrid", type=int, default=512)
    ap.add_argument("--nsample", type=int, default=5_000_000)
    ap.add_argument("--slab", type=float, default=0.02)
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--seed", type=int, default=20260904)
    ap.add_argument("--out", default="figures/where_they_went")
    args = ap.parse_args()

    if len(args.early) != len(args.late):
        raise SystemExit("--early and --late must have the same count")
    labels = args.labels or [f"run{i}" for i in range(len(args.early))]
    if len(labels) != len(args.early):
        raise SystemExit("--labels must have one entry per run")

    rng = np.random.default_rng(args.seed)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    say("=" * 74)
    say("Stage 08 - where the converted particles ended up")
    say("=" * 74)
    say("")
    say("Delta is measured on a CIC grid of ALL the matter at the LATE")
    say("snapshot, in units of its own mean, and sampled at the converted")
    say("particles' positions. The reference is a random sample of")
    say("PartType1 measured on the same grid with the same kernel.")

    de = np.logspace(np.log10(D_LO), np.log10(D_HI), NB + 1)
    dc = np.sqrt(de[:-1] * de[1:])
    C = {}

    for lab, e, l in zip(labels, args.early, args.late):
        for p in (e, l):
            if not os.path.isfile(p):
                raise SystemExit(f"{p} is not a file")
        with h5py.File(e, "r") as f:
            ze = float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])
            ne = np.atleast_1d(f["Header"].attrs["NumPart_Total"])
        with h5py.File(l, "r") as f:
            zl = float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])
            nl = np.atleast_1d(f["Header"].attrs["NumPart_Total"])
        if zl >= ze:
            raise SystemExit(f"{lab}: --late (z={zl}) is not later than "
                             f"--early (z={ze}).")
        if int(ne[0]) + int(ne[1]) != int(nl[0]) + int(nl[1]):
            raise SystemExit(
                f"{lab}: gas + DM is {int(ne[0]) + int(ne[1]):,} early and "
                f"{int(nl[0]) + int(nl[1]):,} late. QLA conversion moves "
                f"particles between types and never creates or destroys "
                f"them, so the doomed set would not mean what it claims.")

        say(f"\n[{lab}]  z = {ze:.4f} -> {zl:.4f}")
        say("       identifying the converted particles ...")
        late_gas = h5ids(l, "PartType0", args.chunk)
        doomed = missing_from(e, late_gas, args.chunk)
        del late_gas
        say(f"       converted between these epochs: {len(doomed):,}")
        if not len(doomed):
            raise SystemExit(f"{lab}: nothing was converted between these two "
                             f"snapshots; there is nothing to locate.")

        say("       locating them in PartType1 ...")
        pos = coords_of(l, "PartType1", doomed, args.chunk)
        found = len(pos)
        say(f"       found {found:,} of {len(doomed):,} in PartType1"
            f"  ({100 * found / len(doomed):.2f}%)")
        if found < 0.99 * len(doomed):
            say("       [!] not all of them are there. Something other than")
            say("           QLA conversion removed gas in this run, and the")
            say("           distribution below is of whatever WAS found.")
        del doomed

        say(f"       building the {args.ngrid}^3 matter field ...")
        grid, box = build_field(l, args.ngrid, args.chunk)
        cell = box / args.ngrid
        say(f"       box {box:.4f}, cell {cell:.5f} (same length units)")

        d_conv = cic_sample(grid, pos, box)
        say("       sampling the reference ...")
        d_ref = sample_ref(l, grid, box, args.nsample, args.chunk, rng)

        # projection slab, converted particles and all matter alike
        half = 0.5 * args.slab * box
        mid = 0.5 * box
        sl = np.abs(pos[:, 2] - mid) < half
        k0 = int((mid - half) / box * args.ngrid)
        k1 = max(k0 + 1, int((mid + half) / box * args.ngrid))

        C[lab] = dict(
            h_conv=np.histogram(d_conv, bins=de)[0].astype(float),
            h_ref=np.histogram(d_ref, bins=de)[0].astype(float),
            d_conv=d_conv, d_ref=d_ref, cell=cell, box=box,
            proj_xy=pos[sl][:, :2],
            proj_mat=grid[:, :, k0:k1].mean(axis=2))
        del grid, pos, d_conv, d_ref

    say("\n" + "-" * 74)
    say("WHERE THE CONVERTED PARTICLES SIT, at the late snapshot")
    say("-" * 74)
    say(f"{'run':8s} {'median Delta':>13s} {'reference':>11s} {'ratio':>8s}"
        f" {'% >100':>8s} {'% >1000':>9s} {'% >1e4':>8s}")
    for lab in labels:
        c = C[lab]
        a, b = c["d_conv"], c["d_ref"]
        ma, mb = float(np.median(a)), float(np.median(b))
        say(f"{lab:8s} {ma:13.3f} {mb:11.3f} {ma / mb if mb else 0:8.3f}"
            f" {100 * (a > 1e2).mean():7.2f}% {100 * (a > 1e3).mean():8.2f}%"
            f" {100 * (a > 1e4).mean():7.2f}%")

    say("")
    say("The 'reference' column is the median Delta of ordinary PartType1 on")
    say("the same grid.  The ratio is the whole result: near 1 means the")
    say("converted particles are dynamically indistinguishable from the rest")
    say("of the matter and asking which haloes hold them is the wrong")
    say("question; well above 1 means they really are locked in collapsed")
    say("structure.")
    say("")
    say("Delta saturates at the cell scale, so nothing above roughly the")
    say("density of a halo filling one cell is resolved.  Cell sizes are")
    say("printed per run above; raise --ngrid to lower that floor and check")
    say("the ratio is stable before quoting it.")

    # The numbers are the expensive part - a 512^3 field over 6.7e8 particles
    # - so they go to disk BEFORE anything touches matplotlib, and a failure
    # in the plotting cannot cost a cluster run.
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"\nwritten -> {args.out}.txt")

    try:
        make_figure(labels, C, dc, args)
    except Exception as exc:                             # noqa: BLE001
        print(f"[!] the figure failed: {exc.__class__.__name__}: {exc}")
        print(f"    every number above is already saved in {args.out}.txt")
    else:
        print(f"written -> {args.out}.png")


def make_figure(labels, C, dc, args):
    nr = len(labels)
    fig = plt.figure(figsize=(12.5, 4.4 + 3.6 * ((nr + 1) // 2)))
    gs = fig.add_gridspec(1 + (nr + 1) // 2, 2, hspace=0.32, wspace=0.24)
    ax = fig.add_subplot(gs[0, 0])
    bx = fig.add_subplot(gs[0, 1])
    PAL = ["#1f4fbf", "#c0271f", "#2ca02c", "#9467bd", "#ff7f0e"]

    for i, lab in enumerate(labels):
        c = C[lab]
        col = PAL[i % len(PAL)]
        for h, ls, tag in ((c["h_conv"], "-", "converted"),
                           (c["h_ref"], ":", "all matter")):
            if h.sum() > 0:
                ax.plot(dc, h / h.sum(), color=col, lw=2 if ls == "-" else 1.3,
                        ls=ls, label=f"{lab}, {tag}")
        if c["h_conv"].sum() > 0:
            bx.plot(dc, np.cumsum(c["h_conv"][::-1])[::-1] / c["h_conv"].sum(),
                    color=col, lw=2, label=lab)
    for a in (ax, bx):
        a.set_xscale("log")
        a.set_xlim(1e-1, 1e5)
        a.grid(alpha=0.25, which="both")
        a.set_xlabel(r"$\Delta$ at the late snapshot")
        a.legend(frameon=False, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("fraction per bin")
    ax.set_title("a  converted particles vs all matter", loc="left",
                 fontsize=11)
    bx.set_ylabel(r"fraction of converted above $\Delta$")
    bx.set_title("b  cumulative", loc="left", fontsize=11)

    for i, lab in enumerate(labels):
        c = C[lab]
        cx = fig.add_subplot(gs[1 + i // 2, i % 2])
        m = c["proj_mat"]
        cx.imshow(np.log10(np.maximum(m, 1e-3)).T, origin="lower",
                  extent=[0, c["box"], 0, c["box"]], cmap="bone_r",
                  aspect="equal")
        p = c["proj_xy"]
        if len(p):
            step = max(1, len(p) // 200_000)
            cx.scatter(p[::step, 0], p[::step, 1], s=0.12, c="#c0271f",
                       alpha=0.35, linewidths=0)
        cx.set_title(f"{lab}: converted particles over the matter slab "
                     f"({100 * args.slab:.0f}% of the box)", loc="left",
                     fontsize=10)
        cx.set_xlabel("x")
        cx.set_ylabel("y")

    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
