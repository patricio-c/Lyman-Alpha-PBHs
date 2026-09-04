#!/usr/bin/env python3
"""
Stage 07 - at what density was the gas when it was still gas?

Stage 06 measures where gas is missing at one epoch.  It cannot say why a
bin is short, because it conflates two different things: gas that was
converted, and gas that moved to another density.  This stage separates
them, using the one piece of information the snapshots do carry.

THE IDEA

A gas particle that exists in `PartType0` at an early snapshot and is
absent from `PartType0` at a later one was converted in between.  QLA
conversion is one-way and particles keep their ID, so

    doomed = IDs(PartType0, early)  \\  IDs(PartType0, late)

is exact, not a statistical match.  Particle splitting only ever ADDS new
IDs at the later time, so it cannot put a particle into this set
spuriously.  Having identified those particles, we read their density and
temperature AT THE EARLY SNAPSHOT, while they were still gas.  That is the
distribution stage 06 cannot reach.

WHAT IT ANSWERS

  - was the doomed gas already dense when we last saw it (Delta ~ 100-1000,
    destined, nothing to do with the forest), or was it ordinary IGM
    (Delta ~ 1-10) that got dragged up afterwards?
  - the second case is the one that would make the threshold a genuine
    forest systematic rather than a bookkeeping detail.

WHAT IT CANNOT ANSWER

It only sees conversions that happen BETWEEN the two snapshots.  Check the
converted fraction at each epoch before trusting the answer: if most of a
run's conversion happened before the early snapshot, this stage
characterises the tail, not the bulk, and it says so at runtime.

Usage
-----
    python stages/07_doomed_gas.py \\
        --early .../murgia-cdm-lyman_0000/murgia-cdm-lyman_0000.hdf5 \\
                .../murgia-M3-lyman_0000/murgia-M3-lyman_0000.hdf5 \\
        --late  .../murgia-cdm-lyman_0002/murgia-cdm-lyman_0002.hdf5 \\
                .../murgia-M3-lyman_0002/murgia-M3-lyman_0002.hdf5 \\
        --labels cdm M3 --out figures/doomed_murgia_z10_to_z5

Options
-------
    --early ...       one snapshot per run, the epoch the densities are
                      read at
    --late ...        one snapshot per run, same order and count
    --labels ...      display names, same order
    --n-ic N          particles per side in the ICs, for the converted
                      fraction (default: inferred from gas+DM totals)
    --check-dm        also verify the doomed IDs turn up in PartType1 at
                      the late snapshot.  Exact but expensive: it streams
                      the whole DM ID array.
    --chunk N         particles read at a time (default 8e6)
    --d-lo X          low edge of the Delta histogram (default 1e-2)
    --d-hi X          high edge (default 1e6)
    --nbins N         log bins between them (default 160)
    --out PREFIX      writes PREFIX.png and PREFIX.txt

A note on picking the early snapshot.  The earliest one maximises coverage
- an IC snapshot sees 100% of the conversion - but at high redshift the
density field is nearly uniform, so the whole distribution collapses into
one or two of the default bins and the percentiles become meaningless.
Use --d-lo / --d-hi / --nbins to zoom on the range the gas occupies; the
run warns when the binning, rather than the data, is setting the answer.
Near the ICs the question also changes: it stops being "was this forest
gas" and becomes "does the material that later converts start out
overdense", which the contrast column answers.
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

RHO_CRIT0_CGS = 1.8788e-29        # g/cm^3, times h^2
D_LO, D_HI = 1e-2, 1e6
NB = 160


def cgs_factor(dset, physical=True):
    """SWIFT tags each dataset with its CGS conversion; match on substring."""
    want = "physical cgs" if physical else "cgs"
    best = None
    for k in dset.attrs:
        kl = k.lower()
        if "conversion factor" not in kl or "cgs" not in kl:
            continue
        if want in kl:
            return float(np.ravel(dset.attrs[k])[0])
        best = float(np.ravel(dset.attrs[k])[0])
    if best is None:
        raise SystemExit(f"{dset.name}: no CGS conversion attribute. "
                         f"Attributes present: {list(dset.attrs)}")
    return best


def cosmology(f):
    """Omega_b, h, z from the file. Fails loudly rather than guessing."""
    at = dict(f["Cosmology"].attrs) if "Cosmology" in f else {}
    at.update({k: v for k, v in f["Header"].attrs.items() if k not in at})

    def grab(names, what):
        for n in names:
            for k in at:
                if k.lower() == n:
                    return float(np.ravel(at[k])[0])
        raise SystemExit(
            f"could not read {what} from Cosmology/Header. Keys: {sorted(at)}")

    return (grab(("omega_b", "omega_baryon", "omega_baryons"), "Omega_b"),
            grab(("h", "hubbleparam", "hubble param"), "h"),
            grab(("redshift",), "Redshift"))


def gas_key(g, *names, required=True):
    for n in names:
        if n in g:
            return n
    if required:
        raise SystemExit(f"PartType0 has none of {names}. Present: {list(g)}")
    return None


def late_gas_ids(path, chunk):
    """Sorted array of the gas IDs that still exist at the late snapshot."""
    with h5py.File(path, "r") as f:
        d = f["PartType0"][gas_key(f["PartType0"], "ParticleIDs")]
        out = np.empty(d.shape[0], dtype=np.uint64)
        for i0 in range(0, d.shape[0], chunk):
            i1 = min(i0 + chunk, d.shape[0])
            out[i0:i1] = d[i0:i1]
    out.sort()
    return out


def doomed_census(early, late_ids, args):
    """
    Histogram, in Delta, the early-time gas split into the particles that
    survive to the late snapshot and the ones that do not.
    """
    with h5py.File(early, "r") as f:
        ob, h, z = cosmology(f)
        g = f["PartType0"]
        ikey = gas_key(g, "ParticleIDs")
        dkey = gas_key(g, "Densities", "Density")
        mkey = gas_key(g, "Masses", "Mass")

        f_rho = cgs_factor(g[dkey])
        f_m = cgs_factor(g[mkey])
        rho_b = ob * RHO_CRIT0_CGS * h * h * (1.0 + z) ** 3

        de = np.logspace(np.log10(args.d_lo), np.log10(args.d_hi),
                         args.nbins + 1)
        h_doom = np.zeros(args.nbins)
        h_surv = np.zeros(args.nbins)
        n = g[ikey].shape[0]
        n_doom = 0
        m_doom = 0.0
        m_all = 0.0
        doomed_ids = []

        for i0 in range(0, n, args.chunk):
            i1 = min(i0 + args.chunk, n)
            ids = np.asarray(g[ikey][i0:i1], dtype=np.uint64)
            # searchsorted membership: exact, and no 3 GB argsort
            pos = np.searchsorted(late_ids, ids)
            np.clip(pos, 0, len(late_ids) - 1, out=pos)
            alive = late_ids[pos] == ids
            gone = ~alive

            m = g[mkey][i0:i1].astype(np.float64) * f_m
            d = g[dkey][i0:i1].astype(np.float64) * f_rho / rho_b
            h_doom += np.histogram(d[gone], bins=de, weights=m[gone])[0]
            h_surv += np.histogram(d[alive], bins=de, weights=m[alive])[0]
            n_doom += int(gone.sum())
            m_doom += float(m[gone].sum())
            m_all += float(m.sum())
            if args.check_dm:
                doomed_ids.append(ids[gone].copy())
            del ids, pos, alive, gone, m, d

    return dict(z=z, n_gas=n, n_doom=n_doom, m_doom=m_doom, m_all=m_all,
                d_edges=de, h_doom=h_doom, h_surv=h_surv,
                doomed_ids=(np.concatenate(doomed_ids) if doomed_ids
                            else None))


def dm_contains(path, ids, chunk):
    """How many of `ids` are in PartType1 at `path`. Streams, never loads."""
    ids = np.sort(ids)
    found = np.zeros(len(ids), dtype=bool)
    with h5py.File(path, "r") as f:
        d = f["PartType1"]["ParticleIDs"]
        for i0 in range(0, d.shape[0], chunk):
            i1 = min(i0 + chunk, d.shape[0])
            c = np.asarray(d[i0:i1], dtype=np.uint64)
            pos = np.searchsorted(ids, c)
            np.clip(pos, 0, len(ids) - 1, out=pos)
            hit = ids[pos] == c
            found[pos[hit]] = True
    return int(found.sum())


def infer_ic_gas(total, given=None):
    """
    Recover the IC gas count from the conserved gas + DM total.

    gas + DM is invariant under QLA conversion, and the ICs lay the gas on a
    cubic grid with an integer number of DM particles per gas particle.  So
    total = n_gas_ic * (1 + ratio), and only one candidate makes n_gas_ic a
    perfect cube.  Refuses to guess if none does - getting this wrong would
    silently rescale every converted fraction printed below.
    """
    if given:
        return int(given), None
    for ratio in (1, 2, 4, 8, 16):
        if total % (1 + ratio):
            continue
        n = total // (1 + ratio)
        side = round(n ** (1.0 / 3.0))
        if side ** 3 == n:
            return int(n), ratio
    raise SystemExit(
        f"cannot infer the IC gas count: gas + DM = {total:,} is not "
        f"n^3 * (1 + r) for any cube n and r in (1, 2, 4, 8, 16). Pass "
        f"--n-ic explicitly rather than letting this stage guess, because "
        f"every converted fraction it prints is scaled by that number.")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--early", nargs="+", required=True)
    ap.add_argument("--late", nargs="+", required=True)
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--n-ic", type=int, default=None)
    ap.add_argument("--check-dm", action="store_true")
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--d-lo", type=float, default=D_LO)
    ap.add_argument("--d-hi", type=float, default=D_HI)
    ap.add_argument("--nbins", type=int, default=NB)
    ap.add_argument("--out", default="figures/doomed_gas")
    args = ap.parse_args()

    if len(args.early) != len(args.late):
        raise SystemExit("--early and --late must have the same count")
    labels = args.labels or [f"run{i}" for i in range(len(args.early))]
    if len(labels) != len(args.early):
        raise SystemExit("--labels must have one entry per run")

    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    say("=" * 74)
    say("Stage 07 - the density of the gas that was about to be converted")
    say("=" * 74)
    say("")
    say("Delta is read at the EARLY snapshot, while the gas was still gas.")
    say("'doomed' = present in PartType0 early, absent from PartType0 late.")

    C = {}
    for lab, e, l in zip(labels, args.early, args.late):
        for p in (e, l):
            if not os.path.isfile(p):
                raise SystemExit(f"{p} is not a file")
        with h5py.File(e, "r") as f:
            ne = np.atleast_1d(f["Header"].attrs["NumPart_Total"])
            ze = float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])
        with h5py.File(l, "r") as f:
            nl = np.atleast_1d(f["Header"].attrs["NumPart_Total"])
            zl = float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])
        if zl >= ze:
            raise SystemExit(
                f"{lab}: --late is at z={zl} and --early at z={ze}. The late "
                f"snapshot has to be the later one or the set difference is "
                f"meaningless.")

        if int(ne[0]) + int(ne[1]) != int(nl[0]) + int(nl[1]):
            raise SystemExit(
                f"{lab}: gas + DM is {int(ne[0]) + int(ne[1]):,} early and "
                f"{int(nl[0]) + int(nl[1]):,} late. Those must be equal - QLA "
                f"conversion moves particles between the two types, it does "
                f"not create or destroy them. Something else is going on and "
                f"the set difference below would not mean what it claims.")
        n_ic, ratio = infer_ic_gas(int(ne[0]) + int(ne[1]), args.n_ic)
        say(f"\n[{lab}]  z = {ze:.4f} -> {zl:.4f}")
        if args.n_ic:
            say(f"       IC gas particles (given): {n_ic:,}")
        else:
            say(f"       IC gas particles: {n_ic:,} = {round(n_ic ** (1/3))}^3,"
                f" with {ratio} DM per gas particle")
        f_e = 1.0 - int(ne[0]) / n_ic
        f_l = 1.0 - int(nl[0]) / n_ic
        say(f"       converted by the early snapshot: {100 * f_e:6.2f}%")
        say(f"       converted by the late  snapshot: {100 * f_l:6.2f}%")
        span = f_l - f_e
        say(f"       this stage therefore sees {100 * span:6.2f} points of"
            f" conversion,")
        say(f"       i.e. {100 * span / f_l if f_l > 0 else 0:.0f}% of all the"
            f" conversion in this run.")
        if f_l > 0 and span / f_l < 0.5:
            say("       [!] most of this run's conversion happened BEFORE the")
            say("           early snapshot. What follows describes the tail,")
            say("           not the bulk. Pick an earlier snapshot if one")
            say("           exists.")

        say("       reading late gas IDs ...")
        li = late_gas_ids(l, args.chunk)
        say("       scanning early snapshot ...")
        C[lab] = doomed_census(e, li, args)
        c = C[lab]
        del li
        say(f"       doomed particles {c['n_doom']:,} of {c['n_gas']:,}"
            f"  ({100 * c['n_doom'] / c['n_gas']:.2f}%)")
        say(f"       doomed gas MASS  {100 * c['m_doom'] / c['m_all']:.2f}%"
            f" of the early gas mass")
        if args.check_dm and c["doomed_ids"] is not None:
            k = dm_contains(l, c["doomed_ids"], args.chunk)
            say(f"       cross-check: {k:,} of {c['n_doom']:,} doomed IDs are"
                f" in PartType1 at the late snapshot"
                f"  ({100 * k / max(c['n_doom'], 1):.2f}%)")
            if k < 0.99 * c["n_doom"]:
                say("       [!] they did not all become PartType1. Something")
                say("           other than QLA conversion removed gas here.")
            c["doomed_ids"] = None

    say("\n" + "-" * 74)
    say("WHERE THE DOOMED GAS WAS, at the early snapshot")
    say("-" * 74)
    say(f"{'run':8s} {'doomed med':>11s} {'surv med':>10s} {'contrast':>9s}"
        f" {'90th pct':>10s} {'% >100':>8s} {'% >1000':>9s} {'% <10':>8s}")
    coarse = []
    for lab in labels:
        c = C[lab]
        dc = np.sqrt(c["d_edges"][:-1] * c["d_edges"][1:])
        w = c["h_doom"]
        tot = w.sum()
        if tot <= 0:
            say(f"{lab:8s}   no doomed gas inside the Delta range")
            continue
        cw = np.cumsum(w) / tot
        med = float(np.interp(0.5, cw, dc))
        p90 = float(np.interp(0.9, cw, dc))
        sw = c["h_surv"]
        if sw.sum() > 0:
            smed = float(np.interp(0.5, np.cumsum(sw) / sw.sum(), dc))
            contrast = f"{med / smed:9.4f}"
        else:
            smed, contrast = float("nan"), "        -"
        above = lambda X: 100 * w[dc >= X].sum() / tot   # noqa: E731
        say(f"{lab:8s} {med:11.4f} {smed:10.4f} {contrast}"
            f" {p90:10.4f} {above(100):7.2f}% {above(1000):8.2f}%"
            f" {100 - above(10):7.2f}%")
        # how many bins actually hold the middle 80% of the doomed mass?
        lo = float(np.interp(0.1, cw, np.arange(len(dc))))
        hi = float(np.interp(0.9, cw, np.arange(len(dc))))
        if hi - lo < 4.0:
            coarse.append(lab)

    say("")
    say("Read it like this: if the doomed gas sat at Delta of order 100-1000")
    say("it was already destined and never was forest gas.  If it sat at")
    say("Delta of a few, ordinary IGM is being removed and the threshold is")
    say("a real forest systematic.")
    say("")
    say("'contrast' is the doomed median over the surviving median at the")
    say("same epoch.  Near the initial conditions that is the only column")
    say("that means anything: Delta is ~1 everywhere and the question is")
    say("whether the material that later converts starts out measurably")
    say("overdense.  Contrast > 1 says yes.")
    if coarse:
        say("")
        say(f"[!] for {', '.join(coarse)} the middle 80% of the doomed mass")
        say("    spans fewer than 4 bins, so these percentiles are set by the")
        say("    binning, not by the data.  Re-run over the range the gas")
        say("    actually occupies, e.g. --d-lo 0.5 --d-hi 2 --nbins 200.")

    # ---------------------------------------------------------------- figure
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    PAL = ["#1f4fbf", "#c0271f", "#2ca02c", "#9467bd", "#ff7f0e"]
    for i, lab in enumerate(labels):
        c = C[lab]
        dc = np.sqrt(c["d_edges"][:-1] * c["d_edges"][1:])
        col = PAL[i % len(PAL)]
        if c["h_doom"].sum() > 0:
            ax.plot(dc, c["h_doom"] / c["h_doom"].sum(), color=col, lw=2,
                    label=f"{lab}, doomed")
        if c["h_surv"].sum() > 0:
            ax.plot(dc, c["h_surv"] / c["h_surv"].sum(), color=col, lw=1.3,
                    ls=":", label=f"{lab}, survives")
        tot = c["h_doom"].sum()
        if tot > 0:
            bx.plot(dc, np.cumsum(c["h_doom"][::-1])[::-1] / tot,
                    color=col, lw=2, label=lab)
    for a in (ax, bx):
        a.axvline(1000.0, color="0.35", ls="--", lw=1.3)
        a.set_xscale("log")
        a.set_xlim(args.d_lo, args.d_hi)
        a.grid(alpha=0.25, which="both")
        a.set_xlabel(r"$\Delta$ at the early snapshot")
        a.legend(frameon=False, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("fraction of the mass per bin")
    ax.set_title("a  density of doomed vs surviving gas", loc="left",
                 fontsize=11)
    bx.set_ylabel("fraction of doomed mass above $\\Delta$")
    bx.set_title("b  cumulative, doomed gas only", loc="left", fontsize=11)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"\nwritten -> {args.out}.png")
    print(f"written -> {args.out}.txt")


if __name__ == "__main__":
    main()
