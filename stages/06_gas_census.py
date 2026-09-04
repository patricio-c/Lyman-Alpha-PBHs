#!/usr/bin/env python3
"""
Stage 06 - where, in the density-temperature plane, is the gas missing?

Joop's question, answered without a halo catalogue.  The star-formation
threshold removes more gas in a model with extra small-scale power; this
stage says *where* it goes, which is what decides whether the threshold
value is load-bearing.

The fork it resolves, and why it is worth doing before re-running SWIFT
with a different threshold:

  - if the extra removal sits just above `Delta = 1000`, the result is
    hypersensitive to a numerical choice and raising the threshold will
    sweep most of the differential away;
  - if it sits far above, the threshold value barely matters, the
    differential is robust, and the re-runs will move nothing.

Those predict opposite outcomes for the threshold sweep and this script
separates them on snapshots that are already on disk.

WHAT IT MEASURES, AND WHAT IT CANNOT

Converted gas is gone from `PartType0`, so nothing here histograms it
directly.  What is measured is the *surviving* gas of each run in the
(Delta, T) plane, and the difference between runs.  Because all runs share
the same Panphasia phases, that difference is not a statistical statement:
cell by cell in the plane it is the gas one run has and the other does not.
It still conflates "converted" with "moved", so read it as missing gas,
not as converted gas.  To follow the converted particles themselves you
need `PartType1` split by mass - run `stages/00_inspect_snapshot.py --deep`
first to find out whether that is possible in your run.

Delta is normalised by the COSMIC mean baryon density from `Omega_b`, not
by the surviving mean gas density of each run.  That matters: a run that
lost 15% of its gas would otherwise get its Delta axis shifted by 15% and
the two runs would not be on the same scale.  The script fails loudly if it
cannot read `Omega_b` rather than guessing.

Usage
-----
    python stages/06_gas_census.py \\
        --runs /data/.../lyman/murgia/cdm /data/.../lyman/murgia/M3 \\
        --labels cdm M3 --z 5.0 --out figures/gas_census_murgia_z5

    # the FCT pair, where the lever is four times larger
    python stages/06_gas_census.py --runs cdm40 fct40 --labels CDM FCT \\
        --z 3.0 --out figures/gas_census_fct_z3

Options
-------
    --runs ...        two or more registry names or directories; the first
                      is the reference every ratio is taken against
    --snaps ...       explicit snapshot paths instead, same order and count.
                      Needed for runs whose snapshots are not named *snap* or
                      *output* - murgia is one (murgia-cdm-lyman_0002.hdf5),
                      and common/runs.py filters on those substrings, so the
                      resolver finds nothing. Locate them with
                        find <run dir> -name '*lyman_*.hdf5' -size +0 | sort
    --labels ...      display names, same order (default: basenames)
    --z Z             pick the snapshot closest to this redshift
    --threshold X     overdensity to mark (default 1000, the QLA value)
    --chunk N         particles read at a time (default 8e6, ~100 MB)
    --nbins-d N       log-Delta bins (default 160, over 1e-2 to 1e6)
    --nbins-t N       log-T bins (default 120, over 1e2 to 1e8)
    --out PREFIX      writes PREFIX.png and PREFIX.txt
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import runs  # noqa: E402

RHO_CRIT0_CGS = 1.8788e-29        # g/cm^3, times h^2
MPC_CM = 3.0856775814913673e24

D_LO, D_HI = 1e-2, 1e6
T_LO, T_HI = 1e2, 1e8


def cgs_factor(dset, physical=True):
    """
    SWIFT tags every dataset with its conversion to CGS.  Match on a
    substring rather than the exact attribute name, which has changed
    between SWIFT versions.
    """
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
    """Omega_b, h and z, read from the file. Fails loudly, never guesses."""
    at = dict(f["Cosmology"].attrs) if "Cosmology" in f else {}
    at.update({k: v for k, v in f["Header"].attrs.items() if k not in at})

    def grab(names, what):
        for n in names:
            for k in at:
                if k.lower() == n:
                    return float(np.ravel(at[k])[0])
        raise SystemExit(
            f"could not read {what} from Cosmology/Header. Keys present: "
            f"{sorted(at)}\nWithout it Delta cannot be put on a common scale "
            f"between runs, so this stage refuses to continue.")

    ob = grab(("omega_b", "omega_baryon", "omega_baryons"), "Omega_b")
    h = grab(("h", "hubbleparam", "hubble param"), "h")
    z = grab(("redshift",), "Redshift")
    return ob, h, z


def census(path, args):
    """One-pass chunked histogram of the gas in (Delta, T). Returns dict."""
    with h5py.File(path, "r") as f:
        ob, h, z = cosmology(f)
        g = f["PartType0"]

        def field(*names, required=True):
            for n in names:
                if n in g:
                    return n
            if required:
                raise SystemExit(f"{path}: PartType0 has none of {names}. "
                                 f"Datasets present: {list(g)}")
            return None

        dkey = field("Densities", "Density")
        mkey = field("Masses", "Mass")
        tkey = field("Temperatures", "Temperature", required=False)

        f_rho = cgs_factor(g[dkey])
        f_m = cgs_factor(g[mkey])
        f_T = cgs_factor(g[tkey]) if tkey else 1.0

        # cosmic mean baryon density, physical, at this redshift
        rho_b = ob * RHO_CRIT0_CGS * h * h * (1.0 + z) ** 3

        n = g[mkey].shape[0]
        de = np.logspace(np.log10(D_LO), np.log10(D_HI), args.nbins_d + 1)
        te = np.logspace(np.log10(T_LO), np.log10(T_HI), args.nbins_t + 1)
        h1 = np.zeros(args.nbins_d)
        h2 = np.zeros((args.nbins_d, args.nbins_t))
        mtot = 0.0

        for i0 in range(0, n, args.chunk):
            i1 = min(i0 + args.chunk, n)
            m = g[mkey][i0:i1].astype(np.float64) * f_m
            d = g[dkey][i0:i1].astype(np.float64) * f_rho / rho_b
            mtot += m.sum()
            h1 += np.histogram(d, bins=de, weights=m)[0]
            if tkey:
                T = g[tkey][i0:i1].astype(np.float64) * f_T
                h2 += np.histogram2d(d, T, bins=[de, te], weights=m)[0]
            del m, d

    return dict(path=path, z=z, omega_b=ob, h=h, rho_b=rho_b, n_gas=n,
                m_gas_total=mtot, d_edges=de, t_edges=te, hist_d=h1,
                hist_dt=h2, has_T=bool(tkey))


def qla_params(run_directory):
    """Whatever SWIFT actually applied, if used_parameters.yml is there."""
    out = []
    for p in glob.glob(os.path.join(run_directory, "**", "used_parameters.yml"),
                       recursive=True)[:1]:
        keep = False
        for line in open(p, errors="replace"):
            s = line.rstrip()
            if s and not s.startswith((" ", "\t", "#")):
                keep = any(t in s.lower() for t in
                           ("qla", "starformation", "entropyfloor", "cooling"))
            if keep and s.strip():
                out.append(s)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--snaps", nargs="*", default=None)
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--z", type=float, required=True)
    ap.add_argument("--threshold", type=float, default=1000.0)
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--nbins-d", type=int, default=160)
    ap.add_argument("--nbins-t", type=int, default=120)
    ap.add_argument("--out", default="figures/gas_census")
    args = ap.parse_args()

    labels = args.labels or [os.path.basename(os.path.normpath(r))
                             for r in args.runs]
    if len(labels) != len(args.runs):
        raise SystemExit("--labels must have one entry per run")
    if args.snaps and len(args.snaps) != len(args.runs):
        raise SystemExit("--snaps must have one entry per run")

    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    say("=" * 74)
    say(f"Stage 06 - gas census in the (Delta, T) plane    z = {args.z}")
    say("=" * 74)
    say("")
    say("Delta here is rho_gas / rho_baryon(cosmic mean), from Omega_b.")
    say(f"The line marked at Delta = {args.threshold:.6g} assumes SWIFT's")
    say("QLAStarFormation:over_density is defined against the mean BARYON")
    say("density. CONFIRM THAT in the SWIFT source or documentation before")
    say("reading the 'load-bearing' verdict below: if it is defined against")
    say("the mean TOTAL MATTER density instead, the line moves by")
    say("Omega_m/Omega_b ~ 6.7 and the marginal-versus-deep answer changes.")

    C = {}
    for i, (r, lab) in enumerate(zip(args.runs, labels)):
        if args.snaps:
            snap = os.path.abspath(args.snaps[i])
            if not os.path.isfile(snap):
                raise SystemExit(f"--snaps: {snap} is not a file")
        else:
            snap = runs.resolve_snapshot(r, z=args.z)
        say(f"\n[{lab}] {snap}")
        C[lab] = census(snap, args)
        c = C[lab]
        say(f"       z = {c['z']:.4f}   Omega_b = {c['omega_b']:.5f}   "
            f"h = {c['h']:.5f}")
        if abs(c["z"] - args.z) > 0.05:
            raise SystemExit(
                f"{lab}: this snapshot is at z={c['z']:.4f} but you asked for "
                f"z={args.z}. Analysing the wrong epoch silently is exactly "
                f"the class of bug this repository exists to catch.")
        say(f"       gas particles {c['n_gas']:,}   "
            f"mean rho_b = {c['rho_b']:.4e} g/cm^3")
        pars = qla_params(runs.run_dir(r))
        if pars:
            say("       used_parameters.yml, star formation block:")
            for line in pars[:24]:
                say(f"         {line}")
        else:
            say("       [!] no used_parameters.yml found - the QLA criterion "
                "is unconfirmed")

    ref = labels[0]
    de = C[ref]["d_edges"]
    dc = np.sqrt(de[:-1] * de[1:])

    # everything is quoted as a fraction of the REFERENCE run's total gas, so
    # the runs share a denominator and the numbers add up to the loss.
    m_ref = C[ref]["m_gas_total"]

    say("\n" + "-" * 74)
    say("SURVIVING GAS MASS, as a fraction of the reference run's total")
    say("-" * 74)
    say(f"{'run':10s} {'M_gas / M_gas(ref)':>20s} {'missing vs ref':>16s}")
    for lab in labels:
        r = C[lab]["m_gas_total"] / m_ref
        say(f"{lab:10s} {r:20.5f} {100 * (1 - r):15.2f}%")

    say("\n" + "-" * 74)
    say(f"CUMULATIVE GAS MASS ABOVE Delta   (fraction of {ref} total)")
    say("-" * 74)
    marks = [1.0, 10.0, 100.0, args.threshold, 3 * args.threshold,
             1e4, 1e5]
    say(f"{'Delta >':>10s} " + " ".join(f"{l:>12s}" for l in labels)
        + f" {'missing':>12s}")
    for D in marks:
        vals = []
        for lab in labels:
            cum = np.cumsum(C[lab]["hist_d"][::-1])[::-1] / m_ref
            vals.append(float(np.interp(np.log10(D), np.log10(dc), cum)))
        miss = vals[0] - vals[-1]
        tag = "  <- threshold" if abs(D - args.threshold) < 1e-9 else ""
        say(f"{D:10.4g} " + " ".join(f"{v:12.5f}" for v in vals)
            + f" {miss:12.5f}{tag}")

    say("\n" + "-" * 74)
    say("WHERE THE MISSING GAS SITS  (reference minus last run, per decade)")
    say("-" * 74)
    dm = C[ref]["hist_d"] - C[labels[-1]]["hist_d"]
    tot_missing = dm[dm > 0].sum()
    say(f"total missing mass = {tot_missing / m_ref:.5f} of {ref}")
    if tot_missing <= 0:
        say("[!] the last run has MORE gas than the reference at every "
            "Delta; nothing to attribute.")
    else:
        say(f"{'decade in Delta':>22s} {'share of the missing mass':>26s}")
        edges = [1e-2, 1e0, 1e1, 1e2, args.threshold, 3 * args.threshold,
                 1e4, 1e5, 1e6]
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (dc >= a) & (dc < b)
            share = dm[sel][dm[sel] > 0].sum() / tot_missing
            say(f"{a:10.4g} - {b:<9.4g} {100 * share:25.1f}%")
        below = dm[(dc < args.threshold)][dm[(dc < args.threshold)] > 0].sum()
        near = dm[(dc >= args.threshold) & (dc < 3 * args.threshold)]
        say("")
        say(f"below the threshold          {100 * below / tot_missing:6.1f}%"
            "   (moved, not converted - the runs differ dynamically too)")
        say(f"within a factor 3 above it   "
            f"{100 * near[near > 0].sum() / tot_missing:6.1f}%"
            "   <- if this is large, the threshold VALUE is load-bearing")

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(12.5, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.24)
    ax = fig.add_subplot(gs[0, 0])
    bx = fig.add_subplot(gs[0, 1])
    cx = fig.add_subplot(gs[1, :])
    PAL = ["#1f4fbf", "#c0271f", "#2ca02c", "#9467bd", "#ff7f0e"]

    for i, lab in enumerate(labels):
        cum = np.cumsum(C[lab]["hist_d"][::-1])[::-1] / m_ref
        ax.plot(dc, cum, color=PAL[i % len(PAL)], lw=2, label=lab)
    ax.axvline(args.threshold, color="0.35", ls="--", lw=1.3)
    ax.text(args.threshold * 1.15, 0.55, r"QLA $\Delta$ threshold",
            rotation=90, fontsize=9, color="0.35", va="center")
    ax.set(xscale="log", yscale="log", xlim=(1e-1, D_HI),
           xlabel=r"$\Delta = \rho_{\rm gas}/\bar\rho_{\rm b}$",
           ylabel=f"gas mass above $\\Delta$  /  total in {ref}")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25, which="both")
    ax.set_title("a  cumulative gas mass", loc="left", fontsize=11)

    for i, lab in enumerate(labels[1:], start=1):
        r = np.divide(C[lab]["hist_d"], C[ref]["hist_d"],
                      out=np.full_like(dc, np.nan),
                      where=C[ref]["hist_d"] > 0)
        bx.plot(dc, r, color=PAL[i % len(PAL)], lw=2, label=f"{lab} / {ref}")
    bx.axhline(1.0, color="k", ls=":", lw=1)
    bx.axvline(args.threshold, color="0.35", ls="--", lw=1.3)
    bx.set(xscale="log", xlim=(1e-1, D_HI), ylim=(0, 1.6),
           xlabel=r"$\Delta$", ylabel="gas mass per bin, ratio")
    bx.legend(frameon=False, fontsize=9)
    bx.grid(alpha=0.25, which="both")
    bx.set_title("b  where the divergence starts", loc="left", fontsize=11)

    if C[ref]["has_T"] and len(labels) > 1:
        te = C[ref]["t_edges"]
        d2 = C[ref]["hist_dt"] - C[labels[-1]]["hist_dt"]
        v = np.nanmax(np.abs(d2)) or 1.0
        im = cx.pcolormesh(de, te, (d2 / v).T, cmap="RdBu_r",
                           vmin=-1, vmax=1, shading="auto")
        cx.axvline(args.threshold, color="k", ls="--", lw=1.3)
        cx.set(xscale="log", yscale="log", xlim=(1e-1, D_HI), ylim=(1e2, 1e8),
               xlabel=r"$\Delta$", ylabel="T  [K]")
        cx.set_title(f"c  gas in {ref} minus gas in {labels[-1]}, by mass "
                     "(red = missing from the model)", loc="left", fontsize=11)
        fig.colorbar(im, ax=cx, label="mass difference / peak")
    else:
        cx.text(0.5, 0.5, "no Temperatures dataset - panel c skipped",
                ha="center", va="center", transform=cx.transAxes)
        cx.axis("off")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=150, bbox_inches="tight")
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"\nwritten -> {args.out}.png")
    print(f"written -> {args.out}.txt")


if __name__ == "__main__":
    main()
