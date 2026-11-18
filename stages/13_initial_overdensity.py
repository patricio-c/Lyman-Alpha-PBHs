#!/usr/bin/env python3
"""
Stage 13 - initial overdensity of the gas that QLA later converts.

Exact mass-weighted medians of Delta at the early snapshot (no histogram, so
no range clipping), for converted and surviving gas in CDM and FCT, and a
figure of the four distributions. Delta, the ID matching and the unit handling
are imported from stage 07, so the numbers are defined exactly as there.

    python stages/13_initial_overdensity.py \\
        --cdm /data/contrib/pad_140/pcolazo/lyman/Old_sim/cdm-box-40-1024 \\
        --fct /data/contrib/pad_140/pcolazo/lyman/Old_sim/2-fct-box-40-1024 \\
        --out figures/qla_initial_overdensity

Writes PREFIX.png (300 dpi), PREFIX.pdf, PREFIX.json and PREFIX.txt.
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
import re
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from common import prov  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "stage07", os.path.join(ROOT, "stages", "07_doomed_gas.py"))
s07 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s07)

COLOR = {"CDM": "#1f4fbf", "FCT": "#c0271f"}
NAME = {"CDM": r"$\Lambda$CDM", "FCT": "FCT"}


def snapshot_at(run_dir, z, tol):
    """
    Files of the snapshot closest to z: the virtual file if SWIFT wrote one,
    otherwise its .N pieces (reading both would count every particle twice).
    """
    found = (glob.glob(os.path.join(run_dir, "*.hdf5"))
             + glob.glob(os.path.join(run_dir, "*", "*.hdf5")))
    snaps = {}
    for p in found:
        if os.path.basename(p).startswith("los_"):
            continue
        stem = re.sub(r"\.\d+\.hdf5$", "", p).removesuffix(".hdf5")
        snaps.setdefault(stem, []).append(p)
    zs = []
    for stem, files in snaps.items():
        virtual = [p for p in files if p == stem + ".hdf5"]
        files = virtual or sorted(
            files, key=lambda p: int(re.search(r"\.(\d+)\.hdf5$", p).group(1)))
        try:
            with h5py.File(files[0], "r") as f:
                if "PartType0" not in f:
                    continue
                zs.append((files, float(np.atleast_1d(
                    f["Header"].attrs["Redshift"])[0])))
        except Exception:
            continue
    if not zs:
        raise SystemExit(f"no snapshots with PartType0 under {run_dir}")
    files, zb = min(zs, key=lambda t: abs(t[1] - z))
    if abs(zb - z) > tol:
        raise SystemExit(
            f"{run_dir}: closest snapshot to z={z} is z={zb:.4f} "
            f"({files[0]}). Available: {sorted(round(t[1], 3) for t in zs)}")
    return files


def numpart(files):
    if len(files) == 1:
        with h5py.File(files[0], "r") as f:
            return np.atleast_1d(f["Header"].attrs["NumPart_Total"]).astype(
                np.int64)
    tot = 0
    for p in files:
        with h5py.File(p, "r") as f:
            tot = tot + np.atleast_1d(
                f["Header"].attrs["NumPart_ThisFile"]).astype(np.int64)
    return tot


def gas_ids(files, chunk):
    parts = []
    for p in files:
        with h5py.File(p, "r") as f:
            d = f["PartType0"][s07.gas_key(f["PartType0"], "ParticleIDs")]
            for i0 in range(0, d.shape[0], chunk):
                parts.append(np.asarray(d[i0:i0 + chunk], dtype=np.uint64))
    out = np.concatenate(parts)
    out.sort()
    return out


def census(early, late, chunk):
    """Delta and mass of every early-time gas particle, and whether it survives."""
    late_ids = gas_ids(late, chunk)
    n = int(numpart(early)[0])
    delta = np.empty(n, dtype=np.float32)
    # internal units: in grams a gas particle (~1e40 g) overflows float32
    mass = np.empty(n, dtype=np.float32)
    alive = np.empty(n, dtype=bool)
    j = 0
    for p in early:
        with h5py.File(p, "r") as f:
            ob, h, z = s07.cosmology(f)
            g = f["PartType0"]
            ikey = s07.gas_key(g, "ParticleIDs")
            dkey = s07.gas_key(g, "Densities", "Density")
            mkey = s07.gas_key(g, "Masses", "Mass")
            f_rho = s07.cgs_factor(g[dkey])
            rho_b = ob * s07.RHO_CRIT0_CGS * h * h * (1.0 + z) ** 3
            nf = g[ikey].shape[0]
            for i0 in range(0, nf, chunk):
                i1 = min(i0 + chunk, nf)
                k0, k1 = j + i0, j + i1
                ids = np.asarray(g[ikey][i0:i1], dtype=np.uint64)
                pos = np.searchsorted(late_ids, ids)
                np.clip(pos, 0, len(late_ids) - 1, out=pos)
                alive[k0:k1] = late_ids[pos] == ids
                delta[k0:k1] = (g[dkey][i0:i1].astype(np.float64)
                                * f_rho / rho_b)
                mass[k0:k1] = g[mkey][i0:i1]
            j += nf
    if j != n:
        raise SystemExit(f"read {j:,} gas particles, header says {n:,}")
    return dict(z=z, delta=delta, mass=mass, alive=alive, n=n)


def wmedian(x, w):
    """Mass-weighted median: sort, accumulate the mass, cut at half."""
    idx = np.argsort(x, kind="stable")
    cw = np.cumsum(w[idx], dtype=np.float64)
    k = int(np.searchsorted(cw, 0.5 * cw[-1]))
    return float(x[idx[k]])


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cdm", required=True,
                    help="CDM run directory (or --cdm-early/--cdm-late)")
    ap.add_argument("--fct", required=True, help="FCT run directory")
    ap.add_argument("--cdm-early", nargs="+", default=None,
                    help="explicit snapshot file(s), overrides --cdm")
    ap.add_argument("--cdm-late", nargs="+", default=None)
    ap.add_argument("--fct-early", nargs="+", default=None)
    ap.add_argument("--fct-late", nargs="+", default=None)
    ap.add_argument("--z-early", type=float, default=198.0)
    ap.add_argument("--z-late", type=float, default=3.0)
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--nbins", type=int, default=400)
    ap.add_argument("--ymin", type=float, default=1e-9)
    ap.add_argument("--fontsize", type=float, default=11.0)
    ap.add_argument("--out", default="figures/qla_initial_overdensity")
    args = ap.parse_args()

    snaps = {
        "CDM": (args.cdm_early or snapshot_at(args.cdm, args.z_early, 2.0),
                args.cdm_late or snapshot_at(args.cdm, args.z_late, 0.05)),
        "FCT": (args.fct_early or snapshot_at(args.fct, args.z_early, 2.0),
                args.fct_late or snapshot_at(args.fct, args.z_late, 0.05)),
    }
    log = list(prov.header(
        inputs=[p for e_l in snaps.values() for fl in e_l for p in fl]))

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    say("=" * 74)
    say("Stage 13 - initial overdensity of converted vs surviving gas")
    say("=" * 74)
    say("Delta = rho_SPH / mean cosmic baryon density, read at the early")
    say("snapshot. Medians are mass-weighted and exact over the particles.")

    R, C = {}, {}
    for lab in ("CDM", "FCT"):
        early, late = snaps[lab]
        ne, nl = numpart(early), numpart(late)
        with h5py.File(late[0], "r") as f:
            zl = float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])
        if int(ne[0]) + int(ne[1]) != int(nl[0]) + int(nl[1]):
            raise SystemExit(f"{lab}: gas + DM not conserved between "
                             f"snapshots; the ID difference would not mean "
                             f"conversion.")
        n_ic, _ = s07.infer_ic_gas(int(ne[0]) + int(ne[1]))
        f_hdr_early = 1.0 - int(ne[0]) / n_ic
        f_hdr_late = 1.0 - int(nl[0]) / n_ic

        say(f"\n[{lab}]")
        say(f"  early  {early}")
        say(f"  late   {late}   z = {zl:.4f}")
        c = census(early, late, args.chunk)
        C[lab] = c
        conv = ~c["alive"]
        n_conv = int(conv.sum())
        m = c["mass"].astype(np.float64)
        say(f"  z_early = {c['z']:.4f},  IC gas particles {n_ic:,}")
        say(f"  converted by z_early (header): {100 * f_hdr_early:.2f}%")
        say(f"  converted by z_late  (header): {100 * f_hdr_late:.2f}%")
        say(f"  converted by ID:               {100 * n_conv / n_ic:.2f}%"
            f"   ({n_conv:,} particles, {100 * m[conv].sum() / m.sum():.2f}%"
            f" of the gas mass)")
        n_bad = int((c["mass"] <= 0).sum())
        if n_bad:
            say(f"  [!] {n_bad:,} gas particles with mass <= 0")

        med_c = wmedian(c["delta"][conv], c["mass"][conv])
        med_s = wmedian(c["delta"][~conv], c["mass"][~conv])
        say(f"  median Delta  converted {med_c:.5f}   surviving {med_s:.5f}"
            f"   separation {med_c - med_s:.5f}")
        R[lab] = dict(
            early_snapshot=early, late_snapshot=late,
            z_early=c["z"], z_late=zl, n_ic_gas=n_ic,
            n_converted=n_conv,
            converted_fraction=n_conv / n_ic,
            converted_mass_fraction=float(m[conv].sum() / m.sum()),
            converted_fraction_header=f_hdr_late,
            median_delta_converted=med_c,
            median_delta_surviving=med_s,
            separation=med_c - med_s,
            nonpositive_masses=n_bad)
        del m, conv

    ratio = R["FCT"]["separation"] / R["CDM"]["separation"]
    say("\n" + "-" * 74)
    say(f"separation FCT / CDM = {R['FCT']['separation']:.5f} / "
        f"{R['CDM']['separation']:.5f} = {ratio:.3f}")

    lo = min(float(C[k]["delta"].min()) for k in C) * 0.999
    hi = max(float(C[k]["delta"].max()) for k in C) * 1.001
    edges = np.logspace(np.log10(lo), np.log10(hi), args.nbins + 1)
    dc = np.sqrt(edges[:-1] * edges[1:])
    H = {}
    for lab in ("CDM", "FCT"):
        c = C[lab]
        conv = ~c["alive"]
        for pop, sel in (("conv", conv), ("surv", ~conv)):
            h = np.histogram(c["delta"][sel], bins=edges,
                             weights=c["mass"][sel].astype(np.float64))[0]
            tot = float(c["mass"][sel].astype(np.float64).sum())
            if abs(h.sum() / tot - 1.0) > 1e-9:
                raise SystemExit(f"{lab} {pop}: histogram lost "
                                 f"{1 - h.sum() / tot:.2e} of the mass")
            H[lab, pop] = h / h.sum()
    say(f"histogram: {args.nbins} log bins over Delta {lo:.4f} - {hi:.4f}, "
        f"the full data range (no mass outside)")

    plt.rcParams.update({"font.size": args.fontsize})
    fig, ax = plt.subplots(figsize=(6.93, 5.35))
    for lab in ("CDM", "FCT"):
        for pop, ls, lw in (("conv", "-", 2.0), ("surv", ":", 1.6)):
            y = np.where(H[lab, pop] > 0, H[lab, pop], np.nan)
            txt = (f"{NAME[lab]}, converted "
                   f"({100 * R[lab]['converted_fraction']:.1f}%)"
                   if pop == "conv" else f"{NAME[lab]}, survives")
            ax.plot(dc, y, color=COLOR[lab], ls=ls, lw=lw, label=txt)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(min(lo, 0.48), max(hi, 2.1))
    ax.set_ylim(bottom=args.ymin)
    ax.set_xticks([0.5, 1, 2])
    ax.set_xticklabels(["0.5", "1", "2"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$\Delta = \rho_{\rm SPH}/\bar\rho_{\rm b}$ at "
                  rf"$z = {R['CDM']['z_early']:.0f}$")
    ax.set_ylabel("fraction of the mass per bin")
    ax.legend(frameon=False, fontsize=args.fontsize - 2, loc="upper right")
    ax.text(0.03, 0.97,
            r"median $\Delta$, converted $-$ surviving:" "\n"
            rf"$\Lambda$CDM {R['CDM']['separation']:.3f}, "
            rf"FCT {R['FCT']['separation']:.3f} ($\times${ratio:.1f})",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=args.fontsize - 2,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7"))
    fig.tight_layout()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    fig.savefig(args.out + ".png", dpi=300)
    fig.savefig(args.out + ".pdf")
    out = dict(
        stage="13_initial_overdensity",
        delta_definition="rho_SPH / (Omega_b rho_crit0 h^2 (1+z)^3) at the "
                         "early snapshot",
        statistic="mass-weighted median, exact over particles (no binning)",
        runs=R,
        ratio_separation_fct_over_cdm=ratio,
        histogram=dict(nbins=args.nbins, delta_lo=lo, delta_hi=hi,
                       mass_outside_range=0.0),
    )
    with open(args.out + ".json", "w") as fh:
        json.dump(out, fh, indent=2)
    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(log) + "\n")
    for ext in ("png", "pdf", "json", "txt"):
        print(f"written -> {args.out}.{ext}")


if __name__ == "__main__":
    main()
