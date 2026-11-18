#!/usr/bin/env python3
"""
Stage 14 - the gas census at every redshift, from SWIFT's own spectra.

WHY

Stage 09 builds P_gas from the snapshots because the OLD pair was run with
only matter-matter and cdm-cdm requested. Every run made since 2026-09
(the new pair, the ten low-resolution legs, the twelve more_power runs)
asked SWIFT for "gas-gas" as well, at 38 redshifts (17 for more_power),
and the new pair and the legs also for "gas-matter". The quantities that
decide the resolution test -- r_gas at matched particle mass -- and the
evolution of the drain with z are therefore already on disk, for every run,
at every redshift, with no snapshot to read.

WHAT IT MEASURES, for each (ref, test) pair and each common redshift

    r_mm(k)  = P_mm,test / P_mm,ref          gravity (+ baryons, small)
    r_gg(k)  = P_gg,test / P_gg,ref          what the forest sees
    D(k)     = r_gg / r_mm                    the drain: gas relative to matter
    b_gm(k)  = P_gm / P_mm   (each run)       gas-matter bias, if gas-matter exists

D = 1 means the gas follows the matter; D < 1 at low k is the large-scale
deficit, measured in 3D, with no tau_eff and no sightlines.

CONTROL

On the old pair stage 09 gave r_gg = 0.8055 at k = 0.15 Mpc^-1, z = 3. On the
NEW pair this stage should give a number close to it: same model, same
phases, same resolution, different run. That is the first check on whether
the new pair carries the same problem.

NORMALISATION -- READ THIS. SWIFT converts EVERY species to a density
contrast with the cosmic mean MATTER density, Omega_m rho_crit
(src/power_spectrum.c, `meanrho`, `invcellmean`), not with the species' own
mean. So for gas

    P_gas,SWIFT(k) = (M_gas / M_tot)^2 P_delta_gas(k),

and a run that converts more gas has a smaller raw P_gas for that reason
alone. Read raw, the new pair gives r_gg(0.15) = 0.285 at z = 3. This stage
divides it out with M_gas(z) / M_tot(z) from the run's own statistics.txt
(columns 3, 4, 5: redshift, total mass, gas mass), interpolated to the
redshift of each spectrum, and reports f_gas for both runs. With that, the
check against stage 09 on the old pair (0.8055) is meaningful. The cross
spectrum gas-matter carries one power of the same factor and is corrected
likewise. Matter needs nothing: M_tot is conserved.

UNITS. SWIFT writes k in internal units, Mpc^-1 comoving (UnitLength = 1 Mpc),
and P with the shot noise already subtracted. Ratios are dimensionless. For
gas in the high-conversion runs the subtracted shot noise is larger; the
column P_noise is carried and the ratio of signal to noise is reported at
the highest k kept (--snr-min sets the cut).

Usage
-----
    python stages/14_swift_spectra.py --pairs cdm40:fct40 \
        lr_zreio75_cdm:lr_zreio75_fct --out figures/swift_spectra
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from common import runs  # noqa: E402

KREPORT = (0.15, 0.3, 1.0, 2.7, 10.0, 27.0)       # Mpc^-1
ZREPORT = (5.0, 4.0, 3.0, 2.2, 2.0)


def load_kind(run_dir, kind):
    """{z: (k, P, Pnoise)} for one spectrum kind ('gas', 'matter', 'gas-matter')."""
    out = {}
    for f in sorted(glob.glob(os.path.join(run_dir, "power_spectra",
                                           f"power_{kind}_[0-9]*.txt"))):
        a = np.loadtxt(f, ndmin=2)
        if a.size == 0:
            continue
        z = round(float(a[0, 0]), 3)
        pn = a[:, 3] if a.shape[1] > 3 else np.zeros(len(a))
        out[z] = (a[:, 1], a[:, 2], pn)
    return out


def gas_fraction(run_dir):
    """f_gas(z) = M_gas / M_tot from statistics.txt, as a function of z."""
    f = os.path.join(run_dir, "statistics.txt")
    t = np.loadtxt(f, usecols=(3, 4, 5), comments="#", ndmin=2)
    z, idx = np.unique(t[:, 0], return_index=True)
    fg = t[idx, 2] / t[idx, 1]
    return lambda zz: float(np.interp(zz, z, fg))


def ratio(num, den):
    k1, p1, n1 = num
    k0, p0, n0 = den
    k = k0[(k0 >= max(k0.min(), k1.min())) & (k0 <= min(k0.max(), k1.max()))]
    a = np.interp(np.log(k), np.log(k1), p1)
    b = np.interp(np.log(k), np.log(k0), p0)
    sn = np.minimum(np.interp(np.log(k), np.log(k1), p1 / np.maximum(n1, 1e-300)),
                    np.interp(np.log(k), np.log(k0), p0 / np.maximum(n0, 1e-300)))
    return k, a / b, sn


def at(k, r, k0):
    if k0 < k.min() or k0 > k.max():
        return np.nan
    return float(np.interp(np.log(k0), np.log(k), r))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--pairs", nargs="+", required=True, help="ref:test, registry names or dirs")
    ap.add_argument("--snr-min", type=float, default=1.0,
                    help="drop k where P/P_noise < this in either run")
    ap.add_argument("--out", default="figures/swift_spectra")
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    rows, store = [], {}
    for pair in a.pairs:
        ref, test = pair.split(":")
        dr, dt = runs.run_dir(ref), runs.run_dir(test)
        S = {n: {kind: load_kind(d, kind) for kind in ("gas", "matter", "gas-matter")}
             for n, d in (("ref", dr), ("test", dt))}
        zs = sorted(set(S["ref"]["gas"]) & set(S["test"]["gas"])
                    & set(S["ref"]["matter"]) & set(S["test"]["matter"]), reverse=True)
        if not zs:
            print(f"!! {pair}: no common redshift with gas and matter spectra")
            continue
        print(f"== {ref} -> {test}: {len(zs)} redshifts, {zs[0]} .. {zs[-1]}")
        fgr, fgt = gas_fraction(dr), gas_fraction(dt)
        for z in zs:
            fr, ft = fgr(z), fgt(z)
            k, rgg, sng = ratio(S["test"]["gas"][z], S["ref"]["gas"][z])
            rgg = rgg * (fr / ft) ** 2          # undo SWIFT's Omega_m normalisation
            km, rmm, _ = ratio(S["test"]["matter"][z], S["ref"]["matter"][z])
            ok = sng >= a.snr_min
            k, rgg = k[ok], rgg[ok]
            rmm_k = np.interp(np.log(k), np.log(km), rmm)
            D = rgg / rmm_k
            key = f"{ref}__{test}__z{z:g}"
            store[key + "__k"] = k
            store[key + "__rgg"] = rgg
            store[key + "__rmm"] = rmm_k
            for n, fn in (("ref", fr), ("test", ft)):
                if z in S[n]["gas-matter"]:
                    kg, bgm, _ = ratio(S[n]["gas-matter"][z], S[n]["matter"][z])
                    store[f"{key}__bgm_{n}"] = np.interp(np.log(k), np.log(kg), bgm) / fn
            row = [ref, test, z, fr, ft] + [at(k, rgg, q) for q in KREPORT] \
                + [at(k, rmm_k, q) for q in KREPORT] + [at(k, D, q) for q in KREPORT]
            rows.append(row)

    np.savez(a.out + ".npz", **store)
    hdr = ["ref", "test", "z", "fgas_ref", "fgas_test"] + [f"rgg_{q:g}" for q in KREPORT] \
        + [f"rmm_{q:g}" for q in KREPORT] + [f"D_{q:g}" for q in KREPORT]
    with open(a.out + ".txt", "w") as fh:
        fh.write("# stage 14: SWIFT shot-noise-subtracted spectra, k in Mpc^-1 comoving\n")
        fh.write("# rgg = P_delta_gas ratio (SWIFT's Omega_m normalisation divided out with\n")
        fh.write("#       fgas = M_gas/M_tot from statistics.txt), rmm = P_matter ratio, D = rgg/rmm\n")
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join([r[0], r[1], f"{r[2]:g}"] + [f"{v:.4f}" for v in r[3:]]) + "\n")
    print(f"wrote {a.out}.txt and .npz ({len(rows)} rows)")

    # the short version on stdout: the drain at the largest scale, vs z
    print("\nD(k = 0.15 Mpc^-1) and r_gg(27) at the report redshifts")
    for r in rows:
        if any(abs(r[2] - zr) < 1e-6 for zr in ZREPORT):
            print(f"  {r[0]:>16s} -> {r[1]:<16s} z={r[2]:<4g} "
                  f"fgas={r[3]:.4f}/{r[4]:.4f} rgg(0.15)={r[5]:.4f} rmm(0.15)={r[11]:.4f} "
                  f"D(0.15)={r[17]:.4f} rgg(27)={r[10]:.4f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pairs = list(dict.fromkeys((r[0], r[1]) for r in rows))
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
        for (ref, test) in pairs:
            zz = np.array([r[2] for r in rows if (r[0], r[1]) == (ref, test)])
            dd = np.array([r[17] for r in rows if (r[0], r[1]) == (ref, test)])
            axs[0].plot(zz, dd, "o-", ms=3, label=f"{test}/{ref}")
            key = f"{ref}__{test}__z3"
            if key + "__k" in store:
                axs[1].semilogx(store[key + "__k"], store[key + "__rgg"] / store[key + "__rmm"],
                                label=f"{test}/{ref}")
        axs[0].axhline(1, c="k", lw=0.6)
        axs[0].set_xlabel("z"); axs[0].set_ylabel("D(k=0.15 Mpc$^{-1}$) = r$_{gg}$/r$_{mm}$")
        axs[0].invert_xaxis()
        axs[1].axhline(1, c="k", lw=0.6)
        axs[1].set_xlabel("k [Mpc$^{-1}$]"); axs[1].set_ylabel("D(k), z = 3")
        axs[0].legend(fontsize=6)
        fig.tight_layout(); fig.savefig(a.out + ".png", dpi=130)
        print(f"wrote {a.out}.png")
    except Exception as e:  # the table is the product; the figure is a convenience
        print("figure skipped:", e)


if __name__ == "__main__":
    main()
