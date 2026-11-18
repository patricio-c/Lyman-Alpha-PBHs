#!/usr/bin/env python3
"""
Stage 15 - order every model by how much extra power it has, and say how
much gas each one converted.

THE QUESTION

The suite has fifteen non-CDM initial power spectra: production FCT, the
six NB (broken spectrum) and six poisson runs of more_power, and the FCT of
the two low-resolution legs. To ask how the forest responds to extra power,
each needs ONE number for "how much extra", computed from the input spectrum
that monofonIC actually used, not from model parameters.

WHAT IT COMPUTES, from IC_dir/input_powerspec.txt (k in h/Mpc, P_dtot at
a_start), against a reference CDM spectrum on the same k grid:

  low_k      P/P_ref at k < 0.01 h/Mpc: must be 1. If it is not, the sigma_8
             normalisation has moved the large scales and that offset enters
             every ratio. Reported first because it is a precondition.
  dA_lnk     integral of (Delta^2 - Delta^2_ref) dln k, divided by the
             same integral for the reference: the extra area in Delta^2.
  dI1D(kp)   (1/2pi) integral_{kp}^{kcut} k (P - P_ref) dk, divided by the
             reference's own value: the extra power that the P1D kernel
             feeds into the flux at k_par = kp. Evaluated at the k_par of
             0.001, 0.003, 0.01, 0.03 s/km at z = 3. This is the quantity the
             forest integrates; the others are for comparison.
  sig_ratio  sigma_R / sigma_R,ref, top hat, R = 0.067 Mpc/h (1e8 Msun/h),
             the scale on which the 2026-09-16 entry found that conversion
             collapses both grid families onto one relation.

kcut is reported twice: the gas Nyquist of THIS run (pi N_gas^(1/3) / L,
from the IC header), and a COMMON cut, the smallest of those over all runs,
so that integrals over different resolutions can be compared.

THE RESPONSE side, from headers only (no particle is read):

  conv(z)    1 - N_gas(snapshot) / N_gas(IC), at the snapshot nearest each z
             in --zconv. By particle count, as in the 2026-09-04 census.

Usage
-----
    python stages/15_extra_power.py --out figures/extra_power
"""
import argparse
import glob
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from common import runs  # noqa: E402

ICB = "/data/contrib/pad_140/pcolazo/IC_lyman_alpha"
H = 0.681
HZ3 = 306.3077          # km/s/Mpc, physical, from Cosmology/H of the z = 3 LOS files
KPAR_SKM = (0.001, 0.003, 0.01, 0.03)
R_SIGMA = 0.067         # Mpc/h

# name: (IC directory, run registry name or dir, reference name)
MODELS = {
    "CDM_b40":   (f"{ICB}/CDM/6_10^6/B40", "cdm40", None),
    "FCT_b40":   (f"{ICB}/6_10^6/B40", "fct40", "CDM_b40"),
    "CDM_lr_c":  (f"{ICB}/low_resolution/CDM/c", "lr_zreio75_cdm", None),
    "FCT_lr_c":  (f"{ICB}/low_resolution/FCT/c", "lr_zreio75_fct", "CDM_lr_c"),
    "CDM_lr_a":  (f"{ICB}/low_resolution/CDM/a", "lr_a_cdm", None),
    "FCT_lr_a":  (f"{ICB}/low_resolution/FCT/a", "lr_a_fct", "CDM_lr_a"),
    # the CDM anchor of more_power: same TF as CDM_b40, the grid's own sc load
    "CDM_sc":    (f"{ICB}/more_power/CDM_sc", "/data/contrib/pad_140/pcolazo/lyman/more_power/CDM_sc", "CDM_b40"),
}
for _f in ("NB", "poisson"):
    for _i in range(1, 7):
        MODELS[f"{_f}_{_i}"] = (f"{ICB}/more_power/{_f}_{_i}", f"{_f}_{_i}", "CDM_b40")


def ic_file(d):
    c = [f for f in glob.glob(os.path.join(d, "IC_*.hdf5")) if "z5" not in os.path.basename(f)]
    return sorted(c)[0] if c else None


def ngas_box(f):
    with h5py.File(f, "r") as h:
        a = h["Header"].attrs
        n = np.array(a["NumPart_Total"], dtype=np.int64)
        if "NumPart_Total_HighWord" in a:
            n = n + (np.array(a["NumPart_Total_HighWord"], dtype=np.int64) << 32)
        L = float(np.atleast_1d(a["BoxSize"])[0])
        z = float(np.atleast_1d(a["Redshift"])[0]) if "Redshift" in a else np.nan
    return int(n[0]), L, z


def snapshots(run_dir):
    """virtual snapshot files only, with their redshift and gas count"""
    out = []
    for f in sorted(glob.glob(os.path.join(run_dir, "*_[0-9][0-9][0-9][0-9]", "*.hdf5"))
                    + glob.glob(os.path.join(run_dir, "*_[0-9][0-9][0-9][0-9].hdf5"))):
        b = os.path.basename(f)
        if not b.endswith(".hdf5") or b.count(".") > 1 or b.startswith("los_"):
            continue
        try:
            n, _, z = ngas_box(f)
            out.append((z, n, f))
        except Exception:
            pass
    return out


trapz = getattr(np, "trapezoid", None) or np.trapz


def tophat(x):
    return 3 * (np.sin(x) - x * np.cos(x)) / x**3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zconv", type=float, nargs="+", default=[7.0, 5.0, 3.0, 2.0])
    ap.add_argument("--out", default="figures/extra_power")
    a = ap.parse_args()

    spec, info = {}, {}
    for m, (icd, run, ref) in MODELS.items():
        p = os.path.join(icd, "input_powerspec.txt")
        f = ic_file(icd)
        if not os.path.exists(p) or f is None:
            print(f"-- {m}: no IC yet in {icd}, skipped")
            continue
        t = np.loadtxt(p)
        try:
            ng, _, _ = ngas_box(f)
        except (KeyError, OSError) as e:      # an IC still being written has no Header
            print(f"-- {m}: {os.path.basename(f)} not readable yet ({e.__class__.__name__}), skipped")
            continue
        spec[m] = (t[:, 0], t[:, 1])
        info[m] = dict(ng=ng, run=run, ref=ref, icfile=f)
    # L from the box: all 40 Mpc/h except leg a (80)
    for m in info:
        info[m]["L"] = 80.0 if m.endswith("_lr_a") else 40.0
        info[m]["kny"] = np.pi * round(info[m]["ng"] ** (1 / 3)) / info[m]["L"]   # h/Mpc
    kcom = min(v["kny"] for v in info.values())
    print(f"common cut: gas Nyquist min = {kcom:.2f} h/Mpc")

    kpar_hmpc = [q * HZ3 / 4.0 / H for q in KPAR_SKM]   # s/km -> 1/Mpc comoving -> h/Mpc

    rows = []
    for m, v in info.items():
        ref = v["ref"] or m
        k, P = spec[m]
        kr, Pr = spec[ref]
        Pr = np.exp(np.interp(np.log(k), np.log(kr), np.log(Pr)))
        R = P / Pr
        lowk = float(np.median(R[k < 0.01]))
        res = dict(model=m, ref=ref, ng=v["ng"], kny=v["kny"], low_k=lowk)
        for cutname, kc in (("own", v["kny"]), ("com", kcom)):
            s = k <= kc
            lk = np.log(k[s])
            D2, D2r = k[s]**3 * P[s] / (2 * np.pi**2), k[s]**3 * Pr[s] / (2 * np.pi**2)
            res[f"dA_lnk_{cutname}"] = trapz(D2 - D2r, lk) / trapz(D2r, lk)
            for q, kp in zip(KPAR_SKM, kpar_hmpc):
                u = s & (k >= kp)
                num = trapz(k[u] * (P[u] - Pr[u]), k[u])
                den = trapz(k[u] * Pr[u], k[u])
                res[f"dI1D_{q:g}_{cutname}"] = num / den
        x = k * R_SIGMA
        w = tophat(x) ** 2 * k**3 / (2 * np.pi**2)
        res["sig_ratio"] = np.sqrt(trapz(w * P, np.log(k)) / trapz(w * Pr, np.log(k)))
        # response: conversion fraction from headers
        try:
            snaps = snapshots(runs.run_dir(v["run"]))
        except SystemExit:
            snaps = []
        for zt in a.zconv:
            if snaps:
                z, n, _ = min(snaps, key=lambda s_: abs(s_[0] - zt))
                res[f"conv_z{zt:g}"] = 1 - n / v["ng"] if abs(z - zt) < 0.05 else np.nan
            else:
                res[f"conv_z{zt:g}"] = np.nan
        rows.append(res)

    keys = list(rows[0].keys())
    with open(a.out + ".txt", "w") as fh:
        fh.write("# stage 15: extra linear power vs reference CDM (input_powerspec.txt, a_start)\n")
        fh.write(f"# k_par(s/km) -> h/Mpc at z=3: "
                 + ", ".join(f"{q:g}->{kp:.3f}" for q, kp in zip(KPAR_SKM, kpar_hmpc)) + "\n")
        fh.write(f"# common cut {kcom:.3f} h/Mpc; conv = 1 - N_gas/N_gas(IC)\n")
        fh.write("\t".join(keys) + "\n")
        for r in sorted(rows, key=lambda r_: r_.get("dI1D_0.003_com", 0)):
            fh.write("\t".join(str(r[q]) if isinstance(r[q], str) else
                               (f"{r[q]:d}" if isinstance(r[q], int) else f"{r[q]:.6g}")
                               for q in keys) + "\n")
    print(open(a.out + ".txt").read())


if __name__ == "__main__":
    main()
