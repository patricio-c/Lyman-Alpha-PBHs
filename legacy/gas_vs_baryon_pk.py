#!/usr/bin/env python
"""
gas_vs_baryon_pk.py - Test directo de la hipotesis "QLA se llevo la potencia
de gran escala".

Mide en 3D, desde el snapshot, el P(k) de:
    gas          -> solo PartType0 (lo que ve el bosque)
    convertidas  -> PartType1 con masa de gas (lo que QLA saco del gas)
    bariones     -> gas + convertidas (el campo bariónico completo)
    materia      -> todo PartType0 + PartType1

y devuelve los cocientes FCT/CDM. La lectura:

  * cociente(materia) = 1 y cociente(bariones) = 1, pero cociente(gas) < 1
        -> el deficit ES la remocion de gas por QLA.
  * cociente(gas) tambien = 1 a k bajo
        -> QLA no toca la gran escala; el deficit viene de la temperatura
           o de la proyeccion no lineal 3D->1D.

No necesita LOS ni extraccion. Uso:

    python gas_vs_baryon_pk.py \\
        --cdm ../lyman/cdm-box-40-1024/cdm-40-m6-lyman_0003.hdf5 \\
        --fct ../lyman/2-fct-box-40-1024/<basename>_0003.hdf5 \\
        --ngrid 256 --out gas_vs_baryon.png
"""
import argparse

import h5py
import numpy as np

CHUNK = 4_000_000


def cic(coords_iter, ngrid, box):
    """Deposicion CIC con lectura por chunks. Devuelve la grilla de masa."""
    grid = np.zeros((ngrid, ngrid, ngrid), dtype=np.float64)
    total = 0.0
    total2 = 0.0
    for pos, m in coords_iter:
        f = (pos / box) * ngrid
        f -= np.floor(f / ngrid) * ngrid           # periodicidad
        i0 = np.floor(f).astype(np.int64)
        d = f - i0
        i0 %= ngrid
        i1 = (i0 + 1) % ngrid
        flat = grid.reshape(-1)
        n2 = ngrid * ngrid
        for a in (0, 1):
            wx = d[:, 0] if a else 1.0 - d[:, 0]
            ix = i1[:, 0] if a else i0[:, 0]
            for b in (0, 1):
                wy = d[:, 1] if b else 1.0 - d[:, 1]
                iy = i1[:, 1] if b else i0[:, 1]
                wxy = wx * wy
                base = ix * n2 + iy * ngrid
                for c in (0, 1):
                    wz = d[:, 2] if c else 1.0 - d[:, 2]
                    iz = i1[:, 2] if c else i0[:, 2]
                    flat += np.bincount(base + iz, weights=m * wxy * wz,
                                        minlength=flat.size)
        total += float(m.sum())
        total2 += float((m ** 2).sum())
    return grid, total, total2


def read_parts(path, ptype, mass_sel=None):
    """Itera (coords, masses) por chunks. mass_sel filtra por masa."""
    with h5py.File(path, "r") as f:
        key = f"PartType{ptype}"
        if key not in f:
            return
        n = f[key]["Masses"].shape[0]
        for s in range(0, n, CHUNK):
            e = min(s + CHUNK, n)
            m = f[key]["Masses"][s:e].astype(np.float64)
            p = f[key]["Coordinates"][s:e].astype(np.float64)
            if mass_sel is not None:
                k = mass_sel(m)
                m, p = m[k], p[k]
            if m.size:
                yield p, m


def mass_populations(path):
    """Separa las dos poblaciones de PartType1 por masa."""
    with h5py.File(path, "r") as f:
        m1 = f["PartType1"]["Masses"][:1_000_000].astype(np.float64)
        m0 = f["PartType0"]["Masses"][:1000].astype(np.float64)
    m_gas = float(np.median(m0))
    # las convertidas conservan la masa de gas; la DM original pesa
    # Omega_cdm/Omega_b veces mas. Umbral robusto a masas de IC variables.
    thr = 2.0 * m_gas
    if (m1 < thr).sum() == 0:
        return m_gas, None                       # no hay convertidas aca
    return m_gas, thr


def pk_from_grid(grid, total, total2, box, ngrid, subtract_shot=True):
    delta = grid / (total / ngrid ** 3) - 1.0
    dk = np.fft.rfftn(delta) * (box / ngrid) ** 3 / box ** 1.5
    p3 = (np.abs(dk) ** 2).astype(np.float64)
    kf = 2.0 * np.pi / box
    kx = np.fft.fftfreq(ngrid, d=1.0 / ngrid) * kf
    kz = np.fft.rfftfreq(ngrid, d=1.0 / ngrid) * kf
    kk = np.sqrt(kx[:, None, None] ** 2 + kx[None, :, None] ** 2
                 + kz[None, None, :] ** 2)
    # deconvolucion de la ventana CIC: P_med = prod_i sinc^4 * P_true
    nx = np.fft.fftfreq(ngrid, d=1.0 / ngrid)
    nz = np.fft.rfftfreq(ngrid, d=1.0 / ngrid)
    w = (np.sinc(nx[:, None, None] / ngrid) ** 4
         * np.sinc(nx[None, :, None] / ngrid) ** 4
         * np.sinc(nz[None, None, :] / ngrid) ** 4)
    p3 /= w
    edges = np.logspace(np.log10(kf * 0.9), np.log10(kf * ngrid / 2), 40)
    nb = len(edges) - 1
    idx = np.digitize(kk.ravel(), edges) - 1
    ok = (idx >= 0) & (idx < nb)
    idx = idx[ok]
    pw, kb, cnt = np.zeros(nb), np.zeros(nb), np.zeros(nb)
    np.add.at(pw, idx, p3.ravel()[ok])
    np.add.at(kb, idx, kk.ravel()[ok])
    np.add.at(cnt, idx, 1.0)
    good = cnt > 0
    # ruido de disparo para trazadores pesados: V * sum(m^2) / (sum m)^2
    shot = box ** 3 * total2 / total ** 2 if subtract_shot else 0.0
    p = pw[good] / cnt[good] - shot
    return kb[good] / cnt[good], p, shot


def measure(path, ngrid, label):
    with h5py.File(path, "r") as f:
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
    m_gas, thr = mass_populations(path)
    out = {}
    fields = {"gas": (0, None)}
    if thr is not None:
        fields["conv"] = (1, lambda m, t=thr: m < t)
        fields["dm"] = (1, lambda m, t=thr: m >= t)
    else:
        fields["dm"] = (1, None)
        with h5py.File(path, "r") as f:          # convertidas como PartType4
            has4 = "PartType4" in f and f["PartType4"]["Masses"].shape[0] > 0
        if has4:
            fields["conv"] = (4, None)
    grids = {}
    for nm, (pt, sel) in fields.items():
        g, tot, tot2 = cic(read_parts(path, pt, sel), ngrid, box)
        grids[nm] = (g, tot, tot2)
        print(f"  {label} {nm:5s}: M = {tot:.4e}")
    combos = {"gas": ["gas"], "matter": list(grids)}
    if "conv" in grids:
        combos["baryon"] = ["gas", "conv"]
    for nm, parts in combos.items():
        g = sum(grids[p][0] for p in parts)
        tot = sum(grids[p][1] for p in parts)
        tot2 = sum(grids[p][2] for p in parts)
        kk, pp, shot = pk_from_grid(g, tot, tot2, box, ngrid)
        print(f"  {label} {nm:6s}: P_shot = {shot:.3e} Mpc^3")
        out[nm] = (kk, pp)
    return out, box


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdm", required=True)
    ap.add_argument("--fct", required=True)
    ap.add_argument("--ngrid", type=int, default=256)
    ap.add_argument("--out", default="gas_vs_baryon.png")
    args = ap.parse_args()

    print("CDM:"); A, box = measure(args.cdm, args.ngrid, "CDM")
    print("FCT:"); B, _ = measure(args.fct, args.ngrid, "FCT")

    common = [c for c in ("matter", "baryon", "gas") if c in A and c in B]
    k = A[common[0]][0]
    print(f"\n{'k [1/Mpc]':>10}" + "".join(f"{c:>10}" for c in common))
    for i in range(0, len(k), max(1, len(k) // 16)):
        row = "".join(f"{B[c][1][i]/A[c][1][i]:10.3f}" for c in common)
        print(f"{k[i]:10.4f}{row}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    for c, col in zip(common, ["k", "C2", "C3"]):
        ax.plot(A[c][0], B[c][1] / A[c][1], color=col, lw=1.8, label=c)
    ax.axhline(1, color="0.5", ls="--", lw=1)
    ax.set(xscale="log", xlabel=r"$k$ [Mpc$^{-1}$]",
           ylabel="P(FCT) / P(CDM)", ylim=(0.7, 2.0),
           title="3D power: gas vs full baryons vs matter, z = 3")
    ax.legend(frameon=False); ax.grid(alpha=0.2, which="both")
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nescrito -> {args.out}")


if __name__ == "__main__":
    main()

