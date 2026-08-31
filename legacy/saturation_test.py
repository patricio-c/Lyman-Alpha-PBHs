#!/usr/bin/env python
"""
saturation_test.py - Tests de saturacion sobre los caches de tau.

Hipotesis a testear
-------------------
QLA convierte en estrellas todo el gas con Delta > 1000 y T < 1e5 K. La FCT
pierde asi el 49.8% de su gas contra el 14.1% de la CDM. Ese gas es el que
produce los sistemas SATURADOS (F ~ 0). Bolton et al. 2017 (Sherwood,
conclusiones) reportan que un modelo de subgrilla que retiene ese gas produce
mas sistemas saturados y AUMENTA el P1D a escalas grandes. Leido al reves:
sacar gas denso baja el P1D a k chico. La FCT saca 3.5x mas -> deficit.

Dos tests, ambos desde el cache, sin re-extraer:

  stats   Cuantos pixeles y cuantos sistemas saturados hay en cada corrida,
          con A = 1 (sin reescalar) y con el reescalado aplicado. Si la
          hipotesis es correcta: FCT tiene MENOS saturacion que CDM con A=1.

  ratio   Recalcula el cociente P1D despues de sacar los sistemas saturados.
          Si el deficit a k chico desaparece, queda confirmado. Y no es un
          parche: DESI DR1 enmascara DLAs y modela HCDs, asi que esto alinea
          el analisis de la simulacion con el de los datos.

Dos formas de sacarlos, ambas implementadas:

  --mode pair-reject   Descarta la LOS ENTERA si tiene un sistema saturado en
                       CUALQUIERA de las dos corridas. Como las LOS son las
                       mismas (mismo seed, mismas posiciones), el descarte es
                       pareado: se comparan exactamente los mismos rayos. No
                       introduce funcion ventana. Es el modo limpio.

  --mode mask          Enmascara solo los pixeles del sistema mas un padding,
                       poniendoles delta_F = 0. Conserva mas estadistica pero
                       mete una ventana que suprime potencia. Como se aplica
                       igual a las dos, se cancela en gran medida en el
                       cociente, pero no exactamente. Es el modo comparable a
                       lo que hace DESI.

Ejemplos
--------
    python saturation_test.py stats --a cache_cdm.npz --b cache_fct.npz \\
        --tau-eff 0.3719

    python saturation_test.py ratio --a cache_cdm.npz --b cache_fct.npz \\
        --tau-eff 0.3719 --mode pair-reject \\
        --tau-thresh 5 10 20 --out ratio_masked.png
"""

import argparse
import sys

import numpy as np
from scipy.ndimage import maximum_filter1d

from sherwood_los import flux_power_1d, rescale_tau, tau_eff

LAMBDA_LYA = 1215.6701
C_KMS = 2.99792458e5


# ---------------------------------------------------------------------------
# carga
# ---------------------------------------------------------------------------

def load_cache(path):
    d = np.load(path, allow_pickle=True)
    return dict(tau=np.asarray(d["tau"], dtype=np.float64),
                dv=float(d["dv"]), z=float(d["z"]),
                box_kms=float(d["box_kms"]), h=float(d["h"]),
                hz=float(d["hz"]), label=str(d["label"]))


def k_conversions(z, hz, h):
    return hz / (1.0 + z), hz / ((1.0 + z) * h)


# ---------------------------------------------------------------------------
# deteccion de sistemas saturados
# ---------------------------------------------------------------------------

def saturated_mask(tau, tau_thresh, pad_px=0):
    """Mascara booleana (nlos, npix) de pixeles saturados, dilatada pad_px.

    La dilatacion se hace con maximum_filter1d en modo 'wrap' porque las LOS
    son periodicas. Es O(N) con ventana deslizante, no O(N*pad).
    """
    m = tau > tau_thresh
    if pad_px > 0:
        m = maximum_filter1d(m, size=2 * pad_px + 1, axis=1, mode="wrap")
    return m


def count_systems(mask):
    """Numero de rachas contiguas por LOS, respetando periodicidad.

    Una racha empieza donde mask pasa de False a True mirando el pixel previo
    con wraparound. Si una LOS esta entera saturada no hay flanco: se cuenta
    como un sistema aparte.
    """
    prev = np.roll(mask, 1, axis=1)
    starts = mask & ~prev
    n = starts.sum(axis=1)
    full = mask.all(axis=1)
    n = np.where(full, 1, n)
    return n


def system_lengths(mask):
    """Longitudes de todas las rachas, en pixeles. Loop por LOS con rachas."""
    out = []
    rows = np.where(mask.any(axis=1))[0]
    for i in rows:
        row = mask[i]
        if row.all():
            out.append(row.size)
            continue
        # rotar para que arranque en un False, asi ninguna racha cruza el borde
        j = np.argmin(row)
        r = np.roll(row, -j)
        d = np.diff(np.concatenate(([0], r.view(np.int8), [0])))
        out.extend(np.flatnonzero(d == -1) - np.flatnonzero(d == 1))
    return np.asarray(out)


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------

def cmd_stats(args):
    A, B = load_cache(args.a), load_cache(args.b)
    la, lb = args.labels

    print(f"z = {A['z']:.4f}   dv = {A['dv']:.4f} km/s   "
          f"{A['tau'].shape[0]} LOS x {A['tau'].shape[1]} px\n")

    print("=" * 68)
    print("TEST 1a. Saturacion SIN reescalar (A = 1). Campo tal cual sale")
    print("         de la simulacion. Aca es donde se ve el efecto fisico.")
    print("=" * 68)
    _sat_table(A, B, la, lb, scale=None, thresholds=args.tau_thresh,
               dv=A["dv"])

    print()
    print("=" * 68)
    print(f"TEST 1b. Saturacion DESPUES de reescalar a "
          f"tau_eff = {args.tau_eff:.4f}")
    print("         Aca el reescalado ya movio todo. Sirve para entender la")
    print("         PDF, no para el argumento fisico.")
    print("=" * 68)
    _sat_table(A, B, la, lb, scale=args.tau_eff, thresholds=args.tau_thresh,
               dv=A["dv"])

    print("\nLectura: si en el TEST 1a la FCT tiene MENOS pixeles y MENOS")
    print("sistemas saturados que la CDM, la hipotesis queda confirmada y el")
    print("mecanismo de Bolton+2017 aplica directo a tu cociente.")


def _sat_table(A, B, la, lb, scale, thresholds, dv):
    for D, lab in [(A, la), (B, lb)]:
        t = D["tau"]
        if scale is not None:
            a_ = rescale_tau(t, scale)
            t = a_ * t
            note = f"  [A = {a_:.5f}]"
        else:
            note = f"  [tau_eff nominal = {tau_eff(t):.5f}]"
        print(f"\n{lab}{note}")
        print(f"  {'tau>':>6}{'F<':>10}{'%px':>10}{'sistemas':>11}"
              f"{'/LOS':>8}{'<len> km/s':>12}")
        for th in thresholds:
            m = t > th
            ns = count_systems(m)
            L = system_lengths(m)
            mlen = L.mean() * dv if L.size else 0.0
            print(f"  {th:6.1f}{np.exp(-th):10.2e}{100*m.mean():10.4f}"
                  f"{ns.sum():11d}{ns.mean():8.3f}{mlen:12.2f}")


# ---------------------------------------------------------------------------
# ratio con enmascarado
# ---------------------------------------------------------------------------

def _p1d_pair_reject(A, B, target, tau_thresh, pad_px):
    """Descarta la LOS entera si tiene sistema saturado en cualquiera de las dos.

    Ojo: el umbral se aplica sobre tau SIN reescalar, porque el reescalado
    depende de que LOS sobrevivan y eso seria circular. Se itera una vez:
    seleccionar, reescalar sobre lo que queda, y listo.
    """
    mA = saturated_mask(A["tau"], tau_thresh, pad_px)
    mB = saturated_mask(B["tau"], tau_thresh, pad_px)
    keep = ~(mA.any(axis=1) | mB.any(axis=1))
    if keep.sum() < 50:
        sys.exit(f"Solo quedan {keep.sum()} LOS con tau>{tau_thresh}. "
                 "Subi el umbral.")
    out = []
    for D in (A, B):
        t = D["tau"][keep]
        a_ = rescale_tau(t, target)
        k, p, e = flux_power_1d(np.exp(-a_ * t), D["dv"])
        out.append((k, p, e, a_))
    info = dict(nkeep=int(keep.sum()), ntot=int(keep.size),
                fkeep=float(keep.mean()), fmask=0.0)
    return out, info


def _p1d_mask(A, B, target, tau_thresh, pad_px):
    """Enmascara pixeles saturados poniendoles el flujo medio (delta_F = 0).

    La misma mascara UNION se aplica a las dos corridas. Usar la union y no
    la mascara propia de cada una es esencial: si cada una enmascara sus
    propios pixeles, se comparan conjuntos distintos de gas y el cociente
    deja de significar algo.
    """
    m = (saturated_mask(A["tau"], tau_thresh, pad_px)
         | saturated_mask(B["tau"], tau_thresh, pad_px))
    good = ~m
    out = []
    for D in (A, B):
        a_ = rescale_tau(D["tau"][good], target)
        F = np.exp(-a_ * D["tau"])
        F[m] = F[good].mean()          # delta_F = 0 en lo enmascarado
        k, p, e = flux_power_1d(F, D["dv"])
        out.append((k, p, e, a_))
    info = dict(nkeep=int(good.shape[0]), ntot=int(good.shape[0]),
                fkeep=1.0, fmask=float(m.mean()))
    return out, info


def _ratio_from(out):
    (ka, pa, ea, Aa), (kb, pb, eb, Ab) = out
    lo, hi = max(ka[1], kb[1]), min(ka[-1], kb[-1]) * 0.5
    sel = (ka >= lo) & (ka <= hi)
    kk = ka[sel]
    pbi, ebi = np.interp(kk, kb, pb), np.interp(kk, kb, eb)
    r = pbi / pa[sel]
    rel = np.sqrt((ebi / pbi) ** 2 + (ea[sel] / pa[sel]) ** 2)
    return kk, r, rel, Aa, Ab


def cmd_ratio(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    A, B = load_cache(args.a), load_cache(args.b)
    if abs(A["z"] - B["z"]) > 0.02:
        sys.exit(f"Redshifts distintos: {A['z']} vs {B['z']}")
    la, lb = args.labels
    z, h, hz = A["z"], A["h"], A["hz"]
    to_mpc, to_hmpc = k_conversions(z, hz, h)
    pad_px = int(round(args.pad_kms / A["dv"]))

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.axhline(1.0, color="k", ls="--", lw=1)
    R_z = C_KMS * 0.8 / ((1.0 + z) * LAMBDA_LYA)
    ax.axvspan(1.0e-3, 0.5 * np.pi / R_z, color="green", alpha=0.07, zorder=0)

    # linea base: sin enmascarar nada
    base = []
    for D in (A, B):
        a_ = rescale_tau(D["tau"], args.tau_eff)
        k, p, e = flux_power_1d(np.exp(-a_ * D["tau"]), D["dv"])
        base.append((k, p, e, a_))
    kk0, r0, rel0, Aa0, Ab0 = _ratio_from(base)
    ax.fill_between(kk0, r0 * (1 - rel0), r0 * (1 + rel0),
                    color="0.6", alpha=0.3)
    ax.plot(kk0, r0, color="k", lw=2.2, label="no masking")
    print(f"padding = {args.pad_kms:.1f} km/s = {pad_px} px    "
          f"modo = {args.mode}\n")
    print(f"{'caso':>18}{'A('+la+')':>10}{'A('+lb+')':>10}"
          f"{'LOS':>8}{'%mask':>8}{'r(0.005)':>10}{'r(0.01)':>9}")
    print(f"{'sin mascara':>18}{Aa0:10.4f}{Ab0:10.4f}"
          f"{A['tau'].shape[0]:8d}{0.0:8.2f}"
          f"{np.interp(0.005, kk0, r0):10.4f}"
          f"{np.interp(0.010, kk0, r0):9.4f}")

    fn = _p1d_pair_reject if args.mode == "pair-reject" else _p1d_mask
    for th, col in zip(args.tau_thresh, ["C0", "C2", "C3", "C4", "C5"]):
        out, info = fn(A, B, args.tau_eff, th, pad_px)
        kk, r, rel, Aa, Ab = _ratio_from(out)
        ax.plot(kk, r, color=col, lw=1.6,
                label=rf"masked $\tau > {th:g}$  "
                      rf"($F < {np.exp(-th):.1e}$)")
        print(f"{'tau > '+format(th,'g'):>18}{Aa:10.4f}{Ab:10.4f}"
              f"{info['nkeep']:8d}{100*info['fmask']:8.2f}"
              f"{np.interp(0.005, kk, r):10.4f}"
              f"{np.interp(0.010, kk, r):9.4f}")

    if args.kt:
        ax.axvline(args.kt / to_mpc, color="C3", ls=":", lw=1.2)
    ax.axvline(2.0 * np.pi / A["box_kms"], color="grey", ls=":", lw=1.2)

    ax.set(xscale="log", xlabel=r"$k$  [s km$^{-1}$]",
           ylabel=rf"$P({lb})\,/\,P({la})$",
           title=(f"z = {z:.1f}, saturated systems removed "
                  f"({args.mode})"))
    if args.xlim:
        ax.set_xlim(*args.xlim)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    ax.grid(alpha=0.2, which="both")
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    tx = ax.secondary_xaxis("top", functions=(lambda x: x * to_hmpc,
                                              lambda x: x / to_hmpc))
    tx.set_xlabel(r"$k$  [$h$ Mpc$^{-1}$]")

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nescrito -> {args.out}")
    print("\nLectura: si r(0.005) sube hacia 1 al subir la agresividad del")
    print("enmascarado, el deficit a k chico viene de los sistemas saturados")
    print("que le faltan a la FCT, o sea de la conversion de gas de QLA.")
    print("Si no se mueve, es otra cosa y hay que seguir buscando.")


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("stats")
    s.add_argument("--a", required=True)
    s.add_argument("--b", required=True)
    s.add_argument("--tau-eff", type=float, required=True)
    s.add_argument("--tau-thresh", nargs="*", type=float,
                   default=[2.0, 5.0, 10.0, 20.0, 50.0])
    s.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    s.set_defaults(func=cmd_stats)

    r = sub.add_parser("ratio")
    r.add_argument("--a", required=True)
    r.add_argument("--b", required=True)
    r.add_argument("--tau-eff", type=float, required=True)
    r.add_argument("--tau-thresh", nargs="*", type=float,
                   default=[5.0, 10.0, 20.0])
    r.add_argument("--pad-kms", type=float, default=0.0,
                   help="ensanchar la mascara este tanto a cada lado, para "
                        "cubrir las alas de la linea. DESI usa algo asi.")
    r.add_argument("--mode", choices=["pair-reject", "mask"],
                   default="pair-reject")
    r.add_argument("--labels", nargs=2, default=["CDM", "FCT"])
    r.add_argument("--kt", type=float, default=None)
    r.add_argument("--xlim", nargs=2, type=float, default=None)
    r.add_argument("--ylim", nargs=2, type=float, default=None)
    r.add_argument("--out", required=True)
    r.set_defaults(func=cmd_ratio)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
