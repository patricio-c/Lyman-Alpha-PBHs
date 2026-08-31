#!/usr/bin/env python
"""
compare_gas.py - Compara el contenido de gas entre dos corridas.

Motivo: los snapshots de CDM y FCT tienen el MISMO numero total de particulas
pero distinto reparto gas/DM (en el ejemplo, 47.976.746 particulas que en CDM
son gas en FCT son materia oscura). Con la carga 'masked' de monofonIC, eso
ocurre donde el esquema estandar daria masa bariónica negativa, o sea en las
celdas mas SUBDENSAS.

Si la corrida con menos gas esta perdiendo justamente los voids, el cociente
de P1D entre las dos no mide fisica: mide cobertura bariónica distinta. Este
script lo verifica sobre los archivos de LOS ya regenerados.

Que mira:
  1. Masa total de gas y su cociente respecto de Omega_b (cuanta materia
     bariónica esta efectivamente resuelta como gas).
  2. Distribucion de masas por particula (deberian ser todas iguales).
  3. PDF de sobredensidad del gas a lo largo de las LOS. Si a la corrida con
     menos gas le falta la cola SUBDENSA, la hipotesis queda confirmada.
  4. Densidad de particulas por LOS y su dispersion.

Uso:
    python compare_gas.py cdm40_z3.0.hdf5 fct40_z3.0.hdf5 \
        --labels CDM FCT --n-los 300
"""

import argparse

import h5py
import numpy as np

MPC_CM = 3.0856775814913673e24
MSUN_G = 1.98841e33


def read_meta(path):
    with h5py.File(path, "r") as f:
        c, hd, u = f["Cosmology"].attrs, f["Header"].attrs, f["Units"].attrs
        return dict(
            z=float(np.atleast_1d(c["Redshift"])[0]),
            h=float(np.atleast_1d(c["h"])[0]),
            omega_b=float(np.atleast_1d(c["Omega_b"])[0]),
            box=float(np.atleast_1d(hd["BoxSize"])[0]),
            u_l=float(np.atleast_1d(u["Unit length in cgs (U_L)"])[0]),
            u_m=float(np.atleast_1d(u["Unit mass in cgs (U_M)"])[0]),
        )


def gather(path, n_los):
    """Junta masas y densidades del gas de las primeras n_los lineas."""
    m = read_meta(path)
    mass, rho, hsml, counts = [], [], [], []
    with h5py.File(path, "r") as f:
        names = sorted(k for k in f if k.startswith("LOS_"))[:n_los]
        rho_c = float(np.atleast_1d(f[names[0]]["Densities"].attrs[
            "Conversion factor to CGS (not including cosmological corrections)"])[0])
        m_c = float(np.atleast_1d(f[names[0]]["Masses"].attrs[
            "Conversion factor to CGS (not including cosmological corrections)"])[0])
        for nm in names:
            # float32 -> float64 ANTES de convertir: 1e43 desborda en float32
            mass.append(f[nm]["Masses"][:].astype(np.float64))
            rho.append(f[nm]["Densities"][:].astype(np.float64))
            hsml.append(f[nm]["SmoothingLengths"][:].astype(np.float64))
            counts.append(f[nm]["Coordinates"].shape[0])
    return (m, np.concatenate(mass) * m_c, np.concatenate(rho) * rho_c,
            np.array(counts), len(names), np.concatenate(hsml))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file_a")
    ap.add_argument("file_b")
    ap.add_argument("--labels", nargs=2, default=["A", "B"])
    ap.add_argument("--n-los", type=int, default=300)
    args = ap.parse_args()

    out = {}
    for path, lab in zip([args.file_a, args.file_b], args.labels):
        out[lab] = gather(path, args.n_los)

    la, lb = args.labels
    ma = out[la][0]
    print(f"\nz = {ma['z']:.3f},  caja = {ma['box']:.4f} Mpc internos, "
          f"h = {ma['h']:.4f}\n")

    # --- 1. masa por particula --------------------------------------------
    print("1) Masa de las particulas de gas [g]:")
    for lab in args.labels:
        m = out[lab][1]
        print(f"   {lab:<6} media={m.mean():.5e}  desvio/media={m.std()/m.mean():.2e}"
              f"  min={m.min():.4e}  max={m.max():.4e}")
    ra = out[la][1].mean() / out[lb][1].mean()
    print(f"   cociente de masa media {la}/{lb} = {ra:.5f}")
    if abs(ra - 1) > 0.01:
        print("   >>> LAS MASAS DIFIEREN: la resolucion NO es la misma y el")
        print("       cociente de P1D no cancela el efecto de resolucion.")
    else:
        print("   ok: misma masa por particula, misma resolucion nominal.")

    # --- 2. densidad media del gas en la caja -----------------------------
    print("\n2) Densidad media comovil del gas [g/cm^3] y fraccion de Omega_b:")
    rho_b = 1.87847e-29 * ma["h"] ** 2 * ma["omega_b"]     # comovil
    print(f"   Omega_b implica rho_b comovil = {rho_b:.5e}")
    for lab in args.labels:
        meta, mass, rho, counts, nlos = out[lab][:5]
        # masa total de gas = N_gas * m_part ; N_gas del header del snapshot
        # aca solo tenemos las LOS, asi que estimamos con la densidad media
        # pesada por volumen SPH: <rho> = sum(m) / sum(m/rho)
        vol = (mass / rho).sum()
        rho_mean = mass.sum() / vol
        print(f"   {lab:<6} <rho>_V = {rho_mean:.5e}   "
              f"cociente vs Omega_b = {rho_mean/rho_b:.4f}")
    print("   (Si un run da muy por debajo de 1, parte de sus bariones no")
    print("    estan resueltos como gas: viven dentro de particulas de DM.)")

    # --- 3. PDF de sobredensidad ------------------------------------------
    print("\n3) Distribucion de sobredensidad del gas a lo largo de las LOS:")
    print(f"   {'percentil':>10}" + "".join(f"{l:>12}" for l in args.labels)
          + f"{'cociente':>11}")
    pcts = [0.1, 1, 5, 25, 50, 75, 95, 99]
    qa = np.percentile(out[la][2] / rho_b, pcts)
    qb = np.percentile(out[lb][2] / rho_b, pcts)
    for p, va, vb in zip(pcts, qa, qb):
        print(f"   {p:>9.1f}%{va:>12.4f}{vb:>12.4f}{va/vb:>11.3f}")
    print("\n   Lectura: si al run con menos gas le FALTA la cola subdensa,")
    print("   sus percentiles bajos van a estar muy por ENCIMA del otro.")
    print("   Eso significa que no tiene particulas en los voids, que es")
    print("   justo donde el forest transmite y donde vive la señal a k chico.")

    fa = (out[la][2] / rho_b < 0.1).mean()
    fb = (out[lb][2] / rho_b < 0.1).mean()
    print(f"\n   fraccion de particulas con Delta < 0.1: "
          f"{la} {100*fa:.2f}%   {lb} {100*fb:.2f}%")
    fa2 = (out[la][2] / rho_b < 0.01).mean()
    fb2 = (out[lb][2] / rho_b < 0.01).mean()
    print(f"   fraccion con Delta < 0.01:              "
          f"{la} {100*fa2:.2f}%   {lb} {100*fb2:.2f}%")

    # --- 4. particulas por LOS --------------------------------------------
    print("\n4) Particulas por linea de visión:")
    for lab in args.labels:
        c = out[lab][3]
        print(f"   {lab:<6} mediana={int(np.median(c))}  "
              f"min={c.min()}  max={c.max()}  ({out[lab][4]} LOS)")
    print(f"   cociente de medianas {la}/{lb} = "
          f"{np.median(out[la][3])/np.median(out[lb][3]):.4f}")

    # --- 5. veredicto ------------------------------------------------------
    # Contexto: la formacion estelar de QLA convierte en colisionales las
    # particulas por encima de over_density (1000 por defecto). Eso saca gas
    # DENSO a proposito y es lo esperado. Lo que NO seria aceptable es perder
    # gas SUBDENSO, porque ahi vive la transmision del forest.
    # --- 4b. longitudes de suavizado --------------------------------------
    print("\n4b) Longitudes de suavizado del gas [Mpc comovil]:")
    print(f"   {'':>8}{'p10':>10}{'p50':>10}{'p90':>10}")
    hs = {}
    for lab in args.labels:
        h_ = out[lab][5]
        hs[lab] = np.percentile(h_, [10, 50, 90])
        print(f"   {lab:<8}{hs[lab][0]:>10.4f}{hs[lab][1]:>10.4f}"
              f"{hs[lab][2]:>10.4f}")
    rh = hs[lb][1] / hs[la][1]
    nratio = np.median(out[la][3]) / np.median(out[lb][3])
    print(f"   cociente de h mediana {lb}/{la} = {rh:.4f}")
    print(f"\n   Prediccion SPH: si {lb} tiene {nratio:.2f}x menos gas en el")
    print(f"   IGM, su h deberia crecer como n^(-1/3) = {nratio**(1/3):.3f}x")
    print(f"   y el conteo por LOS caer solo {nratio**(1/3):.2f}x, no {nratio:.2f}x.")
    if rh < 1.05:
        print(f"   >>> h NO crecio ({rh:.3f}x). El gas removido NO estaba")
        print(f"       concentrado en halos, o el h ya estaba limitado por")
        print(f"       h_max. Revisar h_max del yml y la cola superior de h.")
    else:
        print(f"   >>> h crecio {rh:.3f}x, consistente con menos gas difuso.")

    print("\nVeredicto:")
    lo_a, lo_b = qa[0], qb[0]          # percentil 0.1
    hi_a, hi_b = qa[-1], qb[-1]        # percentil 99
    if abs(ra - 1) > 0.01:
        print("  Las masas de particula difieren -> problema de resolucion,")
        print("  el cociente de P1D no cancelaria ese efecto.")
    else:
        print(f"  cola SUBDENSA  (p0.1): {la} {lo_a:.4f}  {lb} {lo_b:.4f}"
              f"   cociente {lo_a/lo_b:.3f}")
        print(f"  cola DENSA     (p99) : {la} {hi_a:.3f}  {lb} {hi_b:.3f}"
              f"   cociente {hi_a/hi_b:.3f}")
        if lo_b > 2.0 * lo_a:
            print("\n  >>> Al run con menos gas le falta la cola SUBDENSA.")
            print("      Eso NO lo explica la formacion estelar de QLA, que")
            print("      remueve gas denso. El cociente de P1D estaria midiendo")
            print("      cobertura bariónica en los voids, no fisica.")
        elif hi_b < 0.5 * hi_a:
            print("\n  >>> Al run con menos gas le falta la cola DENSA, que es")
            print("      exactamente lo que hace la formacion estelar de QLA")
            print("      (over_density: 1000). Esperado y no problematico para")
            print("      el forest, que vive en Delta ~ 0.1-10.")
        else:
            print("\n  >>> Las dos cubren el mismo rango de densidad en las LOS.")
            print("      La diferencia de conteo no sesga ni voids ni halos.")
        # cobertura en el regimen del forest
        print(f"\n  Fraccion de particulas en el regimen del forest "
              f"(0.1 < Delta < 10):")
        for lab, q in [(la, out[la][2] / rho_b), (lb, out[lb][2] / rho_b)]:
            fr = ((q > 0.1) & (q < 10)).mean()
            print(f"    {lab:<6} {100*fr:.2f}%")
        print("  Si las dos tienen cobertura parecida en ese rango, el")
        print("  cociente de P1D sigue siendo interpretable.")


if __name__ == "__main__":
    main()
