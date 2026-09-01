#!/usr/bin/env python
"""
relos.py - Regenera archivos de LOS de SWIFT desde snapshots completos.

Motivo: los archivos los_*.hdf5 originales fueron escritos con
range_when_shooting_down_* = [0, 40] en una caja de 58.7372 Mpc (unidades
internas), asi que a cada linea de visión le falta el 31.9% de las particulas.
Las POSICIONES de los rayos, en cambio, estan intactas en los atributos de
cada grupo LOS. Este script:

  1. Lee eje y posicion de cada rayo del archivo de LOS viejo.
  2. Recorre el snapshot completo (soporta snapshots distribuidos en varios
     archivos) en bloques, y selecciona para cada rayo las particulas con
     b < gamma*h - el MISMO criterio que usa SWIFT (verificado: en los
     archivos originales max(b/(gamma*h)) = 1.0000).
  3. Escribe un archivo nuevo con el mismo esquema: mismos grupos LOS_XXXX,
     mismos atributos de grupo, mismos datasets (los que existan en el
     snapshot), mismos atributos de conversion de unidades.

Uso:
    python relos.py --old-los los_0004.hdf5 \
        --snapshot 'cdm-40-m6-lyman_0004*.hdf5' \
        --out los_0004_fixed.hdf5

    # con rayos nuevos uniformes sobre toda la caja en vez de los viejos
    # (los viejos solo muestrean [0,40] Mpc transversalmente por el mismo
    #  bug en allowed_los_range):
    python relos.py --old-los los_0004.hdf5 --snapshot '...' \
        --out los_0004_uniform.hdf5 --uniform 2048 --seed 12345

Notas de rendimiento: el costo esta dominado por la lectura del snapshot
(decenas de GB), no por el computo; por eso Python + h5py + cKDTree rinde
igual que C aca. La seleccion usa un arbol 2D periodico por eje sobre las
posiciones de los rayos, asi que escala como O(Npart * log Nrayos).
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

M4_GAMMA_DEFAULT = 1.825742

# Campos que la pipeline de extraccion necesita si o si.
ESSENTIAL = ["Coordinates", "Densities", "Temperatures", "Velocities",
             "SmoothingLengths", "Masses"]

_RAY_ATTR_CANDIDATES = [
    ("Xpos", "Ypos"), ("Xproj", "Yproj"), ("x_proj", "y_proj"),
    ("Position_x", "Position_y"),
]


# ---------------------------------------------------------------------------
# Lectura de los rayos del archivo viejo
# ---------------------------------------------------------------------------

def read_rays(old_los_path: str):
    """
    Devuelve (rays, group_attrs, z_old, kernel_gamma, box).

    rays: lista de (nombre, eje, pos2d) con pos2d en unidades internas sobre
    los DOS ejes transversales en orden creciente de indice de eje.
    El eje se detecta por la dispersion de las particulas (las truncadas
    alcanzan: 40 Mpc de extension axial contra <1 Mpc transversal).
    """
    rays, gattrs = [], {}
    with h5py.File(old_los_path, "r") as f:
        z_old = float(np.atleast_1d(f["Cosmology"].attrs["Redshift"])[0])
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        try:
            kg = float(np.atleast_1d(f["HydroScheme"].attrs["Kernel gamma"])[0])
        except KeyError:
            kg = M4_GAMMA_DEFAULT
        names = sorted(k for k in f.keys() if k.startswith("LOS_"))
        for nm in names:
            g = f[nm]
            c = g["Coordinates"][:].astype(np.float64)
            ang = 2.0 * np.pi * c / box
            R = np.hypot(np.cos(ang).mean(axis=0), np.sin(ang).mean(axis=0))
            axis = int(np.argmin(R))
            if R[axis] > 0.6 or np.delete(R, axis).min() < 0.9:
                raise ValueError(f"{nm}: eje ambiguo, R = {R}")
            tr = [i for i in range(3) if i != axis]
            pos = None
            for kx, ky in _RAY_ATTR_CANDIDATES:
                if kx in g.attrs and ky in g.attrs:
                    pos = np.array([float(np.atleast_1d(g.attrs[kx])[0]),
                                    float(np.atleast_1d(g.attrs[ky])[0])])
                    break
            if pos is None:
                # respaldo: mediana de las particulas mas pegadas al rayo
                h = g["SmoothingLengths"][:]
                idx = np.argsort(h)[: max(10, len(h) // 10)]
                pos = np.median(c[idx][:, tr], axis=0)
            rays.append((nm, axis, pos))
            gattrs[nm] = dict(g.attrs)
    return rays, gattrs, z_old, kg, box


def make_uniform_rays(n_per_axis: int, box: float, seed: int):
    """
    Rayos nuevos uniformes sobre toda la cara de la caja, n por eje.

    Los atributos de grupo replican los que escribe SWIFT, no solo Xpos/Ypos:
    Zaxis es el eje de INTEGRACION y Xaxis/Yaxis son los dos transversales en
    orden creciente de indice, que son a los que se refieren Xpos e Ypos.
    spectWizard los exige y aborta sin ellos; nuestro propio extractor no los
    usa (lee Xpos/Ypos), asi que agregarlos no cambia ningun tau.

    El sorteo NO se toca: rng.uniform sigue siendo la unica llamada al
    generador y en el mismo orden, asi que --uniform N --seed S devuelve
    exactamente las mismas posiciones que antes de este cambio.
    """
    rng = np.random.default_rng(seed)
    rays, gattrs, i = [], {}, 0
    for axis in range(3):
        tr = [j for j in range(3) if j != axis]
        for _ in range(n_per_axis):
            nm = f"LOS_{i:04d}"
            pos = rng.uniform(0.0, box, size=2)
            rays.append((nm, axis, pos))
            gattrs[nm] = {"Xaxis": np.array([tr[0]], dtype=np.int32),
                          "Xpos": np.array([pos[0]]),
                          "Yaxis": np.array([tr[1]], dtype=np.int32),
                          "Ypos": np.array([pos[1]]),
                          "Zaxis": np.array([axis], dtype=np.int32)}
            i += 1
    return rays, gattrs


# ---------------------------------------------------------------------------
# Seleccion desde el snapshot
# ---------------------------------------------------------------------------

def _is_virtual(path):
    """
    True si el archivo es el .hdf5 VIRTUAL que SWIFT escribe con
    distributed:1 para agregar los pedazos .0.hdf5 .. .N.hdf5.

    Incluirlo junto con los pedazos hace que cada particula se lea DOS veces:
    la densidad SPH sale al doble y tau tambien. El sintoma es un factor de
    reescalado A ~ 0.5 y wsum ~ 1.7 en vez de ~0.85.

    Se detecta por el atributo Header/Virtual que escribe SWIFT; como
    respaldo, por NumPart_ThisFile == NumPart_Total con
    NumFilesPerSnapshot > 1.
    """
    try:
        with h5py.File(path, "r") as f:
            a = f["Header"].attrs
            if "Virtual" in a and int(np.atleast_1d(a["Virtual"])[0]) == 1:
                return True
            nf = int(np.atleast_1d(a.get("NumFilesPerSnapshot", [1]))[0])
            if nf > 1 and "NumPart_ThisFile" in a and "NumPart_Total" in a:
                this = np.atleast_1d(a["NumPart_ThisFile"]).astype(np.int64)
                tot = np.atleast_1d(a["NumPart_Total"]).astype(np.int64)
                if this.sum() == tot.sum() and tot.sum() > 0:
                    return True
    except Exception:
        pass
    return False


def resolve_snapshot_files(pattern):
    """Expande el patron y descarta el archivo virtual si hay pedazos."""
    files = sorted(glob.glob(pattern))
    if not files:
        return [], []
    real = [p for p in files if not _is_virtual(p)]
    dropped = [p for p in files if p not in real]
    return (real, dropped) if real else (files, [])


def select_particles(snap_files, rays, box, kernel_gamma, fields,
                     chunk=4_000_000, verbose=True):
    """
    Recorre los snapshots en bloques y devuelve, por rayo, los datos de sus
    particulas.

    Estrategia: el arbol se construye sobre las PARTICULAS del bloque (2D,
    periodico, en los dos ejes transversales) y se consulta una vez por RAYO.
    Asi el bucle a nivel Python corre sobre rayos (cientos) y no sobre
    particulas (cientos de millones), y el test exacto b < gamma*h_i queda
    vectorizado con numpy. La version ingenua (arbol sobre rayos, bucle sobre
    particulas) es ~5 ordenes de magnitud mas lenta a esta escala.
    """
    by_axis = {}
    for i, r in enumerate(rays):
        by_axis.setdefault(r[1], []).append(i)

    store = {i: {f: [] for f in fields} for i in range(len(rays))}
    npart_total = 0

    for path in snap_files:
        with h5py.File(path, "r") as f:
            if "PartType0" not in f:
                continue
            p0 = f["PartType0"]
            n = p0["Coordinates"].shape[0]
            npart_total += n
            for s in range(0, n, chunk):
                e = min(s + chunk, n)
                coords = p0["Coordinates"][s:e].astype(np.float64)
                Hsup = p0["SmoothingLengths"][s:e].astype(np.float64) * kernel_gamma
                rmax = float(Hsup.max())
                data = None  # se lee solo si hay candidatos

                for axis, ids in by_axis.items():
                    tr = [j for j in range(3) if j != axis]
                    pts = np.mod(coords[:, tr], box)
                    tree = cKDTree(pts, boxsize=[box, box])
                    ray_pos = np.array([rays[i][2] for i in ids])
                    cand = tree.query_ball_point(np.mod(ray_pos, box),
                                                 r=rmax, workers=-1)
                    for k, i in enumerate(ids):
                        c = np.asarray(cand[k], dtype=np.int64)
                        if c.size == 0:
                            continue
                        d = pts[c] - rays[i][2]
                        d -= box * np.round(d / box)
                        keep = c[(d[:, 0] ** 2 + d[:, 1] ** 2) <= Hsup[c] ** 2]
                        if keep.size == 0:
                            continue
                        if data is None:
                            data = {fl: p0[fl][s:e] for fl in fields}
                        for fl in fields:
                            store[i][fl].append(data[fl][keep])
                if verbose:
                    print(f"    {path}: {e}/{n}", flush=True)

    # concatenar los bloques por rayo
    ndup = 0
    for i in store:
        for fl in fields:
            store[i][fl] = (np.concatenate(store[i][fl]) if store[i][fl]
                            else np.empty((0,)))
        # red de seguridad: si por cualquier via una particula entro dos veces
        # (archivo virtual junto con sus pedazos, glob solapado), deduplicar.
        if "ParticleIDs" in fields and store[i]["ParticleIDs"].size:
            _, keep = np.unique(store[i]["ParticleIDs"], return_index=True)
            if keep.size != store[i]["ParticleIDs"].size:
                ndup += store[i]["ParticleIDs"].size - keep.size
                keep = np.sort(keep)
                for fl in fields:
                    store[i][fl] = store[i][fl][keep]
    if ndup:
        print(f"  AVISO: se eliminaron {ndup} particulas duplicadas "
              "(revisa el patron de snapshot).")
    return store, npart_total


# ---------------------------------------------------------------------------
# Escritura con el mismo esquema
# ---------------------------------------------------------------------------

def write_output(out_path, snap0_path, old_los_path, rays, gattrs, store,
                 fields):
    with h5py.File(out_path, "w") as out, \
         h5py.File(snap0_path, "r") as snap, \
         h5py.File(old_los_path, "r") as old:

        # grupos raiz: copiar atributos del archivo de LOS viejo (Header con
        # OutputType etc.) y completar con los del snapshot si faltan.
        for grp in ["Header", "Cosmology", "Units", "HydroScheme",
                    "InternalCodeUnits", "GravityScheme", "ICs_parameters",
                    "Code"]:
            src = old if grp in old else (snap if grp in snap else None)
            if src is None:
                continue
            g = out.create_group(grp)
            for k, v in src[grp].attrs.items():
                g.attrs[k] = v

        counts = []
        p0 = snap["PartType0"]
        for i, (nm, axis, pos) in enumerate(rays):
            g = out.create_group(nm)
            for k, v in gattrs.get(nm, {}).items():
                g.attrs[k] = v
            n_i = len(store[i][fields[0]])
            counts.append(n_i)
            # NumParts describe ESTE archivo. Por el camino que reusa los
            # rayos viejos se copia del LOS original y queda con el conteo
            # del original, que no es el que acabamos de escribir; por el
            # camino --uniform directamente no existe. Se fija siempre.
            g.attrs["NumParts"] = np.array([n_i], dtype=np.int32)
            # Y si el LOS viejo no traia los ejes, se completan: sin ellos
            # spectWizard no sabe a lo largo de que eje corre la linea.
            if "Zaxis" not in g.attrs:
                tr = [j for j in range(3) if j != axis]
                g.attrs["Xaxis"] = np.array([tr[0]], dtype=np.int32)
                g.attrs["Yaxis"] = np.array([tr[1]], dtype=np.int32)
                g.attrs["Zaxis"] = np.array([axis], dtype=np.int32)
            for fld in fields:
                arr = np.asarray(store[i][fld])
                d = g.create_dataset(fld, data=arr)
                for k, v in p0[fld].attrs.items():
                    d.attrs[k] = v

        tot = int(np.sum(counts))
        npt = np.zeros(7, dtype=np.int64)
        npt[0] = tot
        out["Header"].attrs["NumPart_ThisFile"] = npt
        out["Header"].attrs["NumPart_Total"] = npt.astype(np.uint32)
        out["Header"].attrs["TotalNumberOfParticles"] = npt
    return counts


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-los", required=True,
                    help="archivo de LOS original (fuente de rayos y esquema)")
    ap.add_argument("--snapshot", required=True,
                    help="snapshot al mismo z; acepta glob para distribuidos")
    ap.add_argument("--out", required=True)
    ap.add_argument("--uniform", type=int, default=None,
                    help="ignorar los rayos viejos y tirar N nuevos por eje, "
                         "uniformes sobre TODA la caja")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-los", type=int, default=None,
                    help="limitar rayos (para pruebas o tandas)")
    ap.add_argument("--skip-los", type=int, default=0,
                    help="saltear los primeros N rayos (para correr por "
                         "tandas: --skip-los 0/2048/4096 con --max-los 2048)")
    ap.add_argument("--chunk", type=int, default=4_000_000)
    args = ap.parse_args()

    snap_files, dropped = resolve_snapshot_files(args.snapshot)
    if not snap_files:
        sys.exit(f"No encontre snapshots con el patron {args.snapshot}")
    if dropped:
        print("Archivos VIRTUALES descartados (evitan doble conteo):")
        for p in dropped:
            print(f"  {p}")

    rays, gattrs, z_old, kg, box = read_rays(args.old_los)
    if args.uniform:
        rays, gattrs = make_uniform_rays(args.uniform, box, args.seed)
        print(f"Rayos nuevos uniformes: {len(rays)}")
    if args.skip_los:
        rays = rays[args.skip_los:]
    if args.max_los:
        rays = rays[: args.max_los]

    with h5py.File(snap_files[0], "r") as f:
        z_snap = float(np.atleast_1d(f["Cosmology"].attrs["Redshift"])[0])
        box_snap = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        avail = set(f["PartType0"].keys())
    if abs(z_snap - z_old) > 0.02:
        sys.exit(f"Redshift no coincide: LOS z={z_old}, snapshot z={z_snap}")
    if abs(box_snap - box) > 1e-6 * box:
        sys.exit(f"BoxSize no coincide: {box} vs {box_snap}")

    with h5py.File(args.old_los, "r") as f:
        first = sorted(k for k in f.keys() if k.startswith("LOS_"))[0]
        wanted = list(f[first].keys())
    fields = [x for x in wanted if x in avail]
    missing_ess = [x for x in ESSENTIAL if x not in fields]
    if missing_ess:
        sys.exit(f"Al snapshot le faltan campos esenciales: {missing_ess}")
    skipped = [x for x in wanted if x not in avail]
    if skipped:
        print(f"Campos del LOS viejo ausentes en el snapshot (se omiten): "
              f"{skipped}")

    print(f"\nCaja = {box:.4f} internas, z = {z_snap:.4f}, "
          f"gamma = {kg:.6f}, rayos = {len(rays)}, "
          f"snapshots = {len(snap_files)} archivo(s)")
    print(f"Campos a copiar: {fields}\n")

    # --- validacion de los datos de entrada -------------------------------
    print("Validando el snapshot (primeras 5e6 particulas por archivo):")
    bad = False
    for p in snap_files:
        with h5py.File(p, "r") as f:
            p0 = f["PartType0"]
            n = min(p0["Coordinates"].shape[0], 5_000_000)
            for fl, lim in [("Masses", 0.0), ("Densities", 0.0),
                            ("Temperatures", 0.0)]:
                if fl not in p0:
                    continue
                v = p0[fl][:n]
                nneg = int((v <= lim).sum())
                if nneg:
                    bad = True
                    print(f"  {os.path.basename(p)} {fl}: {nneg} valores "
                          f"<= {lim} (min = {v.min():.4e})")
    if bad:
        print("  ATENCION: hay masas/densidades/temperaturas no positivas.")
        print("  Eso produce tau NEGATIVO aguas abajo y rompe el reescalado.")
        print("  Es el problema conocido de masas de bariones negativas de las")
        print("  IC de monofonIC con el boost FCT. Hay que resolverlo en las IC,")
        print("  no aca.")
    else:
        print("  OK: masas, densidades y temperaturas positivas.\n")

    store, ntot = select_particles(snap_files, rays, box, kg, fields,
                                   chunk=args.chunk)
    counts = write_output(args.out, snap_files[0], args.old_los, rays,
                          gattrs, store, fields)

    print(f"\nListo -> {args.out}")
    print(f"  particulas del snapshot recorridas: {ntot}")
    print(f"  particulas por LOS: mediana = {int(np.median(counts))}, "
          f"min = {min(counts)}, max = {max(counts)}")
    print("  Verifica con: python diagnose_extraction.py (pixeles vacios "
          "debe dar ~0%)")


if __name__ == "__main__":
    main()
