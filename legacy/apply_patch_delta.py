#!/usr/bin/env python
"""
apply_patch_delta.py - Aplica el flag --delta-max a swift_extract.py y
forest_tools.py. Idempotente: si ya esta aplicado, no toca nada.
Hace backup .bak_delta antes de escribir. Correr en el directorio LOS:

    python apply_patch_delta.py
"""
import shutil
import sys

OK, CHANGED = [], []


def patch(fname, repls):
    src = open(fname).read()
    if "delta_max" in src:
        OK.append(f"{fname}: ya parcheado, no se toca")
        return
    out = src
    for old, new in repls:
        if old not in out:
            sys.exit(f"ANCLA NO ENCONTRADA en {fname}:\n---\n{old}\n---\n"
                     "El archivo difiere de la version auditada; parchear a "
                     "mano con patch_delta_cut.md.")
        if out.count(old) != 1:
            sys.exit(f"ANCLA AMBIGUA ({out.count(old)} veces) en {fname}: {old[:60]}")
        out = out.replace(old, new)
    shutil.copy(fname, fname + ".bak_delta")
    open(fname, "w").write(out)
    CHANGED.append(fname)


# --------------------------- swift_extract.py ------------------------------
SIG_OLD = ("                normalize: bool = True, w_floor: float = 0.05):")
SIG_NEW = ("                normalize: bool = True, w_floor: float = 0.05,\n"
           "                delta_max: float | None = None):")

CUT_ANCHOR = "    # Geometria en cm propios."
CUT_NEW = """    frac_cut = 0.0
    if delta_max is not None:
        # rho ya esta en cgs FISICAS: Delta = rho / rho_b_media(z)
        rho_b = 1.87847e-29 * meta.h ** 2 * meta.omega_b * (1.0 + meta.z) ** 3
        keep = rho / rho_b <= delta_max
        frac_cut = float(1.0 - keep.mean())
        if keep.sum() < 10:
            raise ValueError(f"{los_name}: delta_max={delta_max} deja "
                             f"{int(keep.sum())} particulas.")
        coords_int = coords_int[keep]
        rho, T, vel = rho[keep], T[keep], vel[keep]
        hsml_int, mass = hsml_int[keep], mass[keep]

    # Geometria en cm propios."""

DIAG_OLD = '        "npart": len(x_par), "band_halfwidth": hw,'
DIAG_NEW = ('        "npart": len(x_par), "band_halfwidth": hw,\n'
            '        "frac_delta_cut": frac_cut,')

patch("swift_extract.py", [(SIG_OLD, SIG_NEW),
                           (CUT_ANCHOR, CUT_NEW),
                           (DIAG_OLD, DIAG_NEW)])

# ---------------------------- forest_tools.py ------------------------------
CALL_OLD = """    tau, dv, _, diags = extract_all(args.los, args.npix, G,
                                    max_los=args.max_los,
                                    normalize=not args.no_normalize,
                                    collect_diag=True)"""
CALL_NEW = """    tau, dv, _, diags = extract_all(args.los, args.npix, G,
                                    max_los=args.max_los,
                                    normalize=not args.no_normalize,
                                    delta_max=args.delta_max,
                                    collect_diag=True)"""

SAVE_OLD = "        normalized=not args.no_normalize, wsum_raw_med=wr)"
SAVE_NEW = ("        normalized=not args.no_normalize, wsum_raw_med=wr,\n"
            "        delta_max=args.delta_max if args.delta_max is not None"
            " else 0.0)")

ARG_OLD = '    c.add_argument("--out", required=True)\n    c.set_defaults(func=cmd_cache)'
ARG_NEW = ('    c.add_argument("--delta-max", type=float, default=None,\n'
           '                   help="descartar gas con Delta = rho/rho_b por '
           'encima de este valor antes de depositar (test del esquema QLA).")\n'
           '    c.add_argument("--out", required=True)\n'
           '    c.set_defaults(func=cmd_cache)')

patch("forest_tools.py", [(CALL_OLD, CALL_NEW),
                          (SAVE_OLD, SAVE_NEW),
                          (ARG_OLD, ARG_NEW)])

# --------------------------------- checks ----------------------------------
import py_compile
for f in ("swift_extract.py", "forest_tools.py"):
    py_compile.compile(f, doraise=True)
for m in OK:
    print(m)
for f in CHANGED:
    print(f"parcheado -> {f}  (backup en {f}.bak_delta)")
print("compilacion OK; probar:  python forest_tools.py cache -h | grep delta")
