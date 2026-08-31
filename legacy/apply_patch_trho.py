#!/usr/bin/env python
"""
apply_patch_trho.py - Agrega el flag --impose-trho T0 GAMMA.

Reemplaza la temperatura de cada particula por T = T0 * Delta^(gamma-1).
Corriendo las DOS simulaciones con los MISMOS T0 y gamma, la diferencia
termica queda apagada y lo que sobrevive en el cociente es estructura de
densidad pura.

Requiere apply_patch_delta.py aplicado antes. Idempotente, hace backup.
Correr en el directorio LOS:

    python apply_patch_trho.py
"""
import shutil
import sys

CHANGED, OK = [], []


def patch(fname, repls, marker):
    src = open(fname).read()
    if marker in src:
        OK.append(f"{fname}: ya parcheado")
        return
    out = src
    for old, new in repls:
        if out.count(old) != 1:
            sys.exit(f"ANCLA no unica ({out.count(old)}) en {fname}:\n{old[:70]}\n"
                     "Aplicaste apply_patch_delta.py antes?")
        out = out.replace(old, new)
    shutil.copy(fname, fname + ".bak_trho")
    open(fname, "w").write(out)
    CHANGED.append(fname)


# --------------------------- swift_extract.py ------------------------------
SIG_OLD = """                normalize: bool = True, w_floor: float = 0.05,
                delta_max: float | None = None):"""
SIG_NEW = """                normalize: bool = True, w_floor: float = 0.05,
                delta_max: float | None = None,
                impose_trho: tuple | None = None):"""

T_OLD = "    v_par = vel[:, axis]                                   # km/s peculiar fisica"
T_NEW = """    if impose_trho is not None:
        # T = T0 * Delta^(gamma-1) con Delta = rho / rho_b_media(z).
        # rho ya esta en cgs FISICAS aca.
        rho_b = 1.87847e-29 * meta.h ** 2 * meta.omega_b * (1.0 + meta.z) ** 3
        T = float(impose_trho[0]) * (rho / rho_b) ** (float(impose_trho[1]) - 1.0)

    v_par = vel[:, axis]                                   # km/s peculiar fisica"""

patch("swift_extract.py", [(SIG_OLD, SIG_NEW), (T_OLD, T_NEW)], "impose_trho")

# ---------------------------- forest_tools.py ------------------------------
CALL_OLD = """                                    delta_max=args.delta_max,
                                    collect_diag=True)"""
CALL_NEW = """                                    delta_max=args.delta_max,
                                    impose_trho=args.impose_trho,
                                    collect_diag=True)"""

ARG_OLD = '    c.add_argument("--out", required=True)\n    c.set_defaults(func=cmd_cache)'
ARG_NEW = ('    c.add_argument("--impose-trho", type=float, nargs=2,\n'
           '                   metavar=("T0", "GAMMA"), default=None,\n'
           '                   help="reemplazar T por T0*Delta^(gamma-1). Con los '
           'mismos valores en ambas corridas apaga la diferencia termica.")\n'
           '    c.add_argument("--out", required=True)\n'
           '    c.set_defaults(func=cmd_cache)')

patch("forest_tools.py", [(CALL_OLD, CALL_NEW), (ARG_OLD, ARG_NEW)], "impose_trho")

import py_compile
for f in ("swift_extract.py", "forest_tools.py"):
    py_compile.compile(f, doraise=True)
for m in OK:
    print(m)
for f in CHANGED:
    print(f"parcheado -> {f}  (backup en {f}.bak_trho)")
print("compila OK; probar:  python forest_tools.py cache -h | grep trho")
