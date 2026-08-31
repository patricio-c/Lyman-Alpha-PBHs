#!/usr/bin/env python
"""
check_versions.py - Verifica que los modulos de la pipeline sean las versiones
corregidas y no copias viejas.

Cada correccion que fuimos haciendo dejo una marca en el codigo. Este script
busca esas marcas y te dice, archivo por archivo, si te falta alguna.

Uso:
    python check_versions.py            # en el directorio con los .py
    python check_versions.py /ruta/LOS
"""

import ast
import os
import sys

# archivo -> [(marca a buscar, que correccion representa)]
CHECKS = {
    "relos.py": [
        ("--run-dir", "resolucion por redshift (evita confundir indices)"),
        ("resolve_by_redshift", "buscador de LOS y snapshot por z"),
        ("_is_virtual", "descarte del archivo virtual (doble conteo, factor 2)"),
        ("ParticleIDs", "deduplicacion de seguridad"),
        ("Validando el snapshot", "validacion de masas/densidades/temperaturas"),
        ("source = snap if grp in snap", "metadatos desde el snapshot, no del LOS viejo"),
        ("cKDTree(pts", "arbol sobre particulas (5 ordenes mas rapido)"),
        ("--skip-los", "corridas por tandas"),
    ],
    "swift_extract.py": [
        ("b_perp", "parametro de impacto transversal en la deposicion SPH"),
        ("max_b_violation", "chequeo estructural b <= H"),
        ("concentracion circular", "deteccion de eje robusta al borde periodico"),
        ('"Masses", mass', "validacion de masas no positivas"),
        ("tau negativo", "aborto ante tau negativo"),
        ("mHI_over_mH", "deposicion conservativa de HI"),
    ],
    "sherwood_los.py": [
        ("_expected_size", "deteccion de endianness por tamano de archivo"),
        ("deconvolve_pixel", "deconvolucion de la ventana de pixel"),
        ("rescale_tau", "reescalado a tau_eff por biseccion"),
    ],
    "ionization.py": [
        ("neutral_fraction", "equilibrio de fotoionizacion"),
        ("alpha_A_HII", "recombinacion caso A"),
        ("he_state", "contribucion del helio a n_e"),
    ],
    "sherwood_postprocess.py": [
        ("convolve_lsf", "perfil instrumental"),
        ("rebin_velocity", "rebin conservativo"),
        ("continuum_bias", "sesgo de continuo"),
    ],
    "where_holes.py": [("P2 exige", "reporte de posicion de agujeros")],
    "diagnose_extraction.py": [("tau_eff(A)", "curva de reescalado")],
    "scan_runs.py": [("TRUNCADO", "deteccion de truncamiento")],
    "analyze_holes.py": [("Test 1 - histograma", "test de pileup de h")],
}

OPTIONAL = {"analyze_holes.py", "sherwood_postprocess.py", "where_holes.py"}


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    print(f"\nVerificando modulos en {os.path.abspath(root)}\n")
    problems = 0

    for fname, marks in CHECKS.items():
        path = os.path.join(root, fname)
        if not os.path.exists(path):
            tag = "opcional" if fname in OPTIONAL else "FALTA"
            print(f"  [{tag:>8}] {fname}")
            if fname not in OPTIONAL:
                problems += 1
            continue

        text = open(path, encoding="utf-8", errors="replace").read()

        # 1) compila?
        try:
            ast.parse(text)
        except SyntaxError as ex:
            print(f"  [  ROTO  ] {fname}: error de sintaxis en linea {ex.lineno}")
            print(f"             -> tipico de pegar en vim sin ':set paste'")
            problems += 1
            continue

        # 2) tiene todas las marcas?
        missing = [d for m, d in marks if m and m not in text]
        if missing:
            print(f"  [ VIEJO  ] {fname}: le faltan {len(missing)} correcciones")
            for d in missing:
                print(f"             - {d}")
            problems += 1
        else:
            print(f"  [   OK   ] {fname}")

    print()
    if problems:
        print(f"{problems} archivo(s) con problemas. Volve a bajarlos y pegalos con:")
        print("    cat > ARCHIVO.py << 'EOF'")
        print("    ...contenido...")
        print("    EOF")
        print("  (vim con autoindent rompe la indentacion de Python)")
        sys.exit(1)
    print("Todos los modulos estan en su version corregida.")


if __name__ == "__main__":
    main()
