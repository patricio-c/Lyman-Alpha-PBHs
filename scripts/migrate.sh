#!/usr/bin/env bash
# Bring the working directory on Clementina into the repo.
#
# This COPIES. It deletes nothing. Delete the old directory yourself, by
# hand, only after `scripts/verify.sh` passes and you have pushed.
#
#   bash scripts/migrate.sh /data/contrib/pad_140/pcolazo/LOS
#
set -euo pipefail
SRC="${1:?usage: migrate.sh <old LOS directory>}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "source : $SRC"
echo "repo   : $ROOT"

mkdir -p "$ROOT/legacy" "$ROOT/cache" "$ROOT/data" "$ROOT/figures"

# 1. the extraction and physics code the results actually came from.
#    These stay untouched in legacy/ and are imported, never rewritten.
for f in swift_extract.py forest_tools.py sherwood_los.py ionization.py \
         relos.py relosz.py fix_los_attrs.py plot_p1d.py p1d_three_v2.py \
         gas_vs_baryon_pk.py removal_curve_pk.py mass_match_pk.py \
         saturation_test.py cross_flux.py los_paneles.py where_holes.py \
         diagnose_extraction.py deficit_from_caches.py debug_tau.py \
         compare_gas.py cociente_cdm_fct.py check_versions.py scan_runs.py \
         apply_patch_delta.py apply_patch_trho.py pipeline_sherwood.py \
         sherwood_postprocess.py plot_3d_1d.py; do
  [ -f "$SRC/$f" ] && cp -v "$SRC/$f" "$ROOT/legacy/" || true
done

# 2. the .bak files are the pre-patch versions. Keep them: they are the
#    only record of what the code looked like before --delta-max existed.
mkdir -p "$ROOT/legacy/prepatch"
for f in "$SRC"/*.bak_delta "$SRC"/*.bak_trho; do
  [ -e "$f" ] && cp -v "$f" "$ROOT/legacy/prepatch/" || true
done

# 3. caches. These are the expensive part - hours of extraction each.
for f in "$SRC"/cache_*.npz; do
  [ -e "$f" ] && cp -v "$f" "$ROOT/cache/" || true
done

# 4. external tables the physics depends on.
for f in TREECOOL_HM12_G+Q tauH1_lya_z3.0.dat; do
  [ -f "$SRC/$f" ] && cp -v "$SRC/$f" "$ROOT/data/" || true
done

# 5. figures from the talk, for the record. Not tracked by git (see
#    .gitignore) but kept locally so nothing is lost.
for f in "$SRC"/*.png; do
  [ -e "$f" ] && cp "$f" "$ROOT/figures/" || true
done

echo
echo "copied. Nothing in $SRC was modified or removed."
echo "Next:  bash scripts/verify.sh"
