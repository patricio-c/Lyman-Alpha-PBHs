#!/bin/bash
#SBATCH --job-name=xrand
#SBATCH --partition=cpunode
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=08:00:00
#SBATCH --output=logs/xrand_%j.out
#SBATCH --error=logs/xrand_%j.err

# Muestra al azar de N sightlines (misma semilla en todas las corridas) y tau
# con stage 01, varias corridas en paralelo en un nodo (stage 01 es serial).
# Cada subset es un archivo LOS de SWIFT valido: es EXACTAMENTE lo que corre
# SpecWizard (docs/PROTOCOL_2026-11-18_validation_suite.md, seccion "Track").
#
# Lanzar desde la raiz del repo:
#   RUNS="cdm40 fct40" Z=3.0 sbatch scripts/run_extract_random.sh
#   RUNS="NB_1 NB_2 ... " Z=3.0 N=1500 SEED=1 sbatch scripts/run_extract_random.sh
# RUNS: nombres del registro (common/runs.py) o label=directorio.
#
# Salidas, con TAG = <label>_z<Z>_n<N>_s<SEED>:
#   cache/los_subsets/<TAG>.hdf5      el subset (entrada de SpecWizard)
#   cache/cache_<TAG>.npz             tau de stage 01
#   logs/extract_<TAG>.log            incluye "tau_eff (raw"
#   figures/los_manifest.tsv          una fila por subset, con sha256

set -euo pipefail
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

REPO=/data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
cd "$REPO"
Z=${Z:?Z requerido, p.ej. Z=3.0}
N=${N:-1500}
SEED=${SEED:-1}
RUNS=${RUNS:?RUNS requerido}
mkdir -p cache/los_subsets figures logs
MAN=figures/los_manifest.tsv
[[ -f $MAN ]] || printf "tag\trun_dir\tz_req\tlos_file\tz_file\tsubset\tN\tseed\tsha256\tcache\tgit\n" > $MAN
GIT=$(git rev-parse --short HEAD)

resolve() {  # imprime: run_dir los_file z_file
  python - "$1" "$Z" <<'EOF'
import sys, h5py
sys.path.insert(0, '.')
from common import runs
name, z = sys.argv[1], float(sys.argv[2])
d = runs.run_dir(name)
f = runs.resolve_los_file(d, z=z, verbose=False)
with h5py.File(f, 'r') as h:
    zf = float(h['Header'].attrs['Redshift'])
print(d, f, f"{zf:.4f}")
EOF
}

pids=()
for spec in $RUNS; do
  if [[ "$spec" == *=* ]]; then label=${spec%%=*}; src=${spec#*=}; else label=$spec; src=$spec; fi
  read -r RUN_DIR LOS ZF < <(resolve "$src")
  TAG="${label}_z${Z}_n${N}_s${SEED}"
  SUB=cache/los_subsets/${TAG}.hdf5
  echo "$TAG  <-  $LOS  (z_file=$ZF)"
  (
    python scripts/make_los_subset.py --in "$LOS" --out "$SUB" --n "$N" --random --seed "$SEED"
    python stages/01_extract_los.py --los-file "$SUB" --npix 2048 \
        --treecool data/TREECOOL_HM12_G+Q --progress-every 300 \
        --out "cache/cache_${TAG}.npz" > "logs/extract_${TAG}.log" 2>&1
    SHA=$(sha256sum "$SUB" | cut -c1-16)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$TAG" "$RUN_DIR" "$Z" "$LOS" "$ZF" \
        "$SUB" "$N" "$SEED" "$SHA" "cache/cache_${TAG}.npz" "$GIT" >> $MAN
  ) &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
grep -H "tau_eff (raw" logs/extract_*_z${Z}_n${N}_s${SEED}.log || true
exit $fail
