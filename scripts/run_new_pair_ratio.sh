#!/bin/bash
#SBATCH --job-name=newpair
#SBATCH --partition=cpunode
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=logs/newpair_%j.out
#SBATCH --error=logs/newpair_%j.err

# FCT/CDM de un par (por defecto el de produccion NUEVO, LOS nativas), con N sightlines al azar,
# las mismas en los dos modelos. Lanzar desde la raiz del repo:
#   cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs && mkdir -p logs
#   Z=3.0 sbatch scripts/run_new_pair_ratio.sh
#   Z=2.2 sbatch scripts/run_new_pair_ratio.sh        # otro redshift
# Opcionales: N=1500 (sightlines), SEED=1 (otra semilla = otra muestra),
# CDM_DIR / FCT_DIR / LABEL para otro par, p.ej. la leg 1/8:
#   CDM_DIR=$L/low_resolution/c/z_reio/7.5/CDM FCT_DIR=$L/low_resolution/c/z_reio/7.5/FCT \
#   LABEL=lr75 Z=3.0 sbatch scripts/run_new_pair_ratio.sh

set -euo pipefail
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

REPO=/data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
L=/data/contrib/pad_140/pcolazo/lyman
cd "$REPO"

Z=${Z:-3.3}
N=${N:-1500}
SEED=${SEED:-1}
CDM_DIR=${CDM_DIR:-$L/cdm-box-40-1024}
FCT_DIR=${FCT_DIR:-$L/fct-box-40-1024}
LABEL=${LABEL:-new40}
TAG="${LABEL}_z${Z}_n${N}_s${SEED}"
mkdir -p cache/los_subsets figures logs

# el archivo LOS se elige leyendo el redshift del header, no por indice
los_at() {
  python -c "import sys; sys.path.insert(0, '.'); from common import runs; \
print(runs.resolve_los_file(sys.argv[1], z=float(sys.argv[2]), verbose=False))" "$1" "$Z"
}
CDM_LOS=$(los_at "$CDM_DIR")
FCT_LOS=$(los_at "$FCT_DIR")
echo "z=$Z  N=$N  seed=$SEED"
echo "CDM  $CDM_LOS"
echo "FCT  $FCT_LOS"

# 1) misma muestra al azar en los dos archivos
python scripts/make_los_subset.py --in "$CDM_LOS" \
    --out "cache/los_subsets/cdm_${TAG}.hdf5" --n "$N" --random --seed "$SEED"
python scripts/make_los_subset.py --in "$FCT_LOS" \
    --out "cache/los_subsets/fct_${TAG}.hdf5" --n "$N" --random --seed "$SEED"

# 2) tau de cada uno, en paralelo
for m in cdm fct; do
  python stages/01_extract_los.py --los-file "cache/los_subsets/${m}_${TAG}.hdf5" \
      --npix 2048 --treecool data/TREECOOL_HM12_G+Q --progress-every 300 \
      --out "cache/cache_${m}_${TAG}.npz" > "logs/extract_${m}_${TAG}.log" 2>&1 &
done
wait
grep -H "tau_eff (raw" logs/extract_{cdm,fct}_${TAG}.log

# 3) son las mismas lineas en los dos? si no, bootstrap sin emparejar
python stages/02_check_los_match.py --a "cache/los_subsets/cdm_${TAG}.hdf5" \
    --b "cache/los_subsets/fct_${TAG}.hdf5" | tee "logs/match_${TAG}.log"
PAIR=""
grep -q "matched line by line" "logs/match_${TAG}.log" || PAIR="--unpaired"

# 4) cociente y Delta P1D, normalizado al tau_eff crudo del propio CDM
TE=$(grep "tau_eff (raw" "logs/extract_cdm_${TAG}.log" | awk '{print $NF}')
echo "common tau_eff target = $TE (CDM raw)   ${PAIR:-paired}"
python stages/10_dp1d.py "cache/cache_cdm_${TAG}.npz" "cache/cache_fct_${TAG}.npz" \
    --labels CDM FCT --tau-eff "$TE" $PAIR --out "figures/dp1d_${TAG}"

echo "listo: figures/dp1d_${TAG}.png  y  .txt"
