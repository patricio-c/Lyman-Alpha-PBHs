#!/bin/bash
#SBATCH --job-name=ratios
#SBATCH --partition=cpunode
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=logs/ratios_%j.out
#SBATCH --error=logs/ratios_%j.err

# FCT/CDM (o cualquier test/ref) con stage 10, a partir de los caches de
# run_extract_random.sh, normalizado al tau_eff CRUDO de la referencia.
# Paired si stage 02 dice que son las mismas lineas; si no, --unpaired.
#
#   PAIRS="cdm40:fct40 lr_zreio75_cdm:lr_zreio75_fct" Z=3.0 sbatch scripts/run_pair_ratios.sh
#
# Salida: figures/dp1d_<ref>__<test>_z<Z>_n<N>_s<SEED>.{txt,png}

set -euo pipefail
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
Z=${Z:?}; N=${N:-1500}; SEED=${SEED:-1}; PAIRS=${PAIRS:?}
S="z${Z}_n${N}_s${SEED}"

one() {
  local ref=$1 test=$2 out="${1}__${2}_${S}"
  local TE PAIR=""
  TE=$(grep "tau_eff (raw" "logs/extract_${ref}_${S}.log" | awk '{print $NF}')
  python stages/02_check_los_match.py --a "cache/los_subsets/${ref}_${S}.hdf5" \
      --b "cache/los_subsets/${test}_${S}.hdf5" > "logs/match_${out}.log" 2>&1 || true
  grep -q "matched line by line" "logs/match_${out}.log" || PAIR="--unpaired"
  echo "$out  tau_eff_ref=$TE  ${PAIR:-paired}"
  python stages/10_dp1d.py "cache/cache_${ref}_${S}.npz" "cache/cache_${test}_${S}.npz" \
      --labels "$ref" "$test" --tau-eff "$TE" $PAIR --out "figures/dp1d_${out}" \
      > "logs/dp1d_${out}.log" 2>&1
}
for p in $PAIRS; do one "${p%%:*}" "${p#*:}" & done
wait
