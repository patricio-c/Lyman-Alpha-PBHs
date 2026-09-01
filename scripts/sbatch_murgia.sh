#!/bin/bash
#SBATCH --job-name=lya_murgia
#SBATCH --output=logs/murgia_%j.out
#SBATCH --error=logs/murgia_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=10:00:00
#
# V3 external check: extract the murgia LOS caches for one model.
#
# Usage:
#   sbatch scripts/sbatch_murgia.sh <MODEL> [Z ...]
#
#   sbatch scripts/sbatch_murgia.sh cdm          # 4.6 5.0 5.4, the default
#   sbatch scripts/sbatch_murgia.sh M2
#   sbatch scripts/sbatch_murgia.sh M3
#   sbatch scripts/sbatch_murgia.sh cdm 4.4 4.8 5.2 5.6   # the other four
#
# Why these three redshifts by default: they are three of the four data bins
# in Murgia, Scelfo, Viel & Raccanelli 2019 (arXiv:1903.10509), which used
# MIKE and HIRES/KECK at z = 4.2, 4.6, 5.0, 5.4. We have no z=4.2 LOS file.
# The murgia LOS grid is z = 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6.
#
# The models, confirmed by Pato 2026-09-01: M2 is M_PBH = 10^2 M_sun and M3
# is M_PBH = 10^3 M_sun, both with f_PBH = 1. Those are exactly the two
# curves in the paper's Figure 1, so the ratio measured here has a published
# answer to hit. cdm is the matched LCDM run: all three share the Panphasia
# descriptor [Panph6,L20,(235287,445214,422255),S1,KK1025,CH-999,COLIBRE050]
# and the same Omega_cdm, so they are the same realisation with a different
# input transfer function - the sample variance cancels in the ratio by
# construction, at the level of the initial conditions.
#
# No snapshot is needed: stage 01 reads the SWIFT LOS files directly. Each
# file holds 1536 sightlines and takes about 70 minutes, so three redshifts
# is roughly three and a half hours. Submit the three models in parallel.
set -euo pipefail
MODEL="${1:?model: cdm, M2 or M3}"; shift || true
if [ "$#" -eq 0 ]; then ZS=(4.6 5.0 5.4); else ZS=("$@"); fi

source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate astro
cd "${SLURM_SUBMIT_DIR:-.}"; mkdir -p logs cache figures

RUNDIR="/data/contrib/pad_140/pcolazo/lyman/murgia/${MODEL}"
[ -d "$RUNDIR" ] || { echo "no existe $RUNDIR"; exit 1; }
TC="$(ls data/TREECOOL* 2>/dev/null | head -1)"
echo "model     ${MODEL}"
echo "run dir   ${RUNDIR}"
echo "redshifts ${ZS[*]}"
echo "TREECOOL  ${TC:-NONE}"
echo

for Z in "${ZS[@]}"; do
    OUT="cache/cache_murgia_${MODEL}_z${Z}.npz"
    if [ -f "$OUT" ]; then
        echo "===== z=${Z}: ${OUT} ya existe, salteo ====="
        continue
    fi
    echo "===== z=${Z} -> ${OUT} ====="
    # --run takes a raw directory and --z matches on the redshift each LOS
    # file REPORTS, not on the index in its name. Never match on the index:
    # los_0003 is not the epoch of snap_0003 and the offset is not even
    # consistent between runs.
    python stages/01_extract_los.py --run "$RUNDIR" --z "$Z" --npix 2048 \
        ${TC:+--treecool "$TC"} --out "$OUT"
    echo
done

echo "===== listo. caches de ${MODEL}: ====="
ls -la cache/cache_murgia_${MODEL}_z*.npz
