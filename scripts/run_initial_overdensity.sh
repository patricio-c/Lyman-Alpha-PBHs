#!/bin/bash
#SBATCH --job-name=stage13
#SBATCH --partition=cpunode
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=03:00:00
#SBATCH --output=logs/stage13_%j.out
#SBATCH --error=logs/stage13_%j.err
# Lanzar desde la raiz del repo:
#   cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs && mkdir -p logs && sbatch scripts/run_initial_overdensity.sh
# Por defecto el par VIEJO (el resultado 2.65 del 2026-10-27). Otro par:
#   CDM=$L/cdm-box-40-1024 FCT=$L/fct-box-40-1024 OUT=figures/qla_initial_overdensity_new40 \
#       sbatch scripts/run_initial_overdensity.sh
# EXTRA se pasa tal cual a stage 13 (p.ej. EXTRA="--z-early 15").
set -eu
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro
cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
OLD=/data/contrib/pad_140/pcolazo/lyman/Old_sim
CDM=${CDM:-$OLD/cdm-box-40-1024}
FCT=${FCT:-$OLD/2-fct-box-40-1024}
OUT=${OUT:-figures/qla_initial_overdensity}
python stages/13_initial_overdensity.py \
    --cdm "$CDM" \
    --fct "$FCT" \
    --out "$OUT" ${EXTRA:-}
