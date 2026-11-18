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

set -eu
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro

cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
OLD=/data/contrib/pad_140/pcolazo/lyman/Old_sim

python stages/13_initial_overdensity.py \
    --cdm "$OLD/cdm-box-40-1024" \
    --fct "$OLD/2-fct-box-40-1024" \
    --out figures/qla_initial_overdensity
