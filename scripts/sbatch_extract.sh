#!/bin/bash
#SBATCH --job-name=lya_extract
#SBATCH --output=logs/extract_%j.out
#SBATCH --error=logs/extract_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
# Usage:  sbatch scripts/sbatch_extract.sh fct40 cache/cache_fct.npz
# Never run this with `srun --pty`: the session dies with your ssh
# connection and takes the whole extraction with it.
set -euo pipefail
RUN="${1:?run name}"; OUT="${2:?output cache}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate astro
cd "${SLURM_SUBMIT_DIR:-.}"; mkdir -p logs cache
python stages/01_extract_los.py --run "$RUN" --z 3.0 --npix 2048 \
    --treecool data/TREECOOL_HM12_G+Q --out "$OUT"
