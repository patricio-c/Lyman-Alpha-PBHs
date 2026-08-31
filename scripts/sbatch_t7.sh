#!/bin/bash
#SBATCH --job-name=lya_t7
#SBATCH --output=logs/t7_%j.out
#SBATCH --error=logs/t7_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=24:00:00
# Usage: sbatch scripts/sbatch_t7.sh fct40 2.0e-3
# Builds the augmented snapshots (T0 sweep) and re-extracts each one.
set -euo pipefail
RUN="${1:?run name}"; MMAX="${2:?--conv-mass-max from stage 00 --deep}"
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate astro
cd "${SLURM_SUBMIT_DIR:-.}"; mkdir -p logs cache data

python tests/t7_stars_back.py --run "$RUN" --z 3.0 \
    --conv-parttype 1 --conv-mass-max "$MMAX" \
    --t0-sweep 0.5 1.0 2.0 --out "data/aug_${RUN}.hdf5"

for F in 0.5 1 2; do
  # NOTE: t7 writes an augmented SNAPSHOT. SWIFT has to re-shoot the LOS
  # through it before stage 01 can run. See README, t7 section.
  python stages/01_extract_los.py --los-file "data/aug_${RUN}_T0x${F}_los.hdf5" \
      --npix 2048 --treecool data/TREECOOL_HM12_G+Q \
      --out "cache/cache_${RUN}_t7_T0x${F}.npz"
done

python stages/04_p1d.py cache/cache_cdm.npz "cache/cache_${RUN}.npz" \
    cache/cache_${RUN}_t7_T0x0.5.npz cache/cache_${RUN}_t7_T0x1.npz \
    cache/cache_${RUN}_t7_T0x2.npz \
    --labels CDM FCT "FCT+gas T0x0.5" "FCT+gas T0x1" "FCT+gas T0x2" \
    --out figures/t7
