#!/bin/bash
#SBATCH --job-name=lya_roundtrip
#SBATCH --output=logs/roundtrip_%j.out
#SBATCH --error=logs/roundtrip_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#
# V1: does legacy/relos.py reproduce SWIFT's own LOS output?
#
# Usage:
#   sbatch scripts/sbatch_roundtrip.sh <NATIVE_LOS.hdf5> <'SNAPSHOT*.hdf5'> [NLOS]
#
# NATIVE_LOS must be a file SWIFT wrote ON THE FLY. Not one relos.py
# produced - that compares relos.py against itself - and not one written
# with a truncated range_when_shooting_down_*, which fails by construction
# because a third of every sightline is missing. Check the candidate first:
#
#   python tests/t12_relos_roundtrip.py --a <candidate>
#
# and read the ray-position range it prints. If the rays do not reach the
# far side of the box, that file is truncated: pick another run.
#
# The snapshot glob must be quoted so the shell does not expand it, and it
# must be at the SAME redshift as the LOS file - relos.py refuses otherwise.
# Quote it: 'cdm-40-m6-lyman_0003*.hdf5'.
#
# This goes through the queue rather than an interactive session because
# relos.py reads the entire snapshot, tens of GB, however few rays are
# asked for. The cost is the read, not the geometry.
set -euo pipefail
LOS="${1:?native SWIFT LOS file}"
SNAP="${2:?snapshot glob, quoted}"
NLOS="${3:-100}"

source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate astro
cd "${SLURM_SUBMIT_DIR:-.}"; mkdir -p logs data figures

REGEN="data/regen_$(basename "${LOS%.hdf5}")_n${NLOS}.hdf5"

echo "===== V1 step 1: regenerating ${NLOS} sightlines from the snapshot ====="
echo "  reference : ${LOS}"
echo "  snapshot  : ${SNAP}"
echo "  output    : ${REGEN}"
echo
python legacy/relos.py --old-los "${LOS}" --snapshot "${SNAP}" \
    --out "${REGEN}" --max-los "${NLOS}"

echo
echo "===== V1 step 2: comparing against what SWIFT itself wrote ====="
python tests/t12_relos_roundtrip.py --a "${LOS}" --b "${REGEN}" \
    --max-los "${NLOS}" --out figures/t12_roundtrip.txt
