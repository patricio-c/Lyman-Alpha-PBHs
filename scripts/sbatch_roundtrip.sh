#!/bin/bash
#SBATCH --job-name=lya_roundtrip
#SBATCH --output=logs/roundtrip_%j.out
#SBATCH --error=logs/roundtrip_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#
# V1: does our LOS regeneration reproduce SWIFT's own output?
#
# Usage:
#   sbatch scripts/sbatch_roundtrip.sh <RUN_DIR> <Z> [NLOS]
#
# e.g.
#   sbatch scripts/sbatch_roundtrip.sh \
#       /data/contrib/pad_140/pcolazo/lyman/murgia/cdm 5.6 100
#
# Two arguments, both plain: a run DIRECTORY and a redshift. No globs, so
# nothing for the shell to expand behind your back. legacy/relosz.py does
# the resolving - it finds the LOS file and the snapshot by the redshift
# each file reports, handles snapshots kept in their own subdirectory,
# builds the pieces pattern, and calls relos.py. That is the path Pato
# normally uses and it is the one that knows about the subdirectory layout.
#
# relosz.py looks for snapshots as DIRECTORIES. A run that keeps its
# snapshots as plain files next to the LOS files (the 40 Mpc/h pair does)
# will make it report "No hay snapshot a z=..."; for those, call relos.py
# directly with an explicit --snapshot glob instead.
#
# Choosing the input. The reference has to be a LOS file SWIFT wrote on the
# fly and did NOT truncate. The 40 Mpc/h pair only holds data out to 40 of
# 58.7372 Mpc from the h mix-up, so it cannot serve: relos.py would be
# asked to reproduce a file that is already wrong. murgia is clean - 99.8%
# transverse coverage, verified 2026-09-01. Step 4 of
# scripts/run_validation_A.sh lists every (run, z) that has both a LOS file
# and a snapshot within dz <= 0.02, and step 4b prints each one's ray range,
# so a truncated candidate is visible before any queue time goes into it.
set -euo pipefail
RUNDIR="${1:?run directory, e.g. .../lyman/murgia/cdm}"
Z="${2:?redshift, e.g. 5.6}"
NLOS="${3:-100}"

source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate astro
cd "${SLURM_SUBMIT_DIR:-.}"; mkdir -p logs data figures

TAG="$(basename "$RUNDIR")_z${Z}_n${NLOS}"
REGEN="data/regen_${TAG}.hdf5"
RELOG="logs/relosz_${TAG}.log"

echo "===== V1 step 1: regenerating ${NLOS} sightlines from the snapshot ====="
echo "  run      : ${RUNDIR}"
echo "  redshift : ${Z}"
echo "  output   : ${REGEN}"
echo
python legacy/relosz.py "$RUNDIR" "$Z" --out "$REGEN" --max-los "$NLOS" 2>&1 | tee "$RELOG"

# Compare against exactly the file relos.py used as its template, taken from
# what relosz.py printed. Re-deriving it here with a second copy of the
# matching rule is how the two quietly drift apart.
TEMPLATE="$(sed -n 's/.*plantilla:[[:space:]]*//p' "$RELOG" | head -1)"
if [ -z "$TEMPLATE" ] || [ ! -f "$TEMPLATE" ]; then
    echo
    echo "Could not read the template LOS file back from relosz.py output."
    echo "Nothing to compare against; see ${RELOG}."
    exit 1
fi

echo
echo "===== V1 step 2: comparing against what SWIFT itself wrote ====="
echo "  reference: ${TEMPLATE}"
python tests/t12_relos_roundtrip.py --a "$TEMPLATE" --b "$REGEN" \
    --max-los "$NLOS" --out figures/t12_roundtrip_${TAG}.txt
