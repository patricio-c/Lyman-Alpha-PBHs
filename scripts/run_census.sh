#!/bin/bash
#SBATCH --job-name=census
#SBATCH --partition=cpunode
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --output=logs/census_%j.out
#SBATCH --error=logs/census_%j.err

# Stage 14 (SWIFT gas/matter spectra, every z, every pair of the suite) y
# stage 15 (potencia extra lineal + fraccion convertida por headers).
# Solo lee .txt y headers: minutos.
#   sbatch scripts/run_census.sh

set -euo pipefail
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro
cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
mkdir -p figures logs

# M1 el par nuevo | M2 umbral | M4 z_reio | resolucion | caja grande
python stages/14_swift_spectra.py --out figures/swift_spectra_suite --pairs \
    cdm40:fct40 \
    lr_zreio75_cdm:lr_zreio75_fct \
    lr_zreio6_cdm:lr_zreio6_fct \
    lr_zreio9_cdm:lr_zreio9_fct \
    lr_delta_cdm:lr_delta_fct \
    lr_a_cdm:lr_a_fct \
    lr_zreio75_cdm:lr_delta_cdm   lr_zreio75_fct:lr_delta_fct \
    lr_zreio75_cdm:lr_zreio6_cdm  lr_zreio75_cdm:lr_zreio9_cdm \
    lr_zreio75_fct:lr_zreio6_fct  lr_zreio75_fct:lr_zreio9_fct \
    cdm40:lr_zreio75_cdm          fct40:lr_zreio75_fct

# M3: la grilla contra CDM b40 (masked). NO es el ancla limpia: cuando
# more_power/CDM_sc exista, repetir con ella como referencia.
GRID=""
for f in NB poisson; do for i in 1 2 3 4 5 6; do GRID="$GRID cdm40:${f}_$i"; done; done
python stages/14_swift_spectra.py --out figures/swift_spectra_grid_vs_cdm40 --pairs $GRID

python stages/15_extra_power.py --out figures/extra_power
