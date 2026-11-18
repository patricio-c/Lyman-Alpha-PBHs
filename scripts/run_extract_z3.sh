#!/bin/bash
#SBATCH --job-name=lya_z3
#SBATCH --partition=cpunode
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --output=logs/extract_z3_%j.out
#SBATCH --error=logs/extract_z3_%j.err

# Uso (parado en el repo):
#   cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
#   mkdir -p logs && sbatch scripts/run_extract_z3.sh              # test: 768 sightlines por corrida
#   mkdir -p logs && MAXLOS=0 sbatch scripts/run_extract_z3.sh     # completo: 6144 sightlines
# Si un cache ya existe lo saltea, asi que se puede relanzar sin repetir trabajo.

set -u
source /home/pcolazo/miniconda3/etc/profile.d/conda.sh
conda activate astro
# muchas extracciones en paralelo en un nodo: un hilo cada una
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

REPO=/data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs
LR=/data/contrib/pad_140/pcolazo/lyman/low_resolution
cd "$REPO"

MAXLOS=${MAXLOS:-768}
NPAR=${NPAR:-22}
if [ "$MAXLOS" -gt 0 ]; then TAG="_n${MAXLOS}"; LIM="--max-los $MAXLOS"; else TAG=""; LIM=""; fi
mkdir -p cache logs figures

# label  run (nombre de registro o directorio)
RUNS=(
  "NB_1 NB_1" "NB_2 NB_2" "NB_3 NB_3" "NB_4 NB_4" "NB_5 NB_5" "NB_6 NB_6"
  "poisson_1 poisson_1" "poisson_2 poisson_2" "poisson_3 poisson_3"
  "poisson_4 poisson_4" "poisson_5 poisson_5" "poisson_6 poisson_6"
  "a_cdm $LR/a/CDM"                  "a_fct $LR/a/FCT"
  "delta_cdm $LR/c/delta_dif/CDM"    "delta_fct $LR/c/delta_dif/FCT"
  "zreio6_cdm $LR/c/z_reio/6/CDM"    "zreio6_fct $LR/c/z_reio/6/FCT"
  "zreio75_cdm $LR/c/z_reio/7.5/CDM" "zreio75_fct $LR/c/z_reio/7.5/FCT"
  "zreio9_cdm $LR/c/z_reio/9/CDM"    "zreio9_fct $LR/c/z_reio/9/FCT"
)

run_one() {
  local label=$1 run=$2
  local out=cache/cache_${label}${TAG}.npz
  if [ -f "$out" ]; then echo "skip   $label (ya existe $out)"; return; fi
  python stages/01_extract_los.py --run "$run" --z 3.0 --npix 2048 \
      --treecool data/TREECOOL_HM12_G+Q $LIM --progress-every 200 \
      --out "$out" > "logs/extract_${label}${TAG}.log" 2>&1 \
    && echo "done   $label" || echo "FAILED $label  (ver logs/extract_${label}${TAG}.log)"
}

echo "inicio $(date)   MAXLOS=$MAXLOS  NPAR=$NPAR  nodo=$(hostname)"
for r in "${RUNS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do wait -n; done
  run_one $r &
done
wait
echo "extraccion terminada $(date)"

# ---- plots (stage 04, normalizacion cruda: --norm none) ----
plot() {  # plot <out> <label_ref> <cache_ref> <label> ...
  local out=$1; shift
  local files=() labels=()
  while [ $# -gt 0 ]; do
    local lab=$1 f=$2; shift 2
    [ -f "$f" ] && { files+=("$f"); labels+=("$lab"); }
  done
  [ ${#files[@]} -ge 2 ] || { echo "plot $out: faltan caches"; return; }
  python stages/04_p1d.py "${files[@]}" --labels "${labels[@]}" \
      --norm none --ylim-ratio 0.5 2.0 --out "$out"
}

# 1) los 12 more_power contra el CDM b40 de produccion (cache completo ya hecho)
args=(CDM_b40 cache/cache_cdm40.npz)
for m in NB_1 NB_2 NB_3 NB_4 NB_5 NB_6 poisson_1 poisson_2 poisson_3 poisson_4 poisson_5 poisson_6; do
  args+=("$m" "cache/cache_${m}${TAG}.npz")
done
plot "figures/p1d_z3_more_power${TAG}" "${args[@]}"

# 2) cada leg nueva contra SU propio CDM (misma caja y resolucion)
for leg in a delta zreio6 zreio75 zreio9; do
  plot "figures/p1d_z3_${leg}${TAG}" \
       "${leg}_CDM" "cache/cache_${leg}_cdm${TAG}.npz" \
       "${leg}_FCT" "cache/cache_${leg}_fct${TAG}.npz"
done

echo "todo listo $(date)"
