#!/usr/bin/env bash
#
# Validation block A - everything cheap, in one log.
#
#   bash scripts/run_validation_A.sh 2>&1 | tee ~/validacion_A.log
#
# Deliberately NO `set -e`. A step that fails should not take the rest of
# the block with it: the point is that one run answers several questions,
# and a failure is itself an answer.
#
# What it does not do: t8_single_los.py. Stage 02 showed the sightlines of
# the two runs are unrelated, so t8 compares unrelated lines and would only
# add an invalid figure to the log. It comes back after the LOS are
# re-shot with a shared position list.
#
# The long job (V1, the relos.py round trip) is NOT here - it reads a whole
# snapshot and belongs in the queue. See scripts/sbatch_roundtrip.sh, and
# pick its input from what section 4 below prints.

cd /data/contrib/pad_140/pcolazo/lya-repro || exit 1
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate astro

LYA=/data/contrib/pad_140/pcolazo/lyman
CDM_LOS=$LYA/cdm-box-40-1024/los_0010.hdf5
FCT_LOS=$LYA/2-fct-box-40-1024/los_0010.hdf5
mkdir -p figures cache data logs

b(){ echo; echo "################################################################"; echo "## $*"; echo "################################################################"; }

b "0. repo, environment, data"
git pull --ff-only
git log --oneline -6
python -c "import numpy,scipy,h5py,matplotlib as m;print('numpy',numpy.__version__,'| scipy',scipy.__version__,'| h5py',h5py.__version__,'| mpl',m.__version__)"
ls -la data/ 2>/dev/null | head -15
TC=$(ls data/TREECOOL* 2>/dev/null | head -1); echo "TREECOOL picked: ${TC:-NONE}"

b "1. estimator unit tests (no data, one second)"
python tests/test_estimator.py

b "2. geometry at z=3 (does not read a single particle)"
python stages/01_extract_los.py --run fct40 --z 3.0 --geometry-only

b "3. WHY the pairing failed: transverse ray range, CDM"
python tests/t12_relos_roundtrip.py --a "$CDM_LOS"
b "3b. same for FCT"
python tests/t12_relos_roundtrip.py --a "$FCT_LOS"

b "4. which runs can serve V1: a native LOS file AND a snapshot at the same z"
python - <<'PY'
import glob, os, sys
sys.path.insert(0, '.')
from common import runs
BASE = runs.BASE
dirs = []
for pat in ('*', 'murgia/*', 'more_power/*'):
    dirs += [d for d in sorted(glob.glob(os.path.join(BASE, pat)))
             if os.path.isdir(d)]
print(f"{'run':<34s} {'#LOS':>5s} {'#snap':>6s}  z with BOTH")
for d in dirs:
    try:
        files = runs.list_hdf5(d)
    except Exception as e:
        print(f"{os.path.relpath(d, BASE):<34s}  ERROR {e}")
        continue
    zl, zs = {}, {}
    for p in files:
        try:
            if runs.is_los_file(p):
                z = runs.file_redshift(p)
                if z is not None:
                    zl[round(z, 2)] = p
            else:
                z = runs.snapshot_redshift(p)
                if z is not None:
                    zs[round(z, 2)] = p
        except Exception:
            pass
    both = sorted(set(zl) & set(zs))
    print(f"{os.path.relpath(d, BASE):<34s} {len(zl):>5d} {len(zs):>6d}  {both}")
    for z in both:
        print(f"{'':<34s}   z={z:<6g} LOS {os.path.basename(zl[z])}"
              f"   SNAP {os.path.basename(zs[z])}")
PY

b "5. V3 - t0, Matteo's test: the rescaling does not manufacture the effect"
python tests/t0_rescaling.py --cdm cache/cache_cdm.npz --fct cache/cache_fct.npz --out figures/t0_rescaling

b "6. V3 - P1D and its ratio"
python stages/04_p1d.py cache/cache_cdm.npz cache/cache_fct.npz --labels CDM FCT --out figures/p1d

b "7. V3 - correlation function (a check on the binning, not on the field)"
python stages/05_xi.py cache/cache_cdm.npz cache/cache_fct.npz --labels CDM FCT --out figures/xi

b "8. V3 - t9 without positions: the two bootstraps side by side"
python tests/t9_unpaired_significance.py --cdm cache/cache_cdm.npz --fct cache/cache_fct.npz --out figures/t9_nojk

b "9. V3 - t9 with positions: adds the block jackknife"
echo "(if it aborts saying the positions do not match the cache, THAT is the"
echo " answer to 6144 vs 1536 - step 8 already produced the numbers)"
python tests/t9_unpaired_significance.py --cdm cache/cache_cdm.npz --fct cache/cache_fct.npz \
    --los-a "$CDM_LOS" --los-b "$FCT_LOS" --out figures/t9

b "10. V3 - t10, additive Delta P1D: amplitude or shape"
python tests/t10_delta_p1d.py --pair cache/cache_cdm.npz cache/cache_fct.npz --out figures/t10

b "11. V2 - cut 100 sightlines for Maria"
python scripts/make_los_subset.py --in "$CDM_LOS" --out data/los_cdm_z3_sample100.hdf5 --n 100

b "12. V2 - our own tau on those same 100 sightlines"
python stages/01_extract_los.py --los-file data/los_cdm_z3_sample100.hdf5 --npix 2048 \
    ${TC:+--treecool "$TC"} --out cache/cache_cdm_sample100.npz

b "END. figures and caches produced"
ls -la figures/ cache/ data/*.hdf5 2>/dev/null
