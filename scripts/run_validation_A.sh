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
# snapshot and belongs in the queue. See scripts/sbatch_roundtrip.sh; step
# 4b below prints its exact arguments.

cd /data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs || exit 1
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
echo "These two files are KNOWN to be truncated to 40 of 58.7372 Mpc (the h"
echo "mix-up), so the truncation warning below is expected, not news. What"
echo "is news is whether the two are truncated the SAME way: if both cover"
echo "the same range along the same axes, truncation does NOT explain the"
echo "pairing failure and a different seed does."
python tests/t12_relos_roundtrip.py --a "$CDM_LOS"
b "3b. same for FCT"
python tests/t12_relos_roundtrip.py --a "$FCT_LOS"

b "4. which runs can serve V1: a native LOS file AND a snapshot at the same z"
echo "Snapshots live in a subdirectory in most of these runs and beside the"
echo "LOS files in others, and the indices do not correspond. So nothing is"
echo "matched by name: this reads the redshift out of every file."
python - <<'PY'
import glob, os, sys
sys.path.insert(0, '.')
from common import runs

# relos.py refuses a snapshot more than this far from the LOS redshift, so
# proposing a wider pair would only produce a job that dies on arrival.
ZTOL = 0.02

BASE = runs.BASE
rows = []
dirs = []
for pat in ('*', 'murgia/*', 'more_power/*'):
    dirs += [d for d in sorted(glob.glob(os.path.join(BASE, pat)))
             if os.path.isdir(d)]
# `murgia` and `more_power` are containers, not runs: listing them too would
# report every child's files a second time under the parent's name. Drop any
# directory that is an ancestor of another one in the list.
dirs = sorted(set(dirs))
dirs = [d for d in dirs
        if not any(o != d and o.startswith(d + os.sep) for o in dirs)]

print(f"{'run':<34s} {'#LOS':>5s} {'#snap':>6s}  pairs within dz<={ZTOL}")
for d in dirs:
    try:
        # recursive, so snapshots kept in their own subdirectory are found
        files = runs.list_hdf5(d)
    except Exception as e:
        print(f"{os.path.relpath(d, BASE):<34s}  ERROR {e}")
        continue
    zl, zs = [], []
    for p in files:
        try:
            if runs.is_los_file(p):
                z = runs.file_redshift(p)
                if z is not None:
                    zl.append((z, p))
            else:
                z = runs.snapshot_redshift(p)
                if z is not None:
                    zs.append((z, p))
        except Exception:
            pass
    # Match on the redshift READ FROM EACH FILE, never on the index in the
    # filename: los_0003 and snap_0003 are not the same epoch, and in these
    # runs they routinely are not.
    pairs, misses = [], []
    for zlos, plos in sorted(zl):
        if not zs:
            continue
        zsnap, psnap = min(zs, key=lambda t: abs(t[0] - zlos))
        (pairs if abs(zsnap - zlos) <= ZTOL else misses).append(
            (zlos, plos, zsnap, psnap))
    print(f"{os.path.relpath(d, BASE):<34s} {len(zl):>5d} {len(zs):>6d}  "
          f"{len(pairs)}")
    for zlos, plos, zsnap, psnap in pairs:
        print(f"{'':<34s}   MATCH  LOS {os.path.basename(plos):<22s} "
              f"z={zlos:.4f}   SNAP {os.path.relpath(psnap, d):<44s} "
              f"z={zsnap:.4f}   dz={abs(zsnap - zlos):.4f}")
        # the run directory and the redshift are all relosz.py needs: it
        # redoes this resolution itself, so nothing here has to guess a glob
        row = f"{zlos}\t{plos}\t{d}"
        if row not in rows:
            rows.append(row)
    for zlos, plos, zsnap, psnap in misses:
        print(f"{'':<34s}   no snap for LOS {os.path.basename(plos):<18s} "
              f"z={zlos:.4f}, nearest snapshot z={zsnap:.4f} "
              f"(dz={abs(zsnap - zlos):.4f})")

os.makedirs('logs', exist_ok=True)
with open('logs/v1_candidates.txt', 'w') as fh:
    fh.write("\n".join(rows) + ("\n" if rows else ""))
print(f"\n{len(rows)} candidate(s) written to logs/v1_candidates.txt")
PY

b "4b. ray range of every V1 candidate: a truncated file cannot be the reference"
echo "Appearing here does NOT mean a run is usable. Read the ray range: if"
echo "the rays stop short of the box the file is truncated, and relos.py"
echo "would be asked to reproduce something already wrong. Pick a clean one."
while IFS=$'\t' read -r Z LOS RUNDIR; do
    [ -z "$LOS" ] && continue
    echo
    echo "-------- candidate  z=$Z  $(basename "$LOS")  in $(basename "$RUNDIR")"
    python tests/t12_relos_roundtrip.py --a "$LOS" --max-los 200
    echo "   IF THIS ONE IS CLEAN, the V1 command is exactly:"
    echo "   sbatch scripts/sbatch_roundtrip.sh $RUNDIR $Z 100"
done < logs/v1_candidates.txt

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
echo "This file is truncated at 40 of 58.7372 Mpc. That is fine for V2 and"
echo "has to be said out loud when handing it over: a code-vs-code check"
echo "feeds both codes the same input, and a truncated input is truncated"
echo "identically for both. It is the tau ALGORITHM being compared, not the"
echo "physics of this particular box."
python scripts/make_los_subset.py --in "$CDM_LOS" --out data/los_cdm_z3_sample100.hdf5 --n 100

b "12. V2 - our own tau on those same 100 sightlines"
python stages/01_extract_los.py --los-file data/los_cdm_z3_sample100.hdf5 --npix 2048 \
    ${TC:+--treecool "$TC"} --out cache/cache_cdm_sample100.npz

b "13. V2 - the same cut from a clean (untruncated) run, if there is one"
echo "The production file above is the one whose tau goes in the paper, so it"
echo "is the one worth comparing. This second cut is the same test on a run"
echo "that is not truncated, so a disagreement cannot be blamed on the input."
CLEAN=$(head -1 logs/v1_candidates.txt 2>/dev/null | cut -f2)
if [ -n "$CLEAN" ] && [ -f "$CLEAN" ]; then
    echo "clean candidate: $CLEAN"
    python scripts/make_los_subset.py --in "$CLEAN" --out data/los_clean_sample100.hdf5 --n 100
    python stages/01_extract_los.py --los-file data/los_clean_sample100.hdf5 --npix 2048 \
        ${TC:+--treecool "$TC"} --out cache/cache_clean_sample100.npz
else
    echo "no candidate in logs/v1_candidates.txt - skipping"
fi

b "END. figures and caches produced"
ls -la figures/ cache/ data/*.hdf5 2>/dev/null
