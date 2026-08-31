#!/usr/bin/env bash
# Does the repo actually work on this machine, with these caches?
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
fail=0
echo "== python =="
python -c "import numpy, scipy, h5py, matplotlib; print('  deps ok')" || fail=1
echo "== modules compile =="
python -m compileall -q common stages tests >/dev/null && echo "  ok" || fail=1
echo "== legacy present =="
for f in swift_extract.py forest_tools.py ionization.py; do
  [ -f "legacy/$f" ] && echo "  $f" || { echo "  MISSING $f"; fail=1; }
done
echo "== caches readable =="
for f in cache/cache_cdm.npz cache/cache_fct.npz; do
  [ -f "$f" ] && python -c "
from common import cache
c = cache.load('$f'); print(c.describe())" || { echo "  MISSING $f"; fail=1; }
done
echo
[ $fail -eq 0 ] && echo "VERIFY OK" || echo "VERIFY FAILED"
exit $fail
