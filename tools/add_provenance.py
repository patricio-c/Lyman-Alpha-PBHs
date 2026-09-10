#!/usr/bin/env python3
"""
One-shot: make every stage record what produced its log.

Inserts `from common import prov` and replaces `log = []` with
`log = list(prov.header(...))`, so the first lines of every stage .txt are
the command, the git hash, the host and time, and the hashed inputs.

Idempotent: running it twice changes nothing.  Stages 00, 01 and 02 are
skipped because they write no log; stage 01's outputs are caches, which
already carry the same provenance through cache.save().

Run from the repository root, then commit the result:

    python tools/add_provenance.py .
"""
import sys
from pathlib import Path

# stages whose inputs are worth hashing, and the expression that names them
INPUTS = {
    "04_p1d.py": "args.caches",
    "10_dp1d.py": "[args.cache_ref, args.cache_test]",
    "11_drain_shape.py": "[args.boot_npz, args.pk_ref, args.pk_test]",
}

ANCHOR = ("sys.path.insert(0, os.path.dirname("
          "os.path.dirname(os.path.abspath(__file__))))")

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
touched, skipped = [], []

for f in sorted((root / "stages").glob("*.py")):
    src = f.read_text()
    if "    log = []\n" not in src:
        skipped.append(f"{f.name} (no log to stamp)")
        continue
    if "prov.header" in src:
        skipped.append(f"{f.name} (already done)")
        continue
    if ANCHOR not in src:
        skipped.append(f"{f.name} (no sys.path anchor - patch by hand)")
        continue

    src = src.replace(ANCHOR,
                      ANCHOR + "\nfrom common import prov  # noqa: E402", 1)
    arg = INPUTS.get(f.name)
    call = f"prov.header(inputs={arg})" if arg else "prov.header()"
    src = src.replace("    log = []\n", f"    log = list({call})\n", 1)
    f.write_text(src)
    touched.append(f.name)

print("patched:")
for t in touched:
    print("  ", t)
print("skipped:")
for s in skipped:
    print("  ", s)
