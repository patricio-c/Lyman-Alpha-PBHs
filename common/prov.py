"""
Provenance header for stage logs.

Every stage writes a .txt next to its figure.  Until this module existed
those logs recorded the numbers but not the invocation, so months later a
figure could be read but not reproduced: you could see that a run used
`cache_fct.npz` and `tau_eff = 0.42461`, but not whether it was called with
`--exact-boot 300` or against which commit of the code, and both change
the error bars.

`cache.stamp()` already collects exactly the right fields for the caches.
This is the same information formatted for a human reading a log file, so
that the first thing in every stage output answers "what produced this".

Usage in a stage, replacing `log = []`:

    log = list(prov.header())

and, when the stage knows its input files:

    log = list(prov.header(inputs=[args.cache_ref, args.cache_test]))

Input hashing is the first and last MB plus the size, not the whole file,
because full hashes of 100 GB snapshots are not worth the wall time.  It
is enough to catch "I regenerated the input and forgot".
"""

from __future__ import annotations

import os

from . import cache as _cache


def header(inputs=None, echo=True, width=74):
    """
    Lines identifying this run: command, code version, host, time, inputs.

    Returns the lines and, unless echo is False, prints them too, so that
    the terminal and the log file say the same thing.

    A missing or unreadable input is reported rather than raised: a stage
    that can run without one of its optional files should still record
    that the file was absent.
    """
    st = _cache.stamp()
    lines = ["-" * width,
             f"command : {st['argv']}",
             f"git     : {st['git']}",
             f"host    : {st['host']}    run at {st['time']}",
             f"python  : {st['python']}   numpy {st['numpy']}",
             f"cwd     : {os.getcwd()}"]
    for p in (inputs or []):
        if p is None:
            continue
        if os.path.exists(p):
            lines.append(f"input   : {p}  "
                         f"[{_cache.file_hash(p)}, "
                         f"{os.path.getsize(p)} bytes]")
        else:
            lines.append(f"input   : {p}  [MISSING at run time]")
    lines.append("-" * width)
    if echo:
        for ln in lines:
            print(ln)
    return lines
