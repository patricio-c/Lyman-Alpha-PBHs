"""
Reading and writing the tau caches, with provenance.

A cache is a .npz holding at minimum:

    tau   float32/float64  [n_los, n_pix]   optical depth, NOT rescaled
    dv    float                             pixel width in km/s
    z     float                             redshift of the snapshot

and, when written by this repo, also:

    box_mpc, nlos, npix, run, engine, uvb, hubble_mode,
    prov   a JSON string with git hash, argv, host, time, input hashes

Caches written before this repo existed only have the first three keys.
`load` accepts them and fills the rest with NaN / "unknown", so nothing
breaks when you point a new script at an old file.
"""

from __future__ import annotations

import hashlib
import json
import os
import socket
import subprocess
import sys
import time

import numpy as np

REQUIRED = ("tau", "dv", "z")


# --- provenance ------------------------------------------------------------

def git_hash(default="not-a-git-checkout"):
    try:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        out = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            dirty = subprocess.run(["git", "-C", here, "status", "--porcelain"],
                                   capture_output=True, text=True, timeout=5)
            suffix = "-dirty" if dirty.stdout.strip() else ""
            return out.stdout.strip() + suffix
    except Exception:
        pass
    return default


def file_hash(path, nbytes=1 << 20):
    """
    MD5 of the first and last MB of a file plus its size.  Full hashes of
    100 GB snapshots are not worth the wall time; this is enough to catch
    "I regenerated the snapshot and forgot".
    """
    try:
        size = os.path.getsize(path)
        h = hashlib.md5(str(size).encode())
        with open(path, "rb") as f:
            h.update(f.read(nbytes))
            if size > 2 * nbytes:
                f.seek(-nbytes, os.SEEK_END)
                h.update(f.read(nbytes))
        return h.hexdigest()[:16]
    except Exception as exc:
        return f"unhashable:{exc.__class__.__name__}"


def stamp(inputs=None, extra=None):
    """Provenance dict for whatever is about to be written."""
    d = {
        "git": git_hash(),
        "argv": " ".join(sys.argv),
        "host": socket.gethostname(),
        "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }
    if inputs:
        d["inputs"] = {p: file_hash(p) for p in inputs}
    if extra:
        d.update(extra)
    return d


# --- read / write ----------------------------------------------------------

def save(path, tau, dv, z, inputs=None, **meta):
    payload = {
        "tau": np.asarray(tau, dtype=np.float32),
        "dv": float(dv),
        "z": float(z),
        "nlos": int(np.atleast_2d(tau).shape[0]),
        "npix": int(np.atleast_2d(tau).shape[1]),
        "prov": json.dumps(stamp(inputs=inputs), indent=1),
    }
    for k, v in meta.items():
        payload[k] = v if isinstance(v, np.ndarray) else np.array(v)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    np.savez_compressed(path, **payload)
    return path


class Cache(dict):
    """Attribute access over the npz contents, with sane defaults."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    def describe(self):
        lines = [f"  file      {self.get('_path', '?')}",
                 f"  tau       {self['tau'].shape}  "
                 f"{self['tau'].dtype}  "
                 f"[{np.nanmin(self['tau']):.3e}, "
                 f"{np.nanmax(self['tau']):.3e}]",
                 f"  z         {self['z']:.4f}",
                 f"  dv        {self['dv']:.5f} km/s",
                 f"  run       {self.get('run', 'unknown')}",
                 f"  engine    {self.get('engine', 'unknown')}",
                 f"  uvb       {self.get('uvb', 'unknown')}",
                 f"  hubble    {self.get('hubble_mode', 'unknown')}"]
        return "\n".join(lines)


def load(path, z_expect=None, z_tol=0.02):
    d = np.load(path, allow_pickle=False)
    missing = [k for k in REQUIRED if k not in d]
    if missing:
        raise SystemExit(f"{path}: missing required keys {missing}. "
                         f"Present: {sorted(d.files)}")

    c = Cache()
    c["_path"] = path
    c["tau"] = np.atleast_2d(d["tau"]).astype(np.float64)
    c["dv"] = float(d["dv"])
    c["z"] = float(d["z"])
    for k in d.files:
        if k in REQUIRED:
            continue
        v = d[k]
        c[k] = v.item() if v.ndim == 0 else v
    c.setdefault("run", "unknown")
    c.setdefault("engine", "unknown")
    c.setdefault("uvb", "unknown")
    c.setdefault("hubble_mode", "unknown")

    if not np.isfinite(c["tau"]).all():
        n = int((~np.isfinite(c["tau"])).sum())
        raise SystemExit(f"{path}: {n} non-finite values in tau")
    if (c["tau"] < 0).any():
        raise SystemExit(f"{path}: negative optical depths")
    if z_expect is not None and abs(c["z"] - z_expect) > z_tol:
        raise SystemExit(f"{path}: z = {c['z']:.4f}, expected {z_expect:.4f}")
    return c


def load_pair(path_cdm, path_fct, z_tol=0.02):
    """
    Load two caches and refuse to continue if they are not comparable.

    Every test in this repo compares two runs pixel-for-pixel or LOS-for-LOS,
    which is only meaningful when the two extractions used the same geometry.
    This is where that gets enforced once instead of in each script.
    """
    a = load(path_cdm)
    b = load(path_fct)
    if abs(a.z - b.z) > z_tol:
        raise SystemExit(f"redshift mismatch: {a.z:.4f} vs {b.z:.4f}")
    if abs(a.dv - b.dv) / a.dv > 1e-6:
        raise SystemExit(f"dv mismatch: {a.dv:.6f} vs {b.dv:.6f} km/s. "
                         f"If this is intentional (Hubble-rescaling test) "
                         f"use load() twice and say so explicitly.")
    if a.tau.shape != b.tau.shape:
        raise SystemExit(f"shape mismatch: {a.tau.shape} vs {b.tau.shape}")
    return a, b
