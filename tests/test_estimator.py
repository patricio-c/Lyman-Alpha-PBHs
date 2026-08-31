#!/usr/bin/env python3
"""
Unit tests for common/p1d.py.

These are checks on ARITHMETIC, not on science. The inputs are analytic
(a cosine, a constant) and the expected outputs are known in closed form.
No simulation output is ever synthesised here; every scientific number in
this repository comes from a real snapshot.

    python tests/test_estimator.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.p1d import (flux, p1d_from_flux, solve_A, tau_eff, xi_from_flux,
                        logbin, ratio)

fails = []


def check(name, got, want, rtol=1e-9):
    ok = np.allclose(got, want, rtol=rtol)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {got!r} vs {want!r}")
    if not ok:
        fails.append(name)


n, dv = 2048, 2.19629
L = n * dv
v = np.arange(n) * dv

# 1. a single cosine mode must land in one bin with amplitude L a^2 / 4
for nmode, amp in ((7, 0.3), (31, 0.05)):
    d = amp * np.cos(2 * np.pi * nmode * v / L)
    F = (1.0 + d)[None, :]
    k, P, _ = p1d_from_flux(F, dv, global_mean=True)
    check(f"cosine n={nmode} peak bin", int(np.argmax(P)), nmode)
    check(f"cosine n={nmode} amplitude", P[nmode], L * amp ** 2 / 4, rtol=1e-8)
    check(f"cosine n={nmode} k", k[nmode], 2 * np.pi * nmode / L, rtol=1e-12)

# 2. Parseval: xi(0) must equal the variance of delta_F
d = 0.3 * np.cos(2 * np.pi * 7 * v / L) + 0.05 * np.cos(2 * np.pi * 31 * v / L)
F = (1.0 + d)[None, :]
r, xi = xi_from_flux(F, dv)
check("xi(0) == var(delta_F)", xi[0], np.var(d), rtol=1e-10)

# 3. solve_A must invert tau_eff exactly
tau = np.abs(1.0 + 0.5 * np.cos(2 * np.pi * 3 * v / L))[None, :] * 0.4
for target in (0.2, 0.3719, 0.8):
    A = solve_A(tau, target)
    check(f"solve_A round trip target={target}", tau_eff(tau, A), target,
          rtol=1e-9)

# 4. a constant flux field has zero power everywhere
k, P, _ = p1d_from_flux(np.full((3, n), 0.7), dv)
check("constant field has no power", float(np.abs(P[1:]).max()), 0.0,
      rtol=1e-12)

# 5. the ratio of a field to itself is 1 on every bin
kb, pb = logbin(*p1d_from_flux(F, dv)[:2], nbins=20)
_, rr = ratio(kb, pb, kb, pb, kgrid=kb)
check("ratio of a curve to itself", float(np.abs(rr - 1).max()), 0.0,
      rtol=1e-12)

# 6. rescaling by A=1 must be a no-op
check("flux A=1", flux(tau, 1.0), np.exp(-tau))

print()
if fails:
    print(f"{len(fails)} FAILED: {fails}")
    sys.exit(1)
print("all estimator unit tests pass")
