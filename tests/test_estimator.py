#!/usr/bin/env python3
"""
Unit tests for common/p1d.py and common/boot.py.

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


def check(name, got, want, rtol=1e-9, atol=1e-11):
    ok = np.allclose(got, want, rtol=rtol, atol=atol)
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


# --- common/boot.py -------------------------------------------------------
# The resampling module claims three things that would be silent if wrong:
# that dividing the mean flux out of the periodogram is exact, that binning
# commutes with the average over sightlines, and that a multinomial count
# vector is the same thing as drawing indices with replacement. t9 and t10
# are built entirely on those three, so they get checked here.

from common import boot  # noqa: E402

rng = np.random.default_rng(7)
nl, npx = 240, 256
Fm = (0.7
      + 0.05 * np.cos(2 * np.pi * 3 * np.arange(npx) / npx)[None, :]
      + 0.03 * rng.standard_normal((nl, 1))
      + 0.02 * rng.standard_normal((nl, npx)))

kk, Q, mF = boot.periodogram(Fm, dv, chunk=64)
k_ref, P_ref, _ = p1d_from_flux(Fm, dv)
msk = k_ref > 0
check("boot: Q/Fbar^2 reproduces p1d_from_flux",
      float(np.max(np.abs(Q.mean(axis=0)[msk] / mF.mean() ** 2
                          / P_ref[msk] - 1))), 0.0, rtol=1e-11)
check("boot: chunking changes nothing",
      float(np.abs(boot.periodogram(Fm, dv, chunk=1000)[1] - Q).max()), 0.0)

kb_b, sel, W = boot.bin_operator(kk, nbins=20)
Qb = Q[:, sel] @ W
kb_ref, P_binned_ref = logbin(k_ref, P_ref, nbins=20)
check("boot: bin_operator centres match logbin",
      float(np.abs(kb_b / kb_ref - 1).max()), 0.0)
check("boot: bin-then-average == average-then-bin",
      float(np.abs(boot.full_sample(Qb, mF) / P_binned_ref - 1).max()), 0.0)

cnt = boot.draw_counts(nl, 64, rng)
check("boot: every draw resamples n out of n",
      float(np.abs(cnt.sum(axis=1) - nl).max()), 0.0, rtol=1e-12)
check("boot: all-ones weights give the full sample",
      float(np.abs(boot.apply_counts(np.ones((1, nl), dtype=np.float32),
                                     Qb, mF)[0]
                   / boot.full_sample(Qb, mF) - 1).max()), 0.0)

idx = np.random.default_rng(11).integers(0, nl, size=(200, nl))
brute = np.array([Qb[i].mean(axis=0) / mF[i].mean() ** 2 for i in idx])
cmat = np.array([np.bincount(i, minlength=nl) for i in idx], dtype=np.float64)
check("boot: count matmul == explicit index gather",
      float(np.abs(boot.apply_counts(cmat, Qb, mF) / brute - 1).max()), 0.0)

sub = np.zeros(nl, bool)
sub[::3] = True
check("boot: subsample refits Fbar on the subsample",
      float(np.abs(boot.subsample(Qb, mF, sub)
                   / logbin(*p1d_from_flux(Fm[sub], dv)[:2],
                            nbins=20)[1] - 1).max()), 0.0)

pos = np.array([[0.1, 0.1], [0.9, 0.1], [0.1, 0.9], [0.6, 0.6], [4.1, 0.1]])
lab = boot.tile_labels(pos, 1.0, 4)
check("boot: tiling wraps periodically", int(lab[0]), int(lab[4]))
check("boot: tiling separates corners", len(set(lab[:4].tolist())), 4)

t = np.array([1.0, 2.0, 3.0, 4.0])
check("boot: delete-one jackknife variance",
      boot.jackknife_var(t), 0.75 * np.sum((t - t.mean()) ** 2), rtol=1e-13)

dv_vec = np.array([1.0, 1.2, 0.9, 1.1])
ev = np.full(4, 0.1)
c1, s1, x1, _ = boot.fit_constant(dv_vec, err=ev)
c2, s2, x2, _ = boot.fit_constant(dv_vec, Cinv=np.diag(1.0 / ev ** 2))
check("boot: fit_constant diagonal Cinv == err form", (c1, s1, x1),
      (c2, s2, x2), rtol=1e-12)
check("boot: fit_constant recovers the mean", c1, float(dv_vec.mean()),
      rtol=1e-13)

S = rng.standard_normal((2000, 8)) * 0.3
C, Ci, hf = boot.covariance(S)
check("boot: Hartlap factor", hf, (2000 - 8 - 2) / (2000 - 1), rtol=1e-12)
check("boot: Cinv is the corrected inverse",
      float(np.abs(Ci @ C - hf * np.eye(8)).max()), 0.0)


# --- fit_linear / fit_template --------------------------------------------
# The constant and the (k_max^2 - k^2) template are the same model to within
# (k/k_max)^2, and t11's first result is that the data cannot tell them
# apart. Here the input is exact, so the failure to distinguish must be
# exact too: that is the check.

kf = np.geomspace(1e-3, 2.5e-2, 14)
kmx = 5.0
tmpl = kmx ** 2 - kf ** 2
amp_true = 3.7e-5
d_exact = amp_true * tmpl
ef = np.full(14, 1e-9)
a_t, _, x2_t, dof_t = boot.fit_template(d_exact, tmpl, err=ef)
check("boot: fit_template recovers a known amplitude", a_t, amp_true,
      rtol=1e-12)
check("boot: exact template fit has zero chi2", x2_t, 0.0)
c_t, _, x2_ct, dof_ct = boot.fit_constant(d_exact, err=ef)
check("boot: template and constant fits share dof", dof_t, dof_ct)
check("boot: a k<<k_max template is a constant to (k/k_max)^2",
      float(np.abs(c_t / (amp_true * kmx ** 2) - 1)),
      float((kf ** 2).mean() / kmx ** 2), rtol=0.05)

Mlin = np.column_stack([np.ones(14), np.log(kf)])
tt = np.array([1.7, -1.25])
th_l, cov_l, x2_l, dof_l = boot.fit_linear(Mlin @ tt, Mlin, err=ef)
check("boot: fit_linear recovers a known power law", th_l, tt, rtol=1e-10)
check("boot: fit_linear dof = n - npar", dof_l, 12)


# --- growth factor and the fixed window -----------------------------------
# Einstein-de Sitter has D(a) = a exactly, which is the only closed-form
# anchor available for a numerically integrated growth factor.

from common.units import Cosmology, common_window, desi_window, growth_factor  # noqa: E402

eds = Cosmology(Om=1.0)
for zz in (1.0, 3.0, 9.0):
    check(f"growth: EdS D(z={zz:g}) == 1/(1+z)", growth_factor(zz, eds),
          1.0 / (1.0 + zz), rtol=1e-5)
check("growth: D(0) == 1", growth_factor(0.0), 1.0, rtol=1e-12)
check("growth: D decreases with z",
      bool(np.all(np.diff([growth_factor(z) for z in (0, 1, 2, 3, 4)]) < 0)),
      True)

zlist = [2.2, 2.6, 3.0, 3.4, 4.0]
w1, w2 = common_window(zlist)
check("window: k1 is the largest of the individual k1",
      w1, max(desi_window(z)[0] for z in zlist), rtol=1e-12)
check("window: k2 is the smallest of the individual k2",
      w2, min(desi_window(z)[1] for z in zlist), rtol=1e-12)
check("window: the fixed k2 sits at the lowest-z bin",
      w2, desi_window(min(zlist))[1], rtol=1e-12)

print()
if fails:
    print(f"{len(fails)} FAILED: {fails}")
    sys.exit(1)
print("all estimator unit tests pass")
