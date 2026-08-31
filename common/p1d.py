"""
Flux statistics from an array of optical depths.

This module is the whole downstream half of the pipeline.  It takes
tau[n_los, n_pix] and dv [km/s] and produces everything else: the flux, the
effective optical depth, the rescaling factor A, the 1D power spectrum, and
the correlation function.

Design boundary: nothing in here knows where tau came from.  Swap our
extractor for spectWizard and every number below is computed the same way.

Conventions
-----------
    F      = exp(-A * tau)
    tau_eff = -ln(<F>)                 mean over all pixels of all LOS
    delta_F = F / <F> - 1              <F> global, not per-LOS
    P1D(k) = <|FFT(delta_F) * dv|^2> / L      L = n_pix * dv
    k      = 2 pi * rfftfreq(n_pix, dv)

P1D comes out in km/s, k in s/km.  The quantity usually plotted is
k*P1D/pi, which is dimensionless.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


# --- effective optical depth and the rescaling factor A --------------------

def tau_eff(tau, A=1.0):
    """Effective optical depth of the whole sample for a rescaling A."""
    mF = np.mean(np.exp(-A * np.asarray(tau, dtype=np.float64)))
    return -np.log(max(float(mF), 1e-300))


def solve_A(tau, target, lo=1e-6, hi=1e3):
    """
    Find A such that tau_eff(A*tau) == target.

    tau_eff is monotonically increasing in A, so brentq is safe.  Raises if
    the target is outside what this tau field can reach.
    """
    tau = np.asarray(tau, dtype=np.float64)
    f = lambda a: tau_eff(tau, a) - target
    flo, fhi = f(lo), f(hi)
    if flo > 0:
        raise ValueError(f"target tau_eff={target:.5f} below the floor "
                         f"reachable by this field ({tau_eff(tau, lo):.5f})")
    if fhi < 0:
        raise ValueError(f"target tau_eff={target:.5f} above the ceiling "
                         f"reachable by this field ({tau_eff(tau, hi):.5f})")
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-12)


# --- the estimator ---------------------------------------------------------

def flux(tau, A=1.0):
    return np.exp(-A * np.asarray(tau, dtype=np.float64))


def p1d_from_flux(F, dv, global_mean=True):
    """
    P1D from a flux array [n_los, n_pix].

    Returns (k, P, err) with k in s/km and P in km/s.  err is the error on
    the mean across LOS (std / sqrt(n_los)), which is the sampling error of
    the estimator, not a cosmic-variance error bar.
    """
    F = np.atleast_2d(np.asarray(F, dtype=np.float64))
    n_los, n_pix = F.shape
    Fbar = F.mean() if global_mean else F.mean(axis=1, keepdims=True)
    d = F / Fbar - 1.0

    L = n_pix * dv
    fk = np.fft.rfft(d, axis=-1) * dv
    Pk = (fk.real ** 2 + fk.imag ** 2) / L          # [n_los, n_pix//2+1]

    k = 2.0 * np.pi * np.fft.rfftfreq(n_pix, d=dv)
    P = Pk.mean(axis=0)
    err = Pk.std(axis=0) / np.sqrt(max(n_los, 1))
    return k, P, err


def p1d_from_tau(tau, dv, target=None, A=None, global_mean=True):
    """
    Convenience wrapper.  Exactly one of `target` (rescale to this tau_eff)
    or `A` (use this rescaling) may be given; neither means A=1, i.e. the
    raw field with no rescaling at all.

    Returns (k, P, err, A).
    """
    if target is not None and A is not None:
        raise ValueError("give either target or A, not both")
    if target is not None:
        A = solve_A(tau, target)
    elif A is None:
        A = 1.0
    k, P, err = p1d_from_flux(flux(tau, A), dv, global_mean=global_mean)
    return k, P, err, A


# --- binning and masking ---------------------------------------------------

def logbin(k, P, err=None, nbins=28, kmin=None, kmax=None):
    """Average P into log-spaced k bins.  Drops the k=0 mode."""
    m = k > 0
    if kmin is not None:
        m &= k >= kmin
    if kmax is not None:
        m &= k <= kmax
    kk, pp = k[m], P[m]
    edges = np.geomspace(kk.min() * 0.999, kk.max() * 1.001, nbins + 1)
    idx = np.digitize(kk, edges) - 1
    idx = np.clip(idx, 0, nbins - 1)
    cnt = np.bincount(idx, minlength=nbins)
    ksum = np.bincount(idx, weights=kk, minlength=nbins)
    psum = np.bincount(idx, weights=pp, minlength=nbins)
    ok = cnt > 0
    out_k = ksum[ok] / cnt[ok]
    out_p = psum[ok] / cnt[ok]
    if err is None:
        return out_k, out_p
    e2 = np.bincount(idx, weights=err[m] ** 2, minlength=nbins)
    out_e = np.sqrt(e2[ok]) / cnt[ok]
    return out_k, out_p, out_e


def nyquist_cut(k, dv, frac=0.5):
    """
    Boolean mask keeping k below `frac` times the Nyquist frequency.

    The default frac=0.5 exists because the P1D estimator aliases badly in
    the top octave: above half-Nyquist the ratio between two runs picks up
    structure that is pure grid artefact.
    """
    return (k > 0) & (k <= frac * np.pi / dv)


# --- ratio between two runs on a common grid -------------------------------

def ratio(k_a, P_a, k_b, P_b, kgrid=None):
    """
    P_a / P_b interpolated onto a common grid.

    If kgrid is None the grid is the overlap of the two, taken from A.  Both
    curves are interpolated in log-log, which matters at low k where P is a
    steep power law and linear interpolation biases the ratio low.
    """
    if kgrid is None:
        lo = max(k_a[k_a > 0].min(), k_b[k_b > 0].min())
        hi = min(k_a.max(), k_b.max())
        kgrid = k_a[(k_a >= lo) & (k_a <= hi)]
    la = np.interp(np.log(kgrid), np.log(k_a[k_a > 0]),
                   np.log(np.maximum(P_a[k_a > 0], 1e-300)))
    lb = np.interp(np.log(kgrid), np.log(k_b[k_b > 0]),
                   np.log(np.maximum(P_b[k_b > 0], 1e-300)))
    return kgrid, np.exp(la - lb)


# --- configuration space ---------------------------------------------------

def xi_from_flux(F, dv, rmax_kms=None, global_mean=True):
    """
    1D flux correlation function xi(dv) along the line of sight, by FFT.

    This is the Fourier transform of exactly the same estimator as above, so
    it is not an independent measurement of the *field*; it is an independent
    check of the *estimator and its binning*.  If a feature in the P1D ratio
    survives here without a matching feature in xi, the feature is in the
    binning, not the data.
    """
    F = np.atleast_2d(np.asarray(F, dtype=np.float64))
    n_los, n_pix = F.shape
    Fbar = F.mean() if global_mean else F.mean(axis=1, keepdims=True)
    d = F / Fbar - 1.0

    fk = np.fft.rfft(d, axis=-1)
    xi = np.fft.irfft((fk.real ** 2 + fk.imag ** 2), n=n_pix, axis=-1) / n_pix
    xi = xi.mean(axis=0)

    half = n_pix // 2
    r = np.arange(half) * dv
    xi = xi[:half]
    if rmax_kms is not None:
        m = r <= rmax_kms
        r, xi = r[m], xi[m]
    return r, xi
