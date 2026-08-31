"""
Sightline resampling for P1D error bars.

Why this module exists
----------------------
Every error bar in this repo is a resampling over sightlines, and the naive
implementation puts the FFT inside the loop: 2000 bootstrap draws over a
6144 x 2048 cache is hours of wall time for a number that should take
seconds.  Two facts make it both cheap and exact.

1.  P1D is a *mean over sightlines* of per-sightline periodograms.  Compute
    the periodogram matrix once and every draw becomes a mean over rows.

2.  The global mean flux enters only as an overall factor.  With
    delta_i = F_i / Fbar - 1,

        rfft(delta_i)[k>0] = rfft(F_i)[k>0] / Fbar

    because the constant offset lives entirely in the k=0 mode, which the
    estimator drops anyway.  So define the Fbar-free periodogram

        Q_i(k) = |rfft(F_i)(k) * dv|^2 / L

    and for any subsample S,

        P1D_S(k) = mean_{i in S} Q_i(k) / Fbar_S^2

    exactly, not approximately.  Fbar_S is a mean over the per-sightline
    mean fluxes, so it is a mean over rows as well and costs nothing.  This
    matters: freezing Fbar at its full-sample value inside a bootstrap is a
    real (if small) bias, and here there is no reason to accept it.

Binning in k is a linear map, so it commutes with the mean over sightlines.
Binning the periodogram matrix *before* resampling turns each draw into a
mean over a [n_los, n_bins] array, and the whole bootstrap into a single
matrix product.  A draw is a multinomial count vector, which is what
sampling n items from n with replacement is; `counts @ Qb / n_los` does
every draw at once.

Consistency with common/p1d.py is not assumed.  `check_against_p1d` asserts
that this machinery reproduces `p1d.p1d_from_flux` on the full sample to
machine precision, and both scripts call it before doing anything else.  If
someone changes the estimator, this fails loudly instead of quietly
producing error bars for a different quantity.
"""

from __future__ import annotations

import numpy as np


# --- the Fbar-free periodogram --------------------------------------------

def periodogram(F, dv, chunk=512):
    """
    Per-sightline periodogram of the flux, with the mean flux left out.

    Returns (k [s/km], Q [n_los, n_k], mF [n_los]).  See the module
    docstring for why Q is the periodogram of F and not of the contrast.

    The FFT runs in chunks of `chunk` sightlines: a 6144 x 2048 float64
    cache would otherwise allocate a 100 MB complex intermediate on top of
    everything else, for no reason.
    """
    F = np.atleast_2d(np.asarray(F, dtype=np.float64))
    n_los, n_pix = F.shape
    L = n_pix * dv
    Q = np.empty((n_los, n_pix // 2 + 1), dtype=np.float64)
    mF = np.empty(n_los, dtype=np.float64)
    for s in range(0, n_los, chunk):
        e = min(s + chunk, n_los)
        blk = F[s:e]
        fk = np.fft.rfft(blk, axis=-1) * dv
        Q[s:e] = (fk.real ** 2 + fk.imag ** 2) / L
        mF[s:e] = blk.mean(axis=1)
    k = 2.0 * np.pi * np.fft.rfftfreq(n_pix, d=dv)
    return k, Q, mF


def check_against_p1d(k, Q, mF, dv, F, rtol=1e-9):
    """
    Assert that this module reproduces common.p1d on the full sample.

    Returns the worst fractional deviation over k > 0.  Raises if it is
    above `rtol`: at that point the error bars would belong to a different
    estimator than the measurement, which is exactly the class of silent
    bug this repository exists to catch.
    """
    from common.p1d import p1d_from_flux
    k_ref, P_ref, _ = p1d_from_flux(F, dv)
    P_mine = Q.mean(axis=0) / mF.mean() ** 2
    m = k_ref > 0
    worst = float(np.max(np.abs(P_mine[m] / np.maximum(P_ref[m], 1e-300) - 1.0)))
    if worst > rtol:
        raise SystemExit(
            f"common/boot.py disagrees with common/p1d.py by {worst:.3e} "
            f"(tolerance {rtol:.1e}). The resampling machinery and the "
            f"estimator have drifted apart. Fix that before trusting any "
            f"error bar from this script.")
    return worst


# --- binning as a linear operator -----------------------------------------

def bin_operator(k, nbins=24, kmin=None, kmax=None):
    """
    Log-spaced k binning expressed as a matrix.

    Returns (kb [n_bins], sel [bool over k], W [n_sel, n_bins]) so that

        Qb = Q[:, sel] @ W

    is the per-sightline binned power.  Empty bins are dropped, so kb and
    the columns of W are as long as the bins that actually got modes.

    The edges follow common.p1d.logbin exactly (geomspace over the selected
    range, digitize, clip), so a curve binned here overlays a curve binned
    by stage 04 without a half-bin offset.
    """
    m = k > 0
    if kmin is not None:
        m &= k >= kmin
    if kmax is not None:
        m &= k <= kmax
    kk = k[m]
    if kk.size == 0:
        raise SystemExit(f"no k modes in [{kmin}, {kmax}] s/km")
    edges = np.geomspace(kk.min() * 0.999, kk.max() * 1.001, nbins + 1)
    idx = np.clip(np.digitize(kk, edges) - 1, 0, nbins - 1)
    cnt = np.bincount(idx, minlength=nbins)
    ok = cnt > 0
    kb = np.bincount(idx, weights=kk, minlength=nbins)[ok] / cnt[ok]

    keep = np.flatnonzero(ok)
    W = np.zeros((kk.size, keep.size), dtype=np.float64)
    col = {b: j for j, b in enumerate(keep)}
    W[np.arange(kk.size), [col[b] for b in idx]] = 1.0
    W /= cnt[ok][None, :]
    return kb, m, W


def band_operator(k, kmin, kmax):
    """
    A single wide bin, for a band-averaged statistic.

    Same object as `bin_operator` with nbins=1, kept separate because the
    band is a different thing conceptually: the curve is for looking at,
    the band is the number that gets quoted.
    """
    kb, sel, W = bin_operator(k, nbins=1, kmin=kmin, kmax=kmax)
    return float(kb[0]), sel, W, int(sel.sum())


# --- resampling ------------------------------------------------------------

def draw_counts(n_los, n_boot, rng, dtype=np.float32):
    """
    Bootstrap draws as multinomial count vectors, [n_boot, n_los].

    Sampling n items out of n with replacement *is* a multinomial draw over
    the n indices, and the count form is what lets the whole bootstrap be
    one matrix product.  float32 keeps a 2000 x 6144 draw at 49 MB; counts
    are small integers so nothing is lost.
    """
    p = np.full(n_los, 1.0 / n_los)
    return rng.multinomial(n_los, p, size=n_boot).astype(dtype)


def apply_counts(counts, Qb, mF, chunk=256):
    """
    P1D of every weighted resample.  Returns [n_boot, n_bins].

    `counts` may be a bootstrap count matrix or any weight matrix whose rows
    sum to n_los (jackknife masks are handled by `subsample` instead).
    """
    n_los = Qb.shape[0]
    out = np.empty((counts.shape[0], Qb.shape[1]), dtype=np.float64)
    for s in range(0, counts.shape[0], chunk):
        e = min(s + chunk, counts.shape[0])
        c = counts[s:e].astype(np.float64)
        P = (c @ Qb) / n_los
        Fb = (c @ mF) / n_los
        out[s:e] = P / Fb[:, None] ** 2
    return out


def subsample(Qb, mF, mask):
    """P1D of the sightlines selected by a boolean mask.  [n_bins]."""
    if not mask.any():
        raise SystemExit("empty subsample")
    return Qb[mask].mean(axis=0) / mF[mask].mean() ** 2


def full_sample(Qb, mF):
    return Qb.mean(axis=0) / mF.mean() ** 2


# --- spatial blocks --------------------------------------------------------

def tile_labels(pos, box, nside):
    """
    Assign each sightline to a tile of an nside x nside grid on the
    transverse plane.  Returns integer labels in [0, nside^2).

    Both runs are tiled with the same grid over the same box, so tile b is
    the same physical region in CDM and in FCT.  That is the part of the
    pairing that survives stage 02's verdict: individual lines are not
    matched, but regions are, and the large-scale cosmic variance that
    pairing was supposed to cancel is a property of the region.
    """
    u = np.mod(np.asarray(pos, dtype=np.float64), box) / box
    ij = np.clip((u * nside).astype(np.int64), 0, nside - 1)
    return ij[:, 0] * nside + ij[:, 1]


def jackknife_var(theta_b, n_blocks=None):
    """
    Delete-one variance for a block jackknife.

    var = (B-1)/B * sum_b (theta_b - mean)^2.  Exact for equal block sizes
    and standard practice for roughly equal ones; the caller is expected to
    report the occupancy spread so the reader can judge.
    """
    t = np.asarray(theta_b, dtype=np.float64)
    B = int(n_blocks or t.size)
    return float((B - 1) / B * np.sum((t - t.mean()) ** 2))


# --- covariance and constant fits -----------------------------------------

def covariance(samples, hartlap=True):
    """
    Covariance of a [n_draw, n_bin] set of resampled statistics, and its
    inverse.

    Returns (C, Cinv, factor).  With `hartlap`, the inverse is multiplied by
    (N - p - 2)/(N - 1), the correction for the bias of an inverted noisy
    covariance.  It is derived for a covariance estimated from Gaussian
    samples; applied to a bootstrap covariance it is approximate, and the
    honest use is as a sanity margin rather than a precision correction.
    Keep N/p above ~50 and it barely matters.
    """
    S = np.asarray(samples, dtype=np.float64)
    N, p = S.shape
    C = np.cov(S, rowvar=False)
    C = np.atleast_2d(C)
    if N - p - 2 <= 0:
        raise SystemExit(
            f"{N} draws for {p} bins is not enough to invert a covariance. "
            f"Raise --n-boot or lower --nbins.")
    f = (N - p - 2) / (N - 1) if hartlap else 1.0
    return C, f * np.linalg.inv(C), f


def fit_linear(d, M, Cinv=None, err=None):
    """
    Generalised least squares for a linear model  d = M theta.

    `M` is [n_bin, n_par].  Returns (theta, cov_theta, chi2, dof) with
    dof = n_bin - n_par.  `Cinv` is the *inverse* covariance, as returned by
    `covariance` and already Hartlap-corrected if that was asked for; give
    `err` instead for the diagonal case.

    Everything fitted in t10 and t11 is a one- or two-parameter linear model
    (a constant, an amplitude times a fixed template, a power law in log
    space), so they all go through here and the covariance is handled the
    same way once.
    """
    d = np.asarray(d, dtype=np.float64)
    M = np.atleast_2d(np.asarray(M, dtype=np.float64))
    if M.shape[0] != d.size:
        M = M.T
    n, p = M.shape
    if n <= p:
        raise SystemExit(f"{n} bins for {p} parameters is not a fit")
    if Cinv is not None:
        Ci = np.atleast_2d(Cinv)
        if Ci.shape != (n, n):
            raise SystemExit(f"Cinv is {Ci.shape}, expected ({n}, {n})")
    else:
        Ci = np.diag(1.0 / np.asarray(err, dtype=np.float64) ** 2)
    A = M.T @ Ci @ M
    cov = np.linalg.inv(A)
    theta = cov @ (M.T @ Ci @ d)
    r = d - M @ theta
    return theta, cov, float(r @ Ci @ r), n - p


def fit_constant(d, Cinv=None, err=None):
    """
    Best-fit constant to d.  Returns (c, sigma_c, chi2, dof), dof = n - 1.
    Thin wrapper over `fit_linear` with a column of ones.
    """
    d = np.asarray(d, dtype=np.float64)
    th, cov, chi2, dof = fit_linear(d, np.ones((d.size, 1)), Cinv=Cinv,
                                    err=err)
    return float(th[0]), float(np.sqrt(cov[0, 0])), chi2, dof


def fit_template(d, template, Cinv=None, err=None):
    """
    Best-fit amplitude of a fixed shape:  d = amp * template.

    Returns (amp, sigma_amp, chi2, dof), dof = n - 1, the same dof as
    `fit_constant`, so the two chi2 values compare directly.  That
    comparison is the point: a constant and (k_max^2 - k^2) are the same
    model to within (k/k_max)^2, and if the data cannot tell them apart
    that IS the measurement of k << k_max.
    """
    d = np.asarray(d, dtype=np.float64)
    t = np.asarray(template, dtype=np.float64).reshape(-1, 1)
    th, cov, chi2, dof = fit_linear(d, t, Cinv=Cinv, err=err)
    return float(th[0]), float(np.sqrt(cov[0, 0])), chi2, dof
