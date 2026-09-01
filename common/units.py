"""
Cosmology and unit conversions.

Everything here is derived from a single (h, Omega_m) pair plus the comoving
box size.  Nothing is hard-coded to z=3: the z=3 numbers quoted in the
docstrings are what this module reproduces for the 40 Mpc/h runs, and are
kept only as a regression anchor.

Regression anchors for box = 40 Mpc/h, h = 0.681, Omega_m = 0.3053, z = 3:

    L_com   = 58.7372 Mpc
    H(z)    = 306.31  km/s/Mpc
    v_box   = 4498.0  km/s          (= H(z) * L_com / (1+z))
    dv      = 2.19629 km/s          (2048 pixels)
    k[Mpc^-1] = k[s/km] * 76.578
    k[h/Mpc]  = k[s/km] * 112.45
"""

from __future__ import annotations

import numpy as np

C_KMS = 299792.458            # speed of light, km/s
LAMBDA_LYA = 1215.67          # Lyman-alpha rest wavelength, Angstrom


class Cosmology:
    """Flat LCDM.  Only what the forest needs."""

    def __init__(self, h=0.681, Om=0.3053, Ob=0.0486):
        self.h = float(h)
        self.Om = float(Om)
        self.Ob = float(Ob)
        self.OL = 1.0 - self.Om
        self.H0 = 100.0 * self.h      # km/s/Mpc

    def H(self, z):
        """Hubble parameter in km/s/Mpc.  Matter + Lambda only (fine at z<10)."""
        return self.H0 * np.sqrt(self.Om * (1.0 + z) ** 3 + self.OL)

    def hubble_velocity_scale(self, z):
        """
        Factor that turns a comoving distance [Mpc] into a velocity [km/s]
        along the line of sight:  v = x_com * H(z) / (1+z).
        """
        return self.H(z) / (1.0 + z)

    def box_kms(self, box_mpc, z):
        """Velocity extent of the box at redshift z, in km/s."""
        return box_mpc * self.hubble_velocity_scale(z)

    def dv(self, box_mpc, z, npix):
        return self.box_kms(box_mpc, z) / float(npix)

    def k_kms_to_mpc(self, z):
        """Multiply k[s/km] by this to get k[Mpc^-1]."""
        return self.hubble_velocity_scale(z)

    def k_kms_to_hmpc(self, z):
        """Multiply k[s/km] by this to get k[h/Mpc]."""
        return self.hubble_velocity_scale(z) / self.h

    def summary(self, z, box_mpc=None, npix=None):
        lines = [f"h = {self.h:.4f}   Omega_m = {self.Om:.4f}   z = {z:.3f}",
                 f"H(z)          = {self.H(z):.3f} km/s/Mpc",
                 f"k[s/km] -> k[Mpc^-1]  x {self.k_kms_to_mpc(z):.4f}",
                 f"k[s/km] -> k[h/Mpc]   x {self.k_kms_to_hmpc(z):.4f}"]
        if box_mpc is not None:
            lines.append(f"box           = {box_mpc:.4f} Mpc "
                         f"= {box_mpc * self.h:.4f} Mpc/h")
            lines.append(f"v_box         = {self.box_kms(box_mpc, z):.3f} km/s")
            if npix:
                lines.append(f"dv ({npix} px) = "
                             f"{self.dv(box_mpc, z, npix):.5f} km/s")
                lines.append(f"k_fundamental = "
                             f"{2*np.pi/self.box_kms(box_mpc, z):.6f} s/km")
                lines.append(f"k_nyquist     = "
                             f"{np.pi/self.dv(box_mpc, z, npix):.6f} s/km")
        return "\n".join(lines)


# Default cosmology of the runs.  Both CDM and FCT share it: the FCT model
# modifies the primordial spectrum, not the background expansion.
DEFAULT = Cosmology()


# --- observational anchors -------------------------------------------------

def tau_eff_turner24(z):
    """Turner et al. 2024 effective optical depth.  tau_eff(3.0) = 0.3719."""
    return 2.46e-3 * (1.0 + z) ** 3.62


TAU_EFF_SHERWOOD_Z3 = 0.40839


def desi_window(z, resolution_ang=0.8):
    """
    (k_min, k_max) of the DESI DR1 P1D window in s/km.
    Upper end is the half-Nyquist of the spectrograph resolution.
    """
    R_z = C_KMS * resolution_ang / ((1.0 + z) * LAMBDA_LYA)
    return 1.0e-3, 0.5 * np.pi / R_z


# --- linear growth ---------------------------------------------------------

def _trapz(y, x):
    return (np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(y, x)


def growth_factor(z, cosmo=None, n=4096, a_min=1e-8):
    """
    Linear growth factor D(z), normalised to D(0) = 1.

    Flat LCDM, matter + Lambda, the standard quadrature

        D(a)  =  (5 Om / 2) E(a) int_0^a da' / (a' E(a'))^3

    integrated in ln a, which is where the integrand is smooth.  No scipy:
    this gets called once per redshift, not in a loop.

    Why it is here: the Poisson isocurvature term of the FCT model is fixed
    in comoving scale and grows as D^2(z), which is the one piece of the
    redshift evolution that is calculable rather than calibrated.  Anything
    that separates "what the growth factor already explains" from "what is
    left" needs this.

    Regression anchor: for Omega_m = 1 (Einstein-de Sitter) D(a) = a
    exactly, so D(z=1) = 0.5.  tests/test_estimator.py checks that.
    """
    cosmo = cosmo or DEFAULT

    def E(a):
        return np.sqrt(cosmo.Om * a ** -3 + cosmo.OL)

    def unnorm(a):
        aa = np.geomspace(a_min, a, n)
        # integrand of int da'/(a' E)^3 rewritten in ln a': da' = a' dln a',
        # so the integrand is 1/(a'^2 E^3) and the abscissa is ln a'.
        return (2.5 * cosmo.Om * E(a)
                * _trapz(1.0 / (aa ** 2 * E(aa) ** 3), np.log(aa)))

    a = 1.0 / (1.0 + np.asarray(z, dtype=np.float64))
    d0 = unnorm(1.0)
    if a.ndim == 0:
        return float(unnorm(float(a)) / d0)
    return np.array([unnorm(float(x)) for x in a]) / d0


# --- a window that does not move with redshift -----------------------------

def common_window(zs, resolution_ang=0.8):
    """
    The (k1, k2) intersection of the DESI windows of several redshifts.

    `desi_window` returns a k_max that grows with z, because the resolution
    in s/km improves as (1+z) rises.  Integrating P1D over that window bin
    by bin means integrating over a DIFFERENT window in each bin, and the
    "redshift evolution" that comes out then contains a purely instrumental
    component.  Any integrated observable A(z) must use one fixed window for
    every bin, and this is it.  Fixing it is five lines; not fixing it is a
    referee report.
    """
    zs = np.atleast_1d(np.asarray(zs, dtype=np.float64))
    lo = max(desi_window(float(z), resolution_ang)[0] for z in zs)
    hi = min(desi_window(float(z), resolution_ang)[1] for z in zs)
    if not hi > lo:
        raise ValueError(f"no common window across z = {zs}: "
                         f"k1 = {lo:.5g} >= k2 = {hi:.5g} s/km")
    return lo, hi
