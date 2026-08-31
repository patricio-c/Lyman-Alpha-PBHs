"""
The UV background as an input, not a constant.

A TREECOOL file is seven columns of ASCII:

    log10(1+z)   gH0   gHe0   gHep   hH0   hHe0   hHep

with the photoionisation rates in s^-1 and the photoheating rates in
erg s^-1.  Interpolating gH0 at the redshift of the snapshot gives the
Gamma_HI that sets the neutral fraction, and therefore sets tau up to a
multiplicative constant.

That last point is worth stating plainly: within the optically thin
approximation tau scales as 1/Gamma_HI, so changing the UV background is
*exactly* the same operation as changing the rescaling factor A.  It cannot
change the shape of P1D and it cannot change a ratio between two runs that
share a background.  We keep it switchable anyway, because the first thing
anyone asks is whether the choice of background matters, and the honest
answer is easier to give with a flag than with an argument.

Anchor: HM12 gives Gamma_HI(z=3) = 8.2758e-13 s^-1.
"""

from __future__ import annotations

import os

import numpy as np

HM12_GAMMA_HI_Z3 = 8.2758e-13


def read_treecool(path):
    """Return (z, gH0, gHe0, gHep) arrays from a TREECOOL file."""
    d = np.loadtxt(path)
    if d.ndim != 2 or d.shape[1] < 4:
        raise SystemExit(f"{path}: expected at least 4 columns, got "
                         f"{d.shape}")
    z = 10.0 ** d[:, 0] - 1.0
    return z, d[:, 1], d[:, 2], d[:, 3]


def gamma_hi(z, table=None, value=None):
    """
    Photoionisation rate of HI at redshift z.

    Priority: an explicit --gamma-hi value, then a TREECOOL table, then the
    HM12 value at z=3 with a loud complaint if z is not 3.
    """
    if value is not None:
        return float(value), f"explicit ({value:.4e} s^-1)"
    if table is not None:
        zt, g, _, _ = read_treecool(table)
        order = np.argsort(zt)
        val = float(np.interp(z, zt[order], g[order]))
        return val, f"{os.path.basename(table)} interpolated at z={z:.3f}"
    if abs(z - 3.0) > 0.05:
        raise SystemExit(
            f"no UV background given and z = {z:.3f} is not 3. Pass either "
            f"--treecool <file> or --gamma-hi <value>; the built-in default "
            f"is only valid at z=3.")
    return HM12_GAMMA_HI_Z3, "HM12 built-in default at z=3"
