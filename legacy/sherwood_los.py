"""
Lector del formato binario de espectros de la simulacion Sherwood
(Bolton et al. 2017, MNRAS 464, 897) + estimador del espectro de potencias
1D del flujo transmitido (P1D).

Traduccion directa de read_spectra.pro / get_uvb.pro / get_wavelength.pro
(James Bolton, Sep 2017), con el estimador de P1D agregado.

Layout del binario (stream plano, sin record markers de Fortran):

    npix     int32
    nlos     int32
    ztime    float32   redshift
    omegaM   float32
    omegaL   float32
    omegab   float32
    hubble   float32   h = H0/100
    boxsize  float32   [h^-1 ckpc]
    axis     int32   [nlos]        1=x, 2=y, 3=z
    coord1   float32 [nlos]        [h^-1 ckpc]
    coord2   float32 [nlos]        [h^-1 ckpc]
    pixels   float32 [npix]        posiciones de pixel [h^-1 ckpc]
    tau      float32 [npix*nlos]   profundidad optica HI Lya, C-order (los, pix)
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np

LAMBDA_LYA = 1215.6701  # [Angstrom], reposo
CKMS = 2.99792458e5     # [km/s]

HEADER_BYTES = 2 * 4 + 6 * 4  # npix, nlos + 6 floats


# ---------------------------------------------------------------------------
# Lectura
# ---------------------------------------------------------------------------

@dataclass
class SherwoodSpectra:
    """Contenedor de un archivo de espectros de Sherwood."""

    z: float
    omega_m: float
    omega_l: float
    omega_b: float
    hubble: float
    boxsize: float          # [h^-1 ckpc]
    axis: np.ndarray        # [nlos] int32
    coord1: np.ndarray      # [nlos] float32 [h^-1 ckpc]
    coord2: np.ndarray      # [nlos] float32 [h^-1 ckpc]
    pixels: np.ndarray      # [npix] float32 [h^-1 ckpc]
    tau: np.ndarray         # [nlos, npix] float32

    @property
    def nlos(self) -> int:
        return self.tau.shape[0]

    @property
    def npix(self) -> int:
        return self.tau.shape[1]

    # -- geometria en velocidad ---------------------------------------------

    @property
    def hz(self) -> float:
        """H(z) [km/s/Mpc]. Universo plano materia + Lambda, como en Sherwood."""
        return 100.0 * self.hubble * np.sqrt(
            self.omega_m * (1.0 + self.z) ** 3 + self.omega_l
        )

    @property
    def box_kms(self) -> float:
        """Tamano de la caja en velocidad de Hubble [km/s]."""
        box_pmpc = 1.0e-3 * self.boxsize / (self.hubble * (1.0 + self.z))
        return self.hz * box_pmpc

    @property
    def dv(self) -> float:
        """Ancho de pixel [km/s]."""
        return self.box_kms / self.npix

    @property
    def velocity(self) -> np.ndarray:
        """Eje de velocidad de Hubble por pixel [km/s]."""
        pix_pmpc = 1.0e-3 * self.pixels / (self.hubble * (1.0 + self.z))
        return self.hz * pix_pmpc

    @property
    def wavelength(self) -> np.ndarray:
        """Longitud de onda observada [Angstrom] (relativista, como en el IDL)."""
        w = self.velocity / CKMS
        return LAMBDA_LYA * (1.0 + self.z) * np.sqrt((1.0 + w) / (1.0 - w))

    def los_position(self, i: int) -> str:
        """Equivalente a get_los_position.pro."""
        names = {1: ("x", "y", "z"), 2: ("y", "x", "z"), 3: ("z", "x", "y")}
        along, c1, c2 = names[int(self.axis[i])]
        return (
            f"LOS {i} corre a lo largo del eje {along}: "
            f"{c1} = {self.coord1[i]:.3f}, {c2} = {self.coord2[i]:.3f} [h^-1 ckpc]"
        )


def _expected_size(npix: int, nlos: int) -> int:
    return HEADER_BYTES + 4 * (3 * nlos + npix + npix * nlos)


def read_sherwood_spectra(filename: str, mmap: bool = False) -> SherwoodSpectra:
    """
    Lee un archivo tipo `tauH1_lya_z2.0.dat`.

    El endianness se detecta solo: se prueban los dos y se elige aquel para el
    cual el tamano del archivo coincide exactamente con el implicado por el
    header. Si ninguno cierra, el archivo no tiene este formato y se levanta
    ValueError en vez de devolver basura silenciosamente.

    mmap=True evita cargar `tau` en RAM (util para cajas grandes o para leer
    solo un subconjunto de LOS); el array devuelto es de solo lectura.
    """
    fsize = os.path.getsize(filename)
    if fsize <= HEADER_BYTES:
        raise ValueError(f"{filename}: archivo demasiado chico ({fsize} bytes).")

    with open(filename, "rb") as fh:
        raw = fh.read(8)

    endian = None
    for e in ("<", ">"):
        npix, nlos = np.frombuffer(raw, dtype=np.dtype(f"{e}i4"), count=2)
        npix, nlos = int(npix), int(nlos)
        if npix > 0 and nlos > 0 and _expected_size(npix, nlos) == fsize:
            endian = e
            break

    if endian is None:
        le = np.frombuffer(raw, dtype="<i4", count=2)
        be = np.frombuffer(raw, dtype=">i4", count=2)
        raise ValueError(
            f"{filename}: el tamano ({fsize} bytes) no coincide con ningun header "
            f"valido. Leido little-endian: npix={le[0]}, nlos={le[1]}; "
            f"big-endian: npix={be[0]}, nlos={be[1]}. "
            "Revisa que sea el formato binario de Sherwood."
        )

    i4 = np.dtype(f"{endian}i4")
    f4 = np.dtype(f"{endian}f4")

    with open(filename, "rb") as fh:
        fh.seek(8)
        head = np.fromfile(fh, dtype=f4, count=6)
        z, om, ol, ob, h, box = (float(v) for v in head)

        axis = np.fromfile(fh, dtype=i4, count=nlos)
        coord1 = np.fromfile(fh, dtype=f4, count=nlos)
        coord2 = np.fromfile(fh, dtype=f4, count=nlos)
        pixels = np.fromfile(fh, dtype=f4, count=npix)
        offset_tau = fh.tell()

        if mmap:
            tau = np.memmap(
                filename, dtype=f4, mode="r", offset=offset_tau,
                shape=(nlos, npix),
            )
        else:
            tau = np.fromfile(fh, dtype=f4, count=npix * nlos).reshape(nlos, npix)

    return SherwoodSpectra(
        z=z, omega_m=om, omega_l=ol, omega_b=ob, hubble=h, boxsize=box,
        axis=axis, coord1=coord1, coord2=coord2, pixels=pixels, tau=tau,
    )


# ---------------------------------------------------------------------------
# UVB (get_uvb.pro)
# ---------------------------------------------------------------------------

def gamma_hi(treecool_file: str, z: float) -> float:
    """
    Tasa de fotoionizacion de HI [s^-1] interpolada de una tabla TREECOOL.
    Columna 0 = log10(1+z), columna 1 = Gamma_HI.

    Nota: IDL `interpol` extrapola linealmente fuera del rango; np.interp
    satura en los extremos. Para las tablas HM12 (z hasta ~15) no importa
    en el rango del forest.
    """
    table = np.loadtxt(treecool_file, usecols=(0, 1))
    zz = 10.0 ** table[:, 0] - 1.0
    return float(np.interp(z, zz, table[:, 1]))


# ---------------------------------------------------------------------------
# Flujo y P1D
# ---------------------------------------------------------------------------

def tau_eff(tau: np.ndarray, scale: float = 1.0) -> float:
    """Profundidad optica efectiva: -ln(<exp(-A*tau)>) sobre todos los pixeles."""
    return float(-np.log(np.mean(np.exp(-scale * np.asarray(tau, dtype=np.float64)))))


def rescale_tau(tau: np.ndarray, target_tau_eff: float,
                tol: float = 1e-8, max_iter: int = 100) -> float:
    """
    Devuelve el factor A tal que -ln(<exp(-A*tau)>) = target_tau_eff.

    Reescalar tau por una constante es el procedimiento estandar para absorber
    la incerteza en Gamma_HI (equivale a Gamma_HI -> Gamma_HI / A). Hacelo
    SIEMPRE antes de comparar P1D con datos, porque el P1D es muy sensible al
    flujo medio. tau_eff es monotona creciente en A, asi que biseccion alcanza.
    """
    t = np.asarray(tau, dtype=np.float64)

    lo, hi = 1.0e-3, 1.0e3
    if tau_eff(t, lo) > target_tau_eff or tau_eff(t, hi) < target_tau_eff:
        raise ValueError(
            f"target_tau_eff={target_tau_eff} fuera del rango alcanzable "
            f"[{tau_eff(t, lo):.4f}, {tau_eff(t, hi):.4f}]."
        )

    for _ in range(max_iter):
        mid = np.sqrt(lo * hi)
        if tau_eff(t, mid) < target_tau_eff:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1.0 < tol:
            break
    return float(np.sqrt(lo * hi))


def flux_power_1d(flux: np.ndarray, dv: float, mean_flux: float | None = None,
                  deconvolve_pixel: bool = True, chunk: int = 4096):
    """
    Espectro de potencias 1D del flujo transmitido.

    Parametros
    ----------
    flux : (nlos, npix) o (npix,)
        F = exp(-tau).
    dv : float
        Ancho de pixel [km/s].
    mean_flux : float, opcional
        <F> usado para definir delta_F = F/<F> - 1. Por defecto la media global
        de la muestra (NO por LOS: normalizar LOS por LOS te borra potencia en
        las escalas grandes).
    deconvolve_pixel : bool
        Divide por el cuadrado de la ventana de pixel top-hat, W(k) =
        sinc(k*dv/2). Corrige el suavizado por binning finito; a k pequeno es
        despreciable, cerca de Nyquist llega a ~20%.
    chunk : int
        LOS por bloque de FFT. Mantiene el pico de memoria acotado sin perder
        velocidad; la FFT ya es el cuello de botella.

    Devuelve
    --------
    k : (npix//2 + 1,) [s/km]
    p1d : (npix//2 + 1,) [km/s]
    p1d_err : (npix//2 + 1,) [km/s]
        Error estandar de la media sobre LOS. Ojo: las LOS de una misma caja no
        son independientes, asi que esto subestima la varianza real.
    """
    f = np.asarray(flux)
    if f.ndim == 1:
        f = f[None, :]
    nlos, npix = f.shape

    if mean_flux is None:
        mean_flux = float(np.mean(f, dtype=np.float64))
    if not np.isfinite(mean_flux) or mean_flux <= 0:
        raise ValueError(f"mean_flux invalido: {mean_flux}")

    k = 2.0 * np.pi * np.fft.rfftfreq(npix, d=dv)
    nk = k.size
    length = npix * dv

    # Acumuladores en float64: media y suma de cuadrados por bloque.
    acc = np.zeros(nk, dtype=np.float64)
    acc2 = np.zeros(nk, dtype=np.float64)

    for start in range(0, nlos, chunk):
        block = np.asarray(f[start:start + chunk], dtype=np.float64)
        delta = block / mean_flux - 1.0
        pk = (length / npix ** 2) * np.abs(np.fft.rfft(delta, axis=1)) ** 2
        acc += pk.sum(axis=0)
        acc2 += (pk ** 2).sum(axis=0)

    p1d = acc / nlos
    var = np.maximum(acc2 / nlos - p1d ** 2, 0.0)
    p1d_err = np.sqrt(var / max(nlos - 1, 1))

    if deconvolve_pixel:
        # np.sinc(x) = sin(pi x)/(pi x)  =>  sin(k dv/2)/(k dv/2) = sinc(k dv / 2pi)
        w2 = np.sinc(k * dv / (2.0 * np.pi)) ** 2
        p1d = p1d / w2
        p1d_err = p1d_err / w2

    return k, p1d, p1d_err
