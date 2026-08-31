"""
Post-proceso de los espectros sinteticos siguiendo Bolton et al. (2017),
seccion 3.2. Estos son los pasos que van DESPUES de tener tau y ANTES de
medir el P1D.

Orden correcto (no lo cambies, no conmutan):

    tau  ->  reescalar a tau_eff        (sherwood_los.rescale_tau)
         ->  F = exp(-A*tau)
         ->  convolve_lsf(F, fwhm=7 km/s)
         ->  rebin_velocity(F, dv_out=3 km/s)
         ->  add_noise(F, snr=50)            [opcional]
         ->  apply_continuum_bias(F, z)      [opcional]
         ->  flux_power_1d(...)

Por que no conmutan: convolucionar despues de rebinnear suaviza sobre una
grilla ya degradada y pierde informacion; agregar ruido antes de convolucionar
lo correlaciona artificialmente (el ruido instrumental entra despues de la
optica, no antes).
"""

from __future__ import annotations

import numpy as np

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))  # 0.42466


# ---------------------------------------------------------------------------
# 1. Perfil instrumental
# ---------------------------------------------------------------------------

def convolve_lsf(flux: np.ndarray, dv: float, fwhm_kms: float = 7.0,
                 periodic: bool = True) -> np.ndarray:
    """
    Convoluciona con una gaussiana de ancho a media altura `fwhm_kms`.

    Se hace por FFT y no por np.convolve porque las LOS de caja son periodicas:
    la convolucion circular es la condicion de borde fisicamente correcta, y de
    paso es O(N log N) en vez de O(N*M). Con (5000, 2048) la diferencia se nota.

    periodic=False hace padding con reflexion, para LOS que no son de caja
    periodica (por ejemplo un light-cone armado pegando snapshots).
    """
    f = np.atleast_2d(np.asarray(flux, dtype=np.float64))
    npix = f.shape[1]
    sigma = fwhm_kms * FWHM_TO_SIGMA

    if sigma <= 0:
        return flux

    if sigma < dv / 2:
        raise ValueError(
            f"sigma={sigma:.3f} km/s es menor que medio pixel (dv={dv:.3f}). "
            "El kernel esta submuestreado; la convolucion no significa nada."
        )

    if periodic:
        # Kernel gaussiano en el espacio de Fourier: exp(-k^2 sigma^2 / 2).
        # Analitico, sin errores de discretizacion del kernel.
        k = 2.0 * np.pi * np.fft.rfftfreq(npix, d=dv)
        win = np.exp(-0.5 * (k * sigma) ** 2)
        out = np.fft.irfft(np.fft.rfft(f, axis=1) * win, n=npix, axis=1)
    else:
        half = int(np.ceil(4.0 * sigma / dv))
        x = np.arange(-half, half + 1) * dv
        kern = np.exp(-0.5 * (x / sigma) ** 2)
        kern /= kern.sum()
        pad = np.pad(f, ((0, 0), (half, half)), mode="reflect")
        out = np.empty_like(f)
        for i in range(f.shape[0]):
            out[i] = np.convolve(pad[i], kern, mode="valid")

    return out if np.ndim(flux) > 1 else out[0]


# ---------------------------------------------------------------------------
# 2. Rebinning
# ---------------------------------------------------------------------------

def rebin_velocity(flux: np.ndarray, dv_in: float, dv_out: float = 3.0):
    """
    Rebinnea de una grilla de dv_in a una de dv_out conservando la integral.

    Devuelve (flux_rebin, dv_real). Ojo: dv_out casi nunca divide exacto la
    caja, asi que se usa nbins = round(L/dv_out) y el ancho real devuelto
    (dv_real = L/nbins) difiere un poco del pedido. Usa SIEMPRE dv_real para
    el P1D; si usas el nominal te desplazas la escala de k.

    Implementacion: se interpola la integral acumulada del flujo en los bordes
    de los bins nuevos y se diferencia. Esto conserva exactamente el flujo
    total, cosa que interpolar el flujo directamente no hace.
    """
    f = np.atleast_2d(np.asarray(flux, dtype=np.float64))
    nlos, npix = f.shape
    length = npix * dv_in

    nbins = int(round(length / dv_out))
    if nbins < 8:
        raise ValueError(f"dv_out={dv_out} deja solo {nbins} bins. Demasiado grueso.")
    if nbins > npix:
        raise ValueError(
            f"dv_out={dv_out} km/s es mas fino que la grilla nativa "
            f"({dv_in:.4f} km/s). Rebinnear hacia arriba inventa informacion."
        )
    dv_real = length / nbins

    # Bordes de los bins de entrada: 0, dv_in, 2*dv_in, ...
    edges_in = np.arange(npix + 1) * dv_in
    edges_out = np.arange(nbins + 1) * dv_real

    # Integral acumulada del flujo, con un 0 al principio.
    cum = np.zeros((nlos, npix + 1), dtype=np.float64)
    np.cumsum(f * dv_in, axis=1, out=cum[:, 1:])

    # np.interp no vectoriza sobre filas; se hace por bloque de columnas.
    idx = np.searchsorted(edges_in, edges_out, side="right") - 1
    idx = np.clip(idx, 0, npix - 1)
    frac = (edges_out - edges_in[idx]) / dv_in

    cum_at = cum[:, idx] + frac * (cum[:, idx + 1] - cum[:, idx])
    out = np.diff(cum_at, axis=1) / dv_real

    return (out if np.ndim(flux) > 1 else out[0]), dv_real


# ---------------------------------------------------------------------------
# 3. Ruido
# ---------------------------------------------------------------------------

def add_noise(flux: np.ndarray, snr: float = 50.0, seed: int | None = None):
    """
    Agrega ruido gaussiano uniforme de sigma = 1/snr por pixel.

    Devuelve (flux_ruidoso, sigma).

    ATENCION: esto es forward-modelling de la observacion, no dato simulado.
    Solo tiene sentido si vas a comparar contra un P1D observado que TAMPOCO
    tiene sustraido su ruido. Si compras contra un P1D ya corregido por ruido
    (que es lo normal en DESI y eBOSS), NO agregues ruido aca: agregas un piso
    blanco P_N = sigma^2 * dv / <F>^2 que a k grande domina todo.

    Si lo agregas igual, restale ese piso con noise_power() antes de comparar.
    """
    f = np.asarray(flux, dtype=np.float64)
    sigma = 1.0 / snr
    rng = np.random.default_rng(seed)
    return f + rng.normal(0.0, sigma, size=f.shape), sigma


def noise_power(sigma: float, dv: float, mean_flux: float) -> float:
    """
    Piso de potencia blanca que introduce add_noise, en km/s.

    delta = F/<F> - 1  =>  var(ruido en delta) = sigma^2 / <F>^2
    Ruido blanco: P_N = var * dv, plano en k.
    """
    return (sigma / mean_flux) ** 2 * dv


# ---------------------------------------------------------------------------
# 4. Sesgo de continuo
# ---------------------------------------------------------------------------

def continuum_bias(z: float, a: float = 1.58e-5, b: float = 5.63) -> float:
    """
    C_corr = C_est / C_true, la razon entre el continuo estimado a mano sobre
    los datos y el verdadero (Faucher-Giguere et al. 2008b).

    Forma funcional C_corr = 1 - a*(1+z)^b. Da ~2% a z=2.5 y ~17% a z=4.2,
    que son los valores que cita la literatura para FG08b.

    OJO: verifica los coeficientes contra la Eq. (2) de Bolton et al. (2017)
    antes de usar esto en un paper. Esta parametrizacion la reconstrui de como
    la citan otros trabajos, no de la ecuacion impresa en el paper de Sherwood.
    Ademas FG08b la calibraron para 2 <= z <= 4; a z=4.2 ya estas extrapolando,
    y a z=5 la extrapolacion se va a ~40%, que no es creible.
    """
    c = 1.0 - a * (1.0 + z) ** b
    if c <= 0:
        raise ValueError(f"C_corr = {c:.3f} en z={z}. La extrapolacion se rompio.")
    return c


def apply_continuum_bias(flux: np.ndarray, z: float, **kwargs) -> np.ndarray:
    """
    Emula un continuo mal colocado dividiendo por C_corr.

    Direccion: si el continuo se coloca demasiado bajo (C_est < C_true, o sea
    C_corr < 1), el flujo normalizado que mide el observador es MAS ALTO que el
    verdadero, porque divide por un continuo mas chico. Por eso se divide.

    Esto sube <F>, o sea baja tau_eff. Es el paso que mas mueve el P1D a z>4.
    """
    return np.asarray(flux, dtype=np.float64) / continuum_bias(z, **kwargs)
