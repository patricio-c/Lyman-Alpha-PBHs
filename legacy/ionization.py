"""
Equilibrio de fotoionizacion: densidad total de gas + temperatura -> n_HI.

Este es el paso que falta entre el output de LOS de SWIFT (que da Densities y
Temperatures) y la extraccion de tau. Las corridas sin tablas de especies
quimicas no guardan la fraccion neutra, asi que hay que reconstruirla en
post-proceso a partir de la tasa de fotoionizacion del fondo UV.

Aproximaciones y cuando valen:

  - Caso A para la recombinacion. Correcto para el IGM opticamente delgado,
    que es todo el forest a densidades relevantes. Si te interesan sistemas
    autoblindados (N_HI > 1e17) esto subestima la fraccion neutra y necesitas
    un modelo de self-shielding tipo Rahmati et al. (2013).
  - Equilibrio de ionizacion instantaneo. Falla durante la reionizacion; a
    z < 6 con un fondo uniforme es buena.
  - Helio: solo entra via su contribucion a n_e. Se parametriza con el estado
    de ionizacion del He, que a z ~ 4 es HeII.

Referencias de los coeficientes: Hui & Gnedin (1997) para la recombinacion
caso A, Theuns et al. (1998) para la ionizacion colisional.
"""

from __future__ import annotations

import numpy as np

M_H = 1.67262192e-24     # masa del proton [g]
T_HI = 157807.0          # 13.6 eV / k_B [K]


# ---------------------------------------------------------------------------
# Coeficientes atomicos
# ---------------------------------------------------------------------------

def alpha_A_HII(T: np.ndarray) -> np.ndarray:
    """Recombinacion caso A de HII [cm^3/s]. Hui & Gnedin (1997), ec. A4."""
    T = np.asarray(T, dtype=np.float64)
    lam = 2.0 * T_HI / T
    return (1.269e-13 * lam ** 1.503
            / (1.0 + (lam / 0.522) ** 0.470) ** 1.923)


def alpha_B_HII(T: np.ndarray) -> np.ndarray:
    """Recombinacion caso B de HII [cm^3/s]. Hui & Gnedin (1997), ec. A6."""
    T = np.asarray(T, dtype=np.float64)
    lam = 2.0 * T_HI / T
    return (2.753e-14 * lam ** 1.500
            / (1.0 + (lam / 2.740) ** 0.407) ** 2.242)


def gamma_coll_HI(T: np.ndarray) -> np.ndarray:
    """Ionizacion colisional de HI [cm^3/s]. Theuns et al. (1998), tabla."""
    T = np.asarray(T, dtype=np.float64)
    # exp(-T_HI/T) subdesborda a T bajas; se acota para evitar warnings.
    arg = np.clip(T_HI / T, 0.0, 700.0)
    return (1.17e-10 * np.sqrt(T) * np.exp(-arg)
            / (1.0 + np.sqrt(T / 1.0e5)))


# ---------------------------------------------------------------------------
# Solucion del equilibrio
# ---------------------------------------------------------------------------

def neutral_fraction(n_H: np.ndarray, T: np.ndarray, gamma_HI: float,
                     X_H: float = 0.76, he_state: str = "HeII",
                     case: str = "A", collisional: bool = True,
                     n_iter: int = 8) -> np.ndarray:
    """
    Fraccion neutra de hidrogeno f = n_HI / n_H en equilibrio.

    Parametros
    ----------
    n_H : densidad numerica de hidrogeno TOTAL (neutro + ionizado) [cm^-3]
    T : temperatura [K]
    gamma_HI : tasa de fotoionizacion de HI [s^-1], de la tabla TREECOOL
    X_H : fraccion de masa de hidrogeno. 0.76 es primordial; si tu corrida
          tiene enriquecimiento y guardaste ElementMassFractions, usa esa.
    he_state : 'HeII' o 'HeIII'. Solo afecta n_e. A z ~ 4, antes de la
          reionizacion de HeII, el helio esta mayormente en HeII.
    case : 'A' (opticamente delgado, lo normal para el forest) o 'B'.
    collisional : incluir ionizacion colisional. Importa solo en gas chocado
          a T > 1e4.5 K; en el IGM fotoionizado es despreciable pero barato.

    Balance resuelto, por atomo de hidrogeno:

        n_HI * (gamma_HI + q_coll * n_e)  =  alpha * n_e * n_HII

    Con n_e = n_HII + (contribucion del He), que depende de n_HII. Se resuelve
    por iteracion de punto fijo: se arranca suponiendo n_e = n_H (gas casi
    totalmente ionizado, cierto en el forest) y converge en pocas vueltas.
    """
    n_H = np.asarray(n_H, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)

    if np.any(n_H < 0) or np.any(T <= 0):
        raise ValueError("n_H negativa o T no positiva en la entrada.")

    alpha = alpha_A_HII(T) if case.upper() == "A" else alpha_B_HII(T)
    q = gamma_coll_HI(T) if collisional else np.zeros_like(T)

    # Electrones aportados por el helio, por unidad de n_H.
    # n_He / n_H = (1 - X_H) / (4 X_H)
    y = (1.0 - X_H) / (4.0 * X_H)
    he_e = y * (1.0 if he_state == "HeII" else 2.0)

    # Punto fijo. x = n_HII / n_H.
    x = np.ones_like(n_H)
    for _ in range(n_iter):
        n_e = n_H * (x + he_e)
        # (1 - x) * (gamma + q*n_e) = alpha * n_e * x
        num = gamma_HI + q * n_e
        x_new = num / (num + alpha * n_e)
        # Donde n_e -> 0 el gas queda neutro por falta de recombinaciones;
        # el limite correcto es x -> 1 si gamma > 0.
        x = np.where(n_e > 0, x_new, 1.0)

    return np.clip(1.0 - x, 0.0, 1.0)


def n_HI_from_density(rho_gas: np.ndarray, T: np.ndarray, gamma_HI: float,
                      X_H: float = 0.76, **kwargs) -> np.ndarray:
    """
    Densidad numerica de HI [cm^-3] desde densidad de masa de gas [g/cm^3].

    rho_gas debe estar en unidades FISICAS (proper), no comoviles. SWIFT
    guarda densidades comoviles: convertilas con rho_fisica = rho_com * (1+z)^3
    antes de llamar a esto, o te vas a equivocar por un factor 140 a z = 4.2.
    """
    n_H = np.asarray(rho_gas, dtype=np.float64) * X_H / M_H
    return n_H * neutral_fraction(n_H, T, gamma_HI, X_H=X_H, **kwargs)


def mean_hydrogen_density(z: float, omega_b: float, h: float,
                          X_H: float = 0.76) -> float:
    """
    Densidad numerica media de hidrogeno [cm^-3] a redshift z, en fisicas.
    Util para expresar resultados en terminos de sobredensidad Delta = n_H/<n_H>.
    """
    rho_crit0 = 1.87847e-29 * h ** 2          # [g/cm^3]
    return rho_crit0 * omega_b * (1.0 + z) ** 3 * X_H / M_H
