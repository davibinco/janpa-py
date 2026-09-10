"""Gaussian-type orbital overlap integrals. Replaces OverlapIntegrals.java."""
from __future__ import annotations

import numpy as np
from scipy.special import comb

from .polynom3d import Polynom3D


def gamma_hi(n: int) -> float:
    """Γ(n + 1/2) = √π * (1/2)(3/2)...(n-1/2)."""
    result = np.sqrt(np.pi)
    for i in range(1, n + 1):
        result *= 0.5 + (i - 1)
    return result


def primitive_int_1d(n: int, alpha: float) -> float:
    """∫_{-∞}^{∞} x^n exp(-α x²) dx. Zero for odd n."""
    if n % 2 == 1:
        return 0.0
    return alpha ** (-0.5 - n / 2.0) * gamma_hi(n // 2)


def primitive_int_1d_sphr(n: int, alpha: float) -> float:
    """∫_0^{∞} r^n exp(-α r²) dr (radial integral)."""
    if n % 2 == 0:
        return primitive_int_1d(n, alpha) / 2.0
    result = 0.5 / alpha
    for i in range(1, (n - 1) // 2 + 1):
        result *= i / alpha
    return result


# Precompute C_nk table up to max expected power (L_MAX = 10 -> 2*L_MAX = 20)
_MAX_PWR = 20
_QUICK_CNK = np.zeros((_MAX_PWR + 1, _MAX_PWR + 1))
for _i in range(_MAX_PWR + 1):
    for _j in range(_i + 1):
        _QUICK_CNK[_i, _j] = comb(_i, _j, exact=True)


def bs_bs_overlap(
    zeta1: float,
    ylm1: Polynom3D,
    r01: np.ndarray,
    zeta2: float,
    ylm2: Polynom3D,
    r02: np.ndarray,
) -> float:
    """Overlap of two primitive Gaussian basis functions with angular parts ylm1, ylm2."""
    zeta_eff = zeta1 + zeta2
    r0eff = (zeta1 * r01 + zeta2 * r02) / zeta_eff
    prefactor = np.exp(-zeta1 * zeta2 * np.sum((r01 - r02)**2) / zeta_eff)
    
    nmax = 0
    for i in range(ylm1.n_terms):
        for j in range(ylm2.n_terms):
            for mu in range(3):
                ncurr = ylm1.powers[i, mu] + ylm2.powers[j, mu]
                if ncurr > nmax:
                    nmax = ncurr
                    
    quick_ints = np.array([primitive_int_1d(i, zeta_eff) for i in range(nmax + 1)])
    
    quick_expans_1 = np.ones((nmax + 1, 3))
    quick_expans_2 = np.ones((nmax + 1, 3))
    
    for pwr in range(1, nmax + 1):
        quick_expans_1[pwr] = quick_expans_1[pwr - 1] * (r0eff - r01)
        quick_expans_2[pwr] = quick_expans_2[pwr - 1] * (r0eff - r02)
        
    result = 0.0
    for i in range(ylm1.n_terms):
        for j in range(ylm2.n_terms):
            term = ylm1.coefs[i] * ylm2.coefs[j]
            for mu in range(3):
                tmp = 0.0
                pow1 = ylm1.powers[i, mu]
                pow2 = ylm2.powers[j, mu]
                for pwr1 in range(pow1 + 1):
                    c1 = _QUICK_CNK[pow1, pwr1] * quick_expans_1[pwr1, mu]
                    for pwr2 in range(pow2 + 1):
                        c2 = _QUICK_CNK[pow2, pwr2] * quick_expans_2[pwr2, mu]
                        idx = (pow1 - pwr1) + (pow2 - pwr2)
                        tmp += c1 * c2 * quick_ints[idx]
                term *= tmp
                if tmp == 0:
                    break
            result += term
            
    return float(result * prefactor)


def same_center_radial_overlap(
    coefs1: np.ndarray,
    exponents1: np.ndarray,
    l1: int,
    coefs2: np.ndarray,
    exponents2: np.ndarray,
    l2: int,
) -> float:
    """Same-center radial overlap integral (angular parts integrated analytically)."""
    result = 0.0
    for c1, e1 in zip(coefs1, exponents1):
        for c2, e2 in zip(coefs2, exponents2):
            result += c1 * c2 * primitive_int_1d_sphr(l1 + l2 + 2, e1 + e2)
    return result * 4 * np.pi
