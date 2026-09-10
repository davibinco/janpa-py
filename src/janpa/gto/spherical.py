"""Hardcoded real solid harmonic polynomials up to L=4 (G functions), normalization tables, and Cartesian-to-Pure conversion coefficients. Replaces SphericalHarmonics.java. All data tables are module-level constants. """

from __future__ import annotations

import copy
import numpy as np

from .polynom3d import Polynom3D, YLMSeries
from .overlap import primitive_int_1d, primitive_int_1d_sphr, bs_bs_overlap

L_MAX = 4

def get_quick_ylm() -> list[list[Polynom3D]]:
    """Return Quick_YLM[L][L+m] as Polynom3D objects, up to L=4."""
    S3 = np.sqrt(3); S15 = np.sqrt(15); S10 = np.sqrt(10)
    S6 = np.sqrt(6); S7 = np.sqrt(7); S35 = np.sqrt(35)
    S350 = np.sqrt(350); S2 = np.sqrt(2); S14 = np.sqrt(14)

    Y = [[None] * (2 * L + 1) for L in range(5)]

    Y[0][0] = Polynom3D(np.array([1.0]), np.array([[0, 0, 0]]))

    Y[1][0] = Polynom3D(np.array([1.0]), np.array([[0, 1, 0]]))
    Y[1][1] = Polynom3D(np.array([1.0]), np.array([[0, 0, 1]]))
    Y[1][2] = Polynom3D(np.array([1.0]), np.array([[1, 0, 0]]))

    Y[2][0] = Polynom3D(np.array([1.0]), np.array([[1, 1, 0]]))
    Y[2][1] = Polynom3D(np.array([1.0]), np.array([[0, 1, 1]]))
    Y[2][2] = Polynom3D(
        np.array([1.0 / S3, -0.5 / S3, -0.5 / S3]),
        np.array([[0, 0, 2], [0, 2, 0], [2, 0, 0]]),
    )
    Y[2][3] = Polynom3D(np.array([1.0]), np.array([[1, 0, 1]]))
    Y[2][4] = Polynom3D(
        np.array([0.5, -0.5]),
        np.array([[2, 0, 0], [0, 2, 0]]),
    )

    Y[3][0] = Polynom3D(
        np.array([1.5 / S6, -0.5 / S6]),
        np.array([[2, 1, 0], [0, 3, 0]]),
    )
    Y[3][1] = Polynom3D(np.array([1.0]), np.array([[1, 1, 1]]))
    Y[3][2] = Polynom3D(
        np.array([2.0 / S10, -0.5 / S10, -0.5 / S10]),
        np.array([[0, 1, 2], [0, 3, 0], [2, 1, 0]]),
    )
    Y[3][3] = Polynom3D(
        np.array([1.0 / S15, -1.5 / S15, -1.5 / S15]),
        np.array([[0, 0, 3], [0, 2, 1], [2, 0, 1]]),
    )
    Y[3][4] = Polynom3D(
        np.array([2.0 / S10, -0.5 / S10, -0.5 / S10]),
        np.array([[1, 0, 2], [1, 2, 0], [3, 0, 0]]),
    )
    Y[3][5] = Polynom3D(
        np.array([0.5, -0.5]),
        np.array([[2, 0, 1], [0, 2, 1]]),
    )
    Y[3][6] = Polynom3D(
        np.array([0.5 / S6, -1.5 / S6]),
        np.array([[3, 0, 0], [1, 2, 0]]),
    )

    Y[4][0] = Polynom3D(
        np.array([0.5, -0.5]),
        np.array([[3, 1, 0], [1, 3, 0]]),
    )
    Y[4][1] = Polynom3D(
        np.array([1.5 / S2, -0.5 / S2]),
        np.array([[2, 1, 1], [0, 3, 1]]),
    )
    Y[4][2] = Polynom3D(
        np.array([3.0 / S7, -0.5 / S7, -0.5 / S7]),
        np.array([[1, 1, 2], [1, 3, 0], [3, 1, 0]]),
    )
    Y[4][3] = Polynom3D(
        np.array([10.0 / S350, -7.5 / S350, -7.5 / S350]),
        np.array([[0, 1, 3], [0, 3, 1], [2, 1, 1]]),
    )
    Y[4][4] = Polynom3D(
        np.array([1.0 / S35, -3.0 / S35, -3.0 / S35, 0.375 / S35, 0.75 / S35, 0.375 / S35]),
        np.array([[0, 0, 4], [0, 2, 2], [2, 0, 2], [0, 4, 0], [2, 2, 0], [4, 0, 0]]),
    )
    Y[4][5] = Polynom3D(
        np.array([10.0 / S350, -7.5 / S350, -7.5 / S350]),
        np.array([[1, 0, 3], [1, 2, 1], [3, 0, 1]]),
    )
    Y[4][6] = Polynom3D(
        np.array([1.5 / S7, -1.5 / S7, 0.25 / S7, -0.25 / S7]),
        np.array([[2, 0, 2], [0, 2, 2], [0, 4, 0], [4, 0, 0]]),
    )
    Y[4][7] = Polynom3D(
        np.array([0.5 / S2, -1.5 / S2]),
        np.array([[3, 0, 1], [1, 2, 1]]),
    )
    Y[4][8] = Polynom3D(
        np.array([0.125, -0.750, 0.125]),
        np.array([[4, 0, 0], [2, 2, 0], [0, 4, 0]]),
    )

    return Y

MOLDEN_CART_POWERS = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
    [1, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [3, 0, 0],
    [0, 3, 0],
    [0, 0, 3],
    [1, 2, 0],
    [2, 1, 0],
    [2, 0, 1],
    [1, 0, 2],
    [0, 1, 2],
    [0, 2, 1],
    [1, 1, 1],
    [4, 0, 0],
    [0, 4, 0],
    [0, 0, 4],
    [3, 1, 0],
    [3, 0, 1],
    [1, 3, 0],
    [0, 3, 1],
    [1, 0, 3],
    [0, 1, 3],
    [2, 2, 0],
    [2, 0, 2],
    [0, 2, 2],
    [2, 1, 1],
    [1, 2, 1],
    [1, 1, 2],
], dtype=int)

ALEG_ROOTS_SQUARED = [
    [[]],
    [[], []],
    [[1.0 / 3.0], [], []],
    [[3.0 / 5.0], [1.0 / 5.0], [], []],
    [[0.11558710999704793517, 0.74155574714580920772], [3.0 / 7.0], [1.0 / 7.0], [], []],
    [[0.28994919792569030224, 0.82116191318542080892], [0.08135701799384851519, 0.58530964867281815149], [1.0 / 3.0], [1.0 / 9.0], [], []],
]

def get_aleg_normalizer() -> list[list[float]]:
    """Compute ALegNormalizer[L][L+m] for L=0..5, m=-L..L."""
    n_l = len(ALEG_ROOTS_SQUARED)
    result: list[list[float]] = [[] for _ in range(n_l)]
    for L in range(n_l):
        result[L] = [0.0] * (2 * L + 1)
        tmp = 1.0
        for i in range(1, L + 1):
            tmp *= 2.0 - 1.0 / i
        result[L][L] = tmp * np.sqrt((2.0 * L + 1.0) / (4.0 * np.pi))
        for m in range(1, L + 1):
            result[L][L + m] = result[L][L + m - 1] / np.sqrt((L + m) / (L - m + 1.0))
            if m == 1:
                result[L][L + m] *= np.sqrt(2.0)
            result[L][L - m] = result[L][L + m]
    return result

ALEG_NORMALIZER = get_aleg_normalizer()

def get_quick_ylm_norm2(quick_ylm: list[list[Polynom3D]]) -> list[list[float]]:
    """Compute squared norms of each YLM component (over 4π)."""
    result = []
    for L in range(len(quick_ylm)):
        result_L = []
        radial_norm2 = primitive_int_1d_sphr(2 + 2 * L, 2.0) * 4 * np.pi
        for m in range(len(quick_ylm[L])):
            ylm = quick_ylm[L][m]
            overlap = bs_bs_overlap(
                1.0, ylm, np.array([0, 0, 0], dtype=float),
                1.0, ylm, np.array([0, 0, 0], dtype=float)
            )
            result_L.append(overlap / radial_norm2)
        result.append(result_L)
    return result

def get_quick_ylm_norm4pi() -> list[list[Polynom3D]]:
    """Return Quick_YLM normalized so that each component has ∫Y²dΩ = 4π."""
    result = copy.deepcopy(get_quick_ylm())
    norms2 = get_quick_ylm_norm2(result)
    for L in range(len(result)):
        for m in range(len(result[L])):
            result[L][m].scale_coefs_by(1.0 / np.sqrt(norms2[L][m]))
    return result

def cartesian_to_pure() -> list[YLMSeries]:
    """Return the 35 Cartesian-to-Pure conversion YLM_Series entries.
    Entry index matches MOLDEN_CART_POWERS ordering.
    """
    S2 = np.sqrt(2)
    S3 = np.sqrt(3)
    S7 = np.sqrt(7)
    S10 = np.sqrt(10)
    S14 = np.sqrt(14)
    S15 = np.sqrt(15)
    S35 = np.sqrt(35)
    
    def Y(l, lm_pairs, coefs):
        ls = np.array([p[0] for p in lm_pairs], dtype=int)
        ms = np.array([p[1] for p in lm_pairs], dtype=int)
        return YLMSeries(primary_l=l, ls=ls, ms=ms, coefs=np.array(coefs))
        
    result = [None] * 35
    result[ 0] = Y(0, [(0, 0)], [1.0])
    result[ 1] = Y(1, [(1, 0)], [1.0])
    result[ 2] = Y(1, [(1, 1)], [1.0])
    result[ 3] = Y(1, [(1, -1)], [1.0])
    result[ 4] = Y(2, [(2, 0), (2, 2), (0, 0)], [-1.0/S3, 1.0, 1.0/3.0])
    result[ 5] = Y(2, [(2, 0), (2, 2), (0, 0)], [-1.0/S3, -1.0, 1.0/3.0])
    result[ 6] = Y(2, [(2, 0), (0, 0)], [2.0/S3, 1.0/3.0])
    result[ 7] = Y(2, [(2, -2)], [1.0])
    result[ 8] = Y(2, [(2, 1)], [1.0])
    result[ 9] = Y(2, [(2, -1)], [1.0])
    result[10] = Y(3, [(3, 1), (3, 3), (1, 0)], [-3.0/S10, S3/S2, 0.6])
    result[11] = Y(3, [(3, -3), (3, -1), (1, 1)], [-S3/S2, -3.0/S10, 0.6])
    result[12] = Y(3, [(3, 0), (1, -1)], [6.0/S15, 0.6])
    result[13] = Y(3, [(3, 1), (3, 3), (1, 0)], [-1.0/S10, -S3/S2, 0.2])
    result[14] = Y(3, [(3, -3), (3, -1), (1, 1)], [S3/S2, -1.0/S10, 0.2])
    result[15] = Y(3, [(3, 0), (3, 2), (1, -1)], [-3.0/S15, 1.0, 0.2])
    result[16] = Y(3, [(3, 1), (1, 0)], [4.0/S10, 0.2])
    result[17] = Y(3, [(3, -1), (1, 1)], [4.0/S10, 0.2])
    result[18] = Y(3, [(3, 0), (3, 2), (1, -1)], [-3.0/S15, -1.0, 0.2])
    result[19] = Y(3, [(3, -2)], [1.0])
    result[20] = Y(4, [(4, 0), (4, 2), (4, 4), (2, 0), (2, 2), (0, 0)], [3.0/S35, -2.0/S7, 1.0, -2.0*S3/7.0, 6.0/7.0, 0.2])
    result[21] = Y(4, [(4, 0), (4, 2), (4, 4), (2, 0), (2, 2), (0, 0)], [3.0/S35, 2.0/S7, 1.0, -2.0*S3/7.0, -6.0/7.0, 0.2])
    result[22] = Y(4, [(4, 0), (2, 0), (0, 0)], [8.0/S35, 4.0*S3/7.0, 0.2])
    result[23] = Y(4, [(4, -4), (4, -2), (2, -2)], [1.0, -1.0/S7, 3.0/7.0])
    result[24] = Y(4, [(4, 1), (4, 3), (2, 1)], [-3.0/S14, 1.0/S2, 3.0/7.0])
    result[25] = Y(4, [(4, -4), (4, -2), (2, -2)], [-1.0, -1.0/S7, 3.0/7.0])
    result[26] = Y(4, [(4, -3), (4, -1), (2, -1)], [-1.0/S2, -3.0/S14, 3.0/7.0])
    result[27] = Y(4, [(4, 1), (2, 1)], [4.0/S14, 3.0/7.0])
    result[28] = Y(4, [(4, -1), (2, -1)], [4.0/S14, 3.0/7.0])
    result[29] = Y(4, [(4, 0), (4, 4), (2, 0), (0, 0)], [1.0/S35, -1.0, -2.0*S3/21.0, 1.0/15.0])
    result[30] = Y(4, [(4, 0), (4, 2), (2, 0), (2, 2), (0, 0)], [-4.0/S35, 2.0/S7, S3/21.0, 1.0/7.0, 1.0/15.0])
    result[31] = Y(4, [(4, 0), (4, 2), (2, 0), (2, 2), (0, 0)], [-4.0/S35, -2.0/S7, S3/21.0, -1.0/7.0, 1.0/15.0])
    result[32] = Y(4, [(4, -3), (4, -1), (2, -1)], [1.0/S2, -1.0/S14, 1.0/7.0])
    result[33] = Y(4, [(4, 1), (4, 3), (2, 1)], [-1.0/S14, -1.0/S2, 1.0/7.0])
    result[34] = Y(4, [(4, -2), (2, -2)], [2.0/S7, 1.0/7.0])

    return result

def molden_cart_norms2_over_4pi() -> np.ndarray:
    """Normalization factors for MOLDEN cartesian functions over 4π."""
    result = np.zeros(len(MOLDEN_CART_POWERS))
    for i, p in enumerate(MOLDEN_CART_POWERS):
        res = primitive_int_1d(2 * p[0], 2.0)
        res *= primitive_int_1d(2 * p[1], 2.0)
        res *= primitive_int_1d(2 * p[2], 2.0)
        L = p[0] + p[1] + p[2]
        res /= primitive_int_1d_sphr(2 + 2 * L, 2.0)
        res /= 4 * np.pi
        result[i] = res
    return result