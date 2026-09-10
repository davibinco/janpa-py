"""Precomputed angular integrals for YLM products and nested-function derivatives.

Replaces AngularInts.java and NestedFunctionsDerivs.java.
"""

from __future__ import annotations

import numpy as np

from .polynom3d import ensure_cnk_enough

__all__ = [
    "angular_eval",
    "angular_evaluate",
    "ensure_ank_enough",
    "x_sqrt_n_inv_derivs",
    "AngularInts",
    "NestedFunctionsDerivs",
    "eval",
    "evaluate",
]

_SQRT2 = np.sqrt(2.0)
_SQRT3 = np.sqrt(3.0)
_SQRT5 = np.sqrt(5.0)
_SQRT7 = np.sqrt(7.0)
_SQRT10 = np.sqrt(10.0)
_SQRT_PI = np.sqrt(np.pi)
_SQRT_PI_3 = np.pi ** (3.0 / 2.0)

_SP3 = _SQRT_PI_3
_S2 = _SQRT2
_S3 = _SQRT3
_S5 = _SQRT5
_S7 = _SQRT7
_S10 = _SQRT10


def x_sqrt_n_inv_derivs(max_n: int, max_j: int) -> np.ndarray:
    """Derivatives of x^{-n/2}.

    result[n, j] is the coefficient of x^{-n/2-j} in d^j/dx^j x^{-n/2}.
    """
    max_n = int(max_n)
    max_j = int(max_j)

    result = np.zeros((max_n + 1, max_j + 1), dtype=float)
    for n in range(max_n + 1):
        result[n, 0] = 1.0
        for j in range(1, max_j + 1):
            result[n, j] = -result[n, j - 1] * (n / 2.0 + j - 1.0)
    return result


_ANK: list[list[float]] = []
_ZETA_EFF_COEFS = np.zeros((0, 0), dtype=float)


def ensure_ank_enough(max_n: int) -> list[list[float]]:
    """Ensure the module-level falling-factorial coefficient cache covers 0..max_n.

    Ank[n][mu] = n! / (n-mu)! for 0 <= mu <= n.
    """
    global _ANK

    max_n = int(max_n)
    if max_n < 0:
        return _ANK

    if len(_ANK) <= max_n:
        for n in range(len(_ANK), max_n + 1):
            if n == 0:
                _ANK.append([1.0])
            else:
                row = [0.0] * (n + 1)
                row[n] = _ANK[n - 1][n - 1] * float(n)
                for mu in range(n - 1, -1, -1):
                    row[mu] = row[mu + 1] / float(n - mu)
                _ANK.append(row)

    return _ANK


def _ensure_zeta_eff_coefs(max_n: int, max_j: int) -> np.ndarray:
    """Ensure cached coefficients for derivatives of zeta_eff^{-n/2}."""
    global _ZETA_EFF_COEFS

    max_n = int(max_n)
    max_j = int(max_j)

    current_n = _ZETA_EFF_COEFS.shape[0] - 1
    current_j = _ZETA_EFF_COEFS.shape[1] - 1

    if current_n < max_n or current_j < max_j:
        new_n = max(current_n, max_n)
        new_j = max(current_j, max_j)
        _ZETA_EFF_COEFS = x_sqrt_n_inv_derivs(new_n, new_j)

    return _ZETA_EFF_COEFS


# ----------------------------------------------------------------------
# Hardcoded AngularInts coefficient tables
# ----------------------------------------------------------------------
_ANGULAR_COEF_TABLE = {
    # l1 = 0
    (0, 0, 0): ((0, 0, _SP3),),
    (0, 1, 0): ((0, 1, -_SP3 * _S3),),
    (0, 2, 0): ((0, 2, _SP3 * _S5),),
    (0, 3, 0): ((0, 3, -_SP3 * _S7),),
    (0, 4, 0): ((0, 4, 3.0 * _SP3),),
    # l1 = 1, l2 = 1
    (1, 1, -1): ((0, 0, 1.5 * _SP3),),
    (1, 1, 0): ((1, 1, 3.0 * _SP3), (0, 0, 1.5 * _SP3)),
    (1, 1, 1): ((0, 0, 1.5 * _SP3),),
    # l1 = 1, l2 = 2
    (1, 2, -1): ((0, 1, -1.5 * _SP3 * _S5),),
    (1, 2, 0): ((1, 2, -_SP3 * _S3 * _S5), (0, 1, -_SP3 * _S3 * _S5)),
    (1, 2, 1): ((0, 1, -1.5 * _SP3 * _S5),),
    # l1 = 1, l2 = 3
    (1, 3, -1): ((0, 2, 1.5 * _SP3 * _S2 * _S7),),
    (1, 3, 0): ((1, 3, _SP3 * _S3 * _S7), (0, 2, 1.5 * _SP3 * _S3 * _S7)),
    (1, 3, 1): ((0, 2, 1.5 * _SP3 * _S2 * _S7),),
    # l1 = 1, l2 = 4
    (1, 4, -1): ((0, 3, -1.5 * _SP3 * _S3 * _S10),),
    (1, 4, 0): ((1, 4, -3.0 * _SP3 * _S3), (0, 3, -6.0 * _SP3 * _S3)),
    (1, 4, 1): ((0, 3, -1.5 * _SP3 * _S3 * _S10),),
    # l1 = 2, l2 = 2
    (2, 2, -2): ((0, 0, (15.0 / 4.0) * _SP3),),
    (2, 2, -1): ((1, 1, (15.0 / 2.0) * _SP3), (0, 0, (15.0 / 4.0) * _SP3)),
    (2, 2, 0): (
        (2, 2, 5.0 * _SP3),
        (1, 1, 10.0 * _SP3),
        (0, 0, (15.0 / 4.0) * _SP3),
    ),
    (2, 2, 1): ((1, 1, (15.0 / 2.0) * _SP3), (0, 0, (15.0 / 4.0) * _SP3)),
    (2, 2, 2): ((0, 0, (15.0 / 4.0) * _SP3),),
    # l1 = 2, l2 = 3
    (2, 3, -2): ((0, 1, -(15.0 / 4.0) * _SP3 * _S7),),
    (2, 3, -1): (
        (1, 2, -1.5 * _SP3 * _S5 * _S2 * _S7),
        (0, 1, -1.5 * _SP3 * _S5 * _S2 * _S7),
    ),
    (2, 3, 0): (
        (2, 3, -_SP3 * _S5 * _S7),
        (1, 2, -3.0 * _SP3 * _S5 * _S7),
        (0, 1, -(9.0 / 4.0) * _SP3 * _S5 * _S7),
    ),
    (2, 3, 1): (
        (1, 2, -1.5 * _SP3 * _S5 * _S2 * _S7),
        (0, 1, -1.5 * _SP3 * _S5 * _S2 * _S7),
    ),
    (2, 3, 2): ((0, 1, -(15.0 / 4.0) * _SP3 * _S7),),
    # l1 = 2, l2 = 4
    (2, 4, -2): ((0, 2, (45.0 / 4.0) * _SP3 * _S3),),
    (2, 4, -1): (
        (1, 3, (15.0 / 2.0) * _SP3 * _S3 * _S2),
        (0, 2, (45.0 / 4.0) * _SP3 * _S3 * _S2),
    ),
    (2, 4, 0): (
        (2, 4, 3.0 * _SP3 * _S5),
        (1, 3, 12.0 * _SP3 * _S5),
        (0, 2, (27.0 / 2.0) * _SP3 * _S5),
    ),
    (2, 4, 1): (
        (1, 3, (15.0 / 2.0) * _SP3 * _S3 * _S2),
        (0, 2, (45.0 / 4.0) * _SP3 * _S3 * _S2),
    ),
    (2, 4, 2): ((0, 2, (45.0 / 4.0) * _SP3 * _S3),),
    # l1 = 3, l2 = 3
    (3, 3, -3): ((0, 0, (105.0 / 8.0) * _SP3),),
    (3, 3, -2): ((1, 1, (105.0 / 4.0) * _SP3), (0, 0, (105.0 / 8.0) * _SP3)),
    (3, 3, -1): (
        (2, 2, 21.0 * _SP3),
        (1, 1, 42.0 * _SP3),
        (0, 0, (105.0 / 8.0) * _SP3),
    ),
    (3, 3, 0): (
        (3, 3, 7.0 * _SP3),
        (2, 2, (63.0 / 2.0) * _SP3),
        (1, 1, (189.0 / 4.0) * _SP3),
        (0, 0, (105.0 / 8.0) * _SP3),
    ),
    (3, 3, 1): (
        (2, 2, 21.0 * _SP3),
        (1, 1, 42.0 * _SP3),
        (0, 0, (105.0 / 8.0) * _SP3),
    ),
    (3, 3, 2): ((1, 1, (105.0 / 4.0) * _SP3), (0, 0, (105.0 / 8.0) * _SP3)),
    (3, 3, 3): ((0, 0, (105.0 / 8.0) * _SP3),),
    # l1 = 3, l2 = 4
    (3, 4, -3): ((0, 1, -(315.0 / 8.0) * _SP3),),
    (3, 4, -2): (
        (1, 2, -(45.0 / 4.0) * _SP3 * _S3 * _S7),
        (0, 1, -(45.0 / 4.0) * _SP3 * _S3 * _S7),
    ),
    (3, 4, -1): (
        (2, 3, -3.0 * _SP3 * _S3 * _S7 * _S5),
        (1, 2, -9.0 * _SP3 * _S3 * _S7 * _S5),
        (0, 1, -(45.0 / 8.0) * _SP3 * _S3 * _S7 * _S5),
    ),
    (3, 4, 0): (
        (3, 4, -3.0 * _SP3 * _S7),
        (2, 3, -18.0 * _SP3 * _S7),
        (1, 2, -(81.0 / 2.0) * _SP3 * _S7),
        (0, 1, -(45.0 / 2.0) * _SP3 * _S7),
    ),
    (3, 4, 1): (
        (2, 3, -3.0 * _SP3 * _S3 * _S7 * _S5),
        (1, 2, -9.0 * _SP3 * _S3 * _S7 * _S5),
        (0, 1, -(45.0 / 8.0) * _SP3 * _S3 * _S7 * _S5),
    ),
    (3, 4, 2): (
        (1, 2, -(45.0 / 4.0) * _SP3 * _S3 * _S7),
        (0, 1, -(45.0 / 4.0) * _SP3 * _S3 * _S7),
    ),
    (3, 4, 3): ((0, 1, -(315.0 / 8.0) * _SP3),),
    # l1 = 4, l2 = 4
    (4, 4, -4): ((0, 0, (945.0 / 16.0) * _SP3),),
    (4, 4, -3): ((1, 1, (945.0 / 8.0) * _SP3), (0, 0, (945.0 / 16.0) * _SP3)),
    (4, 4, -2): (
        (2, 2, (405.0 / 4.0) * _SP3),
        (1, 1, (405.0 / 2.0) * _SP3),
        (0, 0, (945.0 / 16.0) * _SP3),
    ),
    (4, 4, -1): (
        (3, 3, 45.0 * _SP3),
        (2, 2, (405.0 / 2.0) * _SP3),
        (1, 1, (2025.0 / 8.0) * _SP3),
        (0, 0, (945.0 / 16.0) * _SP3),
    ),
    (4, 4, 0): (
        (4, 4, 9.0 * _SP3),
        (3, 3, 72.0 * _SP3),
        (2, 2, 243.0 * _SP3),
        (1, 1, 270.0 * _SP3),
        (0, 0, (945.0 / 16.0) * _SP3),
    ),
    (4, 4, 1): (
        (3, 3, 45.0 * _SP3),
        (2, 2, (405.0 / 2.0) * _SP3),
        (1, 1, (2025.0 / 8.0) * _SP3),
        (0, 0, (945.0 / 16.0) * _SP3),
    ),
    (4, 4, 2): (
        (2, 2, (405.0 / 4.0) * _SP3),
        (1, 1, (405.0 / 2.0) * _SP3),
        (0, 0, (945.0 / 16.0) * _SP3),
    ),
    (4, 4, 3): ((1, 1, (945.0 / 8.0) * _SP3), (0, 0, (945.0 / 16.0) * _SP3)),
    (4, 4, 4): ((0, 0, (945.0 / 16.0) * _SP3),),
}


def angular_eval(l1: int, l2: int, m: int) -> np.ndarray:
    """Return the coefs matrix for angular integral of Y(l1,m)*Y(l2,m).

    Shape: (l1+1, l2+1). Hardcoded for L <= 4.

    If l1 > l2, the swapped table entry is returned transposed as a
    convenience. The Java code expects callers to canonicalize l1 <= l2.
    """
    l1 = int(l1)
    l2 = int(l2)
    m = int(m)

    if l1 > l2:
        return angular_eval(l2, l1, m).T.copy()

    key = (l1, l2, m)
    terms = _ANGULAR_COEF_TABLE.get(key)
    if terms is None:
        raise ValueError(
            f"Unable to integrate YLMs for l1={l1}, l2={l2}, m={m} ! "
        )

    coefs = np.zeros((l1 + 1, l2 + 1), dtype=float)
    for i, j, value in terms:
        coefs[i, j] = value
    return coefs


def angular_evaluate(
    l1: int,
    l2: int,
    m: int,
    zeta1_derivs: int,
    zeta2_derivs: int,
    zeta1: float,
    zeta2: float,
    dz: float,
) -> np.ndarray:
    """Full angular integral with zeta-derivative support.

    Returns shape (zeta1_derivs+1, zeta2_derivs+1), in the original
    (pre-swap) derivative ordering.
    """
    l1 = int(l1)
    l2 = int(l2)
    m = int(m)
    zeta1_derivs = int(zeta1_derivs)
    zeta2_derivs = int(zeta2_derivs)
    zeta1 = float(zeta1)
    zeta2 = float(zeta2)
    dz = float(dz)

    transposed = False

    # Canonicalize l1 <= l2, swapping all associated quantities.
    if l2 < l1:
        l1, l2 = l2, l1
        zeta1, zeta2 = zeta2, zeta1
        zeta1_derivs, zeta2_derivs = zeta2_derivs, zeta1_derivs
        dz = -dz
        transposed = True

    coefs = angular_eval(l1, l2, m)

    zeta_eff = zeta1 + zeta2
    inv_zeta_eff = 1.0 / zeta_eff
    sqrt_inv_zeta_eff_dz = dz * np.sqrt(inv_zeta_eff)

    max_pwr_zetas = l1 + l2

    zeta1_pwrs = np.empty(max_pwr_zetas + 1, dtype=float)
    zeta2_pwrs = np.empty(max_pwr_zetas + 1, dtype=float)
    zeta1_pwrs[0] = 1.0
    zeta2_pwrs[0] = 1.0
    for i in range(1, max_pwr_zetas + 1):
        zeta1_pwrs[i] = zeta1_pwrs[i - 1] * zeta1
        zeta2_pwrs[i] = zeta2_pwrs[i - 1] * zeta2

    zeta_eff_inv_pwrs = np.empty(zeta1_derivs + zeta2_derivs + 1, dtype=float)
    zeta_eff_inv_pwrs[0] = 1.0
    for i in range(1, zeta_eff_inv_pwrs.shape[0]):
        zeta_eff_inv_pwrs[i] = zeta_eff_inv_pwrs[i - 1] * inv_zeta_eff

    zeta_eff_md_inv_pwrs05 = np.empty(l1 + l2 + 1, dtype=float)
    zeta_eff_md_inv_pwrs05[0] = inv_zeta_eff ** ((l1 + l2 + 3.0) / 2.0)
    for i in range(1, zeta_eff_md_inv_pwrs05.shape[0]):
        zeta_eff_md_inv_pwrs05[i] = (
            zeta_eff_md_inv_pwrs05[i - 1] * sqrt_inv_zeta_eff_dz
        )

    cnk = ensure_cnk_enough(max(zeta1_derivs, zeta2_derivs))
    ank = ensure_ank_enough(max(l1, l2))

    max_zeta_eff_coef_n = 2 * (l1 + l2) + 3
    max_zeta_eff_coef_j = zeta1_derivs + zeta2_derivs
    zeta_eff_coefs = _ensure_zeta_eff_coefs(
        max_zeta_eff_coef_n, max_zeta_eff_coef_j
    )

    result = np.zeros((zeta1_derivs + 1, zeta2_derivs + 1), dtype=float)

    for i in range(l1 + 1):
        for j in range(l2 + 1):
            coef_ij = coefs[i, j]
            if coef_ij == 0.0:
                continue

            sign = 1.0 if (i % 2 == 0) else -1.0

            result[0, 0] += (
                sign
                * coef_ij
                * zeta2_pwrs[i]
                * zeta1_pwrs[j]
                * zeta_eff_md_inv_pwrs05[i + j]
            )

            for n in range(zeta1_derivs + 1):
                for k in range(zeta2_derivs + 1):
                    if n + k == 0:
                        continue

                    tmp = 0.0
                    max_mu = min(n, j)
                    max_nu = min(k, i)

                    for mu in range(max_mu + 1):
                        for nu in range(max_nu + 1):
                            p = i + j + l1 + l2 + 3
                            r = n - mu + k - nu

                            tmp += (
                                cnk[n][mu]
                                * cnk[k][nu]
                                * ank[j][mu]
                                * zeta1_pwrs[j - mu]
                                * ank[i][nu]
                                * zeta2_pwrs[i - nu]
                                * zeta_eff_coefs[p, r]
                                * zeta_eff_inv_pwrs[n + k - mu - nu]
                            )

                    tmp *= zeta_eff_md_inv_pwrs05[i + j] * coef_ij
                    result[n, k] += sign * tmp

    if transposed:
        return result.T.copy()

    return result


class AngularInts:
    """Java-compatible static facade around angular_eval/angular_evaluate."""

    @staticmethod
    def eval(l1: int, l2: int, m: int) -> np.ndarray:
        """Java AngularInts.eval."""
        return angular_eval(l1, l2, m)

    @staticmethod
    def evaluate(
        l1: int,
        l2: int,
        m: int,
        zeta1_derivs: int,
        zeta2_derivs: int,
        zeta1: float,
        zeta2: float,
        dz: float,
    ) -> np.ndarray:
        """Java AngularInts.evaluate."""
        return angular_evaluate(
            l1,
            l2,
            m,
            zeta1_derivs,
            zeta2_derivs,
            zeta1,
            zeta2,
            dz,
        )


class NestedFunctionsDerivs:
    """Nested-function derivative utilities.

    Replaces NestedFunctionsDerivs.java. The binomial coefficients are
    provided by janpa.gto.polynom3d.ensure_cnk_enough.
    """

    max_deriv = 3
    Cnk_max = 5

    def __init__(self) -> None:
        self._zeta_diff_ci: list[float] = []
        ensure_cnk_enough(self.Cnk_max)
        self._ensure_zeta_diff_ci_enough(self.max_deriv)

    # ------------------------------------------------------------------
    # Cache maintenance
    # ------------------------------------------------------------------
    def _ensure_zeta_diff_ci_enough(self, max_n: int) -> None:
        """Ensure zetaDiffCi coefficients up to order max_n."""
        max_n = int(max_n)
        if len(self._zeta_diff_ci) - 1 >= max_n:
            return

        old_len = len(self._zeta_diff_ci)
        self._zeta_diff_ci.extend([0.0] * (max_n + 1 - old_len))

        for i in range(old_len, max_n + 1):
            if i == 0:
                self._zeta_diff_ci[i] = -1.0
            else:
                self._zeta_diff_ci[i] = (
                    -self._zeta_diff_ci[i - 1] / 2.0 * (2.0 * i - 3.0)
                )

    @staticmethod
    def ensure_Ank_enough(max_n: int) -> list[list[float]]:
        """Java-compatible static ensure_Ank_enough."""
        return ensure_ank_enough(max_n)

    # ------------------------------------------------------------------
    # Derivative formulas
    # ------------------------------------------------------------------
    def product_derivs(
        self,
        fd: np.ndarray,
        gd: np.ndarray,
        k: float,
        max_i: int,
        max_j: int,
    ) -> np.ndarray:
        """Derivatives of k * f * g from derivative tables f and g."""
        fd = np.asarray(fd, dtype=float)
        gd = np.asarray(gd, dtype=float)
        max_i = int(max_i)
        max_j = int(max_j)
        k = float(k)

        cnk = ensure_cnk_enough(max(max_i, max_j))
        result = np.zeros((max_i + 1, max_j + 1), dtype=float)

        for i in range(max_i + 1):
            for j in range(max_j + 1):
                total = 0.0
                for mu in range(i + 1):
                    for nu in range(j + 1):
                        total += (
                            cnk[i][mu]
                            * cnk[j][nu]
                            * fd[mu, nu]
                            * gd[i - mu, j - nu]
                        )
                result[i, j] = total * k

        return result

    def ab_derivs(
        self,
        zeta1: float,
        zeta2: float,
        zab: float,
        max_i: int,
        max_j: int,
        transpose_result: bool,
    ) -> np.ndarray:
        """Derivatives of the A/B center coefficient factor."""
        zeta1 = float(zeta1)
        zeta2 = float(zeta2)
        zab = float(zab)
        max_i = int(max_i)
        max_j = int(max_j)
        transpose_result = bool(transpose_result)

        zeta_eff = zeta1 + zeta2
        self._ensure_zeta_diff_ci_enough(max_i + max_j)

        zeta_eff_pwr = np.empty(max_i + max_j + 1, dtype=float)
        zeta_eff_pwr[0] = 1.0 / np.sqrt(zeta_eff)
        for i in range(1, zeta_eff_pwr.shape[0]):
            zeta_eff_pwr[i] = zeta_eff_pwr[i - 1] / zeta_eff

        if transpose_result:
            result = np.zeros((max_j + 1, max_i + 1), dtype=float)
            for i in range(max_i + 1):
                for j in range(max_j + 1):
                    result[j, i] = (
                        zeta_eff_pwr[i + j]
                        * ((2.0 * j - 1.0) * zeta1 - 2.0 * i * zeta2)
                        * zab
                        * self._zeta_diff_ci[i + j]
                    )
        else:
            result = np.zeros((max_i + 1, max_j + 1), dtype=float)
            for i in range(max_i + 1):
                for j in range(max_j + 1):
                    result[i, j] = (
                        zeta_eff_pwr[i + j]
                        * ((2.0 * j - 1.0) * zeta1 - 2.0 * i * zeta2)
                        * zab
                        * self._zeta_diff_ci[i + j]
                    )

        return result

    def exponential_derivs(
        self,
        fd: np.ndarray,
        max_i: int,
        max_j: int,
    ) -> np.ndarray:
        """Derivatives of exp(f) from derivative table f."""
        fd = np.asarray(fd, dtype=float)
        max_i = int(max_i)
        max_j = int(max_j)

        result = np.zeros((max_i + 1, max_j + 1), dtype=float)
        result[0, 0] = np.exp(fd[0, 0])

        cnk = ensure_cnk_enough(max(max_i, max_j))

        for i in range(1, max_i + 1):
            total = 0.0
            for mu in range(i):
                total += cnk[i - 1][mu] * result[mu, 0] * fd[i - mu, 0]
            result[i, 0] = total

        for i in range(max_i + 1):
            for j in range(max_j):
                total = 0.0
                for mu in range(i + 1):
                    for nu in range(j + 1):
                        total += (
                            cnk[i][mu]
                            * cnk[j][nu]
                            * result[mu, nu]
                            * fd[i - mu, j - nu + 1]
                        )
                result[i, j + 1] = total

        return result

    def prefactor_exponent_derivs(
        self,
        max_i: int,
        max_j: int,
        dz: float,
        zeta1: float,
        zeta2: float,
    ) -> np.ndarray:
        """Derivatives of the Gaussian prefactor exponent."""
        max_i = int(max_i)
        max_j = int(max_j)
        dz = float(dz)
        zeta1 = float(zeta1)
        zeta2 = float(zeta2)

        aux = np.empty(max_i + 1 + max_j + 1, dtype=float)
        minus_zeta_eff_inv = -1.0 / (zeta1 + zeta2)

        aux[0] = 1.0
        aux[1] = minus_zeta_eff_inv
        for i in range(2, aux.shape[0]):
            aux[i] = aux[i - 1] * minus_zeta_eff_inv * float(i - 1)

        result = np.zeros((max_i + 1, max_j + 1), dtype=float)
        dz2 = dz * dz
        zeta12 = zeta1 * zeta1
        zeta22 = zeta2 * zeta2

        result[0, 0] = zeta1 * zeta2 * minus_zeta_eff_inv * dz2

        for i in range(max_i + 1):
            for j in range(max_j + 1):
                if i + j == 0:
                    continue

                tmp = (float(i) * zeta22 + float(j) * zeta12) * minus_zeta_eff_inv
                if i > 0 and j > 0:
                    tmp += (
                        float(i)
                        * float(j)
                        * (zeta1 + zeta2)
                        / float(i + j - 1)
                    )

                result[i, j] = aux[i + j] * tmp * (-dz2)

        return result

    @staticmethod
    def xSqrtNInvDerivs(max_n: int, max_j: int) -> np.ndarray:
        """Java-compatible xSqrtNInvDerivs."""
        return x_sqrt_n_inv_derivs(max_n, max_j)

    # Java-compatible aliases
    productDerivs = product_derivs
    abDerivs = ab_derivs
    exponentialDerivs = exponential_derivs
    prefactorExponentDerivs = prefactor_exponent_derivs


# Java/module-compatible aliases
eval = angular_eval
evaluate = angular_evaluate