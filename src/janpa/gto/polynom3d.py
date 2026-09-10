"""3D polynomial representation for solid harmonics.

This module replaces the Java classes Polynom3D and YLM_Series.

The Java "temporary array" pattern (`coefs_temp`, `powers_temp`,
`temp_num_terms`) is replaced here by Python dictionaries / intermediate
lists inside helper methods.

After construction, `Polynom3D.coefs` and `Polynom3D.powers` are always
NumPy arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.special import comb


__all__ = [
    "Polynom3D",
    "YLMSeries",
    "ensure_cnk_enough",
    "ensure_Cnk_enough",
    "c_nk",
    "C_nk",
]


# ---------------------------------------------------------------------------
# Binomial coefficient cache
# ---------------------------------------------------------------------------

_QUICK_CNK: list[list[float]] = []


def ensure_cnk_enough(max_n: int) -> list[list[float]]:
    """Ensure the module-level binomial coefficient cache covers 0..max_n.

    Returns the cache as a list of rows:
        cache[n][k] == C(n, k), for 0 <= k <= n.
    """
    max_n = int(max_n)
    if max_n < 0:
        return _QUICK_CNK

    if len(_QUICK_CNK) <= max_n:
        for n in range(len(_QUICK_CNK), max_n + 1):
            _QUICK_CNK.append(
                [float(comb(n, k, exact=True)) for k in range(n + 1)]
            )

    return _QUICK_CNK


# Java-style alias.
ensure_Cnk_enough = ensure_cnk_enough


def _binom(n: int, k: int) -> float:
    """Cached binomial coefficient used by Polynom3D.c_nk."""
    n = int(n)
    k = int(k)

    if k < 0 or k > n:
        return 0.0

    ensure_cnk_enough(n)
    return _QUICK_CNK[n][k]


def c_nk(n: int, k: int) -> float:
    """Module-level binomial coefficient C(n, k)."""
    return _binom(n, k)


# Java-style alias.
C_nk = c_nk


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _prepare_powers(powers: Any, n_coefs: int) -> np.ndarray:
    """Convert a powers-like object to an `(n_coefs, 3)` int array."""
    if powers is None:
        return np.zeros((n_coefs, 3), dtype=int)

    arr = np.asarray(powers, dtype=int)

    if arr.size == 0:
        return np.zeros((n_coefs, 3), dtype=int)

    if arr.ndim == 1:
        if arr.size % 3 != 0:
            raise ValueError("powers must contain triples (nx, ny, nz).")
        arr = arr.reshape((-1, 3))
    elif arr.ndim == 2 and arr.shape[1] == 3:
        pass
    else:
        if arr.size % 3 != 0:
            raise ValueError("powers must contain triples (nx, ny, nz).")
        arr = arr.reshape((-1, 3))

    if arr.shape[0] != n_coefs:
        if arr.size == n_coefs * 3:
            arr = arr.reshape((n_coefs, 3))
        else:
            raise ValueError("powers must have one row per coefficient.")

    return arr.copy()


# ---------------------------------------------------------------------------
# Polynom3D
# ---------------------------------------------------------------------------


@dataclass(eq=False)
class Polynom3D:
    """A 3D polynomial: sum of c_i * x^nx * y^ny * z^nz terms.

    Attributes
    ----------
    coefs:
        `(n_terms,)` float array.
    powers:
        `(n_terms, 3)` int array. Each row is `(nx, ny, nz)`.

    Constructor forms
    -----------------
    * `Polynom3D()`
        Zero polynomial.
    * `Polynom3D(num_terms)`
        Java-compatible allocation of `num_terms` zero terms.
    * `Polynom3D(coefs, powers)`
        Create from coefficient and power arrays.
    * `Polynom3D(source)`
        Copy constructor.
    """

    coefs: Any = None
    powers: Any = None

    def __post_init__(self) -> None:
        coefs_in = self.coefs
        powers_in = self.powers

        # Copy constructor.
        if isinstance(coefs_in, Polynom3D):
            source = coefs_in
            self.coefs = source.coefs.copy()
            self.powers = source.powers.copy()
            return

        # Default zero polynomial.
        if coefs_in is None:
            self.set_zero()
            return

        # Java-compatible constructor Polynom3D(int NumTerms).
        # An integer with omitted powers is treated as a term capacity.
        if (
            powers_in is None
            and isinstance(coefs_in, (int, np.integer))
            and not isinstance(coefs_in, bool)
        ):
            num_terms = int(coefs_in)
            if num_terms <= 0:
                self.set_zero()
            else:
                self.coefs = np.zeros(num_terms, dtype=float)
                self.powers = np.zeros((num_terms, 3), dtype=int)
            return

        coefs = np.asarray(coefs_in, dtype=float)
        if coefs.ndim == 0:
            coefs = coefs.reshape(1)
        else:
            coefs = coefs.ravel()

        if coefs.size == 0:
            self.set_zero()
            return

        self.coefs = coefs.copy()
        self.powers = _prepare_powers(powers_in, coefs.size)

    # ------------------------------------------------------------------
    # Basic Java-compatible helpers
    # ------------------------------------------------------------------

    @property
    def n_terms(self) -> int:
        """Number of stored terms."""
        return int(self.coefs.size)

    def set_zero(self) -> Polynom3D:
        """Set this polynomial to the canonical zero polynomial."""
        self.coefs = np.array([0.0], dtype=float)
        self.powers = np.array([[0, 0, 0]], dtype=int)
        return self

    def copy(self) -> Polynom3D:
        """Return a deep copy."""
        return self.__class__(self.coefs.copy(), self.powers.copy())

    def create_from_arrays(
        self,
        num_terms: int,
        coefs_: Any,
        powers_: Any,
    ) -> Polynom3D:
        """Java-compatible CreateFromArrays."""
        num_terms = int(num_terms)
        if num_terms <= 0:
            return self.set_zero()

        if coefs_ is None:
            coefs_arr = np.array([], dtype=float)
        else:
            coefs_arr = np.asarray(coefs_, dtype=float).ravel()

        if coefs_arr.size < num_terms:
            raise ValueError("coefs_ is too short for num_terms.")

        if powers_ is None:
            powers_arr = np.zeros((num_terms, 3), dtype=int)
        else:
            powers_arr = np.asarray(powers_, dtype=int)

            if powers_arr.size == 0:
                powers_arr = np.zeros((num_terms, 3), dtype=int)
            else:
                if powers_arr.ndim == 1:
                    if powers_arr.size % 3 != 0:
                        raise ValueError(
                            "powers_ must contain triples (nx, ny, nz)."
                        )
                    powers_arr = powers_arr.reshape((-1, 3))
                elif powers_arr.ndim == 2 and powers_arr.shape[1] == 3:
                    pass
                else:
                    if powers_arr.size % 3 != 0:
                        raise ValueError(
                            "powers_ must contain triples (nx, ny, nz)."
                        )
                    powers_arr = powers_arr.reshape((-1, 3))

                if powers_arr.shape[0] < num_terms:
                    raise ValueError("powers_ is too short for num_terms.")

                powers_arr = powers_arr[:num_terms]

        self.coefs = coefs_arr[:num_terms].copy()
        self.powers = powers_arr.copy()
        return self

    # ------------------------------------------------------------------
    # Internal term-dictionary conversion
    # ------------------------------------------------------------------

    def _assign_terms(
        self,
        terms: dict[tuple[int, int, int], float],
        epsilon: float = 0.0,
    ) -> Polynom3D:
        """Replace self by the terms in a dictionary.

        Terms with `abs(coef) <= epsilon` are discarded.
        If no terms remain, the canonical zero polynomial is stored.
        """
        coefs: list[float] = []
        powers: list[tuple[int, int, int]] = []

        for key, value in terms.items():
            if abs(value) > epsilon:
                coefs.append(float(value))
                powers.append((int(key[0]), int(key[1]), int(key[2])))

        if coefs:
            self.coefs = np.asarray(coefs, dtype=float)
            self.powers = np.asarray(powers, dtype=int)
        else:
            self.set_zero()

        return self

    @classmethod
    def _from_terms(
        cls,
        terms: dict[tuple[int, int, int], float],
        epsilon: float = 0.0,
    ) -> Polynom3D:
        """Create a new polynomial from a term dictionary."""
        return cls()._assign_terms(terms, epsilon)

    # Compatibility with the older skeleton name.
    _from_dict = _assign_terms

    # ------------------------------------------------------------------
    # Polynomial arithmetic
    # ------------------------------------------------------------------

    def multiply_by(self, other: Polynom3D) -> Polynom3D:
        """self *= other. Mutates self and returns self."""
        terms: dict[tuple[int, int, int], float] = {}

        for c1, p1 in zip(self.coefs, self.powers):
            if c1 == 0.0:
                continue

            x1 = int(p1[0])
            y1 = int(p1[1])
            z1 = int(p1[2])

            for c2, p2 in zip(other.coefs, other.powers):
                if c2 == 0.0:
                    continue

                key = (
                    x1 + int(p2[0]),
                    y1 + int(p2[1]),
                    z1 + int(p2[2]),
                )
                terms[key] = terms.get(key, 0.0) + c1 * c2

        return self._assign_terms(terms)

    def add_poly(self, other: Polynom3D, c: float = 1.0) -> Polynom3D:
        """self += c * other. Mutates self and returns self."""
        terms: dict[tuple[int, int, int], float] = {}

        for cf, p in zip(self.coefs, self.powers):
            if cf == 0.0:
                continue
            key = (int(p[0]), int(p[1]), int(p[2]))
            terms[key] = terms.get(key, 0.0) + cf

        c = float(c)
        if c != 0.0:
            for cf, p in zip(other.coefs, other.powers):
                value = c * cf
                if value == 0.0:
                    continue
                key = (int(p[0]), int(p[1]), int(p[2]))
                terms[key] = terms.get(key, 0.0) + value

        return self._assign_terms(terms)

    def laplacian(self) -> Polynom3D:
        """Return a new polynomial equal to ∇²(self)."""
        terms: dict[tuple[int, int, int], float] = {}

        for cf, p in zip(self.coefs, self.powers):
            if cf == 0.0:
                continue

            nx = int(p[0])
            ny = int(p[1])
            nz = int(p[2])

            if nx >= 2:
                key = (nx - 2, ny, nz)
                terms[key] = terms.get(key, 0.0) + cf * nx * (nx - 1)

            if ny >= 2:
                key = (nx, ny - 2, nz)
                terms[key] = terms.get(key, 0.0) + cf * ny * (ny - 1)

            if nz >= 2:
                key = (nx, ny, nz - 2)
                terms[key] = terms.get(key, 0.0) + cf * nz * (nz - 1)

        return self._from_terms(terms)

    def differentiate_mu(self, mu: int) -> Polynom3D:
        """Return ∂/∂x_mu as a new polynomial.

        mu = 0 -> x
        mu = 1 -> y
        mu = 2 -> z
        """
        mu = int(mu)
        if mu not in (0, 1, 2):
            raise ValueError("mu must be 0, 1, or 2.")

        terms: dict[tuple[int, int, int], float] = {}

        for cf, p in zip(self.coefs, self.powers):
            if cf == 0.0:
                continue

            power = int(p[mu])
            if power == 0:
                continue

            new_power = [int(p[0]), int(p[1]), int(p[2])]
            new_power[mu] = power - 1
            key = (new_power[0], new_power[1], new_power[2])

            terms[key] = terms.get(key, 0.0) + cf * power

        return self._from_terms(terms)

    def scale_coefs_by(self, factor: float) -> Polynom3D:
        """Multiply all coefficients by factor. Mutates self."""
        self.coefs = self.coefs * float(factor)
        return self

    # ------------------------------------------------------------------
    # Evaluation
    # ------------------------------------------------------------------

    def evaluate_at_point(self, r: Any) -> float:
        """Evaluate the polynomial at `r = (x, y, z)`.

        This keeps the Java optimization: powers of x, y, and z are
        precomputed once up to the maximum needed exponent.
        """
        r_arr = np.asarray(r, dtype=float).ravel()
        if r_arr.size != 3:
            raise ValueError("r must be a length-3 point: (x, y, z).")

        if self.coefs.size == 0 or self.powers.size == 0:
            return 0.0

        if np.any(self.powers < 0):
            raise ValueError("Polynomial powers must be non-negative.")

        max_powers = np.max(self.powers, axis=0)

        xyz_powers: list[np.ndarray] = []
        for mu in range(3):
            max_p = max(0, int(max_powers[mu]))
            vals = np.empty(max_p + 1, dtype=float)
            vals[0] = 1.0
            for pwr in range(1, max_p + 1):
                vals[pwr] = vals[pwr - 1] * r_arr[mu]
            xyz_powers.append(vals)

        values = (
            self.coefs
            * xyz_powers[0][self.powers[:, 0]]
            * xyz_powers[1][self.powers[:, 1]]
            * xyz_powers[2][self.powers[:, 2]]
        )

        return float(np.sum(values))

    # ------------------------------------------------------------------
    # Binomial coefficients
    # ------------------------------------------------------------------

    @staticmethod
    def c_nk(n: int, k: int) -> float:
        """Binomial coefficient C(n, k) as a float."""
        return _binom(n, k)

    @staticmethod
    def ensure_cnk_enough(max_n: int) -> list[list[float]]:
        """Ensure the global binomial cache covers up to max_n."""
        return ensure_cnk_enough(max_n)

    # ------------------------------------------------------------------
    # Python operator conveniences
    # ------------------------------------------------------------------

    def __mul__(self, other: Any) -> Polynom3D:
        if isinstance(other, Polynom3D):
            return self.copy().multiply_by(other)
        if np.isscalar(other):
            return self.copy().scale_coefs_by(float(other))
        return NotImplemented

    def __rmul__(self, other: Any) -> Polynom3D:
        return self.__mul__(other)

    def __imul__(self, other: Any) -> Polynom3D:
        if isinstance(other, Polynom3D):
            return self.multiply_by(other)
        if np.isscalar(other):
            return self.scale_coefs_by(float(other))
        return NotImplemented

    def __add__(self, other: Any) -> Polynom3D:
        if isinstance(other, Polynom3D):
            return self.copy().add_poly(other, 1.0)
        return NotImplemented

    def __iadd__(self, other: Any) -> Polynom3D:
        if isinstance(other, Polynom3D):
            return self.add_poly(other, 1.0)
        return NotImplemented

    # ------------------------------------------------------------------
    # Object protocol
    # ------------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Polynom3D):
            return NotImplemented
        return (
            np.array_equal(self.coefs, other.coefs)
            and np.array_equal(self.powers, other.powers)
        )

    __hash__ = None

    def __repr__(self) -> str:
        parts = [
            f"{c:.4f}·x^{p[0]}·y^{p[1]}·z^{p[2]}"
            for c, p in zip(self.coefs, self.powers)
        ]
        return " + ".join(parts) if parts else "0"

    __str__ = __repr__

    # ------------------------------------------------------------------
    # Java-compatible method-name aliases
    # ------------------------------------------------------------------

    MultiplyBy = multiply_by
    AddPoly = add_poly
    Laplacian = laplacian
    Differentiate_mu = differentiate_mu
    EvaluateAtPoint = evaluate_at_point
    ScaleCoefsBy = scale_coefs_by
    SetZero = set_zero
    CreateFromArrays = create_from_arrays
    C_nk = c_nk
    ensure_Cnk_enough = ensure_cnk_enough


# ---------------------------------------------------------------------------
# YLMSeries
# ---------------------------------------------------------------------------


@dataclass(init=False, eq=False)
class YLMSeries:
    """A linear combination of spherical harmonics Y_L^M.

    Java constructors represented:

    * `YLMSeries(LMain, LMs, Coefs)`
        `LMs` is an array-like of `(L, M)` rows.
    * `YLMSeries(source)`
        Copy constructor.
    * `YLMSeries(NTermsMax)`
        Allocate zero-filled arrays of length `NTermsMax`.

    The normal field-style constructor is also supported:
        `YLMSeries(primary_l, ls, ms, coefs, indexes=None)`.
    """

    primary_l: int
    ls: np.ndarray
    ms: np.ndarray
    coefs: np.ndarray
    indexes: np.ndarray

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        # Empty constructor.
        if not args and not kwargs:
            self._set_fields(
                0,
                np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=float),
                np.array([], dtype=int),
            )
            return

        # Copy constructor or capacity constructor.
        if len(args) == 1 and not kwargs:
            obj = args[0]

            if isinstance(obj, YLMSeries):
                self._set_fields(
                    obj.primary_l,
                    obj.ls,
                    obj.ms,
                    obj.coefs,
                    obj.indexes,
                )
                return

            if isinstance(obj, (int, np.integer)):
                n_terms_max = int(obj)
                if n_terms_max < 0:
                    raise ValueError("NTermsMax must be non-negative.")
                self._set_fields(
                    0,
                    np.zeros(n_terms_max, dtype=int),
                    np.zeros(n_terms_max, dtype=int),
                    np.zeros(n_terms_max, dtype=float),
                    np.zeros(n_terms_max, dtype=int),
                )
                return

        # Java constructor: YLM_Series(int LMain, int[][] LMs, double[] Coefs)
        if len(args) == 3 and not kwargs:
            primary_l, lms, coefs = args
            self._set_from_lms(primary_l, lms, coefs)
            return

        field_names = ("primary_l", "ls", "ms", "coefs", "indexes")
        if len(args) > len(field_names):
            raise TypeError(
                "too many positional arguments for YLMSeries "
                "(expected at most 5)"
            )

        params: dict[str, Any] = dict(zip(field_names, args))
        params.update(kwargs)

        # Allow Java-style primary_L keyword.
        if "primary_L" in params:
            params.setdefault("primary_l", params.pop("primary_L"))

        allowed = set(field_names) | {"lms", "primary_L"}
        unknown = set(params) - allowed
        if unknown:
            raise TypeError(
                f"unexpected keyword argument(s) for YLMSeries: "
                f"{', '.join(sorted(unknown))}"
            )

        # Optional LMs keyword constructor.
        if "lms" in params:
            primary_l = params.get("primary_l", 0)
            lms = params.get("lms")
            coefs = params.get("coefs", np.array([], dtype=float))
            self._set_from_lms(primary_l, lms, coefs)
            return

        self._set_fields(
            params.get("primary_l", 0),
            params.get("ls", np.array([], dtype=int)),
            params.get("ms", np.array([], dtype=int)),
            params.get("coefs", np.array([], dtype=float)),
            params.get("indexes", None),
        )

    # ------------------------------------------------------------------
    # Internal constructors / normalizers
    # ------------------------------------------------------------------

    def _set_from_lms(self, primary_l: Any, lms: Any, coefs: Any) -> None:
        if lms is None:
            lms_arr = np.asarray([], dtype=int)
        else:
            lms_arr = np.asarray(lms, dtype=int)

        if lms_arr.size == 0:
            lms_arr = lms_arr.reshape((0, 2))
        elif lms_arr.ndim == 1:
            if lms_arr.size % 2 != 0:
                raise ValueError("LMs must contain pairs (L, M).")
            lms_arr = lms_arr.reshape((-1, 2))
        elif lms_arr.ndim == 2 and lms_arr.shape[1] == 2:
            pass
        else:
            if lms_arr.size % 2 != 0:
                raise ValueError("LMs must contain pairs (L, M).")
            lms_arr = lms_arr.reshape((-1, 2))

        self._set_fields(
            primary_l,
            lms_arr[:, 0],
            lms_arr[:, 1],
            coefs,
            np.zeros(lms_arr.shape[0], dtype=int),
        )

    def _set_fields(
        self,
        primary_l: Any,
        ls: Any,
        ms: Any,
        coefs: Any,
        indexes: Any,
    ) -> None:
        self.primary_l = 0 if primary_l is None else int(primary_l)

        ls_arr = np.asarray([] if ls is None else ls, dtype=int).ravel().copy()
        ms_arr = np.asarray([] if ms is None else ms, dtype=int).ravel().copy()
        coefs_arr = (
            np.asarray([] if coefs is None else coefs, dtype=float)
            .ravel()
            .copy()
        )

        n = ls_arr.size

        if ms_arr.size != n:
            if ms_arr.size == 0:
                ms_arr = np.zeros(n, dtype=int)
            else:
                raise ValueError("ls and ms must have the same length.")

        if coefs_arr.size != n:
            if coefs_arr.size == 0:
                coefs_arr = np.zeros(n, dtype=float)
            else:
                raise ValueError("ls and coefs must have the same length.")

        if indexes is None:
            indexes_arr = np.zeros(n, dtype=int)
        else:
            indexes_arr = (
                np.asarray(indexes, dtype=int).ravel().copy()
            )
            if indexes_arr.size != n:
                if indexes_arr.size == 0:
                    indexes_arr = np.zeros(n, dtype=int)
                else:
                    raise ValueError(
                        "indexes must have the same length as ls."
                    )

        self.ls = ls_arr
        self.ms = ms_arr
        self.coefs = coefs_arr
        self.indexes = indexes_arr

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    @property
    def n_terms(self) -> int:
        """Number of stored series terms."""
        return int(self.coefs.size)

    def copy(self) -> YLMSeries:
        """Return a deep copy."""
        return YLMSeries(self)

    @classmethod
    def from_lms(cls, primary_l: int, lms: Any, coefs: Any) -> YLMSeries:
        """Create from `(L, M)` rows and coefficients."""
        return cls(primary_l, lms, coefs)

    @classmethod
    def zeros(cls, n_terms_max: int) -> YLMSeries:
        """Create a zero-filled series with capacity `n_terms_max`."""
        return cls(n_terms_max)

    # Java-style field alias.
    @property
    def primary_L(self) -> int:
        return self.primary_l

    @primary_L.setter
    def primary_L(self, value: int) -> None:
        self.primary_l = int(value)

    # ------------------------------------------------------------------
    # Object protocol
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self.n_terms

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, YLMSeries):
            return NotImplemented
        return (
            self.primary_l == other.primary_l
            and np.array_equal(self.ls, other.ls)
            and np.array_equal(self.ms, other.ms)
            and np.array_equal(self.indexes, other.indexes)
            and np.array_equal(self.coefs, other.coefs)
        )

    __hash__ = None