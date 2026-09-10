# janpa/gto/basis.py
"""Core data structures for Gaussian-Type Orbital basis sets and molecular data.

Replaces Java classes: AtomicCenter, BasisFunction, RadialPartOfBasisFunction, MO.
All linked-list patterns (_next pointers) are replaced by plain Python lists.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class AtomicCenter:
    """A single atomic center (nucleus) in the molecule."""
    name: str
    id: int                        # 1-based index, matching MOLDEN convention
    z: float                       # nuclear charge
    r0: np.ndarray = field(default_factory=lambda: np.zeros(3))  # (3,) coordinates

    def distance_to(self, r: np.ndarray) -> float:
        return float(np.linalg.norm(self.r0 - r))


@dataclass
class RadialPart:
    """Shared radial part (exponents + contraction coefficients) for one shell.

    Replaces RadialPartOfBasisFunction. Multiple BasisFunction objects with
    the same RadialPart share the same exponents/coefficients but differ in m.
    """
    exponents: np.ndarray          # (n_prim,)
    coefs: np.ndarray              # (n_prim,)
    center_id: int                 # 1-based atom index
    l_used_with: int               # angular momentum L this radial part is used with
    addit_r_power: int = 0         # additional r^n radial power (must be even)

    @property
    def n_prim(self) -> int:
        return len(self.exponents)


@dataclass
class BasisFunction:
    """A single (contracted) Gaussian basis function.

    Replaces BasisFunction. The linked-list _next is gone; collections of
    basis functions are plain Python lists.
    """
    l: int                         # angular momentum quantum number
    m: int                         # magnetic quantum number (MOLDEN ordering)
    center_id: int = -1            # 1-based atom index
    r0: np.ndarray | None = None   # (3,) center coordinates (view into AtomicCenter.r0)
    exponents: np.ndarray | None = None   # (n_prim,)
    coefs: np.ndarray | None = None       # (n_prim,)
    radial_part_id: int = -1       # index into RadialParts list
    additional_r_power: int = 0    # extra r^n factor
    is_spherical: bool = True
    weight: float = 0.0            # NAO weight (set during NPA)
    is_nrb: bool = False           # True if classified as Natural Rydberg Basis

    @property
    def n_prim(self) -> int:
        return 0 if self.exponents is None else len(self.exponents)

    def scale_coefs_by(self, factor: float) -> None:
        if self.coefs is not None:
            self.coefs *= factor

    def remove_small_contraction_coefs(self, threshold: float = 1e-15) -> None:
        if self.coefs is None:
            return
        mask = np.abs(self.coefs) > threshold
        self.coefs = self.coefs[mask].copy()
        self.exponents = self.exponents[mask].copy()


@dataclass
class MolecularOrbital:
    """A single molecular orbital. Replaces MO."""
    energy: float = 0.0
    occupancy: float = 0.0
    spin: int = 0                  # 0=restricted, +1=alpha, -1=beta
    bs_coefs: np.ndarray = field(default_factory=lambda: np.array([]))  # (n_basis,)