# janpa/clpo/lodesc.py
"""Localized orbital description. Replaces LOdescription.java."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

LO_TYPE_RY = 0
LO_TYPE_LP = 1
LO_TYPE_BD = 2
LO_TYPE_NB = 3

@dataclass
class LODescription:
    n_nao: int
    nao_to_hybrids: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    lo_to_hybrids: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    lo_labels: list[str] = field(default_factory=list)
    lo_types: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    hybrid_labels: list[str] = field(default_factory=list)
    host_atom_of_hybrid: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    hybrids_of_lo: list[list[int]] = field(default_factory=list)
    bd_per_atomic_pair: np.ndarray | None = None

    def __post_init__(self):
        n = self.n_nao
        if self.nao_to_hybrids.size == 0:
            self.nao_to_hybrids = np.zeros((n, n))
        if self.lo_to_hybrids.size == 0:
            self.lo_to_hybrids = np.zeros((n, n))
        if not self.lo_labels:
            self.lo_labels = [""] * n
        if self.lo_types.size == 0:
            self.lo_types = np.zeros(n, dtype=int)
        if self.host_atom_of_hybrid.size == 0:
            self.host_atom_of_hybrid = np.zeros(n, dtype=int)
        if not self.hybrids_of_lo:
            self.hybrids_of_lo = [[] for _ in range(n)]