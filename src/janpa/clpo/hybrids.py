# janpa/clpo/hybrids.py
"""Atomic hybrids management. Replaces AtomicHybrids.java."""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from janpa.utils.matrix import symmetric_orthogonalize

@dataclass
class AtomicHybrids:
    nao_indices: np.ndarray
    u: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    u_saved: np.ndarray | None = None
    n_valid_hybrids: int = 0
    friend_atom_index: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    friend_hybrid_index: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    friend_refs_saved: list[np.ndarray] = field(default_factory=list)

    def __post_init__(self):
        n = len(self.nao_indices)
        if self.u.size == 0:
            self.u = np.zeros((n, n))
        if self.friend_atom_index.size == 0:
            self.friend_atom_index = np.full(n, -1, dtype=int)
        if self.friend_hybrid_index.size == 0:
            self.friend_hybrid_index = np.full(n, -1, dtype=int)
        self.friend_refs_saved = [np.array([], dtype=int), np.array([], dtype=int)]

    def backup_u(self):
        self.u_saved = self.u.copy()

    def restore_u(self):
        if self.u_saved is not None:
            self.u = self.u_saved.copy()

    def backup_bonding(self):
        self.friend_refs_saved[0] = self.friend_atom_index.copy()
        self.friend_refs_saved[1] = self.friend_hybrid_index.copy()

    def restore_bonding(self):
        if self.friend_refs_saved[0].size > 0:
            self.friend_atom_index = self.friend_refs_saved[0].copy()
        if self.friend_refs_saved[1].size > 0:
            self.friend_hybrid_index = self.friend_refs_saved[1].copy()

    def get_hybrid(self, i: int) -> np.ndarray:
        return self.u[:, i]

    def hybrid_scalar_mul(self, i: int, v: np.ndarray) -> float:
        return float(np.dot(self.u[:, i], v[:, 0]))

    def append_hybrid(self, vec: np.ndarray) -> bool:
        if self.n_valid_hybrids == len(self.nao_indices):
            return False
        self.u[:, self.n_valid_hybrids] = vec
        self.n_valid_hybrids += 1
        return True

    def test_lin_indep(self, eig_s_threshold: float, eigenvals_s: list[np.ndarray] | None = None) -> bool:
        if self.n_valid_hybrids <= 1 and eigenvals_s is None:
            return True
        u_small = self.u[:, :self.n_valid_hybrids]
        s = u_small.T @ u_small
        if self.n_valid_hybrids <= 1 and eigenvals_s is not None:
            eigenvals_s[0] = np.array([s[0, 0]])
            return True
        vals, _ = np.linalg.eigh(s)
        vals = np.sort(vals)
        if eigenvals_s is not None:
            eigenvals_s[0] = vals
        return vals[-1] > eig_s_threshold

    def append_hybrid_if_lin_indep(self, eig_s_threshold: float, vec: np.ndarray) -> bool:
        if not self.append_hybrid(vec):
            return False
        if self.test_lin_indep(eig_s_threshold, None):
            return True
        else:
            self.n_valid_hybrids -= 1
            return False

    def symmetr_orth(self):
        if self.n_valid_hybrids <= 1:
            return
        u_small = self.u[:, :self.n_valid_hybrids]
        self.u[:, :self.n_valid_hybrids] = symmetric_orthogonalize(u_small)