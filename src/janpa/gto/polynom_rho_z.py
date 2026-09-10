"""2D polynomial in z and rho for cylindrical integration.

Replaces Polynom_rho_z.java.

The Java linked stack based on `stackPointer` is replaced by an explicit
Python list of saved states. Each saved state contains a copy of the
coefficient matrix and the active powers.
"""

from __future__ import annotations

import numpy as np

__all__ = ["Polynom_rho_z"]


class Polynom_rho_z:
    """Polynomial in z and rho with specialized multiplication routines."""

    def __init__(self, max_z_pwr: int, max_rho_pwr: int) -> None:
        max_z_pwr = int(max_z_pwr)
        max_rho_pwr = int(max_rho_pwr)

        if max_z_pwr < 0 or max_rho_pwr < 0:
            raise ValueError("max_z_pwr and max_rho_pwr must be non-negative.")

        self.cf_matrix = np.zeros((max_z_pwr + 1, max_rho_pwr + 1), dtype=float)
        self.tmp_matrix = np.zeros((max_z_pwr + 1, max_rho_pwr + 1), dtype=float)

        self.cf_matrix[0, 0] = 1.0

        self.last_active_rho_power = 0
        self.last_active_z_power = 0

        self._stack: list[tuple[np.ndarray, int, int]] = []

    # ------------------------------------------------------------------
    # Stack mechanism
    # ------------------------------------------------------------------
    def push_to_stack(self) -> None:
        """Save current cf_matrix and active powers on an internal stack."""
        self._stack.append(
            (
                self.cf_matrix.copy(),
                int(self.last_active_rho_power),
                int(self.last_active_z_power),
            )
        )

    def restore_from_stack(self, delete_source_from_stack: bool) -> None:
        """Restore state from stack.

        Parameters
        ----------
        delete_source_from_stack:
            If True, pop the saved state. If False, keep it and copy it.
        """
        if not self._stack:
            return

        if delete_source_from_stack:
            cf_matrix, rho_power, z_power = self._stack.pop()
        else:
            cf_matrix, rho_power, z_power = self._stack[-1]

        self.cf_matrix = cf_matrix.copy()
        self.last_active_rho_power = int(rho_power)
        self.last_active_z_power = int(z_power)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _swap_old_new(self) -> None:
        """Swap cf_matrix and tmp_matrix."""
        self.cf_matrix, self.tmp_matrix = self.tmp_matrix, self.cf_matrix

    # ------------------------------------------------------------------
    # Multiplication routines
    # ------------------------------------------------------------------
    def mul_sq(self, cf_z: float, cf_rho: float, z_a_minus_z0: float) -> None:
        """Multiply by cf_z * (z - z_a_minus_z0)^2 + cf_rho * rho^2.

        This is a direct translation of the Java loop structure, including
        the temporary-matrix initialization pattern.
        """
        lz = self.last_active_z_power
        lr = self.last_active_rho_power

        cf = self.cf_matrix
        tmp = self.tmp_matrix

        # tmp[i+2, j] = cf_z * cf[i, j]
        for i in range(lz + 1):
            for j in range(lr + 1):
                tmp[i + 2, j] = cf_z * cf[i, j]

        # Zero first two rows over the future active rho range.
        for j in range(lr + 3):
            tmp[0, j] = 0.0
            tmp[1, j] = 0.0

        # Zero the two future active rho columns.
        for i in range(lz + 3):
            tmp[i, lr + 1] = 0.0
            tmp[i, lr + 2] = 0.0

        # -2 cf_z z_a_minus_z0 z
        tmp_coeff = cf_z * 2.0 * z_a_minus_z0
        for i in range(lz + 1):
            for j in range(lr + 1):
                tmp[i + 1, j] -= tmp_coeff * cf[i, j]

        # cf_rho rho^2
        for i in range(lz + 1):
            for j in range(lr + 1):
                tmp[i, j + 2] += cf_rho * cf[i, j]

        # cf_z z_a_minus_z0^2
        tmp_coeff = cf_z * z_a_minus_z0 * z_a_minus_z0
        for i in range(lz + 1):
            for j in range(lr + 1):
                tmp[i, j] += tmp_coeff * cf[i, j]

        self.last_active_rho_power += 2
        self.last_active_z_power += 2
        self._swap_old_new()

    def mul_lin(self, z_a_minus_z0: float) -> None:
        """Multiply by (z - z_a_minus_z0)."""
        lz = self.last_active_z_power
        lr = self.last_active_rho_power
        cf = self.cf_matrix

        # Copy top active row to the new top row.
        cf[lz + 1, : lr + 1] = cf[lz, : lr + 1]

        # Descend so that cf[i-1] is still unchanged when used.
        for i in range(lz, 0, -1):
            for j in range(lr + 1):
                cf[i, j] *= -z_a_minus_z0
                cf[i, j] += cf[i - 1, j]

        for j in range(lr + 1):
            cf[0, j] *= -z_a_minus_z0

        self.last_active_z_power += 1

    def mul_assoc_leg(self, L: int, m: int, z_a_minus_z0: float) -> None:
        """Multiply by associated Legendre polynomial factor roots."""
        from .spherical import ALEG_ROOTS_SQUARED

        L = int(L)
        abs_m = abs(int(m))

        squared_roots = ALEG_ROOTS_SQUARED[L][abs_m]
        for root_sq in squared_roots:
            self.mul_sq(1.0 - root_sq, -root_sq, z_a_minus_z0)

        if (L - abs_m) % 2 == 1:
            self.mul_lin(z_a_minus_z0)

    def mul_r2(self, r2_pwr: int, z_a_minus_z0: float) -> None:
        """Multiply by r^2 = rho^2 + (z - z_a_minus_z0)^2, r2_pwr times."""
        r2_pwr = int(r2_pwr)
        for _ in range(r2_pwr):
            self.mul_sq(1.0, 1.0, z_a_minus_z0)

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    def printout(self) -> None:
        """Print non-negligible terms, Java-compatible."""
        for i in range(self.cf_matrix.shape[0]):
            something_printed = False
            for j in range(self.cf_matrix.shape[1]):
                if abs(self.cf_matrix[i, j]) > 1.0e-10:
                    print(
                        f" {self.cf_matrix[i, j]:+.3f} * z^{i} * rho^{j} ",
                        end="",
                    )
                    something_printed = True
            if something_printed:
                print()
        print("---")

    # Java-compatible aliases
    pushToStack = push_to_stack
    restoreFromStack = restore_from_stack
    mul_AssocLeg = mul_assoc_leg