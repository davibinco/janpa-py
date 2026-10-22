"""Simultaneous diagonalization of a set of symmetric matrices by Jacobi rotations. Cardellini & Sauer (1996) algorithm. Replaces CS96_SimultDiag.java. This is the ONE matrix algorithm with no direct scipy equivalent. """
from __future__ import annotations

import math

import numpy as np

def simultaneous_diagonalization(
    matrices: list[np.ndarray],
    max_iters: int = -1,
    u0: np.ndarray | None = None,
    threshold: float = 1e-7,
) -> tuple[np.ndarray, int]:
    """Simultaneously diagonalize a set of symmetric matrices.

    Parameters
    ----------
    matrices : list of (n, n) symmetric arrays (modified in-place)
    max_iters : maximum iterations (-1 for unlimited)
    u0 : initial rotation matrix (n, n), or None for identity
    threshold : convergence threshold on max |sin(theta)|

    Returns
    -------
    U : (n, n) rotation matrix such that U^T A_k U ≈ diagonal for all k
    n_iters : number of iterations performed
    """
    if not matrices:
        raise ValueError("matrices list cannot be empty")

    # All matrices are rotated together as one (K, n, n) stack, so each Jacobi rotation costs a few numpy calls
    # instead of a few per matrix.
    a = np.array(matrices, dtype=float)
    n = a.shape[1]

    if u0 is None:
        U = np.eye(n)
    else:
        # Match Java logic: result = u0.transpose()
        U = u0.T.copy()

    n_iter = 0

    while True:
        s_max = 0.0

        # Single sweep iteration
        for i in range(n):
            for j in range(i + 1, n):
                # Evaluate angles required to simultaneously annihilate off-diagonals
                h0 = a[:, i, i] - a[:, j, j]
                h1 = a[:, i, j] + a[:, j, i]
                t_on = float(np.dot(h0, h0) - np.dot(h1, h1))
                t_off = 2.0 * float(np.dot(h0, h1))
                theta = math.atan2(t_off, t_on) / 4.0
                c = math.cos(theta)
                s = math.sin(theta)

                if abs(s) > s_max:
                    s_max = abs(s)
                if s == 0.0:
                    continue  # identity rotation

                # Apply computed rotation R = [[c, s], [-s, c]] as A -> R A R^T to all matrices
                row_i = a[:, i, :].copy()
                row_j = a[:, j, :].copy()
                a[:, i, :] = c * row_i + s * row_j
                a[:, j, :] = c * row_j - s * row_i
                col_i = a[:, :, i].copy()
                col_j = a[:, :, j].copy()
                a[:, :, i] = c * col_i + s * col_j
                a[:, :, j] = c * col_j - s * col_i

                # Update transformation matrix U (multiplies on the left by rotation)
                U_i = U[i, :].copy()
                U_j = U[j, :].copy()
                U[i, :] = c * U_i + s * U_j
                U[j, :] = c * U_j - s * U_i

        n_iter += 1

        # Break conditions equivalent to Java's do-while
        if s_max <= threshold:
            break
        if max_iters >= 0 and n_iter >= max_iters:
            break

    for k, m in enumerate(matrices):
        m[...] = a[k]

    # Return transpose to match Java `return result.transpose()`
    return U.T, n_iter
