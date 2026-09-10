"""Simultaneous diagonalization of a set of symmetric matrices by Jacobi rotations. Cardellini & Sauer (1996) algorithm. Replaces CS96_SimultDiag.java. This is the ONE matrix algorithm with no direct scipy equivalent. """
from __future__ import annotations

import numpy as np

def _2x2_rotate(a: np.ndarray, i: int, j: int, c: float, s: float) -> None:
    """Apply a 2x2 Jacobi rotation to rows/cols i and j of matrix a (in-place)."""
    aii = a[i, i]
    ajj = a[j, j]
    aij = a[i, j]
    aji = a[j, i]

    # Temporarily copy rows & columns
    row_i = a[i, :].copy()
    row_j = a[j, :].copy()
    
    # Apply rotation over rows
    a[i, :] = c * row_i + s * row_j
    a[j, :] = -s * row_i + c * row_j

    col_i = a[:, i].copy()
    col_j = a[:, j].copy()
    
    # Apply rotation over columns
    a[:, i] = c * col_i + s * col_j
    a[:, j] = -s * col_i + c * col_j

    # Fix the intersecting 2x2 block that was inadvertently modified 
    # by the overlapping row & column transformations.
    a[i, i] = c*c*aii + s*s*ajj + c*s*(aij + aji)
    a[j, j] = c*c*ajj + s*s*aii - c*s*(aij + aji)
    a[i, j] = c*c*aij - s*s*aji + c*s*(ajj - aii)
    a[j, i] = c*c*aji - s*s*aij + c*s*(ajj - aii)

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
        
    n = matrices[0].shape[0]
    
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
                t_on = 0.0
                t_off = 0.0
                
                # Evaluate angles required to simultaneously annihilate off-diagonals
                for k in range(len(matrices)):
                    a = matrices[k]
                    h0 = a[i, i] - a[j, j]
                    h1 = a[i, j] + a[j, i]
                    t_on += h0 * h0 - h1 * h1
                    t_off += h0 * h1
                
                t_off *= 2.0
                theta = np.arctan2(t_off, t_on) / 4.0
                c = np.cos(theta)
                s = np.sin(theta)
                
                if np.abs(s) > s_max:
                    s_max = np.abs(s)
                    
                # Apply computed rotation to all matrices in the list
                for k in range(len(matrices)):
                    _2x2_rotate(matrices[k], i, j, c, s)
                    
                # Update transformation matrix U (multiplies on the left by rotation)
                U_i = U[i, :].copy()
                U_j = U[j, :].copy()
                U[i, :] = c * U_i + s * U_j
                U[j, :] = -s * U_i + c * U_j
                
        n_iter += 1
        
        # Break conditions equivalent to Java's do-while
        if s_max <= threshold:
            break
        if max_iters >= 0 and n_iter >= max_iters:
            break
            
    # Return transpose to match Java `return result.transpose()`
    return U.T, n_iter