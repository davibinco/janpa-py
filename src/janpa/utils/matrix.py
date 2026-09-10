"""Thin wrappers around numpy/scipy replacing EigenEngine.java. Most operations are one-liners. The only non-trivial algorithm is CS96_SimultDiag (simultaneous diagonalization by Jacobi rotations). """
from __future__ import annotations

import numpy as np
from scipy import linalg

import numpy as np

def matrix_sqrt(A: np.ndarray, inverse: bool = False) -> np.ndarray:
    eigvals, eigvecs = np.linalg.eigh(A)
    
    if inverse:
        # Match original Java: absolute threshold of 1e-12
        eigvals_inv_sqrt = np.zeros_like(eigvals)
        mask = eigvals > 1e-12
        eigvals_inv_sqrt[mask] = 1.0 / np.sqrt(eigvals[mask])
        return eigvecs @ np.diag(eigvals_inv_sqrt) @ eigvecs.T
    else:
        eigvals_sqrt = np.sqrt(np.maximum(eigvals, 0.0))
        return eigvecs @ np.diag(eigvals_sqrt) @ eigvecs.T

def transform_symmetric_to_new_basis(a: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Q^T A Q for symmetric A. Exploits symmetry for efficiency."""
    return q.T @ a @ q

def transform_to_new_basis(a: np.ndarray, u: np.ndarray, is_symmetric: bool) -> np.ndarray:
    """Transform matrix A to new basis U.
    If is_symmetric is True, uses Q^T A Q.
    If is_symmetric is False, uses U A U^T.
    """
    if is_symmetric:
        return transform_symmetric_to_new_basis(a, u)
    else:
        return u @ a @ u.T

def generalized_eigh(a: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Solve A v = λ S v for symmetric A, S. Returns (eigenvalues, eigenvectors) with eigenvectors S-orthonormal.
    """
    return linalg.eigh(a, s)

def symmetric_orthogonalize(vecs: np.ndarray) -> np.ndarray:
    """Symmetric orthogonalization: U S V^T → U V^T (from SVD of Vecs)."""
    u, _, vt = linalg.svd(vecs, full_matrices=False)
    return u @ vt

# ---------------------------------------------------------
# Sorter methods (trivially replaced by numpy)
# ---------------------------------------------------------

def merge_ascending(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Merge two arrays and sort ascending."""
    return np.sort(np.concatenate((a, b)))

def array_sort(array: np.ndarray, ascending: bool = True) -> np.ndarray:
    """Return indices that would sort the array."""
    if ascending:
        return np.argsort(array)
    else:
        return np.argsort(array)[::-1]

def array_remap(array: np.ndarray, remap: np.ndarray) -> np.ndarray:
    """Reorder array according to remap indices."""
    return array[remap]

def array_sort_in_place(array: np.ndarray, ascending: bool = True) -> np.ndarray:
    """Sort array in place and return the remap indices."""
    remap = array_sort(array, ascending)
    array[:] = array[remap]
    return remap

def sort_eigenvalues(eigenvalues: np.ndarray, descending: bool = True) -> np.ndarray:
    """Sort eigenvalues and return the sorting indices.
    Defaults to descending order to match the Java Sorter.sort_eigenvalues behavior.
    """
    if descending:
        return np.argsort(eigenvalues)[::-1]
    else:
        return np.argsort(eigenvalues)

def sorted_eigen_data(eigenvalues: np.ndarray, eigenvectors: np.ndarray, descending: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Sort eigenvalues and eigenvectors.
    Defaults to descending order to match the Java Sorter.sorted_eigenData behavior.
    """
    remap = sort_eigenvalues(eigenvalues, descending)
    return eigenvalues[remap], eigenvectors[:, remap]