# janpa/utils/printout.py
"""Formatting and exporting matrices and orbitals to MOLDEN format. Replaces onpa/printout.java."""

from __future__ import annotations
import copy
import numpy as np
from janpa.gto.basis import MolecularOrbital
from janpa.io.molden import MoldenFile

STARS = " * * * "


def print_stars() -> None:
    print()
    print(STARS)
    print()


def export_orbitals(
    ao_info: MoldenFile,
    orbitals_to_ao: np.ndarray,
    occupancies: np.ndarray | None,
    energies: np.ndarray | None,
    file_name: str,
    renorm_bf_coefs: bool,
    title: str | None = None,
) -> None:
    """Exports orbitals into a MOLDEN file."""
    molden2 = copy.deepcopy(ao_info)
    molden2.is_spherical = True
    if title is not None:
        molden2.title = title
    else:
        molden2.title = "Natural Orbitals prepared by JANPA"
        
    n_bfs = orbitals_to_ao.shape[1]
    n_orbitals = orbitals_to_ao.shape[0]
    
    mos = []
    for mo_idx in range(n_orbitals):
        occ = occupancies[mo_idx] if occupancies is not None else 0.0
        ene = energies[mo_idx] if energies is not None else 0.0
        spin = 1
        bs_coefs = orbitals_to_ao[mo_idx, :].copy()
        mos.append(MolecularOrbital(energy=ene, occupancy=occ, spin=spin, bs_coefs=bs_coefs))
        
    molden2.mos = mos
    molden2.save(file_name)


def print_matrix(
    a: np.ndarray,
    comment: str,
    column_names: list[str] | None,
    row_names: list[str] | None,
    fname: str,
    matrix_float_number_format: str = "%12.5f",
    matrix_line_width: int = 0,
) -> None:
    """Formats and exports a matrix to a plain text file."""
    if not fname:
        return
        
    print(f'\nExporting matrix "{comment}" to {fname}')
    
    with open(fname, 'w') as f:
        if comment is not None:
            f.write(comment + "\n")
        else:
            f.write("\n")
            
        f.write(f"{a.shape[0]}\t{a.shape[1]}\n")
        
        if column_names is None:
            f.write("\n")
        else:
            for i, name in enumerate(column_names):
                if matrix_line_width != 0 and i > 0 and (i % matrix_line_width == 0):
                    f.write("\n")
                f.write(f"{name}\t")
            f.write("\n")
            
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                if (j == 0) or (matrix_line_width != 0 and (j % matrix_line_width == 0)):
                    f.write("\n")
                try:
                    val_str = matrix_float_number_format % a[i, j]
                except TypeError:
                    val_str = f"{a[i, j]:{matrix_float_number_format}}"
                f.write(f"{val_str}\t")
                
            if row_names is None:
                f.write("\t")
            else:
                f.write(f"\t{row_names[i]}")
            f.write("\n")
            
    print("done")