# janpa/utils/constants.py
"""Physical constants and periodic table data."""

import numpy as np

BOHR_RADIUS = 0.52917721092  # Angstrom per Bohr
ANGSTROM_TO_BOHR = 1.0 / BOHR_RADIUS

# Periodic table: index 0 → Z=1 (H)
ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn",
]

_ELEMENT_TO_Z = {name: i + 1 for i, name in enumerate(ELEMENTS)}

def name_to_z(atom_name: str) -> int:
    """Convert element name/symbol to atomic number Z. Case-insensitive."""
    z = _ELEMENT_TO_Z.get(atom_name.strip().capitalize())
    if z is None:
        raise ValueError(f"Unknown atom type: '{atom_name}'")
    return z

def z_to_name(z: int) -> str:
    """Convert atomic number Z to element symbol."""
    if z < 1 or z > len(ELEMENTS):
        raise ValueError(f"Z={z} out of range")
    return ELEMENTS[z - 1]

# Angular momentum label helpers
SPDF_LABELS = "spdfgh"