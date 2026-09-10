# janpa/clpo/hybridization.py
"""Hybridization analysis of LHOs. Replaces CLPO/HybridizAn.java."""
from __future__ import annotations
import numpy as np
from .lpo import LODescription
from janpa.gto.basis import BasisFunction, AtomicCenter

def hybr_an(clpo: LODescription, naos: list[BasisFunction], centers: list[AtomicCenter]):
    """Print hybridization analysis of LHOs."""
    print("Hydridization analysis of LHOs:")
    l_max = 4
    hybridization = np.zeros((len(clpo.host_atom_of_hybrid), l_max + 1))
    
    for i in range(len(clpo.host_atom_of_hybrid)):
        for nao_idx in range(len(naos)):
            tmp = clpo.nao_to_hybrids[nao_idx, i]
            hybridization[i, naos[nao_idx].l] += tmp * tmp
            
        a_idx = clpo.host_atom_of_hybrid[i]
        print(f"{i + 1:4d} {centers[a_idx].name}{a_idx + 1}", end="")
        for l in range(l_max + 1):
            print(f"{hybridization[i, l] * 100:7.2f} ", end="")
            
        s_char = hybridization[i, 0]
        p_ratio = hybridization[i, 1] / s_char if s_char > 0 else 0.0
        d_ratio = hybridization[i, 2] / s_char if s_char > 0 else 0.0
        
        print(f"s p ^ ({p_ratio:.1f}) d ^ ({d_ratio:.1f})")
    print()