from __future__ import annotations
import os
import contextlib
import tempfile
import numpy as np
from copy import deepcopy
from warnings import warn
from typing import Tuple, List, Optional, Any

from pyscf import scf
from pyscf.tools import molden as pyscf_molden

# Import the pure Python JANPA modules
from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa
from janpa.clpo.lpo import create_clpos, CLPOOptions
from janpa.clpo.lodesc import LO_TYPE_RY, LO_TYPE_LP, LO_TYPE_BD, LO_TYPE_NB

def extract_clpo_graph(npa_res, clpo_res) -> List[Tuple[int, ...]]:
    """
    Extracts the CLPO graph (edges) directly from janpa-py results.
    Returns: List of tuples: (i,) for lone pairs, (i, i+1) for BD/NB pairs
    """
    Q = clpo_res.clpo.nao_to_hybrids @ clpo_res.clpo.lo_to_hybrids.T
    sds_clpo = Q.T @ npa_res.sds_nao @ Q
    occupancies = np.diag(sds_clpo)
    lo_types = clpo_res.clpo.lo_types
    
    nodes = []
    n_lo = len(occupancies)
    i = 0
    while i < n_lo:
        lo_type = lo_types[i]
        occ = occupancies[i]
        
        if lo_type == LO_TYPE_LP:
            if 0.5 < occ < 1.5:
                warn(f'Lone Pair {i} found with occupation close to 1 => {occ:.5f}, take care.')
            nodes.append((i,))
            i += 1
        elif lo_type == LO_TYPE_RY:
            i += 1
        elif lo_type == LO_TYPE_BD:
            if i + 1 >= n_lo or lo_types[i+1] != LO_TYPE_NB:
                raise ValueError(f"BD orbital {i} without following NB orbital")
            bd_occ = occ
            nb_occ = occupancies[i+1]
            if bd_occ + nb_occ < 1.7:
                warn(f'Bond pair orbitals [{i},{i+1}] population expected under expected 2e-, predicted: {bd_occ + nb_occ:.5f}, take care with predicted edges.')
            nodes.append((i, i + 1))
            i += 2
        else:
            i += 1
    return nodes

def _run_janpa_pipeline(mol, mf, thres=1e-9, silent=True, custom_edges=None):
    """
    Internal helper to run the JANPA pipeline entirely in memory/tempfiles.
    """
    # Use a temporary directory so we don't clutter the user's workspace
    with tempfile.TemporaryDirectory() as tmpdir:
        molden_in = os.path.join(tmpdir, "temp.molden")
        
        # Set up stdout suppression
        devnull = open(os.devnull, 'w')
        ctx = contextlib.redirect_stdout(devnull) if silent else contextlib.nullcontext()
        
        try:
            with ctx:
                # 1. Export initial MOLDEN
                pyscf_molden.from_mo(mol, molden_in, mf.mo_coeff, occ=mf.mo_occ, ene=mf.mo_energy)
                
                # 2. Load and prepare
                molden_file = MoldenFile.load(molden_in)
                molden_file.coords_to_au()
                molden_file.to_unnormalized_primitive_coefs()
                
                # 3. NPA
                npa_res = run_npa(molden_file)
                
                # 4. CLPO
                clpo_opts = CLPOOptions(hybr_opt_conv_thresh=thres, hybr_opt_max_iter=1000)
                if custom_edges is not None:
                    # JANPA's parse_edges expects a string like "[(1, 2), (3, 4)]"
                    clpo_opts.edges = str(custom_edges) 
                
                clpo_res = create_clpos(npa_res.sds_nao, npa_res.nao, molden_file.centers, clpo_opts)
                
                # 5. Extract Graph
                graph = extract_clpo_graph(npa_res, clpo_res)
                
                # 6. Get Transformation Matrices
                aho_to_ao = clpo_res.clpo.nao_to_hybrids.T @ npa_res.nao_to_ao
                clpo_to_ao = clpo_res.clpo.lo_to_hybrids @ aho_to_ao
                
        finally:
            if silent:
                devnull.close()
                
    return aho_to_ao, clpo_to_ao, graph

def generate_HAO_molecule(mol, mf=None, silent=True, thres=1e-9, **kwargs):
    """
    Generates a molecule with Hybrid Atomic Orbitals (HAO/AHO) via janpa-py.
    """
    if mf is None:
        mf = scf.RHF(mol).run()
        
    aho_to_ao, _, _ = _run_janpa_pipeline(mol, mf, thres=thres, silent=silent)
    
    # Compatibility with Tequila/Sunrise Molecule objects
    if hasattr(mol, 'integral_manager'):
        nmol = deepcopy(mol)
        nmol.integral_manager.orbital_coefficients = aho_to_ao
        nmol.integral_manager._orbital_type = "HAO"
        return nmol
        
    # Fallback for pure PySCF
    return aho_to_ao

def generate_CLPO_molecule_edges(mol, mf=None, edges=None, silent=True, thres=1e-12, **kwargs):
    """
    Generates a molecule with CLPO orbitals and returns the molecule and SPA edges.
    """
    if mf is None:
        mf = scf.RHF(mol).run()
        
    _, clpo_to_ao, graph = _run_janpa_pipeline(mol, mf, thres=thres, silent=silent, custom_edges=edges)
    
    if hasattr(mol, 'integral_manager'):
        nmol = deepcopy(mol)
        nmol.integral_manager.orbital_coefficients = clpo_to_ao
        nmol.integral_manager._orbital_type = "CLPO"
        return nmol, graph
        
    return clpo_to_ao, graph

def generate_CLPO_molecule(mol, mf=None, edges=None, silent=True, thres=1e-12, **kwargs):
    """
    Generates a molecule with CLPO orbitals.
    """
    res = generate_CLPO_molecule_edges(mol, mf, edges, silent, thres, **kwargs)
    if isinstance(res, tuple):
        return res[0]
    return res