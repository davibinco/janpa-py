from __future__ import annotations
import os
import contextlib
import tempfile
import numpy as np
from copy import deepcopy

import psi4

# Import the pure Python JANPA modules
from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa
from janpa.clpo.lpo import create_clpos, CLPOOptions
from janpa.clpo.graph import extract_clpo_graph_from_janpa


def _run_scf(mol):
    """Runs the SCF on a psi4 Molecule with the current psi4 options (basis, reference, ...)."""
    _, wfn = psi4.energy("scf", molecule=mol, return_wfn=True)
    return wfn


def _molden_to_psi4_ao_order(wfn) -> np.ndarray:
    """
    perm[k] is the psi4 AO index of the k-th MOLDEN basis function. MOLDEN orders the p shells as x, y, z while the
    pure psi4 shells are p0, p+1, p-1 = z, x, y. The other shells are ordered alike (0, +1, -1, +2, -2, ...).
    """
    bs = wfn.basisset()
    perm = []
    for sh in range(bs.nshell()):
        shell = bs.shell(sh)
        order = [1, 2, 0] if shell.am == 1 and bs.has_puream() else range(shell.nfunction)
        perm += [shell.function_index + i for i in order]
    return np.array(perm)


def _run_janpa_pipeline(wfn, thres=1e-9, silent=True, custom_edges=None):
    """
    Internal helper to run the JANPA pipeline on a psi4 Wavefunction. The returned matrices have one row per orbital and
    one column per psi4 AO.
    """
    bs = wfn.basisset()
    if not bs.has_puream() and bs.max_am() >= 2:
        raise ValueError("JANPA needs spherical basis functions, set the psi4 option 'puream' to True")

    with tempfile.TemporaryDirectory() as tmpdir:
        molden_in = os.path.join(tmpdir, "temp.molden")

        devnull = open(os.devnull, 'w')
        ctx = contextlib.redirect_stdout(devnull) if silent else contextlib.nullcontext()

        try:
            with ctx:
                # 1. Export initial MOLDEN
                psi4.molden(wfn, molden_in)

                # 2. Load and prepare. psi4 writes the raw contraction coefs, which are not normalized for general
                # contractions (e.g. cc-pVXZ), while its MO coefs refer to normalized basis functions
                molden_file = MoldenFile.load(molden_in)
                molden_file.coords_to_au()
                molden_file.to_unnormalized_primitive_coefs()
                molden_file.normalize_contractions()

                # 3. NPA
                npa_res = run_npa(molden_file)

                # 4. CLPO
                clpo_opts = CLPOOptions(hybr_opt_conv_thresh=thres, hybr_opt_max_iter=1000)
                if custom_edges is not None:
                    # Bonds between 0-based atoms, a string like "[(1, 2), (1, 2)]"
                    clpo_opts.edges = str(custom_edges)

                clpo_res = create_clpos(npa_res.sds_nao, npa_res.nao, molden_file.centers, clpo_opts)

                # 5. Extract Graph
                graph = extract_clpo_graph_from_janpa(npa_res, clpo_res)

                # 6. Get Transformation Matrices
                aho_to_ao = clpo_res.clpo.nao_to_hybrids.T @ npa_res.nao_to_ao
                clpo_to_ao = clpo_res.clpo.lo_to_hybrids @ aho_to_ao

        finally:
            devnull.close()

    perm = _molden_to_psi4_ao_order(wfn)
    aho_psi4 = np.zeros_like(aho_to_ao)
    aho_psi4[:, perm] = aho_to_ao
    clpo_psi4 = np.zeros_like(clpo_to_ao)
    clpo_psi4[:, perm] = clpo_to_ao
    return aho_psi4, clpo_psi4, graph


def generate_HAO_molecule(mol=None, wfn=None, silent=True, thres=1e-9, **kwargs):
    """
    Generates a molecule with Hybrid Atomic Orbitals (HAO/AHO) via janpa-py from a psi4 calculation.
    mol is a psi4 Molecule, used to run the SCF with the current psi4 options when wfn (a psi4 Wavefunction) is not given.
    """
    if wfn is None:
        wfn = _run_scf(mol)

    aho_to_ao, _, _ = _run_janpa_pipeline(wfn, thres=thres, silent=silent)

    # Compatibility with Tequila/Sunrise Molecule objects
    if hasattr(mol, 'integral_manager'):
        nmol = deepcopy(mol)
        nmol.integral_manager.orbital_coefficients = aho_to_ao
        nmol.integral_manager._orbital_type = "HAO"
        return nmol

    # Fallback for pure psi4
    return aho_to_ao


def generate_CLPO_molecule_edges(mol=None, wfn=None, edges=None, silent=True, thres=1e-12, **kwargs):
    """
    Generates a molecule with CLPO orbitals from a psi4 calculation and returns the molecule and SPA edges.
    mol is a psi4 Molecule, used to run the SCF with the current psi4 options when wfn (a psi4 Wavefunction) is not given.
    edges optionally imposes bonds on the CLPOs, between 0-based atoms, one pair per bond, e.g. [(1, 2), (1, 2)] for a
    double bond between atoms 1 and 2. The bonds between the other atom pairs are found as usual.
    """
    if wfn is None:
        wfn = _run_scf(mol)

    _, clpo_to_ao, graph = _run_janpa_pipeline(wfn, thres=thres, silent=silent, custom_edges=edges)

    if hasattr(mol, 'integral_manager'):
        nmol = deepcopy(mol)
        nmol.integral_manager.orbital_coefficients = clpo_to_ao
        nmol.integral_manager._orbital_type = "CLPO"
        return nmol, graph

    return clpo_to_ao, graph


def generate_CLPO_molecule(mol=None, wfn=None, edges=None, silent=True, thres=1e-12, **kwargs):
    """
    Generates a molecule with CLPO orbitals from a psi4 calculation.
    """
    res = generate_CLPO_molecule_edges(mol, wfn, edges, silent, thres, **kwargs)
    if isinstance(res, tuple):
        return res[0]
    return res
