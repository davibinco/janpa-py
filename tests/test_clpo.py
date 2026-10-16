import pytest
import numpy as np
from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa
from janpa.clpo.lpo import create_clpos, CLPOOptions
from pyscf import gto, scf
from janpa.clpo.lodesc import LO_TYPE_BD, LO_TYPE_NB, LO_TYPE_LP
from janpa.janpa_interface import generate_CLPO_molecule_edges

def test_h4_clpo_bonds(h4_molden_path):
    """Test that H4 forms two distinct H-H bonds in CLPO."""
    mf = MoldenFile.load(h4_molden_path)
    mf.coords_to_au()
    mf.to_unnormalized_primitive_coefs()
    
    npa_res = run_npa(mf)
    
    opts = CLPOOptions(hybr_opt_conv_thresh=1e-9, hybr_opt_max_iter=1000)
    clpo_res = create_clpos(npa_res.sds_nao, npa_res.nao, mf.centers, opts)
    
    # We expect exactly 2 bonding orbitals (BD) and 2 antibonding (NB) for H4
    bd_count = np.sum(clpo_res.clpo.lo_types == LO_TYPE_BD)
    nb_count = np.sum(clpo_res.clpo.lo_types == LO_TYPE_NB)
    
    assert bd_count == 2, f"Expected 2 BD orbitals, found {bd_count}"
    assert nb_count == 2, f"Expected 2 NB orbitals, found {nb_count}"

def _clpo_bonds(molden_path, edges):
    mf = MoldenFile.load(molden_path)
    mf.coords_to_au()
    mf.to_unnormalized_primitive_coefs()
    npa_res = run_npa(mf)
    opts = CLPOOptions(hybr_opt_conv_thresh=1e-9, hybr_opt_max_iter=1000, edges=edges)
    clpo = create_clpos(npa_res.sds_nao, npa_res.nao, mf.centers, opts).clpo
    return sorted(tuple(sorted(clpo.host_atom_of_hybrid[h] for h in clpo.hybrids_of_lo[i]))
                  for i in range(len(clpo.lo_types)) if clpo.lo_types[i] == LO_TYPE_BD)


def test_h4_imposed_bonds(h4_molden_path):
    """Bonds given on atoms are kept, even when they are not the best ones, and the other atom pairs are found as usual."""
    assert _clpo_bonds(h4_molden_path, "") == [(0, 1), (2, 3)]
    assert _clpo_bonds(h4_molden_path, "[(0, 1)]") == [(0, 1), (2, 3)]
    assert _clpo_bonds(h4_molden_path, "[(1, 2)]") == [(1, 2)]
    assert _clpo_bonds(h4_molden_path, str([(0, 3), (1, 2)])) == [(0, 3), (1, 2)]


def test_lewis_structure_validation(h4_molden_path):
    with pytest.raises(ValueError, match="does not exist"):
        _clpo_bonds(h4_molden_path, "[(0, 4)]")
    with pytest.raises(ValueError, match="only 1 hybrids"):
        _clpo_bonds(h4_molden_path, "[(0, 1), (0, 1)]")
    with pytest.raises(ValueError, match="only bonds are given"):
        _clpo_bonds(h4_molden_path, "[(0,)]")


def test_h2co_clpo(h2co_molden_path):
    """Formaldehyde: C=O double bond, two C-H bonds, core LPs on C and O and two O lone pairs. The reference values are
    from the per-hybrid Python implementation, before the hybrid optimization was batched."""
    mf = MoldenFile.load(h2co_molden_path)
    mf.coords_to_au()
    mf.to_unnormalized_primitive_coefs()
    npa_res = run_npa(mf)
    assert np.allclose(npa_res.npa_charges, [0.22980029, -0.5434533, 0.1568265, 0.1568265], atol=1e-6)

    opts = CLPOOptions(hybr_opt_conv_thresh=1e-9, hybr_opt_max_iter=1000)
    clpo = create_clpos(npa_res.sds_nao, npa_res.nao, mf.centers, opts).clpo
    q = clpo.nao_to_hybrids @ clpo.lo_to_hybrids.T
    occ = np.diag(q.T @ npa_res.sds_nao @ q)
    los = sorted((tuple(sorted(int(clpo.host_atom_of_hybrid[h]) for h in clpo.hybrids_of_lo[i])), round(occ[i], 6))
                 for i in range(len(clpo.lo_types)) if clpo.lo_types[i] in (LO_TYPE_LP, LO_TYPE_BD))
    assert los == [((0,), 1.999734), ((0, 1), 1.999838), ((0, 1), 2.0), ((0, 2), 1.995784), ((0, 3), 1.995784),
                   ((1,), 1.924390), ((1,), 1.987155), ((1,), 1.999969)]


def test_clpo_interface_orthonormal():
    """The CLPOs returned to PySCF/sunrise are orthonormal in the AO basis and come with their BD/LP graph."""
    mol = gto.M(atom="C 0 0 0; O 0 0 1.205; H 0 0.943 -0.587; H 0 -0.943 -0.587", basis="6-31g", verbose=0)
    mf = scf.RHF(mol).run()
    clpo_to_ao, graph = generate_CLPO_molecule_edges(mol, mf)
    s = mol.intor("int1e_ovlp")
    assert np.allclose(clpo_to_ao @ s @ clpo_to_ao.T, np.eye(mol.nao), atol=1e-10)
    assert sum(len(e) == 2 for e in graph) == 4 and sum(len(e) == 1 for e in graph) == 4
