import contextlib
import os
import pytest
import numpy as np
from pyscf import gto, scf
from pyscf.tools import molden as pyscf_molden

psi4 = pytest.importorskip("psi4")

from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa
from janpa.janpa_interface import generate_CLPO_molecule_edges as pyscf_CLPO_edges
from janpa.janpa_psi4_interface import generate_CLPO_molecule_edges, generate_HAO_molecule

H2CO = "C 0 0 0; O 0 0 1.205; H 0 0.943 -0.587; H 0 -0.943 -0.587"


def _psi4_scf(basis, **options):
    psi4.core.clean()
    psi4.core.be_quiet()
    psi4.set_options({"basis": basis, "scf_type": "pk", "e_convergence": 1e-10, "d_convergence": 1e-10, **options})
    mol = psi4.geometry("\n".join(H2CO.split("; ")) + "\nunits angstrom")
    _, wfn = psi4.energy("scf", molecule=mol, return_wfn=True)
    return mol, wfn


def _npa_charges(molden_path, normalize):
    mf = MoldenFile.load(molden_path)
    mf.coords_to_au()
    mf.to_unnormalized_primitive_coefs()
    if normalize:
        mf.normalize_contractions()
    with contextlib.redirect_stdout(open(os.devnull, "w")):
        return run_npa(mf).npa_charges


@pytest.mark.parametrize("basis", ["6-31g", "cc-pvdz"])
def test_psi4_npa_charges_match_pyscf(basis, tmp_path):
    """The psi4 MOLDEN file (no [5D] tags for s/p bases, unnormalized general contractions for cc-pVDZ) gives the same
    NPA charges as the PySCF one."""
    _, wfn = _psi4_scf(basis, puream=True)
    psi4.molden(wfn, str(tmp_path / "psi4.molden"))
    mol = gto.M(atom=H2CO, basis=basis, verbose=0)
    mf = scf.RHF(mol).run(conv_tol=1e-10)
    pyscf_molden.from_mo(mol, str(tmp_path / "pyscf.molden"), mf.mo_coeff, occ=mf.mo_occ, ene=mf.mo_energy)
    assert np.allclose(_npa_charges(tmp_path / "psi4.molden", True), _npa_charges(tmp_path / "pyscf.molden", False),
                       atol=1e-5)


@pytest.mark.parametrize("basis, puream", [("6-31g", False), ("cc-pvdz", True)])
def test_psi4_clpo_orthonormal(basis, puream):
    """The CLPOs and AHOs are orthonormal in the psi4 AO basis, also for the Cartesian p shells of 6-31G, and form the
    C=O double bond, two C-H bonds, core LPs on C and O and two O lone pairs."""
    mol, wfn = _psi4_scf(basis, puream=puream)
    s = np.asarray(psi4.core.MintsHelper(wfn.basisset()).ao_overlap())
    nbf = wfn.basisset().nbf()
    clpo_to_ao, graph = generate_CLPO_molecule_edges(mol, wfn)
    assert np.allclose(clpo_to_ao @ s @ clpo_to_ao.T, np.eye(nbf), atol=1e-8)
    assert sum(len(e) == 2 for e in graph) == 4 and sum(len(e) == 1 for e in graph) == 4
    aho_to_ao = generate_HAO_molecule(mol, wfn)
    assert np.allclose(aho_to_ao @ s @ aho_to_ao.T, np.eye(nbf), atol=1e-8)


def test_psi4_clpo_matches_pyscf():
    """Same 6-31G CLPOs as from PySCF. The two programs order the shells of an atom differently, so the orbitals are
    compared by their sorted absolute AO coefs."""
    mol, wfn = _psi4_scf("6-31g")
    clpo_psi4, graph_psi4 = generate_CLPO_molecule_edges(mol, wfn)
    pmol = gto.M(atom=H2CO, basis="6-31g", verbose=0)
    clpo_pyscf, graph_pyscf = pyscf_CLPO_edges(pmol, scf.RHF(pmol).run(conv_tol=1e-10))
    assert graph_psi4 == graph_pyscf
    assert np.allclose(np.sort(np.abs(clpo_psi4)), np.sort(np.abs(clpo_pyscf)), atol=1e-4)


def test_psi4_runs_scf_and_imposes_edges():
    """Without a Wavefunction the SCF is run with the current psi4 options; imposed bonds are kept."""
    mol, _ = _psi4_scf("cc-pvdz", puream=True)
    _, graph = generate_CLPO_molecule_edges(mol)
    _, imposed = generate_CLPO_molecule_edges(mol, edges=[(0, 1), (0, 1)])
    assert graph == imposed


def test_psi4_cartesian_d_rejected():
    mol, wfn = _psi4_scf("6-31g*", puream=False)
    with pytest.raises(ValueError, match="puream"):
        generate_CLPO_molecule_edges(mol, wfn)
