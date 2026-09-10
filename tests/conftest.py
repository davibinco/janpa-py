import pytest
from pyscf import gto, scf
from pyscf.tools import molden as pyscf_molden

@pytest.fixture(scope="session")
def h2_molden_path(tmp_path_factory):
    """Generate a simple H2 MOLDEN file using PySCF."""
    tmp_dir = tmp_path_factory.mktemp("data")
    molden_file = tmp_dir / "h2.molden"
    
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", unit="Angstrom")
    mf = scf.RHF(mol).run()
    
    pyscf_molden.from_mo(mol, str(molden_file), mf.mo_coeff, occ=mf.mo_occ, ene=mf.mo_energy)
    return str(molden_file)

@pytest.fixture(scope="session")
def h4_molden_path(tmp_path_factory):
    """Generate a linear H4 MOLDEN file using PySCF."""
    tmp_dir = tmp_path_factory.mktemp("data")
    molden_file = tmp_dir / "h4.molden"
    
    mol = gto.M(atom="H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", basis="sto-3g", unit="Angstrom")
    mf = scf.RHF(mol).run()
    
    pyscf_molden.from_mo(mol, str(molden_file), mf.mo_coeff, occ=mf.mo_occ, ene=mf.mo_energy)
    return str(molden_file)