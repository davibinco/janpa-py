import pytest
from janpa.io.molden import MoldenFile

def test_load_molden_h2(h2_molden_path):
    mf = MoldenFile.load(h2_molden_path)
    assert mf.n_atoms == 2
    assert mf.n_basis == 2
    assert len(mf.mos) == 2
    assert mf.is_spherical is True
    # Check that occupations were read correctly (2.0 for occupied, 0.0 for virtual)
    assert mf.mos[0].occupancy == 2.0
    assert mf.mos[1].occupancy == 0.0