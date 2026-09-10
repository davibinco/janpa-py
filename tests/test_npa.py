import pytest
import numpy as np
from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa

def test_h2_npa_charges(h2_molden_path):
    """Test that H2 has symmetric NPA charges (approx 0.0)."""
    mf = MoldenFile.load(h2_molden_path)
    mf.coords_to_au()
    mf.to_unnormalized_primitive_coefs()
    
    res = run_npa(mf)
    
    # H2 should have 2 electrons total -> Sum of NPA charges = 2 (protons) - 2 (electrons) = 0
    assert np.isclose(np.sum(res.npa_charges), 0.0, atol=1e-5)
    
    # Both H atoms should have ~0 charge due to symmetry
    assert np.allclose(res.npa_charges, 0.0, atol=1e-3)
    
    # Wiberg bond order for H-H should be close to 1.0
    assert res.wiberg_indices[0, 1] > 0.8

def test_h4_npa_electrons(h4_molden_path):
    """Test that H4 has 4 electrons and correct total charge."""
    mf = MoldenFile.load(h4_molden_path)
    mf.coords_to_au()
    mf.to_unnormalized_primitive_coefs()
    
    res = run_npa(mf)
    
    # Total NPA charge should be 0 (4 protons - 4 electrons)
    assert np.isclose(np.sum(res.npa_charges), 0.0, atol=1e-5)