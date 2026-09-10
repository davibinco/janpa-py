import pytest
import numpy as np
from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa
from janpa.clpo.lpo import create_clpos, CLPOOptions
from janpa.clpo.lodesc import LO_TYPE_BD, LO_TYPE_NB

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