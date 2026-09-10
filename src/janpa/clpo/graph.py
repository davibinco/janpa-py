import numpy as np
from warnings import warn
from janpa.clpo.lodesc import LO_TYPE_RY, LO_TYPE_LP, LO_TYPE_BD, LO_TYPE_NB

def extract_clpo_graph_from_janpa(npa_res, clpo_res):
    """
    Extracts the CLPO graph (edges) directly from janpa-py results,
    mimicking the exact behavior of sunrise's extract_clpo_graph.
    
    Parameters:
    -----------
    npa_res : NPAResult
        The result object from run_npa().
    clpo_res : CLPOResult
        The result object from create_clpos().
        
    Returns:
    --------
    list of tuples:
        - (i,) for lone pairs (occupancy >= 1.0)
        - (i, i+1) for BD/NB pairs
        where i is the 0-based index of the CLPO orbital.
    """
    # 1. Get the transformation matrix from NAO to CLPO
    # nao_to_hybrids: (n_nao, n_hybrids)
    # lo_to_hybrids: (n_lo, n_hybrids)
    # The columns of Q are the CLPOs expressed in the NAO basis.
    Q = clpo_res.clpo.nao_to_hybrids @ clpo_res.clpo.lo_to_hybrids.T
    
    # 2. Compute the SDS matrix in the CLPO basis
    # sds_nao is (n_nao, n_nao)
    sds_clpo = Q.T @ npa_res.sds_nao @ Q
    
    # 3. Extract occupancies (diagonal elements)
    occupancies = np.diag(sds_clpo)
    lo_types = clpo_res.clpo.lo_types
    
    nodes = []
    n_lo = len(occupancies)
    
    i = 0
    while i < n_lo:
        lo_type = lo_types[i]
        occ = occupancies[i]
        
        if lo_type == LO_TYPE_LP:
            # Replicate sunrise warning for radical-like lone pairs
            if 0.5 < occ < 1.5:
                warn(f'Lone Pair {i} found with occupation close to 1 => {occ:.5f}, take care.')
            nodes.append((i,))
            i += 1
            
        elif lo_type == LO_TYPE_RY:
            # Skip Rydberg orbitals (they appear as "(LP)" with occ < 1.0 in the sunrise graph file)
            i += 1
            
        elif lo_type == LO_TYPE_BD:
            # The next orbital MUST be the corresponding NB
            if i + 1 >= n_lo or lo_types[i+1] != LO_TYPE_NB:
                raise ValueError(f"BD orbital {i} without following NB orbital")
                
            bd_occ = occ
            nb_occ = occupancies[i+1]
            
            # Replicate sunrise warning for depleted bonds
            if bd_occ + nb_occ < 1.7:
                warn(f'Bond pair orbitals [{i},{i+1}] population expected under expected 2e-, predicted: {bd_occ + nb_occ:.5f}, take care with predicted edges.')
                
            nodes.append((i, i + 1))
            i += 2
            
        else:
            # Fallback for unexpected types
            i += 1

    return nodes