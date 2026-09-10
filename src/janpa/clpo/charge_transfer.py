# janpa/clpo/charge_transfer.py
"""Charge transfer estimation between molecular fragments. Replaces CLPO/CT_Estimator.java. """
from __future__ import annotations
import numpy as np
from janpa.utils.matrix import transform_symmetric_to_new_basis
from .lpo import LODescription

class CTEstimator:
    """Approximate charge transfer analysis in CLPO basis."""
    def __init__(self, lo_desc: LODescription, sds_nao: np.ndarray, ct_threshold: float = 0.01):
        self.los = lo_desc
        self.ct_threshold = ct_threshold
        nao2lo = lo_desc.nao_to_hybrids @ lo_desc.lo_to_hybrids.T
        self.sds_los = transform_symmetric_to_new_basis(sds_nao, nao2lo.T)

    def ensure_clpos_intrafragment(self, fragment_ids):
        clpo2frag_id = np.zeros(len(self.los.hybrids_of_lo), dtype=int)
        for i in range(len(self.los.hybrids_of_lo)):
            hybrs = self.los.hybrids_of_lo[i]
            same_frag = True
            frag_prev = fragment_ids[self.los.host_atom_of_hybrid[hybrs[0]]]
            clpo2frag_id[i] = frag_prev
            for h in range(1, len(hybrs)):
                same_frag &= fragment_ids[self.los.host_atom_of_hybrid[hybrs[h]]] == frag_prev
            if not same_frag:
                print(f"ERROR: CLPO {i + 1} has hybrids belonging to different fragments!")
                clpo2frag_id[i] = -1
        return clpo2frag_id

    def perform_ct(self, frag1, frag2, clpo2frag_id):
        tot_ct = 0.0
        lowest_donor_occ = 1.0
        highest_accpt_occ = 1.0
        n = len(clpo2frag_id)
        n_ct_below_threshold = 0
        q_ct_below_threshold = 0.0
        
        print(" orb.num. description occup. --> charge, e --> occup. description orb.num.")
        for i in range(n):
            if clpo2frag_id[i] != frag1: continue
            dii = self.sds_los[i, i]
            for j in range(n):
                if j == i: continue
                if clpo2frag_id[j] != frag2: continue
                djj = self.sds_los[j, j]
                if (dii > lowest_donor_occ) and (djj < dii) and (djj < highest_accpt_occ):
                    dij = self.sds_los[i, j]
                    q_ct = dij * dij / dii
                    if q_ct > self.ct_threshold:
                        print(f" {i + 1:5d} {self.los.lo_labels[i]:>15s} {dii:7.4f} --> {q_ct:7.5f} --> {djj:7.4f} {self.los.lo_labels[j]:>15s} {j + 1:5d}")
                    else:
                        n_ct_below_threshold += 1
                        q_ct_below_threshold += q_ct
                    tot_ct += q_ct
        print()
        print(f"{n_ct_below_threshold} orbital pairs with total charge transfer of {q_ct_below_threshold:.5f} were not printed")
        print()
        return tot_ct

    def compute_ct(self, fragment_ids: np.ndarray) -> None:
        """Print charge transfer analysis between fragments."""
        print("Approximate charge transfer analysis in CLPO basis")
        print("NOTE: this is an experimental feature AND IS SUBJECT TO CHANGE!")
        print(" We expect to have the underlying theory published soon...")
        print("")
        print(f"Threshold for printing: {self.ct_threshold:.5f} ")
        print("")
        
        n_tot_fragments = np.max(fragment_ids)
        fragment_from = np.zeros(n_tot_fragments)
        fragment_to = np.zeros(n_tot_fragments)
        
        clpo2frag_id = self.ensure_clpos_intrafragment(fragment_ids)
        
        print("IntrAfragment charge transfers (conjugation analysis, etc.)")
        for i_frag in range(1, n_tot_fragments + 1):
            print(f"CT within fragment {i_frag}\n")
            self.perform_ct(i_frag, i_frag, clpo2frag_id)
            
        print()
        print("IntErfragment charge transfers (donor-acceptor analysis, etc.)")
        print()
        for i_frag in range(1, n_tot_fragments + 1):
            for j_frag in range(i_frag + 1, n_tot_fragments + 1):
                print(f"CT between from fragment {i_frag} to fragment {j_frag}")
                tot_ct = self.perform_ct(i_frag, j_frag, clpo2frag_id)
                fragment_from[i_frag - 1] += tot_ct
                fragment_to[j_frag - 1] += tot_ct
                
                print(f"CT between from fragment {j_frag} to fragment {i_frag}")
                tot_ct = self.perform_ct(j_frag, i_frag, clpo2frag_id)
                fragment_from[j_frag - 1] += tot_ct
                fragment_to[i_frag - 1] += tot_ct
                
        print("Inter-molecular charge transfer summary (NEW):")
        print("Mol Accepted -Donated = got_electrons")
        for f in range(n_tot_fragments):
            print(f"{f + 1:3d} +{fragment_to[f]:.5f} -{fragment_from[f]:.5f} = {fragment_to[f] - fragment_from[f]:+.5f}")