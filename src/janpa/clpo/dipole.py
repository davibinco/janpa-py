# janpa/clpo/dipole.py
"""Dipole moment analysis in localized orbital bases. Replaces CLPO/DipoleAn.java. """
from __future__ import annotations
import numpy as np
from janpa.utils.matrix import transform_symmetric_to_new_basis

class DipoleAnalysis:
    """Compute and print dipole contributions in NAO/LHO/CLPO bases."""
    def __init__(self, npa_result, clpo_result):
        self.npa = npa_result
        self.lpo_clpo = clpo_result
        self.ao_dipole_matr = [None] * 3
        self.nao_dipole_matr = [None] * 3
        self.lho_dipole_matr = [None] * 3
        self.clpo_dipole_matr = [None] * 3
        self.lpo_dipole_matr = [None] * 3
        self.true_dipole_au = np.zeros(3)
        
        self.bohr_radius = 0.529177210903 

    def print_dipolar_contribs(self, dipoles, sds, lbl, labels, print_occ_threshold):
        d_printed_sum = np.zeros(3)
        d_all_sum = np.zeros(3)
        print(" [ DipX/au DipY/au DipZ/au ] x Occupancy from Contributing Orbital -> [ Charge Center (Angstroms) ]")
        for b in range(len(labels)):
            occ = sds[b, b]
            if occ > print_occ_threshold:
                print(f" += [{dipoles[0][b, b]:10.5f} {dipoles[1][b, b]:10.5f} {dipoles[2][b, b]:10.5f}] x ({occ:7.4f}) from {labels[b]:>22s} ({lbl} {b + 1}) -> [{-dipoles[0][b, b]*self.bohr_radius:10.5f} {-dipoles[1][b, b]*self.bohr_radius:10.5f} {-dipoles[2][b, b]*self.bohr_radius:10.5f}]")
                dip_abs = 0.0
                for mu in range(3):
                    d_printed_sum[mu] += dipoles[mu][b, b] * occ
                    dip_abs += (dipoles[mu][b, b] * occ) ** 2
        for b in range(len(labels)):
            occ = sds[b, b]
            for mu in range(3):
                d_all_sum[mu] += dipoles[mu][b, b] * occ
                
        print(" -------------------------------------------------")
        print(f"{'Localized total:':>25s} [{d_all_sum[0]:10.5f} {d_all_sum[1]:10.5f} {d_all_sum[2]:10.5f}] a.u.")
        print(f"{'Exact:':>25s} [{self.true_dipole_au[0]:10.5f} {self.true_dipole_au[1]:10.5f} {self.true_dipole_au[2]:10.5f}] a.u.")
        print(f"{'Diff = Localized-Exact:':>25s} [{d_all_sum[0] - self.true_dipole_au[0]:10.5f} {d_all_sum[1] - self.true_dipole_au[1]:10.5f} {d_all_sum[2] - self.true_dipole_au[2]:10.5f}] a.u.")
        print(f"{'Locl. in printed:':>25s} [{d_printed_sum[0]:10.5f} {d_printed_sum[1]:10.5f} {d_printed_sum[2]:10.5f}] a.u.")
        print(f"{'Diff = Printed-Exact:':>25s} [{d_printed_sum[0] - self.true_dipole_au[0]:10.5f} {d_printed_sum[1] - self.true_dipole_au[1]:10.5f} {d_printed_sum[2] - self.true_dipole_au[2]:10.5f}] a.u.")

    def lho_dip(self, sds_orb, el_dip_matrices, occ_print_thresh, orb_type):
        print(f"\n{orb_type} dipole moments (neutralized by the equal fraction of nuclear charge):")
        print(f"{'Id':>5s}\t{'Atom':<7s}\t{'occup.':>9s} * [{'dx ':>10s} {'dy ':>10s} {'dz ':>10s}]")
        d_printed_sum = np.zeros(3)
        d_all_sum = np.zeros(3)
        los = self.lpo_clpo.lpo if orb_type == "AHO" else self.lpo_clpo.clpo
        
        for i in range(len(los.host_atom_of_hybrid)):
            a = los.host_atom_of_hybrid[i]
            occ = sds_orb[i, i]
            for mu in range(3):
                tmp = el_dip_matrices[mu][i, i] + self.npa.centers[a].r0[mu] * 1.0
                d_all_sum[mu] += tmp * occ
            if occ > occ_print_thresh:
                dip_abs = 0.0
                print(f"{i + 1:5d}\t{self.npa.centers[a].name + str(a + 1):<7s}\t{occ:9.5f} * [", end="")
                for mu in range(3):
                    tmp = el_dip_matrices[mu][i, i] + self.npa.centers[a].r0[mu] * 1.0
                    if mu != 2:
                        print(f"{tmp:10.5f} ", end="")
                    else:
                        print(f"{tmp:10.5f}", end="")
                    d_printed_sum[mu] += tmp * occ
                    dip_abs += (tmp * occ) ** 2
                dip_abs = np.sqrt(dip_abs)
                print(f"], abs. = {dip_abs:10.5f} a.u.")
                
        print(f"{'Total in printed:':>25s} [{d_printed_sum[0]:10.5f} {d_printed_sum[1]:10.5f} {d_printed_sum[2]:10.5f}] a.u.")
        print(f"{'Total in localized:':>25s} [{d_all_sum[0]:10.5f} {d_all_sum[1]:10.5f} {d_all_sum[2]:10.5f}] a.u.")
        print(f"{'Exact:':>25s} [{self.true_dipole_au[0]:10.5f} {self.true_dipole_au[1]:10.5f} {self.true_dipole_au[2]:10.5f}] a.u.")
        print(f"Note: {occ_print_thresh:.5f} occupancy threshold was used for printing\n")

    def bond_dipoles(self, sds_hybrids, sds_los, los, hybr_dipoles, orb_type):
        tot_dip = np.zeros(3)
        for i in range(len(los.host_atom_of_hybrid)):
            hybrs = los.hybrids_of_lo[i]
            if len(hybrs) == 1:
                occ = sds_los[i, i]
                if occ > 1.0:
                    a = los.host_atom_of_hybrid[hybrs[0]]
                    print(f"{i + 1:5d} {los.lo_labels[i]:>25s}: ", end="")
                    print(f"{'':>10s} occ*<h|d|h> = [", end="")
                    for mu in range(3):
                        tmp = occ * (hybr_dipoles[mu][hybrs[0], hybrs[0]] + self.npa.centers[a].r0[mu])
                        print(f"{tmp:10.5f} ", end="")
                        tot_dip[mu] += tmp
                    print("]")
            if len(hybrs) == 2:
                occ = sds_los[i, i]
                ca = los.lo_to_hybrids[i, hybrs[0]]
                cb = los.lo_to_hybrids[i, hybrs[1]]
                a = los.host_atom_of_hybrid[hybrs[0]]
                b = los.host_atom_of_hybrid[hybrs[1]]
                print(f"{i + 1:5d} {los.lo_labels[i]:>25s}: ", end="")
                print(f"{'':>10s} occ*|c1|^2*<h1|d|h1> = [", end="")
                for mu in range(3):
                    tmp = occ * ca * ca * (hybr_dipoles[mu][hybrs[0], hybrs[0]] + self.npa.centers[a].r0[mu])
                    print(f"{tmp:10.5f} ", end="")
                    tot_dip[mu] += tmp
                print("] ", end="")
                print(f"{'':>10s} occ*|c2|^2*<h2|d|h2> = [", end="")
                for mu in range(3):
                    tmp = occ * cb * cb * (hybr_dipoles[mu][hybrs[1], hybrs[1]] + self.npa.centers[b].r0[mu])
                    print(f"{tmp:10.5f} ", end="")
                    tot_dip[mu] += tmp
                print("] ", end="")
                print(f"{'':>10s} occ*c1*c2*2*<h1|d|h2> = [", end="")
                for mu in range(3):
                    tmp = occ * 2 * ca * cb * hybr_dipoles[mu][hybrs[0], hybrs[1]]
                    print(f"{tmp:10.5f} ", end="")
                    tot_dip[mu] += tmp
                print("] ")
                
        print("total dipole without atomic chagres:")
        for mu in range(3):
            print(f"{tot_dip[mu]:10.5f}", end="")
        print()
        print("total dipole with atomic chagres:")
        for mu in range(3):
            for a in range(len(self.npa.centers)):
                tot_dip[mu] += self.npa.centers[a].r0[mu] * self.npa.npa_charges[a]
                if mu == 0:
                    print(f"{a + 1}: + {self.npa.centers[a].r0[mu] * self.npa.npa_charges[a]:10.5f} ")
            print(f"{tot_dip[mu]:10.5f}", end="")
        print()

    def print_dipoles(self, occ_print_threshold: float = 1e-3) -> None:
        print()
        print("Localized analysis of the dipole moment")
        print()
        print("Forming the dipole matrix elements in AO and NAO bases...")
        
        lho_2_ao = self.lpo_clpo.clpo.nao_to_hybrids.T @ self.npa.nao_2_ao
        clpo_2_ao = self.lpo_clpo.clpo.lo_to_hybrids @ lho_2_ao
        lpo_2_ao = self.lpo_clpo.lpo.lo_to_hybrids @ self.lpo_clpo.lpo.nao_to_hybrids.T @ self.npa.nao_2_ao
        
        for mu in range(3):
            nucl_dip = 0.0
            for c in range(len(self.npa.centers)):
                nucl_dip += self.npa.centers[c].r0[mu] * self.npa.centers[c].z
                
            sz = len(self.npa.basis)
            self.ao_dipole_matr[mu] = np.zeros((sz, sz))
            for i in range(sz):
                for j in range(sz):
                    tmp = - self.npa.bs_integrals.dipole_matrix[i][j][mu]
                    self.ao_dipole_matr[mu][i, j] = tmp
                    
            tot_el_dip = np.trace(self.npa.d_global @ self.ao_dipole_matr[mu])
            print(f"Dipole [ {'XYZ'[mu]} ]: Total = {tot_el_dip + nucl_dip:7.4f} a.u. ({(tot_el_dip + nucl_dip)/0.393430:7.4f} Debye), nuclear = {nucl_dip:8.5f}, electronic = {tot_el_dip:8.5f}, ")
            self.true_dipole_au[mu] = tot_el_dip + nucl_dip
            
            self.nao_dipole_matr[mu] = transform_symmetric_to_new_basis(self.ao_dipole_matr[mu], self.npa.nao_2_ao)
            self.lho_dipole_matr[mu] = transform_symmetric_to_new_basis(self.ao_dipole_matr[mu], lho_2_ao)
            self.clpo_dipole_matr[mu] = transform_symmetric_to_new_basis(self.ao_dipole_matr[mu], clpo_2_ao)
            self.lpo_dipole_matr[mu] = transform_symmetric_to_new_basis(self.ao_dipole_matr[mu], lpo_2_ao)
            
        sds_clpo = transform_symmetric_to_new_basis(self.npa.sds, clpo_2_ao)
        sds_lpo = transform_symmetric_to_new_basis(self.npa.sds, lpo_2_ao)
        sds_lho = transform_symmetric_to_new_basis(self.npa.sds, lho_2_ao)
        
        self.lho_dip(sds_lho, self.lho_dipole_matr, occ_print_threshold, "LHO")
        self.bond_dipoles(sds_lho, sds_clpo, self.lpo_clpo.clpo, self.lho_dipole_matr, "")
        
        print()
        print("Dipole analysis in CLPO basis")
        self.print_dipolar_contribs(self.clpo_dipole_matr, sds_clpo, "CLPO", self.lpo_clpo.clpo.lo_labels, 1.0)