# janpa/npa/npa.py
"""Natural Population Analysis (NPA) — the 7-step algorithm. Replaces onpa/NPA.java."""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from janpa.gto.basis import AtomicCenter, BasisFunction
from janpa.gto.integrals import IntegralEngine
from janpa.io.molden import MoldenFile
from janpa.utils.matrix import matrix_sqrt, generalized_eigh, transform_symmetric_to_new_basis
from janpa.utils.printout import print_matrix, STARS
from janpa.npa.options import NPAOptions


# Natural Minimal Basis size per angular momentum (s, p, d, f, g) per atomic number (Z-1)
NMB_PER_ATOM = [
    [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15, 15, 15, 15],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 15, 15, 10, 10, 10, 10, 10, 15, 10, 10, 10, 10, 10, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]


def create_nao_labels(basis: list[BasisFunction]) -> list[str]:
    """Creates standard labels for NAO/PNAO/AO basis functions."""
    result = []
    for bf in basis:
        fn_id = f"R{bf.radial_part_id + 1}*{'spdfg'[bf.l]}({bf.m})"
        if bf.additional_r_power != 0:
            fn_id += f"*r^{bf.additional_r_power}"
        if not bf.is_nrb:
            res = f"A{bf.center_id}: {fn_id}"
        else:
            res = f"A{bf.center_id}*: {fn_id}"
        result.append(res)
    return result


def intracenter_basis_orthogonalization(
    initial_basis: list[BasisFunction],
    n_atoms: int,
    sds: np.ndarray,
    overlap_matrix: np.ndarray,
    verbose_print: bool = False,
    glb_print: bool = False,
) -> list[BasisFunction]:
    """Diagonalizes SDS blocks per atom and per L to produce natural atomic orbitals."""
    result = [None] * len(initial_basis)
    input_basis_not_empty = len(initial_basis) != 0
    
    l_max_of_atom = np.zeros(n_atoms, dtype=int)
    for bf in initial_basis:
        if bf.l > l_max_of_atom[bf.center_id - 1]:
            l_max_of_atom[bf.center_id - 1] = bf.l
            
    basis_functions_per_atom = np.zeros(n_atoms, dtype=int)
    for bf in initial_basis:
        basis_functions_per_atom[bf.center_id - 1] += 1
        
    sorted_basis_functions = np.zeros(len(initial_basis), dtype=int)
    first_sorted_basis_function_of_atom = np.zeros(n_atoms, dtype=int)
    for i in range(1, n_atoms):
        first_sorted_basis_function_of_atom[i] = first_sorted_basis_function_of_atom[i-1] + basis_functions_per_atom[i-1]
        
    sliding_indexes = first_sorted_basis_function_of_atom.copy()
    for i, bf in enumerate(initial_basis):
        sorted_basis_functions[sliding_indexes[bf.center_id - 1]] = i
        sliding_indexes[bf.center_id - 1] += 1
        
    n_pnaos = 0
    
    if verbose_print:
        print("Producing natural orbitals for each center...")
        
    for cntr in range(n_atoms):
        if not input_basis_not_empty:
            break
        for L in range(l_max_of_atom[cntr] + 1):
            if verbose_print:
                print(f"center = {cntr+1}, L = {L}")
                
            radial_parts_used_for_this_l = []
            basis_functions_used_for_this_l = []
            
            for i in range(first_sorted_basis_function_of_atom[cntr], first_sorted_basis_function_of_atom[cntr] + basis_functions_per_atom[cntr]):
                bf_idx = sorted_basis_functions[i]
                bf = initial_basis[bf_idx]
                if bf.l == L:
                    rid = bf.radial_part_id
                    if rid not in radial_parts_used_for_this_l:
                        radial_parts_used_for_this_l.append(rid)
                    if bf_idx not in basis_functions_used_for_this_l:
                        basis_functions_used_for_this_l.append(bf_idx)
                        
            n_radial_parts = len(radial_parts_used_for_this_l)
            n_basis_functions = len(basis_functions_used_for_this_l)
            
            if n_radial_parts == 0:
                continue
                
            local_sds = np.zeros((n_radial_parts, n_radial_parts))
            local_s = np.zeros((n_radial_parts, n_radial_parts))
            
            for r_idx1 in range(n_radial_parts):
                for r_idx2 in range(r_idx1, n_radial_parts):
                    n_fns = 0
                    for i in range(n_basis_functions):
                        for j in range(i, n_basis_functions):
                            bf_i = initial_basis[basis_functions_used_for_this_l[i]]
                            bf_j = initial_basis[basis_functions_used_for_this_l[j]]
                            
                            same_m = bf_i.m == bf_j.m
                            ith_eq_r1 = bf_i.radial_part_id == radial_parts_used_for_this_l[r_idx1]
                            ith_eq_r2 = bf_i.radial_part_id == radial_parts_used_for_this_l[r_idx2]
                            jth_eq_r1 = bf_j.radial_part_id == radial_parts_used_for_this_l[r_idx1]
                            jth_eq_r2 = bf_j.radial_part_id == radial_parts_used_for_this_l[r_idx2]
                            
                            if same_m and ((ith_eq_r1 and jth_eq_r2) or (ith_eq_r2 and jth_eq_r1)):
                                n_fns += 1
                                local_s[r_idx1, r_idx2] += overlap_matrix[basis_functions_used_for_this_l[i], basis_functions_used_for_this_l[j]]
                                local_sds[r_idx1, r_idx2] += sds[basis_functions_used_for_this_l[i], basis_functions_used_for_this_l[j]]
                                
                    if n_fns > 0:
                        local_s[r_idx1, r_idx2] /= n_fns
                        local_sds[r_idx1, r_idx2] /= n_fns
                        
                    local_s[r_idx2, r_idx1] = local_s[r_idx1, r_idx2]
                    local_sds[r_idx2, r_idx1] = local_sds[r_idx1, r_idx2]
                    
            if glb_print:
                print("Local_SDS")
                print(local_sds)
                print("Local_S")
                print(local_s)
                
            eigvals, eigvecs = generalized_eigh(local_sds, local_s)
            
            for source_bf in range(n_basis_functions):
                bf_orig = initial_basis[basis_functions_used_for_this_l[source_bf]]
                new_bf = BasisFunction(
                    l=L,
                    m=bf_orig.m,
                    center_id=cntr + 1,
                    radial_part_id=bf_orig.radial_part_id,
                    coefs=np.zeros(len(initial_basis)),
                    is_nrb=False
                )
                
                r_idx1 = radial_parts_used_for_this_l.index(bf_orig.radial_part_id)
                new_bf.weight = eigvals[r_idx1]
                
                for j in range(n_basis_functions):
                    bf_j = initial_basis[basis_functions_used_for_this_l[j]]
                    if bf_j.m == new_bf.m:
                        r_idx2 = radial_parts_used_for_this_l.index(bf_j.radial_part_id)
                        new_bf.coefs[basis_functions_used_for_this_l[j]] = eigvecs[r_idx2, r_idx1]
                        
                if glb_print:
                    tmp = 0.0
                    coefs_str = ""
                    for i in range(len(initial_basis)):
                        coefs_str += f"{new_bf.coefs[i]:12.5f}"
                        for j in range(len(initial_basis)):
                            tmp += new_bf.coefs[i] * overlap_matrix[i, j] * new_bf.coefs[j]
                    print(coefs_str)
                    print(f"pnao {n_pnaos}: norm2 = {tmp:.7f}")
                    
                result[n_pnaos] = new_bf
                n_pnaos += 1
                
    print(f"Total number of natural functions produced: {n_pnaos}")
    if verbose_print:
        print("weights of the natural functions produced:")
        occ_total = 0.0
        for i in range(n_pnaos):
            print(f"{result[i].weight:15.7f}", end="")
            occ_total += result[i].weight
        print()
        print(f"sum = {occ_total:12.7f}\n")
        
    return result[:n_pnaos]


@dataclass
class NPAResult:
    """Container for all NPA results."""
    nao: list[BasisFunction] = field(default_factory=list)
    npa_charges: np.ndarray = field(default_factory=lambda: np.array([]))
    overlap_matrix: np.ndarray | None = None
    dipole_matrix: np.ndarray | None = None
    density_matrix: np.ndarray | None = None
    sds: np.ndarray | None = None
    sds_nao: np.ndarray | None = None
    nao_to_ao: np.ndarray | None = None
    overlap_nao: np.ndarray | None = None
    wiberg_indices: np.ndarray | None = None


class NPA:
    def __init__(self, options: NPAOptions):
        self.options = options
        self.basis: list[BasisFunction] = []
        self.centers: list[AtomicCenter] = []
        self._data_loaded = False
        
        self.overlap_matrix: np.ndarray | None = None
        self.d_global: np.ndarray | None = None
        self.s05ds05: np.ndarray | None = None
        self.sds: np.ndarray | None = None
        self.pnao_overlap_matrix: np.ndarray | None = None
        self.sds_pnao: np.ndarray | None = None
        self.nao_to_ao: np.ndarray | None = None
        self.overlap_nao: np.ndarray | None = None
        self.sds_nao: np.ndarray | None = None
        
        self.nao: list[BasisFunction] = []
        self.pnaos: list[BasisFunction] = []
        self.npa_charges: np.ndarray | None = None
        self.ao_names: list[str] = []
        self.pnao_labels: list[str] = []
        
        self.verbose_print = options.VerbosePrint.get_boolean()
        
    def load_mo_from_molden(self, molden: MoldenFile) -> bool:
        if molden is None:
            return False
            
        self._data_loaded = False
        self.basis = molden.basis
        self.centers = molden.centers
        
        n_bas = len(self.basis)
        print("Building overlap S and dipole matrices...")
        
        engine = IntegralEngine(self.basis, molden.radial_parts, self.centers)
        overlap_matrix, dipole_matrix = engine.build_overlap_and_dipole()
        self.overlap_matrix = overlap_matrix
        
        print(f"Overlap & dipole integrals evaluated.")
        
        bs_norms_ok = False
        mo_norms_ok = False
        mo_orth_ok = False
        
        print("Checking whether the basis functions are unity-normalized...")
        bs_norms = np.diag(overlap_matrix)
        max_dev = 0.0
        max_dev_index = -1
        
        for basf in range(n_bas):
            dev = abs(bs_norms[basf] - 1.0)
            if dev > max_dev:
                max_dev = dev
                max_dev_index = basf
                
        print(f"Maximum deviation of the basis function norm2 from unity: BFN {max_dev_index} (0-based num.), max|norm2-1| = {max_dev:.3E}")
        
        if max_dev > self.options.bf_nrm2_dev_threshold and self.options.do_BS_renormalize:
            print("Trying to make basis functions unity-normalized...")
            for basf1 in range(n_bas):
                for basf2 in range(n_bas):
                    overlap_matrix[basf1, basf2] /= bs_norms[basf1] * bs_norms[basf2]
                    
        if max_dev < self.options.bf_nrm2_dev_threshold or self.options.do_BS_renormalize:
            bs_norms_ok = True
            
        print("Checking for the eigenvalues (linear (in)dependency) of the basis function overlap matrix...")
        eigvals = np.linalg.eigvalsh(overlap_matrix)
        s_min_abs = np.min(np.abs(eigvals))
        print(f"The smallest eigenvalue of the basis function overlap matrix: {s_min_abs:7.3E}")
        
        if s_min_abs < self.options.S_AO_lin_dep_thresh:
            print("WARNING: the basis set functions seems to be linearly dependent!")
            if self.options.fail_if_lin_dep:
                return False
                
        overlap_matrix_sqrt = matrix_sqrt(overlap_matrix)
        
        print("Checking whether the orbitals are unity-normalized...")
        max_dev = 0.0
        max_dev_index = -1
        n_mo = len(molden.mos)
        
        mo_coef_matrix = np.zeros((n_mo, n_bas))
        for mo_idx in range(n_mo):
            mo_coef_matrix[mo_idx, :] = molden.mos[mo_idx].bs_coefs
            
        if self.options.do_MO_renormalize:
            mo_norm2s = np.zeros(n_mo)
            for mo_idx in range(n_mo):
                mo_norm2s[mo_idx] = mo_coef_matrix[mo_idx] @ overlap_matrix @ mo_coef_matrix[mo_idx]
                dev = abs(mo_norm2s[mo_idx] - 1.0)
                if dev > max_dev:
                    max_dev = dev
                    max_dev_index = mo_idx
                    
            print(f"Maximum deviation of the orbital norm2 from unity: MO {max_dev_index} (0-based num.), max|norm2-1| = {max_dev:.3E}")
            
            if max_dev > self.options.mo_non1_renorm_thresh:
                print("Trying to make MOs unity-normalized...")
                for mo_idx in range(n_mo):
                    molden.mos[mo_idx].bs_coefs /= np.sqrt(mo_norm2s[mo_idx])
                    mo_coef_matrix[mo_idx, :] = molden.mos[mo_idx].bs_coefs
                    
        print("Checking MO overlap matrix...")
        mo_overlap = transform_symmetric_to_new_basis(overlap_matrix, mo_coef_matrix.T)
        
        off_diag_i = 0
        off_diag_j = 1
        non1_mo = 0
        for i in range(n_mo):
            if abs(mo_overlap[i, i] - 1.0) > abs(mo_overlap[non1_mo, non1_mo] - 1.0):
                non1_mo = i
            for j in range(i + 1, n_mo):
                if abs(mo_overlap[i, j]) > abs(mo_overlap[off_diag_i, off_diag_j]) and molden.mos[i].spin == molden.mos[j].spin:
                    off_diag_i = i
                    off_diag_j = j
                    
        mo_nrm_dev = abs(mo_overlap[non1_mo, non1_mo] - 1.0)
        print(f"Maximum of MO |norm2-1|: {mo_nrm_dev:7.3E} (MO {non1_mo + 1:4d})")
        
        if mo_nrm_dev < self.options.mo_nrm2_dev_threshold:
            mo_norms_ok = True
        else:
            print("WARNING: molecular orbital normalization problem! This may spoil MO occupancies...")
            
        if n_mo > 1:
            mo_ov = abs(mo_overlap[off_diag_i, off_diag_j])
            print(f"Maximum absolute value of off-diagonal MO overlap element: {mo_ov:7.3E} (< MO {off_diag_i + 1:4d} | MO {off_diag_j + 1:4d} >)")
            if mo_ov < self.options.mo_off_diag_threshold:
                mo_orth_ok = True
            else:
                print("WARNING: molecular orbitals seem to be non-orthogonal ?!")
        else:
            mo_orth_ok = True
            
        if not bs_norms_ok or not mo_norms_ok or not mo_orth_ok:
            print("WARNING: input data seems to be improper!")
        else:
            print("First-order reduced density matrix is OK.")
            
        print()
        n_alpha = sum(1 for mo in molden.mos if mo.spin > 0)
        if n_alpha != 0 and n_alpha != n_mo and n_alpha != n_mo // 2:
            print("WARNING: the system seems to be opened-shell!\nNote that this case is current not supported so that the obtained results will be meaningless!")
            
        self.ao_names = create_nao_labels(self.basis)
        print_matrix(overlap_matrix, "The Overlap Matrix in spherical function AO basis:", self.ao_names, self.ao_names, self.options.S_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        print("Building density matrix D... ", end="")
        self.d_global = np.zeros((n_bas, n_bas))
        for basf1 in range(n_bas):
            for basf2 in range(basf1, n_bas):
                tmp = 0.0
                for mo in molden.mos:
                    tmp += mo.bs_coefs[basf1] * mo.bs_coefs[basf2] * mo.occupancy
                self.d_global[basf1, basf2] = tmp
                self.d_global[basf2, basf1] = tmp
        print("done.")
        
        print_matrix(self.d_global, "The D density Matrix in spherical function AO basis:", self.ao_names, self.ao_names, self.options.D_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        print("Building D.S... ", end="")
        ds = self.d_global @ overlap_matrix
        print("done.")
        
        tot_el_num = np.trace(ds)
        print(f"Total number of electrons: {tot_el_num:f}")
        
        sum_nucl_charge = sum(c.z for c in self.centers)
        print(f"Sum of electrons charges and the nuclei charges: {sum_nucl_charge - tot_el_num:7.5f}")
        print()
        
        print("Performing Mulliken and Lowdin population analyses...")
        a_lowdin = overlap_matrix_sqrt @ self.d_global @ overlap_matrix_sqrt
        self.s05ds05 = a_lowdin
        
        mulliken_population = np.zeros(len(self.centers))
        lowdin_population = np.zeros(len(self.centers))
        
        for i in range(n_bas):
            lowdin_population[self.basis[i].center_id - 1] += a_lowdin[i, i]
            mulliken_population[self.basis[i].center_id - 1] += ds[i, i]
            
        print("Atom\tMulliken \tLowdin \tMulliken\tLowdin")
        print("\tPopulation\tPopulation\tCharge \tCharge")
        for c_idx, c in enumerate(self.centers):
            print(f"{c.name}{c_idx + 1:7d}\t{mulliken_population[c_idx]:8.5f}\t{lowdin_population[c_idx]:8.5f}\t{c.z - mulliken_population[c_idx]:8.5f}\t{c.z - lowdin_population[c_idx]:8.5f}")
            
        print()
        print("Building S.D.S... ", end="")
        self.sds = overlap_matrix.T @ self.d_global @ overlap_matrix
        print("done.")
        
        print_matrix(self.sds, "S.D.S matrix in spherical function AO basis:", self.ao_names, self.ao_names, self.options.SDS_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        print()
        
        self._data_loaded = True
        return True
        
    def create_naos(self) -> bool:
        print(f"{STARS}\n\nCreating NAOs\n")
        
        print("\nSTEP 1. Produce PNAOs\n")
        self.pnaos = intracenter_basis_orthogonalization(self.basis, len(self.centers), self.sds, self.overlap_matrix, self.verbose_print, self.options.glbPrint)
        
        print("Sorting PNAOs...")
        sorted_pnaos = np.arange(len(self.pnaos))
        found = True
        while found:
            found = False
            for pnao in range(len(self.pnaos) - 1):
                if self.pnaos[sorted_pnaos[pnao]].weight < self.pnaos[sorted_pnaos[pnao + 1]].weight:
                    found = True
                    sorted_pnaos[pnao], sorted_pnaos[pnao + 1] = sorted_pnaos[pnao + 1], sorted_pnaos[pnao]
                    
        print("\nSTEP 2. Split PNAOs into NMB / NRB sets\n")
        print("Number of basis functions in the Natural Minimal Basis (NMB) set for each center: ")
        lmx = 5
        functions_used = np.zeros(lmx, dtype=int)
        
        for cntr_idx, center in enumerate(self.centers):
            z = int(center.z)
            if z == 0:
                print(f"Warning: center {cntr_idx + 1} has zero nuclear charge; no NMB functions will be assigned to it!")
                
            functions_used[:] = 0
            for i in range(len(self.pnaos)):
                pnao = self.pnaos[sorted_pnaos[i]]
                if pnao.center_id == cntr_idx + 1:
                    pnao.is_nrb = False
                    functions_used[pnao.l] += 1
                    nmb_lz = 0
                    if z > 0:
                        nmb_lz = NMB_PER_ATOM[pnao.l][z - 1]
                    if functions_used[pnao.l] > nmb_lz or pnao.weight < 1.e-7:
                        pnao.is_nrb = True
                        
            functions_used[:] = 0
            for pnao in self.pnaos:
                if pnao.center_id == cntr_idx + 1 and not pnao.is_nrb:
                    functions_used[pnao.l] += 1
                    
            print(f"center {cntr_idx + 1:3d}: ", end="")
            for L in range(lmx):
                print(f" {functions_used[L]:2d} of {'spdfg'[L]} |", end="")
            print()
            
        sorted_pnaos = None
        
        n_nmb = sum(1 for pnao in self.pnaos if not pnao.is_nrb)
        n_nrb = sum(1 for pnao in self.pnaos if pnao.is_nrb)
        print(f"In total: NMB set has {n_nmb} functions, NRB set has {n_nrb} functions;\n")
        
        if self.verbose_print:
            print("Does PNAO belong to NRB?: ")
            for pnao in self.pnaos:
                print(f"{pnao.is_nrb!s:15}", end="")
            print()
            
        original_nmb_numbers = []
        original_nrb_numbers = []
        for bs_idx, pnao in enumerate(self.pnaos):
            if pnao.is_nrb:
                original_nrb_numbers.append(bs_idx)
            else:
                original_nmb_numbers.append(bs_idx)
                
        if self.verbose_print:
            print("Transforming overlap and SDS matrices to PNAO basis...")
            
        pnao_to_ao_matrix = np.array([pnao.coefs for pnao in self.pnaos])
        self.pnao_labels = create_nao_labels(self.pnaos)
        self.pnao_overlap_matrix = pnao_to_ao_matrix @ self.overlap_matrix @ pnao_to_ao_matrix.T
        
        print_matrix(self.pnao_overlap_matrix, "PNAO overlap matrix:", self.pnao_labels, self.pnao_labels, self.options.PNAO_OverlapMatrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        print(f"Trace of the PNAO overlap matrix: {np.trace(self.pnao_overlap_matrix):.7f} (should be equal to {len(self.pnaos)}, the total number of PNAOs)")
        
        self.sds_pnao = pnao_to_ao_matrix @ self.sds @ pnao_to_ao_matrix.T
        print(f"The trace of SDS matrix in PNAO basis = {np.trace(self.sds_pnao)}")
        
        print_matrix(self.sds_pnao, "The S.D.S matrix in PNAO basis:", self.pnao_labels, self.pnao_labels, self.options.PNAO_SDS_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        print("\nSTEP 3. Weighted orthogonalization of NMB PNAOs\n")
        ow1 = np.zeros((len(self.pnaos), len(self.pnaos)))
        for i in range(n_nrb):
            ow1[original_nrb_numbers[i], original_nrb_numbers[i]] = 1.0
            
        nmb_overlap_matrix = self.pnao_overlap_matrix[np.ix_(original_nmb_numbers, original_nmb_numbers)]
        nmb_sds_matrix = self.sds_pnao[np.ix_(original_nmb_numbers, original_nmb_numbers)]
        
        nmb_labels = []
        for i in range(n_nmb):
            pnao = self.pnaos[original_nmb_numbers[i]]
            nmb_labels.append(f"A{pnao.center_id}: R{pnao.radial_part_id + 1}*{'spdfg'[pnao.l]}({pnao.m})")
            
        print_matrix(nmb_overlap_matrix, "NMB_old overlap:", nmb_labels, nmb_labels, self.options.NMB_old_Overlap_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        print_matrix(nmb_sds_matrix, "NMB_old SDS:", nmb_labels, nmb_labels, self.options.NMB_old_SDS_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        wsw_nmb = np.zeros((n_nmb, n_nmb))
        min_weight = self.pnaos[original_nmb_numbers[0]].weight
        for i in range(n_nmb):
            w_i = self.pnaos[original_nmb_numbers[i]].weight
            if w_i <= min_weight:
                min_weight = w_i
            if w_i <= 1.e-7:
                w_i = 1.e-7
            for j in range(i, n_nmb):
                w_j = self.pnaos[original_nmb_numbers[j]].weight
                if w_j <= 1.e-7:
                    w_j = 1.e-7
                wsw_ij = nmb_overlap_matrix[i, j] * w_i * w_j
                wsw_nmb[i, j] = wsw_ij
                wsw_nmb[j, i] = wsw_ij
                
        print(f"min weight of NMB PNAO = {min_weight}")
        wsw_nmb_inv_sqrt = matrix_sqrt(wsw_nmb, inverse=True)
        
        for j in range(n_nmb):
            w_j = self.pnaos[original_nmb_numbers[j]].weight
            for i in range(n_nmb):
                ow1[original_nmb_numbers[i], original_nmb_numbers[j]] = w_j * wsw_nmb_inv_sqrt[j, i]
                
        overlap_new = ow1 @ self.pnao_overlap_matrix @ ow1.T
        sds_new = ow1 @ self.sds_pnao @ ow1.T
        
        if self.verbose_print:
            print("Weights of NMB PNAOs:")
            s = 0.0
            for i, idx in enumerate(original_nmb_numbers):
                s += sds_new[idx, idx]
                print(f"{i:4d}\t{idx:4d}\t{sds_new[idx, idx]:.15f}")
            print(f"sum of NMB weights = {s:.10f}")
            
        print("\nSTEP 4. Schmidt orthogonalization of NRBs to new NMBs\n")
        nrb = [None] * n_nrb
        for nrb_idx in range(n_nrb):
            orig_idx = original_nrb_numbers[nrb_idx]
            pnao = self.pnaos[orig_idx]
            new_bf = BasisFunction(
                l=pnao.l,
                m=pnao.m,
                center_id=pnao.center_id,
                radial_part_id=pnao.radial_part_id,
                coefs=np.zeros(len(self.pnaos)),
                is_nrb=True
            )
            new_bf.coefs[orig_idx] = 1.0
            nrb[nrb_idx] = new_bf
            
        for nrb_idx in range(n_nrb):
            for nmb_idx in range(n_nmb):
                nrb[nrb_idx].coefs[original_nmb_numbers[nmb_idx]] -= overlap_new[original_nmb_numbers[nmb_idx], original_nrb_numbers[nrb_idx]]
                
        os1 = np.zeros((len(self.pnaos), len(self.pnaos)))
        for b in range(n_nrb):
            for c in range(len(self.pnaos)):
                os1[original_nrb_numbers[b], c] = nrb[b].coefs[c]
        for b in range(n_nmb):
            os1[original_nmb_numbers[b], original_nmb_numbers[b]] = 1.0
            
        overlap_new = os1 @ overlap_new @ os1.T
        sds_new = os1 @ sds_new @ os1.T
        
        print("Diagonal elements of SDS matrix after 1-st Schmidt transformation for new NRB functions:")
        s = 0.0
        for i, idx in enumerate(original_nrb_numbers):
            s += sds_new[idx, idx]
            print(f"{i:4d}\t{idx:4d}\t{sds_new[idx, idx]:.15f}")
        print(f"Sum of NRB 'occupancies' = {s:.10f}")
        
        print("\nSTEP 5. Intracenter naturalization of new NRBs\n")
        on2 = np.zeros((len(self.pnaos), len(self.pnaos)))
        for i in range(n_nmb):
            on2[original_nmb_numbers[i], original_nmb_numbers[i]] = 1.0
            
        nrb_labels = create_nao_labels(nrb)
        nrb_overlap_matrix = overlap_new[np.ix_(original_nrb_numbers, original_nrb_numbers)]
        nrb_sds_matrix = sds_new[np.ix_(original_nrb_numbers, original_nrb_numbers)]
        
        print_matrix(nrb_overlap_matrix, "NRB_old overlap:", nrb_labels, nrb_labels, self.options.NRB_old_Overlap_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())

        nrb_new = intracenter_basis_orthogonalization(nrb, len(self.centers), nrb_sds_matrix, nrb_overlap_matrix, self.verbose_print, self.options.glbPrint)
        if n_nrb > 0:
            u = np.array([bf.coefs for bf in nrb_new])
            nrb_new_overlap = u @ nrb_overlap_matrix @ u.T
        else:
            u = np.zeros((0, len(self.pnaos)))
            nrb_new_overlap = np.zeros((0, 0))
        
        print_matrix(nrb_new_overlap, "NRB_new overlap:", nrb_labels, nrb_labels, self.options.NRB_new_Overlap_Matrix_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        for i in range(n_nrb):
            for j in range(n_nrb):
                on2[original_nrb_numbers[i], original_nrb_numbers[j]] = nrb_new[i].coefs[j]
                
        overlap_new = on2 @ overlap_new @ on2.T
        print_matrix(overlap_new, "Overlap_new", self.pnao_labels, self.pnao_labels, self.options.S_Matrix_after_ON2_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        sds_new = on2 @ sds_new @ on2.T
        print_matrix(sds_new, "SDS_new", self.pnao_labels, self.pnao_labels, self.options.SDS_Matrix_after_ON2_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        if self.verbose_print:
            print("NMB and NRB weights and occupancies after intracenter naturalization of NRBs")
            s = 0.0
            print("LocNum.\tGlobNum.\t<phi|SDS|phi>")
            for i, idx in enumerate(original_nmb_numbers):
                s += sds_new[idx, idx]
                print(f"{i:4d}\t{idx:4d}\t{sds_new[idx, idx]:.10f}")
            print(f"Sum of NMB SDS diagonal terms = {s:.10f}")
            
            s = 0.0
            print("LocNum.\tGlobNum.\t<phi|SDS|phi>\tweight\tNRB label")
            for i, idx in enumerate(original_nrb_numbers):
                s += sds_new[idx, idx]
                print(f"{i:4d}\t{idx:4d}\t{sds_new[idx, idx]:.10f}\t{nrb_new[i].weight:.10f}\t{self.pnao_labels[idx]}")
            print(f"Sum of NRB SDS diagonal terms = {s:.10f}")
            
        print("\nSTEP 6. Weighted orthogonalization of naturalized NRBs\n")
        ow2 = np.zeros((len(self.pnaos), len(self.pnaos)))
        for i in range(n_nmb):
            ow2[original_nmb_numbers[i], original_nmb_numbers[i]] = 1.0
        
        heavy_nrb_threshold = self.options.Heavy_NRB_Threshold.get_double()
        if n_nrb == 0:
            ow_nrb = np.zeros((0, 0))
        elif self.options.NRB_OW_Straightforward.get_boolean():

            if self.verbose_print:
                print("Using a direct WSW orthogonalization for all NRB functions.")
            w_nrb = np.zeros((n_nrb, n_nrb))
            if n_nrb > 0:
                min_weight = nrb_new[0].weight
            else:
                min_weight = 0.0
            for i in range(n_nrb):
                w_nrb[i, i] = nrb_new[i].weight
                if nrb_new[i].weight < min_weight:
                    min_weight = nrb_new[i].weight
            print(f"min weight = {min_weight}")
            
            wsw_nrb = np.zeros((n_nrb, n_nrb))
            for i in range(n_nrb):
                for j in range(n_nrb):
                    wsw_nrb[i, j] = nrb_new[i].weight * nrb_new_overlap[i, j] * nrb_new[j].weight
                    
            ow_nrb = w_nrb @ matrix_sqrt(wsw_nrb, inverse=True)
        else:
            if self.verbose_print:
                print(f"Using a WSWSL orthogonalization for NRB functions with occ. Threshold = {heavy_nrb_threshold:.3E}")
            max_weight = 0.0
            for i in range(n_nrb):
                if nrb_new[i].weight > max_weight:
                    max_weight = nrb_new[i].weight
            print(f"Maximum weight of NRB function = {max_weight:.5E}")
            
            n_heavy = 0
            n_light = 0
            if n_nrb > 0:
                min_weight = nrb_new[0].weight
            else:
                min_weight = 0.0
                
            for i in range(n_nrb):
                if nrb_new[i].weight > heavy_nrb_threshold:
                    n_heavy += 1
                else:
                    n_light += 1
                if nrb_new[i].weight < min_weight:
                    min_weight = nrb_new[i].weight
                    
            print(f"min weight = {min_weight}")
            if self.verbose_print:
                print(f"'Heavily occupied NRB set' has {n_heavy} of {n_nrb} functions")
                
            heavy_nrbs = []
            light_nrbs = []
            w_nrb_heavy = np.zeros((n_heavy, n_heavy))
            
            for i in range(n_nrb):
                if nrb_new[i].weight > heavy_nrb_threshold:
                    heavy_nrbs.append(i)
                    w_nrb_heavy[len(heavy_nrbs)-1, len(heavy_nrbs)-1] = nrb_new[i].weight
                else:
                    light_nrbs.append(i)
                    
            heavy_nrb_overlap = nrb_new_overlap[np.ix_(heavy_nrbs, heavy_nrbs)]
            wsw_nrb_heavy = np.zeros((n_heavy, n_heavy))
            for i in range(n_heavy):
                for j in range(i, n_heavy):
                    tmp = heavy_nrb_overlap[i, j] * nrb_new[heavy_nrbs[i]].weight * nrb_new[heavy_nrbs[j]].weight
                    wsw_nrb_heavy[i, j] = tmp
                    wsw_nrb_heavy[j, i] = tmp
                    
            ow_nrb_heavy = w_nrb_heavy @ matrix_sqrt(wsw_nrb_heavy, inverse=True)
            
            ow_nrb = np.zeros((n_nrb, n_nrb))
            ow_nrb[np.ix_(heavy_nrbs, heavy_nrbs)] = ow_nrb_heavy
            for i in range(n_light):
                ow_nrb[light_nrbs[i], light_nrbs[i]] = 1.0
                
            new_nrb_overlap = ow_nrb @ nrb_new_overlap @ ow_nrb.T
            print_matrix(new_nrb_overlap, "new_NRB_Overlap before Schmidt", None, None, self.options.NRB_Overlap_after_OW_heavy_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
            
            os = np.eye(n_nrb)
            for i in range(n_light):
                for j in range(n_heavy):
                    os[light_nrbs[i], heavy_nrbs[j]] -= new_nrb_overlap[light_nrbs[i], heavy_nrbs[j]]
                    
            ow_nrb = ow_nrb @ os.T
            new_nrb_overlap = ow_nrb @ nrb_new_overlap @ ow_nrb.T
            print_matrix(new_nrb_overlap, "new_NRB_Overlap after Schmidt", None, None, self.options.NRB_Overlap_after_OS2_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
            
            if n_light > 0:
                light_overlap = new_nrb_overlap[np.ix_(light_nrbs, light_nrbs)]
                ol = matrix_sqrt(light_overlap, inverse=True)
                olnrb = np.eye(n_nrb)
                olnrb[np.ix_(light_nrbs, light_nrbs)] = ol
                ow_nrb = ow_nrb @ olnrb.T
                
        nrb_overlap_test = ow_nrb @ nrb_new_overlap @ ow_nrb.T
        print_matrix(nrb_overlap_test, "NRB new overlap", None, None, self.options.NRB_Overlap_after_OW2_final_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        diff_norm2 = 0.0
        for i in range(n_nrb):
            for j in range(n_nrb):
                val = nrb_overlap_test[i, j]
                if i == j:
                    diff_norm2 += (val - 1.0)**2
                else:
                    diff_norm2 += val**2
        print(f"|S_NRB - 1| = {np.sqrt(diff_norm2):.5E} (Should be VERY close to zero!)")
        
        ow2[np.ix_(original_nrb_numbers, original_nrb_numbers)] = ow_nrb.T
        
        ao_names = create_nao_labels(self.basis)
        print_matrix(ow2, "2-nd WSW transformation matrix:", None, ao_names, self.options.OW2_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        overlap_new = ow2 @ overlap_new @ ow2.T
        print_matrix(overlap_new, "Overlap matrix in almost-NAO basis (obtained after OW2 step):", self.pnao_labels, self.pnao_labels, self.options.S_Matrix_after_OW2_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        sds_new = ow2 @ sds_new @ ow2.T
        print_matrix(sds_new, "SDS matrix in almost-NAO basis (obtained after OW2 step):", self.pnao_labels, self.pnao_labels, self.options.SDS_Matrix_after_OW2_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        if self.verbose_print:
            print("Diagonal elements of SDS before the final intracenter naturalization step:")
            for cntr_idx in range(len(self.centers)):
                cntr_occ = sum(sds_new[i, i] for i in range(len(self.pnaos)) if self.pnaos[i].center_id - 1 == cntr_idx)
                print(f"center {cntr_idx + 1}: {cntr_occ:.4f}")
                
        print("\nSTEP 7. Final Intracenter Natural Transformation withing the full set of functions\n")
        self.options.dont_print_matrices = False
        self.nao = intracenter_basis_orthogonalization(self.pnaos, len(self.centers), sds_new, overlap_new, self.verbose_print, self.options.glbPrint)
        
        bf_new_to_nao = np.array([bf.coefs for bf in self.nao])
        nao_to_pnao = bf_new_to_nao @ ow2 @ on2 @ os1 @ ow1
        self.nao_to_ao = nao_to_pnao @ pnao_to_ao_matrix
        
        self.overlap_nao = bf_new_to_nao @ overlap_new @ bf_new_to_nao.T
        diff_norm2 = 0.0
        for i in range(len(self.nao)):
            for j in range(len(self.nao)):
                val = self.overlap_nao[i, j]
                if i == j:
                    diff_norm2 += (val - 1.0)**2
                else:
                    diff_norm2 += val**2
        print(f"SQRT{{ SUM[(NaoOverlap_ij - delta_ij)^2] }} = {np.sqrt(diff_norm2):.2e} (should be VERY close to zero)")
        
        max_offdiag = 0.0
        for i in range(len(self.nao)):
            for j in range(i + 1, len(self.nao)):
                if abs(self.overlap_nao[i, j]) > max_offdiag:
                    max_offdiag = abs(self.overlap_nao[i, j])
        print(f"max_offdiag = {max_offdiag:.2e} (should be VERY close to zero)")
        
        if np.sqrt(diff_norm2) > 1.0e-10:
            print("WARNING: NAOs seem to be not strictly orthogonalized!")
            
        self.sds_nao = bf_new_to_nao @ sds_new @ bf_new_to_nao.T
        
        print(STARS)
        print("\nFinal NAO occupancies and leading AO terms:\n")
        old_weights = np.zeros(len(self.nao))
        print(f"{'NAO #':>5s} {'Name':20s} {'Occupancy':>10s} {'Leading term'}")
        for i in range(len(self.nao)):
            old_weights[i] = self.nao[i].weight
            self.nao[i].weight = self.sds_nao[i, i]
            self.nao[i].is_nrb = self.pnaos[i].is_nrb
            
            print(f"{i + 1:5d} {self.pnao_labels[i]:20s} {self.nao[i].weight:10.7f} ", end="")
            
            ao_main = np.argmax(np.abs(self.nao_to_ao[i, :]))
            tmp = self.nao_to_ao[i, ao_main]
            print(f"({tmp:.2f})*BF[{ao_main + 1} = {ao_names[ao_main]}]")
            
        print_matrix(self.sds_nao, "S.D.S in NAO basis:", self.pnao_labels, self.pnao_labels, self.options.SDS_NAO_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        print_matrix(self.nao_to_ao, "AO-to-NAO transformation matrix:", self.pnao_labels, ao_names, self.options.NAO2AO_File.get_string(), self.options.MatrixFloatNumberFormat.get_string(), self.options.MatrixLineWidth.get_int())
        
        print(f"trace = {np.trace(self.sds_nao)}")
        
        print("\nFinal electron populations and NPA charges:\n")
        print("Center\tNuclear\t Electron \t NMB \tNPA ")
        print(" \t charge\t population\t population\tcharge")
        
        self.npa_charges = np.zeros(len(self.centers))
        for cntr_idx, center in enumerate(self.centers):
            cntr_occ = 0.0
            nmb_contrib = 0.0
            for i in range(len(self.nao)):
                if self.nao[i].center_id - 1 == cntr_idx:
                    cntr_occ += self.nao[i].weight
                    if not self.nao[i].is_nrb:
                        nmb_contrib += self.nao[i].weight
            self.npa_charges[cntr_idx] = center.z - cntr_occ
            print(f"{center.name}{cntr_idx + 1:5d}\t{center.z:7.1f}\t{cntr_occ:11.7f}\t{nmb_contrib:12.7f}\t{self.npa_charges[cntr_idx]:13.10f}")

        print()
        print("Angular momentum contributions of the total atomic population:\n")
        print(" Cntr", end="")
        l_contrib = np.zeros(6)
        for l in range(len(l_contrib)):
            print(f"{'spdfgh'[l]:12s}", end="")
        print()

        for cntr_idx, center in enumerate(self.centers):
            l_contrib[:] = 0.0
            for i in range(len(self.nao)):
                if self.nao[i].center_id - 1 == cntr_idx:
                    l_contrib[self.nao[i].l] += self.nao[i].weight
            print(f"{center.name}{cntr_idx + 1:7d}", end="")

            for l in range(len(l_contrib)):
                print(f"{l_contrib[l]:12.7f}", end="")
            print()
            
        print()
        return True


def bond_indexes_generic(n_atoms: int, basis_functions: list[BasisFunction], s05ds05: np.ndarray) -> np.ndarray:
    """Calculates generic Wiberg bond indices."""
    result = np.zeros((n_atoms, n_atoms))
    for i in range(len(basis_functions)):
        for j in range(i + 1, len(basis_functions)):
            result[basis_functions[i].center_id - 1, basis_functions[j].center_id - 1] += s05ds05[i, j] * s05ds05[j, i]
            
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            result[j, i] = result[i, j]
            
    for i in range(n_atoms):
        v = 0.0
        result[i, i] = 0.0
        for j in range(n_atoms):
            v += result[i, j]
        result[i, i] = v
        
    return result


def run_npa(molden: MoldenFile, options: NPAOptions | None = None) -> NPAResult:
    """Execute the full 7-step NPA procedure.

    Steps:
    1. Produce PNAOs via intracenter diagonalization of SDS
    2. Split PNAOs into NMB / NRB
    3. Weighted orthogonalization of NMB PNAOs
    4. Schmidt orthogonalization of NRBs to new NMBs
    5. Intracenter naturalization of new NRBs
    6. Weighted orthogonalization of naturalized NRBs
    7. Final intracenter natural transformation
    
    Returns NPAResult with all matrices and charges.
    """
    if options is None:
        options = NPAOptions()
        
    npa = NPA(options)
    if not npa.load_mo_from_molden(molden):
        raise RuntimeError("Failed to load MOs from MOLDEN file.")
    if not npa.create_naos():
        raise RuntimeError("Failed to create NAOs.")
        
    res = NPAResult(
        nao=npa.nao,
        npa_charges=npa.npa_charges,
        overlap_matrix=npa.overlap_matrix,
        density_matrix=npa.d_global,
        sds=npa.sds,
        sds_nao=npa.sds_nao,
        nao_to_ao=npa.nao_to_ao,
        overlap_nao=npa.overlap_nao,
        wiberg_indices=bond_indexes_generic(len(npa.centers), npa.nao, npa.sds_nao)
    )
    return res