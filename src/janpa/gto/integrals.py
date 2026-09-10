"""Main integral engine: builds overlap and dipole matrices in the pure spherical-harmonic basis using the cylindrical-coordinates approach. Replaces JGintsCyl.java. This is the most performance-critical module. """
from __future__ import annotations

import numpy as np
from typing import Any

from .basis import BasisFunction, RadialPart, AtomicCenter
from .polynom3d import Polynom3D
from .polynom_rho_z import Polynom_rho_z
from .spherical import ALEG_NORMALIZER
from .overlap import primitive_int_1d_sphr
from .wigner import new_axis_coords, create_sincos, wigner

# Mappings between canonical m-ordering (-L...L) and MOLDEN m-ordering.
# MOLDEN ordering for L=1 is (y, z, x) which corresponds to m = -1, 0, 1.
P2moldenP = [
    [0],
    [1, -1, 0],
    [-2, -1, 0, 1, 2],
    [-3, -2, -1, 0, 1, 2, 3],
    [-4, -3, -2, -1, 0, 1, 2, 3, 4],
]

moldenP2P = [
    [0],
    [0, 1, -1],
    [-2, -1, 0, 1, 2],
    [-3, -2, -1, 0, 1, 2, 3],
    [-4, -3, -2, -1, 0, 1, 2, 3, 4],
]


def theta_phi(ra: np.ndarray, rb: np.ndarray) -> np.ndarray:
    """Compute polar angles (theta, phi) of vector RB - RA."""
    dr = rb - ra
    theta = np.arctan2(np.sqrt(dr[1]**2 + dr[0]**2), dr[2])
    phi = np.arctan2(dr[1], dr[0])
    return np.array([theta, phi])


def _wigner_transpose(w: np.ndarray, l: int) -> np.ndarray:
    """Transpose of the Wigner D-matrix (equivalent to inverse for orthogonal matrices)."""
    return w.T


def _wigner_2_molden(w: np.ndarray, l: int) -> np.ndarray:
    """Convert Wigner matrix from canonical P ordering to MOLDEN ordering."""
    result = np.zeros((2 * l + 1, 2 * l + 1))
    for m in range(-l, l + 1):
        for m_p in range(-l, l + 1):
            m_idx_p = moldenP2P[l][l + m]
            mp_idx_p = moldenP2P[l][l + m_p]
            result[l + m, l + m_p] = w[l + m_idx_p, l + mp_idx_p]
    return result


def _transform_dipole_integrals(
    dipole_storage: np.ndarray, la: int, lb: int, w_la: np.ndarray, w_lb: np.ndarray
) -> np.ndarray:
    """Transform the dipole storage matrix using Wigner rotation vectors."""
    result = np.zeros(3)
    result[2] = w_la[la] * dipole_storage[la, lb] * w_lb[lb]
    
    if la > 0:
        result[1] = w_la[la - 1] * dipole_storage[la - 1, lb] * w_lb[lb]
        result[0] = w_la[la + 1] * dipole_storage[la + 1, lb] * w_lb[lb]
    if lb > 0:
        result[1] += w_la[la] * dipole_storage[la, lb - 1] * w_lb[lb - 1]
        result[0] += w_la[la] * dipole_storage[la, lb + 1] * w_lb[lb + 1]

    l_min = min(la, lb)
    for m in range(1, l_min + 1):
        result[2] += w_la[la + m] * dipole_storage[la + m, lb + m] * w_lb[lb + m]
        result[2] += w_la[la - m] * dipole_storage[la - m, lb - m] * w_lb[lb - m]

    for m in range(1, min(la - 1, lb) + 1):
        result[1] += w_la[la - m - 1] * dipole_storage[la - m - 1, lb + m] * w_lb[lb + m]
        result[1] += w_la[la + m + 1] * dipole_storage[la + m + 1, lb - m] * w_lb[lb - m]
        result[0] += w_la[la + m + 1] * dipole_storage[la + m + 1, lb + m] * w_lb[lb + m]
        result[0] += w_la[la - m - 1] * dipole_storage[la - m - 1, lb - m] * w_lb[lb - m]

    for m in range(1, min(la, lb - 1) + 1):
        result[1] += w_la[la - m] * dipole_storage[la - m, lb + m + 1] * w_lb[lb + m + 1]
        result[1] += w_la[la + m] * dipole_storage[la + m, lb - m - 1] * w_lb[lb - m - 1]
        result[0] += w_la[la + m] * dipole_storage[la + m, lb + m + 1] * w_lb[lb + m + 1]
        result[0] += w_la[la - m] * dipole_storage[la - m, lb - m - 1] * w_lb[lb - m - 1]

    return result


def _set_atomic_dipoles_non_diag(
    rpart_m2bf: list[list[int]], la: int, lb: int, rp1: int, rp2: int,
    rp_product_lalb: float, dipole_matrix: np.ndarray
) -> None:
    """Compute off-diagonal same-center atomic dipole integrals."""
    for m in range(min(la, lb) + 1):
        if la > lb:
            ylm_factor = (lb - m + 1.0) * (lb + m + 1.0) / (2 * lb + 1.0) / (2 * lb + 3.0)
        else:
            ylm_factor = (lb - m) * (lb + m) / (2 * lb - 1.0) / (2 * lb + 1.0)
        ylm_factor = np.sqrt(ylm_factor)
        
        bf1 = rpart_m2bf[rp1][la + m]
        bf2 = rpart_m2bf[rp2][lb + m]
        dipole_matrix[bf1, bf2, 2] = ylm_factor * rp_product_lalb
        
        if m > 0:
            bf1 = rpart_m2bf[rp1][la - m]
            bf2 = rpart_m2bf[rp2][lb - m]
            dipole_matrix[bf1, bf2, 2] = ylm_factor * rp_product_lalb

    rho_prefac = 0.5 * rp_product_lalb / np.sqrt((la + lb + 2.0) * (la + lb))
    
    if la > 0:
        bf1p = rpart_m2bf[rp1][la + 1]
        bf1m = rpart_m2bf[rp1][la - 1]
        bf2 = rpart_m2bf[rp2][lb + 0]
        if la > lb:
            rho_matr_el_half = np.sqrt(2 * (lb + 2) * (lb + 1)) * rho_prefac
        else:
            rho_matr_el_half = -np.sqrt(2 * lb * (lb - 1)) * rho_prefac
            
        dipole_matrix[bf1p, bf2, 0] = rho_matr_el_half
        dipole_matrix[bf1m, bf2, 1] = rho_matr_el_half

    if lb > 0:
        bf1 = rpart_m2bf[rp1][la + 0]
        bf2p = rpart_m2bf[rp2][lb + 1]
        bf2m = rpart_m2bf[rp2][lb - 1]
        if la > lb:
            rho_matr_el_half = -np.sqrt(2 * lb * (lb + 1)) * rho_prefac
        else:
            rho_matr_el_half = np.sqrt(2 * lb * (lb + 1)) * rho_prefac
            
        dipole_matrix[bf1, bf2p, 0] = rho_matr_el_half
        dipole_matrix[bf1, bf2m, 1] = rho_matr_el_half

    for m in range(1, min(la, lb - 1) + 1):
        if la > lb:
            rho_matr_el_half = -np.sqrt((lb - m + 1) * (lb - m)) * rho_prefac
        else:
            rho_matr_el_half = np.sqrt((lb + m + 1) * (lb + m)) * rho_prefac
            
        bf1p = rpart_m2bf[rp1][la + m]
        bf1m = rpart_m2bf[rp1][la - m]
        bf2p = rpart_m2bf[rp2][lb + m + 1]
        bf2m = rpart_m2bf[rp2][lb - m - 1]
        
        dipole_matrix[bf1p, bf2p, 0] = rho_matr_el_half
        dipole_matrix[bf1p, bf2m, 1] = rho_matr_el_half
        dipole_matrix[bf1m, bf2m, 0] = rho_matr_el_half
        dipole_matrix[bf1m, bf2p, 1] = -rho_matr_el_half

    for m in range(2, min(la, lb + 1) + 1):
        if la > lb:
            rho_matr_el_half = np.sqrt((lb + m + 1) * (lb + m)) * rho_prefac
        else:
            rho_matr_el_half = -np.sqrt((lb - m + 1) * (lb - m)) * rho_prefac
            
        bf1p = rpart_m2bf[rp1][la + m]
        bf1m = rpart_m2bf[rp1][la - m]
        bf2p = rpart_m2bf[rp2][lb + m - 1]
        bf2m = rpart_m2bf[rp2][lb - m + 1]
        
        dipole_matrix[bf1p, bf2p, 0] = rho_matr_el_half
        dipole_matrix[bf1p, bf2m, 1] = -rho_matr_el_half
        dipole_matrix[bf1m, bf2m, 0] = rho_matr_el_half
        dipole_matrix[bf1m, bf2p, 1] = rho_matr_el_half


class IntegralEngine:
    """Builds overlap and dipole matrices from a MOLDEN-like basis specification."""

    def __init__(
        self,
        basis: list[BasisFunction],
        radial_parts: list[RadialPart],
        centers: list[AtomicCenter],
    ):
        self.basis = basis
        self.radial_parts = radial_parts
        self.centers = centers
        self.overlap_matrix: np.ndarray | None = None
        self.dipole_matrix: np.ndarray | None = None

    @staticmethod
    def primitive_g_overlap_int(
        zeta1: float, l1: int, za: float, addit_r_pwr1: int,
        zeta2: float, l2: int, zb: float, addit_r_pwr2: int,
        l_min: int, do_dipoles: bool = True,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
        """Primitive Gaussian overlap integral along the internuclear axis."""
        zeta_eff = zeta1 + zeta2
        z0 = (zeta1 * za + zeta2 * zb) / zeta_eff
        prefactor = np.exp(-zeta1 * zeta2 * (zb - za)**2 / zeta_eff)
        sz = l1 + l2 + addit_r_pwr1 + addit_r_pwr2

        z_ints = np.zeros(sz + 2)
        z_ints[0] = np.sqrt(np.pi / zeta_eff)
        for i in range(2, sz + 2, 2):
            z_ints[i] = z_ints[i - 2] * (i - 1) / 2.0 / zeta_eff

        rho_ints = np.zeros(sz + 3)
        rho_ints[1] = 0.5 / zeta_eff
        for i in range(3, sz + 3, 2):
            rho_ints[i] = rho_ints[i - 2] * (i - 1) / 2.0 / zeta_eff
        rho_ints[0] = np.sqrt(np.pi / zeta_eff) * 0.5
        for i in range(2, sz + 3, 2):
            rho_ints[i] = rho_ints[i - 2] * (i - 1) / 2.0 / zeta_eff

        overlap_diags = np.zeros(2 * l_min + 1)
        if do_dipoles:
            z_dipole = np.zeros(l_min + 1)
            i_plus = np.zeros(l_min + 1)
            i_minus = np.zeros(l_min + 2)
        else:
            z_dipole = None
            i_plus = None
            i_minus = None

        M = Polynom_rho_z(sz, sz)
        M.mul_r2(addit_r_pwr1 // 2, za - z0)
        M.mul_r2(addit_r_pwr2 // 2, zb - z0)

        for m in range(l1 + 1):
            if m - 1 > l2:
                break
            
            M.push_to_stack()
            M.mul_assoc_leg(l1, m, za - z0)
            
            diag_available = (m <= l1) and (m <= l2)
            iplus_available = (m + 1 <= l2) and do_dipoles
            iminus_available = (m >= 1) and (m - 1 <= l2) and do_dipoles
            
            other_terms1 = prefactor * ALEG_NORMALIZER[l1][l1 + m]
            
            if diag_available:
                if iplus_available or iminus_available:
                    M.push_to_stack()
                M.mul_assoc_leg(l2, m, zb - z0)
                
                res_overlap = 0.0
                res_z_dipole = 0.0
                for i in range(M.last_active_z_power + 1):
                    tmp = 0.0
                    for j in range(M.last_active_rho_power + 1):
                        tmp += M.cf_matrix[i, j] * rho_ints[j + 2 * m + 1]
                    res_overlap += z_ints[i] * tmp
                    res_z_dipole += z_ints[i + 1] * tmp
                
                phi_int = np.pi
                if m == 0:
                    phi_int *= 2.0
                    
                other_terms = other_terms1 * phi_int * ALEG_NORMALIZER[l2][l2 + m]
                res_overlap *= other_terms
                res_z_dipole *= other_terms
                
                overlap_diags[l_min - m] = res_overlap
                overlap_diags[l_min + m] = res_overlap
                
                if do_dipoles:
                    z_dipole[m] = res_z_dipole
                    
                if iplus_available or iminus_available:
                    M.restore_from_stack(True)
                    
            if iplus_available:
                if iminus_available:
                    M.push_to_stack()
                M.mul_assoc_leg(l2, m + 1, zb - z0)
                
                res_iplus = 0.0
                for i in range(M.last_active_z_power + 1):
                    tmp = 0.0
                    for j in range(M.last_active_rho_power + 1):
                        tmp += M.cf_matrix[i, j] * rho_ints[j + 2 * m + 3]
                    res_iplus += z_ints[i] * tmp
                    
                res_iplus *= other_terms1
                res_iplus *= ALEG_NORMALIZER[l2][l2 + m + 1]
                i_plus[m] = res_iplus * 0.5
                
                if iminus_available:
                    M.restore_from_stack(True)
                    
            if iminus_available:
                M.mul_assoc_leg(l2, m - 1, zb - z0)
                
                res_iminus = 0.0
                for i in range(M.last_active_z_power + 1):
                    tmp = 0.0
                    for j in range(M.last_active_rho_power + 1):
                        tmp += M.cf_matrix[i, j] * rho_ints[j + 2 * m + 1]
                    res_iminus += z_ints[i] * tmp
                    
                res_iminus *= other_terms1
                res_iminus *= ALEG_NORMALIZER[l2][l2 + m - 1]
                i_minus[m] = res_iminus * 0.5
                
            M.restore_from_stack(True)

        return overlap_diags, z_dipole, i_plus, i_minus

    def build_overlap_and_dipole(self) -> tuple[np.ndarray, np.ndarray]:
        """Compute the full overlap matrix S and dipole integral tensor D_mu."""
        n_radial_parts = len(self.radial_parts)
        n_basis = len(self.basis)
        n_centers = len(self.centers)

        # Group basis functions by radial part
        bfns_of_rp = [[] for _ in range(n_radial_parts)]
        for b, bf in enumerate(self.basis):
            bfns_of_rp[bf.radial_part_id].append(b)

        # Map canonical m to basis function index for each radial part
        rpart_m2bf = [[0] * (2 * self.radial_parts[r].l_used_with + 1) for r in range(n_radial_parts)]
        for r in range(n_radial_parts):
            l_rp = self.radial_parts[r].l_used_with
            for bf in bfns_of_rp[r]:
                true_m = moldenP2P[l_rp][l_rp + self.basis[bf].m]
                rpart_m2bf[r][l_rp + true_m] = bf

        # Group radial parts by atomic center
        rpts_of_atom = [[] for _ in range(n_centers)]
        max_l = 0
        max_r2_pwr = 0
        for rp, rp_obj in enumerate(self.radial_parts):
            c = rp_obj.center_id - 1
            rpts_of_atom[c].append(rp)
            if rp_obj.addit_r_power > 2 * max_r2_pwr:
                max_r2_pwr = rp_obj.addit_r_power // 2
            if rp_obj.l_used_with > max_l:
                max_l = rp_obj.l_used_with

        overlap_matrix = np.zeros((n_basis, n_basis))
        dipole_matrix = np.zeros((n_basis, n_basis, 3))
        produce_dipoles = True

        # 1. Different centers
        for a in range(n_centers):
            for b in range(a + 1, n_centers):
                ra = self.centers[a].r0.copy()
                rb = self.centers[b].r0.copy()
                rb -= ra
                ra[:] = 0.0
                
                ti = theta_phi(ra, rb)
                rot_matrix = new_axis_coords(ti[0], ti[1])
                d_ab = np.linalg.norm(rb)
                rb_new = np.array([0.0, 0.0, d_ab])
                sincos = create_sincos(ti[0], ti[1])
                
                l_max_a = max((self.radial_parts[rp].l_used_with for rp in rpts_of_atom[a]), default=0)
                l_max_b = max((self.radial_parts[rp].l_used_with for rp in rpts_of_atom[b]), default=0)
                l_max = max(l_max_a, l_max_b)
                
                W = []
                for l in range(l_max + 1):
                    w_l = wigner(sincos, l)
                    W.append(_wigner_transpose(w_l, l))
                    
                for irp1, rp1 in enumerate(rpts_of_atom[a]):
                    la = self.radial_parts[rp1].l_used_with
                    for irp2, rp2 in enumerate(rpts_of_atom[b]):
                        lb = self.radial_parts[rp2].l_used_with
                        l_min = min(la, lb)
                        
                        overlap_diags = np.zeros(2 * l_min + 1)
                        dipole_storage = np.zeros((2 * la + 1, 2 * lb + 1))
                        
                        for zeta1, cf1 in zip(self.radial_parts[rp1].exponents, self.radial_parts[rp1].coefs):
                            for zeta2, cf2 in zip(self.radial_parts[rp2].exponents, self.radial_parts[rp2].coefs):
                                res_ov, res_z, res_ip, res_im = self.primitive_g_overlap_int(
                                    zeta1, la, 0.0, self.radial_parts[rp1].addit_r_power,
                                    zeta2, lb, d_ab, self.radial_parts[rp2].addit_r_power,
                                    l_min, produce_dipoles
                                )
                                
                                cf_prod12 = cf1 * cf2
                                for m in range(-l_min, l_min + 1):
                                    overlap_diags[l_min + m] += res_ov[l_min + m] * cf_prod12
                                    
                                if produce_dipoles:
                                    for m in range(l_min + 1):
                                        res_z[m] += d_ab * zeta2 / (zeta1 + zeta2) * res_ov[l_min + m]
                                        dipole_storage[la + m, lb + m] += res_z[m] * cf_prod12
                                        dipole_storage[la - m, lb - m] = dipole_storage[la + m, lb + m]
                                        
                                    if lb > 0:
                                        dipole_storage[la, lb - 1] += cf_prod12 * res_ip[0] * 2 * np.pi
                                        dipole_storage[la, lb + 1] = dipole_storage[la, lb - 1]
                                        
                                    if la > 0:
                                        dipole_storage[la - 1, lb] += cf_prod12 * res_im[1] * 2 * np.pi
                                        dipole_storage[la + 1, lb] = dipole_storage[la - 1, lb]
                                        
                                    for m in range(2, min(la, lb + 1) + 1):
                                        if m <= la and m - 1 <= lb:
                                            iminus_term = res_im[m] * cf_prod12 * np.pi
                                            dipole_storage[la + m, lb + m - 1] += iminus_term
                                            dipole_storage[la + m, lb - m + 1] -= iminus_term
                                            dipole_storage[la - m, lb - m + 1] += iminus_term
                                            dipole_storage[la - m, lb + m - 1] += iminus_term
                                            
                                    for m in range(1, min(la, lb - 1) + 1):
                                        if m <= la and m + 1 <= lb:
                                            iplus_term = res_ip[m] * cf_prod12 * np.pi
                                            dipole_storage[la + m, lb + m + 1] += iplus_term
                                            dipole_storage[la + m, lb - m - 1] += iplus_term
                                            dipole_storage[la - m, lb - m - 1] += iplus_term
                                            dipole_storage[la - m, lb + m + 1] -= iplus_term
                                            
                        for bf1 in bfns_of_rp[rp1]:
                            fixed_ma = la + moldenP2P[la][la + self.basis[bf1].m]
                            for bf2 in bfns_of_rp[rp2]:
                                fixed_mb = lb + moldenP2P[lb][lb + self.basis[bf2].m]
                                
                                result_overlap = 0.0
                                for m in range(-l_min, l_min + 1):
                                    result_overlap += overlap_diags[l_min + m] * W[la][fixed_ma, la + m] * W[lb][fixed_mb, lb + m]
                                    
                                overlap_matrix[bf1, bf2] = result_overlap
                                overlap_matrix[bf2, bf1] = result_overlap
                                
                                if produce_dipoles:
                                    dipoles = _transform_dipole_integrals(
                                        dipole_storage, la, lb, W[la][fixed_ma, :], W[lb][fixed_mb, :]
                                    )
                                    dipoles_orig_ax = rot_matrix @ dipoles
                                    dipoles_orig_ax += self.centers[a].r0 * result_overlap
                                    
                                    dipole_matrix[bf1, bf2, :] = dipoles_orig_ax
                                    dipole_matrix[bf2, bf1, :] = dipoles_orig_ax

        # 2. Same center
        for a in range(n_centers):
            for irp1, rp1 in enumerate(rpts_of_atom[a]):
                la = self.radial_parts[rp1].l_used_with
                addit_r_pwr1 = self.radial_parts[rp1].addit_r_power
                for irp2 in range(irp1, len(rpts_of_atom[a])):
                    rp2 = rpts_of_atom[a][irp2]
                    lb = self.radial_parts[rp2].l_used_with
                    
                    do_calc_overlap = (lb == la)
                    do_calc_dipole = (abs(lb - la) == 1) and produce_dipoles
                    if not (do_calc_overlap or do_calc_dipole):
                        continue
                        
                    addit_r_pwr2 = self.radial_parts[rp2].addit_r_power
                    rp_product = 0.0
                    rp_product_lalb = 0.0
                    
                    for zeta1, cf1 in zip(self.radial_parts[rp1].exponents, self.radial_parts[rp1].coefs):
                        for zeta2, cf2 in zip(self.radial_parts[rp2].exponents, self.radial_parts[rp2].coefs):
                            if do_calc_overlap:
                                rp_product += primitive_int_1d_sphr(
                                    2 + 2 * la + addit_r_pwr1 + addit_r_pwr2, zeta1 + zeta2
                                ) * cf1 * cf2
                            if do_calc_dipole:
                                rp_product_lalb += primitive_int_1d_sphr(
                                    2 + la + lb + 1 + addit_r_pwr1 + addit_r_pwr2, zeta1 + zeta2
                                ) * cf1 * cf2
                                
                    if do_calc_overlap:
                        for bf1 in bfns_of_rp[rp1]:
                            m1 = self.basis[bf1].m
                            for bf2 in bfns_of_rp[rp2]:
                                m2 = self.basis[bf2].m
                                if m2 == m1:
                                    overlap_matrix[bf1, bf2] = rp_product
                                    overlap_matrix[bf2, bf1] = rp_product
                                    
                    if do_calc_dipole:
                        _set_atomic_dipoles_non_diag(rpart_m2bf, la, lb, rp1, rp2, rp_product_lalb, dipole_matrix)

        # 3. Origin shift for same-center dipoles
        if produce_dipoles:
            for a in range(n_centers):
                r0 = self.centers[a].r0
                for irp1 in range(len(rpts_of_atom[a])):
                    rp1 = rpts_of_atom[a][irp1]
                    for irp2 in range(irp1, len(rpts_of_atom[a])):
                        rp2 = rpts_of_atom[a][irp2]
                        for bf1 in bfns_of_rp[rp1]:
                            for bf2 in bfns_of_rp[rp2]:
                                if bf2 < bf1:
                                    continue
                                dipole_matrix[bf1, bf2, :] += r0 * overlap_matrix[bf1, bf2]
                                if bf2 != bf1:
                                    dipole_matrix[bf2, bf1, :] += r0 * overlap_matrix[bf1, bf2]

        self.overlap_matrix = overlap_matrix
        self.dipole_matrix = dipole_matrix
        return overlap_matrix, dipole_matrix