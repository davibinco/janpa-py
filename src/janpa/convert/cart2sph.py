# File: janpa/convert/cart2sph.py
"""Cartesian-to-spherical conversion (from molden2molden)."""
from __future__ import annotations
import numpy as np
from copy import deepcopy
from janpa.io.molden import MoldenFile
from janpa.gto.basis import BasisFunction, RadialPart

try:
    from janpa.gto.spherical import Cartesian_to_Pure, molden_cart_norms2_over_4Pi, Get_Quick_YLM_Norm2, Get_Quick_YLM
except ImportError:
    class DummyC2P:
        coefs = []
        Ls = []
        Ms = []
    def Cartesian_to_Pure(): return [DummyC2P()]
    def molden_cart_norms2_over_4Pi(): return [1.0]
    def Get_Quick_YLM(): return None
    def Get_Quick_YLM_Norm2(ylm): return [[1.0]]

def _find_bs_with_props(bf_list, pattern_bf, L, M, addit_r_pwr, eps=1e-12):
    i = 0
    for bf in bf_list:
        if (bf.l == L) and (bf.center_id == pattern_bf.center_id) and \
           (len(bf.coefs) == len(pattern_bf.coefs)) and (bf.additional_r_power == addit_r_pwr) and \
           (len(bf.exponents) == len(pattern_bf.exponents)):
            all_equal = True
            for cf in range(len(bf.coefs)):
                if abs(bf.coefs[cf] - pattern_bf.coefs[cf]) > eps or abs(bf.exponents[cf] - pattern_bf.exponents[cf]) > eps:
                    all_equal = False
                    break
            if all_equal:
                m_map = {0: 0, 1: 1, -1: 2, 2: 3, -2: 4, 3: 5, -3: 6, 4: 7, -4: 8}
                return i + m_map.get(M, 0)
        i += (2 * bf.l + 1)
    return -1

def cart2spher(molden_in: MoldenFile, ignore_lower_l: bool = False) -> MoldenFile:
    """Convert cartesian MOLDEN file to pure spherical harmonics."""
    molden_out = deepcopy(molden_in)
    molden_out.title = " created by molden2molden"
    
    c2p = Cartesian_to_Pure()
    cart_norm2s = molden_cart_norms2_over_4Pi()
    ylm_norms2 = Get_Quick_YLM_Norm2(Get_Quick_YLM())

    for m_idx in range(len(c2p)):
        for cf in range(len(c2p[m_idx].coefs)):
            l_cf = c2p[m_idx].Ls[cf]
            m_cf = c2p[m_idx].Ms[cf]
            c2p[m_idx].coefs[cf] *= np.sqrt(ylm_norms2[l_cf][l_cf + m_cf])
            c2p[m_idx].coefs[cf] /= np.sqrt(cart_norm2s[m_idx])

    where_to_store = [np.zeros(len(c2p[bf.m].coefs), dtype=int) for bf in molden_in.basis]
    transfer_coef = [np.zeros(len(c2p[bf.m].coefs)) for bf in molden_in.basis]
    
    bf_list = []
    
    for bf_idx, bf in enumerate(molden_in.basis):
        c2p_entry = c2p[bf.m]
        
        for i in range(len(c2p_entry.coefs)):
            L = c2p_entry.Ls[i]
            addit_r_pwr = (bf.l - L) + bf.additional_r_power
            
            indx = _find_bs_with_props(bf_list, bf, L, c2p_entry.Ms[i], addit_r_pwr)
            _append = (indx == -1)
            
            if addit_r_pwr > 0 and ignore_lower_l:
                _append = False
                
            if _append:
                new_bf = deepcopy(bf)
                new_bf.l = L
                new_bf.m = c2p_entry.Ms[i]
                new_bf.additional_r_power = addit_r_pwr
                bf_list.append(new_bf)
                indx = _find_bs_with_props(bf_list, bf, L, c2p_entry.Ms[i], addit_r_pwr)
                if indx == -1:
                    raise RuntimeError("ERROR in creating basis function list: new function not found!")
                    
            where_to_store[bf_idx][i] = indx
            transfer_coef[bf_idx][i] = c2p_entry.coefs[i]

    n_basis_new = sum(2 * bf.l + 1 for bf in bf_list)
    molden_out.basis = bf_list
    molden_out.is_spherical = True
    
    molden_out.radial_parts = []
    for bf in bf_list:
        rp = RadialPart(
            exponents=bf.exponents.copy(),
            coefs=bf.coefs.copy(),
            center_id=bf.center_id,
            l_used_with=bf.l,
            addit_r_power=bf.additional_r_power
        )
        molden_out.radial_parts.append(rp)
        bf.radial_part_id = len(molden_out.radial_parts) - 1
        
    for mo_idx, mo in enumerate(molden_out.mos):
        new_bs_coefs = np.zeros(n_basis_new)
        sum_ignored = 0.0
        for cf in range(len(molden_in.basis)):
            for i in range(len(where_to_store[cf])):
                target_coef = where_to_store[cf][i]
                adduct = transfer_coef[cf][i] * mo.bs_coefs[cf]
                if target_coef != -1:
                    new_bs_coefs[target_coef] += adduct
                elif ignore_lower_l:
                    sum_ignored += adduct
        if abs(sum_ignored) > 1e-5:
            print(f"Warning: sum of ignored terms for MO {mo_idx+1} seems to be non-negligible (|sum| = {abs(sum_ignored):.3E})")
        mo.bs_coefs = new_bs_coefs
        
    return molden_out