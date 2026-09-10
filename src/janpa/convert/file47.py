# File: janpa/convert/file47.py
"""NBO .47 file import."""
from __future__ import annotations
import re
import numpy as np
from dataclasses import dataclass, field
import logging

logger = logging.getLogger("janpa")

@dataclass
class File47Data:
    natoms: int = 0
    nbas: int = 0
    bodm: bool = False
    upper: bool = False
    bohr: bool = False
    nshell: int = 0
    nexp: int = 0
    density: np.ndarray | None = None
    overlap: np.ndarray | None = None
    basis_center: np.ndarray | None = None
    basis_label: np.ndarray | None = None
    ncomp: np.ndarray | None = None
    nprim: np.ndarray | None = None
    nptr: np.ndarray | None = None
    c_exp: np.ndarray | None = None
    c_CS: np.ndarray | None = None
    c_CP: np.ndarray | None = None
    c_CD: np.ndarray | None = None
    c_CF: np.ndarray | None = None
    basis_remap: np.ndarray | None = None

def import_47(fname: str) -> File47Data:
    """Parse an NBO3-compatible .47 file."""
    data = File47Data()
    
    with open(fname, 'r') as f:
        lines = f.readlines()
        
    section = None
    ind = {
        'D_i': 0, 'D_j': 0,
        'S_i': 0, 'S_j': 0,
        'bas_cntr': 0, 'bas_lbl': 0,
        'ncomp': 0, 'nprim': 0, 'nptr': 0,
        'exp': 0, 'CS': 0, 'CP': 0, 'CD': 0, 'CF': 0
    }
    
    basis_intarr = None
    basis_intarr_IND = None
    contract_IND = None
    
    def add_to_matrix(matrix, index_I, index_J, value):
        i = ind[index_I]
        j = ind[index_J]
        matrix[i, j] = value
        if data.upper:
            matrix[j, i] = value
        if (data.upper and j == i) or (j == data.nbas - 1):
            ind[index_I] += 1
            ind[index_J] = 0
        else:
            ind[index_J] += 1
            
    def add_to_int_arr(arr, index_I, value):
        i = ind[index_I]
        arr[i] = value
        ind[index_I] += 1
        
    def add_to_dbl_arr(arr, index_I, value):
        i = ind[index_I]
        arr[i] = value
        ind[index_I] += 1
        
    nDens = 0
    nOvrlp = 0
    
    for line in lines:
        line = line.strip()
        if not line: continue
            
        tokens = re.split(r'[ =]+', line)
        
        i = 0
        while i < len(tokens):
            s = tokens[i]
            increment = 0
            
            if s == "$END":
                if section == "GENNBO":
                    data.density = np.zeros((data.nbas, data.nbas))
                    data.overlap = np.zeros((data.nbas, data.nbas))
                    data.basis_center = np.zeros(data.nbas, dtype=int)
                    data.basis_label = np.zeros(data.nbas, dtype=int)
                    data.basis_remap = np.zeros(data.nbas, dtype=int)
                    ind['D_i'] = 0; ind['D_j'] = 0
                    ind['S_i'] = 0; ind['S_j'] = 0
                    ind['bas_cntr'] = 0; ind['bas_lbl'] = 0
                section = None
                increment = 1
                
            elif section == "OVERLAP":
                add_to_matrix(data.overlap, 'S_i', 'S_j', float(s))
                increment = 1
                nOvrlp += 1
                
            elif section == "DENSITY":
                add_to_matrix(data.density, 'D_i', 'D_j', float(s))
                increment = 1
                nDens += 1
                
            elif section == "BASIS":
                if s.lower() == "center":
                    basis_intarr = data.basis_center
                    basis_intarr_IND = 'bas_cntr'
                    increment = 1
                elif s.lower() == "label":
                    basis_intarr = data.basis_label
                    basis_intarr_IND = 'bas_lbl'
                    increment = 1
                else:
                    if re.match(r"^[0-9]+$", s):
                        add_to_int_arr(basis_intarr, basis_intarr_IND, int(s))
                        increment = 1
                        
            elif section == "GENNBO":
                sl = s.lower()
                if sl == "bodm": data.bodm = True; increment = 1
                elif sl == "upper": data.upper = True; increment = 1
                elif sl == "bohr": data.bohr = True; increment = 1
                elif sl == "nbas":
                    data.nbas = int(tokens[i+1])
                    increment = 2
                elif sl == "natoms":
                    data.natoms = int(tokens[i+1])
                    increment = 2
                else:
                    if increment == 0:
                        logger.warning(f"UNKNOWN KEYWORD in GENNBO: {s}")
                        increment = 1
                        
            elif section == "CONTRACT":
                sl = s.lower()
                if sl in ("nshell", "nshells"):
                    data.nshell = int(tokens[i+1])
                    increment = 2
                    data.ncomp = np.zeros(data.nshell, dtype=int)
                    data.nprim = np.zeros(data.nshell, dtype=int)
                    data.nptr = np.zeros(data.nshell, dtype=int)
                    ind['ncomp'] = 0; ind['nprim'] = 0; ind['nptr'] = 0
                elif sl == "nexp":
                    data.nexp = int(tokens[i+1])
                    increment = 2
                    data.c_exp = np.zeros(data.nexp)
                    data.c_CS = np.zeros(data.nexp)
                    data.c_CP = np.zeros(data.nexp)
                    data.c_CD = np.zeros(data.nexp)
                    data.c_CF = np.zeros(data.nexp)
                    ind['exp'] = 0; ind['CS'] = 0; ind['CP'] = 0; ind['CD'] = 0; ind['CF'] = 0
                elif sl == "ncomp": contract_IND = 'ncomp'; increment = 1
                elif sl == "nprim": contract_IND = 'nprim'; increment = 1
                elif sl == "nptr": contract_IND = 'nptr'; increment = 1
                elif sl == "exp": contract_IND = 'exp'; increment = 1
                elif sl == "cs": contract_IND = 'CS'; increment = 1
                elif sl == "cp": contract_IND = 'CP'; increment = 1
                elif sl == "cd": contract_IND = 'CD'; increment = 1
                elif sl == "cf": contract_IND = 'CF'; increment = 1
                else:
                    if re.match(r"^[0-9\.DE\+\-]+$", s):
                        if contract_IND in ('exp', 'CS', 'CP', 'CD', 'CF'):
                            attr_name = 'c_exp' if contract_IND == 'exp' else f'c_{contract_IND}'
                            add_to_dbl_arr(getattr(data, attr_name), contract_IND, float(s))
                        else:
                            add_to_int_arr(getattr(data, contract_IND), contract_IND, int(s))
                        increment = 1
                        
            if increment == 0:
                if s == "$NBO": section = "NBO"; increment = 1
                elif s == "$GENNBO": section = "GENNBO"; increment = 1
                elif s == "$COORD": section = "COORD"; increment = 1
                elif s == "$BASIS": section = "BASIS"; increment = 1
                elif s == "$OVERLAP": section = "OVERLAP"; increment = 1
                elif s == "$DENSITY": section = "DENSITY"; increment = 1
                elif s == "$CONTRACT": section = "CONTRACT"; increment = 1
                else: increment = 1
                
            i += increment
            
    logger.info(f"nbas = {data.nbas}")
    logger.info(f"{nDens} elements of the density matrix read")
    logger.info(f"{nOvrlp} elements of the overlap matrix read")
    
    return data

def check_bf_order_and_create_remap(data: File47Data, molden) -> bool:
    ok = True
    for i in range(data.nbas):
        if molden.basis[i].center_id != data.basis_center[i]:
            logger.error(f"ERROR: Basis function {i+1} center ID mistach: {data.basis_center[i]} in .47 file v.s. {molden.basis[i].center_id} required by molden file")
            ok = False
        if molden.basis[i].l != (data.basis_label[i] // 100):
            logger.error(f"ERROR: Basis function {i+1} center L mistach: {data.basis_label[i] // 100} in .47 file v.s. {molden.basis[i].l} required by molden file")
            ok = False
    if not ok:
        return False
        
    proper_ordering_SPH = [
        [51], [151, 152, 153], [255, 252, 253, 254, 251],
        [351, 352, 353, 354, 355, 356, 357], [451, 452, 453, 454, 455, 456, 457, 458, 459],
        [551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561]
    ]
    proper_ordering_CART = [
        [51], [151, 152, 153], [201, 204, 206, 202, 203, 205],
        [301, 307, 310, 304, 302, 303, 306, 309, 308, 305],
        [401, 411, 415, 402, 403, 407, 412, 410, 414, 404, 406, 413, 405, 408, 409]
    ]
    cart_compon_count = [1, 3, 6, 10, 15]
    
    i = 0
    while i < data.nbas:
        L = molden.basis[i].l
        if molden.is_spherical:
            component_count = 2 * L + 1
            ordering = proper_ordering_SPH
        else:
            component_count = cart_compon_count[L] if L < len(cart_compon_count) else 1
            ordering = proper_ordering_CART
            
        for M in range(component_count):
            label_to_find = data.basis_label[i+M]
            if label_to_find < 154:
                label_to_find = (label_to_find // 100) * 100 + 50 + (label_to_find % 10)
                
            found = False
            offset = 0
            if L < len(ordering):
                while (not found) and (offset < len(ordering[L])):
                    found = (ordering[L][offset] == label_to_find)
                    if not found: offset += 1
            if not found:
                logger.error(f"ERROR: function label {data.basis_label[i+M]} is unknown for L = {L}")
                return False
                
            data.basis_remap[i+M] = i + offset
        i += component_count
    return True

def produce_nos(data: File47Data, molden):
    """Creating natural orbitals from the .47 file density/overlap matrices data..."""
    D = data.density
    S = data.overlap
    print(f" tr(S) = {np.trace(S):.10f}")
    print(f" tr(D) = {np.trace(D):.10f}")
    
    if not data.bodm:
        logger.error("ERROR: non-BODM density matrices are not supported yet!")
        return
        
    w, v = np.linalg.eigh(S)
    w = np.maximum(w, 0) 
    S05 = v @ np.diag(np.sqrt(w)) @ v.T
    
    DD = S05 @ D @ S05
    eig_vals, eig_vecs = np.linalg.eigh(DD)
    
    print(f"tr(EIG) = {np.sum(eig_vals)}")
    
    w_inv = np.zeros_like(w)
    mask = w > 1e-12
    w_inv[mask] = 1.0 / np.sqrt(w[mask])
    Sm05 = v @ np.diag(w_inv) @ v.T
    
    U = Sm05 @ eig_vecs
    
    from janpa.gto.basis import MolecularOrbital
    molden.mos = []
    for i in range(data.nbas):
        bs_coefs = np.zeros(data.nbas)
        for cf in range(data.nbas):
            bs_coefs[data.basis_remap[cf]] = U[cf, i]
            
        molden.mos.append(MolecularOrbital(
            energy=0.0,
            occupancy=eig_vals[i],
            spin=1,
            bs_coefs=bs_coefs
        ))

def build_basis_set(data: File47Data, molden) -> bool:
    nRadialParts = 0
    for sh in range(data.nshell):
        if data.ncomp[sh] != 4: nRadialParts += 1
        else: nRadialParts += 2
            
    from janpa.gto.basis import RadialPart
    molden.radial_parts = []
    for i in range(nRadialParts):
        molden.radial_parts.append(RadialPart(
            exponents=np.array([]), coefs=np.array([]), center_id=0, l_used_with=0, addit_r_power=0
        ))
        
    nRP_Components = np.zeros(nRadialParts, dtype=int)
    L_max = 0
    Cartesian_D_found = False
    Cartesian_F_found = False
    Cartesian_G_found = False
    
    rp = 0
    for sh in range(data.nshell):
        molden.radial_parts[rp].exponents = data.c_exp[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
        molden.radial_parts[rp].coefs = np.zeros(data.nprim[sh])
        L = -2
        
        nc = data.ncomp[sh]
        if nc == 1: L = 0
        elif nc == 3: L = 1
        elif nc == 4: L = -1
        elif nc in (5, 6): L = 2
        elif nc in (7, 10): L = 3
        elif nc in (9, 15): L = 4
        
        if L == -2:
            logger.error(f"Unrecognized ncomp value: ncomp[{sh}] = {nc}")
            return False
            
        if L == -1:
            molden.radial_parts[rp].exponents = data.c_exp[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
            molden.radial_parts[rp].coefs = data.c_CS[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
            molden.radial_parts[rp].l_used_with = 0
            nRP_Components[rp] = 1
            rp += 1
            L = 1
            
        if L == 0: molden.radial_parts[rp].coefs = data.c_CS[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
        elif L == 1: molden.radial_parts[rp].coefs = data.c_CP[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
        elif L == 2: molden.radial_parts[rp].coefs = data.c_CD[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
        elif L == 3: molden.radial_parts[rp].coefs = data.c_CF[data.nptr[sh]-1 : data.nptr[sh]-1 + data.nprim[sh]]
            
        molden.radial_parts[rp].l_used_with = L
        if nc != 4: nRP_Components[rp] = nc
        else: nRP_Components[rp] = 3
            
        if L == 2 and nRP_Components[rp] != 5: Cartesian_D_found = True
        if L == 3 and nRP_Components[rp] != 7: Cartesian_F_found = True
        if L == 4 and nRP_Components[rp] != 9: Cartesian_G_found = True
        
        if L > L_max: L_max = L
        rp += 1
        
    logger.info(f"{rp} radial parts have been assembled, L_MAX = {L_max}")
    molden.is_spherical = (not Cartesian_D_found) and (not Cartesian_F_found) and (not Cartesian_G_found)
    
    if not molden.is_spherical:
        if L_max >= 2 and not Cartesian_D_found:
            logger.error("Mixing of spherical D / cartesian others basis is not supported!")
            return False
        if L_max >= 3 and not Cartesian_F_found:
            logger.error("Mixing of spherical F / cartesian others basis is not supported!")
            return False
        if L_max >= 4 and not Cartesian_G_found:
            logger.error("Mixing of spherical G / cartesian others basis is not supported!")
            return False
            
    return True