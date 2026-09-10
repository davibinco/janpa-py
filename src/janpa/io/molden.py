"""MOLDEN file reader/writer. Replaces MOLDEN_IO.java.

Supports:
- [Molden Format], [Title], [Atoms], [GTO], [MO], [5D], [9G] sections
- Cartesian vs spherical detection
- Bohr ↔ Angstrom conversion
- Binary .bmolden format (optional, can be deferred)
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from janpa.gto.basis import AtomicCenter, BasisFunction, RadialPart, MolecularOrbital
from janpa.gto.overlap import primitive_int_1d_sphr
from janpa.utils.constants import BOHR_RADIUS, SPDF_LABELS


@dataclass
class MoldenFile:
    """In-memory representation of a MOLDEN file."""
    title: str = ""
    is_spherical: bool = False
    coords_in_au: bool = False
    centers: list[AtomicCenter] = field(default_factory=list)
    basis: list[BasisFunction] = field(default_factory=list)
    radial_parts: list[RadialPart] = field(default_factory=list)
    mos: list[MolecularOrbital] = field(default_factory=list)
    highest_l: int = 0
    allow_additional_r_power: bool = False

    @property
    def n_basis(self) -> int:
        return len(self.basis)

    @property
    def n_atoms(self) -> int:
        return len(self.centers)

    def coords_to_au(self) -> None:
        if self.coords_in_au:
            return
        for c in self.centers:
            c.r0 = c.r0 / BOHR_RADIUS
        self.coords_in_au = True
        self._sync_basis_coords()

    def coords_to_angstroms(self) -> None:
        if not self.coords_in_au:
            return
        for c in self.centers:
            c.r0 = c.r0 * BOHR_RADIUS
        self.coords_in_au = False
        self._sync_basis_coords()

    def _sync_basis_coords(self) -> None:
        for bf in self.basis:
            if 1 <= bf.center_id <= len(self.centers):
                bf.r0 = self.centers[bf.center_id - 1].r0

    def to_unnormalized_primitive_coefs(self) -> None:
        """Divide contraction coefs by primitive Gaussian norms."""
        for rp in self.radial_parts:
            l = rp.l_used_with
            for cf in range(len(rp.coefs)):
                if rp.addit_r_power != 0:
                    norm2 = primitive_int_1d_sphr(2 * l + 2 + 2 * rp.addit_r_power, 2 * rp.exponents[cf])
                else:
                    norm2 = primitive_int_1d_sphr(2 * l + 2, 2 * rp.exponents[cf])
                rp.coefs[cf] /= np.sqrt(norm2)

    def to_normalized_primitive_coefs(self) -> None:
        """Multiply contraction coefs by primitive Gaussian norms."""
        for rp in self.radial_parts:
            l = rp.l_used_with
            for cf in range(len(rp.coefs)):
                if rp.addit_r_power != 0:
                    norm2 = primitive_int_1d_sphr(2 * l + 2 + 2 * rp.addit_r_power, 2 * rp.exponents[cf])
                else:
                    norm2 = primitive_int_1d_sphr(2 * l + 2, 2 * rp.exponents[cf])
                rp.coefs[cf] *= np.sqrt(norm2)

    def un_normalize_primitives(self, ylm_norms2_over4pi: np.ndarray) -> None:
        if not self.is_spherical:
            raise ValueError("UnnormalizePrimitives() currently works with spherical basis sets only!")
        
        for rp in self.radial_parts:
            l = rp.l_used_with
            for cf in range(len(rp.coefs)):
                if rp.addit_r_power != 0:
                    rnorm2 = primitive_int_1d_sphr(2 * l + 2 + 2 * rp.addit_r_power, 2 * rp.exponents[cf])
                else:
                    rnorm2 = primitive_int_1d_sphr(2 * l + 2, 2 * rp.exponents[cf])
                rnorm2 *= 4 * np.pi
                rp.coefs[cf] /= np.sqrt(rnorm2 * ylm_norms2_over4pi[l, 0])
                
        for bf in self.basis:
            l = bf.l
            for cf in range(len(bf.coefs)):
                rnorm2 = primitive_int_1d_sphr(2 * l + 2 + 2 * bf.additional_r_power, 2 * bf.exponents[cf])
                rnorm2 *= 4 * np.pi
                bf.coefs[cf] /= np.sqrt(rnorm2 * ylm_norms2_over4pi[l, l + bf.m])

    @classmethod
    def load(cls, path: str | Path, allow_additional_r_power: bool = True) -> MoldenFile:
        """Parse a MOLDEN-format file."""
        path = Path(path)
        with open(path, 'r') as f:
            lines = f.readlines()
            
        if not lines or not "[MOLDEN FORMAT]" in lines[0].upper():
            raise ValueError("Not a valid MOLDEN file: missing [Molden Format]")
            
        mf = cls(allow_additional_r_power=allow_additional_r_power)
        
        section = None
        title_lines = []
        
        use_5d = False
        use_7f = False
        use_9g = False
        coords_in_au = False
        
        centers_dict = {}
        radial_parts = []
        mos = []
        
        i = 1
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
                
            if line.startswith('['):
                match = re.match(r"^\[([^\]]+)\](.*)$", line)
                if match:
                    sec_name = match.group(1).strip().upper()
                    extra = match.group(2).strip().upper()
                    
                    if sec_name == "TITLE":
                        section = "TITLE"
                        title_lines = []
                    elif sec_name.startswith("ATOMS"):
                        section = "ATOMS"
                        if "AU" in extra:
                            coords_in_au = True
                    elif sec_name == "GTO":
                        section = "GTO"
                    elif sec_name == "MO":
                        section = "MO"
                    elif sec_name == "5D" or sec_name == "5D7F":
                        use_5d = True
                        use_7f = True
                        section = None
                    elif sec_name == "9G":
                        use_9g = True
                        section = None
                    else:
                        section = None
                    i += 1
                    continue
                
            if section == "TITLE":
                title_lines.append(line)
                i += 1
                continue
                
            if section == "ATOMS":
                if line.startswith('['):
                    section = None
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    name = parts[0]
                    cid = int(parts[1])
                    z = float(parts[2])
                    r0 = np.array([float(parts[3]), float(parts[4]), float(parts[5])])
                    centers_dict[cid] = AtomicCenter(name=name, id=cid, z=z, r0=r0)
                i += 1
                continue
                
            if section == "GTO":
                if line.startswith('['):
                    section = None
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        current_center_id = int(parts[0])
                    except ValueError:
                        i += 1
                        continue
                    
                    i += 1
                    while i < len(lines):
                        shell_line = lines[i].strip()
                        if not shell_line or shell_line.startswith('['):
                            break
                        shell_parts = shell_line.split()
                        if len(shell_parts) >= 2:
                            shell_type = shell_parts[0].lower()
                            add_r = 0
                            match = re.match(r"^([a-z]+)([0-9]+)$", shell_type)
                            if match:
                                if not allow_additional_r_power:
                                    raise ValueError(f"Extended molden command found ({shell_type}) but allow_additional_r_power is false.")
                                add_r = int(match.group(2))
                                shell_type = match.group(1)
                                
                            if shell_type == "sp":
                                raise ValueError("SP shells are not supported!")
                                
                            if shell_type == "s": l = 0
                            elif shell_type == "p": l = 1
                            elif shell_type == "d": l = 2
                            elif shell_type == "f": l = 3
                            elif shell_type == "g": l = 4
                            elif shell_type == "h": l = 5
                            else:
                                i += 1
                                continue
                                
                            n_coef = int(shell_parts[1])
                            exponents = []
                            coefs = []
                            i += 1
                            for _ in range(n_coef):
                                if i >= len(lines): break
                                coef_line = lines[i].strip().upper().replace('D', 'E')
                                cp = coef_line.split()
                                if len(cp) >= 2:
                                    exponents.append(float(cp[0]))
                                    coefs.append(float(cp[1]))
                                i += 1
                                
                            radial_parts.append({
                                'center_id': current_center_id,
                                'l': l,
                                'addit_r_power': add_r,
                                'exponents': np.array(exponents),
                                'coefs': np.array(coefs)
                            })
                        else:
                            i += 1
                    continue
                i += 1
                continue
                
            if section == "MO":
                if line.startswith('['):
                    section = None
                    continue
                
                energy = 0.0
                occupancy = 0.0
                spin = 0
                
                while i < len(lines):
                    prop_line = lines[i].strip()
                    if not prop_line or prop_line.startswith('[') or '=' not in prop_line:
                        break
                    if '=' in prop_line:
                        k, v = prop_line.split('=', 1)
                        k = k.strip().lower()
                        v = v.strip()
                        if k == "ene":
                            energy = float(v)
                        elif k == "occup":
                            occupancy = float(v)
                        elif k == "spin":
                            if "alpha" in v.lower():
                                spin = 1
                            elif "beta" in v.lower():
                                spin = -1
                    i += 1

                bs_coefs_dict = {}
                while i < len(lines):
                    coef_line = lines[i].strip()
                    if not coef_line:
                        i += 1
                        continue
                    if coef_line.startswith('['):
                        break
                    
                    if '=' in coef_line:
                        break
                    
                    cp = coef_line.split()
                    if len(cp) >= 2:
                        try:
                            idx = int(cp[0]) - 1
                            val = float(cp[1])
                            bs_coefs_dict[idx] = val
                        except ValueError:
                            # Fallback break if the line doesn't start with an integer
                            break
                    i += 1
                    
                mos.append({
                    'energy': energy,
                    'occupancy': occupancy,
                    'spin': spin,
                    'bs_coefs_dict': bs_coefs_dict
                })
                continue
                
            i += 1
            
        max_id = max(centers_dict.keys()) if centers_dict else 0
        mf.centers = [centers_dict[k] for k in sorted(centers_dict.keys())]
        mf.title = "\n".join(title_lines).strip()
        mf.coords_in_au = coords_in_au
        
        mf.radial_parts = []
        raw_basis = []
        mf.highest_l = 0
        
        cartesian_counts = [1, 3, 6, 10, 15, 21]
        first_cart_m = [0, 1, 4, 10, 20, 35]
        
        for rp_idx, rp in enumerate(radial_parts):
            mf.radial_parts.append(RadialPart(
                exponents=rp['exponents'],
                coefs=rp['coefs'],
                center_id=rp['center_id'],
                l_used_with=rp['l'],
                addit_r_power=rp['addit_r_power']
            ))
            if rp['l'] > mf.highest_l:
                mf.highest_l = rp['l']
                
            center = mf.centers[rp['center_id'] - 1]
            base_bf = BasisFunction(
                l=rp['l'],
                m=0,
                center_id=rp['center_id'],
                r0=center.r0,
                exponents=rp['exponents'],
                coefs=rp['coefs'].copy(),
                radial_part_id=rp_idx,
                additional_r_power=rp['addit_r_power']
            )
            raw_basis.append(base_bf)
            
        is_spherical = use_5d and use_7f
        if mf.highest_l >= 4:
            is_spherical = is_spherical and use_9g
        mf.is_spherical = is_spherical
        
        mf.basis = []
        for base_bf in raw_basis:
            l = base_bf.l
            if is_spherical:
                bf0 = BasisFunction(
                    l=l, m=0, center_id=base_bf.center_id, r0=base_bf.r0,
                    exponents=base_bf.exponents, coefs=base_bf.coefs.copy(),
                    radial_part_id=base_bf.radial_part_id, additional_r_power=base_bf.additional_r_power
                )
                mf.basis.append(bf0)
                for ii in range(1, l + 1):
                    bf_plus = BasisFunction(
                        l=l, m=ii, center_id=base_bf.center_id, r0=base_bf.r0,
                        exponents=base_bf.exponents, coefs=base_bf.coefs.copy(),
                        radial_part_id=base_bf.radial_part_id, additional_r_power=base_bf.additional_r_power
                    )
                    mf.basis.append(bf_plus)
                    bf_minus = BasisFunction(
                        l=l, m=-ii, center_id=base_bf.center_id, r0=base_bf.r0,
                        exponents=base_bf.exponents, coefs=base_bf.coefs.copy(),
                        radial_part_id=base_bf.radial_part_id, additional_r_power=base_bf.additional_r_power
                    )
                    mf.basis.append(bf_minus)
            else:
                start_m = first_cart_m[l] if l < len(first_cart_m) else 0
                count = cartesian_counts[l] if l < len(cartesian_counts) else 1
                for m_idx in range(count):
                    bf_cart = BasisFunction(
                        l=l, m=start_m + m_idx, center_id=base_bf.center_id, r0=base_bf.r0,
                        exponents=base_bf.exponents, coefs=base_bf.coefs.copy(),
                        radial_part_id=base_bf.radial_part_id, additional_r_power=base_bf.additional_r_power
                    )
                    mf.basis.append(bf_cart)
                    
        n_basis = len(mf.basis)
        for mo_dict in mos:
            bs_coefs = np.zeros(n_basis)
            for idx, val in mo_dict['bs_coefs_dict'].items():
                if 0 <= idx < n_basis:
                    bs_coefs[idx] = val
            mf.mos.append(MolecularOrbital(
                energy=mo_dict['energy'],
                occupancy=mo_dict['occupancy'],
                spin=mo_dict['spin'],
                bs_coefs=bs_coefs
            ))
            
        return mf

    def save(self, path: str | Path) -> None:
        """Write MOLDEN-format file."""
        path = Path(path)
        with open(path, 'w') as f:
            f.write("[Molden Format]\n")
            f.write("[Title]\n")
            f.write(f"{self.title}\n\n")
            
            if self.coords_in_au:
                f.write("[Atoms] AU\n")
            else:
                f.write("[Atoms] Angs\n")
                
            for idx, c in enumerate(self.centers):
                f.write(f"{c.name:<2}{idx+1:4d} {c.z:3.0f} {c.r0[0]:20.10f} {c.r0[1]:20.10f} {c.r0[2]:20.10f}\n")
                
            f.write("[GTO]\n")
            for idx, c in enumerate(self.centers):
                cid = idx + 1
                f.write(f"{cid:3d} 0\n")
                for rp in self.radial_parts:
                    if rp.center_id == cid:
                        shell_char = SPDF_LABELS[rp.l_used_with] if rp.l_used_with < len(SPDF_LABELS) else 'h'
                        if rp.addit_r_power == 0:
                            f.write(f"{shell_char}{rp.n_prim:4d} 1.0\n")
                        else:
                            if not self.allow_additional_r_power:
                                raise ValueError(f"Radial part attempts using additional_r_power>0 but allow_additional_r_power is false!")
                            f.write(f"{shell_char}{rp.addit_r_power}{rp.n_prim:4d} 1.0\n")
                        for e, coef in zip(rp.exponents, rp.coefs):
                            f.write(f"{e:20.10f} {coef:20.10f}\n")
                f.write("\n")
                
            if self.is_spherical:
                f.write("[5D]\n")
                f.write("[9G]\n")
                
            if self.mos:
                f.write("[MO]\n")
                for mo in self.mos:
                    f.write(" Sym= 1a\n")
                    f.write(f" Ene={mo.energy:22.14E}\n")
                    if mo.spin == 1:
                        f.write(" Spin= Alpha\n")
                    else:
                        f.write(" Spin= Beta\n")
                    f.write(f" Occup={mo.occupancy:9.6f}\n")
                    for cf_idx, val in enumerate(mo.bs_coefs):
                        f.write(f"{cf_idx+1:3d} {val:20.12f}\n")

    def get_mo_coef_matrix(self) -> np.ndarray:
        """Return (n_mo, n_basis) matrix of MO coefficients."""
        if not self.mos:
            return np.array([])
        return np.array([mo.bs_coefs for mo in self.mos])