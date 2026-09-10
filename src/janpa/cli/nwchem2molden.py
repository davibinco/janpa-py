# File: janpa/cli/nwchem2molden.py
"""nwchem2molden CLI: NWChem output to MOLDEN converter."""
from __future__ import annotations
import argparse
import re
import sys
import numpy as np

from janpa.io.molden import MoldenFile
from janpa.gto.basis import AtomicCenter, MolecularOrbital

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="nwchem2molden",
        description="NWChem output to MOLDEN converter"
    )
    p.add_argument("logfile", help="NWChem output (.out) file")
    p.add_argument("ascfile", nargs="?", help="NWChem orbitals (.asc) file")
    p.add_argument("-o", "--output", default="nwchem.molden", help="Output MOLDEN file")
    return p

def _matches_fmts(lines, start_idx, fmts):
    result = []
    idx = start_idx
    for fmt in fmts:
        if idx >= len(lines): return None, idx
        s = lines[idx]
        if re.match(fmt, s):
            result.append(s)
            idx += 1
        else:
            return None, idx
    return result, idx

def parse_out_file(fname: str, molden: MoldenFile) -> bool:
    with open(fname, 'r') as f:
        lines = f.readlines()
        
    molden.is_spherical = False
    atoms = []
    basis_title = ""
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        if re.search(r'Basis\s+"([^"]+)"\s+->\s+"([^"]+)"\s+\(cartesian\)', line):
            molden.is_spherical = False
            basis_title = re.search(r'Basis\s+"([^"]+)"', line).group(1)
            print(f"CARTESIAN basis set information found: {basis_title}")
        elif re.search(r'Basis\s+"([^"]+)"\s+->\s+"([^"]+)"\s+\(spherical\)', line):
            molden.is_spherical = True
            basis_title = re.search(r'Basis\s+"([^"]+)"', line).group(1)
            print(f"SPHERICAL basis set information found: {basis_title}")
            
        geom_match = re.search(r'Geometry\s+"[^"]+"\s+->\s+"([^"]*)"', line)
        if geom_match:
            fmts = [r'\s+\-+', r'', r'\s+Output coordinates in.+', r'', r'\s+No.+\s+Tag.+\s+Charge.+\s+X.+\s+Y.+\s+Z', r'\s+\-+\s+\-+\s+\-+\s+\-+\s+\-+\s+\-+']
            res, next_i = _matches_fmts(lines, i+1, fmts)
            if res:
                print("Geometry information found")
                units_line = res[2]
                coords_to_au = 1.0
                scale_match = re.search(r'\(scale by\s+([0-9\.E\-\+]+)\s+to convert to a\.u\.\)', units_line, re.IGNORECASE)
                if scale_match:
                    try: coords_to_au = float(scale_match.group(1))
                    except: pass
                        
                i = next_i
                atoms = []
                while i < len(lines) and lines[i].strip():
                    props = lines[i].split()
                    if len(props) >= 7:
                        try:
                            cid = int(props[1])
                            name = props[2]
                            z = float(props[3])
                            x, y, z_coord = float(props[4]), float(props[5]), float(props[6])
                            atoms.append(AtomicCenter(name=name, id=cid, z=z, r0=np.array([x, y, z_coord]) * coords_to_au))
                        except ValueError:
                            break
                    i += 1
                molden.centers = atoms
                molden.coords_in_au = True
                continue
                
        i += 1
        
    if not atoms:
        print("ERROR: no information about molecular geometry has been found!")
        return False
        
    molden.title = basis_title
    return True

def parse_asc_file(fname: str, molden: MoldenFile) -> bool:
    with open(fname, 'r') as f:
        lines = f.readlines()
        
    print(f"Reading orbitals from {fname}")
    i = 4
    if i < len(lines): print(f'"scftype20": {lines[i].strip()}'); i += 1
    if i < len(lines): print(f'Date: {lines[i].strip()}'); i += 1
    if i < len(lines): print(f'job type: {lines[i].strip()}'); i += 1
    i += 3
    if i < len(lines): 
        bs_title = lines[i].strip()
        print(f"basis set name: {bs_title}")
        i += 1
        
    nsets = int(lines[i].strip()); i += 1
    if nsets != 1:
        print("ERROR: asc files with more than one set of orbitals are not supported!")
        return False
        
    nBF = int(lines[i].strip()); i += 1
    nMO = int(lines[i].strip()); i += 1
    print(f"There are {nBF} basis functions and {nMO} orbitals")
    
    def read_doubles(n):
        nonlocal i
        res = []
        while len(res) < n and i < len(lines):
            vals = lines[i].split()
            for v in vals: res.append(float(v))
            i += 1
        return res
        
    occupations = read_doubles(nBF)
    energies = read_doubles(nBF)
    
    all_positive = all(e > 0 for e in energies)
    
    molden.mos = []
    for mo_idx in range(nMO):
        bs_coefs = read_doubles(nBF)
        ene = 0.0 if all_positive else energies[mo_idx]
        molden.mos.append(MolecularOrbital(
            energy=ene,
            occupancy=occupations[mo_idx],
            spin=1,
            bs_coefs=np.array(bs_coefs)
        ))
        
    bf = 0
    ylm_norms_cart = [1.0] * 35 
    
    for rp in molden.radial_parts:
        L = rp.l_used_with
        if molden.is_spherical:
            if L == 2:
                for mo in molden.mos:
                    c = mo.bs_coefs
                    molden_fns = [c[bf+2], -c[bf+3], c[bf+1], c[bf+4], c[bf+0]]
                    c[bf:bf+5] = molden_fns
            elif L == 3:
                for mo in molden.mos:
                    c = mo.bs_coefs
                    molden_fns = [c[bf+3], -c[bf+4], c[bf+2], c[bf+5], c[bf+1], -c[bf+6], c[bf+0]]
                    c[bf:bf+7] = molden_fns
            elif L == 4:
                for mo in molden.mos:
                    c = mo.bs_coefs
                    molden_fns = [c[bf+4], -c[bf+5], c[bf+3], c[bf+6], c[bf+2], -c[bf+7], c[bf+1], c[bf+8], c[bf+0]]
                    c[bf:bf+9] = molden_fns
            bf += (2 * L + 1)
        else:
            if L == 0: bf += 1
            elif L == 1: bf += 3
            elif L == 2:
                for mo in molden.mos:
                    c = mo.bs_coefs
                    molden_fns = [c[bf+0], c[bf+3], c[bf+5], c[bf+1]*np.sqrt(ylm_norms_cart[7]*5), c[bf+2]*np.sqrt(ylm_norms_cart[8]*5), c[bf+4]*np.sqrt(ylm_norms_cart[9]*5)]
                    c[bf:bf+6] = molden_fns
                bf += 6
            elif L == 3:
                for mo in molden.mos:
                    c = mo.bs_coefs
                    molden_fns = [c[bf+0], c[bf+6], c[bf+9], c[bf+3], c[bf+1], c[bf+2], c[bf+5], c[bf+8], c[bf+7], c[bf+4]]
                    for k in range(10): molden_fns[k] *= np.sqrt(ylm_norms_cart[10+k]*7)
                    c[bf:bf+10] = molden_fns
                bf += 10
            elif L == 4:
                for mo in molden.mos:
                    c = mo.bs_coefs
                    molden_fns = [c[bf+0], c[bf+10], c[bf+14], c[bf+1], c[bf+2], c[bf+6], c[bf+11], c[bf+9], c[bf+13], c[bf+3], c[bf+5], c[bf+12], c[bf+4], c[bf+7], c[bf+8]]
                    for k in range(15): molden_fns[k] *= np.sqrt(ylm_norms_cart[20+k]*9)
                    c[bf:bf+15] = molden_fns
                bf += 15
                
    return True

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    
    molden = MoldenFile()
    if not parse_out_file(args.logfile, molden):
        sys.exit(1)
        
    if args.ascfile:
        if not parse_asc_file(args.ascfile, molden):
            sys.exit(1)
            
    molden.save(args.output)
    print(f"Saved to {args.output}")

if __name__ == "__main__":
    main()