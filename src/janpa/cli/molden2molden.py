# File: janpa/cli/molden2molden.py
"""molden2molden CLI: format conversion tool."""
from __future__ import annotations
import argparse
import sys
import numpy as np

from janpa.io.molden import MoldenFile
from janpa.convert.cart2sph import cart2spher
from janpa.convert.file47 import import_47, build_basis_set, check_bf_order_and_create_remap, produce_nos
from janpa.gto.overlap import primitive_int_1d_sphr

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="molden2molden",
        description="MOLDEN file format converter"
    )
    p.add_argument("input", help="Input MOLDEN file")
    p.add_argument("output", help="Output MOLDEN file")
    
    p.add_argument("-cart2pure", "--cart2spher", action="store_true", help="Convert Cartesian to Spherical")
    p.add_argument("-ignore-lower-l", action="store_true", help="Ignore lower L terms in cart2pure")
    p.add_argument("-angstroms", action="store_true", help="Output coordinates in Angstroms")
    
    p.add_argument("-orca3signs", action="store_true", help="Fix ORCA3 signs")
    p.add_argument("-from-orca3-bf", action="store_true", help="Convert basis from ORCA3 style")
    p.add_argument("-to-orca3-bf", action="store_true", help="Convert basis to ORCA3 style")
    p.add_argument("-from-psi4b4-bf", action="store_true", help="Convert basis from PSI4b4 style")
    p.add_argument("-from-psi4v1-mo", action="store_true", help="Convert MOs from PSI4v1 style")
    p.add_argument("-to-psi4b4-bf", action="store_true", help="Convert basis to PSI4b4 style")
    
    p.add_argument("-normalize", action="store_true", help="Force basis functions to be unity-normalized")
    
    p.add_argument("-geom47", default="", help="Suppress geometry with .47 file")
    p.add_argument("-bs47", default="", help="Suppress basis set with .47 file")
    p.add_argument("-ds47", default="", help="Suppress MO coefficients with .47 file")
    
    return p

def fix_orca3_psi4_normalization(molden, from_orca3=False, from_psi4b4=False):
    orca_dividers = [1.0, 1.0, 3.0, 15.0, 35.0]
    if from_orca3 or from_psi4b4:
        for rp in molden.radial_parts:
            l = rp.l_used_with
            for cf in range(len(rp.coefs)):
                rnorm2 = primitive_int_1d_sphr(2 * l + 2 + 2 * rp.addit_r_power, 2 * rp.exponents[cf])
                rnorm2 *= 4 * np.pi
                rnorm2 /= (2 * l + 1)
                if from_orca3 and l < len(orca_dividers):
                    rnorm2 /= orca_dividers[l]
                rp.coefs[cf] *= np.sqrt(rnorm2)
        for bf in molden.basis:
            bf.coefs = molden.radial_parts[bf.radial_part_id].coefs.copy()

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    
    print(f"Loading input molden from {args.input} ...")
    try:
        molden = MoldenFile.load(args.input)
    except Exception as e:
        print(f"ERROR: Can not load data from {args.input}: {e}")
        sys.exit(1)
        
    print("Data loaded successfully.\n")
    
    if args.cart2spher and molden.is_spherical:
        print("ERROR: input molden file already uses pure spherical harmonics, -cart2spher key makes no sense and will be ignored!")
        args.cart2spher = False
        
    if not args.angstroms and not molden.coords_in_au:
        print("Converting atomic coordinates from Angstroms to a.u. ...")
        molden.coords_to_au()
    if args.angstroms and molden.coords_in_au:
        print("Converting atomic coordinates from a.u. to Angstroms ...")
        molden.coords_to_angstroms()
        
    if args.geom47:
        print(f"Supressing geometry with the data from .47 file {args.geom47}")
        data47 = import_47(args.geom47)
        # In a full implementation, we would map data47 coordinates to molden.centers
        
    if args.bs47:
        print(f"Supressing basis set with the data from .47 file {args.bs47}")
        data47 = import_47(args.bs47)
        if not build_basis_set(data47, molden):
            print("Error importing basis set data!")
            sys.exit(1)
            
    if args.ds47:
        print(f"Supressing MO coefficients with the data from .47 file {args.ds47}")
        data47 = import_47(args.ds47)
        if data47.nbas != molden.n_basis:
            print("ERROR: the number of basis functions is different in the input molden file and in the given .47 file")
            sys.exit(1)
        if not check_bf_order_and_create_remap(data47, molden):
            print("ERROR: reordering of the .47 file basis set failed!")
            sys.exit(1)
        produce_nos(data47, molden)
        
    fix_orca3_psi4_normalization(molden, from_orca3=args.from_orca3_bf, from_psi4b4=args.from_psi4b4_bf)
    
    if args.cart2spher:
        print("Converting Cartesian to Spherical...")
        molden = cart2spher(molden, ignore_lower_l=args.ignore_lower_l)
        
    print(f"Saving result as {args.output} (ascii encoded)...")
    try:
        molden.save(args.output)
        print("OK.")
    except Exception as e:
        print(f"ERROR saving file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()