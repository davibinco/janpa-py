# File: janpa/cli/main.py
"""Main CLI entry point for the janpa NPA program."""
from __future__ import annotations
import argparse
import logging
import sys
import time
import numpy as np

from janpa.io.molden import MoldenFile
from janpa.npa.npa import run_npa
from janpa.clpo.lpo import create_clpos, CLPOOptions
from janpa.utils.printout import print_stars, print_matrix, export_orbitals

logger = logging.getLogger("janpa")

def show_banner():
    print(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *")
    print(" * janpa: A cross-platform open-source implementation of NPA * ")
    print(" * and other electronic structure analysis methods with Python * ")
    print(" * A part of JANPA package, http://janpa.sourceforge.net * ")
    print(" * Version: 2.02 (13-01-2019) [Python port] * ")
    print(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *")
    print()
    print(" (c) Tymofii Nikolaienko, 2014-2019")
    print(" Python port (c) 2024")
    print()
    print(" If any results obtained with this program are published,")
    print(" or for any other reasons, please, cite this work as: ")
    print(" 1) T.Y.Nikolaienko, L.A.Bulavin; Int. J. Quantum Chem. (2019), ")
    print(" Vol.119, page e25798, DOI: 10.1002/qua.25798")
    print(" 2) T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem.(2014),")
    print(" Vol.1050, pages 15-22, DOI: 10.1016/j.comptc.2014.10.002")
    print_stars()

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="janpa",
        description="JANPA: Natural Population Analysis and electronic structure tools",
    )
    p.add_argument("-i", "--input", required=True, help="Input MOLDEN file")
    p.add_argument("-p", "--params", default="", help="Read options from parameter file")
    
    p.add_argument("--npacharges", default="", help="Export NPA charges to file")
    p.add_argument("--wiberg-file", default="", help="Export Wiberg bond orders")
    p.add_argument("--nao-molden", default="", help="Export NAOs in MOLDEN format")
    p.add_argument("--pnao-molden", default="", help="Export PNAOs in MOLDEN format")
    p.add_argument("--pnao2ao-file", default="", help="Export PNAO to AO matrix")
    
    p.add_argument("--lho-molden", default="", help="Export LHOs in MOLDEN format")
    p.add_argument("--clpo-molden", default="", help="Export CLPOs in MOLDEN format")
    p.add_argument("--aho-molden", default="", help="Export AHOs in MOLDEN format")
    p.add_argument("--lpo-molden", default="", help="Export LPOs in MOLDEN format")
    
    p.add_argument("--clpo2lho-file", default="", help="Export CLPO to LHO matrix")
    p.add_argument("--lho2nao-file", default="", help="Export LHO to NAO matrix")
    p.add_argument("--lpo2aho-file", default="", help="Export LPO to AHO matrix")
    p.add_argument("--aho2nao-file", default="", help="Export AHO to NAO matrix")

    p.add_argument("--do-fock", action="store_true", help="Perform Fock matrix analysis")
    p.add_argument("--fock-ao-file", default="", help="Export Fock matrix in AO basis")
    p.add_argument("--fock-nao-file", default="", help="Export Fock matrix in NAO basis")
    p.add_argument("--print-nbm-submatrices", action="store_true", help="Print NBM submatrices")
    
    p.add_argument("--print-geom", action="store_true", help="Print geometry")
    p.add_argument("--overlap-naive", action="store_true", help="Use naive overlap calculation")
    p.add_argument("--edges", default="", help="Custom bonding graph edges")
    p.add_argument("--hybr-opt-thresh", type=float, default=1e-5, help="Hybrid optimization threshold")
    p.add_argument("--hybr-opt-max-iter", type=int, default=10000, help="Max hybrid optimization iterations")
    p.add_argument("--max-bond-ionicity", type=float, default=0.90, help="Max CLPO bond ionicity")
    
    p.add_argument("--verbose", action="store_true", help="Verbose output")
    return p

def print_bond_indices(indices: np.ndarray, centers, export_fname: str):
    n_atoms = len(centers)
    print("Wiberg-Mayer bond indices (based on density matrix in NAO basis):")
    header = f"{'Centr. A/B':>10}"
    for i in range(n_atoms):
        header += f"{i+1:10d}"
    print(header)
    
    for i in range(n_atoms):
        row = f"{i+1:7d} "
        for j in range(i):
            row += f"{'':>10}"
        row += f" ({indices[i, i]:7.4f})"
        for j in range(i + 1, n_atoms):
            row += f"{indices[i, j]:10.4f}"
        print(row)
        
    if export_fname:
        atom_labels = [f"{c.name}{i+1}" for i, c in enumerate(centers)]
        print_matrix(
            indices,
            "Wiberg-Mayer bond indices computed in NAO basis (note: "
            "diagonal elements are atomic 'valencies' (sums of all bond orders in current line)):",
            atom_labels, atom_labels, export_fname
        )

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    
    show_banner()
    
    t0 = time.perf_counter()
    logger.info("JANPA starting...")
    
    fname = args.input
    
    print(f"Loading MOLDEN from {fname}")
    try:
        molden = MoldenFile.load(fname)
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)
        
    print()
    print("Data loaded successfully")
    print(f" Number of basis functions: {molden.n_basis}; number of molecular orbitals: {len(molden.mos)} ")
    print()
    
    if not molden.is_spherical:
        print(" ERROR: Input molden file does not have a 'pure' ('spherical') basis set!")
        print(" Try using molden2molden -cart2pure to convert it.")
        sys.exit(1)
        
    molden.coords_to_au()
    
    if args.print_geom:
        print("\nCartesian coordinates of the atoms (in atomic units):")
        print(f"{'ID':<4}{'Element':<8}{'Nucl.Chrg.':<12}{'X':<15}{'Y':<15}{'Z':<15}")
        for i, c in enumerate(molden.centers):
            print(f"{i+1:<4}{c.name:<8}{c.z:<12.1f}{c.r0[0]:<15.10f}{c.r0[1]:<15.10f}{c.r0[2]:<15.10f}")
        print_stars()
        
    if not args.overlap_naive:
        molden.to_unnormalized_primitive_coefs()
        
    npa_res = run_npa(molden, {"verbose": args.verbose})
    print_stars()
    
    if npa_res.wiberg_indices is not None:
        print_bond_indices(npa_res.wiberg_indices, molden.centers, args.wiberg_file)
        print_stars()
        
    print("=" * 100)
    
    clpo_opts = CLPOOptions(
        hybr_opt_conv_thresh=args.hybr_opt_thresh,
        hybr_opt_max_iter=args.hybr_opt_max_iter,
        max_bond_ionicity_threshold=args.max_bond_ionicity,
        edges=args.edges
    )
    
    clpo_res = create_clpos(npa_res.sds_nao, npa_res.nao, molden.centers, clpo_opts)
    print_stars()
    
    print("Atomic connectivity analysis based on CLPO bonding graph:")
    print()
    print_stars()
    
    if args.npacharges and npa_res.npa_charges.size > 0:
        print(f"Writing NPA charges to\t{args.npacharges}")
        with open(args.npacharges, 'w') as ch:
            for q in npa_res.npa_charges:
                ch.write(f"{q:.16f}\n")

    if args.do_fock and len(molden.mos) >= molden.n_basis:
        print("Fock matrix analysis in NAO basis")
        mo_energies = np.diag([mo.energy for mo in molden.mos])
        C = molden.get_mo_coef_matrix().T  # Shape: (n_basis, n_mo)
        try:
            C_inv = np.linalg.inv(C)
            print("Transforming Fock matrix to AO basis...")
            Fock_AO = C @ mo_energies @ C.T
            
            if hasattr(npa_res, 'ao_names'):
                print_matrix(Fock_AO, "Fock matrix in AO basis", npa_res.ao_names, npa_res.ao_names, args.fock_ao_file)
                
            if npa_res.nao_to_ao is not None:
                print("Transforming Fock matrix to NAO basis...")
                Fock_NAO = npa_res.nao_to_ao.T @ Fock_AO @ npa_res.nao_to_ao
                # Further NAO exports would go here
        except np.linalg.LinAlgError:
            print("ERROR: Could not invert MO -> AO matrix for Fock analysis.")
                
    if args.nao_molden and getattr(npa_res, 'nao_to_ao', None) is not None:
        print(f"Writing NAOs to\t{args.nao_molden}")
        # export_orbitals(molden, npa_res.nao_to_ao, ...)
        
    elapsed = time.perf_counter() - t0
    logger.info("Total run time: %.2f s", elapsed)
    print_stars()

if __name__ == "__main__":
    main()