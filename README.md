# janpa-py

A fast, pure-Python port of [JANPA](https://janpa.sourceforge.net/) for Natural Population Analysis (NPA) and localized orbitals (CLPO/AHO). This library provides native Python implementations of the 7-step NPA algorithm and Chemist's Localized Property-optimized Orbitals (CLPO) using `numpy` and `scipy`, completely eliminating the need for Java binaries or temporary disk I/O.

## Key Features

- **Pure Python:** Built entirely on standard scientific Python libraries.
- **In-Memory Execution:** No messy `.molden` or `graph` text files left on your hard drive.
- **Seamless PySCF Integration:** Designed to work directly with PySCF molecule objects and mean-field results.
- **Chemical Graph Extraction:** Automatically extracts bonding graphs (edges) from the CLPO density matrix for downstream quantum chemistry workflows.

## Installation

Clone the repository and install it into your Python environment:

    git clone https://github.com/yourusername/janpa.git
    cd janpa
    pip install .

## Quick Start with PySCF

Here is a minimal example demonstrating how to run NPA, generate CLPO orbitals, and extract the chemical bonding graph directly in memory using PySCF and the `janpa_interface` module:

    import numpy as np
    from pyscf import gto, scf
    from janpa.janpa_interface import generate_CLPO_molecule_edges

    # 1. Build molecule and run SCF
    mol = gto.M(atom="H 0 0 0; H 0 0 1", basis="sto-3g", unit="Angstrom")
    mf = scf.RHF(mol).run()

    # 2. Generate CLPO orbitals and extract the bonding graph
    # silent=True suppresses the verbose JANPA terminal output
    clpo_coeff, edges = generate_CLPO_molecule_edges(mol, mf, silent=True)

    print("CLPO Coefficient matrix shape:", clpo_coeff.shape)
    print("Extracted CLPO Edges:", edges)

The `edges` list contains tuples representing the localized orbital pairs:
- `(i,)` represents a Lone Pair (LP) on a single center.
- `(i, i+1)` represents a Bonding (BD) and Antibonding (NB) pair between two centers.

## References and Attribution

This Python port is based on the original JANPA software developed by Tymofii Nikolaienko, L. A. Bulavin, and D. M. Hovorun. All theoretical foundations and original algorithms belong to the respective authors. Please refer to the [original JANPA webpage](https://janpa.sourceforge.net/) for full references, copyright information, and the original Java binaries.

If you use this software or the underlying methods in your research, please cite the original articles:

1. T. Yu. Nikolaienko, L. A. Bulavin; *Localized orbitals for optimal decomposition of molecular properties*, Int. J. Quantum Chem. (2019), Vol.119, page e25798. [DOI: 10.1002/qua.25798](https://doi.org/10.1002/qua.25798)
2. T. Yu. Nikolaienko, L. A. Bulavin, D. M. Hovorun; *JANPA: an open source cross-platform implementation of the Natural Population Analysis on the Java platform*, Comput. Theor. Chem. (2014), V.1050, P.15-22. [DOI: 10.1016/j.comptc.2014.10.002](https://doi.org/10.1016/j.comptc.2014.10.002)