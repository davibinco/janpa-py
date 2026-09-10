# janpa-py: Pure Python Natural Population Analysis

A fast, dependency-free Python port of [JANPA](https://janpa.sourceforge.net/) (Natural Population Analysis and electronic structure tools). This library provides native Python implementations of the 7-step NPA algorithm, Atomic Hybrid Orbitals (AHO/LHO), and Chemist's Localized Property-optimized Orbitals (CLPO/LPO) using `numpy` and `scipy`, completely eliminating the need for Java binaries, subprocess wrappers, or temporary disk I/O.

This repository contains the Python conversion of our modified Java code for JANPA. See the `janpa-java` branch for the Java code, and refer to the [original JANPA project](https://janpa.sourceforge.net/) for the foundational work.

## Key Features

- **Pure Python:** Built entirely on `numpy`, `scipy`, and `networkx`.
- **In-Memory Execution:** No messy `.molden` or `graph` text files left on your hard drive.
- **Seamless Integration:** Designed as a drop-in replacement for PySCF and Tequila/Sunrise workflows.
- **Custom Graphs:** Supports both automatic bonding graph extraction and custom edge definitions for quantum ansatz generation (e.g., HCB-SPA).

## Installation

Clone the repository and install it into your Python environment:

    git clone https://github.com/yourusername/janpa.git
    cd janpa
    pip install .

## Quick Usage

Here is a minimal example demonstrating how to run NPA and CLPO analysis directly in memory using PySCF:

    import numpy as np
    from pyscf import gto, scf
    from janpa.io.molden import MoldenFile
    from janpa.npa.npa import run_npa
    from janpa.clpo.lpo import create_clpos, CLPOOptions
    
    # 1. Build molecule and run SCF
    mol = gto.M(atom="H 0 0 0; H 0 0 1", basis="sto-3g", unit="Angstrom")
    mf = scf.RHF(mol).run()
    
    # 2. Pass the molecule to the JANPA pipeline
    # (See the Tequila/Sunrise integration scripts in the examples folder 
    # for full active-space and custom-edge workflows)

## Tequila / Sunrise Integration

This library acts as a direct replacement for the Java-based `sunrise.CLPO` module. It exposes `generate_CLPO_molecule_edges` and `generate_HAO_molecule` to flawlessly map localized orbitals into frozen-core active spaces, allowing you to feed chemical bonding graphs directly into quantum computing ansatze.

## References and Attribution

This Python port is based on the original JANPA software developed by Tymofii Nikolaienko, L. A. Bulavin, and D. M. Hovorun. All theoretical foundations and original algorithms belong to the respective authors. Please refer to the [original JANPA webpage](https://janpa.sourceforge.net/) for full references, copyright information, and the original Java binaries.

If you use this software or the underlying methods in your research, please cite the original articles:

1. T. Yu. Nikolaienko, L. A. Bulavin; *Localized orbitals for optimal decomposition of molecular properties*, Int. J. Quantum Chem. (2019), Vol.119, page e25798. [DOI: 10.1002/qua.25798](https://doi.org/10.1002/qua.25798)
2. T. Yu. Nikolaienko, L. A. Bulavin, D. M. Hovorun; *JANPA: an open source cross-platform implementation of the Natural Population Analysis on the Java platform*, Comput. Theor. Chem. (2014), V.1050, P.15-22. [DOI: 10.1016/j.comptc.2014.10.002](https://doi.org/10.1016/j.comptc.2014.10.002)

## License

This Python port is released under the MIT License. The original JANPA Java software is subject to its own licensing and copyright as detailed on the original project's website.