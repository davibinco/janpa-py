package molden2molden;

import ProgramOptions.*;

/**
 * A command-line parsing construction for molden2molden program
 * (c) Tymofii Nikolaienko, 2014
 * Revision: 30.Mar.2014
 *
 * This file is a part of the JANPA project. 
 *
 * Copyright (c) 2014, Tymofii Nikolaienko
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 
 * 3. All advertising and/or published materials mentioning features or use 
 *    of this software must display the following acknowledgement:
 * 
 *       This product includes components from JANPA package of programs
 *       ( http://janpa.sourceforge.net/ ) developed by Tymofii Nikolaienko
 * 
 * 4. Neither the name of the developer, Tymofii Nikolaienko,  nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * 
 * 5. In case if the code of JANPA package code and/or its parts and/or any data 
 *    produced with JANPA package of programs are published, the following citations
 *    for the JANPA package of programs should be given:
 *
 *     1) T.Y.Nikolaienko, L.A.Bulavin; Int. J. Quantum Chem. (2019), 
 *        Vol.119, page e25798, DOI: 10.1002/qua.25798
 *     2) T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem. (2014),
 *        V. 1050, P. 15-22, DOI: 10.1016/j.comptc.2014.10.002
 * 
 * THIS SOFTWARE IS PROVIDED BY ''AS IS'' AND ANY EXPRESS OR IMPLIED 
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 * IN NO EVENT SHALL Tymofii Nikolaienko BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
public class m2mOptions {
        // possible operations
        public OptionParameter op_coordsAngstroms = new OptionParameter(null, "-CoordsAngstrom", /*new Boolean(*/false/*)*/,
                "forces coordinates to be written in Angstroms\n\t\t(default is to convert to Bohrs);");
        public OptionParameter _first = op_coordsAngstroms;
        public OptionParameter op_normalizeBF = new OptionParameter(_first, "-NormalizeBF", false,
                "forces all basis functions to be unity-normalized\n\t\tby correcting contraction coefficients;");
        public OptionParameter op_cart2spher = new OptionParameter(_first, "-cart2pure", false,
                "convert cartesian harmonics to pure spherical\n\t\tones;");
        public OptionParameter IgnoreLowerLTerms = new OptionParameter(_first, "-IgnoreLowerL", false,
                "ignore terms with L<L_Main during cartesian->pure\n\t\tconversion;" +
                "\n\nORCA-specific operations:\n");

        // ORCA-specific operations
        public OptionParameter op_orca3signs = new OptionParameter(_first, "-orca3signs", false,
                "make MO coefficients compatible with the spherical harmonics\n\t\tdefinitions used by MOLDEN;");
        public OptionParameter op_from_orca3_BF = new OptionParameter(_first, "-fromorca3bf", false,
                "make basis set contraction coefficients produced by orca_2mkl\n\t\tcompatible with MOLDEN definition;");
        public OptionParameter op_to_orca3_BF = new OptionParameter(_first, "-toorca3bf", false,
                "convert from MOLDEN-style contraction coefficients to the\n\t\tORCA-style (inverse to -FromOrca3BF" +
                "\n\nPSI4-specific operations:\n");
        
        // PSI4-specific operations
        public OptionParameter op_from_psi4v1_MO = new OptionParameter(_first, "-frompsi4v1mo", false,
                "make MO coefficients produced by PSI4 1.0 release\n\t\tcompatible with MOLDEN definition;");                
        public OptionParameter op_from_psi4b4_BF = new OptionParameter(_first, "-frompsi4b4bf", false,
                "make basis set contraction _and_ MO coefficients produced by PSI4-beta4/beta5\n\t\tcompatible with MOLDEN definition;");        
        public OptionParameter op_to_psi4b4_BF = new OptionParameter(_first, "-topsi4b4bf", false,
                "convert from MOLDEN-style contraction coefficients to the\n\t\tPSI4-style (inverse to -FromPsi4b4BF);"+
                "\n\nInput and output file names:\n");

        // I/O file names
        public OptionParameter fname_in = new OptionParameter(_first, "-i", "",
                "INPUT molden file name;");        // INPUT molden file
        public OptionParameter fname_in_BINAR = new OptionParameter(_first, "-bi", "",
                "BINARY INPUT molden file name;");        // BINARY INPUT molden file
        public OptionParameter fname_out = new OptionParameter(_first, "-o", "",
                "OUTPUT molden file name;");       // OUTPUT molden file
        public OptionParameter fname_out_BINAR = new OptionParameter(_first, "-bo", "",
                "BINARY OUTPUT molden file name;");       // OUTPUT molden file
        
        public OptionParameter File_GEOM47 = new OptionParameter(_first, "-geom47", "",
                "a nbo3-compatible .47 file with the geometry\n\t"+
                "to be used for replacing the geometry in the INPUT molden.");          // a nbocompatible .47 file with MO data

        public OptionParameter File_BS47 = new OptionParameter(_first, "-bs47", "",
                "a nbo3-compatible .47 file with the basis set contraction coefficients\n\t"+
                "to be used for replacing the basis set in the INPUT molden.");          // a nbocompatible .47 file with MO data
        public OptionParameter File_DS47 = new OptionParameter(_first, "-ds47", "",
                "a nbo3-compatible .47 file with AO overlap and MO density matrices\n\t"+
                "to be used for replacing orbitals from INPUT molden.");          // a nbocompatible .47 file with MO data

        
}
