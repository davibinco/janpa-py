package onpa;

import ProgramOptions.*;

/**
 * A class for storing various program output options
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
 *
 * @author (c) Tymofii Nikolaienko, 2014
 */
public class ono_options {

    public boolean glbPrint = false; // 'global' print
    public boolean dont_print_matrices = true;

    //--------------------------------------------------------------------------
    // INPUT / OUTPUT OPTIONS
    //--------------------------------------------------------------------------
    // filenames for matrix and/or other data export
    public OptionParameter Input_Molden_File =
            new OptionParameter(null, "-i", "", 
                    "input molden file");

    public OptionParameter _FirstOption = Input_Molden_File; // begins the option parameter list
    //--------------------------------------------------------------------------

    // some of the 'most useful' options:
    public OptionParameter Charges_File = 
            new OptionParameter(_FirstOption, "-npacharges", "", 
                    "file for exporting plain list of NPA charges");
    
    public OptionParameter NAO_Molden_File = 
            new OptionParameter(_FirstOption, "-NAO_Molden_File", "", 
                    "Export NAOs in MOLDEN format (AO basis)");
    public OptionParameter LHO_Molden_File = 
            new OptionParameter(_FirstOption, "-LHO_Molden_File", "", 
                    "Export LHOs in MOLDEN format (AO basis)");
    public OptionParameter CLPO_Molden_File = 
            new OptionParameter(_FirstOption, "-CLPO_Molden_File", "", 
                    "Export CLPOs in MOLDEN format (AO basis)");
    
    public OptionParameter PrintGeometry = 
            new OptionParameter(_FirstOption, "-PrintGeom", false, 
                    "print geometry of the system.");

    public OptionParameter WiebergBondOrders_File = 
            new OptionParameter(_FirstOption, "-WibergBondOrders_File", "", 
                    "file for exporting a matrix of Wiberg bond orders \n calculated in NAO basis");

    //--------------------------------------------------------------------------
    // 'Advanced' options:
    
    // Density / overlap matrices:
    public OptionParameter D_Matrix_File =
            new OptionParameter(_FirstOption, "-D_Matrix_File", "", 
                    "a file name for D (density) matrix D export");
    public OptionParameter S_Matrix_File =
            new OptionParameter(_FirstOption, "-S_Matrix_File", "", 
                    "a file name for S (overlap) matrix export");
    public OptionParameter SDS_Matrix_File =
            new OptionParameter(_FirstOption, "-SDS_Matrix_File", "", 
                    "a file name for SDS matrix export");
    // After step 1:
    public OptionParameter PNAO_OverlapMatrix_File =
            new OptionParameter(_FirstOption, "-PNAO_OverlapMatrix_File", "", 
                    "a file name for PNAO overlap matrix \n export");
    public OptionParameter PNAO_SDS_Matrix_File =
            new OptionParameter(_FirstOption, "-PNAO_SDS_Matrix_File", "", 
                    "a file name for exporting SDS matrix in \n PNAO basis");

    // after steps 2&3 (division to NMB/NRB & weighted symmetrical orthogonalization of NMB):
    public OptionParameter NMB_old_Overlap_Matrix_File =
            new OptionParameter(_FirstOption, "-NMB_old_Overlap_Matrix_File", "", 
                    "");
    public OptionParameter NMB_old_SDS_Matrix_File =
            new OptionParameter(_FirstOption, "-NMB_old_SDS_Matrix_File", "", 
                    "");

    // after step 4 (Schmidt orthogonalization of each NRB to all NMB)
    public OptionParameter NRB_old_Overlap_Matrix_File =
            new OptionParameter(_FirstOption, "-NRB_old_Overlap_Matrix_File", "", 
                    "");
    public OptionParameter NRB_new_Overlap_Matrix_File =
            new OptionParameter(_FirstOption, "-NRB_new_Overlap_Matrix_File", "", 
                    "");

    // after step 5 (â€?intraatomic naturalizationâ€™ transformation within the new NRB set):
    public OptionParameter S_Matrix_after_ON2_File =
            new OptionParameter(_FirstOption, "-S_Matrix_after_ON2_File", "", 
                    "orbital overlap matrix after 2-nd intracenter \n naturalizatoin transformation");
    public OptionParameter SDS_Matrix_after_ON2_File =
            new OptionParameter(_FirstOption, "-SDS_Matrix_after_ON2_File", "", 
                    "");

    // during & after step 6 (weighted symmetrical orthogonalization of functions of the new NRB set):
    public OptionParameter NRB_Overlap_after_OW_heavy_File =
            new OptionParameter(_FirstOption, "-NRB_Overlap_after_OW_heavy_File", "", 
                    "");
    public OptionParameter NRB_Overlap_after_OS2_File =
            new OptionParameter(_FirstOption, "-NRB_Overlap_after_Schmidt2_File", "", 
                    "");
    public OptionParameter NRB_Overlap_after_OW2_final_File =
            new OptionParameter(_FirstOption, "-NRB_Overlap_after_OW2_final_File", "", 
                    "");
    // final results of step 6:
    public OptionParameter OW2_File = 
            new OptionParameter(_FirstOption, "-OW2_File", "", 
                    "2-nd WSW transformation matrix (within _full_ \n NRB subspace)");
    public OptionParameter S_Matrix_after_OW2_File =
            new OptionParameter(_FirstOption, "-S_Matrix_after_OW2_File", "", 
                    "");
    public OptionParameter SDS_Matrix_after_OW2_File =
            new OptionParameter(_FirstOption, "-SDS_Matrix_after_OW2_File", "", 
                    "");

    // after the final step 7:
    public OptionParameter SDS_NAO_File =
            new OptionParameter(_FirstOption, "-SDS_NAO_File", "", 
                    "a file for exporting SDS matrix in NAO basis");
    // S matrix should be unitary matrix => no need to export

    // Some other transformation matrices:
    public OptionParameter NAO2AO_File = 
            new OptionParameter(_FirstOption, "-NAO2AO_File", "", 
                    "a AO->NAO transformation matrix");
    public OptionParameter PNAO2AO_File = 
            new OptionParameter(_FirstOption, "-PNAO2AO_File", "", 
                    "a AO->PNAO transformation matrix");                // 22.07.2016

    public OptionParameter CLPO2LHO_File = 
            new OptionParameter(_FirstOption, "-CLPO2LHO_File", "", 
                    "a LHO->CLPO transformation matrix");
    public OptionParameter LHO2NAO_File = 
            new OptionParameter(_FirstOption, "-LHO2NAO_File", "", 
                    "a NAO->LHO transformation matrix");
    
    public OptionParameter LPO2AHO_File = 
        new OptionParameter(_FirstOption, "-LPO2AHO_File", "", 
                "a AHO->LPO transformation matrix");
    public OptionParameter AHO2NAO_File = 
            new OptionParameter(_FirstOption, "-AHO2NAO_File", "", 
                    "a NAO->AHO transformation matrix");
    
    
    //TODO: public OptionParameter PNAO2NAO_File = new OptionParameter(_FirstOption, "-PNAO2NAO_File", "", "a PNAO->NAO transformation matrix"); // TODO

    // Options related to Fock matrix calculation/export/analysis       // 06.09.2014
    public OptionParameter do_Fock = 
            new OptionParameter(_FirstOption, "-doFock", false, 
                    "whether to perform Fock matrix processing");        
    public OptionParameter Fock_AO_File  = 
            new OptionParameter(_FirstOption, "-Fock_AO_File", "", "Fock matrix in AO basis");
    public OptionParameter Fock_NAO_File = 
            new OptionParameter(_FirstOption, "-Fock_NAO_File", "", "Fock matrix in NAO basis");
    public OptionParameter PrintNBMSubmatrices = 
            new OptionParameter(_FirstOption, "-printnmbfock", false, 
                    "prints Fock and overlap matrices in NMB PNAO/NAO \n bases (works only when Fock AO matrix is available)");
    

    // Molden file export
    public OptionParameter PNAO_Molden_File = 
            new OptionParameter(_FirstOption, "-PNAO_Molden_File", "", 
                    "Export PNAOs in MOLDEN format (AO basis)");
    
    public OptionParameter AHO_Molden_File = 
            new OptionParameter(_FirstOption, "-AHO_Molden_File", "",
                    "Export AHOs in MOLDEN format (AO basis)");
    
    public OptionParameter LPO_Molden_File = 
            new OptionParameter(_FirstOption, "-LPO_Molden_File", "",
                    "Export LPOs in MOLDEN format (AO basis)");

    // matrix output format settings
    public OptionParameter MatrixLineWidth = 
            new OptionParameter(_FirstOption, "-MatrixMaxValuesPerLine", new Integer(0), 
                    "a maximum number of values per line for matrix \n export (0 = no limit)"); // TODO
    public OptionParameter MatrixFloatNumberFormat =
            new OptionParameter(_FirstOption, "-MatrixFloatNumberFormat", "%.5f", 
                    "a printf-compatible format for matrix elements \n printing"); // overrides printout.MatrixFloatNumberFormat
    
    public OptionParameter maxClpoBondIonicityThreshold =
            new OptionParameter(_FirstOption, "-maximumBondIonicity", new Double(0.90), 
                    "maximum allowable ionicity of CLPO bonding (BD) orbitals");
    public OptionParameter RyOccPrintThreshold =
            new OptionParameter(_FirstOption, "-RyOccPrintThreshold", new Double(1e-3), 
                    "occupancy threshold for RY orbital printing\n set to negative value to print all RY orbitals");
    public OptionParameter HybrOptOccConvThresh =
            new OptionParameter(_FirstOption, "-HybrOptOccConvThresh", new Double(1e-5), 
                    "numerical convergence threshold for hybrid optimization");
    public OptionParameter HybrOptMaxIter =
            new OptionParameter(_FirstOption, "-HybrOptMaxIter", new Integer(1000), 
                    "maximum number of iterations during hybrids optimization");
    
    
    // general output options
    public OptionParameter VerbosePrint = 
            new OptionParameter(_FirstOption, "-verboseprint", false, 
                    "print some additional information");

    //--------------------------------------------------------------------------
    // CALCULATION OPTIONS
    //--------------------------------------------------------------------------
    // options for basis set and MO orthonormality checking
    public double bf_nrm2_dev_threshold = 1.0E-07;    // maximum allowable value of  | <BasisFunction, BasisFunction> - 1 |
    public double mo_nrm2_dev_threshold = 1.0E-04;    // maximum allowable value of  | <MO_i, MO_i> - 1
    public double mo_off_diag_threshold = 1.0E-04;    // maximum allowable value of  | <MO_i, MO_j> | with i =/= j
    public double S_AO_lin_dep_thresh = 1.0E-9;            // minimum allowable value for eigenvalue of AO overlap matrix
    public boolean fail_if_lin_dep = false;     // whether to fail if (nearly) linearly dependent basis functions have been found
    // A threshold for an attempt to renormalize MO"
    public double mo_non1_renorm_thresh = mo_nrm2_dev_threshold;  //public double mo_non1_renorm_thresh = 1.0E-05;
    //------------------------------------
    public boolean do_BS_renormalize = true;
    public boolean do_MO_renormalize = false;
    //--------------------------------------------------------------------------
    // some VERY specific options
    // nao creation:
    public OptionParameter NRB_OW_Straightforward = 
            new OptionParameter(_FirstOption, "-directNRBwsw", false,
                    "whether not to use Schmidt-Lowdin scheme for (WSW)^(1/2) \n calculation at spet 6"); // whether to calculate (WSW)^(-1/2) for NRB transformation directly, rather than by Weight-Schmidt-Lowdin method
    public OptionParameter Heavy_NRB_Threshold = 
            new OptionParameter(_FirstOption, "-heavyNRBthreshold", new Double(1.0E-4),
                    "a threshold value for 'heavily' occupied NRBs (valid only \n without -directNRBwsw)");
    // matrix operations:
    /*public OptionParameter No_Enhanced_Matrix_Transform =
            new OptionParameter(_FirstOption, "-DontUseLDL", false, 
                    "turns off using of LDLT transformation in matrix \n transformations"); // whether to use LDL decomposition for transforming symmetric matrices
    
    public OptionParameter MatrixTransform_Force_Symmetric =
            new OptionParameter(_FirstOption, "-forcesymmetric", false, 
                    "whether to apply Matrix_Symmetrize_Check() after U.A.U^T \n"+
                    " calculated by a 'direct' method (valid only with -DontUseLDL; _strongly_ \n"+
                    " recommended value: true)"); // this parameter
            // is valid only if Enhanced_Matrix_Transform == false;
            // It determines whether to apply Matrix_Symmetrize_Check() after U.A.U^T calculated by a 'direct' method
    */
    public OptionParameter Overlap_Naive =
            new OptionParameter(_FirstOption, "-overlap_straightforward", false, 
                    "whether to apply slow-but-reliable method to compute \n"+
                    " an overlap matrix (not recommended)"); // this parameter

    //--------------------------------------------------------------------------
    
    public OptionParameter Parameter_File =
            new OptionParameter(_FirstOption, "-p", "", 
                    "a file with parameters (having lower priority than the \n command line)");
    
    
}
