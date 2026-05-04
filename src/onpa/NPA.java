package onpa;

/**
 * A file containing NPA class which implemets AO->NAO transformation
 * and the Natural Population Analysis 
 *
 * Version: 15.Feb.2014 / 26.Oct.2013
 * Copyright (c) Tymofii Nikolaienko
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
 */

//import MatrixHelper.LDLtransform;
import JGints.*;
import Jama.*;
import Polynom3D.*;
import moldenio.*;
import java.io.*;
import java.util.Arrays;

import ProgramOptions.*;
import MatrixHelper.*;

//==============================================================================
/**
 * A class performing NAO/NPA analysis
 * @author timn
 * @version 08.Feb.2014 / 23.Oct.2013
 */
public class NPA {
    public BasisFunction[] Basis = null;
    public AtomicCenter[] Centers = null;
    private boolean _data_loaded = false;
    public JGintsCyl bsIntegrals;               // may appear useful since it typically contains more integrals than just an overlap matrix
    
    public Matrix D_Global = null;
    public Matrix OverlapMatrix = null;        // overlap matrix for atomic basis functions 
    public Matrix S05DS05 = null;              // S^(1/2).D.S^(1/2), used for Lowdin population analysis
    public Matrix SDS = null;
    private Matrix OverlapMatrix_SQRT = null;   // S^(1/2) (needed for Lowdin population analysis)
    //-----------------------------------
    // The following three matrices are being created by Create_NAOs():
    public Matrix Overlap_NAO = null;          // Overlap matrix in NAO basis (should be a unitary matrix)
    public Matrix SDS_NAO = null;              // P = S.D.S matrix in NAO basis
    public Matrix NAO_2_AO = null;             // the NAO<-AO basis functions transformation matrix
    public Matrix PNAO_Overlap_Matrix = null;
    public BasisFunction[] NAO;
    public double[] NPA_charges = null;        // NPA charges
    BasisFunction[] PNAOs;
    public String[] AO_Names;
    public String[] PNAO_Labels;
    //-----------------------------------        
    Polynom3D[][] Quick_YLM ;// = SphericalHarmonics.Get_Quick_YLM();
    
    private ono_options options; // misc. work options
    public static PrintStream out = System.out;
    public boolean verbose_print = false; // whether to print some additional messages
    
    
    private double[][] AORadialPartsOverlap; // ovrlap matrix for radial parts only (not used!)

    private EigenEngine ee = null;
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // a constructor
    public NPA(ono_options options) {
        if (options == null) options = new ono_options(); // take all settings by default
        this.options = options; // save the options
        try {
            // pass some options to our 'colleagues'
            printout.MatrixFloatNumberFormat = options.MatrixFloatNumberFormat.get_String();
            printout.MatrixLineWidth = options.MatrixLineWidth.get_int();
            verbose_print = options.VerbosePrint.get_boolean();
        }  catch (Exception e) {}

        // Create an object for matrix/eigenvalue problem manipulations
        ee = new EigenEngine();
        //ee.MatrixTransform_Force_Symmetric = options.MatrixTransform_Force_Symmetric; 
        //ee.No_Enhanced_Matrix_Transform = options.No_Enhanced_Matrix_Transform;
        ee.glbPrint = options.glbPrint;
    }
    //--------------------------------------------------------------------------    
    /** Loads MO data from MOLDEN file, sets Basis and Centers arrays,
     * builds basis function overlap matrix (OverlapMatrix field),
     * checks the normalization of basis function and for their linear (in)dependence,
     * checks for normalization of molecular orbitals, and their orthogonality;
     * Creates the density matrix (Global_D) and the S.D.S matrix
     * 
     * Note that src_molden is passed "by reference" and its AO basis set
     * and MO coefficients are being modified!!!
     *
     */
    boolean Load_MO_From_MOLDEN(MOLDEN_IO src_molden) throws Exception {
        
        if (src_molden == null) return false;
        _data_loaded = false;

        Basis = src_molden.Basis; // TODO: make a copy!
        Centers = src_molden.Centers; // TODO: make a copy!

        final boolean do_magic = false;
        if (do_magic) {
            // Adapt coefficient to normalized spherical harmonics
            // chi[i] = SUM( c[i,j]·r[j]·YLM[j] ) => SUM( c[i,j]·sqrt(YLM.YLM)·r[j]·1/sqrt(YLM.YLM)·YLM[j] ) where r[j] is j-th radial part
            // For that we will need a norm of Y_LM
            double [][] YLM_Norm2 = SphericalHarmonics.Get_Quick_YLM_Norm2(Quick_YLM);
            // corrent radial parts
            // first - the copies
            for (int bs=0; bs<Basis.length; bs++)
                for (int rpart=0; rpart<Basis[bs].coefs.length; rpart++)
                    Basis[bs].coefs[rpart] *= Math.sqrt(YLM_Norm2[Basis[bs].L][Basis[bs].L+Basis[bs].m]);
            // and for their copies:
            // I know that this twice the same! :)
            for (int rpart=0; rpart<src_molden.RadialParts.length; rpart++) {
                double tmp = Math.sqrt(YLM_Norm2[ src_molden.RadialParts[rpart].LUsedWith ][ 0 ]);
                for (int c=0; c<src_molden.RadialParts[rpart].Coefs.length; c++)
                    src_molden.RadialParts[rpart].Coefs[c] *= tmp;
            }

            // now correct YLMs
            for (int L=0; L<Quick_YLM.length; L++)
                for (int m=-L; m<=L; m++)
                    for (int mu=0; mu<Quick_YLM[L][L+m].coefs.length; mu++)
                        Quick_YLM[L][L+m].coefs[mu]  /= Math.sqrt(YLM_Norm2[L][L+m]); // divide each coefficient of the polynom
            //done!
            YLM_Norm2 = null;
        }
        //double[][] X = SphericalHarmonics.Get_Quick_YLM_Norm2(Quick_YLM);
        //X[1][1] += 1;
        int NBas = Basis.length;

/*
        // now create a density matrix
        out.println("Building density matrix D...");
        D_Global = new Matrix(NBas, NBas, 0.0);
        double tmp;
        // set terms
        for (int basf1=0; basf1<NBas; basf1++)
           for (int basf2=basf1; basf2<NBas; basf2++) {
                // loop pver orbitals
                tmp = 0;
                for (int mo=0; mo<molden.MOs.length; mo++)
                    tmp += molden.MOs[mo].BS_Coefs[basf1] * molden.MOs[mo].BS_Coefs[basf2] * molden.MOs[mo].Occupancy;
                D_Global.set(basf1, basf2, tmp);
                D_Global.set(basf2, basf1, tmp);
            }
        //
        if (print_D_Matrix) 
            _Print_Matrix_With_Header(D_Global, "The D density Matrix in spherical function AO basis is:");
 *
 */
        //-------------------------------------------------
        final boolean _do_radial_overlaps = false;
        //-------------------------------------------------
        if (_do_radial_overlaps) {
            out.println("Building same-center-radial overlap matrix...");
            int N = src_molden.RadialParts.length;
            AORadialPartsOverlap = new double[N][N];
            OverlapIntegrals OIR = new OverlapIntegrals();
            // loop over basis functions (!) - we can not loop over radial parts since we need to have Ls for those radial parts!
            int RID1,RID2;
            for (int bs1=0; bs1<Basis.length; bs1++)
                for (int bs2=bs1; bs2<Basis.length; bs2++) {
                    // calculate integrals for radial parts on the same center only!
                    if (Basis[bs1].Center_ID == Basis[bs2].Center_ID) {
                        // get IDs of their radial parts
                        RID1 = Basis[bs1].RadialPart_ID;
                        RID2 = Basis[bs2].RadialPart_ID;
                        //
                        AORadialPartsOverlap[RID1][RID2] = OIR.SameCenterRadialOverlapIntegral(
                            src_molden.RadialParts[RID1].Coefs,
                            src_molden.RadialParts[RID1].Exponents,
                            Basis[bs1].L,
                            src_molden.RadialParts[RID2].Coefs,
                            src_molden.RadialParts[RID2].Exponents,
                            Basis[bs2].L);
                        //
                        AORadialPartsOverlap[RID2][RID1] = AORadialPartsOverlap[RID1][RID2];
                    }
                }
            //Matrix AORadialPartsOverlapMatrix = new Matrix(AORadialPartsOverlap);
            //AORadialPartsOverlapMatrix.print(13, 7);
        }
        //-------------------------------------------------
        //
        MO[] MOs = src_molden.MOs;
        out.println("Building overlap S and dipole matrices...");

        if (!options.Overlap_Naive.get_boolean()) {
            // Use quick overlap builder from JGints
            
            if (true) {
                // use a newer version of JGints -- JGintsCyl
                bsIntegrals = new JGintsCyl();
                bsIntegrals.ImportBasisFromMolden(src_molden);

                ProgTimer t = new ProgTimer();
                t.Start();

                bsIntegrals.Build_Ovarlap_Matrix_PURE();

                out.printf("Overlap & dipole integrals evaluated in %.3f seconds%n",  t.Stop()/1000);

                OverlapMatrix = new Matrix(bsIntegrals.OverlapMatrix);
                
            } else {
                // using an older version of JGints
                JGints jg = new JGints();
                jg.ImportBasisFromMolden(src_molden);
                jg.ylm = Quick_YLM;
                jg.ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(Quick_YLM);

                ProgTimer t = new ProgTimer();
                t.Start();

                jg.Build_Ovarlap_Matrix_PURE();

                out.printf("Overlap integrals evaluated in %.3f seconds%n",  t.Stop()/1000);

                OverlapMatrix = new Matrix(jg.OverlapMatrix);
            }
            //jg._test();
        } else {
            // use .OverlapWith method for each basis function (~3...5 times slower)
            OverlapMatrix = new Matrix(NBas, NBas, 0.0);
            // set Quick_YLM table to all basis functions of MOLDEN
            for (int i=0; i<NBas; i++) Basis[i].Quick_YLM = Quick_YLM;
            // calculate overlaps
            for (int basf1=0; basf1<NBas; basf1++)
               for (int basf2=basf1; basf2<NBas; basf2++) {
                    // loop pver orbitals
                    double tmp = Basis[basf1].OverlapWith( Basis[basf2] );
                    OverlapMatrix.set(basf1, basf2, tmp);
                    OverlapMatrix.set(basf2, basf1, tmp);
                }
        }


        // check normalizations for basis functions and MOs
        boolean bs_norms_ok = false;
        boolean mo_norms_ok = false;
        boolean mo_orth_ok  = false;

        // It might be usefull to re-normalize all basis functions
        out.println("Checking whether the basis functions are unity-normalized...");
        double[] bs_norms = new double[NBas];
        double max_dev = 0;
        int max_dev_index = -1;
        for (int basf=0; basf<NBas; basf++) {
                bs_norms[basf] = OverlapMatrix.get(basf, basf);
                // store maximum non-unity-normalized basis function number
                if (Math.abs( bs_norms[basf] - 1 ) > max_dev) {
                       max_dev = Math.abs( bs_norms[basf] - 1 );
                       max_dev_index = basf;
                }
                //bs_norms[basf] = Math.sqrt( bs_norms[basf] );
        }
        out.printf(" Maximum deviation of the basis function norm2 from unity: BFN %d (0-based num.), max|norm2-1| = %.3E %n", max_dev_index, max_dev );
        if (( max_dev > options.bf_nrm2_dev_threshold ) && (options.do_BS_renormalize)) {
            out.println("Trying to make basis functions unity-normalized...");
            for (int basf1=0; basf1<NBas; basf1++)
               for (int basf2=0; basf2<NBas; basf2++)
                   OverlapMatrix.set(basf1, basf2, OverlapMatrix.get(basf1,basf2) / bs_norms[basf1] / bs_norms[basf2] );
        }
        if (( max_dev < options.bf_nrm2_dev_threshold ) || (options.do_BS_renormalize))
                bs_norms_ok = true;

        //---------------------------------------------------------------------
        out.println("Checking for the eigenvalues (linear (in)dependency) of the basis function overlap matrix...");
        EigenvalueDecomposition OverlapMatrix_eig = OverlapMatrix.eig();
        double S_min_abs = 0; // it should be acceptable without 'abs' as well...
        // find the smallest eiganvalue of the overlap matrix
        double[][] s_eig_arr = OverlapMatrix_eig.getD().getArray();
        S_min_abs = Math.abs(s_eig_arr[0][0]);
        for (int i=1; i<NBas; ++i)
            if (Math.abs(s_eig_arr[i][i]) < S_min_abs) S_min_abs = Math.abs(s_eig_arr[i][i]);
        out.printf("The smallest eigenvalue of the basis function overlap matrix: %7.3E %n", S_min_abs);
        if (S_min_abs < options.S_AO_lin_dep_thresh) {
            WarningManager.warning_printf("WARNING: the basis set functions seems to be linearly dependent!");
            if (options.fail_if_lin_dep) return false;
        }
        
        // replace eigenvalues with their square roots
        for (int i=0; i<s_eig_arr.length; ++i)
            s_eig_arr[i][i] = Math.sqrt(s_eig_arr[i][i]);
        // calculate S^(1/2)
        OverlapMatrix_SQRT = OverlapMatrix_eig.getV().times(new Matrix(s_eig_arr)).times(OverlapMatrix_eig.getV().transpose());
        
        //---------------------------------------------------------------------

        out.println("Checking whether the orbitals are unity-normalized...");
        // we'll use the same bs_norms, max_dev and max_dev_bf vars. as above
        max_dev = 0;
        max_dev_index = -1;
        int N_MO = MOs.length;
        Matrix MO_coef_Matrix = new Matrix(N_MO, NBas);
        for (int mo=0; mo<N_MO; mo++)
            for (int basf1=0; basf1<NBas; basf1++)
                MO_coef_Matrix.set(mo, basf1, MOs[mo].BS_Coefs[basf1]);

        // Below we use MO_coef_Matrix to compute full MO overlap matrix;
        // However, if MO renormalization is allowed, it is neccessary to
        // perform this renormalization BEFORE building full MO overlap matrix
        if (options.do_MO_renormalize) {
            double[] MO_norm2s = new double[N_MO];
            for (int mo=0; mo<N_MO; mo++) {
                // D_Global can NOT be used for calculating MO norms...
                MO_norm2s[mo] = 0;
                for (int basf1=0; basf1<NBas; basf1++) {
                    MO_coef_Matrix.set(mo, basf1, MOs[mo].BS_Coefs[basf1]);

                    for (int basf2=0/*(basf1+1)*/; basf2<NBas; basf2++)
                       MO_norm2s[mo] += /*2 */ MOs[mo].BS_Coefs[basf1] * MOs[mo].BS_Coefs[basf2] * OverlapMatrix.get(basf1, basf2);
                   }
                // store maximum non-unity-normalized basis function number
                if (Math.abs( MO_norm2s[mo] - 1 ) > max_dev) {
                       max_dev = Math.abs( bs_norms[mo] - 1 );
                       max_dev_index = mo;
                }
            }
            out.printf(" Maximum deviation of the orbital norm2 from unity: MO %d (0-based num.), max|norm2-1| = %.3E %n", max_dev_index, max_dev );
            if (( max_dev > options.mo_non1_renorm_thresh ) /*&& (do_MO_renormalize)*/) {
                out.println("Trying to make MOs unity-normalized...");
                for (int mo=0; mo<NBas; mo++)
                    for (int basf1=0; basf1<NBas; basf1++)
                        MOs[mo].BS_Coefs[basf1] /= bs_norms[mo];
            }
        }

        out.println("Checking MO overlap matrix...");
        double[][] MO_Overlap = ee.TransformMatrixToNewBasis(OverlapMatrix, MO_coef_Matrix, true).getArray(); //  MO_coef_Matrix.times(OverlapMatrix).times(MO_coef_Matrix.transpose()).getArray();
        // check for maximum deviation of MO norm from unity (MO might have been renormalized!)
        // and checm for off-diagonal elements
        int off_diag_i = 0, off_diag_j = off_diag_i+1, non1_mo = 0;
        for (int i=0; i<N_MO; ++i) {
            if (Math.abs(MO_Overlap[i][i] - 1) > Math.abs(MO_Overlap[non1_mo][non1_mo] - 1)) non1_mo = i;
            // check for off-diagonal elements
            for (int j=(i+1); j<N_MO; ++j)
                if ((Math.abs(MO_Overlap[i][j]) > Math.abs(MO_Overlap[off_diag_i][off_diag_j])) &&
                    (MOs[i].Spin == MOs[j].Spin)) { // PSI4 writes the same orbitals for alpha- and beta-electrons: it is ok,
                                                                  // but such an orbitals have their spetial part non-orthogonal =>
                                                                  // => check the orbitals of the same spin only! => no use of this.max_offdiag()
                    off_diag_i = i;
                    off_diag_j = j;
                }
        }
        double mo_nrm_dev = Math.abs(MO_Overlap[non1_mo][non1_mo]-1);
        out.printf(" Maximum of MO |norm2-1|: %7.3E (MO %4d)%n", mo_nrm_dev, non1_mo+1);
        if (mo_nrm_dev < options.mo_nrm2_dev_threshold)
            mo_norms_ok = true;
        else
            WarningManager.warning_printf("WARNING: molecular orbital normalization problem! This may spoil MO occupancies...");

        if (N_MO > 1) { // if there is only 1 MO, off_diag_j==1 will attemtp to access non-existing array element...
            double mo_ov = Math.abs(MO_Overlap[off_diag_i][off_diag_j]);
            out.printf(" Maximum absolute value of off-diagonal MO overlap element: %7.3E (< MO %4d | MO %4d >)%n",
                    mo_ov, off_diag_i+1, off_diag_j+1);
            if (mo_ov < options.mo_off_diag_threshold) // if a molden file with 5 digits in MO coefs has been used
                mo_orth_ok  = true;
            else
                WarningManager.warning_printf("WARNING: molecular orbitals seem to be non-orthogonal ?!");
        } else
            mo_orth_ok  = true;

        // Print verification summary
        if ((!bs_norms_ok) || (!mo_norms_ok) || (!mo_orth_ok))
            WarningManager.warning_printf("WARNING: input data seems to be improper!%n");
        else
            out.println("First-order reduced density matrix is OK.");
        out.println();

        // check whether we have a closed-shell system
        int nAlpha = 0;
        for (int i=0; i<src_molden.MOs.length; ++i)
            if (src_molden.MOs[i].Spin > 0) ++nAlpha;
        if ((nAlpha !=0) && (nAlpha != src_molden.MOs.length) &&
                (nAlpha != src_molden.MOs.length/2)) {
            WarningManager.warning_printf("WARNING: the system seems to be opened-shell!%n"+
            "Note that this case is current not supported so that "
                    +"the obtained results will be meaningless!");
        }

        //----------------------------------------------------------------------

        AO_Names = Create_NAO_Labels(Basis);
        // overlap matrix is ready now!
        printout.Print_Matrix(OverlapMatrix, "The Overlap Matrix in spherical function AO basis:", 
                AO_Names, AO_Names, options.S_Matrix_File.get_String());

        // now create a density matrix
        out.print("Building density matrix D... ");
        D_Global = new Matrix(NBas, NBas, 0.0);
        double tmp;
        // set terms
        for (int basf1=0; basf1<NBas; basf1++)
           for (int basf2=basf1; basf2<NBas; basf2++) {
                // loop pver orbitals
                tmp = 0;
                for (int mo=0; mo<MOs.length; mo++)
                    tmp += MOs[mo].BS_Coefs[basf1] * MOs[mo].BS_Coefs[basf2] * MOs[mo].Occupancy;
                D_Global.set(basf1, basf2, tmp);
                D_Global.set(basf2, basf1, tmp);
            }
        out.println("done.");
        //
        printout.Print_Matrix(D_Global, "The D density Matrix in spherical function AO basis:",
                AO_Names, AO_Names, options.D_Matrix_File.get_String());


        out.print("Building D.S... ");
        Matrix DS = D_Global.times(OverlapMatrix);
        out.println("done.");

        double Tot_El_Num = DS.trace();
        out.printf("Total number of electrons: %f\n", Tot_El_Num);
        double Sum_Nucl_Charge = 0; // the sum of nuclei charges
        for (int c=0; c<Centers.length; ++c)
            Sum_Nucl_Charge += Centers[c].Z;
        out.printf("Sum of electrons charges and the nuclei charges: %7.5f %n", Sum_Nucl_Charge - Tot_El_Num);
        out.println();


        out.println("Performing Mulliken and Lowdin population analyses...");
        Matrix A_Lowdin = OverlapMatrix_SQRT.times(D_Global).times(OverlapMatrix_SQRT); // no need in explicit symmetrization here
        S05DS05 = A_Lowdin;
        double[] MullikenPopulation = new double[Centers.length];
        double[] LowdinPopulation = new double[Centers.length];
        for (int i=0; i<NBas; ++i) {
            LowdinPopulation[ Basis[i].Center_ID-1 ] += A_Lowdin.get(i, i);            
            MullikenPopulation[ Basis[i].Center_ID-1 ] += DS.get(i, i);
        }
        String LowdinChargeFmt = "%7s\t%8.5f\t%8.5f\t%8.5f\t%8.5f%n";
        out.printf("%7s\t%s\t%s\t%s\t%s%n", "Atom", "Mulliken  ", "Lowdin    ", "Mulliken", "Lowdin");
        out.printf("%7s\t%s\t%s\t%s\t%s%n",     "", "Population", "Population", "Charge  ", "Charge");
        for (int c=0; c<Centers.length; ++c)
            out.printf(LowdinChargeFmt, String.format("%s%d",Centers[c].Name, c+1),
                    MullikenPopulation[c], LowdinPopulation[c],
                    Centers[c].Z - MullikenPopulation[c], Centers[c].Z - LowdinPopulation[c]);
        out.println();
        
 
        out.print("Building S.D.S... ");
        SDS = ee.TransformMatrixToNewBasis(D_Global, OverlapMatrix, true);   // this is OK since OverlapMatrix is symmetrical (S^T = S)
        

        out.println("done.");
        //
        printout.Print_Matrix(SDS, "S.D.S matrix in spherical function AO basis:",
                AO_Names, AO_Names, options.SDS_Matrix_File.get_String());
        

        out.println();

        //
        _data_loaded = true;
        return true;
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // The "heart" of NAO procedure - the "intracenter naturalizer".
    /** Performs intracenter orthogonalization of a given basis.
     * (corresponds to steps "2b" and "2c" Appendix A in [J.Chem.Phys.(1985), 83, 735-746]
     * <!-- INPUT: -->
     * @param InitialBasis[] SHOULD have Center_ID and RadialPart_ID field set with a proper values
     * @param OverlapMatrix should contain proper overlap matrix of InitialBasis[] functions and
     * @param SDS is used to get sub-blocks for eigenvalue problem.
     * <!-- OUTPUT: -->
     * @return
     * The function returns linear combinations of InitialBasis functions; their Occupancies are set to weights
     * obtained as eigenvalues of a sub-matrix problems.
     * 
     * <p><b><u>NOTES:</u></b><br>
     *        1. InitialBasis.coefs and InitialBasis.exponents are not used, so that they can be set to null)<br>
     *        2. Center_ID is assumed to be 1-based !
     */
    BasisFunction[] IntracenterBasisOrthogonalization(BasisFunction[] InitialBasis, int nAtoms, Matrix SDS, Matrix OverlapMatrix) {
        if (!_data_loaded) return null;
        
        // allocate memory
        BasisFunction[] result = new BasisFunction[InitialBasis.length];

        // It CAN happen, if there are no NRBs, for example - typically, if STO-nG minimal basis set is used / Bugfix 14.Oct.2014
        boolean InputBasisNotAmpty = (InitialBasis.length != 0);
        // get Lmax for each atom
        int[] LMaxOfAtom = new int[ nAtoms ]; // default values in "any" array are zeros (see 4.12.5. of http://docs.oracle.com/javase/specs/jls/se7/html/jls-4.html)
        for (int i=0; i<InitialBasis.length; i++)
            if (InitialBasis[i].L > LMaxOfAtom[ InitialBasis[i].Center_ID-1 ])
                LMaxOfAtom[ InitialBasis[i].Center_ID-1 ] = InitialBasis[i].L;


        // Get the number of basis functions with a given L per each center
        int[] BasisFunctionsPerAtom = new int[ nAtoms ];
        for (int i=0; i<InitialBasis.length; i++)
            BasisFunctionsPerAtom[ InitialBasis[i].Center_ID-1 ] ++;
        int MaxBasisFunctionsPerAtom = 0;
        for (int i=0; i<nAtoms; i++)
            if (BasisFunctionsPerAtom[i]>MaxBasisFunctionsPerAtom) MaxBasisFunctionsPerAtom=BasisFunctionsPerAtom[i];
        //
        // Make a list of basis functions making them ordered as:
        // [Atom1.BS_a, Atom1.BS_b, ...] [Atom2.BS_c, Atom2.BS_d,...] [Atom3.BS_e, ...] ...
        //
        int[] SortedBasisFunctions = new int[ InitialBasis.length ];
        // To be able to surf the SortedBasisFunctions array quickly, prepare an array with the beginning indexes of BS of i-th atom
        // example: BasisFunctionsPerAtom = 1, 3, 5, 2, ... => FirstBasisFunctionOfAtom = 0, 1, 4, 9, 11, ...
        int[] FirstSortedBasisFunctionOfAtom = new int[ nAtoms ];
        FirstSortedBasisFunctionOfAtom[0] = 0;
        for (int i=1; i<nAtoms; i++) FirstSortedBasisFunctionOfAtom[i] = FirstSortedBasisFunctionOfAtom[i-1] + BasisFunctionsPerAtom[i-1];
        // now fill in SortedBasisFunctions array
        // prepare auxiliary array of sliding indexes initialized with the values of FirstBasisFunctionOfAtom
        int[] _sliding_indexes = FirstSortedBasisFunctionOfAtom.clone();
        // finally, fill SortedBasisFunctions array in
        for (int i=0; i<InitialBasis.length; i++) {
            SortedBasisFunctions[ _sliding_indexes[InitialBasis[i].Center_ID-1]  ] = i;
            _sliding_indexes[InitialBasis[i].Center_ID-1]++;
        }

        _sliding_indexes = null;
        /* Ready.
         * It is almost for sure that at the first use of this transofrmation
         * (when Initial basis coinsides with atomic basis set functions)
         * SortedBasisFunctions will contain 0,1,2,3,4,5,....,Basis.length-1.
         * However, in subsequent calls it is expected to be quite different
         */
        
        
        /* Now create radial parts;
         * Loop over all centers,
         *  for each center loop over all of its L values,
         *      for each L collect indexes of all basis functions of that center with a given L
         *      make a list of such a basis functions and a list of their radial parts
         */
        
        int[] RadialPartsUsedForThisL = new int[MaxBasisFunctionsPerAtom/*100*/]; // array length = maximum numebr of radial parts with fixed L per atom <= MaxBasisFunctionsPerAtom
        int NRadialPartsUsedForThisL = 0; // number of valid indexes in RadialPartsUsedForThisL
        
        int[] BasisFunctionsUsedForThisL = new int[MaxBasisFunctionsPerAtom/*200*/]; // array length <= MaxBasisFunctionsPerAtom //~= maxLen(RadialPartsUsedForThisL)*(2Lmax+1)
        int NBasisFunctionsUsedForThisL = 0; // number of valid indexes in BasisFunctionsUsedForThisL
        
        // some auxiliary variables
        boolean was_found = false;
        int RID = 0;
        
        double[] PNAO_weights = new double[InitialBasis.length];
        
        int NPNAOs = 0;
        if (verbose_print)
            out.println("Producing natural orbitals for each center...");
        
        for (int cntr=0; cntr < nAtoms && InputBasisNotAmpty; cntr++) {
            for(int L=0; L <= LMaxOfAtom[cntr]; L++) {
                if (verbose_print) out.println("center = "+(cntr+1)+", L = "+L);
                // find all functions of this atom with a given L
                NRadialPartsUsedForThisL = 0;
                NBasisFunctionsUsedForThisL = 0;
                for (int i=FirstSortedBasisFunctionOfAtom[cntr]; i<(FirstSortedBasisFunctionOfAtom[cntr]+BasisFunctionsPerAtom[cntr]); i++)
                    if (InitialBasis[ SortedBasisFunctions[i] ].L == L) {
                        // A desired basis funciton has been found!
                        // Have we already met the radial part of this basis function?
                        was_found = false;
                        RID = InitialBasis[ SortedBasisFunctions[i] ].RadialPart_ID;
                        for (int j=0; (!was_found)&&(j<NRadialPartsUsedForThisL); j++)
                            was_found |= (RadialPartsUsedForThisL[j] == RID);
                        if (!was_found) {
                            RadialPartsUsedForThisL[NRadialPartsUsedForThisL] = RID;
                            NRadialPartsUsedForThisL++;
                        }
                        // simmilar for BasisFunctionsUsedForThisL
                        was_found = false;
                        for (int j=0; (!was_found)&&(j<NBasisFunctionsUsedForThisL); j++)
                            was_found |= (BasisFunctionsUsedForThisL[j] == SortedBasisFunctions[i]);
                        if (!was_found) {
                            BasisFunctionsUsedForThisL[NBasisFunctionsUsedForThisL] = SortedBasisFunctions[i];
                            NBasisFunctionsUsedForThisL++;
                        }
                    }
                //
                /* Find all functions of this atom with a given L once again;
                 * Now - build up a matrix of NRadialPartsUsedForThisL x NRadialPartsUsedForThisL size and
                 * fill it with the data from the basis functions found above
                 */
                if ((NRadialPartsUsedForThisL==0 )) { // MODIFIED
                    continue;
                }
                double[][] Local_SDS = new double[NRadialPartsUsedForThisL][NRadialPartsUsedForThisL];
                double[][] Local_S = new double[NRadialPartsUsedForThisL][NRadialPartsUsedForThisL];
                // fill in the matrices
                for (int RIndex1=0; RIndex1<NRadialPartsUsedForThisL; RIndex1++)
                    for (int RIndex2=RIndex1; RIndex2<NRadialPartsUsedForThisL; RIndex2++) {
                        Local_S[RIndex1][RIndex2] = 0;
                        Local_SDS[RIndex1][RIndex2] = 0;
                        // Loop over all pairs of basis functions used for this center and this L,
                        // and pick up those where both functions have the same m, while 1-st has RIndex1, and 2-nd has RIndex2
                        int NFns = 0;
                        for (int i=0; i<NBasisFunctionsUsedForThisL; i++)
                            for (int j=i; j<NBasisFunctionsUsedForThisL; j++) {
                                boolean same_m = ( InitialBasis[ BasisFunctionsUsedForThisL[i] ].m == InitialBasis[ BasisFunctionsUsedForThisL[j] ].m  );
                                boolean Ith_eq_R1 = ( InitialBasis[ BasisFunctionsUsedForThisL[i] ].RadialPart_ID == RadialPartsUsedForThisL[RIndex1] );
                                boolean Ith_eq_R2 = ( InitialBasis[ BasisFunctionsUsedForThisL[i] ].RadialPart_ID == RadialPartsUsedForThisL[RIndex2] );
                                boolean Jth_eq_R1 = ( InitialBasis[ BasisFunctionsUsedForThisL[j] ].RadialPart_ID == RadialPartsUsedForThisL[RIndex1] );
                                boolean Jth_eq_R2 = ( InitialBasis[ BasisFunctionsUsedForThisL[j] ].RadialPart_ID == RadialPartsUsedForThisL[RIndex2] );                                
                                if (same_m && ((Ith_eq_R1 && Jth_eq_R2) || (Ith_eq_R2 && Jth_eq_R1))) {
                                        NFns++;
                                        Local_S[RIndex1][RIndex2] += OverlapMatrix.get(BasisFunctionsUsedForThisL[i], BasisFunctionsUsedForThisL[j]);
                                        Local_SDS[RIndex1][RIndex2] += SDS.get(BasisFunctionsUsedForThisL[i], BasisFunctionsUsedForThisL[j]);
                                }
                            }
                        // NFns == (2L+1)
                        if (NFns > 0) {
                            Local_S[RIndex1][RIndex2] /= NFns;
                            Local_SDS[RIndex1][RIndex2] /= NFns;
                        }
                        // if InitialBasis coicides with atomic basis set, NFns will be equal to 2L+1
                        Local_S[RIndex2][RIndex1] = Local_S[RIndex1][RIndex2];
                        Local_SDS[RIndex2][RIndex1] = Local_SDS[RIndex1][RIndex2];
                    }

                /* Now Local_S has the same elements as AORadialPartsOverlap/DIV,
                 * where DIV is something like DIV=1(s),1(p),3(d),15(f) and depends
                 * on normalization conventions (and the number of basis functions?)
                 */
                
                //
                final  boolean GlobalSettings_debug = true;
                if (GlobalSettings_debug) {
                    out.println("Local_SDS");
                    new Matrix(Local_SDS).print(13, 7);
                    out.println("Local_S");
                    new Matrix(Local_S).print(13, 7);
                }

                /* Its time now to solve a generalized eigenvalue problem in order to
                 * find mutually normalized radial parts.
                 * The most logical way for this is to construct overlap matrix for radial directly, i.e.,
                 *
                 * // Create overlap (sub-)matrix for radial parts used to cunstruct PNAOs
                 * double[][] local_radial_S = new double[NRadialPartsUsedForThisL][NRadialPartsUsedForThisL];
                 * for (int i=0; i<NRadialPartsUsedForThisL; i++)
                 *     for (int j=i; j<NRadialPartsUsedForThisL; j++) {
                 *         local_radial_S[i][j] = AORadialPartsOverlap[ RadialPartsUsedForThisL[i] ] [ RadialPartsUsedForThisL[j]];
                 *         local_radial_S[j][i] = local_radial_S[i][j];
                 *     }
                 *
                 * and then use the following eigenvalue problem:
                 *  Matrix[] EIV = Generalized_EVD_SymmMatr(new Matrix(Local_SDS), new Matrix(local_radial_S));
                 * But the use of a separate AORadialPartsOverlap matrix can be avoided since
                 * the same overlap matrix can be obtained as:
                 *  Matrix loc_S = new Matrix(Local_S).times(1/YLM_Norm2[L][0]/(2*L+1)); // "[0]": YLM_Norm2[L][m] does not depend on m
                 */
                Matrix loc_S = new Matrix(Local_S);
                /* If YLMs are already taken to be normalized to unity, 1/YLM_Norm2[L][0] can be ommited;
                 * If not, we will need a norm of Y_LM:
                 * double [][] YLM_Norm2 = SphericalHarmonics.Get_Quick_YLM_Norm2(Quick_YLM);
                 * The origin of 1/(2L+1) factor is as follows:
                 * Each radial part enters the basis set exaclty 2L+1 times (by "cloning" methods - see MOLDEN_IO),
                 * where L is the L of basis the function, which provided us with that radial part;
                 * therefore, while summing up local_S we will meet each radial part exctly 2L+1 times
                 * loc_S = loc_S.times( 1.0/(2*L+1) );
                 * such a loc_S whould be equal to local_radial_S described above.
                 *
                 * Note that multiplication of the overlap matrix by a constant will not change eigenvectors,
                 * but WILL change their lengths (see S^(-1/2) in Generalized_EVD);
                 * in order to avoid making additional normalizations we simply modify
                 * the overlap matrix by 1/(2L+1) factor in order to get new radial
                 * parts already normalized to unity.
                 * It might be a good idea to multiply Local_SDS by 1/(2L+1) as well - this will not influence
                 * eigenvectors in any way, but WILL influence eigenvalues, and hence - weights
                 * of new basis functions. However, we choose to make 1/(2L+1) with final weights rather than
                 * with the whole matrix (see below)
                 */
                Matrix loc_SDS = new Matrix(Local_SDS);

                // Now solve the eigenvalue problem
                Matrix[] EIV = ee.Generalized_EVD_SymmMatr(loc_SDS, loc_S );
                /* Here: EIV[0] = eigenvalues
                 *       EIV[1] = eigenvectors
                 */
                double[][] eigenvecs = EIV[1].getArray();
                /* Now EIV[1] contains coefficients of expansion of radial parts
                 * of new, mutually orthogonal, basis functions over RadialPartsUsedForThisL[i]
                 * Now EACH radial part will be mulptiplied by YLM to generate 2*L+1 basis function
                 */
                // the number of newly created functions should be equal to NBasisFunctionsUsedForThisL

                /* It is important to allocate memory for ALL newly created basis functions for this center
                 * BEFORE their coefficients are calculated.
                 */
                for (int SourceBF=0; SourceBF < NBasisFunctionsUsedForThisL; SourceBF++) {
                    // Let the new funciton have the same m and radial label as the source function
                    // (however, the functions standing behing the formal radial part numbers will be different now!)
                    result[NPNAOs] = new BasisFunction(L, InitialBasis[ BasisFunctionsUsedForThisL[SourceBF] ].m, null, InitialBasis.length);
                    result[NPNAOs].Center_ID = cntr+1;//InitialBasis[ BasisFunctionsUsedForThisL[SourceBF] ].Center_ID;
                    result[NPNAOs].RadialPart_ID = InitialBasis[ BasisFunctionsUsedForThisL[SourceBF] ].RadialPart_ID; //RadialPartsUsedForThisL[RIndex1];
                    // Now we need to use an eigenvector component which corresponds to the same
                    // radial part as one used by sourceBF-th basis function =>
                    // find that "local" (i.e., valid as RadialPartsUsedForThisL[] index) number of the radial part
                    was_found = false;
                    int RIndex1 = 0;
                    while ((!was_found)&&(RIndex1<NRadialPartsUsedForThisL)) {
                        was_found |= (RadialPartsUsedForThisL[RIndex1] == InitialBasis[ BasisFunctionsUsedForThisL[SourceBF] ].RadialPart_ID);
                        if (!was_found) RIndex1++;
                    }
                    result[NPNAOs].weight = EIV[0].get(RIndex1, RIndex1);
                    // set up expansion coefs now.
                    // find all basis functions belonging to this center and having m == result[NPNAOs].m and
                    // identify their radial part number
                    for (int j=0; j<NBasisFunctionsUsedForThisL; j++)
                        if (InitialBasis[ BasisFunctionsUsedForThisL[j] ].m == result[NPNAOs].m) {
                            // determine its radial part ID (=RIndex2) and pick RIndex2-th component of RIndex1-th eigenvector
                            // as the BasisFunctionsUsedForThisL[j]-th expansion coefficient
                            was_found = false;
                            int RIndex2 = 0;
                            while ((!was_found)&&(RIndex2<NRadialPartsUsedForThisL)) {
                                was_found |= (RadialPartsUsedForThisL[RIndex2] == InitialBasis[ BasisFunctionsUsedForThisL[j] ].RadialPart_ID);
                                if (!was_found) RIndex2++;
                            }
                            //
                            result[NPNAOs].coefs[ BasisFunctionsUsedForThisL[j] ] = eigenvecs[RIndex2][RIndex1];
                            // where eigenvecs[RIndex1][A] is RIndex1-th component of A-th eigenvector
                        }
                    // Check norm
                    if (options.glbPrint) {
                        double tmp=0; // a norm
                        for (int i=0; i<InitialBasis.length; i++) {
                            out.printf("%12.5f", result[NPNAOs].coefs[i]);
                            for (int j=0; j<InitialBasis.length; j++)
                                tmp += result[NPNAOs].coefs[i] * OverlapMatrix.get(i, j) * result[NPNAOs].coefs[j];
                        }
                        out.printf("\n pnao %d: norm2 = %.7f \n", NPNAOs, tmp);
                    }
                    NPNAOs++;
                }
            } 
        }

        out.println(" Total number of natural functions produced: "+NPNAOs);

        if (verbose_print) {
            out.println(" weights of the natural functions produced:");
            double OccTotal = 0;
            for (int i=0; i<NPNAOs; i++) {
                out.printf("%15.7f", result[i].weight);
                OccTotal += result[i].weight;
            }
            out.println();
            out.printf(" sum = %12.7f%n%n",OccTotal);
        }
        
        return result;
    }
    //--------------------------------------------------------------------------
    // Generates user-readable labels
    static public String[] Create_NAO_Labels(BasisFunction[] Basis) {
        String[] result = new String[Basis.length];
        for (int i=0; i<Basis.length; i++) {
            String FnID =  String.format("R%d*%s(%d)", Basis[i].RadialPart_ID+1,
                    "spdfg".charAt(Basis[i].L), Basis[i].m);
            // special label for extended molden definitions
            if (Basis[i].additional_r_power != 0)
                FnID = FnID + String.format("*r^%d", Basis[i].additional_r_power);
            
            if (!Basis[i].NRB)
                result[i] = String.format("A%d: %s", Basis[i].Center_ID, FnID);
            else
                result[i] = String.format("A%d*: %s", Basis[i].Center_ID, FnID);
        }

        return result;
    }
    //--------------------------------------------------------------------------
    // returns an expansion coefficients of orbitals as a matrix;
    // result[m][c] is c-th coefficient of m-th orbital
    static public Matrix BasisFunctionsToMatrix(BasisFunction[] bs, boolean clone_data) {
        double[][] A = new double[bs.length][];
        for (int b=0; b<bs.length; b++)
            if (bs[b] != null) {
                if (clone_data)
                    A[b] = bs[b].coefs.clone();
                else
                    A[b] = bs[b].coefs;  // this prevents allocating memort twice when not requested explicitly
            } else
                WarningManager.warning_printf("warning: uncomplete matrix!");

        Matrix result;
        if (bs.length>0)        // Note: it CAN be a 'minimal' AO basis => NO NRB functions at all!
            result = new Matrix(A);
        else
            result = new Matrix(0,0);

        //Matrix result = new Matrix(bs.length, bs.length, 0.0);
        //for (int b=0; b<bs.length; b++)
        //    if (bs[b] != null)
        //        for (int c=0; c<bs[b].coefs.length; c++)
        //            result.set(b/*row*/, c/*column*/, bs[b].coefs[c]);
        //    else
        //        out.println("warning: uncomplete matrix!");
        return result;
    }
    //--------------------------------------------------------------------------
    // returns an array with weights of orbitals as a matrix;
    // (This routine is NOT used internally within NPA class)
    static public double[] BasisFunctionsToOccupancies(BasisFunction[] bs) {
        double[] result = new double[bs.length];
        for (int b=0; b<bs.length; b++)
            if (bs[b] != null)
                result[b] = bs[b].weight;
            else
                WarningManager.warning_printf("warning: uncomplete matrix!");
        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns off-diagonal element of the matrix A with the largest absolute value
     */
    public static double max_offdiag(Matrix A) {
        double max_offdiag = 0.0;
        double[][] arr = A.getArray();
        for (int i=0; i<arr.length; i++)
            for (int j=(i+1); j<arr[i].length; j++) {                
                if (arr[i][j] > max_offdiag) max_offdiag = arr[i][j];
                if (arr[j][i] > max_offdiag) max_offdiag = arr[j][i];
            }
        return max_offdiag;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns || A - E ||^2 where E is a unitary matrix
     */
    double NonUnitary_Norm2(Matrix A) {
        double diff_norm2 = 0;
        for (int i=0; i<A.getRowDimension(); ++i)
            for (int j=0; j<A.getColumnDimension(); j++) {
                double aij = A.get(i, j);
                if (i==j) 
                    diff_norm2 += (aij-1)*(aij-1);
                else 
                    diff_norm2 += aij*aij;
            }
        return diff_norm2;
    }
    //--------------------------------------------------------------------------
    // An array for (charge) -> (number of functions in Natural Minimal Basis) transformation
    // These data is based on a ground state electron configurations of atoms from
    // Gary L. Miessler, Donald A. Tarr - Inorganic Chemistry (Prentice Hall, 2003, p. 39)
    // First index: value of L, second index: (Z of atom-1)
    final int[][] NMB_per_Atom = new int[][]{ 

        // Number of s-type functions:
        {1, 1,  2, 2 ,2, 2, 2, 2, 2, 2,  3, 3 ,3, 3, 3, 3, 3, 3,  4, 4 ,4, 4, 4, 4, 4, 4, 4, 4, 4, 4 ,4, 4, 4, 4, 4, 4,
      /* H He| Li Be :B  C  N  O  F Ne|  Na Mg:Al Si P  S Cl Ar|  K Ca:Sc Ti  V Cr Mn Fe Co Ni Cu Zn:Ga Ge As Se Br Kr    <-- rows 1-4
                     :1-st p shell fills      :2-nd p shell fills     :1-st d shell fills           :3-rd p shell fills       */
         5, 5 ,5,  5,  5,  5,  5,  5,  5,  4,  5,  5 , 5,  5,  5,  5,  5,  5,
      /*Rb Sr:Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd :In  Sn  Sb  Te   I  Xe |    <-- row 5
             : 2-nd d shell fills (except for Pd !)  : 4-th p shell fills    |        */
         6,  6 ,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6 ,  6,  6,  6,  6,  6,  6,  6,  6,  6 , 6,  6,  6,  6,  6,  6
      /*Cs  Ba : La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu:  Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg :Tl  Pb  Bi  Po  At  Rn |    <-- row 6
               : Lanthanoids: 1-st f / 3-rd d shell fill                   :  3-rd d shell fills                 : 5-th p-shell fills    |     */
        }, 
        
        // Number of p-type functions:
        {0, 0,  0, 0 ,3, 3, 3, 3, 3, 3,  3, 3 ,6, 6, 6, 6, 6, 6,  6, 6 ,6, 6, 6, 6, 6, 6, 6, 6, 6, 6 ,9, 9, 9, 9, 9, 9,
      /* H He| Li Be :B  C  N  O  F Ne|  Na Mg:Al Si P  S Cl Ar|  K Ca:Sc Ti  V Cr Mn Fe Co Ni Cu Zn:Ga Ge As Se Br Kr    <-- rows 1-4
                     :1-st p shell fills      :2-nd p shell fills     :1-st d shell fills           :3-rd p shell fills       */
         9, 9 ,9,  9,  9,  9,  9,  9,  9,  9,  9,  9 ,12, 12, 12, 12, 12, 12,
      /*Rb Sr:Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd :In  Sn  Sb  Te   I  Xe |    <-- row 5
             : 2-nd d shell fills                    : 4-th p shell fills    |        */
        12, 12 , 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 , 12, 12, 12, 12, 12, 12, 12, 12, 12 ,15, 15, 15, 15, 15, 15
      /*Cs  Ba : La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu:  Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg :Tl  Pb  Bi  Po  At  Rn |    <-- row 6
               : Lanthanoids: 1-st f / 3-rd d shell fill                   :  3-rd d shell fills                 : 5-th p-shell fills    |     */
        }, 

        // Number of d-type functions:
        {0, 0,  0, 0 ,0, 0, 0, 0, 0, 0,  0, 0 ,0, 0, 0, 0, 0, 0,  0, 0 ,5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ,5, 5, 5, 5, 5, 5,
      /* H He| Li Be :B  C  N  O  F Ne|  Na Mg:Al Si P  S Cl Ar|  K Ca:Sc Ti  V Cr Mn Fe Co Ni Cu Zn:Ga Ge As Se Br Kr    <-- rows 1-4
                     :1-st p shell fills      :2-nd p shell fills     :1-st d shell fills           :3-rd p shell fills       */
         5, 5 ,10,10, 10, 10, 10, 10, 10, 10, 10, 10 ,10, 10, 10, 10, 10, 10,
      /*Rb Sr:Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd :In  Sn  Sb  Te   I  Xe |    <-- row 5
             : 2-nd d shell fills                    : 4-th p shell fills    |        */
        10, 10 , 15, 15, 10, 10, 10, 10, 10, 15, 10, 10, 10, 10, 10, 10, 15 , 15, 15, 15, 15, 15, 15, 15, 15, 15 ,15, 15, 15, 15, 15, 15
      /*Cs  Ba : La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu:  Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg :Tl  Pb  Bi  Po  At  Rn |    <-- row 6
               : Lanthanoids: 1-st f / 3-rd d shell fill                   :  3-rd d shell fills                 : 5-th p-shell fills    |     */
        },

        // Number of f-type functions:
        {0, 0,  0, 0 ,0, 0, 0, 0, 0, 0,  0, 0 ,0, 0, 0, 0, 0, 0,  0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0,
      /* H He| Li Be :B  C  N  O  F Ne|  Na Mg:Al Si P  S Cl Ar|  K Ca:Sc Ti  V Cr Mn Fe Co Ni Cu Zn:Ga Ge As Se Br Kr    <-- rows 1-4
                     :1-st p shell fills      :2-nd p shell fills     :1-st d shell fills           :3-rd p shell fills       */
         0, 0 ,0,  0,  0,  0,  0,  0,  0,  0,  0,  0 , 0,  0,  0,  0,  0,  0,
      /*Rb Sr:Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd :In  Sn  Sb  Te   I  Xe |    <-- row 5
             : 2-nd d shell fills                    : 4-th p shell fills    |        */
         0,  0 ,  0,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7 ,  7,  7,  7,  7,  7,  7,  7,  7,  7 , 7,  7,  7,  7,  7,  7
      /*Cs  Ba : La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu:  Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg :Tl  Pb  Bi  Po  At  Rn |    <-- row 6
               : Lanthanoids: 1-st f / 3-rd d shell fill                   :  3-rd d shell fills                 : 5-th p-shell fills    |     */
        },

        // Number of g-type functions:
        {0, 0,  0, 0 ,0, 0, 0, 0, 0, 0,  0, 0 ,0, 0, 0, 0, 0, 0,  0, 0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0,
      /* H He| Li Be :B  C  N  O  F Ne|  Na Mg:Al Si P  S Cl Ar|  K Ca:Sc Ti  V Cr Mn Fe Co Ni Cu Zn:Ga Ge As Se Br Kr    <-- rows 1-4
                     :1-st p shell fills      :2-nd p shell fills     :1-st d shell fills           :3-rd p shell fills       */
         0, 0 ,0,  0,  0,  0,  0,  0,  0,  0,  0,  0 , 0,  0,  0,  0,  0,  0,
      /*Rb Sr:Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd :In  Sn  Sb  Te   I  Xe |    <-- row 5
             : 2-nd d shell fills                    : 4-th p shell fills    |        */
         0,  0 ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ,  0,  0,  0,  0,  0,  0,  0,  0,  0 , 0,  0,  0,  0,  0,  0
      /*Cs  Ba : La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu:  Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg :Tl  Pb  Bi  Po  At  Rn |    <-- row 6
               : Lanthanoids: 1-st f / 3-rd d shell fill                   :  3-rd d shell fills                 : 5-th p-shell fills    |     */
        }, // number of g-type functions

/*
And here is some data for the next elements (Fr(87) - Usb(172)):
s    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9
p   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 15, 15, 15, 15, 15, 15, 15, 15, 15, 18, 18, 18, 18, 18, 18, 18, 18, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 24, 24
d   15, 15, 20, 20, 20, 20, 20, 15, 15, 20, 15, 15, 15, 15, 15, 15, 15, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 25, 25, 20, 20, 25, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25
f    7,  7,  7,  7, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21
g    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9  
 */
    };
    //--------------------------------------------------------------------------
    /** makes NAOs from PNAOs
     *
     */
    boolean Create_NAOs() throws Exception {

        //----------------------------------------------------------------------
        //
        // PNAO formation by intracenter naturalization of atomic basis functions
        //
        //----------------------------------------------------------------------
        out.printf(printout.Stars+"%n%n");
        out.printf("Creating NAOs%n%n");
        out.printf("%nSTEP 1. Produce PNAOs%n%n");
        
        // Create PNAOs
        PNAOs = IntracenterBasisOrthogonalization(Basis, Centers.length, SDS, OverlapMatrix);

        // Sort PNAOs by their weights to divide them into NMB and NRB subsets
        out.println("Sorting PNAOs...");

        int strted_PNAOs[] = new int[PNAOs.length]; // indexes of PNAO sorted by descending their weights
        // Sort PNAOs by weights using bubble sort
        boolean done, found;
        for (int i=0; i<PNAOs.length; i++) strted_PNAOs[i]=i;
        found = true;
        while (found) {
            found = false;
            for (int pnao=0; pnao<(PNAOs.length-1); pnao++)
                if ( PNAOs[strted_PNAOs[pnao]].weight < PNAOs[strted_PNAOs[pnao+1]].weight ) {
                    found = true;
                    // swap pnao-th and (pnao+1)-th functions
                    int tmp = strted_PNAOs[pnao];
                    strted_PNAOs[pnao] = strted_PNAOs[pnao+1];
                    strted_PNAOs[pnao+1] = tmp;
            }
        }
        
        //----------------------------------------------------------------------
        //
        // Split PNAOs into two sets: Natural Minimal Basis (NMB) and Natural Rydberg Basis (NRB)
        //
        //----------------------------------------------------------------------
        out.printf("%nSTEP 2. Split PNAOs into NMB / NRB sets%n%n");
        
        // Label eacha of those weight-sorted PNAOs as "NMB" or "NRB"
        // Use NMB_per_Atom array in that process
        out.println(" Number of basis functions in the Natural Minimal Basis (NMB) set for each center: ");
        int LMX = 5; //'G'+1
        int[] functions_used = new int[ LMX ]; // index = value of L
        // loop over all atoms
        for(int cntr=0; cntr<Centers.length; cntr++){
            // Z of this center
            int Z = new Double(Centers[cntr].Z).intValue();
            if (Z == 0) {
                WarningManager.warning_printf("Warning: center %d has zero nuclear charge; no NMB functions will be assigned to it!%n", cntr+1);
            }
            // loop over all Ls
            for (int L=0; L<LMX; L++) functions_used[L] = 0;
            // loop over all PNAOs centered on atom cntr, and determine their role: NMB / NRB
            for (int i=0; i<PNAOs.length; i++)
                if ( PNAOs[strted_PNAOs[i]].Center_ID == (cntr+1) ) {
                    PNAOs[strted_PNAOs[i]].NRB = false; // assume that this function belongs to NRB
                    functions_used[ PNAOs[strted_PNAOs[i]].L ] ++;                    
                    int nmb_LZ = 0;
                    if (Z > 0) // 13.Jan.2019: a work-around for Ghost atoms (which have Z==0 in molden files)
                        nmb_LZ = NMB_per_Atom[ PNAOs[strted_PNAOs[i]].L ][Z-1];
                    if (functions_used[ PNAOs[strted_PNAOs[i]].L ] > nmb_LZ )
                        PNAOs[strted_PNAOs[i]].NRB = true;
                }

            // make NMB / NRB statistics - get number of NMB functions for this center
            for (int L=0; L<LMX; L++) functions_used[L] = 0;
            for (int i=0; i<PNAOs.length; i++)
                if (( PNAOs[i].Center_ID == (cntr+1) ) && (!PNAOs[i].NRB))
                    functions_used[ PNAOs[i].L ] ++;
            
            // print NMB summary
            out.printf("center %3d: ",cntr+1);
            for (int L=0; L<LMX; L++)
                out.printf(" %2d of %c |", functions_used[L], "spdfg".charAt(L));
            out.println();

        }
        strted_PNAOs = null;
        // done!
        // Get total count of functions in NMB and NRB sets
        int N_NMB = 0;
        int N_NRB = 0;
        for(int i=0; i<PNAOs.length; i++)
            if (PNAOs[i].NRB) N_NRB++; else N_NMB++;
        // done!
        out.printf("In total: NMB set has %d functions, NRB set has %d functions; %n%n", N_NMB, N_NRB);
       
        // print NRB-belonging labels:
        if (verbose_print) {
            out.println(" Does PNAO belong to NRB?: ");
            for (int i=0; i<PNAOs.length; i++)
                    out.printf("%15b", PNAOs[i].NRB);
            out.println();
            /*
            // and corresponding weights
            out.println(" weights:");
            double OccTotal = 0;
            for (int i=0; i<PNAOs.length; i++) {
                out.printf("%15.7f", PNAOs[i].weight);
                OccTotal += PNAOs[i].weight;
            }
            out.println();
            out.printf(" sum = %12.7f\n",OccTotal);
             * 
             */
        }


        // Create separate lists of NMB and NRB indexes of PNAOs
        int[] original_NMB_numbers = new int[N_NMB];
        int[] original_NRB_numbers = new int[N_NRB];
        { // begin local block to make its _i and _j unvisible from the outside
            int _i=0; int _j=0;
            for(int bs=0; bs<PNAOs.length; ++bs)
                if (PNAOs[bs].NRB) {
                    original_NRB_numbers[_i] = bs;
                    ++_i;
                } else {
                    original_NMB_numbers[_j] = bs;
                    ++_j;
                }
        } // end local block

        //----------------------------------------------------------------------
        //
        // Create the Overlap and S.D.S matrices in PNAO basis
        //
        //----------------------------------------------------------------------
        if (verbose_print)
            out.println("Transforming overlap and SDS matrices to PNAO basis...");
        
        // Get AO -> PNAO transformation coefficients:
        // PNAO[i] = PNAO_to_AO[i][j]_sum(j)_AO[j]
        Matrix PNAO_to_AO_Matrix = BasisFunctionsToMatrix(PNAOs, false);


        // Create user-readable labels of PNAOs: includes PNAO_cntr_IDs, PNAO_RadialPart_GlobalIDs, L and m
        PNAO_Labels = Create_NAO_Labels(PNAOs);  

        // Transform AO overlap matrix to the PNAO basis:
        // (PNAO_mu, PNAO_nu) = SUM(i,k)( c[mu][i]·c[nu][k]·S[i,k] ) = (c.S.c^T)[mu,nu], where c[][] is PNAO_to_AO array
        PNAO_Overlap_Matrix = ee.TransformMatrixToNewBasis(OverlapMatrix, PNAO_to_AO_Matrix, true);//  (PNAO_to_AO_Matrix.times(OverlapMatrix)).times(PNAO_to_AO_Matrix.transpose());
        //        

        //debug

        // Print that matrix
        printout.Print_Matrix(PNAO_Overlap_Matrix, "PNAO overlap matrix:",
                PNAO_Labels, PNAO_Labels, options.PNAO_OverlapMatrix_File.get_String());
        
        out.printf("Trace of the PNAO overlap matrix: %.7f (should be equal to %d, the total number of PNAOs)\n",
                PNAO_Overlap_Matrix.trace(), PNAOs.length);

        // Transform global SDS matrix to PNAO basis
        Matrix SDS_PNAO = ee.TransformMatrixToNewBasis(SDS, PNAO_to_AO_Matrix, true); //(PNAO_to_AO_Matrix.times(SDS)).times(PNAO_to_AO_Matrix.transpose());
        out.println(" The trace of SDS matrix in PNAO basis = "+SDS_PNAO.trace());

        printout.Print_Matrix(SDS_PNAO, "The S.D.S matrix in PNAO basis:",
                PNAO_Labels, PNAO_Labels, options.PNAO_SDS_Matrix_File.get_String());
        



        //----------------------------------------------------------------------
        //
        // Weighted orthogonalization of NMB PNAOs
        //
        //----------------------------------------------------------------------
        out.printf("%nSTEP 3. Weighted orthogonalization of NMB PNAOs%n%n");

        double minWeight;

        // This is the transfromation matrix we are about to create;
        Matrix OW1 = new Matrix( PNAOs.length, PNAOs.length, 0.0); // W-orthogonalized NMB + unchanged NRB to PNAO transition matrix
        // It will remain NRBs unchanged:
        for (int i=0; i<N_NRB; i++)
            OW1.set(original_NRB_numbers[i], original_NRB_numbers[i], 1.0);

        // Create a transfromation for NMB:

        // Create "local" overlap and SDS matrices for new NMB functions
        Matrix NMB_Overlap_Matrix = PNAO_Overlap_Matrix.getMatrix(original_NMB_numbers, original_NMB_numbers);
        Matrix NMB_SDS_Matrix = SDS_PNAO.getMatrix(original_NMB_numbers, original_NMB_numbers);
        /*
        Matrix NMB_Overlap_Matrix = new Matrix(N_NMB, N_NMB, 0.0);
        Matrix NMB_SDS_Matrix = new Matrix(N_NMB, N_NMB, 0.0);
        for (int i=0; i<N_NMB; i++)
            for (int j=i; j<N_NMB; j++) {
                double Sij = PNAO_Overlap_Matrix.get(original_NMB_numbers[i], original_NMB_numbers[j]);
                NMB_Overlap_Matrix.set(i, j, Sij);
                NMB_Overlap_Matrix.set(j, i, Sij);
                double sdsij = SDS_PNAO.get(original_NMB_numbers[i], original_NMB_numbers[j]);
                NMB_SDS_Matrix.set(i, j, sdsij);
                NMB_SDS_Matrix.set(j, i, sdsij);
            }
         */

        // Create "local" labels
        String[] NMB_Labels = new String[N_NMB];
        for (int i=0; i<N_NMB; i++) NMB_Labels[i] = String.format("A%d: R%d*%s(%d)", PNAOs[original_NMB_numbers[i]].Center_ID,
                PNAOs[original_NMB_numbers[i]].RadialPart_ID+1,
                    "spdfg".charAt(PNAOs[original_NMB_numbers[i]].L), PNAOs[original_NMB_numbers[i]].m);

        printout.Print_Matrix(NMB_Overlap_Matrix, "NMB_old overlap:", NMB_Labels, NMB_Labels, 
                options.NMB_old_Overlap_Matrix_File.get_String());
        printout.Print_Matrix(NMB_SDS_Matrix, "NMB_old SDS:", NMB_Labels, NMB_Labels, 
                options.NMB_old_SDS_Matrix_File.get_String());
        

        // Create W.S.W^T matrix
        // due to floating-point 'features' it is better NOT to use Matrix.times method to avoid a slightly-non-symmetric result
        Matrix WSW_NMB = new Matrix(N_NMB, N_NMB, 0.0);
        minWeight = PNAOs[original_NMB_numbers[0]].weight;
        for (int i=0; i<N_NMB; i++) {
            double w_i = PNAOs[original_NMB_numbers[i]].weight;
            //
            //out.printf("%3d\t%.15f%n",original_NMB_numbers[i],w_i);
            //
            if (w_i <= minWeight) minWeight = w_i;
            if (w_i <= 1.e-7) w_i = 1.e-7; // Just safety for the inverse
            for (int j=i; j<N_NMB; j++) {
                double w_j = PNAOs[original_NMB_numbers[j]].weight;
                if (w_j <= 1.e-7) w_j = 1.e-7; // Just safety for the inverse
                double wsw_ij = NMB_Overlap_Matrix.get(i, j) * w_i * w_j;
                WSW_NMB.set(i, j, wsw_ij );
                WSW_NMB.set(j, i, wsw_ij );
            }
        }
        out.println(" min weight of NMB PNAO = "+minWeight);


        // Do a W-orthogonalization step:
        // create OW = W(WSW)^(-1/2), where S is PNAO overlap matrix
        // Begin with (WSW)^(-1/2):        
        WSW_NMB = ee.Matrix_minus05(WSW_NMB);
        //Now do: Matrix OW_NMB = W_NMB.times(Matrix_minus05(WSW_NMB)); => OW_NMB = W_NMB.WSW_NMB^(-1/2);
        // OW_NMB[i,j] = W_NMB[i,k].WSW_NMB^(-1/2)[k,j] = W_NMB[i,i].WSW_NMB^(-1/2)[i,j]
        // For good transformation OW_NMB^T should be used:
        // OW_NMB[j,i] = W_NMB[j,j].WSW_NMB^(-1/2)[j,i]
        for (int j=0; j<N_NMB; j++) {
            double w_j = PNAOs[original_NMB_numbers[j]].weight;
            for (int i=0; i<N_NMB; i++)
                //In fact we are doing something like: OW_NMB.set(i,j, w_i * WSW_NMB.get(i,j));
                OW1.set(original_NMB_numbers[i], original_NMB_numbers[j], w_j * WSW_NMB.get(j, i));
        }
        
        /*      
        // create diagonal matrix of weights        
        Matrix W_NMB = new Matrix(N_NMB, N_NMB, 0.0);
        minWeight = 2.0;
        for (int i=0; i<N_NMB; i++) {
                W_NMB.set(i, i, PNAOs[original_NMB_numbers[i]].weight);
                if ( PNAOs[original_NMB_numbers[i]].weight < minWeight) minWeight = PNAOs[original_NMB_numbers[i]].weight;
        }
        // Do a W-orthogonalization step:
        // create OW = W(WSW)^(-1/2), where S is PNAO overlap matrix
        Matrix OW_NMB = W_NMB.times(Matrix_minus05(W_NMB.times(NMB_Overlap_Matrix).times(W_NMB)));
        // coefficients for NMB -> PNAO
        for (int i=0; i<N_NMB; i++)
            for (int j=0; j<N_NMB; j++) OW1.set(original_NMB_numbers[i], original_NMB_numbers[j], OW_NMB.get(j, i));
        // coefficients for NRB -> PNAO (unitary sub matrix)
         */

        // Use OW1 to transform SDS and Overlap matrices into (NRB + updated-NMB) basis              
        Matrix Overlap_new = ee.TransformMatrixToNewBasis(PNAO_Overlap_Matrix, OW1, true); // OW1.times(PNAO_Overlap_Matrix).times(OW1.transpose());
        Matrix SDS_new = ee.TransformMatrixToNewBasis(SDS_PNAO, OW1, true); //OW1.times(SDS_PNAO).times(OW1.transpose());


        /*
        fortPrint.PrintMatrixSimple(Overlap_new,"Overlap_new","aaa",null,null);                     
        out.println("Overlap_new");
        Overlap_new.getMatrix(ordr, ordr).print(13, 7);
        this.options.dont_print_matrices = false;
        this.PrintMatrix(Overlap_new, "Overlap_new", null);
        this.PrintMatrix(SDS_new, "SDS_new", null);
         * 
         */

        //debug
        if (verbose_print) {
            out.println("Weights of NMB PNAOs:");
            double sum=0;
            for(int i=0; i<original_NMB_numbers.length; ++i) {
                sum+=SDS_new.get(original_NMB_numbers[i], original_NMB_numbers[i]);
                out.printf("%4d\t%4d\t%.15f%n", i, original_NMB_numbers[i], SDS_new.get(original_NMB_numbers[i], original_NMB_numbers[i]));
            }
            out.printf("sum of NMB weights = %.10f %n",sum);
        }

        //----------------------------------------------------------------------
        //
        // Schmidt orthogonalization of NRBs to (new) NMBs
        //
        //----------------------------------------------------------------------
        out.printf("%nSTEP 4. Schmidt orthogonalization of NRBs to new NMBs%n%n");

        // new NRB to NRB transformation
        BasisFunction[] NRB = new BasisFunction[N_NRB];


        // NRB[nrb].coefs will contain an expansion overl original NRB PNAOs
        for (int nrb=0; nrb<NRB.length; nrb++) {
            NRB[nrb] = new BasisFunction(PNAOs[ original_NRB_numbers[nrb] ].L, PNAOs[ original_NRB_numbers[nrb] ].m, null, PNAOs.length );
            NRB[nrb].Center_ID = PNAOs[ original_NRB_numbers[nrb] ].Center_ID;
            NRB[nrb].RadialPart_ID = PNAOs[ original_NRB_numbers[nrb] ].RadialPart_ID;
            NRB[nrb].coefs[ original_NRB_numbers[nrb] ] = 1.0; // initialize with unitary transformation
            NRB[nrb].NRB = true;
            //
        }

        //PNAO_Overlap_Matrix = PNAO_Overlap_Matrix.times(PNAO_Overlap_Matrix.inverse()).times(2);
        // Schmidt orthogonalize NRBs to new_NMBs
        for (int nrb=0; nrb<NRB.length; nrb++) 
                for (int nmb=0; nmb<N_NMB; nmb++)
                    NRB[nrb].coefs[ original_NMB_numbers[nmb] ] -= Overlap_new.get(original_NMB_numbers[nmb], original_NRB_numbers[nrb]);

        
        Matrix OS1 = new Matrix(PNAOs.length, PNAOs.length, 0.0);
        for (int b=0; b<N_NRB; b++)
                for (int c=0; c<PNAOs.length; c++)
                    OS1.set(original_NRB_numbers[b], c, NRB[b].coefs[c]);
        for (int b=0; b<N_NMB; b++) OS1.set(original_NMB_numbers[b], original_NMB_numbers[b],  1.0);


        Overlap_new = ee.TransformMatrixToNewBasis(Overlap_new, OS1, true);// OS1.times(Overlap_new).times(OS1.transpose());
        
/*
        Matrix _Overlap_new = OS1.times(Overlap_new).times(OS1.transpose());
        fortPrint.PrintMatrixSimple(PNAO_to_AO_Matrix, "_Overlap_new", null, ordr1, ordr1);

        //////DEBUG
        double[] n2 = new double[N_NRB];
        for (int i=0; i<N_NRB; ++i) n2[i] = _Overlap_new.get( original_NRB_numbers[i], original_NRB_numbers[i] );
        for (int b=0; b<N_NRB; b++)
                for (int c=0; c<PNAOs.length; c++)
                    OS1.set(original_NRB_numbers[b], c, NRB[b].coefs[c] /Math.sqrt(n2[b]) );
        for (int b=0; b<N_NMB; b++) OS1.set(original_NMB_numbers[b], original_NMB_numbers[b],  1.0);
        Overlap_new = OS1.times(Overlap_new).times(OS1.transpose());
        //Overlap_new = _Overlap_new;
        for (int i=0; i<N_NRB; ++i) {
            n2[i] = Overlap_new.get( original_NRB_numbers[i], original_NRB_numbers[i] );
        }
 *
 */

        /*this.options.dont_print_matrices = false;
        PrintMatrix(Overlap_new , "", PNAO_Labels);*/
        SDS_new = ee.TransformMatrixToNewBasis(SDS_new, OS1, true); //OS1.times(SDS_new).times(OS1.transpose());

        //-----------------------------------------------
        if (true) {
            out.println("Diagonal elements of SDS matrix after 1-st Schmidt transformation for new NRB functions:");
            double sum = 0;
            for(int i=0; i<original_NRB_numbers.length; ++i) {
                sum+=SDS_new.get(original_NRB_numbers[i], original_NRB_numbers[i]);
                out.printf("%4d\t%4d\t%.15f%n", i, original_NRB_numbers[i], SDS_new.get(original_NRB_numbers[i], original_NRB_numbers[i]));
            }
            out.printf("Sum of NRB 'occupancies' = %.10f %n",sum);
        }
        //----------------------------------------------------------------------
        //
        // Intracenter naturalization of new NRBs
        //
        //----------------------------------------------------------------------
        out.printf("%nSTEP 5. Intracenter naturalization of new NRBs%n%n");
        
        // The transformation we are going to create
        Matrix ON2 = new Matrix(PNAOs.length, PNAOs.length, 0.0);
        // It will not change NMBs; so, set coefficients for NMB -> PNAO to unitary sub matrix
        for (int i=0; i<N_NMB; i++)
            ON2.set(original_NMB_numbers[i], original_NMB_numbers[i], 1.0);


        String[] NRB_Labels = Create_NAO_Labels(NRB);  // user-readable labels of PNAOs: includes PNAO_cntr_IDs, PNAO_RadialPart_GlobalIDs, L and m
        
        // Create overlap and SDS matrices for new NRB functions
        Matrix NRB_Overlap_Matrix = Overlap_new.getMatrix(original_NRB_numbers, original_NRB_numbers);
        Matrix NRB_SDS_Matrix = SDS_new.getMatrix(original_NRB_numbers, original_NRB_numbers);
        /*
        Matrix NRB_Overlap_Matrix = new Matrix(N_NRB, N_NRB, 0.0);
        for (int i=0; i<N_NRB; i++)
            for (int j=0; j<N_NRB; j++) NRB_Overlap_Matrix.set(i, j, Overlap_new.get(original_NRB_numbers[i], original_NRB_numbers[j]));

        Matrix NRB_SDS_Matrix = new Matrix(N_NRB, N_NRB, 0.0); //new Matrix(SDS_NRB);
        for (int i=0; i<N_NRB; i++)
            for (int j=0; j<N_NRB; j++) NRB_SDS_Matrix.set(i, j, SDS_new.get(original_NRB_numbers[i], original_NRB_numbers[j]));
         */


        printout.Print_Matrix(NRB_Overlap_Matrix, "NRB_old overlap:", NRB_Labels, NRB_Labels, 
                options.NRB_old_Overlap_Matrix_File.get_String());

        
        // Intraatomic natural transformation within NRB set
        BasisFunction[] NRB_new = IntracenterBasisOrthogonalization(NRB, Centers.length, NRB_SDS_Matrix, NRB_Overlap_Matrix);
        // NRB_new is N_NRB x N_NRB in size
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // a small test
        String Lbl[]=Create_NAO_Labels(NRB_new);
        //for (int i=0; i<Lbl.length; i++) out.printf("%15s", Lbl[i]);
        Matrix U = BasisFunctionsToMatrix(NRB_new, false);
        
        Matrix NRB_new_Oevrlap = ee.TransformMatrixToNewBasis(NRB_Overlap_Matrix, U, true); //U.times(NRB_Overlap_Matrix).times(U.transpose());


        printout.Print_Matrix(NRB_new_Oevrlap, "NRB_new overlap:", NRB_Labels, NRB_Labels, 
                options.NRB_new_Overlap_Matrix_File.get_String());


        // coefficients for naturalized NRB -> PNAO
        for (int i=0; i<N_NRB; i++)
            for (int j=0; j<N_NRB; j++) ON2.set(original_NRB_numbers[i], original_NRB_numbers[j], NRB_new[i].coefs[j]);
        
        Overlap_new = ee.TransformMatrixToNewBasis(Overlap_new, ON2, true); //ON2.times(Overlap_new).times(ON2.transpose());
        printout.Print_Matrix(Overlap_new, "Overlap_new", PNAO_Labels, PNAO_Labels, 
                options.S_Matrix_after_ON2_File.get_String());
        
        SDS_new = ee.TransformMatrixToNewBasis(SDS_new, ON2, true); //ON2.times(SDS_new).times(ON2.transpose());
        printout.Print_Matrix(SDS_new, "SDS_new", PNAO_Labels, PNAO_Labels, 
                options.SDS_Matrix_after_ON2_File.get_String());


        if (verbose_print) {
            out.println(" NMB and NRB weights and occupancies after intracenter naturalization of NRBs");
            double sum=0;
            out.println("LocNum.\tGlobNum.\t<phi|SDS|phi>");
            for(int i=0; i<original_NMB_numbers.length; ++i) {
                sum += SDS_new.get(original_NMB_numbers[i], original_NMB_numbers[i]);
                out.printf("%4d\t%4d\t%.10f%n",
                        i,
                        original_NMB_numbers[i],
                        SDS_new.get(original_NMB_numbers[i], original_NMB_numbers[i]));
            }
            out.printf("Sum of NMB SDS diagonal terms = %.10f %n",sum);

            sum = 0;
            out.println("LocNum.\tGlobNum.\t<phi|SDS|phi>\tweight\tNRB label");
            for(int i=0; i<original_NRB_numbers.length; ++i) {
                sum += SDS_new.get(original_NRB_numbers[i], original_NRB_numbers[i]);
                //NRB_new[i].weight = SDS_new.get(original_NRB_numbers[i], original_NRB_numbers[i]);
                out.printf("%4d\t%4d\t%.10f\t%.10f\t%s%n",
                        i,
                        original_NRB_numbers[i],
                        SDS_new.get(original_NRB_numbers[i], original_NRB_numbers[i]),
                        NRB_new[i].weight,
                        PNAO_Labels[original_NRB_numbers[i]]);
            }
            out.printf("Sum of NRB SDS diagonal terms = %.10f %n",sum);
        }

        //----------------------------------------------------------------------
        //
        // Weighted orthogonalization of naturalized NRBs
        //
        //----------------------------------------------------------------------
        out.printf("%nSTEP 6. Weighted orthogonalization of naturalized NRBs%n%n");

        // The transformation matrix we are about to create:
        Matrix OW2 = new Matrix(PNAOs.length, PNAOs.length, 0.0);
        // It will not change NMBs. So, set the coefficients for NMB -> PNAO to unity
        for (int i=0; i<N_NMB; i++)
            OW2.set(original_NMB_numbers[i], original_NMB_numbers[i], 1.0);
        //
        Matrix OW_NRB = null;
        
        boolean NRB_OW_Straightforward = options.NRB_OW_Straightforward.get_boolean();
        double Heavy_NRB_Threshold = options.Heavy_NRB_Threshold.get_double();

        if (NRB_OW_Straightforward) {
            if (verbose_print) out.println("Using a direct WSW orthogonalization for all NRB functions.");
            // Algorithm 1: Direct occupancy-weighted orthogonalization of ALL NRB orbitals
            Matrix W_NRB = new Matrix(N_NRB, N_NRB, 0.0);
            if (N_NRB > 0) // it IS possible that there are no NRB functions!
                minWeight = NRB_new[0].weight;
            else
                minWeight = 0.0;
            
            for (int i=0; i<N_NRB; ++i) {
                    W_NRB.set(i, i, NRB_new[i].weight); // NRB_new2[original_NRB_numbers[i]].weight);
                    if ( NRB_new[i].weight < minWeight) minWeight = NRB_new[i].weight;
            }
            out.println(" min weight = "+minWeight);

            // create OW = W(WSW)^(-1/2), where S is PNAO overlap matrix
            //OW_NRB = W_NRB.times(Matrix_minus05(W_NRB.times(NRB_new_Oevrlap).times(W_NRB)));
            Matrix WSW_NRB = new Matrix(N_NRB, N_NRB, 0.0);
            for (int i=0; i<N_NRB; ++i)
                for (int j=0; j<N_NRB; ++j) {
                    double tmp = NRB_new[i].weight * NRB_new_Oevrlap.get(i, j) * NRB_new[j].weight;
                    WSW_NRB.set(i, j, tmp);
                }
            OW_NRB = W_NRB.times(ee.Matrix_minus05(WSW_NRB)); // yes, it is not very efficient,
                // but this algorithm is not to be used very often at all!
        } else {
            // Algorithm 2
            if (verbose_print) 
                out.printf("Using a WSWSL orthogonalization for NRB functions with "+
                        "occ. Threshold = %.3E%n",Heavy_NRB_Threshold);
            
            //   0. Get max weight
            double maxWeight = 0;
            for (int i=0; i<N_NRB; ++i)
                if (NRB_new[i].weight > maxWeight) maxWeight = NRB_new[i].weight;
            out.printf("Maximum weight of NRB function = %.5E%n", maxWeight);
            
            //   1. Direct occupancy-weighted orthogonalization for 'heavily' occupied NRBs
            int nHeavy = 0;
            int nLight = 0;
            // bugfix for N_NRB == 0
            if (N_NRB > 0)
                minWeight = NRB_new[0].weight;
            else
                minWeight = 0.0;
            for (int i=0; i<N_NRB; ++i) {
                if (NRB_new[i].weight > Heavy_NRB_Threshold) 
                    ++nHeavy;
                else
                    ++nLight;
                //
                if (NRB_new[i].weight < minWeight) minWeight = NRB_new[i].weight;
            }
            out.println(" min weight = "+minWeight);
            if (verbose_print) out.printf("'Heavily occupied NRB set' has %d of %d functions%n", nHeavy, N_NRB);
            
            int[] Heavy_NRBs = new int[nHeavy]; // Note: these arrays contain NOT a global indexes of orbitals,
            int[] Light_NRBs = new int[nLight]; // but the indexes valid within the NRB range only!
            Matrix W_NRB_Heavy = new Matrix(nHeavy, nHeavy, 0.0);
            nHeavy = 0;
            nLight = 0;
            for (int i=0; i<N_NRB; ++i)
                if (NRB_new[i].weight > Heavy_NRB_Threshold) {
                    Heavy_NRBs[nHeavy] = i;
                    W_NRB_Heavy.set(nHeavy, nHeavy, NRB_new[i].weight);
                    ++nHeavy;
                } else Light_NRBs[nLight++] = i;
            
            
            Matrix HeavyNRB_Overlap = NRB_new_Oevrlap.getMatrix(Heavy_NRBs, Heavy_NRBs);
            // W.(W.S.W)^(-1/2)
            Matrix WSW_nrb_heavy = new Matrix(nHeavy, nHeavy, 0.0);//W_NRB_Heavy.times(HeavyNRB_Overlap).times(W_NRB_Heavy);
            // 'manual' matrix multiplication to avoid non-symmetry caused by round-off errors
            for (int i=0; i<nHeavy; ++i)
                for (int j=i; j<nHeavy; ++j) {
                    double tmp = NRB_new_Oevrlap.get(Heavy_NRBs[i], Heavy_NRBs[j]);
                    tmp *= NRB_new[Heavy_NRBs[i]].weight * NRB_new[Heavy_NRBs[j]].weight;
                    WSW_nrb_heavy.set(i, j, tmp);
                    WSW_nrb_heavy.set(j, i, tmp);
                }            
            Matrix OW_NRB_Heavy = W_NRB_Heavy.times(ee.Matrix_minus05(WSW_nrb_heavy));
            /*
            // new overlap matrix of 'heavy' NRBs - to be used in Schmidt orthogonalization process
            Matrix HeavyNRB_Overlap_new = OW_NRB_Heavy.transpose().times(HeavyNRB_Overlap.times(OW_NRB_Heavy));
             */
            // insert OW_NRB_Heavy into the final OW_NRB transformation matrix
            OW_NRB = new Matrix(N_NRB, N_NRB, 0);
            OW_NRB.setMatrix(Heavy_NRBs, Heavy_NRBs, OW_NRB_Heavy);
            // initialize OW_NRB diagonal positions for 'light' NRBs with ones
            for (int i=0; i<nLight; ++i) OW_NRB.set(Light_NRBs[i], Light_NRBs[i], 1.0);
            // update full(!) NRB overlap matrix
            Matrix new_NRB_Overlap = ee.TransformMatrixToNewBasis(NRB_new_Oevrlap, OW_NRB.transpose(), true); // OW_NRB.transpose().times(NRB_new_Oevrlap.times(OW_NRB));
            printout.Print_Matrix(new_NRB_Overlap, "new_NRB_Overlap before Schmidt", null, null, 
                    options.NRB_Overlap_after_OW_heavy_File.get_String());
                        
                        
            //   2. Schmidt orthogonalization for the rest of NRBs w.r.t. the 'heavily' occupied ones
            // substract projections onto NEW(!) 'heavily occupied' NRBs
            Matrix OS = new Matrix(N_NRB, N_NRB, 0); // Schmidt transformation matrix
            for (int i=0; i<N_NRB; ++i) OS.set(i, i, 1.0); // initialize with 1-matrix

            for (int i=0; i<nLight; ++i) // loop over 'light' orbitals
                for (int j=0; j<nHeavy; ++j) {
                    double tmp = OS.get(Light_NRBs[i], Heavy_NRBs[j]);
                    tmp -= new_NRB_Overlap.get(Light_NRBs[i], Heavy_NRBs[j]);
                    OS.set(Light_NRBs[i], Heavy_NRBs[j], tmp);
                }
            // Note that OS is NOT symmetric!
            OW_NRB = OW_NRB.times(OS.transpose());
                    // OW^T = OS.OW^T => OW = OW.OS^T
            // update NRB overlap matrix
            new_NRB_Overlap = ee.TransformMatrixToNewBasis(NRB_new_Oevrlap, OW_NRB.transpose(), true) ; //OW_NRB.transpose().times(NRB_new_Oevrlap.times(OW_NRB));
            
            printout.Print_Matrix(new_NRB_Overlap, "new_NRB_Overlap after Schmidt", null, null, 
                    options.NRB_Overlap_after_OS2_File.get_String());
            
            //   3. Lowdin orthogonalization for these 'rest' NRBs within themselves
            if (nLight > 0) {
                Matrix Light_Overlap = new_NRB_Overlap.getMatrix(Light_NRBs, Light_NRBs);                
                Matrix OL = ee.Matrix_minus05(Light_Overlap);
                Matrix OLNRB = new Matrix(N_NRB, N_NRB, 0); // unchanged 'healy' NRBs and Lowdin transformation for 'light' NRBs
                for (int i=0; i<N_NRB; ++i) OLNRB.set(i, i, 1.0);
                OLNRB.setMatrix(Light_NRBs, Light_NRBs, OL);
                
                // insert Lowdin step into global NRB transformation
                OW_NRB = OW_NRB.times(OLNRB.transpose()); 
            }
        }


        // test the orthonormality of new NRBs
        Matrix nrb_overlap_test = ee.TransformMatrixToNewBasis(NRB_new_Oevrlap, OW_NRB.transpose(), true); //OW_NRB.transpose().times(NRB_new_Oevrlap.times(OW_NRB));

        printout.Print_Matrix(nrb_overlap_test, "NRB new overlap", null, null, 
                options.NRB_Overlap_after_OW2_final_File.get_String());
        
        out.printf(" |S_NRB - 1| = %.5E (Should be VERY close to zero!)%n", Math.sqrt(NonUnitary_Norm2(nrb_overlap_test )));


        OW2.setMatrix(original_NRB_numbers, original_NRB_numbers, OW_NRB.transpose() );


        //Matrix A = OW2.times(ON2).times(OS1).times(OW1).times(PNAO_to_AO_Matrix);

        String[] AO_names = Create_NAO_Labels(Basis);
        printout.Print_Matrix(OW2, "2-nd WSW transformation matrix:", null, AO_names, options.OW2_File.get_String());

        Overlap_new = ee.TransformMatrixToNewBasis(Overlap_new, OW2, true); //OW2.times(Overlap_new).times(OW2.transpose());
        printout.Print_Matrix(Overlap_new, "Overlap matrix in almost-NAO basis (obtained after OW2 step):",
                PNAO_Labels, PNAO_Labels, options.S_Matrix_after_OW2_File.get_String());
        
        SDS_new = ee.TransformMatrixToNewBasis(SDS_new, OW2, true); // OW2.times(SDS_new).times(OW2.transpose());
        printout.Print_Matrix(SDS_new, "SDS matrix in almost-NAO basis (obtained after OW2 step):",
                PNAO_Labels, PNAO_Labels, options.SDS_Matrix_after_OW2_File.get_String());
        


        if (verbose_print) {
            out.println("Diagonal elements of SDS before the final intracenter naturalization step:");
            for (int cntr=0; cntr<Centers.length; cntr++) {
                double cntr_occ = 0;
                for (int i=0; i<Basis.length; i++)
                    if ((PNAOs[i].Center_ID-1) == cntr) cntr_occ+=SDS_new.get(i, i);
                out.printf("center %d: %.4f\n", cntr+1, cntr_occ);
            }
        }

        //----------------------------------------------------------------------
        //
        // Final Intracenter Natural Transformation
        //
        //----------------------------------------------------------------------
        out.printf("%nSTEP 7. Final Intracenter Natural Transformation withing the full set of functions%n%n");

        options.dont_print_matrices = false;

        // Perform final transformation step
        NAO = IntracenterBasisOrthogonalization(PNAOs, Centers.length, SDS_new, Overlap_new);

        // Create corresponding transformation matrix
        Matrix BF_new_2_NAO = BasisFunctionsToMatrix(NAO, false);
        
        // Create the PNAO->NAO transformation matrix
        Matrix NAO_2_PNAO = BF_new_2_NAO.times(OW2).times(ON2).times(OS1).times(OW1);
        // Create the NAO<-AO basis functions transformation matrix
        NAO_2_AO = NAO_2_PNAO.times(PNAO_to_AO_Matrix);

        // Transform ovarlap matrix to NAO basis
        Overlap_NAO = ee.TransformMatrixToNewBasis(Overlap_new, BF_new_2_NAO, true); //BF_new_2_NAO.times(Overlap_new).times(BF_new_2_NAO.transpose());
        
        //PrintMatrix(Overlap_NAO, " Ovrlap in NAO basis:", PNAO_Labels );

        // Compare this matrix with a diagonal one-matrix (it should be equal!)
        double diff_norm2 = NonUnitary_Norm2(Overlap_NAO);
        out.printf(" SQRT{ SUM[(NaoOverlap_ij - delta_ij)^2] } = %.2e (should be VERY close to zero) %n",Math.sqrt(diff_norm2));
        out.printf(" max_offdiag = %.2e (should be VERY close to zero) %n", max_offdiag(Overlap_NAO) );
        if (Math.sqrt(diff_norm2) > 1.0e-10 )
            WarningManager.warning_printf(" WARNING: NAOs seem to be not strictly orthogonalized!");

        // Transform S.D.S matrix to NAO basis
        SDS_NAO = ee.TransformMatrixToNewBasis(SDS_new, BF_new_2_NAO, true); //BF_new_2_NAO.times(SDS_new).times(BF_new_2_NAO.transpose());
        // Set NAO weights to diagonal terms of S.D.S matrix in NAO basis
        out.println(printout.Stars);
        out.printf("%n Final NAO occupancies and leading AO terms:%n%n");
        double[] _old_weights = new double[NAO.length];
        out.printf("%5s %20s %10s %40s%n", "NAO #", "Name", "Occupancy", "Leading term");
        for (int i=0; i<NAO.length; i++) {
            _old_weights[i] = NAO[i].weight;
            NAO[i].weight = SDS_NAO.get(i, i);
            NAO[i].NRB = PNAOs[i].NRB;              // WARNING: this should be clarified !!! - perhaps, we should give this label by taking top-N of NAOs for this center...
            out.printf("%5d %20s %10.7f ", i+1, PNAO_Labels[i], NAO[i].weight);
            // find a major contribution AO
            int AO_main = 0;
            for(int j=0; j<Basis.length; ++j)
                if (Math.abs(NAO_2_AO.get(i, j)) > Math.abs(NAO_2_AO.get(i, AO_main))) AO_main = j;
            double tmp = NAO_2_AO.get(i, AO_main);
            out.printf("%40s", String.format("(%.2f)*BF[%d = %s]", tmp, AO_main+1, AO_Names[AO_main]));
            out.printf("%n");
        }
      
        printout.Print_Matrix(SDS_NAO, " S.D.S in NAO basis:", PNAO_Labels, PNAO_Labels, options.SDS_NAO_File.get_String());
        printout.Print_Matrix(NAO_2_AO, "AO-to-NAO transformation matrix:", PNAO_Labels, AO_Names, options.NAO2AO_File.get_String());

        

        out.println(" trace = " + SDS_NAO.trace() );

        // Perform population analysis
        out.printf("%nFinal electron populations and NPA charges:%n%n");
        out.println("Center\tNuclear\t Electron  \t  NMB       \tNPA   ");
        out.println("      \t charge\t population\t  population\tcharge");
        NPA_charges = new double[Centers.length];
        //double[] NMB_contribs = new double[Centers.length];
        for (int cntr=0; cntr<Centers.length; cntr++) {
            double cntr_occ = 0;
            double nmb_contrib = 0;
            for (int i=0; i<Basis.length; i++)
                if ((NAO[i].Center_ID-1) == cntr) {
                    cntr_occ += NAO[i].weight;
                    if (!NAO[i].NRB) nmb_contrib += NAO[i].weight;
                }//  SDS_NAO.get(i, i);
            // save charge to an array
            NPA_charges[cntr] = Centers[cntr].Z - cntr_occ;
            out.printf("%5s\t%7.1f\t%11.7f\t%12.7f\t%13.10f\n",
                    String.format("%s%d", Centers[cntr].Name,cntr+1),
                    Centers[cntr].Z,
                    cntr_occ,
                    nmb_contrib,
                    Centers[cntr].Z - cntr_occ);
        }
        out.println();

        // Some additional statistics: angular momentum contributions into the population:
        out.printf("Angular momentum contributions of the total atomic population:%n%n");
        out.print("   Cntr");
        double[] L_contrib = new double[SphericalHarmonics.L_MAX+1];
        for (int l=0; l<L_contrib.length; ++l) out.printf("%12s","spdfghijklmnopqrs".charAt(l));
        out.println();
        
        for (int cntr=0; cntr<Centers.length; cntr++) {
            for (int l=0; l<L_contrib.length; ++l) L_contrib[l] = 0;
            // collect contributions with a given L for this center from all possible NAOs
            for (int i=0; i<Basis.length; i++)
                if ((NAO[i].Center_ID-1) == cntr)
                    L_contrib[ NAO[i].L ] += NAO[i].weight;
            // now print the data
            out.printf("%7s", String.format("%s%d", Centers[cntr].Name,cntr+1));
            for (int l=0; l<L_contrib.length; ++l) out.printf("%12.7f", L_contrib[l]);
            out.println();
        }
        out.println();

        // that's all, folks! :)
        return true;
    }
    //--------------------------------------------------------------------------
    /**
     * Computes Wiberg-Mayer bond indexes (see [K. Wiberg, Tetrahedron 1968, V.24, P.1083-1096]
     * and [ I. Mayer "Bond orders and valences from ab initio wave functions" /
     * Int. J. Quant. Chem. 1986, V. 29, P. 477–483] for more details).
     * 
     * main input is: S05DS05 = S^(1/2) . D . S^(1/2)
     * 
     * @return a matrix (nAtoms x nAtoms) containing bond indexes between the pairs of atoms
     * result[A][B] = SUM(S05DS05[i,j]*S05DS05[j,i] for: i-th NAOs on A and j-th NAO on B)
     * Diagonal elements are set to the "Actual total valence" of i-th atom = SUM[j]( result[i][j], j =/= i )
     *
     */
    // Rev. 16.May.2016: old BondIndexes() generalized into BondIndexes_Generic()
    
    public static Matrix BondIndexes_Generic(int Natoms, BasisFunction[] basisFunctions, Matrix S05DS05) {        
        Matrix result = new Matrix(Natoms, Natoms, 0.0);
        double[][] arr = result.getArray();
        for (int i=0; i<basisFunctions.length; i++)
            for (int j=(i+1); j<basisFunctions.length; j++)
                arr[basisFunctions[i].Center_ID-1]
                        [basisFunctions[j].Center_ID-1] += S05DS05.get(i, j)*S05DS05.get(j, i);
        // Make the matrix symmetric
        for (int i=0; i<Natoms; i++)
            for (int j=(i+1); j<Natoms; j++) arr[j][i] = arr[i][j];        
        // Fill the diagonal with the sum of off-diagonal elenemts for this atom:
        for (int i=0; i<Natoms; i++) {
            double v = 0;
            arr[i][i] = 0;
            for (int j=0; j<Natoms; j++) v += arr[i][j]; // arr[i][i] == 0, so it will not spoil the sum
            arr[i][i] = v;
        }

        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Redirects to BondIndexes_Generic() using the data taken from this.NAO, this.Centers.length, and  this.SDS_NAO
     * @returns the Wiberg bond indices (as defined by Wiberg) using the density matrix in NAO basis
     */
    public Matrix BondIndexes() {
        return BondIndexes_Generic(Centers.length, NAO, SDS_NAO);
    }
    //--------------------------------------------------------------------------
    public void MutualDensities() {

        out.printf("NAO-atoms contributions to the densities @ nuclei:%n");
        for (int c=0; c<Centers.length; ++c) {            
            double[] r = Centers[c].R0;
            // evaluate all basis functions @this point
            double[] BS_values = new double[Basis.length];
            for (int bf=0; bf<Basis.length; bf++) {
                Basis[bf].Quick_YLM  = Quick_YLM;
                BS_values[bf] = Basis[bf].EvaluateAtPoint(r, Basis[bf].R0);
            }
            double[] nao_values = new double[NAO.length];
            for (int nao=0; nao<NAO.length; nao++) {
                nao_values[nao] = 0;
                for (int bf=0; bf<Basis.length; bf++)
                    nao_values[nao] += NAO_2_AO.get(nao, bf) * BS_values[bf];                            
            }
            
            
            out.printf("Nucleus %3d:", c+1);
            //
            double rho_sum = 0;
            
            // loop over all atoms and evaluate their densities @r
            for (int c2=0; c2<Centers.length; c2++) {
                double rho2 = 0;
                // loop over all NAOs of c2-th atom                
                for (int nao=0; nao<NAO.length; nao++) {
                    if (NAO[nao].Center_ID == (c2+1)) {
                        /*double naoValue = 0;                        
                        for (int bf=0; bf<Basis.length; bf++)
                            naoValue += NAO_2_AO.get(nao, bf) * BS_values[bf];                            
                        rho2 += NAO[nao].weight * naoValue*naoValue;*/
                        for (int bf=0; bf<Basis.length; bf++)
                            if (NAO[bf].Center_ID == (c2+1)) 
                                rho2 += SDS_NAO.get(nao, bf) * nao_values[nao] * nao_values[bf];
                    }
                }

                // the whole true density @r:
                /*
                rho2 = 0;
                for (int nao=0; nao<NAO.length; nao++) {
                    for (int bf=0; bf<Basis.length; bf++)
                            rho2 += SDS_NAO.get(nao, bf) * nao_values[nao] * nao_values[bf] ;
                }*/
                out.printf("\t%.5f", rho2);
                rho_sum += rho2;
            }
            out.printf("\t|\t%.5f%n",rho_sum);
            
        }
        
    }

}
