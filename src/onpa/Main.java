package onpa;

/**
 * The JANPA program - an independent open-source implementation
 * of Natural Orbitals Analysis algorithms with Java platform
 * supplemented with advanced density-partitioning facilities
 *
 * Version: 13.Jan.2019
 * Created: 26.Oct.2013
 * Copyright (c) T.Yu. Nikolaienko, L.A.Bulavin, D.M.Hovorun
 *
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


import CLPO.*;
import JGints.*;
import Jama.*;
import java.math.*;
import java.util.*;
import java.io.*;
import moldenio.*;
import ProgramOptions.*;
//
import Polynom3D.*;
import MatrixHelper.*;

//==============================================================================
/*class GlobalSettings {
    static final boolean debug = false;
    final static int L_MAX = 10;
    final static boolean no_normalization = true;
}*/
//==============================================================================
//==============================================================================
/**
 * A main class
 * @author timn
 */
public class Main {    

    static String Version = "2.02 (13-01-2019)";
    static PrintStream out = System.out;
    //--------------------------------------------------------------------------
    static void Print_NMB_Submatrix(NPA npa, Matrix A) {
        // header
        out.printf("%15s", "orbital#");
        for (int i=0; i<npa.NAO.length; i++)
            if (!npa.PNAOs[i].NRB)
                out.printf("%15s", npa.PNAO_Labels[i]);
        out.println();
        for (int i=0; i<npa.NAO.length; i++) {
            if (!npa.PNAOs[i].NRB) {
                out.printf("%4d %11s", i+1, npa.PNAO_Labels[i]);
                for (int j=0; j<npa.NAO.length; j++)
                    if (!npa.PNAOs[j].NRB)
                        out.printf("%15.7f", A.get(i, j));
                out.println();
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Prints banner
     */
    private static void ShowBanner() {
        out.println(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        out.println(" *    janpa: A cross-platform open-source implementation of NPA        * ");
        out.println(" *    and other electronic structure analysis methods with Java        * ");
        out.println(" *        A part of JANPA package,   http://janpa.sourceforge.net      * ");
        out.printf (" *                    Version: %20s                    * %n",Version);
        out.println(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        out.println();
        out.println(" (c) Tymofii Nikolaienko, 2014-2019");
        out.println();
        out.println(" If any results obtained with this program are published,");
        out.println(" or for any other reasons, please, cite this work as: ");
        out.println(" 1) T.Y.Nikolaienko, L.A.Bulavin; Int. J. Quantum Chem. (2019), ");
        out.println("    Vol.119, page e25798, DOI: 10.1002/qua.25798");
        out.println(" 2) T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem.(2014),");
        out.println("    Vol.1050, pages 15-22, DOI: 10.1016/j.comptc.2014.10.002");
        printout.PrintStars();
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    private static void print_bond_indices(Matrix indices, AtomicCenter[] centers, String exportFName) throws Exception {
        int nAtoms = centers.length;
        out.println("Wiberg-Mayer bond indices (based on density matrix in NAO basis):");
        out.printf("%10s","Centr. A/B");
        for (int i=0; i<nAtoms; ++i) out.printf("%10d",i+1);
        out.println();
        for (int i=0; i<nAtoms; ++i) {
            out.printf("%7d   ",i+1); // row number
            for (int j=0; j<i; ++j) out.printf("%10s","");  // empty space
            out.printf(" (%7.4f)",indices.get(i, i));  // diagonal elements
            for (int j=(i+1); j<nAtoms; ++j)
                out.printf("%10.4f",indices.get(i, j));   // off-diaginal bond indices
            out.println();
        }        

        // Export these bond orders if requested:        
        if (!exportFName.isEmpty()) {
            String[] atomLabels = new String[nAtoms];
            for (int c=0; c<nAtoms; c++)
                atomLabels[c] = String.format("%s%d", centers[c].Name, c+1);
            printout.Print_Matrix(indices, 
                    "Wiberg-Mayer bond indices computed in NAO basis (note: "+
                    "diagonal elements are atomic "+
                    "'valencies' (sums of all bond orders in current line)):", 
                    atomLabels, atomLabels, exportFName);            
        }
        
    }
    //--------------------------------------------------------------------------
    // A field holding all adjustable program options
    private static ono_options Options = new ono_options();    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {        
        ProgTimer timer = new ProgTimer();


        Locale.setDefault(Locale.US);
        String fname = null;

        // show banner
        ShowBanner();

        // some manual correction of parameters
        Options.Heavy_NRB_Threshold.default_double_format = "%.2E";

        
        

        
        
        //JGints x = new JGints();
        

        if (args.length < 1) {
            out.println("Usage: ");
            out.println(" java -jar janpa.jar input-molden-file");
            out.println("or ");
            out.println(" java -jar janpa.jar -i input-molden-file [other possible options]");
            out.println();
            out.println("where possible options are: ");
            out.println();
            Options._FirstOption.PrintDescriptions(null, " %s:  %s (default: %s)%n%n");
            printout.PrintStars();
            return;
        }
        if (args.length == 1) {
            // there is exactly one command line argument => this is molden file name
            fname = args[0];
            // store this fname for proper printing an option list
            Options.Input_Molden_File.ReadDataFromString(fname);
        } else {
            // Command line is complicated and needs to be parsed
            Options._FirstOption.LoadOptionsFromCommandLine(args, null, true, true);
            // A parameter file could be set => get its name
            String param_file = Options.Parameter_File.get_String();
            if (!param_file.isEmpty()) 
                Options._FirstOption.LoadOptionsFromFile(new BufferedReader(new FileReader(param_file)), null, true);
            // some options could have been changed as compared to the command line => reparse the command line
            Options._FirstOption.LoadOptionsFromCommandLine(args, null, true, true);
            // finally, get the file name
            fname = Options.Input_Molden_File.get_String();
        }
        if ((fname == null) || (fname.isEmpty())) {
            out.println("Error: no file name of the input molden file provided.");
            return;
        }
        if (!(new File(fname).isFile())) {
            out.printf("ERROR: \"%s\" can not be opened or is not a file.%n", fname);
            return;
        }
        
        // Ok, all parameters are fine!
        // Print current program options:
        out.printf("Settings used in this run:%n%n");
        Options._FirstOption.Print_as_List(null, " %s: %s%n");
        printout.PrintStars();
        

        out.println("Loading MOLDEN from "+fname);

        MOLDEN_IO molden = new MOLDEN_IO();
        molden = new MOLDEN_IO();
        molden.Allow_additional_r_power = true;

        if (molden.Load_From_MOLDEN(fname)) {
            out.println();
            out.println("Data loaded successfully");            
            out.printf(" Number of basis functions: %d; number of molecular orbitals: %d %n", molden.Basis.length, molden.MOs.length);
            out.println();
        } else {
            out.println("Error loading data");
            return;
        }
        if (!molden.IsSpherical) {
            out.println(" ERROR: Input molden file does not have a 'pure' ('spherical') basis set!");
            out.println(" Try using molden2molden -cart2pure to convert it.");
            return;
        }

        //molden.Compress();

        if (!Options.Edges.get_String().isEmpty()){
            out.println();
            out.println("User defined Edges : " + Options.Edges.get_String());
        }

        NPA npa = new NPA(Options); // create NPA 'workhorse' and transfer the options to this class

        EigenEngine ee = new EigenEngine();
        //ee.MatrixTransform_Force_Symmetric = Options.MatrixTransform_Force_Symmetric;
        //ee.No_Enhanced_Matrix_Transform = Options.No_Enhanced_Matrix_Transform;
        ee.glbPrint = Options.glbPrint;
        
        //npa.dont_print_matrices = false;
        molden.CoordsToAU();

        if (Options.PrintGeometry.get_boolean()) {
            out.printf("%nCartesian coordinates of the atoms (in atomic units):%n");
            out.printf("%s\t%s\t%s\t%s\t%s\t%s%n","ID","Element","Nucl.Chrg.","     X     ","     Y     ","     Z     ");
            for (int c=0; c<molden.Centers.length; ++c) {
                out.printf("%d\t%s\t%.1f\t%.10f\t%.10f\t%.10f%n",
                        c+1,
                        molden.Centers[c].Name,
                        molden.Centers[c].Z,
                        molden.Centers[c].R0[0], molden.Centers[c].R0[1], molden.Centers[c].R0[2]);
            }
            printout.PrintStars();
        }
        
        // convert to a convenient spherical harmonics
        Polynom3D[][] ylm = SphericalHarmonics.Get_Quick_YLM_Norm4PI();// .Get_Quick_YLM();

        npa.Quick_YLM = ylm;

        if (!Options.Overlap_Naive.get_boolean()) {
            molden.ToUnnormalizedPrimitiveCoefs(); // include 'normalizing constants' implicitly present
            // before each individual primitive Gaussian into the corresponding conefficient within the radial part

            // TODO: check if further molden exports work Ok now !
            // TODO: check if options.Overlap_Naive  works properly now !            
        } else {
            double[][] ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm);
            molden.UnNormalizePrimitives(ylm_norms2);
        }

        // Import MO data from the loaded molden file, and compute the necessary integrals
        boolean molden_accepted = npa.Load_MO_From_MOLDEN(molden);
        if (!molden_accepted) {
            out.println(" MOLDEN data is not acceptable!");
            return;
        }

        

        /*JGints a = new JGints();
        a.ImportBasisFromMolden(molden);
        a.ylm_norms2 = ylm_norms2;
        a.Build_Ovarlap_Matrix_PURE();
        npa.OverlapMatrix = new Matrix(a.OverlapMatrix);
        npa.Basis = molden.Basis;
        npa.Centers = molden.Centers;*/


        
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        // THE KEY POINT: GENERATE AO->NAO TRANSFORMATION
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        npa.Create_NAOs();
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *

        printout.PrintStars();
        
        // Computes Mayer bond indexes [ I. Mayer, Int. J. Quant. Chem. 1986, V. 29, P. 477–483]
        //Matrix MayerIndices = npa.BondIndexes_Generic(npa.Centers.length, npa.Basis, npa.S05DS05);   
        //print_bond_indices(MayerIndices, npa.Centers, "123.txt");
 
        // use NAOs to get Wiberg-Mayer bond indexes
        Matrix WibergIndices = npa.BondIndexes();
        print_bond_indices(WibergIndices, npa.Centers, Options.WiebergBondOrders_File.get_String());
                
        printout.PrintStars();
        out.println("==========================================================================================================");
        PropertyOptimizedOrbitals LPOs = new PropertyOptimizedOrbitals();        
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        // THE 2-ND KEY POINT: CREATE LPOs AND CLPOs
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *

        LPOs.createCLPOs(npa.SDS_NAO, npa.NAO, npa.Centers, Options);
        
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        
        printout.PrintStars();

        out.println("Atomic connectivity analysis based on CLPO bonding graph:");
        out.println();
        int[] clpoFragmentIds = LPOs.LOconnectivity(LPOs.CLPO_descript, npa.NPA_charges);
        printout.PrintStars();

        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        
        
        // perform charge-transfer analysis in CLPO basis:
        (new CT_Estimator(LPOs.CLPO_descript, npa.SDS_NAO)).new_CT(clpoFragmentIds); ;//  npa.SDS_NAO, LPOs.CLPO_descript, clpoFragmentIds);
        
        

        
        
        
                
        
        



              
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *
        // EXPORT RESULTS
        //* * * * * * * * * * * * * * * * * * * * * * * * * * *

        // export plain list of NPA charges if needed
        if (!Options.Charges_File.get_String().isEmpty()) {
            String npach = Options.Charges_File.get_String();
            out.println("Writing NPA charges to\t"+npach);
            PrintWriter ch = new PrintWriter(new FileOutputStream(npach));
            for (int i=0; i<npa.Centers.length; ++i) ch.printf("%.16f%n", npa.NPA_charges[i]);
            ch.close();
        }

        
        // Make RADIAL PARTS Molnde-style-normalized before export!!!
        
         molden.ToNormalizedPrimitiveCoefs();// .UnNormalizePrimitives(ylm_norms2); 
         // That was OK with radial parts, but eralier we've 'tuned' these radial
         // parts to be used with a _specifically_ normalized spherical
         // harmonics. So, now we have to cancel this back
/* // older code (was used with JGints, not JGintsCyl)
         Polynom3D[][] ylm_15 = SphericalHarmonics.Get_Quick_YLM();
         double[][] ylm_norms_15 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm_15);

         for (int rp=0; rp<molden.RadialParts.length; ++rp)
         	for (int cf=0; cf<molden.RadialParts[rp].Coefs.length; cf++)
         		molden.RadialParts[rp].Coefs[cf] *= Math.sqrt(4*Math.PI * ylm_norms_15[molden.RadialParts[rp].LUsedWith][0]);
         // Ready.
         // NOTE however that we have changed only radial parts, but NOT molden.Basis
*/

         // for now, we postpone PNAO/NAO export until Fock matrix MIGHT be
         // available so that we could use its diagonal elements as the
         // orbital energies

        // prepare PNAOs
        Matrix PNAO_2_AO = npa.BasisFunctionsToMatrix(npa.PNAOs, false);
        double[] PNAO_weights = npa.BasisFunctionsToOccupancies(npa.PNAOs);
        
        //*******************************************************
        boolean do_Fock = false;

        if (molden.MOs.length < molden.Basis.length) {
            out.printf(" No information about virtual molecular orbitals was found%n"
                    + " Fock matrix will not be analyzed.%n");
        } ;//else do_Fock = true;
        //
        Matrix Fock_NAO = null;
        Matrix Fock_AO = null;
        Matrix Fock_PNAO = null;
        
        if (Options.do_Fock.get_boolean()) do_Fock = true;

        printout.PrintStars();
        
        if (do_Fock) {
            String[] NAO_Labels = npa.Create_NAO_Labels(npa.NAO);

            out.println("Fock matrix analysis in NAO basis");
            out.println();
            Matrix MO_Energies = new Matrix(molden.MOs.length, molden.MOs.length, 0.0);
            for (int i=0; i<molden.MOs.length; i++) MO_Energies.set(i, i, molden.MOs[i].Energy);
            // MO_Energies = Fock matrix in MO basis;
            // Convert it to AO basis

            out.println("Building MO -> AO matrix...");/*
            Matrix MO_to_AO = new Matrix(molden.MOs.length, molden.Basis.length, 0.0);
            for (int mo=0; mo<molden.MOs.length; mo++)
                for (int bs=0; bs<molden.Basis.length; bs++) MO_to_AO.set(mo, bs, molden.MOs[mo].BS_Coefs[bs]); */
            Matrix MO_to_AO = new Matrix(molden.get_MO2AO_array(false));
            
            out.println("Inverting MO -> AO matrix...");
            Matrix AO_to_MO = MO_to_AO.inverse();
            out.println("Transforming Fock matrix to AO basis...");
            
            Fock_AO = ee.TransformMatrixToNewBasis(MO_Energies, AO_to_MO, true); //AO_to_MO.times(MO_Energies).times(AO_to_MO.transpose());
            printout.Print_Matrix(Fock_AO, "Fock matrix in AO basis", npa.AO_Names, npa.AO_Names, Options.Fock_AO_File.get_String());
            
            
            out.println("Transforming Fock matrix to NAO basis...");
            Fock_NAO = ee.TransformMatrixToNewBasis(Fock_AO, npa.NAO_2_AO, true); // npa.NAO_2_AO.times(Fock_AO).times(npa.NAO_2_AO.transpose());

            // export this matrix if requested            
            printout.Print_Matrix(Fock_NAO, "Fock matrix in NAO basis", NAO_Labels, NAO_Labels, Options.Fock_NAO_File.get_String());
            
            //
            out.printf("NAO\t          Label\tHost\toccup.\t<NAO|F|NAO>%n");
            out.printf("   \t               \tatom\t      \t           %n");
            for (int i=0; i<npa.NAO.length; i++)
                out.printf("%d\t%15s\t%d\t%.7f\t%.7f%n",
                        i+1, NAO_Labels[i],
                        npa.NAO[i].Center_ID, npa.NAO[i].weight, Fock_NAO.get(i, i));
            //
            
            Fock_PNAO = ee.TransformMatrixToNewBasis(Fock_AO, PNAO_2_AO, true);//PNAO_2_AO.times(Fock_AO).times(PNAO_2_AO.transpose());

            if (Options.PrintNBMSubmatrices.get_boolean()) {
                // NMB elements of the Fock Matrix
                out.println("Fock submatrix in NMB PNAO basis");
                Print_NMB_Submatrix(npa, Fock_PNAO);

                out.println("PNAO overlap submatrix in NMB PNAO basis");
                Print_NMB_Submatrix(npa, npa.PNAO_Overlap_Matrix);

                // NMB elements of the Fock Matrix
                out.println("Fock submatrix in NMB NAO basis");
                Print_NMB_Submatrix(npa, Fock_NAO);
            }
        }
        //*******************************************************


        // save PNAOs
        if (!Options.PNAO_Molden_File.get_String().isEmpty()) {
            out.println("Writing PNAOs to\t" + Options.PNAO_Molden_File.get_String());
            // use Fock_PNAO matrix (if available) to set PNAO energies
            double[] PNAO_energies = null;
            if (Fock_PNAO != null) {
                out.println("Note: diagonal elements of the Fock matrix in PNAO "
                        +"baiss will be used as PNAO energies");
                PNAO_energies = new double[npa.PNAOs.length];                
                for (int pnao = 0; pnao<npa.PNAOs.length; ++pnao)
                    PNAO_energies[pnao] = Fock_PNAO.get(pnao, pnao);
            }
            printout.Export_Orbitals(molden, PNAO_2_AO, 
                    PNAO_weights, PNAO_energies, 
                    Options.PNAO_Molden_File.get_String(), false, "PNAOs" );
        }
        
        // Note: Options.NAO2AO_File is processed inside NPA class
        if (!Options.PNAO2AO_File.get_String().isEmpty()) {
            printout.Print_Matrix(PNAO_2_AO, "Pre-NAOs in AO basis", npa.PNAO_Labels, npa.AO_Names, Options.PNAO2AO_File.get_String());
        }
        
        // save NAOs
        double[] NAO_weights = npa.BasisFunctionsToOccupancies(npa.NAO);
        if (!Options.NAO_Molden_File.get_String().isEmpty()) {
            out.println("Writing NAOs to\t" + Options.NAO_Molden_File.get_String());
            // use Fock_NAO matrix (if available) to set NAO energies
            double[] NAO_energies = null;
            if (Fock_NAO != null) {
                out.println("Note: diagonal elements of the Fock matrix in NAO "
                        +"baiss will be used as NAO energies");
                NAO_energies = new double[npa.NAO.length];
                for (int nao = 0; nao<npa.NAO.length; ++nao)
                    NAO_energies[nao] = Fock_NAO.get(nao, nao);
            }

            printout.Export_Orbitals(molden, npa.NAO_2_AO, 
                    NAO_weights, NAO_energies, 
                    Options.NAO_Molden_File.get_String(), false, "NAOs" );
        }

        //*******************************************************
        
        // Hybridized/localized orbital exports (matrix format):        
        if (! Options.CLPO2LHO_File.get_String().isEmpty()) {
            // Note: in CLPO2HO: CLPO2HO[CLPO_index, LHO_index]            
            //printout.Print_Matrix(LPOs.CLPO_to_LHO, "CLPOs in LHO basis", LPOs.LHO_labels, LPOs.CLPO_labels, Options.CLPO2LHO_File.get_String());            
            printout.Print_Matrix(LPOs.CLPO_descript.LO_to_Hybrids, "CLPOs in LHO basis", LPOs.CLPO_descript.Hybrid_labels, LPOs.CLPO_descript.LO_labels, Options.CLPO2LHO_File.get_String());
        }
        if (! Options.LHO2NAO_File.get_String().isEmpty()) {
            //LPOs.NAO2LHO.transpose() [ LHO, NAO ]
            //printout.Print_Matrix(LPOs.NAO_to_LHO.transpose(), "LHOs in NAO basis",  npa.PNAO_Labels, LPOs.LHO_labels, Options.LHO2NAO_File.get_String());            
            printout.Print_Matrix(LPOs.CLPO_descript.NAO_to_Hybrids.transpose(), "LHOs in NAO basis",  npa.PNAO_Labels, LPOs.CLPO_descript.Hybrid_labels, Options.LHO2NAO_File.get_String());
        }

        // and the similar things for LPOs and AHOs
        if (! Options.LPO2AHO_File.get_String().isEmpty()) {
            //printout.Print_Matrix(LPOs.LPO_to_AHO, "LPOs in AHO basis", LPOs.LHO_labels, LPOs.LPO_labels, Options.LPO2AHO_File.get_String());            
            printout.Print_Matrix(LPOs.LPO_descript.LO_to_Hybrids, "LPOs in AHO basis", LPOs.LPO_descript.Hybrid_labels, LPOs.LPO_descript.LO_labels, Options.LPO2AHO_File.get_String());
        }
        if (! Options.AHO2NAO_File.get_String().isEmpty()) {
            //printout.Print_Matrix(LPOs.NAO_to_AHO.transpose(), "AHOs in NAO basis",  npa.PNAO_Labels, LPOs.LHO_labels, Options.AHO2NAO_File.get_String());            
            printout.Print_Matrix(LPOs.LPO_descript.NAO_to_Hybrids.transpose(), "AHOs in NAO basis",  npa.PNAO_Labels, LPOs.LPO_descript.Hybrid_labels, Options.AHO2NAO_File.get_String());
        }


        
        //*******************************************************

        
        
        // Hybridized/localized orbital exports (molden format):
        Matrix[] orbs2nao         = new Matrix[4]; 
        String[] orbs2ao_fnames   = new String[4];
        String[] orbs2ao_comments = new String[4];
        // [0]=LHO, [1]=CLPO, [2] = AHO, [3] = LPO
        final int orbs_LHO_id  = 0;
        final int orbs_CLPO_id = 1;
        final int orbs_AHO_id  = 2;
        final int orbs_LPO_id  = 3;
        
        String lho_molden_fname = Options.LHO_Molden_File.get_String();
        if (! lho_molden_fname.isEmpty() ) {
            orbs2nao[orbs_LHO_id] = LPOs.CLPO_descript.NAO_to_Hybrids.transpose();
            orbs2ao_fnames[orbs_LHO_id] = lho_molden_fname;
            orbs2ao_comments[orbs_LHO_id] = "LHOs";
        }
        
        String clpo_molden_fname = Options.CLPO_Molden_File.get_String();
        if (! clpo_molden_fname.isEmpty() ) {
            // form CLPO->NAO transformation matrix
            orbs2nao[orbs_CLPO_id] = LPOs.CLPO_descript.LO_to_Hybrids.times(LPOs.CLPO_descript.NAO_to_Hybrids.transpose());
            orbs2ao_fnames[orbs_CLPO_id] = clpo_molden_fname;
            orbs2ao_comments[orbs_CLPO_id] = "CLPOs";
        }
        
        String aho_molden_fname = Options.AHO_Molden_File.get_String();
        if (! aho_molden_fname.isEmpty() ) {
            orbs2nao[orbs_AHO_id] = LPOs.LPO_descript.NAO_to_Hybrids.transpose();
            orbs2ao_fnames[orbs_AHO_id] = aho_molden_fname;
            orbs2ao_comments[orbs_AHO_id] = "AHOs";
        }
        
        String lpo_molden_fname = Options.LPO_Molden_File.get_String();
        if (! lpo_molden_fname.isEmpty() ) {
            orbs2nao[orbs_LPO_id] = LPOs.LPO_descript.LO_to_Hybrids.times(LPOs.LPO_descript.NAO_to_Hybrids.transpose());
            orbs2ao_fnames[orbs_LPO_id] = lpo_molden_fname;
            orbs2ao_comments[orbs_LPO_id] = "LPOs";
        }

        // perform the molden export in 'bulk'
        
        for (int i=0; i<orbs2nao.length; i++) {
            if (orbs2nao[i] == null)
                continue;
            // else:     
            // convert X->NAO matrix into X->AO matrix
            Matrix orbs2ao = orbs2nao[i].times(npa.NAO_2_AO);            
            Matrix orbs_d = orbs2nao[i].times(npa.SDS_NAO);
            // d_orbs = orbs2nao[i].times(d_orbs);
            int nNAOs = npa.NAO.length;
            double[] occs = new double[nNAOs];            
            for (int j=0; j<nNAOs; j++) {
                // orbs2nao[i][j,k] * npa.SDS_NAO[k,q] * orbs2nao[i][j,q] == orbs_d[j,q] * orbs2nao[i][j,q]
                for (int k=0; k<nNAOs; k++)
                    occs[j] += orbs_d.get(j, k) * orbs2nao[i].get(j, k);
            }
            
            out.printf("Writing %s to\t%s%n",  orbs2ao_comments[i], orbs2ao_fnames[i] );
            out.printf("Please, note that 'Energy' of the orbitals in MOLDEN file will simply be set to their sequential numbers!%n%n");
            
            double[] ene  = new double[nNAOs];
            for (int j=0; j<nNAOs; j++) {
                //occs[j] =  d_orbs.get(j, j);
                ene[j] = j;
            }
            printout.Export_Orbitals(molden, orbs2ao, occs, ene,  orbs2ao_fnames[i], true, orbs2ao_comments[i] );
        }
        
                
        
        //*******************************************************
        
        
        timer.Stop();
        // summarize warnings        
        WarningManager.summarizeWarnings();
        // Report a total run time
        timer.Print();
        
        printout.PrintStars();

    }

}
