package onpa;

import JGints.MO;
import JGints.SphericalHarmonics;
import java.io.*;
import moldenio.*;
import Jama.*;
import ProgramOptions.WarningManager;


/**
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
 * (c) Tymofii Nikolaienko, 2014
 *
 */
public class printout extends PrintStream {

    public printout(OutputStream user_console) {
        super(user_console);
    }
    
    //--------------------------------------------------------------------------
    public static String Stars="                          * * *                                         ";

    public static void PrintStars() {
        System.out.println();
        System.out.println(Stars);
        System.out.println();
    }    
    //--------------------------------------------------------------------------
    /** Converts Orbital-to-AO matrix into a molden file
     *  // 04.09.2018: @param renormBFcoefs  added
     *  // Note that renormBFcoefs should typically be set to False if JGintsCyl integrals are used
     *  // 13.10.2018: @param title added
    */
    public static void Export_Orbitals(MOLDEN_IO AO_info, Matrix Orbitals_to_AO, 
            double[] occupancies, double[] Energies /*can be null*/,
            String FileName, boolean renormBFcoefs,
            String title) throws Exception {
        // export orbitals to the new MOLDEN
        MOLDEN_IO molden2 = new MOLDEN_IO();
        molden2.CopyAllFrom(AO_info);

        molden2.IsSpherical = true;
        if (title != null)
            molden2.Title = title;
        else
            molden2.Title = "Natural Orbitals prepared by JANPA";
        //
        if (renormBFcoefs) {
            // renormalize to spherical functions typical for MOLDEN
            double[][] YLM_Norm2 = SphericalHarmonics.Get_Quick_YLM_Norm2( SphericalHarmonics.Get_Quick_YLM() );
            double[] L_remormalizers = new double[ YLM_Norm2.length ];
            for (int L=0; L<L_remormalizers.length; L++) L_remormalizers[L] = Math.sqrt(YLM_Norm2[L][0]);
            // correct AO basis coefs
            for (int rp=0; rp<molden2.RadialParts.length; rp++)
                for (int cf=0; cf<molden2.RadialParts[rp].Coefs.length; cf++)
                    molden2.RadialParts[rp].Coefs[cf] /= L_remormalizers[ molden2.RadialParts[rp].LUsedWith ];
            // the correction for MO coefs is implemented below
        }
        
        // re-create molden2.MOs according to the data in Orbitals_to_AO matrix
        
        // bugfix 20.05.2015: the total number of orbitals and basis function 
        // coefficients is now defined by dimensionality of Orbitals_to_AO
        // matrix
        // Bugfix 04.09.2018: nBFs and nOrbitals intercganged (was not important when they are equal, which is most typically the case)
        int nBFs = Orbitals_to_AO.getColumnDimension();    // 1-st index of Orbitals_to_AO
        int nOrbitals = Orbitals_to_AO.getRowDimension();  // 2-nd index of Orbitals_to_AO
        molden2.MOs = new MO[nOrbitals];
        
        // perform some checks for the wrong number of orbitals/basis functions
        if (occupancies != null && nOrbitals != occupancies.length) {            
            WarningManager.warning_printf("WARNING - internal fail in printout.Export_Orbitals: "
                    + "Orbitals_to_AO.getColumnDimension() != occupancies.length");
        }
        if (Energies != null && nOrbitals != Energies.length) {
            WarningManager.warning_printf("WARNING - internal fail in printout.Export_Orbitals: "
                    + "Orbitals_to_AO.getColumnDimension() != Energies.length");
        }
        if (molden2.Basis.length != nBFs) {
            WarningManager.warning_printf("WARNING - internal fail in printout.Export_Orbitals: "
                    + "molden2.Basis.length != nBFs");
        }

        
        for (int mo=0; mo<nOrbitals; mo++) {
            molden2.MOs[mo] = new MO();
            molden2.MOs[mo].BS_Coefs = new double[nBFs];

            if ( occupancies != null)
                molden2.MOs[mo].Occupancy = occupancies[mo];
            //if (Fock_NAO != null)  molden2.MOs[mo].Energy = Fock_NAO.get(mo, mo);

            molden2.MOs[mo].Energy = 0;
            if (Energies != null) molden2.MOs[mo].Energy = Energies[mo];    //  added Aug.2014
            
            // copy coefs from Orbitals_to_AO matrix
            for (int cf=0; cf<nBFs; cf++) {
                //double correction = L_remormalizers[ molden2.Basis[cf].L ];
                molden2.MOs[mo].BS_Coefs[cf] = Orbitals_to_AO.get(mo, cf) ;/// correction;
            }
        }
        //
        molden2.Save_As_MOLDEN(FileName);
    }
    //--------------------------------------------------------------------------
    
    // prints the header of the matrix in the form of:
    // BasisFunctionID ( AtomName AtomNumber, BasisFunctionL / BasisFunctionM )
    static private void _Print_Matrix_Title_AO_basis(int width) {
        /*
        for (int i=0; i<NBas; i++)
            System.out.printf("%"+width+"s",
                String.format("%d[%s%d, %c(%d)]",
                    i+1, molden.Centers[molden.Basis[i].Center_ID-1].Name, molden.Basis[i].Center_ID,
                    "spdfg".charAt(molden.Basis[i].L), molden.Basis[i].m)
            );
        System.out.println();
         *
         */
        System.out.println(" <<<< _Print_Matrix_Title_AO_basis >>>>");
    }
    //--------------------------------------------------------------------------
    static public void _Print_Matrix_With_Header(Matrix A, String Title) {
            // print title
            System.out.println(Title);
            // print header (colums(=row) names)
            _Print_Matrix_Title_AO_basis(15);
            // print matrix
            A.print(13, 7);
    }
    //--------------------------------------------------------------------------
    static String MatrixFloatNumberFormat = "%12.5f"; // default value
    static int MatrixLineWidth = 0; // a maximum number of values per line for matrix export (0 = no limit)
    //--------------------------------------------------------------------------
    /** Exports the elements of matrix A ti the file @fname
     * 
     * if fname != null, prints the matrix into those file
     * 
     */
    static public void Print_Matrix(Matrix A, String Comment, String[] ColumnNames, String[] RowNames, String fname) throws Exception {
        if ((fname == null) ||(fname.isEmpty())) return;

        System.out.println();
        System.out.println("Exporting matrix \""+Comment+"\" to "+fname);

        PrintWriter out = new PrintWriter(fname);

        if (Comment != null) out.println(Comment); else out.println("");        // 1-st line: comment/title (or empty)
        out.printf("%d\t%d%n", A.getRowDimension(), A.getColumnDimension());    // 2-nd line: matrix size
        if (ColumnNames == null) out.println();                                 // 3-rd line: column names (or empty)        
        else {
            for (int i=0; i<ColumnNames.length; ++i) {
                if ((MatrixLineWidth != 0) && i > 0 && (i % MatrixLineWidth == 0))
                    out.println();
                out.printf("%s\t",ColumnNames[i]);
            }
            out.println();
        }
                                                                                // 4-th line: empty always!
        // now print the matrix with row name (or "") at the end
        double[][] a = A.getArray();
        for (int i=0; i<a.length; ++i) {
            for (int j=0; j<a[i].length; ++j) {
                if ((j == 0) || ((MatrixLineWidth != 0) && ((j % MatrixLineWidth) == 0))) 
                    out.println(); // line break before the 1-st line
                    // as well as before each MatrixLineWidth-th element
                out.printf(MatrixFloatNumberFormat, a[i][j]);
                out.print("\t");    // always separate columns with a tab
            }
            if (RowNames == null) 
                out.print("\t"); 
            else 
                out.printf("\t%s",RowNames[i]); // upd.: 06.09.2014
        }
        
        //A.print(out, 13, 7);
        out.close();
        System.out.println("done");
    }
    //--------------------------------------------------------------------------


}
