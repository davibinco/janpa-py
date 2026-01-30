package nwchem2molden;

import java.io.*;
import moldenio.*;
import JGints.*;


/**
 * A class for parsing NwChem basis set library file
 * 
 * Ver.: 14.Feb.2014
 * Created 05.Jan.2014
 * @author timn
 * (c) Tymofii Nikolaienko, 2014
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
public class BasisSetLoader {

    // first index = Z (1-based!)
    // second index = # of basis function for this element
    public BasisFunction[][] TheBasis = new BasisFunction[ PeriodicTableData.all_elements.length +1][];

    public boolean print_loaded_basis = false; // whether to print what has been loaded
    //--------------------------------------------------------------------------
    // loads basis set data
    // and stores it in TheBasis array
    // if do_sort is true, basis functions are sorted as S,P,D,...
    public void LoadFromNwFile(String fname, boolean do_sort) throws Exception{
        FileReader fr1 = new FileReader(fname);
        BufferedReader bsfile = new BufferedReader(fr1);
        //
        String s;

        // 1-based indexing for simplicity
        BasisFunction[] first_bf_of_element = new BasisFunction[ PeriodicTableData.all_elements.length +1]; // java initializes it with NULLs
        BasisFunction[] last_bf_of_element = new BasisFunction[ PeriodicTableData.all_elements.length +1];
        int current_Z = -1;
        int current_L = -2;
        final int MAX_BS_COEFS = 1000;

        s = bsfile.readLine();
        if (s == null) return;

        boolean finished = false;
        
        while (s != null /*!finished*/) {
            // skip this line if it is empty, contains no data or starts with a '#' mark or 'basis', 'end' words
            s = s.trim();
            boolean do_skip = s.isEmpty();
            if (!do_skip) do_skip = s.startsWith("#");
            if (!do_skip) do_skip = s.startsWith("basis");
            if (!do_skip) do_skip = s.startsWith("end");
            if (!do_skip) {
                // new radial part begins
                final String bf_hdr_fmt = "[a-zA-Z]+ +[SPDFGH]+ *[ 0-9\\.]*"; // ' *[10\\.]*' will skip something like '3 1.00       0.000000000000' for gaussian-like style
                if (s.matches(bf_hdr_fmt)) {
                    String[] prms = s.split(" +");
                    int Z;
                    current_Z = PeriodicTableData.Name_to_Z( prms[0] ); // get Z of this element
                    // get L
                    current_L = "SPDFGH".indexOf(prms[1]);
                    if (prms[1].equals("SP")) current_L = -1; // special case of "SP" shells
                    //
                    double[][] bf_data = new double[MAX_BS_COEFS][];
                    int nLines = 0;
                    //
                    // read in all data about this BF
                    boolean local_finished = false;
                    local_finished = ((s = bsfile.readLine()) == null);
                    while (!local_finished ) {
                        s = s.trim();
                        local_finished = !s.matches("[DE0-9\\.\\-\\+\\,]+ +[DE0-9\\.\\-\\+\\,]+.*"); //at least two numbers
                        if (!local_finished) {
                            // this line still belongs to the basis function definition
                            String[] vals = s.split(" +");
                            bf_data[nLines] = new double[vals.length];
                            for (int i=0; i<vals.length; i++) {                                
                                bf_data[nLines][i] = Double.parseDouble(vals[i].toUpperCase().replace("D", "E") ); // 0.458878D-03 -> 0.458878E-03
                            }
                            nLines++;
                            local_finished = ((s = bsfile.readLine()) == null);
                        }
                    }
                    // reading block finished!
                    // add it to the 'database'
                    int L;
                    if (current_L>0)  L = current_L; else L = 0; // 'S' first

                    // allocate memory 
                    if (first_bf_of_element[current_Z] == null)
                        // create new 'first' function if needed
                        last_bf_of_element[current_Z] = ( first_bf_of_element[current_Z] = new BasisFunction(L, 0, null, nLines) );
                    else
                        // update the 'last' pointer
                        last_bf_of_element[current_Z] = ( last_bf_of_element[current_Z]._next = new BasisFunction(L, 0, null, nLines) );

                    BasisFunction bf = last_bf_of_element[current_Z]; // a short reference
                    bf.Center_ID = current_Z;
                    // copy coefs and exponents
                    for (int i=0; i<nLines; i++) {
                        bf.exponents[i] = bf_data[i][0];
                        bf.coefs[i] = bf_data[i][1];
                    }
                    // special cases:
                    if (current_L == (-1)) {
                        // 1) SP
                        bf = (bf._next = new BasisFunction(bf)); // copies exponents and their number
                        bf._next = null;
                        last_bf_of_element[current_Z] = bf;
                        bf.L = 1; // correct L;
                        // copy coefs.
                        for (int i=0; i<nLines; i++)
                            bf.coefs[i] = bf_data[i][2];
                    }
                    // non-SP case with many sets of coefs.
                    if ((current_L != (-1)) && (bf_data[0].length > 2)) {
                        // 2) many coefs for the same exponents
                        // create many an array of clones
                        BasisFunction[] bf_clones = new BasisFunction[bf_data[0].length-2];
                        // create objects
                        for (int i=0; i<bf_clones.length; i++) bf_clones[i] = new BasisFunction(bf); // clone all data
                        // rewrite coefs.
                        for (int i=0; i<nLines; i++)
                            for (int j=0; j<bf_clones.length; j++)
                                bf_clones[j].coefs[i] = bf_data[i][j+2];
                        // correct their ._next fields
                        bf._next = bf_clones[0];
                        for (int j=0; j<(bf_clones.length-1); j++)
                            bf_clones[j]._next = bf_clones[j+1];
                        bf_clones[bf_clones.length-1]._next = null;
                        last_bf_of_element[current_Z] = bf_clones[bf_clones.length-1];
                    }

                    //System.out.println(nLines);
                    //finished = (s == null);
                } else {
                    System.out.printf(" Unknown type of line: '%s'%n",s);
                    System.out.println("skipping...");
                    do_skip = true;
                }
            }
            if (do_skip)
                s = bsfile.readLine();
            
        }        
        //
        bsfile.close();
        fr1.close();
        //
        System.out.println("reading finished.");
        //if (!do_sort) return;  // ready!

       
        
        //
        System.out.println("Sorting...");

        // 1-based numbers!
        for (int Z=1; Z<PeriodicTableData.all_elements.length; Z++) {
            BasisFunction first_bf = first_bf_of_element[Z];
            int nBFs = 0;
            int LMAX = 0;
            BasisFunction bf = first_bf;
            while (bf != null) {
                // Added Dec.25, 2017:
                // Some NwChem's basis set files have coefficient lines like following
                // "       0.785710E+04           0.568250E-03           0.000000E+00           
                //      ....
                //         0.822020E+00           0.000000E+00           0.216890E+00      "
                // i.e., have zero contraction coefs.
                // In this case NwChem itself prints a warning and ignores such primitive Gaussians
                // Thus, let's 'clean-up' the created basis functions and remove such
                // primitives too
                bf.removeSmallContractionCoefs();
                
                
                if (bf.L > LMAX) LMAX = bf.L;                
                nBFs++;
                bf = bf._next;
            }

            // create it even if length == 0  !
            TheBasis[Z] = new BasisFunction[nBFs];

            if (first_bf != null) {
                int b = 0;
                if (do_sort) {
                        // re-sort functions
                        for (int L = 0; L<=LMAX; L++) {
                            // go through list and pick up only functions with a given L
                            bf = first_bf;
                            while (bf != null) {
                                if (bf.L == L) {
                                    TheBasis[Z][b] = bf;
                                    b++;
                                }
                                bf = bf._next;
                            }
                        }
                } else {
                    // if no sorting is required
                    // just copy everything to the array
                    bf = first_bf;
                    while (bf != null) {
                        TheBasis[Z][b] = bf;
                        b++;
                        bf = bf._next;
                    }
                }
            } // if (first_bf != null)
            
        } // for all Z
        //

        // print the result
        if (print_loaded_basis) {
            System.out.println("Here is what has been loaded:");
            for (int Z=1; Z<PeriodicTableData.all_elements.length; Z++) {
                System.out.printf("Element: %2s (Z = %3d); %4d basis functions%n", PeriodicTableData.all_elements[Z-1], Z, TheBasis[Z].length);
                for (int b=0; b<TheBasis[Z].length; b++) {
                    for (int c=0; c<TheBasis[Z][b].coefs.length; c++)
                        System.out.printf("%3d %c  %15.7f  %12.8f%n", b+1, "SPDFGH".charAt(TheBasis[Z][b].L),
                                TheBasis[Z][b].exponents[c], TheBasis[Z][b].coefs[c]);
                    System.out.println();
                }
                System.out.println();
            }
        }
        


    }
    //--------------------------------------------------------------------------

}
