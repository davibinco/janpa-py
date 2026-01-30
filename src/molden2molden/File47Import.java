package molden2molden;

import java.io.*;
import Jama.*;
import MatrixHelper.*;
import moldenio.*;
import JGints.*;
import ProgramOptions.WarningManager;


/**
 * A class for producing the set of Lodwin natural orbitals from the data taken
 * from the nbo3-compatible .47 file
 *
 * (c) Tymofii Nikolaienko, 2014-2015
 * Rev. 24.Jan.2015
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
public class File47Import {

    //--------------------------------------------------------------------------
    // returns result[0] = extracted data; result[1] = rest of the string
    private String[] _Extract_Token(String s, String regexp) {
        if (!s.matches(regexp))
            return new String[]{s};
        else {
            return new String[]{  };
        }
    }
    //--------------------------------------------------------------------------
    public int natoms = 0;
    public int nbas = 0;
    // some possible keywords of the $GENNBO section
    public boolean bodm = false;
    public boolean upper = false;
    public boolean bohr = false;
    // from $CONTRACT section
    public int nshell = 0;
    public int nexp = 0;
    //--------------------------------------------------------------------------
    // from $DENSITY - obligatory!
    public double[][] Density;
    // from $OVERLAP - obligatory!
    public double[][] Overlap;
    // from $BASIS - obligatory!
    public int[] basis_center;
    public int[] basis_label;
    // from $CONTRACT - OPTIONAL!!!
    public int[] ncomp, nprim, nptr;
    public double[] c_exp, c_CS, c_CP, c_CD, c_CF;
    //--------------------------------------------------------------------------
    // global indices for Density/Overlap/... matrix/array elements:
    private int[] _MAT_ind  = new int[14]; // [0] = Density_i, [1] = Density_j, [2] = Overlap_i, [3] = Overlap_j, [4] = bais_center_i, [5] = basis_label_i, ...
    private final static int _D_i = 0;
    private final static int _D_j = 1;
    private final static int _S_i = 2;
    private final static int _S_j = 3;
    // arrays filled by $BASIS section
    private final static int _bas_cntr = 4;
    private final static int _bas_lbl = 5;
    // arrays filled by $CONTRACT section
    private final static int _ncomp_arr =  6;
    private final static int _nprim_arr =  7;
    private final static int _nptr_arr  =  8;
    private final static int _exp_arr   =  9;
    private final static int _CS_arr    = 10;
    private final static int _CP_arr    = 11;
    private final static int _CD_arr    = 12;
    private final static int _CF_arr    = 13;
    //--------------------------------------------------------------------------
    // basis function remap
    int[] basis_remap;  // index = 0..(nbas-1), value = proper position in a MOLDEN-ordering
    //--------------------------------------------------------------------------
    // section IDs
    final static int section_unknown = -1;
    final static int section_NBO     = 0; //  '$NBO'
    final static int section_GENNBO  = 1; //  '$GENNBO'
    final static int section_COORD   = 2; //  '$COORD '         - geometry
    final static int section_BASIS   = 3; //  '$BASIS'          - center IDs and (L,M)-types
    final static int section_OVERLAP = 4; //  '$OVERLAP'        - AO overlap matrix
    final static int section_DENSITY = 5; //  '$DENSITY'        - density matrix
    //final static int section_FOCK    = 6; //  '$FOCK'        - Fock matrix - TODO!
    final static int section_CONTRACT= 7; //  '$CONTRACT'       - expansion of basis functions over primitive gaussians
    //final static int section_LCAOMO  = 8; //  '$LCAOMO'       - expansion of MO over AO - TODO!
    //final static int section_DIPOLE= 9; //  '$DIPOLE'       - some info on dipole integrals - to be ignored!!!

    //--------------------------------------------------------------------------
    /** Writes a value to a 'current position' (indexed by indices from the global array _MAT_ind) of a given int[]
     * and increases this 'position' 
     */
    private void _add_to_int_arr(int[] arr, int index_I, int value) {
        int i = _MAT_ind[index_I];
        arr[i] = value;
        ++_MAT_ind[index_I];
    }
    //--------------------------------------------------------------------------
    /** Writes a value to a 'current position' (indexed by indices from the global array _MAT_ind) of a given double[]
     * and increases this 'position' 
     */
    private void _add_to_dbl_arr(double[] arr, int index_I, double value) {
        int i = _MAT_ind[index_I];
        arr[i] = value;
        ++_MAT_ind[index_I];
    }
    //--------------------------------------------------------------------------
    /** Writes a value to a 'current position' (indexed by indices from the global array _MAT_ind) of a given matrix[][]
     * and increases this 'position' according to the current value of upper flag
     */     
    private void _add_to_matrix(double[][] matrix, int index_I, int index_J, double value) {
        // get current index values
        int i = _MAT_ind[index_I];
        int j = _MAT_ind[index_J];
        // save the value
        matrix[i][j] = value;
        if (upper) matrix[j][i] = value;    // complementary element to symmetric matrices

        // go to the next element
        if ((upper && (j == i)) || ( j == (nbas-1))) {
            // go to the next line if i==j in the 'upper' mode or
            // if the rightmost element has been reached
            ++_MAT_ind[index_I];
            _MAT_ind[index_J] = 0;
        } else {
            // go to the next element to the right otherwise
            ++_MAT_ind[index_J];
            
        }
    }
    //--------------------------------------------------------------------------
    // do some post-processing of the section data
    private void _Finalize_section(int section) {
        switch (section) {
            case section_GENNBO:
                // allocate memory for arrays (matrices and 1D-arrays)

                // $DENSITY section
                Density = new double[nbas][nbas];
                _MAT_ind[_D_i] = 0; // a 'global' index within a Density array
                _MAT_ind[_D_j] = 0;

                // $OVERLAP section
                Overlap = new double[nbas][nbas];
                _MAT_ind[_S_i] = 0;
                _MAT_ind[_S_j] = 0;

                // 'center =' of $BASIS section
                basis_center = new int[nbas];
                _MAT_ind[_bas_cntr] = 0;

                // 'label  =' of $BASIS section
                basis_label = new int[nbas];
                _MAT_ind[_bas_lbl] = 0;

                // basis function remap (any ordering) => (molden ordering) within each L/center combination
                basis_remap = new int[nbas];
                
                break;
        }
    }
    //--------------------------------------------------------------------------
    PrintStream out = System.out;
    //--------------------------------------------------------------------------
    /**
     * Parses the .47 file
     */
    public void import_from_file(String fname) throws Exception{

        boolean _CONTRACT_allocated = false;
        
        out.println("Reading data from the "+fname+" nbo3-compatible .47 file...");

        FileReader fr = new FileReader(fname);    // reads a file as a stream of chars
        BufferedReader in = new BufferedReader(fr);  // collects chars to strings

        int nDens = 0;
        int nOvrlp = 0;
        
        String s, line;
        int[] basis_intarr = null;   // a pointer to a 'current' array of int
        int basis_intarr_IND = -1;   // array index identifier
        int contract_IND = -1;       // array index identifier for $CONTRACT section arrays

        int section = section_unknown;
        while ((line = in.readLine()) != null) {

            line = line.trim();
            if (line.isEmpty())
                continue;           // skip the next instructions and go to the while condition

            String[] data = line.split("[ =]+");
            int i=0;
            while (i < data.length) {
                s = data[i];
                int increment = 0;
                
                // does this data still belong to the current section?
                if (s.equals("$END")) {
                    _Finalize_section(section);
                    section = section_unknown;
                    increment = 1;
                }

                switch (section) {
                    //----------------------------------------------------------
                    case section_OVERLAP:
                        // read in the (Overlap_i, Overlap_j)-th element of the Overlap matrix
                        _add_to_matrix(Overlap, _S_i, _S_j, Double.parseDouble(s));
                        increment = 1;
                        ++nOvrlp;
                        break;
                    //----------------------------------------------------------
                    case section_DENSITY:
                        // read in the (Density_i, Density_j)-th element of the density matrix
                        _add_to_matrix(Density, _D_i, _D_j, Double.parseDouble(s));
                        increment = 1;
                        ++nDens;                        
                        break;
                    //----------------------------------------------------------
                    case section_BASIS:
                        // set a basis_intarr pointer and _MAT_ind variable number to a proper array attributes
                        if (s.toLowerCase().equals("center")) {
                            basis_intarr = basis_center;
                            basis_intarr_IND = _bas_cntr;
                            increment = 1;
                        } else
                        if (s.toLowerCase().equals("label")) {
                            basis_intarr = basis_label;
                            basis_intarr_IND = _bas_lbl;
                            increment = 1;
                        } else
                        {
                            // s probably contains some int value:
                            if (s.matches("[0-9]+")) {
                                _add_to_int_arr(basis_intarr, basis_intarr_IND, Integer.parseInt(s));
                                increment = 1; // one value parsed successfully
                            }
                        }
                        break;
                    //----------------------------------------------------------
                    case section_GENNBO:
                        // parse $GENNBO section lines - some general info about this file
                        if (s.toLowerCase().equals("bodm")) { bodm = true; increment = 1;}
                        if (s.toLowerCase().equals("upper")) { upper = true; increment = 1;}
                        if (s.toLowerCase().equals("bohr")) { bohr = true; increment = 1; }
                        if (s.toLowerCase().equals("nbas")) {
                            nbas = Integer.parseInt( data[i+1] );
                            increment = 2;
                        }
                        if (s.toLowerCase().equals("natoms")) {
                            natoms = Integer.parseInt( data[i+1] );
                            increment = 2;
                        }
                        
                        if (increment == 0) {
                            WarningManager.warning_printf("WARNING: UNKNOWN KEYWORD: "+s);
                        }
                        break;
                    //----------------------------------------------------------                        
                    case section_CONTRACT:
                        // parse $CONTRACT section
                        if (s.toLowerCase().equals("nshell") || s.toLowerCase().equals("nshells"/*NwChem says so...*/)) {
                            nshell = Integer.parseInt(data[i+1]);   // parse value
                            increment = 2; // skip next parameter
                            // allocate memory for  int[] ncomp, nprim and nptr arrays
                            ncomp = new int[nshell];
                            nprim = new int[nshell];
                            nptr =  new int[nshell];
                            // set current 'positions' to zero
                            _MAT_ind[_ncomp_arr] = 0;
                            _MAT_ind[_nprim_arr] = 0;
                            _MAT_ind[_nptr_arr] = 0;
                        } else
                        if (s.toLowerCase().equals("nexp")) { 
                            nexp = Integer.parseInt(data[i+1]);  // parse value
                            increment = 2; // skip next parameter
                            // allocate memory for double[] c_exp, c_CS, c_CP, c_CD, c_CF;
                            c_exp = new double[nexp];
                            c_CS = new double[nexp];
                            c_CP = new double[nexp];
                            c_CD = new double[nexp];
                            c_CF = new double[nexp];
                            // set current 'positions' to zero
                            _MAT_ind[_exp_arr] = 0;
                            _MAT_ind[_CS_arr] = 0;
                            _MAT_ind[_CP_arr] = 0;
                            _MAT_ind[_CD_arr] = 0;
                            _MAT_ind[_CF_arr] = 0;
                        } else
                        // parsing of array names inside the section:
                        if (s.toLowerCase().equals("ncomp")) { contract_IND = _ncomp_arr; increment = 1;} else
                        if (s.toLowerCase().equals("nprim")) { contract_IND = _nprim_arr; increment = 1;} else
                        if (s.toLowerCase().equals("nptr")) { contract_IND = _nptr_arr; increment = 1;} else
                        if (s.toLowerCase().equals("exp")) { contract_IND = _exp_arr; increment = 1;} else
                        if (s.toLowerCase().equals("cs")) { contract_IND = _CS_arr; increment = 1;} else
                        if (s.toLowerCase().equals("cp")) { contract_IND = _CP_arr; increment = 1;} else
                        if (s.toLowerCase().equals("cd")) { contract_IND = _CD_arr; increment = 1;} else
                        if (s.toLowerCase().equals("cf")) { contract_IND = _CF_arr; increment = 1;} else
                        { // if none of the above, probably s contains some int/double value: ?
                            if (s.matches("[0-9\\.DE\\+\\-]+")) {
                                boolean _is_double = false;
                                int[] arrI = null;  // 'pointers'
                                double[] arrD = null;
                                // some 'individual' parsing for each of the arrays
                                switch (contract_IND) {
                                    case _ncomp_arr: _is_double = false; arrI = ncomp; break;
                                    case _nprim_arr: _is_double = false; arrI = nprim; break;
                                    case _nptr_arr: _is_double = false; arrI = nptr; break;
                                    case _exp_arr: _is_double = true; arrD = c_exp; break;
                                    case _CS_arr: _is_double = true; arrD = c_CS; break;
                                    case _CP_arr: _is_double = true; arrD = c_CP; break;
                                    case _CD_arr: _is_double = true; arrD = c_CD; break;
                                    case _CF_arr: _is_double = true; arrD = c_CF; break;
                                }
                                if (_is_double)
                                    _add_to_dbl_arr(arrD, contract_IND, Double.parseDouble(s));
                                else
                                    _add_to_int_arr(arrI, contract_IND, Integer.parseInt(s));
                                
                                increment = 1; // one value parsed successfully
                            }
                        }

                        break;
                }
                if (increment == 0 ) {
                    // if the section data has not been parsed
                    // section header parser
                    if (s.equals("$NBO")) { section = section_NBO; increment = 1; } else
                    if (s.equals("$GENNBO")) { section = section_GENNBO; increment = 1; } else
                    if (s.equals("$COORD")) { section = section_COORD; increment = 1; } else
                    if (s.equals("$BASIS")) { section = section_BASIS; increment = 1; } else
                    if (s.equals("$OVERLAP")) { section = section_OVERLAP; increment = 1; } else
                    if (s.equals("$DENSITY")) { section = section_DENSITY; increment = 1;} else
                    if (s.equals("$CONTRACT")) { section = section_CONTRACT; increment = 1;} else
                        {
                            //out.println("WARNING: data \""+s+"\" skipped");
                            increment = 1;
                        }
                }
                i += increment;
            }

        }
        out.println("Final array index values: ");
        for (int i=0; i<_MAT_ind.length; ++i) out.printf("%d\t",_MAT_ind[i]);
        out.println();

        out.println(" nbas = "+nbas);
        out.println(" "+nDens+" elements of the density matrix read");
        out.println(" "+nOvrlp+" elements of the overlap matrix read");
        in.close();
        fr.close();
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /*
    // searches proper_ordering_CART[][] array for ID47 and returns its absolute index
    int _proper_cart_ID_from_47ID(int ID47) {
        boolean found = false;
        int result = 0;
        for (int L=0; L<proper_ordering_CART.length; ++L) {
            for (int k=0; k<proper_ordering_CART[L].length; ++k) {
                found = (proper_ordering_CART[L][k] == ID47);
                if (found) break;
                ++result;
            }
            if (found) break;
        }
        if (found)
            return result;
        else  return -1;
    }
     * 
     */
    //--------------------------------------------------------------------------
    /**
     * Creates a remap from a .47 file basis function ordering to a molden-compatible one
     * and stores it in basis_remap array
     * @param molden is used to check for the proper basis function ordering
     */
    public boolean Check_BF_order_and_create_remap(MOLDEN_IO molden) {
        // check center/L values/ordering of basis functions in the .47 file
        boolean ok = true;
        for (int i=0; i<nbas; ++i) {            
            if (molden.Basis[i].Center_ID != basis_center[i]) {
                out.printf("ERROR: Basis function %d center ID mistach: %d in .47 file v.s. %d required by molden file%n",
                        i+1, basis_center[i], molden.Basis[i].Center_ID);
                ok = false;
            }            
            if (molden.Basis[i].L != (basis_label[i] / 100)) {
                out.printf("ERROR: Basis function %d center L mistach: %d in .47 file v.s. %d required by molden file%n",
                        i+1, (basis_label[i] / 100), molden.Basis[i].L);
                ok = false;
            }
        }
        if (!ok)
            return false;
        //----------------------------------------------------------------------
        /*if (!molden.IsSpherical) {
            out.print("SORRY, THIS IS NOT SUPORTED FOR CARTESIAN BASIS FUNCTIONS !!!");
            return false;
        }
         * 
         */
        //----------------------------------------------------------------------
        // if we are here, everything is OK with center/L values!
        // Some change in the function ordering might be neccessary since a quantum
        // chemistry software can reorder m-components in an unpredictable way...
        // (as the ORCA does for example...)
        
        // int[] basis_remap array: index = 0..(nbas-1), value = proper position in a MOLDEN-ordering
        int i=0, L;
        /* .47 basis function labels:
        S:
            001 = 051 = s
        P:
            101 = 151 = px
            102 = 152 = py
            103 = 153 = pz
        D:                                      MOLDEN ordering: D 0, D+1, D-1, D+2, D-2
            251 = dxy =     D(-2) = XY
            252 = dxz =     D(+1) = XZ
            253 = dyz =     D(-1) = YZ
            254 = dx2-y2 =  D(+2) = (X^2 - Y^2)/2
            255 = dz2 = d3z2-r2 --------V.S.----------- D( 0) = (Z^2 - 0.5*X^2 - 0.5*Y^2)/SQRT(3)
         NOTE:
           	Pure f "cubic" set  IS NOT SUPPORTED yet
	F: ("standard" Pure f set )             MOLDEN ordering: F 0, F+1, F-1, F+2, F-2, F+3, F-3
            351 = f(0): z(5z2-3r2)= F( 0) = z*(z^2-1.5*y^2-1.5*x^2)/SQRT(15)
            352 = f(c1): x(5z2-r2)= F(+1) = x*(2*z^2-0.5*y^2-0.5*x^2)/SQRT(10)
            353 = f(s1): y(5z2-r2)= F(-1) = y*(2*z^2-0.5*y^2-0.5*x^2)/SQRT(10)
            354 = f(c2): z(x2-y2) = F(+2) = z*(0.5*x^2-0.5*y^2)
            355 = f(s2): xyz =      F(-2) = x*y*z
            356 = f(c3): x(x2-3y2)= F(+3) = x*(0.5*x^2-1.5*y^2)/SQRT(6)
            357 = f(s3): y(3x2-y2)= F(-3) = y*(1.5*x^2-0.5*y^2)/SQRT(6)

        G:
         g (451-459) 0, c1, s1, c2, s2, c3, s3, c4, s4
        G( 0) = (z^4-3*y^2*z^2-3*x^2*z^2+3/8*y^4+3/4*x^2*y^2+3/8*x^4)/SQRT(35) ~ C40
        G(+1) = x*z*(10*z^2-7.5*y^2-7.5*x^2)/SQRT(350) ~ C41
        G(-1) = y*z*(10*z^2-7.5*y^2-7.5*x^2)/SQRT(350) ~ S41
        G(+2) = (1.5*x^2*z^2-1.5*y^2*z^2 + 0.25*y^4-0.25*x^4)/SQRT(7) ~ C42
        G(-2) = x*y*(3*z^2-0.5*y^2-0.5*x^2)/SQRT(7) ~ S42
        G(+3) = x*z*(0.5*x^2-1.5*y^2)/SQRT(2) ~ C43
        G(-3) = 0.5*y*z*(3*x^2-y^2)/SQRT(2) ~ S43
        G(+4) = 1/8*x^4-3/4*x^2*y^2+1/8*y^4 ~ C44
        G(-4) = 0.5*x*y*(x^2-y^2) ~ S44
         It seems, that no true reordering is neccessary for the functions except for D-ones
         HOWEVER, some reordering might still be neccessary due to uncertainty in the 'internal' ordering
           used in the quantum chemistry package while writing the .47 file

         */
        final int[][] proper_ordering_SPH = new int[][]{
            { 51 } /*S*/, // note 051 = 41 decimal...
            {151,152,153}/*P*/,
            {255, 252, 253, 254, 251}/*D*/,
            {351, 352, 353, 354, 355, 356, 357 }/*F*/ ,
            {451, 452, 453, 454, 455, 456, 457, 458, 459} /*G*/,
            {551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561}/*H*/};
        /*
        CARTESIAN:
            s (001) = s
         
            p (101-103) = x, y, z

                          201 202 203 204 205 206
            d (201-206) = xx, xy, xz, yy, yz, zz
                 MOLDEN = xx, yy, zz, xy, xz, yz
         => reordered =>  201 204 206 202 203 205

                          301  302  303  304  305  306  307  308  309  310
            f (301-310) = xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
                 MOLDEN = xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
         => reordered =>  301  307  310  304  302  303  306  309  308  305
         => reordered =>  xxx,  yyy,  zzz  xyy,  xxy,  xxz,  xzz,  yzz,  yyz,  xyz,

                          401,  402,  403,  404,  405,  406,  407,  408,  409   410   411   412   413   414   415
            g (401-415) = xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz
                 MOLDEN = xxxx  yyyy  zzzz  xxxy  xxxz  yyyx  yyyz  zzzx  zzzy  xxyy  xxzz  yyzz  xxyz  yyxz  zzxy
         => reordered =>  xxxx  yyyy  zzzz  xxxy  xxxz  xyyy  yyyz  xzzz  yzzz  xxyy  xxzz  yyzz  xxyz  xyyz  xyzz
         => reordered =>  401,  411,  415,  402,  403,  407,  412,  410,  414,  404,  406,  413,  405,  408,  409

         */
        // A 1-purpose array: 1) re-ordering of MOLDEN basis functions in accordance to their sequence specified by D/S matrices .47 file
        // /////////////and 2) assigning proper absolute values of m for cartesian functions with proper .47-ID
        final int[][] proper_ordering_CART = new int[][]{
            { 51 }/**/, //{ 1 }, // note: in Java 051 = 41 decimal...
            {151, 152, 153}/*P: a 'hack' - we use 'spherical' P functions always */, //{101, 102, 103}/*P*/,
            {201, 204, 206, 202, 203, 205 }/*D*/,
            {301, 307, 310, 304, 302, 303, 306, 309, 308, 305 }/*F*/ ,
            {401, 411, 415, 402, 403, 407, 412, 410, 414, 404, 406, 413, 405, 408, 409 } /*G*/};

        // simply the number of cartesian components for each L
        final int[] cart_compon_count = new int[]{1/*S*/, 3/*P*/, 6/*D*/, 10/*F*/, 15/*G*/}; // equals the lengths of proper_ordering_CART[]

        // final int[][] proper_ordering_CART is declared above
        // final int[] cart_compon_count is declared above

        // loop over all basis functions once again
        while (i<nbas) {
            L = molden.Basis[i].L;

            int component_count;
            if (molden.IsSpherical)
                component_count = (2*L+1);
            else
                component_count = cart_compon_count[L];
            
            for (int M=0; M<component_count; ++M) {
                // find suitable offset for this function
                int label_to_find = basis_label[i+M];
                // a small 'hack': make 051/001(s), 101/151(px), 102/152(py), 103/153(pz) labels indistinguishable:
                if (label_to_find < 154) label_to_find = (label_to_find / 100)*100 + 50 + (label_to_find % 10);

                // find label_to_find in a proper L-subarray of proper_ordering_SPH
                boolean found = false;
                int offset = 0;
                if (molden.IsSpherical)
                    // spherical case
                    while ((!found) && (offset < proper_ordering_SPH[L].length)) {
                        found = (proper_ordering_SPH[L][offset] == label_to_find);
                        if (!found) ++offset;
                    }
                else
                    // cartesian case
                    while ((!found) && (offset < proper_ordering_CART[L].length)) {
                        found = (proper_ordering_CART[L][offset] == label_to_find);
                        if (!found) ++offset;
                    }
                    
                if (!found) {
                    out.printf("ERROR: function label %d is unknown for L = %d %n", basis_label[i+M], L);
                    return false;
                }
                // ok, we've found a position of the M-th component of (2L+1)-component function family
                // (the family sterted with the index i)
                // Now use this offset to calculate a proper position of the M-th component in the final MOLDEN file:
                basis_remap[i+M] = i+offset;
            }
            i += component_count/*(2*L+1)*/; // go to the next (2L+1)-component function family
        }

        // some debug
        /*
        for ( i=0; i<nbas; ++i)
            out.printf("%.10E%n",Overlap[nbas-1][basis_remap[i]]);
         */

        return true;
    }
    //--------------------------------------------------------------------------
    /**
     * Replaces the GEOMETRY in targetMolden witg the one created from the information from $GEOM section
     */
    public boolean Insert_Geom_Into(MOLDEN_IO targetMolden) {
        return true;
    }
    //--------------------------------------------------------------------------
    /**
     * Replaces the basis in targetMolden witg the one created from the information from $CONTRACT section
     */
    public boolean Build_Basis_Set(MOLDEN_IO targetMolden) {
        // Process 'shells' (radial parts)
        // nprim[]: number of primitive gaussians per basin function
        // nptr[i]: 1-based index of the first exponent/contraction coefficient for the i-th radial part
        // ncomp[]: number of function components: 1=S, 3=P, 4=SP, 5=5 of D/6=6 of D, ...
        
        // Get the number of radial parts; this is not very trivial due to the use of SP shells
        int nRadialParts = 0;
        for (int sh=0; sh<nshell; ++sh)
            if (ncomp[sh] != 4) 
                ++nRadialParts; /*'normal' shells use one radial part*/
            else nRadialParts+=2; /*SP shells needs to use 2 radial parts*/
        
        targetMolden._RadialParts_Alloc(nRadialParts);    // allocate memory (and clear previous values)
        // bugfix 03.Sep.2016: create objects in the targetMolden.RadialParts[] array 
        // allocated eralier by targetMolden._RadialParts_Alloc(nRadialParts)
        for (int i=0; i<targetMolden.RadialParts.length; ++i)
            targetMolden.RadialParts[i] = new RadialPartOfBasisFunction();            
        
        int[] nRP_Components = new int[nRadialParts]; // true number of components where the i-th radial part is really used
        int L_max = 0;

        boolean Cartesian_D_found = false;
        boolean Cartesian_F_found = false;
        boolean Cartesian_G_found = false;

        int rp = 0; // index of the current radial part
        for (int sh=0; sh<nshell; ++sh) {            
            // loop over all 'shells' (=radial parts)
            targetMolden.RadialParts[rp].Coefs = new double[ nprim[sh] ];
            targetMolden.RadialParts[rp].Exponents = new double[ nprim[sh] ];
            // determine L
            int L = -2; // unknown
            switch(ncomp[sh]) {
                case 1: L = 0; break; // S
                case 3: L = 1; break; // P
                case 4: L = -1; break; // SP
                case 5: case 6: L = 2; break; // spherical D / cartesian D
                case 7: case 10: L = 3; break; // spherical F / cartesian F
                // non-documented cases:
                case 9: case 15: L = 4; break; // spherical G / cartesian G
            }            
            if (L == -2) {
                out.printf("Unrecognized ncomp value: ncomp[%d] = %d %n", sh+1, ncomp[sh]);
                return false;
            }
            // special treatment of SP shell case: create an S-function and then set L to 1
            if (L == -1) {
                for (int g=0; g<nprim[sh]; ++g) { // loop over primitive gaussians
                    targetMolden.RadialParts[rp].Exponents[g] = c_exp[ nptr[sh]-1 + g ];
                    // .47 files use contraction coefficients with NORMALIZED primitives => no conversion factor is needed
                    // pick the contraction coefficient from the proper array
                    targetMolden.RadialParts[rp].Coefs[g] = c_CS[ nptr[sh]-1 + g ];
                }
                targetMolden.RadialParts[rp].LUsedWith = 0;
                nRP_Components[rp] = 1;
                ++rp; // go to the next function
                L = 1;
            }
            // non-SP shells (or P-component of SP-ones)
            for (int g=0; g<nprim[sh]; ++g) { // loop over primitive gaussians
                targetMolden.RadialParts[rp].Exponents[g] = c_exp[ nptr[sh]-1 + g ];
                // .47 files use contraction coefficients with NORMALIZED primitives => no conversion factor is needed
                // pick the contraction coefficient from the proper array
                switch(L) {
                    case 0: targetMolden.RadialParts[rp].Coefs[g] = c_CS[ nptr[sh]-1 + g ]; break;
                    case 1: targetMolden.RadialParts[rp].Coefs[g] = c_CP[ nptr[sh]-1 + g ] ;break;
                    case 2: targetMolden.RadialParts[rp].Coefs[g] = c_CD[ nptr[sh]-1 + g ] ;break;
                    case 3: targetMolden.RadialParts[rp].Coefs[g] = c_CF[ nptr[sh]-1 + g ]; break;
                }                
            }
            targetMolden.RadialParts[rp].LUsedWith = L;                
            
            // write a correct number of basis functions where this radial part should be used
            if (ncomp[sh] != 4)
                nRP_Components[rp] = ncomp[sh]; // for non-SP shells ncomp[sh] contains an adequate value of the number components where this radial part is to be used
            else nRP_Components[rp] = 3; // since we treat it as a P-function now;
            // was it a cartesian function? (inessential for S/P/SP ones)
            if ((L == 2) && (nRP_Components[rp] != 5)) Cartesian_D_found = true;
            if ((L == 3) && (nRP_Components[rp] != 7)) Cartesian_F_found = true;
            if ((L == 4) && (nRP_Components[rp] != 9)) Cartesian_G_found = true;

            // get a statistics
            if (L>L_max) L_max = L;

            // go to the next radial part index
            ++rp;
        }
        out.printf("%d radial parts have been assembled, L_MAX = %d%n", rp, L_max);

        // set proper cartesian/spherical flags
        targetMolden.IsSpherical = (!Cartesian_D_found) & (!Cartesian_F_found) & (!Cartesian_G_found);
        if (!targetMolden.IsSpherical) {
            if ((L_max>=2) && (!Cartesian_D_found)) {
                out.println("Mixing of spherical D / cartesian others basis is not supported!");
                return false;
            }
            if ((L_max>=3) && (!Cartesian_F_found)) {
                out.println("Mixing of spherical F / cartesian others basis is not supported!");
                return false;
            }
            if ((L_max>=4) && (!Cartesian_G_found)) {
                out.println("Mixing of spherical G / cartesian others basis is not supported!");
                return false;
            }
        }
        if (targetMolden.IsSpherical)
            out.println("Building SPHERICAL basis set...");
        else
            out.println("Building CARTESIAN basis set...");

        // values of m for cartesian functions
        final int[][] cart_compon_m = new int[][]{
            {0}, /*S*/
            {1, 2, 3}, /*P*/
            {4, 5, 6, 7, 8, 9}, /*D*/
            {10, 11, 12, 13, 14, 15, 16, 17, 18, 19}, /*F*/
            {20,21,22,23,24,25,26,27,28,29,30,31,32,33,34} /*G*/};

        
        // Now use these radial parts to build the basis functions
        // In fact, the only really essintial step here is setting proper center_IDs; however, we build true functions to be consistent with other possible uses of this method
        targetMolden.Basis = new BasisFunction[nbas];
        int bf = 0; // basis function index
        for (/*int*/rp=0; rp<nRadialParts; ++rp) { // loop over true radial parts
            int m = 0; // for a spherical basis case this should be so
            if (!targetMolden.IsSpherical) 
                m = cart_compon_m[targetMolden.RadialParts[rp].LUsedWith][0];
            
            // create new basis function (or just its 'first component', if L>0)
            targetMolden.Basis[bf] = new BasisFunction(targetMolden.RadialParts[rp].LUsedWith,
                    m, null, targetMolden.RadialParts[rp].Exponents.length );
            
            // set proper center ID for the radial part (!!!) on the basis of the information from $BASIS section
            targetMolden.RadialParts[rp].CenterID = basis_center[bf];
            targetMolden.Basis[bf].Center_ID = basis_center[bf];
            targetMolden.Basis[bf].RadialPart_ID = rp; // bugfix 03.Sep.2016
            // set atom coords.:
            targetMolden.Basis[bf].R0 = targetMolden.Centers[targetMolden.Basis[bf].Center_ID-1].R0;
            // copy coefs and exponents from the radial part
            targetMolden.Basis[bf].exponents = targetMolden.RadialParts[rp].Exponents.clone();
            targetMolden.Basis[bf].coefs = targetMolden.RadialParts[rp].Coefs.clone();
            
            ++bf; // one function has been created

            // create all the other components
            for (int c=0; c<(nRP_Components[rp]-1/*one function has already been created*/); ++c) {
                targetMolden.Basis[bf] = new BasisFunction(targetMolden.Basis[bf-c-1]); // copy all from the previous basis function
                // correct the value of m
                if (targetMolden.IsSpherical)
                    // proper MOLDEN numbering is:  G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
                    //                                   c=0  c=1  c=2  c=3  c=4  c=5 ...
                    if ((c%2)==0) targetMolden.Basis[bf].m = (c/2)+1; else targetMolden.Basis[bf].m = -((c/2)+1);
                else
                    targetMolden.Basis[bf].m = cart_compon_m[targetMolden.RadialParts[rp].LUsedWith][c+1];

                // one more basis function has been created!
                ++bf;
            }
        }
        // Note that MO coefs are NOT affected now, but they might need to be affected (re-ordered) in future
        // To leave this indicated we create basis functions in a proper MOLDEN, rather than original .47, sequence!

        out.printf("%d basis functions have been built (%d were expected)%n", bf, nbas);

//        try{
        /*            
        }catch (Exception e) {}
         * 
         */
        return true;        
    }
    //--------------------------------------------------------------------------
    /**
     * Diagonalizes the density matrix so as to create Lowdin natural orbitals
     */
    public void Procude_NOs(MOLDEN_IO targetMolden) {
        out.println("Creating natural orbitals from the .47 file density/overlap matrices data...");
        // a wrapper matrices for the two arrays
        Matrix D = new Matrix(Density);
        Matrix S = new Matrix(Overlap);
        //S.print(13, 7);
        out.printf(" tr(S) = %.10f%n", S.trace());
        out.printf(" tr(D) = %.10f%n", D.trace());

        if (!bodm) {
            out.println("ERROR: non-BODM density matrices are not supported yet!");
            return;
        }

        // The basic eigenproblem is
        //  D.S . y = lambda . y
        // here S = (S^(1/2))^T . (S^(1/2)), so
        //  (S^(1/2)).D.(S^(1/2))^T . (S^(1/2)).y = lambda . (S^(1/2)).y
        // Introduce  x = (S^(1/2)).y
        // Then
        //  (S^(1/2)).D.(S^(1/2))^T . x = lambda . x
        // Introcude DD = (S^(1/2)).D.(S^(1/2))^T
        // Then
        //  DD . x = lambda . x
        // is the eigenproblem to solve containing the symmetrical matrix DD 
        Matrix S05 = EigenEngine.Matrix_SQRT(S, false);
        //Matrix DD = LDLtransform.TransformSymmetricMatrixToNewBasis_LDLT(D, S05);        
        Matrix DD = EigenEngine.TransformSymmetricMatrixToNewBasis(D, S05);
        EigenvalueDecomposition eig =  DD.eig();
        Matrix Occupancies = eig.getD();
        //Occupancies.print(13, 7);
        out.println("tr(EIG) = "+ eig.getD().trace() );
        // transform eigenvectors
        Matrix Sm05 = EigenEngine.Matrix_SQRT(S, true);
        Matrix U = Sm05.times(eig.getV());
        //U.print(13, 7);

        // convert this data to molden file
        targetMolden.MOs = new MO[nbas]; // here nbas = dimension of the matrices => the number of eigenvalues
        for (int i=0; i<nbas; ++i) {
            targetMolden.MOs[i] = new MO();
            targetMolden.MOs[i].Energy = 0.0;
            targetMolden.MOs[i].Occupancy = Occupancies.get(i, i);
            targetMolden.MOs[i].BS_Coefs = new double[nbas];
            for (int cf=0; cf<nbas; ++cf)
                targetMolden.MOs[i].BS_Coefs[ basis_remap[cf] ] = U.get(cf, i);
        }
        //
        //U.transpose().times(S).times(U).print(13,7);
    }
    //--------------------------------------------------------------------------

}
