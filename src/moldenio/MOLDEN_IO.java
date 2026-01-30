package moldenio;

/**
 * A file implementing a MOLDEN_IO class for reading and writing
 * text files in a MOLDEN-compatible format (as described at
 * http://www.cmbi.ru.nl/molden/molden_format.html )
 *
 * Version: 12.Jan.2014
 * Created: 26.Oct.2013
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
 *  Some features:
 *   * floating-point values with 'D' instead of 'E' are not fully supported
 *
 */

import JGints.*;
import java.io.*;
import ProgramOptions.WarningManager;

/**
 * A class representing for importing/exporting orbitals from/to text files
 * @author timn
 */
public class MOLDEN_IO {
    public PrintStream out = System.out;
    
    private final boolean debugPrint = false;
    public boolean IsSpherical = false; // default for MOLDEN
    public BasisFunction[] Basis; // basis functions with proper L and m values
    private BasisFunction FirstBasisFunction = null;
    private int NBasisFunctions = 0;
    public int HighestL = 0; // highest value of angular momentum in the basis; is being set by Load_From_MOLDEN()
    final public String spdf = "spdfg";
    //
    // Due to many number of "clones" of basis functions (which differ only with the value of m, but share the same L),
    // and due to computational convenience (in particular, for correct treatment of SP shells) in NAO procedure,
    // we store radial parts of basis functions in a separate array.
    // This array is being constructed in the _BS_List_To_Array() method
    /*
    public double[][] RadialPartsOfBasisFunctions_Exponents = null;
    public double[][] RadialPartsOfBasisFunctions_Coefs = null;
    public int[] RadialPartsOfBasisFunctions_CenterIDs = null; // 1-based identifier of the nuclei !!!
    public int[] RadialPartsOfBasisFunctions_LUsedWith = null;
    public int[] RadialPartsOfBasisFunctions_Addit_r_power = null;
     * 
     */
    public RadialPartOfBasisFunction[] RadialParts = null;
    //--------------------------------------------------------------------------
    /** Allocates memory for radial parts arrays (1-st level only!) */
    // Created: 18.Jan.2014
    public void _RadialParts_Alloc(int nRadParts) {
        RadialParts = new RadialPartOfBasisFunction[nRadParts];
        /*
        RadialPartsOfBasisFunctions_Exponents = new double[nRadParts][];
        RadialPartsOfBasisFunctions_Coefs = new double[nRadParts][];
        RadialPartsOfBasisFunctions_CenterIDs = new int[nRadParts];
        RadialPartsOfBasisFunctions_LUsedWith = new int[nRadParts];
        RadialPartsOfBasisFunctions_Addit_r_power = new int[nRadParts];
         * 
         */
    }
    //--------------------------------------------------------------------------

    // Note: these arrays are the source of basis set information used by Save_As_MOLDEN
    public boolean Allow_additional_r_power = false; // whether the "s2 1 1.00" is acceptable
    
    public String Title=null;
    //
    public MO[] MOs;
    private MO FirstMO = null;

    public boolean Coords_in_AU = false;       // default is Angstroms
    //
    public AtomicCenter[] Centers;
    private AtomicCenter FirstAtom = null;
    //
    private boolean use_5D = false; // default is cartesian
    private boolean use_7F = false; // default is cartesian
    private boolean use_9G = false; // default is cartesian
    //--------------------------------------------------------------------------
    private final int section_unknown = -1;
    private final int section_MOLDEN  = 0;
    private final int section_Title   = 1;
    private final int section_Atoms   = 2;
    private final int section_GTO     = 3;
    private final int section_MO      = 4;
    private final int section_5D      = 5;
    private final int section_5D7F    = 6;
    private final int section_9G      = 7;
    /**
     *  translates "[GTO]", ... to one of pre-defined IDs
     */
    private int _Parse_Section_ID(String s) {
        s = s.toUpperCase();
        if (s.contains("[MOLDEN FORMAT]")) return  section_MOLDEN; // not really used
        if (s.contains("[TITLE]")) return section_Title;
        if (s.contains("[ATOMS]")) return section_Atoms;
        if (s.contains("[GTO]")) return section_GTO;
        if (s.contains("[MO]")) return section_MO;
        if (s.contains("[5D]")) return section_5D;
        if (s.contains("[5D7F]")) return section_5D7F;
        if (s.contains("[9G]")) return section_9G;
        return section_unknown;
    }
    //--------------------------------------------------------------------------
    final public int[] CartesianHarmonicsCount = new int[]{ 1/*S*/, 3/*P*/, 6/*D*/, 10/*F*/, 15/*G*/}; 
        // where index = L, value = number of harmonics
    /* The following values of m identify a cartesian primitive spherical harmonic type:
       1S (L=0): 1(m=0)
       3P (L=1): x(m=1), y(2), z(3)
       6D (L=2): xx(4), yy(5), zz(6), xy(7), xz(8), yz(9)
       10F(L=3): xxx(10), yyy(11), zzz(12), xyy(13), xxy(14), xxz(15), xzz(16), yzz(17), yyz(18), xyz(19)
       15G(L=4): xxxx(20), yyyy(21), zzzz(22), xxxy(23), xxxz(24), yyyx(25), yyyz(26), zzzx(27),
                 zzzy(28), xxyy(29), xxzz(30), yyzz(31), xxyz(32), yyxz(33), zzxy(m=34)           */
    final public int[] First_Cart_m_for_L = new int[]{ 0/*S*/, 1/*P*/, 4/*D*/, 10/*F*/, 20/*G*/}; 
        // where index = L, value = m-index of the first cartesian harmonic with this L
    //--------------------------------------------------------------------------
    // Convert basis function list to an array and make clones m=(-L+1)...+L
    // Version: 14.Jan.2014
    // Created:    Oct.2013
    private void _BS_List_To_Array() throws Exception {
        // get their total number
        BasisFunction bs_tmp, LastBasisFunction;
        //
        int NRadialParts = 0;
        //
        bs_tmp = FirstBasisFunction;
        while (bs_tmp != null) {
            if (bs_tmp.L != -1) {
                if (IsSpherical)
                    NBasisFunctions += (2*bs_tmp.L + 1);
                else
                    NBasisFunctions += CartesianHarmonicsCount[bs_tmp.L];
            } else
                throw new Exception("SP shells are not supported!"); // TODO
            //
            bs_tmp.RadialPart_ID = NRadialParts; // set proper radial part ID number
            bs_tmp = bs_tmp._next;
            NRadialParts++; // number of different basis functions (number of contracted linear combinations in the basis set specification)
        }
        out.printf("Total number of basis functions: %d\n",NBasisFunctions);
        out.printf("Total number of different contractions (radial parts): %d\n",NRadialParts);
        out.printf("Highest angular momentum: %s (Lmax = %d)\n",new String[]{"s","p","d","f","g"}[HighestL], HighestL);
        // allocate arrays
        Basis = new BasisFunction[NBasisFunctions];
        _RadialParts_Alloc(NRadialParts);
        // allocate sub-array for radial parts according to their acual length
        bs_tmp = FirstBasisFunction;
        int i=0;
        while (bs_tmp != null) {
            RadialParts[i] = new RadialPartOfBasisFunction();
            RadialParts[i].Exponents = bs_tmp.exponents.clone();
            RadialParts[i].Coefs = bs_tmp.coefs.clone();
            RadialParts[i].CenterID = bs_tmp.Center_ID;
            RadialParts[i].LUsedWith = bs_tmp.L;
            RadialParts[i].Addit_r_power = bs_tmp.additional_r_power;
            i++;
            bs_tmp = bs_tmp._next;
        }        
        // Move to array and create clones
        /* The following order of D, F en G functions is expected:
           5D: D 0, D+1, D-1, D+2, D-2  /   6D: xx, yy, zz, xy, xz, yz
           7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3    /   10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
           9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4  /   15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy
         For the sake of similarity we assume that: 3P: x=P(0), y=P(+1), z=P(-1)
         */
        bs_tmp = FirstBasisFunction;
        i = 0;
        while (bs_tmp != null) {
            Basis[i] = bs_tmp; // m = 0;
            i++;
            // make clones (i.e., functions with the same L but with different Ms);
            if (bs_tmp.L != (-1)) {
                if (IsSpherical)
                    // Spherical harmonics cloning: create new functions with m=-L...(but skip_m=0)...+L
                    for (int ii=1; ii<=bs_tmp.L; ii++) {
                        Basis[i] = new BasisFunction(bs_tmp);
                        Basis[i].m = ii;
                        /*
                        Basis[i] = new BasisFunction(bs_tmp.L, ii, bs_tmp.R0, bs_tmp.coefs.length);
                        Basis[i].Center_ID = bs_tmp.Center_ID; // clones share the same center_ID
                        Basis[i].RadialPart_ID = bs_tmp.RadialPart_ID; // clones share the same RadialPart_ID
                        // clone exponents and coefficients
                        System.arraycopy(bs_tmp.exponents, 0, Basis[i].exponents, 0, bs_tmp.exponents.length);
                        System.arraycopy(bs_tmp.coefs, 0, Basis[i].coefs, 0, bs_tmp.coefs.length);
                         *
                         */
                        i++;

                        Basis[i] = new BasisFunction(bs_tmp);
                        Basis[i].m = -ii;
                        /*
                        Basis[i] = new BasisFunction(bs_tmp.L, -ii, bs_tmp.R0, bs_tmp.coefs.length);
                        Basis[i].Center_ID = bs_tmp.Center_ID;
                        Basis[i].RadialPart_ID = bs_tmp.RadialPart_ID;
                        // clone exponents and coefficients
                        System.arraycopy(bs_tmp.exponents, 0, Basis[i].exponents, 0, bs_tmp.exponents.length);
                        System.arraycopy(bs_tmp.coefs, 0, Basis[i].coefs, 0, bs_tmp.coefs.length);
                         *
                         */
                        i++;
                    }
                if (!IsSpherical) {
                    // one harmonic 0 has already been created!
                    /*    Basis[i] = new BasisFunction(bs_tmp); */
                    // but by default its m has been set to 0. set m to proper cartesian harmonic index:
                    bs_tmp.m = First_Cart_m_for_L[bs_tmp.L];
                    // ++i has already been done above !!!
                    
                    // Cartesian harmonics cloning
                    for (int ii=(First_Cart_m_for_L[bs_tmp.L]+1); ii<(First_Cart_m_for_L[bs_tmp.L]+CartesianHarmonicsCount[bs_tmp.L]); ++ii) {
                        Basis[i] = new BasisFunction(bs_tmp);
                        Basis[i].m = ii;
                        ++i;
                    }
                }
            } else
                throw new Exception("sp shells are not supported!"); // TODO
            // go to the next function
            LastBasisFunction = bs_tmp._next;
            bs_tmp._next = null; // to prevent subsequent missuse of this pointer
            bs_tmp = LastBasisFunction;
        }
    }

    //--------------------------------------------------------------------------
    // some output flags for Load_From_MOLDEN() method
    public boolean print_MO_info = false;
    //--------------------------------------------------------------------------
    /**
     * Loads basis functions and MO coefs from MOLDEN file
     * follows a specification at: http://www.cmbi.ru.nl/molden/molden_format.html
     * Supported sections are: [Molden Format], [Title], [ATOMS],[GTO] and [MO]
     * Non-supported sections are:  {[GEOCONV], [GEOMETRIES]}, {[FREQ],[FR-COORD],[FR-NORM-COORD], [INT]}
     *
     * Rev.: 08.Feb.2014 / 16.Jan.2014
     */
    public boolean Load_From_MOLDEN(String filename) throws Exception {
        FileReader fr = new FileReader(filename);    // reads a file as a stream of chars
        BufferedReader in = new BufferedReader(fr, 32768);  // collects chars to strings
        String s;
        String[] arr;
        //
        Title = new String();
        //
        AtomicCenter LastAtom = null;
        BasisFunction LastBasisFunction = null;
        BasisFunction bs_tmp = null;
        int section_NOW = section_unknown;
        boolean NewSectionDetected;
        int Current_Center_Number = 0;
        MO LastMO = null;
        //
        //double[] expontnes_tmp = new double[800];
        //double[] coefs_tmp = new double[800];
        int CoefNum = 0;
        int L=0; // -1 for 'sp'
        // read all lines from molden file
        if ( ! in.readLine().toUpperCase().contains("[MOLDEN FORMAT]")) return false; // check for molden signature
        //
        while ((s = in.readLine()) != null) {            
            NewSectionDetected = false;
            // parse a string according to current section type
            switch (section_NOW) {
                /*
                case section_unknown: case section_MOLDEN:
                    // wait until next section begins and parse its type
                    NewSectionDetected = s.contains("[");
                    break;
                 *
                 */
                //----------------------------------------------
                case section_Title: //  "[Title]" section
                    /* Section format: plain text until next section
                     */
                    NewSectionDetected = s.contains("[");
                    if (!NewSectionDetected) {
                        Title += s;
                    }
                    break;
                //----------------------------------------------
                case section_Atoms: // "[Atoms]" section
                    /* Section format:
                      [Atoms] (Angs|AU)
                      element_name number atomic_number x y z
                      ...
                     */
                    NewSectionDetected = s.contains("[");
                    if (!NewSectionDetected) {
                        // parse atomic coords data
                        arr = s.trim().split(" +"); // one or more spaces as a separator
                        // add new atom
                        if (FirstAtom == null)
                            LastAtom = FirstAtom = new AtomicCenter();
                        else
                            LastAtom = (LastAtom._next = new AtomicCenter());
                        // parse atomic properties
                        LastAtom.Name = arr[0]; // it's ok - strings are Immutable Objects
                        LastAtom.ID = Integer.parseInt(arr[1]);
                        LastAtom.Z = Double.parseDouble(arr[2]);
                        for (int mu=0; mu<3; mu++) LastAtom.R0[mu] = Double.parseDouble(arr[3+mu]);
                    }
                    break;
                //----------------------------------------------
                case section_GTO:
                    /* Section Format:
                      [GTO]
                      atom_sequence_number1 0
                      shell_label number_of_primitives 1.00
                      exponent_primitive_1 contraction_coefficient_1 (contraction_coefficient_1)  - 's','p','d','f','sp','g'
                      ...
                      empty line

                      atom_sequence__number2 0
                      shell_label number_of_primitives 1.00
                      exponent_primitive_1 contraction_coefficient_1 (contraction_coefficient_1)
                      ...
                      empty line
                     * NOTE: molden4.4 The 0 on the shell_number line and the 1.00 on the shell_label line are no longer functional and can be left out.
                     */
                    NewSectionDetected = s.contains("[");
                    if ((!NewSectionDetected) && (!s.trim().equals(""))) {
                        // read until an empty line
                        arr = s.trim().split(" +"); //  trim() Returns a copy of the string, with leading and trailing whitespace omitted.
                        // first line of a block should contain "  1 0" - center ID and a zero
                        Current_Center_Number = Integer.parseInt(arr[0]); // 1-based!!!
                        // read all lines before
                        //ready = false;
                        while (((s = in.readLine()) != null) & (!s.trim().equals(""))) { // The & and | operators, when used as logical operators, always evaluate both sides.
                            // if empty line, next non-empty one will be read on the next iteration of parent 'while'
                            //nTempCoefs = 0;
                            arr = s.trim().split(" +");
                            // first new line should be of type
                            // "s   8 1.0 " - like line  => a new basis function starts here
                                // Check for r-extended notation
                                int additional_r_power = 0;
                                if (arr[0].toLowerCase().matches("[a-z]+[0-9]+")) {
                                    if (!Allow_additional_r_power)
                                        throw new Exception("Extended molden command found ("+arr[0]+") but Allow_additional_r_power is false.");
                                    String[] prms = arr[0].toLowerCase().split("[a-z]+");
                                    additional_r_power = Integer.parseInt( prms[1] ); // prms[0] equals ""
                                    arr[0] = arr[0].split("[0-9]+")[0]; // cut the digits after s/p/d/f/...
                                }
                                
                                // determine value of L
                                L = -2; // default: unknown
                                if (arr[0].toLowerCase().equals("s")) L = 0; else
                                if (arr[0].toLowerCase().equals("p")) L = 1; else
                                if (arr[0].toLowerCase().equals("d")) L = 2; else
                                if (arr[0].toLowerCase().equals("f")) L = 3; else
                                if (arr[0].toLowerCase().equals("g")) L = 4; else
                                if (arr[0].toLowerCase().equals("sp")) L = -1;
                                //
                                
                                CoefNum = Integer.parseInt(arr[1]);
                                //
                                /* The following order of D, F en G functions is expected:
                                   5D: D 0, D+1, D-1, D+2, D-2  /   6D: xx, yy, zz, xy, xz, yz
                                   7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3    /   10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
                                   9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4  /   15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy
                                For the sake of similarity we assume that: 3P: x=P(0), y=P(+1), z=P(-1)
                                 */
                                bs_tmp = new BasisFunction(L, 0, Centers[Current_Center_Number-1].R0, CoefNum);
                                bs_tmp.Center_ID = Current_Center_Number;
                                bs_tmp.additional_r_power = additional_r_power;
                                if (LastBasisFunction == null)
                                    LastBasisFunction = (FirstBasisFunction = bs_tmp);
                                else
                                    LastBasisFunction = (LastBasisFunction._next = bs_tmp);
                                //
                                if (L == -1) {
                                    bs_tmp.L = 0; // the just created function is s-type
                                    bs_tmp.m = 0; // by default; to be corrected during the cloning process
                                    // create one more function, p-type
                                    
                                    // !!!: bs_tmp, pointing to an s-function, should be saved !!!
                                    LastBasisFunction = (LastBasisFunction._next = new BasisFunction(bs_tmp));
                                    LastBasisFunction.L = 1; // correct L from 0 to 1 for the last, p-type, function
                                    // !!!: L==(-1) should be saved for the subsequent code!
                                    if (HighestL < 1) HighestL = 1; // since the 'else'-part will not be executed

                                    //System.out.println(" sp shells are not supported yet!"); // TODO
                                    //return false;
                                    /*
                                    // create one more function
                                    bs_tmp.L = 0; // the created one is of s-type, and its reference is still to be stored in bs_tmp
                                    bs_tmp.m = 0;
                                    // create new function of p-type and shift the "LastBasisFunction" pointer
                                    LastBasisFunction = (LastBasisFunction._next = new BasisFunction(1, 0, Centers[Current_Center_Number-1].R0, CoefNum));
                                    // update HighestL
                                    if (HighestL < 1) HighestL = 1;
                                     *
                                     */
                                } else
                                    if (HighestL < L) HighestL = L;
                            // now read CoefNum lines
                            for (int i=0; i<CoefNum; i++) {
                                s = in.readLine();
                                arr = s.trim().toUpperCase().replace("D","E").split(" +"); // convert '0.1171811300D+08' -> '0.1171811300E+08'                                
                                // this is basis function data
                                if (L == -1) {
                                    // a common exponents
                                    LastBasisFunction.exponents[i] = (bs_tmp.exponents[i] = Double.parseDouble(arr[0]));
                                    // but two coefs are allowed
                                    // s-type first
                                    bs_tmp.coefs[i] = Double.parseDouble(arr[1]);
                                    // and now a p- one
                                    LastBasisFunction.coefs[i] = Double.parseDouble(arr[2]);
                                } else {
                                    // for 'normal' functions
                                    LastBasisFunction.exponents[i] = Double.parseDouble(arr[0]);
                                    LastBasisFunction.coefs[i] = Double.parseDouble(arr[1]);
                                }
                            }
                        }
                    }
                    break;
                //----------------------------------------------
                case section_MO:
                    /* Section format
                     [MO]
                     Sym= symmetry_label_1
                     Ene= mo_energy_1
                     Spin= (Alpha|Beta)
                     Occup= mo_occupation_number_1
                     ao_number_1 mo_coefficient_1
                     ...
                     ao_number_n mo_coefficient_n
                     ....
                     Sym= symmetry_label_N
                     Ene= mo_energy_N
                     Spin= (Alpha|Beta)
                     Occup= mo_occupation_number_N
                     ao_number_1 mo_coefficient_1
                     ...
                     ao_number_n mo_coefficient_n
                    */
                    int Current_MO = 0;
                    NewSectionDetected = s.contains("[");
                    if (!NewSectionDetected) {
                        boolean ready = false;
                        while (!ready) {
                            // allocate memory
                            if (FirstMO == null)
                                LastMO = (FirstMO = new MO());
                            else
                                LastMO = (LastMO._next = new MO());
                            Current_MO++;
                            // allocate arrays
                            LastMO.BS_Coefs = new double[ Basis.length ]; // this is why _BS_List_To_Array should be called BEFORE "[MO]" section
                            // parse data
                            // 1) "prop=value" subsection
                            while ((s != null) && (s.contains("="))) {
                                String s2 = s.toLowerCase();
                                // parse orbital property
                                if (s2.contains("ene=")) LastMO.Energy = Double.parseDouble( s2.split("=")[1].trim() ); else
                                if (s2.contains("occup=")) LastMO.Occupancy = Double.parseDouble( s2.split("=")[1].trim() );
                                if (s2.contains("spin=")) {
                                    if (s2.contains("alpha")) LastMO.Spin = +1;
                                    if (s2.contains("beta")) LastMO.Spin = -1;
                                }
                                /*if (s.contains("Sym="))*/
                                // get next line
                                s = in.readLine();
                            }
                            if (print_MO_info)
                                out.printf("MO %d: Energy = %.7f, Occupancy = %.7f\n",Current_MO,LastMO.Energy,LastMO.Occupancy);
                            // data section
                            while ((s!=null) && (!s.contains("="))) {
                                if((s.trim().isEmpty())) { // MODIFIED this hole if statement
                                    s = in.readLine();
                                    ready = (s == null);
                                    if (!ready) {
                                        ready = NewSectionDetected = s.contains("[");}
                                    continue;}
                                arr = s.trim().split(" +");
                                LastMO.BS_Coefs[Integer.parseInt(arr[0])-1] = Double.parseDouble(arr[1]);
                                // get next line
                                s = in.readLine();
                                ready = (s == null);
                                if (!ready) {
                                    ready = NewSectionDetected = s.contains("[");
                                }
                            }
                        }
                        //System.out.println();
                    }
                    break;
                //----------------------------------------------
                default:
                    // wait until next section begins and parse its type
                    NewSectionDetected = s.contains("[");
                    break;
            } // end of switch
            //
            if (NewSectionDetected) {
                // post-processing of the previous section
                switch (section_NOW) {
                    case section_Atoms:
                        // Convert the list of atoms to an array
                        int NAtoms = 0;
                        // get number of atoms
                        LastAtom = FirstAtom;
                        while (LastAtom != null) {
                            NAtoms++;
                            LastAtom = LastAtom._next;
                        }
                        if (debugPrint) out.printf("There are %d atomic centers \n", NAtoms);
                        Centers = new AtomicCenter[NAtoms];
                        LastAtom = FirstAtom;
                        for (int i=0; i<NAtoms; i++) LastAtom = (Centers[i] = LastAtom)._next;
                        break;
                }
                // parse next section ID
                section_NOW = _Parse_Section_ID(s);
                // pre processing of the NEW section
                // some initializations might be neccessary BEFORE actual parsing of data
                /* typical examples include one-line sections like
                   [5D] (use of 'spherical' D and F functions: 5 D and 7 F), the default is 6 D and 10 F ('cartesian' functions )
                   [9G] to specify spherical G functions.
                   [5D10F] ,[7F] (6D en 7F),[5D7F] - the use of mixed spherical and cartesian function
                 */
                switch (section_NOW) {
                    case section_5D: case section_5D7F: // the same code for both cases!
                        use_5D = true;
                        use_7F = true; // [5D] sets both this keys!
                        section_NOW = section_unknown; //wait for next section
                        break;                    
                    case section_9G:
                        use_9G = true;
                        section_NOW = section_unknown; //wait for next section
                        break;
                    case section_Atoms:
                        // parse units
                        Coords_in_AU = s.toUpperCase().contains("AU"); // may contain "ANGS" as well
                        break;
                    case section_MO:
                        // if present, [5D] should go BEFORE "[MO]" section;
                        // "[GTO]" should go BEFORE "[MO]" section;

                        // Check whether all basis functions were denoted as "spherical" - this should be done BEFORE
                        // "cloning" the basis functions
                        IsSpherical = use_5D & use_7F;
                        if (HighestL >= 4) /* if G functions were used*/ IsSpherical &= use_9G;
                        out.println("Basis functions are spherical: "+IsSpherical);
                        
                        _BS_List_To_Array(); // Convert basis function list to an array and make clones m=(-L+1)...+L
                        break;
                    case section_MOLDEN: break; // do nothing
                    case section_unknown:
                        WarningManager.warning_printf("WARNING: section \"%s\" is not supported!%n",s);
                        break;
                }
            }
        }
        // close file
        in.close();
        fr.close();
        

        if (FirstMO == null) return false;   // this is not very good style, though...
        if (Basis == null) return false;  // this is not very good style, though...

        // Convert MO list to an array
        int MO_Count = 0;
        LastMO = FirstMO;
        while (LastMO != null) {
            MO_Count++;
            LastMO = LastMO._next;
        }
        out.printf("Total number of MO: %d\n",MO_Count);
        MOs = new MO[MO_Count];
        // convert a list to an array
        LastMO = FirstMO;
        for (int i=0; i<MO_Count; i++) LastMO = (MOs[i] = LastMO)._next;
        //

        return true;
    }
    //--------------------------------------------------------------------------
    
    // copies all fields from source
    public boolean CopyAllFrom(MOLDEN_IO source) {

        IsSpherical = source.IsSpherical;

        // re-create basis
        Basis = new BasisFunction[source.Basis.length];
        for (int bs=0; bs<Basis.length; bs++) Basis[bs] = new BasisFunction(source.Basis[bs]);
        // no need in copying FirstBasisFunction
        NBasisFunctions = source.NBasisFunctions;
        HighestL = source.HighestL;

        Allow_additional_r_power = source.Allow_additional_r_power;         // bugfix, Aug.2014

        RadialParts = new RadialPartOfBasisFunction[source.RadialParts.length];

        for (int i=0; i<RadialParts.length; i++)
            RadialParts[i] = source.RadialParts[i].CreateClone();
        /*
        // copy radial parts data
        RadialPartsOfBasisFunctions_Exponents = new double[source.RadialParts.length][];
         // .clone() does not work for 2D arrays
        RadialPartsOfBasisFunctions_Coefs = new double[source.RadialParts.length][];
        for (int i=0; i<RadialParts.length; i++)
            RadialParts[i].Coefs = source.RadialParts[i].Coefs.clone();
        RadialPartsOfBasisFunctions_CenterIDs = source.RadialPartsOfBasisFunctions_CenterIDs.clone();
        RadialPartsOfBasisFunctions_LUsedWith = source.RadialPartsOfBasisFunctions_LUsedWith.clone();
        RadialPartsOfBasisFunctions_Addit_r_power = source.RadialPartsOfBasisFunctions_Addit_r_power.clone();
         * 
         */

        Title = new String(source.Title);

        // copy MOs
        MOs = new MO[source.MOs.length];
        for (int mo=0; mo<MOs.length; mo++)
            MOs[mo] = new MO( source.MOs[mo] ); // re-create by copying everithing from source's MOs
        // no need in copying FirstMO 

        Coords_in_AU = source.Coords_in_AU;

        // re-create centers array        
        Centers = new AtomicCenter[source.Centers.length];
        // copy objects
        for (int cntr=0; cntr<Centers.length; cntr++) Centers[cntr] = new AtomicCenter(source.Centers[cntr]);
        // no need in copying FirstAtom
        
        use_5D = source.use_5D;
        use_7F = source.use_7F;
        use_9G = source.use_9G;

        return true;
    }
    //--------------------------------------------------------------------------
    /**
     * Converts basis function coefficients from normalized primitive gaussians to non-normalized
     * assuming that ||Y_LM||^2 == 1
     * 
     * E.g.,
     * { EXP = zeta, COEF = C }   => BF = C*exp(-zeta*r^2)/(Pi/zeta)^(3/4)   { EXP = zeta, COEF = C/(Pi/zeta)^(3/4) }
     * In fact, makes coef /= sqrt( exp(...)_norm2 )
     * 
     * Version: 15.Jan.2014
     *
     *
     */
    public void ToUnnormalizedPrimitiveCoefs() {
        for (int b=0; b<RadialParts.length; b++) {
            int L = RadialParts[b].LUsedWith;
            for (int cf=0; cf<RadialParts[b].Coefs.length; cf++)
                if (RadialParts[b].Addit_r_power != 0)
                    RadialParts[b].Coefs[cf] /= Math.sqrt(
                        OverlapIntegrals.primitive_int_1D_Sphr(2*L+2 + 2*RadialParts[b].Addit_r_power,
                            2*RadialParts[b].Exponents[cf]) );
                else
                    RadialParts[b].Coefs[cf] /= Math.sqrt(
                        OverlapIntegrals.primitive_int_1D_Sphr(2*L+2, 2*RadialParts[b].Exponents[cf]) );
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Converts basis function coefficients from un-normalized primitive gaussians to normalized ones
     * assuming that ||Y_LM||^2 == 1
     * 
     * E.g.,
     * { EXP = zeta, COEF = C }   => BF = C*exp(-zeta*r^2) = C*(Pi/zeta)^(3/4)*exp(-zeta*r^2)/(Pi/zeta)^(3/4)   { EXP = zeta, COEF = C*(Pi/zeta)^(3/4) }
     * In fact, makes coef *= sqrt( exp(...)_norm2 )
     *
     * Version: 15.Jan.2014
     *
     */
    public void ToNormalizedPrimitiveCoefs() {
        for (int b=0; b<RadialParts.length; b++) {
            int L = RadialParts[b].LUsedWith;
            for (int cf=0; cf<RadialParts[b].Coefs.length; cf++)
                if (RadialParts[b].Addit_r_power != 0)
                    RadialParts[b].Coefs[cf] *= Math.sqrt(
                        OverlapIntegrals.primitive_int_1D_Sphr(2*L+2+2*RadialParts[b].Addit_r_power,
                            2*RadialParts[b].Exponents[cf]) );
                else
                    RadialParts[b].Coefs[cf] *= Math.sqrt(
                        OverlapIntegrals.primitive_int_1D_Sphr(2*L+2, 2*RadialParts[b].Exponents[cf]) );
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the values of the orbitals selected by a flag in DoEval array at the point r
     * If DoEval == null, all orbitals are evaluated
     * 
     * Version: 08.Jan.2014
     *
     */
    public double[] EvaluateOrbitalsAtPoint(boolean[] DoEval, double[] r) {
        double[] result = new double[MOs.length];
        // evaluate all basis functions first
        double[] BS_values = new double[Basis.length];
        for (int bf=0; bf<Basis.length; bf++)
            BS_values[bf] = Basis[bf].EvaluateAtPoint(r, Basis[bf].R0);
        
        for (int mo=0; mo<MOs.length; mo++) {
            result[mo] = 0;
            if ( (DoEval == null) || ((DoEval != null) && (DoEval[mo])) )
                for (int bf=0; bf<Basis.length; bf++)
                    result[mo] += BS_values[bf] * MOs[mo].BS_Coefs[bf];
        }
        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Uses EvaluateOrbitalsAtPoint to get the electron density at a given point r
     *
     * Version: 08.Jan.2014
     *
     */
    public double EvaluateDensityAtPoint(double[] r) {
        // evaluate all MOs
        double[] MO_values = EvaluateOrbitalsAtPoint(null, r);
        // get the density
        double result = 0;
        for (int mo=0; mo<MOs.length; mo++)
            result += MOs[mo].Occupancy * MO_values[mo]*MO_values[mo];

        return result;
    }

    //--------------------------------------------------------------------------
    // used for printing MO coefs
    public String MOCoefLineFormat = "%3d %20.12f%n";
    //--------------------------------------------------------------------------
    /**
     * Saves basis functions and MO coefs from MOLDEN file
     * Rev.: 15.Jan.2014
     *
     */
    public boolean Save_As_MOLDEN(String filename /*, double[] L_renormalizers*/) throws Exception {        
        PrintWriter out = new PrintWriter(filename);
        //PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename)), 1024*1024));
        out.println("[Molden Format]");

        // [Title] section
        out.println("[Title]");
        out.println(Title);
        out.println();

        // [Atoms] section
        if (Coords_in_AU) out.println("[Atoms] AU");
            else out.println("[Atoms] Angs");
        // print centers info
        for (int cntr=0; cntr<Centers.length; cntr++)
            out.printf("%-2s%4d %3.0f  %20.10f %20.10f %20.10f %n",
                    Centers[cntr].Name, cntr+1, Centers[cntr].Z,
                    Centers[cntr].R0[0], Centers[cntr].R0[1], Centers[cntr].R0[2]);

        // [GTO] section
        out.println("[GTO]");
        // loop over all centers
        for (int cntr=0; cntr<Centers.length; cntr++) {
            out.printf("%3d 0 %n",cntr+1);
            // loop over all radial parts and pick up those belonging to this center
            for (int rp=0; rp<RadialParts.length; rp++)
                if (RadialParts[rp].CenterID == (cntr+1)) {
                    // print this radial part
                    // header
                    if (RadialParts[rp].Addit_r_power == 0)
                        out.printf("%s%4d 1.0  %n",
                                "spdfgh".charAt(RadialParts[rp].LUsedWith),
                                RadialParts[rp].Coefs.length);
                    else {
                        if (!Allow_additional_r_power)
                            throw new Exception(" Radial part #"+rp+" attempts using additional_r_power>0 but Allow_additional_r_power is false!");
                        out.printf("%s%d%4d 1.0  %n",
                                "spdfgh".charAt(RadialParts[rp].LUsedWith),
                                RadialParts[rp].Addit_r_power,
                                RadialParts[rp].Coefs.length);
                    }

                    
                    // print exponents and coefs
                    for (int cf=0; cf<RadialParts[rp].Coefs.length; cf++)
                        out.printf("%20.10f %20.10f %n",
                                RadialParts[rp].Exponents[cf],
                                RadialParts[rp].Coefs[cf] /* * correction */);
                }
           out.println(); // end with a blank line
        }
        
        if (IsSpherical) {
            out.println("[5D]");
            out.println("[9G]");
        }

        // [MO] section
        if (MOs != null) {
            out.println("[MO]");
            for (int mo=0; mo<MOs.length; mo++) {
                // print MO header
                out.printf(" Sym=     1a%n");
                out.printf(" Ene=%22.14E%n", MOs[mo].Energy); // MOs[mo].Energy
                if (MOs[mo].Spin == +1)
                    out.printf(" Spin= Alpha%n");
                else
                    out.printf(" Spin= Beta%n");

                out.printf(" Occup=%9.6f%n", MOs[mo].Occupancy);
                // print MO coefs
                for (int cf=0; cf<MOs[mo].BS_Coefs.length; cf++)
                    out.printf(MOCoefLineFormat, cf+1, MOs[mo].BS_Coefs[cf]);
            }
        }
        
        out.close();
        //
        return true;
    }
    //--------------------------------------------------------------------------
    final int bmolden_section_FLAGS            = 1;
    // int NumberOfFlags | [FlagsID(4 bytes)+ FlagValue(4 bytes)=8 bytes] ...
    final int bmolden_section_GEOMETRY         = 2;
    // int NumbrOfCenters | Center1 ...
    final int bmolden_section_BASIS            = 3;
    // int NumberOfRadialParts | RadialPart1 ...
    final int bmolden_section_MO               = 4;
    // int NumberOfMO | MO_1 ....
    
    final int bmolden_flagID_CoordsAUAngstrs   = 11;
    final int bmolden_flagID_5D                = 23;
    final int bmolden_flagID_7F                = 24;
    final int bmolden_flagID_9G                = 25;

    final static String bmoldeng_signature = "BMOLDENG"; // where G stands for 'G'aussian
    //--------------------------------------------------------------------------
    /**
     * Saves the data as a 'binary' MOLDEN for space economy
     * 
     * @return true if success
     * @version 20.Jul.2016
     */
    public boolean Save_as_BMOLDEN(String filename) throws Exception {
        FileOutputStream fout = new FileOutputStream(filename);
        //DataInputStream x;        
        DataOutputStream out = new DataOutputStream(fout);
        // file signature
        out.writeChars(bmoldeng_signature); // "BMOLDENG"

        // TITLE
        out.writeInt(Title.length());
        out.writeChars(Title);
        
        // FLAGS
        out.writeInt(bmolden_section_FLAGS);
        out.writeInt(4); // there will be 4 flags
        // Bohrs/Angstroms
        out.writeInt(bmolden_flagID_CoordsAUAngstrs);        
        if (Coords_in_AU) out.writeInt(1); else out.writeInt(0);
        // whether to use 5D
        out.writeInt(bmolden_flagID_5D);
        if (use_5D) out.writeInt(1); else out.writeInt(0);
        // whether to use 7F
        out.writeInt(bmolden_flagID_7F);
        if (use_7F) out.writeInt(1); else out.writeInt(0);
        // whether to use 9G
        out.writeInt(bmolden_flagID_9G);
        if (use_9G) out.writeInt(1); else out.writeInt(0);

        // GEOMETRY
        out.writeInt(bmolden_section_GEOMETRY);
        out.writeInt(Centers.length);
        for (int cntr=0; cntr<Centers.length; cntr++)
            Centers[cntr].SaveToDataStream(out);

        // BASIS: save all radial parts
        out.writeInt(bmolden_section_BASIS);
        out.writeInt(RadialParts.length);
        for (int rp=0; rp<RadialParts.length; rp++)
            RadialParts[rp].SaveToDataStream(out);

        // MO
        out.writeInt(bmolden_section_MO);
        out.writeInt(MOs.length);
        for (int mo=0; mo<MOs.length; mo++)
            MOs[mo].SaveToDataStream(out);
        

        out.close();
        fout.close();

        return true;
    }
    //--------------------------------------------------------------------------
    /**
     * READS the data from a 'binary' MOLDEN
     *
     * @return true if success
     * @version 17.05.2014
     */
    public boolean Load_from_BMOLDEN(String filename) throws Exception {
        FileInputStream fin = new FileInputStream(filename);
        DataInputStream in = new DataInputStream(fin);


        // check signature
        char[] signature = new char[ in.readInt() ];
        for (int i=0; i<signature.length; ++i) signature[i] = in.readChar();
        if (!bmoldeng_signature.equals(signature))
            return false;

        // read TITLE
        char[] title = new char[ in.readInt() ];
        for (int i=0; i<title.length; ++i) title[i] = in.readChar();
        Title = String.valueOf(title);

        // Parse the data for all the other sections
        
        // DataInputStream reports about the end of file by EOFException exception
        try {
            while (true) {
                int section_code = in.readInt();
                int ElementNumber = in.readInt();
                switch(section_code) {
                    case bmolden_section_FLAGS:
                        // we expect ElementNumber pairs of (int FlagID, int FlagValue)
                        for (int i=0; i<ElementNumber; ++i) {
                            // set appropriate flags
                            int flagID = in.readInt();
                            int flagValue = in.readInt();
                            switch(flagID) {
                                case bmolden_flagID_CoordsAUAngstrs:
                                    Coords_in_AU = (flagValue != 0);
                                    break;
                                case bmolden_flagID_5D:
                                    use_5D = (flagValue != 0);
                                    break;
                                case bmolden_flagID_7F:
                                    use_7F = (flagValue != 0);
                                    break;
                                case bmolden_flagID_9G:
                                    use_9G = (flagValue != 0);
                                    break;
                                default:
                                    WarningManager.warning_printf("WARNING: flag ID %d unknown! Ignored!%n",flagID);
                            }
                        }
                        break;
                    case bmolden_section_GEOMETRY:
                        // we expect ElementNumber of AtomicCenter objects
                        Centers = new AtomicCenter[ElementNumber];
                        for (int i=0; i<ElementNumber; ++i)
                            ( Centers[i] = new AtomicCenter() ).ReadFromDataStream(in);
                        break;
                    case bmolden_section_BASIS:
                        // we expect ElementNumber of RadialPartOfBasisFunction objects
                        RadialParts = new RadialPartOfBasisFunction[ElementNumber];
                        for (int i=0; i<ElementNumber; ++i)
                            ( RadialParts[i] = new RadialPartOfBasisFunction() ).ReadFromDataStream(in);
                        
                        // Do we have any radial part with Addit_r_power>0 ?
                        // and what is the largest L ?
                        boolean a_r_found = false;
                        HighestL = 0;
                        for (int j=0; (j<ElementNumber); ++j) {
                            a_r_found |= (RadialParts[j].Addit_r_power > 0);
                            if (RadialParts[j].LUsedWith > HighestL) HighestL = RadialParts[j].LUsedWith;
                        }
                        if ((!Allow_additional_r_power) && (a_r_found))
                                throw new Exception("Addit_r_power > 0 found for one of radial parts, but Allow_additional_r_power is false.");
                        
                        // 'compute' IsSpherical flag
                        IsSpherical = use_5D & use_7F;
                        if (HighestL >= 4) /* if G functions were used*/ IsSpherical &= use_9G;
                        out.println("Basis functions are spherical: "+IsSpherical);

                        // convert radial part list to 'faky' basis function list, and use this._BS_List_To_Array() to
                        // build an actual list of basis functions                        
                        FirstBasisFunction = null; // free previous allocation
                        BasisFunction tmp = null, new_bf; // last element of the list
                        for (int rp=0; rp<RadialParts.length; ++rp) {
                            new_bf = new BasisFunction(
                                    RadialParts[rp].LUsedWith,
                                    0,
                                    Centers[ RadialParts[rp].CenterID-1 ].R0,
                                    RadialParts[rp].Exponents.length);
                            
                            if (FirstBasisFunction == null)
                                // create list
                                tmp = (FirstBasisFunction = new_bf);
                            else
                                // append to list
                                tmp = (tmp._next =  new_bf);
                            
                        }
                        _BS_List_To_Array();

                        break;
                    case bmolden_section_MO:
                        // we expect ElementNumber of MO objects
                        MOs = new MO[ElementNumber];
                        for (int i=0; i<ElementNumber; ++i)
                            ( MOs[i] = new MO() ).ReadFromDataStream(in);
                        break;
                }
            }
        } catch (EOFException e) {
            // end of file reached
            in.close();
            fin.close();
            return true;
        }

        //in.close();
        //fin.close();
        //return false;
    }
    //--------------------------------------------------------------------------

    /** Multiplies coordinates to by a given factor
     * Does NOT affect exponents and coefs. of basis functions!
     * Version: 12.Jan.2014
     */
    public void ScaleCoords(double factor)  {
        //        
        for (int cntr=0; cntr<Centers.length; cntr++)
            for (int mu=0; mu<3; mu++) Centers[cntr].R0[mu] *= factor ;
        // we do not change basis functions, but DO change their R0s
        if (Basis != null)
            for (int bs=0; bs<Basis.length; bs++)
                Basis[bs].R0 = Centers[Basis[bs].Center_ID-1].R0; // this is needed since the constructor of BasisFunction _clones_ (not simply assigns) R0
        // there is nothing to be done with radial parts of basis functions
    }
    //--------------------------------------------------------------------------
    public final static double BohrRadius =  0.52917721092;
    //--------------------------------------------------------------------------
    /** Converts coordinates to atomic units; does nothing if(Coords_in_AU)
     * Does NOT affect exponents and coefs. of basis functions!
     * Version: 12.Jan.2014
     */
    public void CoordsToAU()  {
        if (Coords_in_AU) return;
        ScaleCoords(1.0/BohrRadius); // 1 A = 1/0.529... a.u.
        Coords_in_AU = true;
    }
    //--------------------------------------------------------------------------
    /** Converts coordinates to Angstroms; does nothing if(!Coords_in_AU)
     * Does NOT affect exponents and coefs. of basis functions!
     * Version: 12.Jan.2014
     */
    public void CoordsToAngstroms()  {
        if (!Coords_in_AU) return;
        ScaleCoords(BohrRadius); // 1 a.u. = 0.529... A
        Coords_in_AU = false;
    }    
    //--------------------------------------------------------------------------
    /**
     * Converts contraction coefficients from the ones suitable for a
     * conventional MOLDEN basis function definitions, i.e.,
     * 
     *           (                    Ylm*exp(-zeta*(r-R0)^2)           )
     * BF = SUM_j(  C_j * --------------------------------------------- )
     *           (         SQRT( INT( {Ylm*exp(-zeta*r^2)}^2 * d3r ) )  )
     *
     * to the contraction coefs suitable for a simpler definition
     *
     * BF = SUM_j(  A_j * Ylm*exp(-zeta*(r-R0)^2) )
     *
     * where int( Ylm^2, d_Omega ) = 4*Pi * Ylm_norms2[l][m].
     *
     * Interconnection formula:
     *                              C_j
     * A_j =  ---------------------------------------------
     *         SQRT( INT( {Ylm*exp(-zeta*r^2)}^2 * d3r ) )
     *
     * Affects both radial parts and basis functions.
     * Note that for radial parts 4*Pi*Ylm_norms2_over4Pi[L][0] is taken as a square of the spherical function norm
     *
     * Version: 15.Jan.2014
     */
    public void UnNormalizePrimitives(double[][] Ylm_norms2_over4Pi) {
        if (!this.IsSpherical) {
            out.println("ERROR: MOLDEN_IO.UnnormalizePrimitives() currently works with spherical basis sets only!");
            return;
        }

        // Convert radial parts
        double Rnorm2 = 0;
        for (int rp=0; rp<RadialParts.length; ++rp) {
            int L = RadialParts[rp].LUsedWith;
            // loop over primitive gaussians in this basis function radial part
            for (int cf=0; cf<RadialParts[rp].Coefs.length; ++cf) {
                //
                //  Rnorm2 = INT(  {exp(-zeta*r^2) * r^L}^2 * r^2, r=0..+infinity  )
                //
                if (RadialParts[rp].Addit_r_power != 0) {                    
                    Rnorm2 = OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2 + 2*RadialParts[rp].Addit_r_power,
                            2 * RadialParts[rp].Exponents[cf]);
                } else
                    Rnorm2 = OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2, 2 * RadialParts[rp].Exponents[cf]);

                Rnorm2 *= 4*Math.PI; // since 4*Pi in NOT included into Ylm_norms2_over4Pi
                
                // Now Rnorm2 = 4*Pi / (2*L+1) * INT( {exp(-zeta*r^2) * r^L}^2 * r^2, r=0..+infinity  )
                // Finally, correct the contraction coefficient:
                RadialParts[rp].Coefs[cf] /= Math.sqrt(Rnorm2 *  Ylm_norms2_over4Pi[L][0]);
            }
        }
        
        // Convert basis functions
        for (int bf=0; bf<Basis.length; ++bf) {
            int L = Basis[bf].L;
            // loop over primitive gaussians in this basis function radial part
            for (int cf=0; cf<Basis[bf].coefs.length; ++cf) {
                Rnorm2 = OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2 + 2*Basis[bf].additional_r_power, 2 * Basis[bf].exponents[cf]);
                Rnorm2 *= 4*Math.PI; // since 4*Pi in NOT included into Ylm_norms2_over4Pi
                // Now Rnorm2 = 4*Pi / (2*L+1) * INT( {exp(-zeta*r^2) * r^L}^2 * r^2, r=0..+infinity  )
                // Finally, correct the contraction coefficient:
                Basis[bf].coefs[cf] /= Math.sqrt(Rnorm2 *  Ylm_norms2_over4Pi[L][ L+Basis[bf].m ]);
            }
        }        
    }
    //--------------------------------------------------------------------------
    /**
     * @return an array[mo][bf] of references to/copies of MOs[mo].BS_coefs
     * Note that if createBSCoefsCopy==false, BS_coefs[] are inserted as references!
     *
     * Ver.: 24.Jan.2015
     */
    public double[][] get_MO2AO_array(boolean createBSCoefsCopy) {
        double[][] result = new double[MOs.length][];        
        for (int mo=0; mo<MOs.length; mo++)
            if (createBSCoefsCopy) {
                result[mo] = new double[Basis.length];
                System.arraycopy(MOs[mo].BS_Coefs, 0, result[mo], 0, Basis.length);
            } else
                result[mo] = MOs[mo].BS_Coefs;
        return result;
    }
    //--------------------------------------------------------------------------

    /**
     * Converts MO from Cartesian basis to Spherical (if not already in it)
     */
    //void Ensure_Spherical() {
    //}
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //
    /*
    public void Compress() throws Exception {
        //int mo = 25;
        Double d;
        FileOutputStream f1 = new FileOutputStream("K:\\abInitio\\nucleosides_Roman_PSI4\\Reopt_dC_6-31Gdp\\single_points\\cc-pvtz\\zip\\zip\\A_TEST");
        FileOutputStream f2 = new FileOutputStream("K:\\abInitio\\nucleosides_Roman_PSI4\\Reopt_dC_6-31Gdp\\single_points\\cc-pvtz\\zip\\zip\\A_TEST2");

        //FileStream f = new FileStream()
        //OutputStreamWriter o = new OutputStreamWriter()
        //OutputStream o = new OutputStream()
        //PrintStream x = new PrintStream("K:\\abInitio\\nucleosides_Roman_PSI4\\Reopt_dC_6-31Gdp\\single_points\\cc-pvtz\\zip\\zip\\A_TEST");
        long L;
        int counter = 0;
        byte buffer_M=0, buffer_E=0, b;
        for (int mo=0; mo<MOs.length; ++mo)
         for (int bs=0; bs<Basis.length; ++bs){
            L = Double.doubleToRawLongBits(MOs[mo].BS_Coefs[bs]); // MSB [sign(1 bit)][exponent(11 bits)][mantissa(52 bit)] LSB
                                                                   // Seeeeeee eeeeMMMM MMMMMMMM MMMMMMMM MMMMMMMM MMMMMMMM MMMMMMMM MMMMMMMM
             // save first 6 bytes (mantissa)
             for (int i=0; i<6; ++i) {
                 f1.write( (byte)L );
                 L >>= 8;
             }
             // split 7-th to buffer_E / buffer_M

            if (false) {
                 b = (byte)L;
                 buffer_M <<= 4;
                 buffer_M |= (b & 0x0F);
                 buffer_E >>= 4;
                 buffer_E |= (b & 0xF0);
                 ++counter;
                 if ((counter %2)==0) {
                     f1.write(buffer_M);
                     f2.write(buffer_E);
                 }
                 
                 L >>= 8;
                 f2.write( (byte)L );
             } else {
                 f1.write( (byte)L );
                 L >>= 8;
                 f1.write((byte)L);
             }
             //f1.write( (byte)L );
             //x.println( Double.toHexString(MOs[mo].BS_Coefs[bs]) );
         }
        
         //f1.write(buffer_M);
         f2.write(buffer_E);

        f1.close();
        f2.close();
        return;
    }
     * 
     */
    //--------------------------------------------------------------------------

}
