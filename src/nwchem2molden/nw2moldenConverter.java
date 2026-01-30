package nwchem2molden;

import moldenio.*;
import java.io.*;
import Polynom3D.*;
import JGints.*;

/**
 * The workhorse class for nwchem2molden program
 * @author (c) Tymofii Nikolaienko, 2014
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
public class nw2moldenConverter {

    public MOLDEN_IO molden = new MOLDEN_IO();
    public boolean do_basis_comparison = true; // whether to compare basis set coefs/exponents loaded from -nwbasis file with the ones from nwchem log file

    //--------------------------------------------------------------------------
    // Reads in the numebr of lines equal to fmts.length from logfile and
    // returns all the strings read if all lines match regexps given in fmts,
    // or null if the strings did nto match specified formats
    private String[] _matches_fmts(BufferedReader logfile, String[] fmts) throws Exception{
        //
        String[] result = new String[fmts.length];
        String s;
        for (int i=0; i<fmts.length; i++) {
            s = logfile.readLine();
            if (s == null) return null;     // if no line is avaibable
            if (s.matches(fmts[i]))
                result[i] = s;
            else
                return null; // if s has wrong format
        }
        // if we got here, everuthing is fine!
        return result;
    }
    //--------------------------------------------------------------------------
    /* Expects reading
     * "  H (Hydrogen)
     *    ------------
     *              Exponent  Coefficients
     *         -------------- --------------------------------------------------------- "
     * Returns short eletemt name if successfull
     *
     */
    private String FourLinesHeaderParse(BufferedReader logfile, String FirstLine) throws Exception{
        boolean good = true;
        String result;
        String s;
        if (FirstLine == null) {
            s = logfile.readLine();
            if (s == null) return null;     // if no line is avaibable
        } else
            s = FirstLine;      // first line may be available from the previous reads
        //
        if (!s.matches(" +[a-zA-Z]+ \\(.+\\)")) return null; // if s has wrong format (different from " X(...)")
        result = s.split("\\(")[0].trim(); // short name of the element
        // check whether next three lines are of expected type:
        String[] fmts = new String[]{" +\\-+" /* some "---"s*/,
            " +Exponent +Coefficients +",
            " +\\-+ +\\-+" /* some "---"s, " " and some "---"s */}; // regexps for expected formats
        if (_matches_fmts(logfile, fmts) == null) return null;
        // if we got here, everything was fine!
        return result;
    }
    //--------------------------------------------------------------------------
    /* Reads and parses information like
     * "  1 S  1.30100000E+01  0.019685
          1 S  1.96200000E+00  0.137977
          1 S  4.44600000E-01  0.478148

          2 S  1.22000000E-01  1.000000

          3 P  7.27000000E-01  1.000000 "
     *
     * Reads until a line not matching "3 P  7.27000000E-01  1.000000"-like format is read;
     * returns such a "strange" line (the last one read).
     * In fact, more than one basis function is being read in a single call; instead,
     * the whole set of basis functions for a current center is to be read.
     * All basis functions read are stored in bs, and marked with Z as a center_ID
     *
     * Rev.: 05.02.2014
     */
    private String ParseElementBasis(BufferedReader logfile, int Z, BasisFunction bs) throws Exception{
        String s = null;
        final int MAX_BS_LENGTH = 1000; // maximum length of a radial part expansion
        double[] coefs = new double[MAX_BS_LENGTH]; // working arrays
        double[] expons = new double[MAX_BS_LENGTH];
        int NTerms = 0;

        boolean finished = false;
        String BS_FMT = " +[0-9]+ +[SPDFGHI]+ +[E0-9\\.\\-\\+\\,]+ +[0-9\\.\\,\\-\\+]+";
        boolean create_new_BS = false;
        while (!finished) {
            // read line
            finished = ((s = logfile.readLine()) == null); // read the first line
            if (!finished) {
                // if not end of file
                if (s.equals("")) { // it might be better to do s.matches(" +")
                    finished = create_new_BS; // a trick to "catch" two blank lines one after another:
                        // a previous blank line would set create_new_BS to true, while a non-blank
                        // line belonging to basis would set it false;
                    if (!finished) {
                        // this line is empty => some set of exponents has finished;
                        // save previously read data into basis function object
                        bs.coefs = new double[NTerms];
                        bs.exponents = new double[NTerms];
                        System.arraycopy(expons, 0, bs.exponents, 0, NTerms);
                        System.arraycopy(coefs, 0, bs.coefs, 0, NTerms);
                        NTerms = 0; // !!! 'clear' temporary buffers!
                        create_new_BS = true; // create new object before using 'bs' variable; not needed for the very first use!!!
                    }
                } else {
                    // this line is not empty;
                    // it may or may not belong to the basis set data
                    finished = !s.matches(BS_FMT);
                    if (!finished) {
                        // this line belongs to the basis set data
                        // creation of a new object in bs variable might be needed
                        if (create_new_BS) {
                            // trick: this will not be executed if s is already outside teh basis
                            bs._next = new BasisFunction(null);
                            bs._next.IsSpherical = bs.IsSpherical; // inherit this property from the previous function in the list
                            bs = bs._next; // make bs point to the newly created function
                        }
                        create_new_BS = false; // do not create anything up to an empty line occasion

                        // parse this line and store in bs
                        bs.Center_ID = Z;
                        String params[] = s.split(" +");
                        bs.L = "SPDFGHI".indexOf(params[2].toUpperCase());
                        if (params[2].length() > 1) {
                            System.out.printf(" ERROR: \"%s\" type of basis functions is unknown!%n", params[2]); // for "SP" and "SPD" cases
                        }
                        expons[NTerms] = Double.parseDouble(params[3]);
                        coefs[NTerms] = Double.parseDouble(params[4]);
                        NTerms++;
                    }
                }

            }
            /*
            // analyze this line
            if ((!finished) && (!s.equals(""))) {
            } else {
                // line s is empty
                // Read the next line
                finished = ((s = logfile.readLine()) == null);
                if (!finished) {
                    // does this line still belong to a basis set ?
                    finished = (!s.matches(BS_FMT));
                }
                if (!finished) {
                    // we're still inside a basis set section;
                    // save previous data on exps and coefs
                }
            }
             *
             */
        }
        return s; // return the last line read!
    }
    //--------------------------------------------------------------------------
    public BasisSetLoader Extended_Basis_Info = null;
    //--------------------------------------------------------------------------
    /**
     * Tries to extract basis set and geometry information from the NWChem logfile;
     * Stores results into molden object
     * Sets molden 'title' field to the name of the basis set used
     *
     * @param fname
     * Rev.: 05.02.2014
     */

    public boolean _Parse_OUT_File(String fname, MOLDEN_IO molden) throws Exception {
        FileReader fr1 = new FileReader(fname);
        BufferedReader logfile = new BufferedReader(fr1);
        //
        String s;
        String basis_title = null;

        boolean basis_cart = false;
        boolean basis_spher = false;
        int NAtoms = 0;
        AtomicCenter First_Atom = null; // first center in the list
        AtomicCenter _atom = null;
        BasisFunction First_BS = null; // first basis function in the list
        BasisFunction bs = null; // temporary variable
        //boolean coords_in_au = false; // default: a.u.
        double coords_to_au = 1.0; // a factor for converting coords. to a.u. (printed by NwChem)

        // which algorithm to use to get geometry record:
        final boolean XYZ_geom_label = false;                   // this will take an INITIAL geometry only!!!! -> DO NOT USE THIS!
        final boolean GeomGeom_geom_label = true;


        /*
         * What should we do reading file is to find and extract the basis set information
         * and to find and extract the geometry information;
         * Only the most last basis set and geometry data are stored.
         *
         */
        while ( (s = logfile.readLine()) != null) {
            // TODO: detect "ecp_print" and/or "so_print"

            // Detect basis set information by ' Basis "ao basis" -> "ao basis" (cartesian)'-like line

            basis_cart = s.matches(" +Basis \\\".+\\\" -> \\\".+\\\" \\(cartesian\\)");  // regexp: some spaces, "Basis", some name, "->", some name, "(cartesian)"
            basis_spher = s.matches(" +Basis \\\".+\\\" -> \\\".+\\\" \\(spherical\\)");  // regexp: some spaces, "Basis", some name, "->", some name, "(spherical)"
            // in both cases next line should be "                      -----"-like
            if (basis_cart || basis_spher) {
                basis_title = s;
                s = logfile.readLine();
                if (!s.matches(" +\\-+")) basis_cart  = ( basis_spher = false ); // if no "----"-like line found
            }

            // when the basis is found:
            if (basis_cart) {
                System.out.println(" CARTESIAN basis set information found:");
                molden.IsSpherical = false;
            }
            if (basis_spher) {
                System.out.println(" SPHERICAL ('PURE') basis set information found:");
                molden.IsSpherical = true;
            }
            //
            if (basis_cart || basis_spher) {
                System.out.println(basis_title); // cite basis set information

                System.out.println(" reading basis set information...");

                // release previous basis set
                First_BS = null;

                String first_line = null;
                bs = null; // temporary variable
                String at_name;
                do {
                    at_name = FourLinesHeaderParse(logfile, first_line); // try reading basis header for atom
                    // if at_name == null, those lines are already outside the basis set section
                    if (at_name != null) {
                        System.out.printf(" Basis for %s%n", at_name);

                        // allocate memory
                        if (bs == null)
                            First_BS = (bs = new BasisFunction(null));
                        else
                            bs = (bs._next = new BasisFunction(null)); // create new object in ._next field and save it to bs
                        // set bs.IsSpherical to a proper value
                        bs.IsSpherical = molden.IsSpherical;

                        // parse BS information
                        first_line = ParseElementBasis(logfile, PeriodicTableData.Name_to_Z(at_name), bs);
                        // now bs._next field has been modified, but bs itself has been not;
                        // so, go to the end of the list
                        while ((bs._next) != null) bs = bs._next;
                    }
                } while (at_name != null);

                s = first_line; // update s with the last line read before calling FourLinesHeaderParse()

            }
            /////////////////////////////////////////////////////////////////////
            // are we about to find geometry information?

            ///////////////////////////////////////////////////////////////////
            // the most simple way
            if (s.trim().matches("XYZ format geometry") && XYZ_geom_label) {
                s = logfile.readLine();
                if ((s != null) && s.matches(" +\\-+")) {
                    // we are in geometry section
                    NAtoms = Integer.parseInt( logfile.readLine().trim() );
                    molden.Centers = new AtomicCenter[NAtoms];
                    System.out.printf(" Reading geometry \"%s\"%n", logfile.readLine());
                    // write data directly to molden onject
                    for (int i=0; i<NAtoms; i++) {
                        // line format: " C                    -0.00209924     0.38028059    -0.32508001 "
                        String[] vals = logfile.readLine().trim().split(" +");
                        molden.Centers[i] = new AtomicCenter();
                        molden.Centers[i].ID = i+1;
                        // name and Z
                        molden.Centers[i].Name = vals[0];
                        molden.Centers[i].Z = PeriodicTableData.Name_to_Z(vals[0]);
                        // coords
                        for (int mu=0; mu<3; mu++) molden.Centers[i].R0[mu] = Double.parseDouble(vals[mu+1]);
                        // some tricks for compatibility between the algorithms
                        if (i==0) First_Atom = molden.Centers[i];
                            else  molden.Centers[i-1]._next = molden.Centers[i];
                    }
                    molden.Coords_in_AU = false;
                }
            }
            ///////////////////////////////////////////////////////////////////

            // more complicated way
            if (s.matches(" +Geometry \\\".+\\\" -> \\\".*\\\"") && (GeomGeom_geom_label)) { // something like (Geometry "geometry" -> "") is also acceptable
                // is it realy 'Geometry' record ?
                /* We expect finding something like this:
                 "                             -------------------------

                 Output coordinates in angstroms (scale by  1.889725989 to convert to a.u.)

                  No.       Tag          Charge          X              Y              Z
                 ---- ---------------- ---------- -------------- -------------- -------------- "
                 *
                 */
                String[] fmts = new String[]{" +\\-+" /* "---"s */,
                    "" /*emty line*/,
                    " +Output coordinates in.+",
                    "" /*emty line*/,
                    " +No.+ +Tag.+ +Charge.+ +X.+ +Y.+ +Z",
                    " +\\-+ +\\-+ +\\-+ +\\-+ +\\-+ +\\-+"};
                String[] gm_header = _matches_fmts(logfile, fmts);
                if (gm_header != null) {
                    System.out.println(" Geometry information found ");
                    // get geometry units from 2-nd line of gm_header
                    System.out.printf(" geomrtey units: \"%s\"%n", gm_header[2]);
                    // are coords. in "angstroms" ?
                    //coords_in_au = true;
                    try {
                        String conv_factor = null;
                        conv_factor = gm_header[2].toLowerCase().split("\\(scale by")[1].split("to convert to a\\.u\\.")[0];
                        conv_factor = conv_factor.trim();
                        coords_to_au = Double.parseDouble(conv_factor);
                    } catch (Exception e) {
                        coords_to_au = 1.0;
                        //coords_in_au = !gm_header[2].toLowerCase().contains("Output coordinates in angstroms");
                        System.out.println("WARNING: failed to read a coefficient for transformation of coords. to a.u.");
                        System.out.println(" a value of 1.0 is assumed");
                    }

                    /* Possible choices are:
                        angstroms or an -- Angstroms , the default (converts to A.U. using the Angstrom to A.U. conversion factor)
                        au or atomic or bohr -- Atomic units (A.U.)
                        nm or nanometers -- nanometers (converts to A.U. using a conversion factor computed as 10.0 times the Angstrom to A.U. conversion factor)
                        pm or picometers -- picometers (converts to A.U. using a conversion factor computed as 0.01 times the Angstrom to A.U. conversion factor)
                     */

                    // release previous list of atoms
                    First_Atom = null;
                    NAtoms = 0;

                    // a temporary variable:
                    _atom = null;

                    // read until an empty line
                    while (!(s = logfile.readLine()).equals("")) {
                        // Create an object to store information about the center
                        if (_atom == null)
                            // if this is the first atom:
                            First_Atom = (_atom = new AtomicCenter());
                        else
                            _atom = ( _atom._next = new AtomicCenter()); // create new object in ._next field and assign it to _atom
                        //
                        // parse line
                        // each line has "    1 H                    1.0000     0.25448866     1.32107362     0.00000000"-like format
                        String[] props = s.split(" +"); // be careful: props[0]==""
                        _atom.ID = Integer.parseInt(props[1]);
                        _atom.Name = props[2];
                        _atom.Z = Double.parseDouble(props[3]);
                        for (int mu=0; mu<3; mu++) _atom.R0[mu] = Double.parseDouble(props[4+mu]); // parse coords.
                        NAtoms++;
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////
        }
        //
        logfile.close();
        fr1.close();

        System.out.println(" Reading finished.");

        if (First_Atom == null) {
            System.out.println("ERROR: no information about molecular geometry has been found!");
            return false;
        }
        System.out.printf(" %d atomic centers found%n", NAtoms);

        if (First_BS == null) {
            System.out.println("ERROR: no information about the basis set has been found!");
            return false;
        }

        // extract basis set name from basis_title string
        String bs_name = basis_title.split("\\\"")[1].split("\\\"")[0].trim();
        molden.Title = bs_name;

        if (GeomGeom_geom_label) {
            // trasfer the data to molden object:
            // a list of atoms
            molden.Centers = new AtomicCenter[NAtoms];
            molden.Coords_in_AU = true;
            _atom = First_Atom;
            for (int i=0; i<NAtoms; i++) {
                molden.Centers[i] = _atom;
                // convert coords. to a.u.
                for (int mu=0; mu<3; mu++) molden.Centers[i].R0[mu] *= coords_to_au;
                _atom = _atom._next;
            }
        }
        // molden already contains gemoetry data if XYZ-simple algorithms has been used

        //  molden.IsSpherical = basis_spher;   - already done above
        // TODO: check for 7F key...


        if ( Extended_Basis_Info != null ) {
            System.out.println("Using extended information about the basis set from an external source.");
            // replace coefs. of BS with those from Extended_Basis_Info
            bs = First_BS;
            int[] indxs = new int [ PeriodicTableData.all_elements.length +1]; // last function used for this Z
            while (bs != null) {
                // check the difference
                //System.out.println("element: "+PeriodicTableData.all_elements[bs.Center_ID-1]);
                BasisFunction bs_ex = Extended_Basis_Info.TheBasis[ bs.Center_ID ][ indxs[bs.Center_ID] ];
                boolean do_copy = true;
                if (bs.L != bs_ex.L ) {
                    System.out.printf(" ERROR: different L for Z = %2d, bs # %d%n", bs.Center_ID, indxs[bs.Center_ID]);
                    do_copy = false;
                }
                if (bs.coefs.length != bs_ex.coefs.length) {
                    System.out.printf(" ERROR: different number of primitives for Z = %2d, bs # %d%n", bs.Center_ID, indxs[bs.Center_ID]);
                    do_copy = false;
                }
                if (do_copy) {
                    double diff = 0;
                    for (int c=0; c<bs.coefs.length; c++)
                        diff += Math.abs(bs.exponents[c] - bs_ex.exponents[c]);
                    //System.out.printf("expon. diff = %.10f%n",diff);
                    for (int c=0; c<bs.coefs.length; c++)
                        diff += Math.abs(bs.coefs[c] - bs_ex.coefs[c]);
                    //System.out.printf("Coef. diff = %.10f%n",diff);
                    //
                    if ((do_basis_comparison) && (diff > 1.0E-4)) {
                        System.out.printf("ERROR: coefs./exponents look very different for Z = %2d, bs # %d%n", bs.Center_ID, indxs[bs.Center_ID]);
                        do_copy = false;
                    }
                }

                if (do_copy) {
                    // important to preserve _next field!
                    bs.coefs = bs_ex.coefs;
                    bs.exponents = bs_ex.exponents;
                }
                indxs[bs.Center_ID] ++;
                bs = bs._next;
            }
        }

        // Each function from First_BS list may be used much more times than once! (for different centers!)
        // Get the total number of basis functions
        int NBasisTotal = 0;
        for (int i=0; i<molden.Centers.length; i++){
            // find all basis functions with this Z
            bs = First_BS;
            while (bs != null) {
                if (bs.Center_ID == molden.Centers[i].Z) NBasisTotal++;
                bs = bs._next;
            }
        }
        System.out.printf("Total number of radial parts of basis functions: %d%n",NBasisTotal);
        // what we have in First_BS in fact are radial parts (!) of basis functions
        // allocate memory for basis functions
        //molden.Basis = new BasisFunction[NBasisTotal]; // it is not used by MOLDEN_IO in fact...
        molden.RadialParts = new RadialPartOfBasisFunction[NBasisTotal];
        /*] OfBasisFunctions_CenterIDs = new int[NBasisTotal];
        molden.RadialPartsOfBasisFunctions_Coefs = new double[NBasisTotal][];
        molden.RadialPartsOfBasisFunctions_Exponents = new double[NBasisTotal][];
        molden.RadialPartsOfBasisFunctions_LUsedWith = new int[NBasisTotal];*/
        // fill the array in


        // Molden uses normalized basis functions, but
        // nwchem uses normalized (primitive) gaussians!
        //
        SphericalHarmonics ylm = new SphericalHarmonics();
        Polynom3D[][] YLM = ylm.Get_Quick_YLM_Norm4PI();  //  ||YLM||^2 / (4*Pi) = 1
        double[][] ylm_norms = ylm.Get_Quick_YLM_Norm2(ylm.Get_Quick_YLM());
        double[] ylm_norms_cart = ylm.molden_cart_norms2_over_4Pi();
        OverlapIntegrals oi = new OverlapIntegrals();
        bs = First_BS;
        //
        int k=0;
        //BasisFunction[] bfs = new BasisFunction[NBasisTotal]; int kk=0;
        while (bs != null) {
            // now ensure that basis functions are unity-normalized;
            // NwChem's basis set notation (as well as MOLDEN's one) assumes that
            // contraction coefficients are to be applied for unity-normalized primitive gaussians (with some angular part)
            // But for norm calculation we need to have contraction coefficients valid for UNnormalized gaussians
            // => create a copy of the basis function
            BasisFunction bs_new = new BasisFunction(bs);
            for (int j1=0; j1<bs_new.coefs.length; j1++) {
                double primitive_gaussian_norm2 = 4 * Math.PI * oi.primitive_int_1D_Sphr(2*bs_new.L+2, 2*bs_new.exponents[j1]);
                // spherical harmonics are assumed to be unity-normalized => we need to normalize radial part only!
                bs_new.coefs[j1] /= Math.sqrt( primitive_gaussian_norm2 ); 
            }
            // It may happen now that the basis function itself is not unity-normalized; so, the second step is
            // 2) make a basis function unity-normalized
            // Our SphericalHarmonics uses the same spherical harmonics as MOLDEN does, so there is no need to remormalize spherical harmonics in any way.
            bs_new.Quick_YLM = YLM; // needs to be set for calling OverlapWith
            bs_new.OI = oi;
            bs_new.R0 = new double[]{0.0, 0.0, 0.0}; // its value is not important since we're to calculate a norm of this function only
            double bs_norm2 = bs_new.OverlapWith(bs_new);

            // make true basis function unity-normalized
            for (int j1=0; j1<bs.coefs.length; j1++)
                bs.coefs[j1] /= Math.sqrt( bs_norm2 );

            bs = bs._next;
        }
        /*
        // some DEBUG:
        for (int i=0; i<kk; i++)
                System.out.printf("%13d", bfs[i].Center_ID);
        System.out.println();

        for (int i=0; i<kk; i++) {
            for (int j=0; j<kk; j++) {
                System.out.printf("%13.7f", bfs[i].OverlapWith(bfs[j]));
            }
            System.out.println();
        }
         *
         */


        //int
        k=0;

        for (int i=0; i<molden.Centers.length; i++){
            // find all basis functions with this Z
            bs = First_BS;
            //
            while (bs != null) {
                if (bs.Center_ID == molden.Centers[i].Z) {
                    //molden.Basis[k] = bs;
                    molden.RadialParts[k] = new RadialPartOfBasisFunction();
                    molden.RadialParts[k].CenterID = i+1;
                    molden.RadialParts[k].Exponents = bs.exponents;
                    molden.RadialParts[k].LUsedWith = bs.L;
                    molden.RadialParts[k].Coefs = bs.coefs;
                 /* molden.RadialPartsOfBasisFunctions_CenterIDs[k] = i+1;
                    molden.RadialPartsOfBasisFunctions_Exponents[k] = bs.exponents;
                    molden.RadialPartsOfBasisFunctions_LUsedWith[k] = bs.L;
                    molden.RadialPartsOfBasisFunctions_Coefs[k] = bs.coefs; */
                    k++;
                }
                bs = bs._next;
            }
        }
        // done!
        return true;
    }
    //--------------------------------------------------------------------------
    // reads nDoublesToRead doubles from a given file
    private double[] _read_doubles(BufferedReader ascfile, int nDoublesToRead) throws Exception{
        double[] result = new double[nDoublesToRead];
        int nReads = 0;
        String s;
        while (nReads < nDoublesToRead) {
            if ((s = ascfile.readLine()) == null) return null;
            String[] vals = s.trim().split(" +");
            for (int i=0; i<vals.length; i++) {
                result[nReads] = Double.parseDouble(vals[i]);
                nReads++;
            }
        }
        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts MO or NO information from nwchem ascii file created by mov2asc
     * @param fname
     */
    public boolean _Parse_ASC_File(String fname, MOLDEN_IO molden) throws Exception{
        FileReader fr1 = new FileReader(fname);
        BufferedReader ascfile = new BufferedReader(fr1);

        String s;

        String bs_title = molden.Title;
        System.out.println(" Reading orbitals from "+fname);
        for (int i=0; i<4; i++) ascfile.readLine(); // skip 4 lines
        System.out.println("\"scftype20\": "+ascfile.readLine());
        System.out.println("Date: "+ascfile.readLine());
        System.out.println("job type: "+ascfile.readLine());
        for (int i=0; i<3; i++) ascfile.readLine(); // skip 3 lines
        System.out.println("basis set name: "+( s = ascfile.readLine() ));
        if (!s.equals(bs_title))
            System.out.printf("WARNING: basis set name in asc file (\"%s\") and log file (\"%s\") differ%n", s, bs_title);
        int nsets = Integer.parseInt( ascfile.readLine().trim() );
        if (nsets != 1) {
            System.out.println("ERROR: asc files with more than one set of orbitals are not supported!");
            return false;
        }
        //"NBF" - number of basis functions
        int nBF = Integer.parseInt( ascfile.readLine().trim() );
        int nMO = Integer.parseInt( ascfile.readLine().trim() );
        System.out.printf(" There are %d basis functions and %d orbitals%n", nBF, nMO);
        if (nMO < nBF) {
            System.out.println("WARNING: the number of orbitals is less than the number of basis functions%n");
        }
/*
 * mov2asc.F :
      do jset = 1, Nsets
         read(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! Occupation numbers
         Write(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
C
         read(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! Eigenvalues
         Write(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
C
         do i = 1, nmo(jset)
            read(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! An eigenvector
            Write(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
         enddo
      enddo
 * Hence, there are nBF occupation numbers and nBF eigenvalues,
 * but nMO molecular orbitals...
 */
        // create MO array
        molden.MOs = new MO[nMO];

        // read in MO occupations
        double[] occupations = _read_doubles(ascfile, nBF /* see mov2asc.F */ );
        if (occupations == null) {
            System.out.println("ERROR: can not read orbital occupations.");
            return false;
        }
        for (int i=0; i<nMO; i++) {
            molden.MOs[i] = new MO();
            molden.MOs[i].Occupancy = occupations[i];
        }

        // read in MO energies
        double[] energies = _read_doubles(ascfile, nBF /* see mov2asc.F */);
        if (occupations == null) {
            System.out.println("ERROR: can not read orbital energies.");
            return false;
        }
        boolean all_positive = true;
        for (int i=0; (i<nMO) && (all_positive) ; i++) all_positive = (energies[i] > 0);
        if (all_positive) {
            System.out.println("WARNING: There are no negative values in the set of orbital energies");
            System.out.println(" (this is typical for non-SCF type calculations).");
            System.out.println(" These values will not be written to molden file!");
            for (int i=0; i<nMO; i++) molden.MOs[i].Energy = 0.0;
        } else {
            // set MO energies
            for (int i=0; i<nMO; i++) molden.MOs[i].Energy = energies[i];
        }

        // read in basis function coefs for each MO
        for (int mo=0; mo<nMO; mo++) {
            molden.MOs[mo].BS_Coefs = _read_doubles(ascfile, nBF);
        }

        double[] ylm_norms_cart = SphericalHarmonics.molden_cart_norms2_over_4Pi();

        // loop over all basis functions and swap ordering of basis function
        // coefficients if neccessary
        int bf=0;
        if (molden.IsSpherical) {
                // swaps for 'pure' functions
                for(int rp=0; rp < molden.RadialParts.length; rp++) {
                    double[] molden_fns ;
                    // no swaps for s- and p-functions
                    /*
                     if (molden.RadialPartsOfBasisFunctions_LUsedWith[rp] == 1) {
                        // change the order of p functions:
                        // MOLDEN(spherical): py     px     pz
                        // nwchem(spherical): px     py     pz
                        //             index: bf+0   +1     +2
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+1],
                                molden.MOs[mo].BS_Coefs[bf+0],
                                molden.MOs[mo].BS_Coefs[bf+2]};
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                    }
                     */
                    if (molden.RadialParts[rp].LUsedWith == 2) {
                        // change the order of d functions:
                        // MOLDEN(spherical): D 0,   D+1,   D-1,  D+2,   D-2;
                        // nwchem(spherical): D(-2), D(-1), D(0), D(+1), D(+2)
                        //             index: bf+0     +1     +2    +3     +4
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+2],
                                -molden.MOs[mo].BS_Coefs[bf+3],     // NwChem uses D(+1) with a different sign as compared to MOLDEN ! - see 'ang(n,Ad4) =-cd(1)*x*z' in xc_eval_basis.F
                                molden.MOs[mo].BS_Coefs[bf+1],
                                molden.MOs[mo].BS_Coefs[bf+4],
                                molden.MOs[mo].BS_Coefs[bf+0]};
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                    }
                    if (molden.RadialParts[rp].LUsedWith == 3) {
                        // change the order of f functions:
                        // MOLDEN(spherical): F 0,   F+1,   F-1,   F+2,  F-2,   F+3,   F-3
                        // nwchem(spherical): F(-3), F(-2), F(-1), F(0), F(+1), F(+2), F(+3)
                        //             index: bf+0     +1     +2    +3     +4      +5    +6
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+3], // F 0
                                -molden.MOs[mo].BS_Coefs[bf+4], // F+1  - NwChem uses F(+1) function with a different sign!
                                molden.MOs[mo].BS_Coefs[bf+2], // F-1
                                molden.MOs[mo].BS_Coefs[bf+5], // F+2
                                molden.MOs[mo].BS_Coefs[bf+1], // F-2
                                -molden.MOs[mo].BS_Coefs[bf+6], // F+3 - NwChem uses F(+3) function with a different sign!
                                molden.MOs[mo].BS_Coefs[bf+0]}; // F-3
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                    }
                    if (molden.RadialParts[rp].LUsedWith == 4) {
                        // change the order of g functions:
                        // MOLDEN(spherical): G 0,   G+1,   G-1,   G+2,   G-2,  G+3,   G-3,   G+4,   G-4
                        // nwchem(spherical): G(-4), G(-3), G(-2), G(-1), G(0), G(+1), G(+2), G(+3), G(+4)
                        //             index: bf+0     +1     +2    +3     +4      +5    +6      +7    +8
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+4], // G 0
                                -molden.MOs[mo].BS_Coefs[bf+5], // G+1 - NwChem uses G(+1) function with a different sign!
                                molden.MOs[mo].BS_Coefs[bf+3], // G-1
                                molden.MOs[mo].BS_Coefs[bf+6], // G+2
                                molden.MOs[mo].BS_Coefs[bf+2], // G-2
                                -molden.MOs[mo].BS_Coefs[bf+7], // G+3 - NwChem uses G(+3) function with a different sign!
                                molden.MOs[mo].BS_Coefs[bf+1], // G-3
                                molden.MOs[mo].BS_Coefs[bf+8], // G+4
                                molden.MOs[mo].BS_Coefs[bf+0]}; // G-4
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                    }

                    bf += (2*molden.RadialParts[rp].LUsedWith + 1);
                }
        } else {
            // swaps for cartesian functions
            /* Nwchem logic in cartesian function ordering:
                       for (int i=0; i<3; i++)
                            for (int j=i; j<3; j++)
                                for (int k=j; k<3; k++)
                                    System.out.println(" "+("XYZ".charAt(i))+"XYZ".charAt(j)+"XYZ".charAt(k)+", ");
             */

                for(int rp=0; rp < molden.RadialParts.length; rp++) {
                    double[] molden_fns ;
                    // no correction / order change is needed for for p- and s-functions
                    if (molden.RadialParts[rp].LUsedWith == 0) bf += 1;
                    if (molden.RadialParts[rp].LUsedWith == 1) bf += 3;
                    // but it IS for d,f,g,... ones:
                    if (molden.RadialParts[rp].LUsedWith == 2) {
                        // change the order of d functions:
                        // MOLDEN(cartesian):  xx,  yy,  zz,  xy,  xz,  yz
                        // nwchem(cartesian):  XX,  XY,  XZ,  YY,  YZ,  ZZ
                        //             index: bf+0  +1   +2   +3   +4   +5
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+0], // Math.sqrt(ylm_norms_cart[4] * 5) == 1 // xx
                                molden.MOs[mo].BS_Coefs[bf+3], // Math.sqrt(ylm_norms_cart[5] * 5) == 1 // yy
                                molden.MOs[mo].BS_Coefs[bf+5], // Math.sqrt(ylm_norms_cart[6] * 5) == 1 // zz
                                molden.MOs[mo].BS_Coefs[bf+1] * Math.sqrt(ylm_norms_cart[7] * 5),  // xy, here 5=1/(1/(2*L+1))=2*L+1
                                molden.MOs[mo].BS_Coefs[bf+2] * Math.sqrt(ylm_norms_cart[8] * 5),  // xz
                                molden.MOs[mo].BS_Coefs[bf+4] * Math.sqrt(ylm_norms_cart[9] * 5)}; // yz
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                        bf += 6;
                    }
                    if (molden.RadialParts[rp].LUsedWith == 3) {
                        // change the order of f functions:
                        // MOLDEN(cartesian): xxx,  yyy,  zzz,  xyy,  xxy,  xxz,  xzz,  yzz,  yyz,  xyz
                        // nwchem(cartesian): XXX,  XXY,  XXZ,  XYY,  XYZ,  XZZ,  YYY,  YYZ,  YZZ,  ZZZ
                        //             index: bf+0   +1    +2    +3    +4    +5    +6    +7    +8    +9
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+0],  // xxx
                                molden.MOs[mo].BS_Coefs[bf+6],  // yyy
                                molden.MOs[mo].BS_Coefs[bf+9],  // zzz
                                molden.MOs[mo].BS_Coefs[bf+3],  // xyy
                                molden.MOs[mo].BS_Coefs[bf+1],  // xxy
                                molden.MOs[mo].BS_Coefs[bf+2],  // xxz
                                molden.MOs[mo].BS_Coefs[bf+5],  // xzz
                                molden.MOs[mo].BS_Coefs[bf+8],  // yzz
                                molden.MOs[mo].BS_Coefs[bf+7],  // yyz
                                molden.MOs[mo].BS_Coefs[bf+4]}; // xyz
                            // some correction is still needed since different cartesian functions have different norms!
                            for (int k=0; k<molden_fns.length; ++k) molden_fns[k] *= Math.sqrt(ylm_norms_cart[10+k] * 7); // 7 = 1/(1/(2*L=1))
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                        bf += 10;
                    }
                    if (molden.RadialParts[rp].LUsedWith == 4) {
                        // change the order of g functions:
                        // MOLDEN(cartesian): xxxx  yyyy  zzzz  xxxy  xxxz  yyyx  yyyz  zzzx  zzzy  xxyy  xxzz  yyzz  xxyz  yyxz  zzxy
                        // nwchem(cartesian): XXXX, XXXY, XXXZ, XXYY, XXYZ, XXZZ, XYYY, XYYZ, XYZZ, XZZZ, YYYY, YYYZ, YYZZ, YZZZ, ZZZZ
                        //             index: bf+0    +1    +2    +3    +4    +5    +6    +7    +8    +9   +10   +11   +12   +13   +14
                        for (int mo=0; mo<nMO; mo++) {
                            molden_fns = new double[]{
                                molden.MOs[mo].BS_Coefs[bf+ 0],  // xxxx
                                molden.MOs[mo].BS_Coefs[bf+10],  // yyyy
                                molden.MOs[mo].BS_Coefs[bf+14],  // zzzz
                                molden.MOs[mo].BS_Coefs[bf+ 1],  // xxxy
                                molden.MOs[mo].BS_Coefs[bf+ 2],  // xxxz
                                molden.MOs[mo].BS_Coefs[bf+ 6],  // yyyx
                                molden.MOs[mo].BS_Coefs[bf+11],  // yyyz
                                molden.MOs[mo].BS_Coefs[bf+ 9],  // zzzx
                                molden.MOs[mo].BS_Coefs[bf+13],  // zzzy
                                molden.MOs[mo].BS_Coefs[bf+ 3],  // xxyy
                                molden.MOs[mo].BS_Coefs[bf+ 5],  // xxzz
                                molden.MOs[mo].BS_Coefs[bf+12],  // yyzz
                                molden.MOs[mo].BS_Coefs[bf+ 4],  // xxyz
                                molden.MOs[mo].BS_Coefs[bf+ 7],  // yyxz
                                molden.MOs[mo].BS_Coefs[bf+ 8]}; // zzxy
                            for (int k=0; k<molden_fns.length; ++k) molden_fns[k] *= Math.sqrt(ylm_norms_cart[20+k] * 9); // 9 = 1/(1/(2*L=1))
                            System.arraycopy(molden_fns, 0, molden.MOs[mo].BS_Coefs, bf, molden_fns.length);
                        }
                        bf += 15;
                    }
                }
        }

        //
        ascfile.close();
        fr1.close();

        return true;

    }

}
