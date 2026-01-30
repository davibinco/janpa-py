package molden2molden;

/*
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
 */

import moldenio.*;
import Polynom3D.*;
import JGints.*;
import java.io.*;
import ProgramOptions.WarningManager;

/**
 * A class for performing basis set / MO coefs. transformation within molden files
 * Version: 21.Jan.2014
 * @author timn
 */
public class MoldenConvertor {
    //--------------------------------------------------------------------------
    PrintStream out = System.out;
    //--------------------------------------------------------------------------
    /**
     * Searches for the first basis function with L==patternL, m==patternM, additional_r_power == pattern_AdRPow,
     * being centered at the same center as pattern, and having
     * the same exponents and coefficients (last two - within +/- eps) as the
     * @pattern has. The search is being performed the list of radial parts (!) beginning with @first_bf .
     * @Returns the index of the m=0-th ('pure'!) component of the basis function found,
     * or -1 if none was found
     *
     * Rev.: 15.Jan.2014
     */
    private static int _Find_BS_with_Props(BasisFunction first_bf, BasisFunction pattern, int patternL, int patternM,
                                           int pattern_AdRPow, double eps) {
        int i=0;
        boolean found = false;
        BasisFunction bf = first_bf;
        while ((!found) && (bf != null)) {
            boolean all_equal = false;
            if ((bf.L == patternL) &&
                /*(bf.m == patternM) &&*/ // this should NOT be checked!!!
                (bf.Center_ID == pattern.Center_ID) && // very important!!!
                (bf.coefs.length == pattern.coefs.length) &&
                (bf.additional_r_power == pattern_AdRPow) &&
                (bf.exponents.length == pattern.exponents.length)) {
                    all_equal = true; // an assumption
                    for(int cf=0; cf<bf.coefs.length; ++cf) {
                        all_equal = all_equal & ( Math.abs(bf.coefs[cf]-pattern.coefs[cf]) <= eps );
                        all_equal = all_equal & ( Math.abs(bf.exponents[cf]-pattern.exponents[cf]) <= eps );
                    }
            }
            found = all_equal;
            if (!found) {
                i += (2*bf.L + 1);
                bf = bf._next;
            } else {
                // adjust i according to MOLDEN m-component ordering
                switch (patternM) {
                    case 0:  return i;
                    case +1: return (i+1);
                    case -1: return (i+2);
                    case +2: return (i+3);
                    case -2: return (i+4);
                    case +3: return (i+5);
                    case -3: return (i+6);
                    case +4: return (i+7);
                    case -4: return (i+8);
                    // *** And so on for higher Ls ***
                }
            }
        }
        return -1; // none found!
    }
    //--------------------------------------------------------------------------
        //public String fname_in = null;
        //public String fname_out = null;

        MOLDEN_IO molden_in = new MOLDEN_IO();
        MOLDEN_IO molden_out = null;

        m2mOptions options = new m2mOptions();

    //--------------------------------------------------------------------------
    // check for non-controversy of the commands:
    public boolean CommandsCheck() {
        if ((options.op_from_orca3_BF.get_boolean() && options.op_from_psi4b4_BF.get_boolean()))
                WarningManager.warning_printf(" WARNING: using and -FromOrca3BF and -FromPsi4b4BF keys together might give senseless results!");
        if ((options.op_to_orca3_BF.get_boolean() && options.op_to_psi4b4_BF.get_boolean()))
                WarningManager.warning_printf(" WARNING: using and -ToOrca3BF and -ToPsi4b4BF keys together might give senseless results!");

        // check whether we have in- and out- filenames
        if (options.fname_in.get_String().isEmpty() && options.fname_in_BINAR.get_String().isEmpty()) {
            out.println("Cann't continue without an input file name!");
            return false;
        }
        if (options.fname_out.get_String().isEmpty() && options.fname_out_BINAR.get_String().isEmpty()) {
            out.println("Cann't continue without an output file name!");
            return false;
        }
        return true;            
    }
    //--------------------------------------------------------------------------
    // does the main job
    public boolean CONVERT() throws Exception {
        
        boolean load_ok = false;
        if ( !options.fname_in_BINAR.get_String().isEmpty() ) {
            out.println(" Loading input molden from "+options.fname_in_BINAR.get_String()+" ...");
            load_ok = molden_in.Load_from_BMOLDEN( options.fname_in_BINAR.get_String() );
        } else {
            out.println(" Loading input molden from "+options.fname_in.get_String()+" ...");
            load_ok = molden_in.Load_From_MOLDEN(options.fname_in.get_String());
        }
               
        if (load_ok)
            out.printf("Data loaded successfully.%n%n");
        else {
            out.println(" ERROR: Can not load data from "+options.fname_in.get_String());
            return false;
        }

        if (options.op_cart2spher.get_boolean() && molden_in.IsSpherical) {
            out.println(" ERROR: input molden file already uses pure spherical harmonics,");
            out.println("   -cart2spher key makes no sense and will be ignored!");
            options.op_cart2spher.set_data( false );
        }

        //
        //=================================================================================================================

        // Convert coordinates if neccessary
        if ((!options.op_coordsAngstroms.get_boolean()) && (!molden_in.Coords_in_AU)) {
            out.println("Converting atomic coordinates from Angstroms to a.u. ...");
            molden_in.CoordsToAU(); // default is to convert coords. to a.u.
        }
        if ((options.op_coordsAngstroms.get_boolean()) && (molden_in.Coords_in_AU)) {
            out.println("Converting atomic coordinates from a.u. to Angstroms ...");
            molden_in.CoordsToAngstroms(); // convert to Angstroms if neccessary
        }
        //=================================================================================================================

        // SOME PERPARATIONS FOR BASIS FUNCTION MANIPULATIONS

        // STEP 0: check whether the highest angular momentum is not greater than G (L=4)
        Polynom3D[][] ylm = SphericalHarmonics.Get_Quick_YLM();
        boolean L_to_high = (molden_in.HighestL > SphericalHarmonics.L_MAX);
        if (L_to_high)
            WarningManager.warning_printf(" WARNING: input molden file contains %c-functions, which are not fully supported now! %n",
                     "spdfghijklmnopqr".charAt(molden_in.HighestL));

        // Check whether ORCA input file uses purely spherical harmonics
        if ((options.op_orca3signs.get_boolean() || options.op_from_orca3_BF.get_boolean() ||
                options.op_to_orca3_BF.get_boolean()) && (!molden_in.IsSpherical)) {
            // ORCA file should be in the pure spherical basis
            out.println(" ERROR: Orca is said to use pure spherical functions only, but the input molden file does not!");
            return false;
        }
        //=================================================================================================================

        // MAIN CONVERSIONS START HERE

        //=================================================================================================================
        // WORKMODE1: BASIS SET / MO SIGNS CONVERSION: PACKAGES -> MOLDEN
        
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // a) supress MO data with the one from .47 file
        //    it is better to do it BEFORE MO coefs sign change since ORCA uses its own basis funtion for writing .47 file...
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // a.1) Are we supposed to replace the GEOMETRY with the one from the .47 file ?
        if (!options.File_GEOM47.get_String().isEmpty()) {
            String fname47 = options.File_GEOM47.get_String();
            out.println("Supressing geometry with the data from .47 file "+fname47);
            File47Import a = new File47Import();
            a.import_from_file(fname47);
            if (!a.Insert_Geom_Into( molden_in )) {
                out.println("Error importing basis set data!");
                return false;
            }
        }

        // a.2) Are we supposed to replace the BASIS SET with the one from the .47 file ?
        if (!options.File_BS47.get_String().isEmpty()) {
            String fname47 = options.File_BS47.get_String();
            out.println("Supressing basis set with the data from .47 file "+fname47);
            File47Import a = new File47Import();
            a.import_from_file(fname47);
            if (!a.Build_Basis_Set( molden_in )) {
                out.println("Error importing basis set data!");
                return false;
            }
        }

        // a.3) Are we supposed to replace MOLECULAR ORBITALS with the ones from the .47 file ?
        if (!options.File_DS47.get_String().isEmpty()) {
            String fname47 = options.File_DS47.get_String();
            out.println("Supressing MO coefficients with the data from .47 file "+fname47);
            File47Import a = new File47Import();
            a.import_from_file(fname47);
            if (a.nbas != molden_in.Basis.length) {
                out.printf("ERROR: the number of basis functions is different in %n"+
                        "the input molden file and in the given .47 file");
                return false;
            }
            if (!a.Check_BF_order_and_create_remap(molden_in))  {
                out.println("ERROR: reordering of the .47 file basis set failed!");
                return false;
            }
            a.Procude_NOs( molden_in );
        }
        
        // Note[03.Sep.2016]: importing data from .47 file moved before basis set and MO manipulations!
        
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-        
        // b) basis set
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        
        double Rnorm2 = 0;
        double[] orca_dividers = new double[]{1/*S*/, 1/*P*/, 3/*D*/, 15/*F*/, 35/*G*/}; // additional denominators for 4*Pi/(2*L+1)

        // Basis set contraction coefficient corrections should be done within RadialParts... arrays,
        // but not with BasisFunctions, since the former are used in creating new MOLDEN file by Save_As_MOLDEN()
        // Note that molden_in.Basis should also be affected by this conversion since they might
        // be further used as a source of information for cartesian->pure conversion
        
        // Execult this for inputs from: ORCA3, PSI4 beta4/beta5; but NOT PSI4 1.0
        if (options.op_from_orca3_BF.get_boolean() || options.op_from_psi4b4_BF.get_boolean()) {
            if (options.op_from_orca3_BF.get_boolean())
                out.println("Converting the basis set from the ORCA 3-style to a conventional MOLDEN style...");
            if (options.op_from_psi4b4_BF.get_boolean())
                out.println("Converting the basis set from the PSI4.b4-style to a conventional MOLDEN style...");

            // loop over all radial parts and re-calculate radial parts contractions coefs by making them
            // suitable for SUM_i{ C[i]*NormalizedExp(...) } rather than for SUM_i{ C[i]*exp(...) } -
            // a conventional MOLDEN file stores C[i] in the former form, while Orca3/Psi4(beta4/beta5) did it in the latter form.
            // However, Psi4 1.0 release stores C[i] in a conventional fashion!
            for (int rp=0; rp<molden_in.RadialParts.length; ++rp) {
                int L = molden_in.RadialParts[rp].LUsedWith;
                // loop over primitive gaussians in this basis function radial part
                for (int cf=0; cf<molden_in.RadialParts[rp].Coefs.length; ++cf) {
                    //
                    //  Rnorm2 = INT(  {exp(-zeta*r^2) * r^L * r^addit.}^2 * r^2, r=0..+infinity  )
                    //
                    Rnorm2 = OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2 + 2*molden_in.RadialParts[rp].Addit_r_power ,
                            2 * molden_in.RadialParts[rp].Exponents[cf]);
                    Rnorm2 *= 4*Math.PI;
                    Rnorm2 /= (2*L + 1);
                    // now
                    //   Rnorm2 = 4*Pi / (2*L+1) * INT( {exp(-zeta*r^2) * r^L * r^addit.}^2 * r^2, r=0..+infinity  )
                    // Additional divider if ORCA file is used:
                    if (options.op_from_orca3_BF.get_boolean())
                        Rnorm2 /= orca_dividers[ L ];
                    // Now Rnorm2 contains the proper norm2 of the basis function
                    // Finally, correct the contraction coefficient:
                    molden_in.RadialParts[rp].Coefs[cf] *= Math.sqrt(Rnorm2);
                }
            }
            // Now we should copy updated contraction coefs. from radial parts back to all the basis functions
            for (int bf=0; bf<molden_in.Basis.length; ++bf)
                molden_in.Basis[bf].coefs = molden_in.RadialParts[ molden_in.Basis[bf].RadialPart_ID ].Coefs.clone();
        }
        
        // Conversion for MO coefs.
        // Execult this for: PSI4 beta4/beta5/1.0 (cartesian only!), but NOT ORCA
        if (options.op_from_psi4b4_BF.get_boolean() || options.op_from_psi4v1_MO.get_boolean()) {            
            if (!molden_in.IsSpherical) {
                double[] cart_norm2s = SphericalHarmonics.molden_cart_norms2_over_4Pi();
                double[] scale_factors = new double[molden_in.Basis.length]; // scale factors for individual MO coefs.
                for (int i=0; i<scale_factors.length; ++i) {
                    scale_factors[i] = Math.sqrt( cart_norm2s[molden_in.Basis[i].m] ) * Math.sqrt(2*molden_in.Basis[i].L+1);
                }
                // scale MO coefs
                for (int mo=0; mo<molden_in.MOs.length; ++mo)
                    for (int i=0; i<scale_factors.length; ++i)
                        molden_in.MOs[mo].BS_Coefs[i] *= scale_factors[i];
            }

        }


        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // c) MO coefs. conversion from packages to conventional MOLDEN
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        if (options.op_orca3signs.get_boolean()) {
            out.println("Converting MO coefficients signs to make them compatible with a conventional MOLDEN spherical harmonics...");

            // Orca uses different signs for:
            // F(+3), F(-3), G(-3), G(+3), G(-4), G(+4)
            // get the total number of these basis functions
            int nSingChanges = 0;
            for (int bf=0; bf<molden_in.Basis.length; ++bf) {
                if ( molden_in.Basis[bf].L == 3) nSingChanges += 2; // two F functions should have their sign changed
                if ( molden_in.Basis[bf].L == 4) nSingChanges += 4; // four G functions should have their sign changed
            }
            // there are two possibilities: to change the sign of contraction coefs or to change
            // the signs of MO coefs. We follow the second way.
            int[] sings_inverse = new int[nSingChanges]; // numbers of F(+3), F(-3), G(-3), G(+3), G(-4), G(+4) components of basis functions
            //int component_index = 0; // the index of the 1-st MO coef of the basis function
            int i = 0;
            for (int bf=0; bf<molden_in.Basis.length; ++bf) {
                // MOLDEN ordering of components of F and G functions is: F 0, F+1, F-1, F+2, F-2, F+3, F-3 and
                //  G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
                if (( molden_in.Basis[bf].L == 3) && (molden_in.Basis[bf].m == -3))     // F(-3)
                    sings_inverse[i++] = bf;
                if (( molden_in.Basis[bf].L == 3) && (molden_in.Basis[bf].m ==  3))     // F(+3)
                    sings_inverse[i++] = bf;

                if (( molden_in.Basis[bf].L == 4) && (molden_in.Basis[bf].m == -4))     // G(-4)
                    sings_inverse[i++] = bf;
                if (( molden_in.Basis[bf].L == 4) && (molden_in.Basis[bf].m ==  4))     // G(+4)
                    sings_inverse[i++] = bf;
                if (( molden_in.Basis[bf].L == 4) && (molden_in.Basis[bf].m == -3))     // G(-3)
                    sings_inverse[i++] = bf;
                if (( molden_in.Basis[bf].L == 4) && (molden_in.Basis[bf].m ==  3))     // G(+3)
                    sings_inverse[i++] = bf;

            }

            // now change the signs of MO coefs
            for (int mo=0; mo<molden_in.MOs.length; ++mo)
                for (/*int*/ i=0; i<sings_inverse.length; ++i)
                    molden_in.MOs[mo].BS_Coefs[ sings_inverse[i] ] = -molden_in.MOs[mo].BS_Coefs[ sings_inverse[i] ];
        }
        //=================================================================================================================

        // OK, now mondel_in is guaranteed to be in a conventional MOLDEN format

        // WORKMODE2
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // a) NORMALIZATION OF BFs
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        if (options.op_normalizeBF.get_boolean()) {
            // work with radial parts only!
            out.println("Forcing basis functions to be unity-normalized...");

            for (int rp=0; rp<molden_in.RadialParts.length; ++rp) {
                int L = molden_in.RadialParts[rp].LUsedWith;

                int nPrimitiveGaussians = molden_in.RadialParts[rp].Coefs.length;
                double[] gnorms2 = new double[nPrimitiveGaussians];

                // Fill gnorms2 with
                // gnorms2[i] = 4*Pi/(2*L+1) * INT(  {exp(-zeta_i*r^2) * r^L}^2 * r^2, r=0..+infinity  )
                for (int i=0; i<nPrimitiveGaussians; ++i)
                    gnorms2[i] = /* 4.0*Math.PI/(2*L+1) * */
                                OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2 + 2*molden_in.RadialParts[rp].Addit_r_power,
                                    2 * molden_in.RadialParts[rp].Exponents[i] );

                // Get the squared norm of the basis function
                double norm2 = 0;
                // Loop over primitive gaussians in this basis function radial part and compute the sum
                //
                //  norm2 = SUM_i_j( c_i*c_j/SQRT(gnorms2[i]*gnorms2[j]) * INT(  {exp(-zeta_i*r^2) * r^L}*{exp(-zeta_j*r^2) * r^L} * r^2, r=0..+infinity  ) )
                //
                for (int i=0; i<nPrimitiveGaussians; ++i) {
                    norm2 += molden_in.RadialParts[rp].Coefs[i] *
                            molden_in.RadialParts[rp].Coefs[i] / gnorms2[i] *
                            OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2 + 2*molden_in.RadialParts[rp].Addit_r_power,
                                2 * molden_in.RadialParts[rp].Exponents[i] );
                    for (int j=(i+1); j<nPrimitiveGaussians; ++j)
                        norm2 += 2 * molden_in.RadialParts[rp].Coefs[i] / Math.sqrt(gnorms2[i]) *
                                    molden_in.RadialParts[rp].Coefs[j] / Math.sqrt(gnorms2[j])*
                                    OverlapIntegrals.primitive_int_1D_Sphr( 2* L + 2 + 2*molden_in.RadialParts[rp].Addit_r_power,
                                            molden_in.RadialParts[rp].Exponents[i] +
                                            molden_in.RadialParts[rp].Exponents[j]); // each integral is being encountered twice
                }
                // Muptiply by the squared norm of a spherical harmonic
                /* norm2 *= (4.0*Math.PI / (2*L + 1) ); */  // this can be skipped as soon as 4.0*Math.PI/(2*L+1) is absent in gnorms2
                //out.printf("radial part: %3d, norm2 = %17.15f %n", rp+1, norm2);
                // Finally, correct the contraction coefficient:
                for (int cf=0; cf<molden_in.RadialParts[rp].Coefs.length; ++cf)
                    molden_in.RadialParts[rp].Coefs[cf] /= Math.sqrt(norm2);
            }
        }


        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        // b) cartesian <-> spherical interconversion of basis MO coefs and basis functions
        // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
        if (options.op_cart2spher.get_boolean()) {

            molden_out = new MOLDEN_IO();
            molden_out.Allow_additional_r_power = true;
            molden_out.Centers = molden_in.Centers;
            molden_out.Coords_in_AU = molden_in.Coords_in_AU;
            molden_out.Title = " created by molden2molden from "+options.fname_in.get_String();

            // Get a transformation matrix
            YLM_Series[] c2p = SphericalHarmonics.Cartesian_to_Pure();
            // these expansion coefs. are suitable if pure functions defined in SphericalHarmonics.Get_Quick_YLM() are used.
            // We should use pure harmonics with a slightly different norms, however, to be consistent with MOLDEN
            // Therefore, coefs in c2p need some corrections:
            //                                 SQRT(||Y_lm_ours||^2)                                         Y_lm_ours                  1
            //  c_ours * Y_lm_ours = c_ours * ---------------------- * Y_lm_ours = c_molden * SQRT(2*L+1) * ----------------------- * --------------
            //                                 SQRT(||Y_lm_ours||^2)                                         SQRT(||Y_lm_ours||^2)    SQRT( 2*L+1 )
            // So, c_molden  = c_ours  * SQRT(||Y_lm_ours||^2) / SQRT( ||x^nx*y^ny*z^nz||^2 )
            //
            double[] cart_norm2s = SphericalHarmonics.molden_cart_norms2_over_4Pi();
            double[][] ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(SphericalHarmonics.Get_Quick_YLM());
            //BasisFunction[][] X = new BasisFunction[c2p.length][10];
            for (int m=0; m<c2p.length; ++m)
                for (int cf=0; cf<c2p[m].coefs.length; ++cf) {
                    c2p[m].coefs[cf] *= Math.sqrt(ylm_norms2[c2p[m].Ls[cf]][c2p[m].Ls[cf] + c2p[m].Ms[cf]]);
                    c2p[m].coefs[cf] /= Math.sqrt(cart_norm2s[m]);
                }

            int[][] where_to_store; // [b][i]-th index element contains the index of pure function
                         // needed by the (c2p[basis[b].m].Ls[i] && c2p[basis[b].m].Ms[i])-th component
                         // of the expansion of b-th cartesian basis function of molden_in
            double[][] transfer_coef; // a factor for the MO coef to be multiplied by when transforming
                         // from cartesian to pure basis set
            // initialize a where_to_store[][] array
            where_to_store = new int[molden_in.Basis.length][];
            transfer_coef = new double[molden_in.Basis.length][];

            for (int b=0; b<molden_in.Basis.length; ++b) {
                where_to_store[b] = new int[ c2p[ molden_in.Basis[b].m ].coefs.length ];
                transfer_coef[b] = new double[ c2p[ molden_in.Basis[b].m ].coefs.length ];
            }

            // Get a total number of pure functions needed to represent all the cartesian ones
            // At this stage, bs_first and bf are being used as a list for storing RADIAL PARTS of
            // the pure basis functions
            BasisFunction bf_first = null;
            BasisFunction bf_last = null;
            out.println("Creating a list of pure spherical harmonics...");
            int nIgnorable = 0; // number of basis functions to ignore (if IgnoreLowerLTerms==true)
            // loop over all basis functions of the input file
            for (int bf=0; bf<molden_in.Basis.length; ++bf) {
                // loop over all pure harmonics needed to represent bf-th cartesian function
                for (int i=0; i<c2p[ molden_in.Basis[bf].m ].coefs.length; ++i) {
                    // try finding this a basis function with the given L, primitive exponents and contraction coefs.
                    int L = c2p[molden_in.Basis[bf].m ].Ls[i]; // desired L
                    int addit_r_pwr = (molden_in.Basis[bf].L - L) + molden_in.Basis[bf].additional_r_power; // the power of r multiplier required by the expansion
                    int indx = _Find_BS_with_Props(bf_first,
                                molden_in.Basis[bf], // desired primitive gaussian exponents and contraction coefs.
                                L, // desired L
                                c2p[molden_in.Basis[bf].m ].Ms[i], // desired M
                                addit_r_pwr,
                                0.0);
                    boolean _append = (indx == -1); // whether to create new basis function
                    if ((addit_r_pwr > 0) && options.IgnoreLowerLTerms.get_boolean()) {
                        _append = false;
                        ++nIgnorable;
                    }

                    if (_append) {
                        // append (a list of) the new function with desired L, and m=0 only(!)
                        if (bf_first == null)
                            bf_last = (bf_first = new BasisFunction(molden_in.Basis[bf]));
                        else
                            bf_last = (bf_last._next = new BasisFunction(molden_in.Basis[bf]));
                        bf_last._next = null; // since molden_in.Basis[bf] could have it set to non-null
                        bf_last.L = L; // L of i-th term of expansion CAN be different from molden_in.Basis[bf].L
                        // The value of m of the newly created function will be ignored anyway!
                        bf_last.additional_r_power = addit_r_pwr;

                        // Repeat the search to get proper indx value (its not very 'clever',
                        // but should be easily understandable; in particular, this should give proper m-dependent value
                        indx = _Find_BS_with_Props(bf_first,
                                    molden_in.Basis[bf], // desired primitive gaussian exponents and contraction coefs.
                                    L, // desired L
                                    c2p[molden_in.Basis[bf].m ].Ms[i], // desired M
                                    addit_r_pwr,
                                    0.0);
                        if (indx == -1) {
                            // in fact, this is really an impossible error
                            out.println("ERROR in creating basis function list: new function not found!");
                            return false;
                        }
                    }
                    // Now store the index in where_to_store array
                    where_to_store[bf][i] = indx;
                    transfer_coef[bf][i] = c2p[molden_in.Basis[bf].m ].coefs[i];
                }
            }

            // get the number of radial parts and the number of basis functions
            int nRadParts = 0;
            int nBasisFns = 0;
            BasisFunction tmp = bf_first;
            if (tmp == null) {
                out.println("OOOPS! bf_first == null: no basis functions have been created...");
                return false;
            }
            while (tmp != null) {
                //out.printf("BS: %4d, CNTR %2d, %c, 1st-exp: %7.5f %n", nRadParts, tmp.Center_ID, "spdfg".charAt(tmp.L), tmp.exponents[0]);
                ++nRadParts;
                nBasisFns += (2*tmp.L + 1);
                tmp = tmp._next;
            }
            out.printf("New number of basis functions: %4d, new number of radial parts: %4d%n",nBasisFns, nRadParts);
            // fill in molden_out.*RadialParts* arrays with the information from the bf_first list
            // (which actually contains the radial parts stored in the BasisFunction objects for the sake of convenience!)
            // Allocate mem
            molden_out.RadialParts = new RadialPartOfBasisFunction[nRadParts];
         /* molden_out.RadialPartsOfBasisFunctions_CenterIDs = new int[nRadParts];
            molden_out.RadialPartsOfBasisFunctions_Coefs = new double[nRadParts][];
            molden_out.RadialPartsOfBasisFunctions_Exponents = new double[nRadParts][];
            molden_out.RadialPartsOfBasisFunctions_LUsedWith = new int[nRadParts];
            molden_out.RadialPartsOfBasisFunctions_Addit_r_power = new int[nRadParts]; */
            // fill arrays in
            tmp = bf_first;
            for (int i=0; i<nRadParts; ++i) {
                molden_out.RadialParts[i] = new RadialPartOfBasisFunction();
                molden_out.RadialParts[i].CenterID = tmp.Center_ID;
                molden_out.RadialParts[i].LUsedWith = tmp.L;
                molden_out.RadialParts[i].Coefs = tmp.coefs.clone();
                molden_out.RadialParts[i].Exponents = tmp.exponents.clone();
                molden_out.RadialParts[i].Addit_r_power = tmp.additional_r_power;
                tmp = tmp._next;
            }
            // Basis functions can now be creating by proper cloning radial parts
            molden_out.IsSpherical = true;

            out.println(" Transforming MO coefficients...");
            int target_coef;
            double adduct;
            // Now create MOs arrays and fill them in
            molden_out.MOs = new MO[molden_in.MOs.length];
            for (int mo=0; mo<molden_out.MOs.length; ++mo) {
                molden_out.MOs[mo] = new MO( molden_in.MOs[mo] ); // copy energy, occupancy, ...
                molden_out.MOs[mo].BS_Coefs = null; // distroy the reference to original array
                // and alocate a new one, with a different number of basis functions coefs.
                molden_out.MOs[mo].BS_Coefs = new double[nBasisFns];
                double sumIgnored = 0;
                // Calculate new MO coefs.
                // loop over all coefs of the mo-th orbital in cartesian representation
                for (int cf=0; cf<molden_in.Basis.length; ++cf)
                    // loop over all pure conponents of cf-th cartesian function
                    for (int i=0; i<where_to_store[cf].length; ++i) {
                        target_coef = where_to_store[cf][i];
                        adduct = transfer_coef[cf][i] * molden_in.MOs[mo].BS_Coefs[cf];
                        if (target_coef != (-1))
                            molden_out.MOs[mo].BS_Coefs[ target_coef ] += adduct;
                        else                            
                            if (options.IgnoreLowerLTerms.get_boolean()) sumIgnored += adduct; // if target_coef==(-1)
                    }
                    if (Math.abs(sumIgnored) > 1.0E-5)
                        WarningManager.warning_printf("Warning: sum of ignored terms (~r^n with n>0) for MO %3d  %n"+
                                "   seems to be non-negligible (|sum| = %.3E)%n",
                                mo+1, sumIgnored);

            }
            out.println(" Cartesian to pure transformation finished!");

            molden_in = molden_out; // there is no need in saving an 'old' file

        } else {
            molden_out = molden_in; // if no cartesian->pure conversion was required
        }

        //=================================================================================================================
        // WORKMODE3: BASIS SET CONVERSION: MOLDEN -> PACKAGES
        // MO coefs remain unchanged!
        if (options.op_to_orca3_BF.get_boolean() || options.op_to_psi4b4_BF.get_boolean()) {
            if (options.op_to_orca3_BF.get_boolean()) {
                out.println("Converting the basis set from a conventional MOLDEN style to the ORCA 3-style...");
                if (!molden_in.IsSpherical)
                    WarningManager.warning_printf(" WARNING: ORCA 3-style molden file with cartesian basis%n"+
                                      " functions seems to be useless...");
            }
            if (options.op_to_psi4b4_BF.get_boolean()) {
                out.println("Converting the basis set from a conventional MOLDEN style to the PSI4.b4-style...");
                if (!molden_in.IsSpherical)  {
                    out.println("Sorry, not implemented yet!"); // since a conversion for MO coefs is required
                    return false;
                }
            }

            // the simplest way is to use molden_in.UnNormalizePrimitives() in a 'reverse' fasion
            // Create an array of PSI4-compatible norms of spherical functions:
            double[][] ylm_norm2s = new double[ molden_in.HighestL+1 ][];
            for (int L=0; L<=molden_in.HighestL; ++L) {
                ylm_norm2s[L] = new double[ 2*L + 1 ];
                double nrm2 = 1.0 / (2*L + 1); // without 4*Pi !
                // Additional divider if ORCA file is used:
                if (options.op_to_orca3_BF.get_boolean())
                    nrm2 /= orca_dividers[ L ];
                // Set these values
                for (int m=-L; m<=L; ++m)
                    ylm_norm2s[L][L+m] = nrm2;  // 'reverse' fashion
            }
            // now convert the basis set
            molden_in.UnNormalizePrimitives(ylm_norm2s);
            //throw new Exception(" check UnnormalizePrimitives for 2*2*Additional power !");
            molden_out = molden_in;
        }
        //=================================================================================================================


        // SAVE RESULT
        //molden_out = molden_in;

        out.print("Saving result as "+options.fname_out.get_String()+" ");

        boolean save_result = false;

        if (!options.fname_out_BINAR.get_String().isEmpty()) {
            out.print(" (binary encoded)...");
            save_result = molden_out.Save_as_BMOLDEN(options.fname_out_BINAR.get_String());
        } else {
            out.print(" (ascii encoded)...");
            save_result = molden_out.Save_As_MOLDEN(options.fname_out.get_String());
        }

        if (save_result)
            out.println("OK.");
        else {
            out.println("ERROR saving file.");
            return false;
        }

        return true;
    }
    //--------------------------------------------------------------------------


}
