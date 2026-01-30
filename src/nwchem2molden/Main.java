package nwchem2molden;

/**
 * An application to convert nwchem ascii files produced by mov2asc into molden-compatible format
 * Copyright (c) Tymofii Nikolaienko, 2013
 * Version 05.05.2014
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


import java.io.*;
import moldenio.*;
import java.util.*;
import Polynom3D.*;


//==============================================================================
/**
 * Main class
 * @author Tymofii Nikolaienko
 */
public class Main {
    //--------------------------------------------------------------------------
    final static String Version = "13-01-2019"; // 05.May.2014
    //--------------------------------------------------------------------------
    static nw2molden_options Options = new nw2molden_options();
    //--------------------------------------------------------------------------
    /** An entry point; args is the command line arguments */
    public static void main(String[] args) throws Exception {

        Locale.setDefault(Locale.US);
        System.out.println(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        System.out.println(" * nwchem2molden: produces MOLDEN-compatible files from NwChem results * ");
        System.out.println(" *        A part of JANPA package,   http://janpa.sourceforge.net      * ");
        System.out.printf (" *                      VERSION: %14s                        * %n",Version);
        System.out.println(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        System.out.println();
        System.out.println("         (c) Tymofii Nikolaienko, 2014");
        System.out.println();
        System.out.println(" If any results obtained with this program are published,");
        System.out.println(" or for any other reasons, please, cite this work as: ");
        System.out.println("  T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem.(2014),");
        System.out.println("  V.1050, P.15-22, DOI: 10.1016/j.comptc.2014.10.002");
        System.out.println("                       *  *  *                                   ");
        System.out.println();

        Options = new nw2molden_options();
        
        Options._FirstOption.LoadOptionsFromCommandLine(args, null, true, true);
        
        // are there all neccessary parameters available?
        if (Options.LogFile.get_String().isEmpty() || Options.MO_File.get_String().isEmpty() ||
        Options.Molden_File.get_String().isEmpty()) {
            // not enough parameters => print usage information
            System.out.println("Usage:");
            System.out.println("java -jar nwchem2molden.jar -option1 value1 -option2 value2 ...");
            System.out.println();
            System.out.println("where possible options are");
            Options._FirstOption.PrintDescriptions(null, " %s: %s%n");
            System.out.println("                       *  *  *                               ");
            return;
        }

        // create a class for actual conversion jobs (a 'workhorse')
        nw2moldenConverter workhorse = new nw2moldenConverter();
        workhorse.molden.MOCoefLineFormat = "%3d %22.15f%n"; // 15 decimal digits are available by nwchem


        // is there an additional basis set data available ?
        if (!Options.BasisSet_File.get_String().isEmpty()) {
            System.out.println("Loading external basis set information from " + Options.BasisSet_File.get_String());
            workhorse.Extended_Basis_Info = new BasisSetLoader();
            workhorse.Extended_Basis_Info.print_loaded_basis = Options.PrintNwBasis.get_boolean();
            workhorse.Extended_Basis_Info.LoadFromNwFile(Options.BasisSet_File.get_String(), false);
            workhorse.do_basis_comparison = !Options.CheckNwBasis.get_boolean();
        }


        
        // 1) get the geometry (and basis set) data from the log file
        System.out.println("Analyzing the log file "+Options.LogFile.get_String());
        if (workhorse._Parse_OUT_File(Options.LogFile.get_String(), workhorse.molden))
            System.out.println(" Basis set and geometry data loaded successfully.");
        else {
            System.out.println(" There was an error loading Basis set and/or geometry data.");
            return ;
        }

        // 2) extract MO/NO data
        System.out.println("Getting the MO/NO information from "+Options.MO_File.get_String());
        workhorse._Parse_ASC_File(Options.MO_File.get_String(), workhorse.molden);

        // write a comment to the target molden file
        workhorse.molden.Title = String.format("Created by nwchem2molden from\tlogfile=%s\tand\tascfile=%s", 
                Options.LogFile.get_String(), Options.MO_File.get_String());
        // convert coords. to A.U.
        workhorse.molden.CoordsToAU();

        // 3) save molden file
        if (workhorse.molden.Save_As_MOLDEN(Options.Molden_File.get_String()))
            System.out.println(" Molden file saved to "+Options.Molden_File.get_String());
        else {
            System.out.println(" There was an error saving molden file.");
            return ;
        }
    }

}
