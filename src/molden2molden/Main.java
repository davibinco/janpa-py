package molden2molden;

/**
 *  A molden2molden program
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


import moldenio.*;
import Polynom3D.*;
import ProgramOptions.WarningManager;
import java.util.*;

/**
 * A 'running' class of the molden2molden program
 * @author timn
 */
public class Main {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception{
        //
        Locale.setDefault(Locale.US);

        

        // Welcome banner

        System.out.println(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        System.out.println(" * molden2molden: a conversion tool for MOLDEN-(semi)compatible files  * ");
        System.out.println(" *        A part of JANPA package,   http://janpa.sourceforge.net      * ");
        System.out.printf (" *                      VERSION: %12s                          * %n","13-01-2019");
        System.out.println(" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *");
        System.out.println();

        System.out.println(" Copyright (c) Tymofii Nikolaienko, 2014");
        System.out.println();


        MoldenConvertor mc = new MoldenConvertor ();
        m2mOptions opts = mc.options;
        
        // Usage info
        if (args.length < 5) {
            System.out.println("Usage:");
            System.out.println();
            System.out.println("  molden2molden -operation1 -operation2 ... -i input-molden -o output-molden [-47 nbofile]");
            System.out.println();
            System.out.println("where input-molden and output-molden are the names of input and output");
            System.out.println("molden-format files, and possible oprtation(s) are:");
            System.out.println();
            System.out.println("General-purpose operations:\n");

            opts._first.PrintDescriptions(null, "%s\t%s%n");
            return;
        }
        
        // Citation notice
        System.out.println(" If any results obtained with this program are used in a publication or ");
        System.out.println(" are used in any other way, please cite this software as following: ");
        System.out.println("   molden2molden (a part of JANPA package, http://janpa.sourceforge.net ),");
        System.out.println("   version 13-01-2019, T.Yu.Nikolaienko");
        System.out.println(" and cite the following paper:");
        System.out.println("  T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem.(2014),");
        System.out.println("  V.1050, P.15-22, DOI: 10.1016/j.comptc.2014.10.002");
        System.out.println();


        
        //----------------------------------------------------------------------
        // parse command line
        int nJobs = 0;
        
        opts._first.LoadOptionsFromCommandLine(args, null, true, true);
        // get job count
        if (opts.op_cart2spher.get_boolean()) ++nJobs;
        if (opts.op_coordsAngstroms.get_boolean()) ++nJobs;
        if (opts.op_from_orca3_BF.get_boolean()) ++nJobs;
        if (opts.op_from_psi4b4_BF.get_boolean()) ++nJobs;
        if (opts.op_normalizeBF.get_boolean()) ++nJobs;
        if (opts.op_orca3signs.get_boolean()) ++nJobs;
        if (opts.op_to_orca3_BF.get_boolean()) ++nJobs;
        if (opts.op_to_psi4b4_BF.get_boolean()) ++nJobs;


        // check for non-controversy of the commands:
        mc.CommandsCheck();

        //
        //=================================================================================================================
        System.out.printf("Total number of jobs to be done: %d%n",nJobs);

        mc.molden_in.Allow_additional_r_power = true;
        mc.CONVERT();

        
        WarningManager.summarizeWarnings();

        System.out.println("All done!");

        // DONE!
        

    }

}
