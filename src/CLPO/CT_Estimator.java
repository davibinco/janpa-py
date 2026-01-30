package CLPO;

/**
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

import Jama.Matrix;
import MatrixHelper.EigenEngine;
import java.io.*;

/**
 *
 * @author timAdmin
 */
public class CT_Estimator {
    static PrintStream out = System.out;
    
    private LOdescription LOs;
    public Matrix SDS_LOs;
    public double ct_threshold = 0.01; // electrons
    
    public CT_Estimator(LOdescription LOs, Matrix sds_NAO) {
        this.LOs = LOs;
        Matrix NAO2LO = LOs.NAO_to_Hybrids.times(LOs.LO_to_Hybrids.transpose());
        SDS_LOs = EigenEngine.TransformSymmetricMatrixToNewBasis(sds_NAO, NAO2LO.transpose() );
        //Matrix SDS_LOs_new = SDS_LOs.copy();
    }
    //--------------------------------------------------------------------------
    
    /**
     * Checks if all hybrids of each a single CLPO belongs to
     * one and the same atomic fragment (should always be this way
     * if fragment identification works properly)
     * @returns fragmentIds of each CLPO
     */
    // TODO: import warning manager
    private int[] ensureCLPOsIntrafragment(int[] fragmentIds) {
        int[] clpo2fragId = new int[LOs.hybridsOfLO.length];
        
        for (int i=0; i<LOs.hybridsOfLO.length; i++) {
            int[] hybrs = LOs.hybridsOfLO[i];
            boolean sameFrag = true;
            int fragPrev = fragmentIds[LOs.hostAtomOfHybrid[hybrs[0]]];
            clpo2fragId[i] = fragPrev; // an optimistic assumption )
            for (int h=1; h<hybrs.length; h++)
                sameFrag &= fragmentIds[LOs.hostAtomOfHybrid[hybrs[h]]] == fragPrev;
            if (!sameFrag) {
                out.printf("ERROR: CLPO %d has hybrids belonging to different fragments!%n", i+1);
                clpo2fragId[i] = -1; // sign of an error
            }
        }
        return clpo2fragId;
    }
    //--------------------------------------------------------------------------
    
    /**
     * 
     * @returns result[] = {frag1Lost, frag2Lost}
     */
    private double perform_CT(int frag1, int frag2, int[] clpo2fragId) {
        //double[] result = new double[2];
        double totCT = 0.0;
        
        double lowest_donor_occ = 1.0;
        double highest_accpt_occ = 1.0;
        int n = clpo2fragId.length;

        int n_CT_belowThreshold = 0;
        double qCT_belowThreshold = 0;

        out.println(" orb.num.   description  occup.  -->  charge, e  -->   occup.   description  orb.num.");
        
        for (int i=0; i<n; i++) {            
            if (clpo2fragId[i] != frag1)
                continue;
            double Dii = SDS_LOs.get(i, i);
            for (int j=0; j<n; j++) {
                if (j==i)
                    continue;
                if (clpo2fragId[j] != frag2)
                    continue;
                double Djj = SDS_LOs.get(j, j);    
                // TODO: count how much were skipped due to  lowest_donor_occ / highest_accpt_occ limitations
                // TODO: sum up Dij^2 and report them at the end
                if ((Dii > lowest_donor_occ) && (Djj < Dii) && (Djj < highest_accpt_occ)) {
                    double Dij = SDS_LOs.get(i, j);
                    double qCT = Dij*Dij / Dii;
                    
                    if (qCT > ct_threshold) {
                        out.printf(" %5d %15s  %7.4f  -->   %7.5f   -->  %7.4f  %15s %5d%n",
                                i+1,  LOs.LO_labels[i], Dii,
                                qCT,
                                Djj,  LOs.LO_labels[j], j+1 );
                    } else {
                        // if not printing
                        n_CT_belowThreshold++;
                        qCT_belowThreshold += qCT;
                    }
                    totCT += qCT;                    
                }
            }
        }
        out.println();
        out.printf("%d orbital pairs with total charge transfer of %.5f were not printed%n",
                n_CT_belowThreshold, 
                qCT_belowThreshold);
        out.println();        
        return totCT;
    }
    //--------------------------------------------------------------------------
    
    public void new_CT(int[] fragmentIds){
        
        out.println("Approximate charge transfer analysis in CLPO basis");
        out.println("NOTE: this is an experimental feature AND IS SUBJECT TO CHANGE!");
        out.println(" We expect to have the underlying theory published soon...");
        out.println("");
        //out.printf("Threshold for printing: %.5f (use %s to adjust)%n", ct_threshold, "<todo>");
        out.printf("Threshold for printing: %.5f %n", ct_threshold);
        out.println("");

        int nTotFramgents = 0;
        for (int i = 0; i < fragmentIds.length; i++)
            if (fragmentIds[i] > nTotFramgents)
                nTotFramgents = fragmentIds[i];
        double[] fragmentFrom = new double[nTotFramgents];
        double[] fragmentTo = new double[nTotFramgents];
        
        
        int n = SDS_LOs.getRowDimension(); //todo: get from LOs.LO_to_Hybrids

        int[] clpo2fragId = ensureCLPOsIntrafragment(fragmentIds);
        // fragmentIds are 1-based !!!

        
        out.println("IntrAfragment charge transfers (conjugation analysis, etc.)");
        for (int iFrag=1; iFrag<=nTotFramgents; ++iFrag) {            
            out.printf("CT within fragment %d%n%n", iFrag);
            perform_CT(iFrag, iFrag, clpo2fragId);            
        }
        
        out.println();
        
        out.println("IntErfragment charge transfers (donor-acceptor analysis, etc.)");
        out.println();
        
        for (int iFrag=1; iFrag<=nTotFramgents; ++iFrag) {
            for (int jFrag=iFrag+1; jFrag<=nTotFramgents; ++jFrag) {
                out.printf("CT between from fragment %d to fragment %d%n", iFrag, jFrag);
                double totCT = perform_CT(iFrag, jFrag, clpo2fragId);
                fragmentFrom[iFrag-1] += totCT;
                fragmentTo[jFrag-1] += totCT;
                
                out.printf("CT between from fragment %d to fragment %d%n", jFrag, iFrag);
                totCT = perform_CT(jFrag, iFrag, clpo2fragId);                
                fragmentFrom[jFrag-1] += totCT;
                fragmentTo[iFrag-1] += totCT;                
            }
        }
                        

        out.println("Inter-molecular charge transfer summary (NEW):");
        out.println("Mol  Accepted -Donated  =  got_electrons");
        //TODO: change sign in the output!
        for (int f=0; f<nTotFramgents; f++) {
            out.printf("%3d  +%.5f -%.5f  =   %+.5f%n", f+1, fragmentTo[f], fragmentFrom[f], fragmentTo[f] - fragmentFrom[f]);
        }

/*        double[] acceptor_occups = new double[n];
        double[] donor_occups = new double[n];
        
        out.println("Occupancy vs total CT:");
        for (int j=0; j<n; j++) {
            double Djj = SDS_LOs.get(j, j);
            if ((Djj > ct_threshold) || ( accaptor_occups[j] > ct_threshold)) {
                out.printf(" %5d %15s  %7.5f  > %8.5f = %7.5f - %7.5f  %s%n", 
                        j+1, LOs.LO_labels[j], Djj, accaptor_occups[j] - donor_occups[j],
                        accaptor_occups[j], donor_occups[j],
                        (Djj > accaptor_occups[j] - donor_occups[j])?("ok"):("strange") );
            }
        }
        out.println();
     */
        
    }
    //--------------------------------------------------------------------------
    
    
    /**
     * 
     * @param sds_NAO
     * @param LOs 
     */
    
}
