/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CLPO;

import JGints.AtomicCenter;
import JGints.BasisFunction;
import java.io.PrintStream;

/**
 *
 * @author timAdmin
 */
public class HybridizAn {
    public static void hybr_an(LOdescription clpo, BasisFunction[] NAOs, AtomicCenter[] centers) {
        PrintStream out = System.out;
        out.println("Hydridization analysis of LHOs:");
        int L_max = 4; // G
        double[][] hybridization = new double[ clpo.hostAtomOfHybrid.length ][L_max+1];
        for (int i=0; i<clpo.hostAtomOfHybrid.length; i++) {
            for (int nao=0; nao<NAOs.length; nao++) {
                double tmp =  clpo.NAO_to_Hybrids.get(nao, i);
                hybridization[i][ NAOs[nao].L ] += tmp * tmp;
            }
            out.printf("%4d %s%d", i+1, 
                    centers[  clpo.hostAtomOfHybrid[i] ].Name,
                    clpo.hostAtomOfHybrid[i]+1 
            );
            for (int L=0; L<=L_max; L++) {
                out.printf("%7.2f ", hybridization[i][L] * 100);
            }
            out.printf("s p ^ (%.1f) d ^ (%.1f)%n", 
                    hybridization[i][1]/hybridization[i][0],
                    hybridization[i][2]/hybridization[i][0]
            );          
        }
        out.println();
/*        
        for (int i=0; i<clpo.hybridsOfLO.length; i++) {
            int[] hybrs = clpo.hybridsOfLO[i];
            if (hybrs.length == 2) {
                
            }
        }*/
    }
    //--------------------------------------------------------------------------
        

}
