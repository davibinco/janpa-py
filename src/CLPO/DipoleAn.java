package CLPO;

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
 */

import Jama.*;
import java.io.PrintStream;
import onpa.NPA;
import MatrixHelper.*;
import moldenio.*;


/**
 *
 * @author timAdmin
 */
public class DipoleAn {
    
    public Matrix[] AO_dipoleMatr = new Matrix[3]; //[x,y,z], incl. nuclear contribution
    public Matrix[] NAO_dipoleMatr = new Matrix[3]; //[x,y,z], incl. nuclear contribution
    public Matrix[] LHO_dipoleMatr = new Matrix[3]; //[x,y,z], incl. nuclear contribution
    public Matrix[] CLPO_dipoleMatr = new Matrix[3]; //[x,y,z], incl. nuclear contribution
    public Matrix[] LPO_dipoleMatr = new Matrix[3]; //[x,y,z], incl. nuclear contribution
    public double[] true_Dipole_au = new double[3];
    
    private NPA npa = null;
    private PropertyOptimizedOrbitals lpoClpo = null;
    private PrintStream out = System.out;            
    //--------------------------------------------------------------------------
    public DipoleAn(NPA npa, PropertyOptimizedOrbitals lpoClpo) {
        this.npa = npa;
        this.lpoClpo = lpoClpo;
    }
    //--------------------------------------------------------------------------
    /*private void dipCharges() {

        double[] dCharges = new double[ npa.Centers.length ];
        double[] dNPA = new double[ npa.Centers.length ];
        
        for (int b=0; b<npa.NAO.length; b++) {
            int c = npa.NAO[b].Center_ID-1;
            double[] R0 = npa.Centers[c].R0;
            double tmp1 = 0.0;
            double tmp2 = 0.0;
            for (int mu=0; mu<3; mu++) {
                tmp1 += NAO_dipoleMatr[mu].get(b, b)*R0[mu];
                tmp2 += R0[mu]*R0[mu];
            }                
            //out.printf("%d  %.7f%n", npa.NAO[b].Center_ID, -tmp1/tmp2);
            dCharges[c] += -tmp1/tmp2 * npa.SDS_NAO.get(b, b) ;
            dNPA[c] += npa.SDS_NAO.get(b, b) ;
        }
        out.println("NPA dipole-optimized charges:");
        for (int c=0; c<npa.Centers.length; c++) {
            out.printf("%3d  %10.7f  %10.7f%n", c+1, -dNPA[c] + npa.Centers[c].Z ,  -dNPA[c]*2 + npa.Centers[c].Z + dCharges[c]);
        }
        
    }*/
    //--------------------------------------------------------------------------
    private void printDipolarContribs(Matrix[] dipoles, Matrix sds,  String lbl, String[] labels, double printOccThreshold) {
        
        double[] dPrintedSum = new double[3];
        double[] dAllSum = new double[3];
        out.println("    [  DipX/au    DipY/au    DipZ/au ] x Occupancy  from  Contributing  Orbital ->  [    Charge Center (Angstroms)    ]");
        for (int b=0; b<labels.length; b++) {
            //double occ = npa.SDS_NAO.get(b, b);
            double occ = sds.get(b, b);
            if (occ > printOccThreshold) {
                out.printf(" += [%10.5f %10.5f %10.5f] x ( %7.4f) from %22s ->  [%10.5f  %10.5f %10.5f]%n",
                        dipoles[0].get(b, b), dipoles[1].get(b, b), dipoles[2].get(b, b),
                        occ,
                        String.format("%s (%s %d)", labels[b], lbl, b+1),
                        -dipoles[0].get(b, b)*MOLDEN_IO.BohrRadius, -dipoles[1].get(b, b)*MOLDEN_IO.BohrRadius, -dipoles[2].get(b, b)*MOLDEN_IO.BohrRadius
                         );
                double dipAbs = 0.0;                
                for(int mu=0; mu<3; mu++) {
                    dPrintedSum[mu] += dipoles[mu].get(b, b)*occ;
                    dipAbs += dipoles[mu].get(b, b)*occ * dipoles[mu].get(b, b)*occ;
                }
                //out.printf("abs = %.4f %n", Math.sqrt(dipAbs));
            }
            for(int mu=0; mu<3; mu++) {
                dAllSum[mu] += dipoles[mu].get(b, b)*occ;
            }                
            
        }
        out.println("   -------------------------------------------------");
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Localized total:", dAllSum[0], dAllSum[1], dAllSum[2]);        
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Exact:", true_Dipole_au[0], true_Dipole_au[1], true_Dipole_au[2]);
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Diff = Localized-Exact:", dAllSum[0] - true_Dipole_au[0],
                dAllSum[1] - true_Dipole_au[1], dAllSum[2] - true_Dipole_au[2]);
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Locl. in printed:", dPrintedSum[0], dPrintedSum[1], dPrintedSum[2]);
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Diff = Printed-Exact:", dPrintedSum[0] - true_Dipole_au[0],
                dPrintedSum[1] - true_Dipole_au[1], dPrintedSum[2] - true_Dipole_au[2]);
        
        
    }
    //--------------------------------------------------------------------------
    private void LHO_dip(Matrix sds_orb, Matrix[] el_dip_matrices, double occ_print_Thresh, String orbType) {
        out.printf("%n%s dipole moments (neutralized by the equal fraction of nuclear charge):%n", orbType);
        
        out.printf("%5s\t%-7s\t%9s * [%10s  %10s  %10s]%n", "Id", "Atom",
                "occup.",  "dx  ", "dy  ", "dz  " );
        double[] dPrintedSum = new double[3];
        double[] dAllSum = new double[3];

        LOdescription LOs ;
        if (orbType.equals("AHO"))
            LOs = lpoClpo.LPO_descript;
        else
            LOs = lpoClpo.CLPO_descript;
        
        for(int i=0; i<LOs.hostAtomOfHybrid.length; i++) {
            int a = LOs.hostAtomOfHybrid[i];
            double occ = sds_orb.get(i, i);

            for(int mu=0; mu<3; mu++) {
                double tmp = el_dip_matrices[mu].get(i, i) + npa.Centers[a].R0[mu] * 1.0; // '1.0', since this dipole is _to be multiplied_ by occ., as well, as the electronic one);
                dAllSum[mu] += tmp * occ;
            }

            if (occ > occ_print_Thresh) {
                double dip_abs = 0.0;
                out.printf("%5d\t%-7s\t%9.5f * [",
                        i+1,
                        String.format("%s%d",  npa.Centers[a].Name, a+1),
                        occ
                );
                for(int mu=0; mu<3; mu++) {
                    double tmp = el_dip_matrices[mu].get(i, i) + npa.Centers[a].R0[mu] * 1.0; // '1.0', since this dipole is _to be multiplied_ by occ., as well, as the electronic one);
                    if (mu!=2)
                        out.printf("%10.5f  ", tmp); 
                    else
                        out.printf("%10.5f", tmp); 

                    dPrintedSum[mu] += tmp * occ;

                    dip_abs += tmp*tmp * occ*occ;
                }
                dip_abs = Math.sqrt(dip_abs);
                out.printf("], abs. = %10.5f a.u.%n", dip_abs);
            }                
            //out.println();
        }
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Total in printed:", dPrintedSum[0], dPrintedSum[1], dPrintedSum[2]);
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Total in localized:", dAllSum[0], dAllSum[1], dAllSum[2]);        
        out.printf("%25s [%10.5f  %10.5f %10.5f] a.u.%n", "Exact:", true_Dipole_au[0], true_Dipole_au[1], true_Dipole_au[2]);
        
        out.printf("Note: %.5f occupancy threshold was used for printing%n%n", occ_print_Thresh);
    }
    //--------------------------------------------------------------------------
    private void bond_dipoles(Matrix sds_Hybrids, Matrix sds_LOs,
            LOdescription LOs, Matrix[] hybr_dipoles, String orbType) {
        double[] tot_dip = new double[3];
        
        for(int i=0; i<LOs.hostAtomOfHybrid.length; i++) {
            int[] hybrs = LOs.hybridsOfLO[i];
            if (hybrs.length == 1) {
                double occ = sds_LOs.get(i, i);
                if (occ > 1.0) {
                    int a = LOs.hostAtomOfHybrid[hybrs[0]];
                    out.printf("%5d %25s: ", i+1, LOs.LO_labels[i]);
                    
                    out.printf("%10s occ*<h|d|h> = [", "");
                    for(int mu=0; mu<3; mu++) {
                        double tmp = occ * (hybr_dipoles[mu].get(hybrs[0], hybrs[0]) + npa.Centers[a].R0[mu] );
                        out.printf("%10.5f ", tmp );
                        tot_dip[mu] += tmp;
                    }
                    out.printf("]%n");
                }
            }
            if (hybrs.length == 2) {                
                double occ = sds_LOs.get(i, i);
                //if (occ > 1.0) {
                    double ca = LOs.LO_to_Hybrids.get(i, hybrs[0]);
                    double cb = LOs.LO_to_Hybrids.get(i, hybrs[1]);
                    int a = LOs.hostAtomOfHybrid[hybrs[0]];
                    int b = LOs.hostAtomOfHybrid[hybrs[1]];
                    out.printf("%5d %25s: ", i+1, LOs.LO_labels[i]);
                    
                    out.printf("%10s occ*|c1|^2*<h1|d|h1> = [", "");
                    for(int mu=0; mu<3; mu++) {
                        double tmp = occ * ca*ca*(hybr_dipoles[mu].get(hybrs[0], hybrs[0]) + npa.Centers[a].R0[mu] );
                        out.printf("%10.5f ", tmp );
                        tot_dip[mu] += tmp;
                    }
                    out.printf("] ");
                    
                    out.printf("%10s occ*|c2|^2*<h2|d|h2> = [", "");
                    for(int mu=0; mu<3; mu++) {
                        double tmp = occ * cb*cb*(hybr_dipoles[mu].get(hybrs[1], hybrs[1]) + npa.Centers[b].R0[mu] ) ;
                        out.printf("%10.5f ", tmp );                    
                        tot_dip[mu] += tmp;
                    }
                    out.printf("] ");
                    
                    out.printf("%10s occ*c1*c2*2*<h1|d|h2> = [", "");
                    for(int mu=0; mu<3; mu++) {
                        double tmp = occ * 2*ca*cb*hybr_dipoles[mu].get(hybrs[0], hybrs[1]) ;
                        out.printf("%10.5f ", tmp );
                        tot_dip[mu] += tmp;
                    }
                    //out.printf("] ");

                    out.printf("] %n");
                    
                //}                
            }
        }
        //
        out.printf("total dipole without atomic chagres:%n");
        for(int mu=0; mu<3; mu++) {
            out.printf("%10.5f", tot_dip[mu]);
        }
        out.println();
        out.printf("total dipole with    atomic chagres:%n");
        for(int mu=0; mu<3; mu++) {
            for(int a=0; a<npa.Centers.length; a++) {                
                tot_dip[mu] += npa.Centers[a].R0[mu] * npa.NPA_charges[a];
                if (mu==0)
                    out.printf("%d: + %10.5f %n", a+1, npa.Centers[a].R0[mu] * npa.NPA_charges[a]);
            }
            out.printf("%10.5f", tot_dip[mu]);
        }
        out.println();
        
    }
    //--------------------------------------------------------------------------
    public void printDipoles(double occ_print_Thresh) {
    
        out.println();
        out.println("Localized analysis of the dipole moment");
        out.println();        
        out.println("Forming the dipole matrix elements in AO and NAO bases...");
        
        Matrix LHO_2_AO =  lpoClpo.CLPO_descript.NAO_to_Hybrids/* .NAO_to_LHO*/.transpose() .times( npa.NAO_2_AO );        
        Matrix CLPO_2_AO = lpoClpo.CLPO_descript.LO_to_Hybrids/*.CLPO_to_LHO*/.times( LHO_2_AO);//lpoClpo.NAO_to_LHO.transpose() ).times( npa.NAO_2_AO );        
        //Matrix LPO_2_AO = lpoClpo.LPO_to_AHO.times( lpoClpo.NAO_to_AHO.transpose() ).times( npa.NAO_2_AO );
        Matrix LPO_2_AO = lpoClpo.LPO_descript.LO_to_Hybrids.times( lpoClpo.LPO_descript.NAO_to_Hybrids.transpose() ).times( npa.NAO_2_AO );
        
        
        double numEl = npa.D_Global.times(npa.OverlapMatrix).trace();
        for(int mu=0; mu<3; mu++) {
            double nuclDip = 0.0;
            for (int c=0; c<npa.Centers.length; c++) {                        
                nuclDip += npa.Centers[c].R0[mu] * npa.Centers[c].Z;
            }
            
            int sz = npa.Basis.length;
            AO_dipoleMatr[mu] = new Matrix(sz, sz);
            for (int i=0; i<sz; i++) {
                for (int j=0; j<sz; j++) {
                    double tmp = - npa.bsIntegrals.DipoleMatrix[i][j][mu];
                    //tmp += nuclDip / numEl * npa.OverlapMatrix.get(i,j);
                    AO_dipoleMatr[mu].set(i, j, tmp  );
                }
            }
            
            //nuclDip *= 0;
            
            System.out.printf("Dipole [ %s ]: ", "XYZ".charAt(mu) );
            //PrintWriter x = new PrintWriter("xyz".charAt(mu) + "_dMu.txt");
            //x.close();
            double tot_el_dip = npa.D_Global.times(AO_dipoleMatr[mu]).trace();
            out.printf("Total = %7.4f a.u. (%7.4f Debye), nuclear = %8.5f, electronic = %8.5f, %n",
                            tot_el_dip + nuclDip, (tot_el_dip + nuclDip)/0.393430, nuclDip, tot_el_dip);
            true_Dipole_au[mu] = tot_el_dip + nuclDip;
            
            // TODO: execute LDL^T only once, not 3 times...
            
            NAO_dipoleMatr[mu] = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(AO_dipoleMatr[mu], npa.NAO_2_AO );
            LHO_dipoleMatr[mu] = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(AO_dipoleMatr[mu], LHO_2_AO );
            CLPO_dipoleMatr[mu] = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(AO_dipoleMatr[mu], CLPO_2_AO);
            LPO_dipoleMatr[mu] = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(AO_dipoleMatr[mu], LPO_2_AO);
        }

        Matrix sds_clpo = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(npa.SDS, CLPO_2_AO  );
        Matrix sds_lpo = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(npa.SDS, LPO_2_AO  );
        Matrix sds_lho = MatrixHelper.EigenEngine.TransformSymmetricMatrixToNewBasis(npa.SDS, LHO_2_AO  );
        
        // todo: option for printing the same for AHOs
        LHO_dip(sds_lho, LHO_dipoleMatr, occ_print_Thresh, "LHO");
        //out.printf("%.5f", sds_lho.get(121-1, 136-1));
        /*
        for(int mu=0; mu<3; mu++)
            out.printf("%.5f ", sds_clpo.get(123-1, 250-1) * CLPO_dipoleMatr[mu].get(250-1, 123-1));
        out.println();
        for(int mu=0; mu<3; mu++)
            out.printf("%.5f ", sds_clpo.get(250-1, 123-1) * CLPO_dipoleMatr[mu].get(123-1, 250-1));
        out.println();
        */
        
        bond_dipoles(sds_lho, sds_clpo, lpoClpo.CLPO_descript, LHO_dipoleMatr, "");

        out.println();
        out.println("Dipole analysis in CLPO basis");        
        printDipolarContribs(CLPO_dipoleMatr, sds_clpo, "CLPO", lpoClpo.CLPO_descript.LO_labels, 1.0);
        
        /*
        
        out.println();
        out.println("Dipole analysis in NAO basis");
        printDipolarContribs(NAO_dipoleMatr, npa.SDS_NAO , "NAO", npa.PNAO_Labels , 0.1);

        //dipCharges();
                
        out.println();
        out.println("Dipole analysis in LPO basis");
        printDipolarContribs(LPO_dipoleMatr, sds_lpo, "LPO", lpoClpo.LPO_labels, 1.0);

        out.println();
        out.println("Dipole analysis in CLPO basis");        
        printDipolarContribs(CLPO_dipoleMatr, sds_clpo, "CLPO", lpoClpo.CLPO_labels, 1.0);
        */
        
        
    }
    
}
