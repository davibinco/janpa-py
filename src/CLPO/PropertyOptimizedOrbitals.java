package CLPO;

import JGints.AtomicCenter;
import JGints.BasisFunction;
import Jama.*;
import MatrixHelper.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import onpa.ono_options;
import onpa.printout;

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
public class PropertyOptimizedOrbitals {
    private Matrix SDS_NAO;
    private int nAtoms = 0;
    private int[][] NAOs_at_center;
    private int nNAOs = 0;
    private AtomicHybrids[] hybridsOfAtoms;
    private PrintStream out = System.out;
    private Matrix[][] DAB;
    
    public LOdescription LPO_descript;
    public LOdescription CLPO_descript;
    /*
    // LPO results
    public Matrix NAO_to_AHO; // 1-st index: NAO, 2-nd index: hybrid
    public Matrix  LPO_to_AHO; // 1-st index: LPO, 2-nd index: hybrid // LPO
    public String[] LPO_labels;
    //public String[] AHO_labels; // same as LHO_labels defined below
    public int[][] atomic_AHOs_to_globalHybridIds; // [atomId][0...nHybr 'local' hybrid id]
            
    // CLPO results
    public Matrix NAO_to_LHO; // 1-st index: NAO, 2-nd index: hybrid
    public Matrix CLPO_to_LHO; // 1-st index: LPO, 2-nd index: hybrid // Chemist's LPOs (Lewis-optimized)    
    public String[] CLPO_labels;
    public String[] LHO_labels;
    public int[][] atomic_LHOs_to_globalHybridIds; // [atomId][0...nHybr 'local' hybrid id]
    */
    
    private ono_options options; // misc. work options    
    
    private AtomicCenter[] Centers;
    //--------------------------------------------------------------------------
    /**
     * Uses information from nAtoms and NAOs and @returns an array:
     * result[0-based center ID] = array of NAO IDs centered at this center
     */
    private int[][] get_NAOs_at_centers(BasisFunction[] NAOs) {
        int[][] result = new int[nAtoms][1]; // pre-init
        
        // compute the total number of NAOs at each center
        for (int b=0; b<NAOs.length; b++)
            result[ NAOs[b].Center_ID-1 ] [0] ++;
        
        // allocate memory
        for (int c=0; c<nAtoms; c++)
            result[c] = new int[ result[c][0] ];
        
        // loop over NAOs and fill in the result[][] array
        for (int b=0; b<NAOs.length; b++) {
            int c = NAOs[b].Center_ID-1;
            int lastInd = result[c][ result[c].length - 1  ] ++ ; // use the last element as a counter
            // same as: lastInd <- result[c][-1]; result[c][-1]++;
            result[c][lastInd] = b;
        }
        
        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * 
     * @returns an array[i][j] of D[NAOs_at[i], NAOs_at[j]], valid for: j <= i
     */
    private Matrix[][] diatomic_Submatrices() {
        Matrix[][] result = new Matrix[nAtoms][];
        for (int i=0; i<nAtoms; i++) {
            result[i] = new Matrix[i+1];
            for (int j=0; j<=i; j++) {
                result[i][j] = SDS_NAO.getMatrix(NAOs_at_center[i], NAOs_at_center[j]);
            }
        }
        return result;
    }
    //--------------------------------------------------------------------------    
    /**
     * Returns 2x2 sub-matrix of D in the basis of given hybrids
     * Note that the hybrids MUST be added into the atomic U matrices _before_
     * calling this method !
     */
    private Matrix dihybridMatrix(int atomA, int iA, int atomB, int iB) {
        Matrix DAAh = DAB_dot_h(atomA, atomA, iA); // DAA|hA>
        Matrix DBBh = DAB_dot_h(atomB, atomB, iB); // DBB|hB>
        Matrix DABh = DAB_dot_h(atomA, atomB, iB); // DAB|hB>
        double d11 = hybridsOfAtoms[atomA].hybridScalarMul(iA, DAAh); // <hA| DAA|hA>
        double d22 = hybridsOfAtoms[atomB].hybridScalarMul(iB, DBBh); // <hB| DBB|hB>
        double d12 = hybridsOfAtoms[atomA].hybridScalarMul(iA, DABh); // <hA| DAB|hB>
        // form a matrix from d11, d12 and d22 and return it to a caller
        Matrix result = new Matrix(2,2);
        result.set(0, 0, d11);
        result.set(0, 1, d12);
        result.set(1, 0, d12);
        result.set(1, 1, d22);
        return result;        
    }
    //--------------------------------------------------------------------------    
    
    
    private void _printHybridStatTable(boolean skipPairs, boolean skipSingles, double occThreshold) {
        out.printf(" AtomA (VecId) AtomB (VecId) occupancy     Dij%n");
        double occSumAboveThresh = 0.0;
        double occSumBelowThresh = 0.0;
        int nItemsAboveThresh = 0;
        int nItemsBelowThresh = 0;
        for(int a=0; a<nAtoms; a++) {
            // compute and print occupancies of the hybrids
            for (int h=0; h<hybridsOfAtoms[a].nValidHybrids; h++) {
                
                Matrix hA = hybridsOfAtoms[a].getHybrid(h);
                int b = hybridsOfAtoms[a].friendAtomIndex[h];
                int iB = hybridsOfAtoms[a].friendHybridIndex[h];
                if (skipPairs && (b != -1))
                    continue; // print nothing
                if (skipSingles && (b == -1))
                    continue; // print nothing                
                
                double dii = hA.transpose().times(DAB[a][a]).times(hA).get(0, 0); // <hA|D|hA>
                if (dii < occThreshold) {
                    occSumBelowThresh += dii;
                    nItemsBelowThresh++;
                    continue;
                } else {
                    nItemsAboveThresh++;
                    occSumAboveThresh += dii;
                }
                out.printf(" %5d (%5d) %5d (%5d) %9.5f ", a+1, h+1, // use 1-based numbering for both the atom and its hybrid
                        (b == -1) ? -1 : (b+1), // print -1 or 1-based friend atom number
                        (iB == -1)? -1 : (iB+1), // print -1 or 1-based friend hybrid number
                        dii);

                // Append an off-diagonal element of the density matrix between
                // the two hybrids, or print 'n/a' if this is a non-bonding hybrid
                if ( b != -1 ) {
                    out.printf("%7.3f%n", hA.transpose().times( DAB_dot_h(a, b, iB) ).get(0, 0) );
                    //hybridsOfAtoms[a].hybridScalarMul(h, DAB_dot_h(a, b, iB))
                } else {
                    out.printf("%7s%n", "n/a");
                }
            }
        }
        //out.printf("Sum of occupancies: %.3f in %d items%n", occSum, nItems);
        if (occThreshold > 0){
            out.printf("Total: %.3f electrons = %.3f electrons in %d printed items %n" +
                    "and %.3f electrons in %d items (each having occupancy "+
                    "below threshold = %.4f)%n", 
                        occSumAboveThresh + occSumBelowThresh,
                        occSumAboveThresh, nItemsAboveThresh, 
                        occSumBelowThresh, nItemsBelowThresh,
                        occThreshold);
        } else {
            out.printf("Total: %.3f electrons in %d printed items %n", 
                    occSumAboveThresh, nItemsAboveThresh);
        }
        out.println();
        
    }
    //--------------------------------------------------------------------------
    /**
     * @returns the column-vector equal to 
     *  DAB[@param a][@param b].(@param iB-th hybrid of b-th atom)
     * Note that this method is applicable for both a > b, a==b and a < b
     * Since this method is expected to be called frequently, it avoids
     * calling .transpose() and other 'heavy' matrix operations/computations
     */
    private Matrix DAB_dot_h(int a, int b, int iB) {
        boolean dij_transposed = false;
        int nResComponents;
        Matrix _dij;
        if (b <= a) {
            // Note that DAB[I][J] is valid only for J<=I
            _dij = DAB[a][b];
            nResComponents = _dij.getRowDimension();
        } else {
            _dij = DAB[b][a];
            dij_transposed = true; 
            nResComponents = _dij.getColumnDimension();
        }
        Matrix result = new Matrix(nResComponents, 1);
        int bHybridLength = hybridsOfAtoms[b].U.getRowDimension();
        // compute _dij.(iB-th vec of b) 'manually'
        for (int i=0; i<nResComponents; i++) {
            double tmp = 0;            
            if (! dij_transposed) {
                for (int j=0; j<bHybridLength; j++)
                    tmp += _dij.get(i, j) * hybridsOfAtoms[b].U.get(j, iB);
            } else {
                for (int j=0; j<bHybridLength; j++)
                    tmp += _dij.get(j, i) * hybridsOfAtoms[b].U.get(j, iB);                
            }
            result.set(i, 0, tmp);
        }
        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Computes and @returns the gradient of the DNLO target function w.r.t. 
     * the components of the U matrix of the @param a-th atom
     * If @param pwin != null, pwin[0] receives the 'win' value (i.e., a 
     * contribution of the a-th atom into the target function)
     * 
     */
    private Matrix grad_DNLO_sq(int a, double[] pwin) {
        Matrix Ua = hybridsOfAtoms[a].U; 
        int nVecs = Ua.getColumnDimension();
        int vecLen = Ua.getRowDimension();
        Matrix result = new Matrix(vecLen, nVecs);
        double win = 0.0;
        for (int i=0; i<nVecs; i++) {
            Matrix Daa_ha = DAB_dot_h(a, a, i);            
            double wii = hybridsOfAtoms[a].hybridScalarMul(i, Daa_ha);
            int b = hybridsOfAtoms[a].friendAtomIndex[i];
            
            if (b == -1) {
                // grad[:, i] += Daa_ha * wii
                for (int k=0; k<vecLen; k++)
                    result.set(k, i,  Daa_ha.get(k, 0) * wii  );
                win += wii*wii;                
            } else {                
                Matrix Dab_hb = DAB_dot_h(a, b, 
                        hybridsOfAtoms[a].friendHybridIndex[i] );                
                double wij = hybridsOfAtoms[a].hybridScalarMul(i, Dab_hb);
                // grad[:, i] += Daa_ha * wii + Dab_hb * wij
                for (int k=0; k<vecLen; k++)
                    result.set(k, i,
                            Daa_ha.get(k, 0) * wii + Dab_hb.get(k, 0) * wij );
                win += wii*wii + wij*wij;                
            }
        }
        if (pwin != null)
            pwin[0] = win;
        return result;
    }
    //
    
    private Matrix grad_DNLO_occ2(int a, double[] pwin) { // TODO: correct it
        Matrix Ua = hybridsOfAtoms[a].U; 
        int nVecs = Ua.getColumnDimension();
        int vecLen = Ua.getRowDimension();
        Matrix result = new Matrix(vecLen, nVecs);
        double win = 0.0;
        for (int i=0; i<nVecs; i++) {
            Matrix Daa_ha = DAB_dot_h(a, a, i);            
            double dii = hybridsOfAtoms[a].hybridScalarMul(i, Daa_ha);
            int b = hybridsOfAtoms[a].friendAtomIndex[i];
            
            if (b == -1) {
                // grad[:, i] += Daa_ha * wii
                double wii = dii;
                for (int k=0; k<vecLen; k++)
                    result.set(k, i, result.get(k,i) + Daa_ha.get(k, 0) * wii  );
                    //result.set(k, i, result.get(k,i) + Daa_ha.get(k, 0)   );
                win += dii*dii;
            } else {    
                int hb = hybridsOfAtoms[a].friendHybridIndex[i];
                Matrix Dab_hb = DAB_dot_h(a, b, hb );       
                double dij = hybridsOfAtoms[a].hybridScalarMul(i, Dab_hb);
                Matrix Dbb_hb = DAB_dot_h(b, b, hb);            
                double djj = hybridsOfAtoms[b].hybridScalarMul(hb, Dbb_hb);
                double sqrt = Math.sqrt( (dii-djj)*(dii-djj) + 4*dij*dij );
                double bd_occ = 0.5 * (dii + djj + sqrt);
                double wii = bd_occ * 0.5 * ( 1 + (dii-djj)/sqrt );
                double wij = bd_occ * dij/sqrt;
                
                // grad[:, i] += Daa_ha * wii + Dab_hb * wij
                for (int k=0; k<vecLen; k++)
                    result.set(k, i,  result.get(k, i) + 
                            Daa_ha.get(k, 0) * wii + Dab_hb.get(k, 0) * wij );
                win += bd_occ * bd_occ / 2.0; // since each pair A-B will be 
                // encountered twice: once for atom A and once more for B
            }
        }
        if (pwin != null)
            pwin[0] = win;
        return result;
                
    }
    //
    private boolean opt_Lewis_mode = false; // 
    
    private Matrix grad_DNLO(int a, double[] pwin) {
        if (!opt_Lewis_mode)
            return grad_DNLO_sq(a, pwin); // default mode 
        else
            return grad_DNLO_occ2(a, pwin); // Lewis mode 
    }
    //--------------------------------------------------------------------------
    private int optimizeHybrids_nItersDone = -1;
    /**
     * The key procedure of DNLO which optimizes the hybrids
     * @returns true if iterations converged
     * @param win_prev is the value of the target function at the current 
     * atomic U matrices and the current hybrid pairing. Can be set to an
     * arbitrary negative value if unknown
     */

    //boolean stop = false;
    //boolean converged = false; They were inside optimizeHybrids but annoying there
    private double optimizeHybrids(){//double win_prev) {
        int iter = 0;
        
        Matrix[] G = new Matrix[nAtoms]; // gradients w.r.t. U[a]
        for (int a=0; a<nAtoms; a++)
            hybridsOfAtoms[a].backup_U(); // initialize storage
                
        Matrix[] DAA_hBasis = new Matrix[nAtoms]; // contains < hybrids@A | D | hybrids@A > and is updated at each iteration;
        // used further for computing DOLOs and for controlling optimization convergence
        for(int a=0; a<nAtoms; a++) {
            DAA_hBasis[a] = hybridsOfAtoms[a].U.transpose().times(DAB[a][a]).times(hybridsOfAtoms[a].U);
            // within this method it is used to control convergence only => symmetry is not critical
        }
        double[][] changes = new double[2][2]; // [criteria][i=0:sum, i=1:max]
        
        double[] pwin = new double[]{ 0.0 };
        
        double win_prev = -1; // an arbitrary negative number
        
        
        //t.printf("  Itr.   opt.funct   new-old   ||Unew-U||^2  max||unew-u||^2 ||Dnew-D||^2  max||dnew-d||^2   lambda%n");
        out.println("Iterative optimization of atomic hybrids");
        out.printf("           Target                 Total        max.atomic       Total        max.atomic     Nxt.step%n");
        out.printf("  Itr.    function    new-old  ||Unew-U||^2   ||unew-u||^2   ||Dnew-D||^2   ||dnew-d||^2     lambda %n");
        
        double invlambda = 0.0; // 1/lambda
        boolean dont_scale_stepsize = true; // means that lambda = +infinite and no stepsize scaling is required
        //double prev_dfdlam = 0.0;
        while(true) {
            //TODO?: check if all dij >0
            iter++;
            optimizeHybrids_nItersDone = iter;            
            double win = 0.0;
            
            
            // re-compute G and the target function value for the values of U matrices
            // currently stored in hybridsOfAtoms[a].U
            double invlambda0 = 0.0;
            double dfdlam = 0;
            double g2_tot = 0;
            double g2_max = 0;
            
            for (int a=0; a<nAtoms; a++) {
                // compute gradient
                G[a] = grad_DNLO(a, pwin);
                win += pwin[0];
                
                invlambda0 += Math.abs ( G[a].transpose().times(hybridsOfAtoms[a].U).trace() );
                /*int nDim = hybridsOfAtoms[a].NAO_indices.length;
                double gtu0_ndim = Math.abs( G[a].transpose().times(hybridsOfAtoms[a].U).trace() )/ nDim;
                if (gtu0_ndim > lambda0)
                    lambda0 = gtu0_ndim;
                */
                // tr( G^T.G - G^T.U0.G^T.U0 )
                dfdlam += G[a].times(G[a].transpose()).trace();
                Matrix tmp = G[a].transpose().times(hybridsOfAtoms[a].U);
                dfdlam -= tmp.times(tmp).trace();                
                
            }
            invlambda0 = invlambda0 / nNAOs ; // nNAOs == nDim
            /*dfdlam /= 2.0; 
            dfdlam *= 16;// note also that we're computing gradient which is 4 times smaller than really
            */

            out.printf("%5d %12.4f   %+8.1e", iter, win, (iter==1)?(win):(win - win_prev));
            //out.printf("%5d %12.4f   %+8.1e", iter, win, win - win_prev);
            
            
            // are these current values of hybridsOfAtoms[a].U better than those
            // backuped previously?
            if (win > win_prev) {
                // all ok, the target function has increased!
                invlambda = 1e-16 * invlambda0; // effectively zero
                dont_scale_stepsize = true;
                //prev_dfdlam = dfdlam;
                // notice that win was calculated at hybridsOfAtoms[a].U values
                // while win_prev comes from their older values; hence since
                // win > win_prev, the current hybridsOfAtoms[a].U values is
                // the best approximation we currently have and we must save it
                // since we don't know yet, whether the newer U obtained from 
                // G[a] below are better -- in fact, we will find it out only
                // at the next iteration
                for (int a=0; a<nAtoms; a++) {
                    hybridsOfAtoms[a].backup_U();
                }
            
                // check if converged
                //if (maxUDiff < 0.01 ) { // works badly!!!
                if (win - win_prev < target_function_conv_thresh) {
                    out.printf("    << Converged! (threshold = %.2e) >>%n", target_function_conv_thresh);
                    out.println("Optimization of hybrids finished");
                    //stop = true;
                    //converged = true;
                    return win;
                }
                win_prev = win;         // save the current target function 
                // value -- it is the highest one we've ever encountered!                
            } else {
                // IF win < win_prev:
                // the values currently saved in hybridsOfAtoms[a].U values gave 
                // worse results than those backuped previously => we can safely 
                // delete the current ('bad') values of hybridsOfAtoms[a].U and 
                // replace them with our previously saved best guess
                for (int a=0; a<nAtoms; a++)
                    hybridsOfAtoms[a].restore_U();
                                
                if (win_prev - win < target_function_conv_thresh) {            
                    out.println("    << No step can be taken >>");
                    out.println("Optimization of hybrids finished");
                    return win_prev;
                }
                //converged = false;
                //return false;
                
                double tmp;
                if (dont_scale_stepsize) {
                    // ha! this is the first unsucessfull iteration! => we MUST
                    // initialize lambda with some reasonable value
                    invlambda = invlambda0;
                    //tmp = 1.0;
                } else {
                    // hm... lambda already had some value, but it was too big in 
                    // order to restore convergence => we must reduce the stepsize
                    // (reduce lambda and increase 1/lambda)
                    invlambda *= 2.0;
                    //tmp = 1/invlambda;
                }                
                //tmp = -0.5*prev_dfdlam*tmp*tmp / (win - win_prev - prev_dfdlam*tmp) ;
                //invlambda = 1/tmp;
                
                dont_scale_stepsize = false;
                //return false;
            }
            
            if (iter >= Opt_Max_Iter) {
                out.println("    << Maximum number of iterations reached >>");
                out.println("Optimization of hybrids finished");
                for (int a=0; a<nAtoms; a++)
                    hybridsOfAtoms[a].restore_U();                
                return win_prev; // this value corresponds to lastly saved U matrices
                //break;
            }
            
            for (int i=0; i<changes.length; i++)
                for(int j=0; j<changes[i].length; j++) changes[i][j] = 0.0;

            
            
            // calculate new values for hybridsOfAtoms[a].U using the stepsize scaling
            // obtained above
            for (int a=0; a<nAtoms; a++) {            
                // Ua_new = g.(g^T.g)^(-0.5) where g = G + invlambda * U0, which is 
                // essentially the same as with g = G*(1/invlambda) + U0, as in the paper
                Matrix Ua_new = EigenEngine.symmetrOrth( G[a].plus(  hybridsOfAtoms[a].U.times(invlambda) ) );
                /*
                // g = (1-lam)*U0 + lam*G 
                Matrix Ua_new ;
                if (dont_scale_stepsize)
                    Ua_new = EigenEngine.symmetrOrth( G[a] );
                else {
                    Matrix g = hybridsOfAtoms[a].U.times(1-1/invlambda).plus( G[a].times(1/invlambda)  );
                    Ua_new = EigenEngine.symmetrOrth( g );
                }*/
                
                // compute || Unew - Uold ||^2 for the a-th atom
                double du2 = Ua_new.minus( hybridsOfAtoms[a].U ).normF();
                changes[0][0] += du2;                
                int nDim = hybridsOfAtoms[a].NAO_indices.length;
                if (du2/nDim > changes[0][1])
                    changes[0][1] = du2/nDim;
                
                // compute || Unew^T.Daa.Unew - Uold^T.Daa.Uold ||^2 for the a-th atom
                Matrix daa_new = Ua_new.transpose().times(DAB[a][a]).times(Ua_new);
                /*double */du2 = DAA_hBasis[a].minus( daa_new ).normF(); // || Dnew - Dold ||^2
                DAA_hBasis[a] = daa_new; // update the old value
                changes[1][0] += du2;                
                if (du2 > changes[1][1])
                    changes[1][1] = du2;

                // transfer the elements of Unew into this.hybridsOfAtoms.U
                // These values will be assessed in comparison with the saved ones at the next iteration
                hybridsOfAtoms[a].U = Ua_new;
            }
            //out.printf("%5d %12s   %10s %12.7f  %15.7f", iter, "", "next:", changes[0][0], changes[0][1]); //  Itr.    opt.funct   ||Unew-U||^2  max||unew-u||^2 
            
            out.printf("%12.6f  %13.7f %12.5f    %13.7f   ",//%12.6f", 
                        changes[0][0], changes[0][1], // ||Unew-U||^2  max||unew-u||^2 
                        changes[1][0], changes[1][1]  // ||Dnew-D||^2  max||dnew-d||^2 
                        //,dfdlam
                    );
            if (dont_scale_stepsize)
                out.printf("%10s%n", "(none)");
            else
                out.printf("%10.2e%n", 1/invlambda);                       
        }
        //return false; // never executed
        //return converged;
    }
    //--------------------------------------------------------------------------
    /**
     * @Returns 'global' indices for all hybrids; result[atomId][localHybridId]
     * Note: only .nValidHybrids are used for each atom!
     * 
     * --if @param pDuplicate != null, a copy of result is created in pDuplicate[0]
     * 
     */
    private int[][] getGlobalHybridIndices(){//LOdescription hybridInfo_holder){//int[][][] pDuplicate) {
        int[][] addr2iH = new int[nAtoms][]; // [atom][hybridId] -> global hybrid Id                
        int iH = 0; // global hybrid number 
        for (int a=0; a<nAtoms; a++) {
            int[] a_naos = hybridsOfAtoms[a].NAO_indices;            
            addr2iH[a] = new int[a_naos.length];
            for (int h=0; h<hybridsOfAtoms[a].nValidHybrids; h++) {
                addr2iH[a][h] = iH;
                iH++;
            }
        }
        
        /*
        // create a copy of addr2iH and place it into pDuplicate[0], if pDuplicate != null:
        if (pDuplicate != null) {
            pDuplicate[0] = new int[nAtoms][];
            for (int a=0; a<nAtoms; a++) {
                pDuplicate[0][a] = addr2iH[a].clone();
            }
        }
        */
        
        return addr2iH;
    }
    //--------------------------------------------------------------------------
    /**
     * Transforms and @returns SDS matrix as U^T.SDS_NAO.U where U is constructed
     * from hybridsOfAtoms[a].U 'diagonal' blocks
     * @param hybrAddresses[atomId][hybridId] is used to determine the 'place'
     * of the specific hybrid in 'global' SDS matrix
     * 
     * --if @param NAO2HO != null, NAO2HO[0] gets the 'global' NAO->HO transfromation
     * --matrix: NAO2HO[0][nao_index][hybrid_index]; in this the matrix to be placed 
     * --into the [0]-th element of the array is NOT allocated by this method itself
     * --but must be pre-allocated by the caller !
     * 
     * if @param NAO2HO_holder != null, its .NAO_to_Hybrids matrix gets the 'global' NAO->HO transfromation,
     * that is: [nao_index][hybrid_index]; 
     * Note that .NAO_to_Hybrids matrix must be allocated _before_ calling this method !
     * 
     */
    private Matrix SDS_in_hybrid_basis(int[][] hybrAddresses, LOdescription NAO2HO_holder){ //Matrix[] NAO2HO) {
        /** determine the size of SDS matrix in hybrid basos (nTotHybrids x nTotHybrids)**/
        int nTotHybrids = 0;
        for (int a=0; a<nAtoms; a++) nTotHybrids += hybridsOfAtoms[a].nValidHybrids;
        
        Matrix D_OHO = new Matrix(nTotHybrids, nTotHybrids);
        /*if (NAO2HO != null)
            NAO2HO[0] = new Matrix(nNAOs, nNAOs ); // 1-st index: NAO, 2-nd index: hybrid*/
                
        for (int a=0; a<nAtoms; a++) {
            int[] a_naos = hybridsOfAtoms[a].NAO_indices;            
            //if (NAO2HO != null) {
            if (NAO2HO_holder != null) {
                for (int h=0; h<hybridsOfAtoms[a].nValidHybrids; h++) {
                    Matrix vec = hybridsOfAtoms[a].U.getMatrix(0, a_naos.length-1, h, h );
                    int iH = hybrAddresses[a][h];
                    //NAO2HO[0].setMatrix(a_naos , iH, iH, vec );
                    NAO2HO_holder.NAO_to_Hybrids.setMatrix(a_naos , iH, iH, vec );
                }
            }
            // now
            // transform D_AA to Hybrid basis            
            Matrix daa = hybridsOfAtoms[a].U.transpose().times(DAB[a][a]).times(hybridsOfAtoms[a].U);
            D_OHO.setMatrix(a_naos, a_naos,  daa.times(0.5)  ); // we're going to do D_OHO += D_OHO^T, 
                                                                // hence scale diagonal blocks by 0.5
            // compute one half of off-diagonal blocks
            for (int b=0; b<a; b++) {
                int[] b_naos = hybridsOfAtoms[b].NAO_indices;
                Matrix dab = hybridsOfAtoms[a].U.transpose().times(DAB[a][b]).times(hybridsOfAtoms[b].U);
                D_OHO.setMatrix(a_naos, b_naos,  dab);
            }
        }
        D_OHO = D_OHO.plus(D_OHO.transpose()); // symmetrizes D_AA and moves D_AB^T into D_BA sub-blocks
        return D_OHO;
    }
    //--------------------------------------------------------------------------
    private void printBondMatrix(int[][] numBonds) {
        out.println("Number of two-center(2C) BD orbitals for each pair of atoms");
        out.printf("%10s","Centr. A/B");
        for (int i=0; i<nAtoms; ++i) out.printf("%10d",i+1);
        out.println();
        for (int i=0; i<nAtoms; ++i) {
            out.printf("%7d   ",i+1); // row number
            for (int j=0; j<i; ++j) out.printf("%10s","");  // empty space
            //int totBondsHere = 0;
            //for (int j=0; j<nAtoms; ++j) totBondsHere += numBonds[i][j];
            //out.printf(" (%7d)", totBondsHere );  // diagonal elements
            out.printf(" (%7d)", numBonds[i][i] );  // diagonal elements
            //numBonds[i][i] = totBondsHere;
            for (int j=(i+1); j<nAtoms; ++j)
                out.printf("%10d", numBonds[i][j]);   // off-diaginal bond indices
            out.println();
        }                
    }
    //--------------------------------------------------------------------------
    /*private void printBondMatrix(double[][] vals, String comment) {
        out.println(comment);
        out.printf("%10s","Centr. A/B");
        for (int i=0; i<nAtoms; ++i) out.printf("%10d",i+1);
        out.println();
        for (int i=0; i<nAtoms; ++i) {
            out.printf("%7d   ",i+1); // row number
            for (int j=0; j<i; ++j) out.printf("%10s","");  // empty space
            out.printf(" (%7.5f)", vals[i][i]);  // diagonal elements
            for (int j=(i+1); j<nAtoms; ++j)
                out.printf("%10.5f", vals[i][j]);   // off-diaginal bond indices
            out.println();
        }                
    }*/
    //--------------------------------------------------------------------------
    /**
     * Counts the number of 2C-hybrids belonging to each atomic pair and
     * @returns it
     */
    private int[][] numBonds() {
        int[][] result = new int[nAtoms][nAtoms];
        
        for (int a=0; a<nAtoms; a++) {
            for (int ha=0; ha<hybridsOfAtoms[a].nValidHybrids; ha++) {
                int b = hybridsOfAtoms[a].friendAtomIndex[ha];
                if (b != -1)
                    result[a][b]++; // note that we rely on the fact that
                // the hybridsOfAtoms[].friendAtomIndex[] is filled in
                // in a 'symmetrical' manned, i.e., a.ha -> b.hb if and
                // only if a.ha <- b.hb
            }
        }
        // fill in diagonal elements with the total number of bonds for each atom
        for (int a=0; a<nAtoms; a++) {
            result[a][a] = 0; // well, it must had been 0 anyway...
            for (int b=0; b<nAtoms; b++)
                if (b != a)
                    result[a][a] += result[a][b];
        }
        
        return result;
                
    }
    //--------------------------------------------------------------------------
    /**
     * @returns true if any of the elements in the @param bOld and @param bNew
     * array differs
     */
    /*private boolean numBonds_differs(int[][] bOld, int[][] bNew) {        
        for(int a=0; a<nAtoms; a++) {
            for(int b=0; b<nAtoms; b++)
                if (bOld[a][b] != bNew[a][b])
                    return true;
        }
        return false;
    }*/
    //--------------------------------------------------------------------------
    
    private void CS_Guess() {
        Matrix[] Ui = new Matrix[nAtoms];
        for (int a=0; a<nAtoms; a++) {
            // prepare matrices
            Matrix[] dgs = new Matrix[nAtoms];            
            for (int b=0; b<nAtoms; b++) {
                if (b<=a)
                    dgs[b] = DAB[a][b].times(DAB[a][b].transpose());
                else
                    dgs[b] = DAB[b][a].transpose().times(DAB[b][a]);
            }            
            // and diag. them
            int nMaxIters = NAOs_at_center[a].length * 100;
            int[] itrs = new int[]{ nMaxIters };
            Ui[a] = CS96_SimultDiag.simult_diag(dgs, itrs, null);
            if (itrs[0] >= nMaxIters) {
                out.printf("Warning (non-critical): simultaneous diag. did not converge when "+
                        "generating initial AO-to-AHO guess for atom %d%n", a+1);
            }
            hybridsOfAtoms[a].U = Ui[a];
        }
        
        for (int a=0; a<nAtoms; a++) {
            hybridsOfAtoms[a].nValidHybrids = hybridsOfAtoms[a].NAO_indices.length;
        }        
    }
    /**
     * Computes an initial guess for the Atomic Hybrid Orbitals (AHOs) by
     * essentially using a simultaneous diag method from [ Cardoso, J. F., 
     * Souloumiac, A. (1996), SIAM journal on matrix analysis and applications, 
     * 17(1), 161-164.]. The desired transformation matrix is initialize
     * with the McWeeny hybrids (eigenvalues of Daa) for each atom
     * The found unitary transformation matrices are written directly
     * into hybridsOfAtoms[a].U and the .nValidHybrids fields are
     * set accordingly
     */
    private void CS_Guess2() {
        Matrix[] Ui = new Matrix[nAtoms];
        for (int a=0; a<nAtoms; a++) {
            // prepare matrices
            Matrix[] dgs = new Matrix[nAtoms]; // array of matrices to be 
            // diagonalized computed as [b] = Dab.Dba in McWeeny basis of th
            // a-th atom, i.e., [b] = Ua^T.Dab.Dab^T.Ua
            Matrix Ua = DAB[a][a].eig().getV();
            for (int b=0; b<nAtoms; b++) {
                Matrix tmp;
                if (b<=a) {
                    // b<=a => ok to use DAB[a][b]
                    tmp = Ua.transpose().times( DAB[a][b] );
                    dgs[b] = tmp.times(tmp.transpose()); // (U^T.Dab).(U^T.Dab)^T
                } else {
                    // b>a => must use DAB[b][a] since in DAB the 2-nd index must 
                    // not exceed the first one
                    tmp = DAB[b][a].times( Ua );
                    dgs[b] = tmp.transpose().times(tmp); // (Dba.Ua)^T.(Dba.Ua)
                }
            }            
            // and diag. them
            Ui[a] = CS96_SimultDiag.simult_diag(dgs, null, Ua );
            hybridsOfAtoms[a].U = Ui[a];
        }
        
        for (int a=0; a<nAtoms; a++) {
            hybridsOfAtoms[a].nValidHybrids = hybridsOfAtoms[a].NAO_indices.length;
        }        
    }    
    //--------------------------------------------------------------------------
    private void CS_Guess3() {
        Matrix[] Ui = new Matrix[nAtoms];
        for (int a=0; a<nAtoms; a++) {
            // prepare matrices
            Matrix[] dgs = new Matrix[nAtoms]; // array of matrices to be diagonalized
            for (int b=0; b<nAtoms; b++) {
                if (b<=a)
                    dgs[b] = DAB[a][b].times(DAB[a][b].transpose());
                else
                    dgs[b] = DAB[b][a].transpose().times(DAB[b][a]);
            }                        
            // create an initial guess            
            int nDim_a = NAOs_at_center[a].length;
            
            Matrix g = new Matrix(nDim_a, nDim_a);
            // (g)_ij = sum_k( (M_k)_ij * (M_k)_jj  )
            for (int b=0; b<nAtoms; b++) {                
                for(int mu=0; mu<nDim_a; mu++)
                    for(int nu=0; nu<nDim_a; nu++)
                        g.set(mu, nu,  g.get(mu, nu) + dgs[b].get(mu, nu) * dgs[b].get(nu, nu) );
            }
            g.plusEquals( Matrix.identity(nDim_a, nDim_a).times(1) );
            Matrix u0 = EigenEngine.symmetrOrth(g);
            // transform all matrices into new basis, ensuring that the resulting matrices are all symmetrical
            for (int b=0; b<nAtoms; b++) {
                Matrix tmp;
                if (b<=a) {
                    // b<=a => ok to use DAB[a][b]
                    tmp = u0.transpose().times( DAB[a][b] );
                    dgs[b] = tmp.times(tmp.transpose()); // (U^T.Dab).(U^T.Dab)^T
                } else {
                    // b>a => must use DAB[b][a] since in DAB the 2-nd index must 
                    // not exceed the first one
                    tmp = DAB[b][a].times( u0 );
                    dgs[b] = tmp.transpose().times(tmp); // (Dba.Ua)^T.(Dba.Ua)
                }
                //dgs[b].print(5, 2);
            }            
            
            // and diag. them
            Ui[a] = CS96_SimultDiag.simult_diag(dgs, null, u0 );
            hybridsOfAtoms[a].U = Ui[a];
        }
        
        for (int a=0; a<nAtoms; a++) {
            hybridsOfAtoms[a].nValidHybrids = hybridsOfAtoms[a].NAO_indices.length;
        }        
    }
    
    //--------------------------------------------------------------------------

    double reconnectHybrids(Matrix D_in_hybrid_basis, int[][] hybrAddresses, boolean allowall)
    {
        double maxBondIonicityThreshold = options.maxClpoBondIonicityThreshold.get_double();

        out.println("Finding an optimal hybrid pairing...");
        boolean do_print = false;

        int[][] nao_owners = new int[nNAOs][];
        List<String> graphTable = new ArrayList<>();

        graphTable.add(String.format(
                "%-4s %-30s %-10s %s",
                "ID", "Description", "Occupancy", "Composition"
        ));

        int nEdgesMax = nNAOs * (2 * nNAOs - 1) / 2;
        if (opt_Lewis_mode)
            nEdgesMax += nNAOs;

        PairEdge[] edges_new = new PairEdge[nEdgesMax];
        int k = 0;

        // -------- Build edges --------
        for (int a = 0; a < nAtoms; a++) {
            for (int ha = 0; ha < hybridsOfAtoms[a].NAO_indices.length; ha++) {

                hybridsOfAtoms[a].friendAtomIndex[ha] = -1;
                hybridsOfAtoms[a].friendHybridIndex[ha] = -1;

                int iA = hybrAddresses[a][ha];
                nao_owners[iA] = new int[]{a, ha};

                double daa = D_in_hybrid_basis.get(iA, iA);

                if (opt_Lewis_mode)
                    edges_new[k++] = new PairEdge(iA, nNAOs + iA, daa * daa);

                for (int b = 0; b < nAtoms; b++) {
                    if (b == a) continue;

                    for (int hb = 0; hb < hybridsOfAtoms[b].NAO_indices.length; hb++) {

                        int iB = hybrAddresses[b][hb];
                        if (iB > iA) continue;

                        double dbb = D_in_hybrid_basis.get(iB, iB);
                        double dab = D_in_hybrid_basis.get(iA, iB);

                        double bd_occ = 0.5 * (daa + dbb +
                                Math.sqrt((daa - dbb) * (daa - dbb) + 4 * dab * dab));
                        double nb_occ = daa + dbb - bd_occ;

                        boolean f = (bd_occ > 1.0) && (nb_occ < 1.0) &&
                                Math.cos(Math.atan2(2 * dab, Math.abs(daa - dbb)))
                                        < maxBondIonicityThreshold;

                        if (allowall || f) {
                            edges_new[k++] = opt_Lewis_mode
                                    ? new PairEdge(iA, iB, bd_occ * bd_occ)
                                    : new PairEdge(iA, iB, dab * dab);
                        }
                    }
                }
            }
        }

        if (k == 0) return 0.0;

        PairEdge[] edges_new2 = new PairEdge[k];
        System.arraycopy(edges_new, 0, edges_new2, 0, k);
        edges_new = edges_new2;

        Blossom x = new Blossom(edges_new);
        int[] remap = x.maxWeightMatching();

        double result = 0.0;
        int rowId = 1;

        // -------- Build CLPO table --------
        for (int i = 0; i < remap.length && i < nNAOs; i++) {

            int j = remap[i];

            if (j >= nNAOs || j == -1) {
                double occ = D_in_hybrid_basis.get(i, i);

                int a = nao_owners[i][0];
                int ha = nao_owners[i][1];

                graphTable.add(String.format(
                        "%4d (LP) %-26s %10.5f   1.0 * h%d@%s%d",
                        rowId++,
                        Centers[a].Name + (a + 1),
                        occ,
                        ha + 1,
                        Centers[a].Name,
                        a + 1
                ));

                result += occ * occ;
                continue;
            }

            if (i < j) {
                double d11 = D_in_hybrid_basis.get(i, i);
                double d12 = D_in_hybrid_basis.get(i, j);
                double d22 = D_in_hybrid_basis.get(j, j);

                int a = nao_owners[i][0];
                int ha = nao_owners[i][1];
                int b = nao_owners[j][0];
                int hb = nao_owners[j][1];

                hybridsOfAtoms[a].friendAtomIndex[ha] = b;
                hybridsOfAtoms[a].friendHybridIndex[ha] = hb;
                hybridsOfAtoms[b].friendAtomIndex[hb] = a;
                hybridsOfAtoms[b].friendHybridIndex[hb] = ha;

                double bd_occ = 0.5 * (d11 + d22 +
                        Math.sqrt((d11 - d22) * (d11 - d22) + 4 * d12 * d12));
                double nb_occ = d11 + d22 - bd_occ;

                double io = Math.cos(Math.atan2(2 * d12, Math.abs(d11 - d22)));
                double theta = 0.5 * Math.atan2(2 * d12, d11 - d22);

                double c = Math.cos(theta);
                double s = Math.sin(theta);

                String atomA = Centers[a].Name + (a + 1);
                String atomB = Centers[b].Name + (b + 1);

                graphTable.add(String.format(
                        "%4d (BD) %-26s %10.5f   h%d@%s * (%7.4f) + h%d@%s * (%7.4f)",
                        rowId++,
                        atomA + "-" + atomB + String.format(", Io = %.4f", io),
                        bd_occ,
                        ha + 1, atomA,  c,
                        hb + 1, atomB, -s
                ));

                graphTable.add(String.format(
                        "%4d     %-26s %10.5f   h%d@%s * (%7.4f) + h%d@%s * (%7.4f)",
                        rowId++,
                        atomA + "-" + atomB + ", antibonding (NB)",
                        nb_occ,
                        ha + 1, atomA, -s,
                        hb + 1, atomB, -c
                ));

                result += opt_Lewis_mode
                        ? bd_occ * bd_occ
                        : d11 * d11 + 2 * d12 * d12 + d22 * d22;
            }
        }

        // -------- Write table to file --------
        try (PrintWriter pw = new PrintWriter(new FileWriter("graph"))) {
            for (String line : graphTable)
                pw.println(line);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return result;
    }

    //--------------------------------------------------------------------------
    /**
     * Creates a set of localized one-/two-center orbitals using the 
     * previuosly optimized atomic hybrids
     * 
     * --TEST: @returns false if the structure has not been validated (e.g., if some 
     * --strange bonds have been found) and re-optimization might be needed
     * 
     * --if @param LO2HYBR is not null, @param LO2HYBR[0] (which MUSH already be allocated!)
     * --receives the LO -> HYBRIDS transfromation matrix
     * --if @LO_labels is not null, its elements (which should already exist!) get
     * --short-hand labels of the created localized orbitals
     * 
     * if @param lo_info != null, the fields of this object are filled in with
     * appropriate 'public' information about the created LOs
     * 
     * @returns D in AHO/LHO basis
     */
    private Matrix createLOs( String loNameStr, LOdescription lo_info){ //Matrix[] pLO2hAO, Matrix[] pNAO2hAO, String[] LO_labels ) {
                
        int[][] addr2iH = getGlobalHybridIndices();//null);
        //Matrix[] pNAO2HO = new Matrix[] { null };
        // compute NAO->HO transformation and D in AHO basis

        
        Matrix D_AHO  = SDS_in_hybrid_basis(addr2iH, lo_info);//pNAO2hAO);
                
        int iLPO = 0;
                
        boolean[] isLewis = new boolean[nNAOs];

        double Ry_print_threshold = options.RyOccPrintThreshold.get_double();
                
        out.println();
        out.printf("*** Summary of %s results%n", loNameStr);
        out.println();
        out.printf(" %5s\t%35s\t%9s\t%s%n", loNameStr, " D e s c r i p t i o n     ", "Occupancy","Composition");

        
        // compute LP + BD populations
        double L2_sum = 0.0;
        double LP_sum = 0.0;  int num_lp = 0;
        double RY_sum = 0.0;  int num_ry = 0;
        double BD_sum = 0.0;  int num_bd = 0;
        double NB_sum = 0.0; 
        double ry_below_thresh_occs = 0.0;
        int num_ry_below_thresh = 0;
        
        // loop over all hybrids and see if they have a partner
        for (int a=0; a<nAtoms; a++) {
            for (int ha=0; ha<hybridsOfAtoms[a].nValidHybrids; ha++) {
                int iHa = addr2iH[a][ha];
                if (lo_info != null) {
                    lo_info.hostAtomOfHybrid[iHa] = a;
                }
                
                int b = hybridsOfAtoms[a].friendAtomIndex[ha];
                if (b == -1) {
                    //if (lo_info != null) {
                        //lo_info.partnerHybridGlobalIndex[iHa] = -1;
                    //}
                    
                    // no partner => LP or RY
                    double occ = D_AHO.get(iHa, iHa);
                    isLewis[iLPO] = occ > 1.0;
                    
                    if (lo_info != null) {                         
                        lo_info.LO_to_Hybrids.set(iLPO, iHa, 1.0); 
                        lo_info.hybridsOfLO[iLPO] = new int[]{iHa};                        
                        //lo_info.hostLOindex[iHa] = iLPO;
                        // 02.07.2018: fixed: iLPO++ was here; moved below! (bugfix)
                    }                    
                    
                    iLPO++;
                                        
                    String host_info = String.format("%s%d", Centers[a].Name, a+1);
                    String composition = String.format("1.0 * h%d@%s", iHa+1, host_info);

                    if (occ > 1.0) {                        
                        LP_sum += occ;                        
                        L2_sum += occ*occ;
                        num_lp++;
                        out.printf("%5d\t%-35s\t%-9.5f\t%s%n", 
                                iLPO, "(LP)  "+host_info, occ, composition);
                        
                        if (lo_info != null) {
                            lo_info.LO_labels[iLPO-1] = host_info+":LP"; // '-1' due to above-made 'iLPO++'
                            lo_info.LO_types [iLPO-1] = LOdescription.LO_type_LP; // '-1' due to above-made 'iLPO++'
                        }
                    } else {
                        RY_sum += occ;
                        num_ry++;
                        
                        if (lo_info != null) {
                            lo_info.LO_labels[iLPO-1] = host_info+":RY"; // '-1' due to above-made 'iLPO++'
                            lo_info.LO_types [iLPO-1] = LOdescription.LO_type_RY; // '-1' due to above-made 'iLPO++'
                        }
                                                
                        if (occ > Ry_print_threshold) {
                            out.printf("%5d\t%35s\t%9.5f\t%s%n", 
                                    iLPO, host_info+" (RY)", occ, composition);
                        } else {
                            ry_below_thresh_occs += occ;
                            num_ry_below_thresh ++;
                        }
                        
                    }
                    //out.printf("(%4s %5d) RY: %3d [ho%5d]            : occ = %7.4f%n", loNameStr, iLPO, a+1, iHa+1, occ); // TODO: ionicity, coefs.                    
                    
                } else {
                    // has partner => bonding
                    int hb = hybridsOfAtoms[a].friendHybridIndex[ha];
                    int iHb = addr2iH[b][ hb ];
                    //if (lo_info != null) {
                        //lo_info.partnerHybridGlobalIndex[iHa] = iHb; // bijective: [iHa]->iHb & [iHb]->iHa => no need to test if iHb>iHa
                    //}
                    
                    if (iHb > iHa) {
                        // note that dihybridMatrix() is inapplicable here since
                        // it works on DAB array which contains SDS_NAO sub-
                        // blocks rather than D_OHO elements !
                        double d11 = D_AHO.get(iHa, iHa);
                        double d12 = D_AHO.get(iHa, iHb);
                        double d22 = D_AHO.get(iHb, iHb);
                        Matrix d = new Matrix(2,2);
                        d.set(0, 0, d11);
                        d.set(0, 1, d12);
                        d.set(1, 0, d12);
                        d.set(1, 1, d22);
                        Object[] e = Sorter.sorted_eigenData(d.eig(), -1);
                        double[] vals = (double[])e[0]; // vals[] are sorted in sort in descending order => [0] = BD, [1] = NB
                        BD_sum += vals[0];
                        NB_sum += vals[1];
                        Matrix cf = (Matrix)e[1];

                        String host_atoms = String.format("%s%d-%s%d", 
                                Centers[a].Name, a+1, 
                                Centers[b].Name, b+1);
                        
                        // bonding (BD) orbital
                        if (lo_info != null) {
                            lo_info.LO_to_Hybrids.set(iLPO, iHa , cf.get(0, 0) );
                            lo_info.LO_to_Hybrids.set(iLPO, iHb , cf.get(1, 0) );
                            lo_info.hybridsOfLO[iLPO] = new int[]{iHa, iHb};
                            lo_info.LO_types[iLPO] = LOdescription.LO_type_BD;
                            lo_info.LO_labels[iLPO] = host_atoms+":BD"; 
                        }
                        
                        isLewis[iLPO] = true;
                        iLPO++;
                        // we don't print 'local' hybrid numbers; these are the 'global'
                        // hybrid numbers which are exported as .molden and hence
                        // can be visualized and analyzed, - not the 'local' ones                        

                        // antibonding (NB) orbital
                        if (lo_info != null) {
                            lo_info.LO_to_Hybrids.set(iLPO, iHa , cf.get(0, 1) );
                            lo_info.LO_to_Hybrids.set(iLPO, iHb , cf.get(1, 1) );
                            lo_info.hybridsOfLO[iLPO] = new int[]{iHa, iHb};
                            lo_info.LO_types[iLPO] = LOdescription.LO_type_NB;
                            lo_info.LO_labels[iLPO] = host_atoms+":NB";
                        }
                        
                        iLPO++;                        
                        
                        
                        String info1 = String.format("%s * (%7.4f) + %s * (%7.4f)", 
                                "h" + String.format("%d@%s%d", iHa+1, Centers[a].Name, a+1),
                                 cf.get(0, 0),
                                 "h" + String.format("%d@%s%d", iHb+1, Centers[b].Name, b+1),
                                 cf.get(1, 0) );
                        
                        String info2 = String.format("%s * (%7.4f) + %s * (%7.4f)", 
                                "h" + String.format("%d@%s%d", iHa+1, Centers[a].Name, a+1),
                                 cf.get(0, 1),
                                 "h" + String.format("%d@%s%d", iHb+1, Centers[b].Name, b+1),
                                 cf.get(1, 1) );
                        
                        //String info2 = String.format("d11 = %.4f, d12 = %.4f, d22 = %.4f",  d11, d12, d22);                        
                        
                        out.printf("%5d\t%-35s\t%-9.5f\t%s\t%n", 
                                iLPO-1,  
                                String.format("(BD) %s, Io = %.4f", 
                                        host_atoms, 
                                        Math.cos(Math.atan2(2*d12, Math.abs(d11-d22)))
                                ),
                                vals[0],
                                info1
                            );
                        out.printf("%5d\t%-35s\t%9.5f\t%s\t%n", 
                                iLPO,
                                "     "+host_atoms+", antibonding (NB)",
                                 vals[1], 
                                info2);
                        
                        L2_sum += vals[0]*vals[0];
                        num_bd++;
                    }
                }
            }
        }

        if ( Ry_print_threshold > 0) {
            out.println();
            if (num_ry_below_thresh > 0) {
                out.printf("Note: %d one-center RY orbitals, each having occupancy below %.2e,%n"+
                           "      were not printed (use the program option %s to change this behavior)%n",
                           num_ry_below_thresh, Ry_print_threshold,
                           options.RyOccPrintThreshold.Name );
                out.printf("      Total occupancy of these RY orbitals is %.5f%n", ry_below_thresh_occs);            
            }
        }
        
        out.println();
        out.println();
        
        
        lo_info.BDperAtomicPair = numBonds();
        printBondMatrix(lo_info.BDperAtomicPair);
        
        
        // (debug-like): print 'atomic valences' if LOs are CLPOs
        if (loNameStr.equals("CLPO"))  {
            for(int a=0; a<nAtoms; a++) {
                out.printf("VAL: \t %2s \t %d %n", Centers[a].Name, lo_info.BDperAtomicPair[a][a]);
            }
        }
        
        out.println();
        out.printf(">> %s occupancy summary >>%n", loNameStr);
        out.printf("       bonding (BD): %12.5f in %4d oribtals%n", BD_sum, num_bd);
        out.printf("  anti-bonding (NB): %12.5f in %4d oribtals%n", NB_sum, num_bd);
        out.printf(" 1c-lone pairs (LP): %12.5f in %4d oribtals%n", LP_sum, num_lp);
        out.printf(" 1c-unoccupied (RY): %12.5f in %4d oribtals%n", RY_sum, num_ry);        
        out.println();

        out.printf("Method        BD+LP....in      NB+RY     BD+NB+LP  BD+NB+LP+RY      trace(D)   Sum[Bd^2+Lp^2]    ||D||^2%n");
        double tmp = D_AHO.normF();
        out.printf(" %5s %12.5f %5d %10.5f %12.5f %12.5f %12.3f   %12.4f  %12.4f", 
                        //bd+lp  nb+ry  b+n+l total  tr(d)  sum_sq  ||d||^^2
                loNameStr,                                
                BD_sum + LP_sum, num_bd + num_lp,
                NB_sum + RY_sum, 
                BD_sum + LP_sum + NB_sum, 
                BD_sum + LP_sum + NB_sum + RY_sum,
                D_AHO.trace(),
                L2_sum,
                tmp*tmp
                );
        out.println();
                
        
     
        return D_AHO;
    }
    
    //--------------------------------------------------------------------------
    private double targ_func1(int[][] hybr2iH) {
        Matrix Dnew = SDS_in_hybrid_basis(hybr2iH, null);
        double tmp = 0;
        for(int a=0; a<nAtoms; a++) {
            for (int ha=0; ha<hybridsOfAtoms[a].nValidHybrids; ha++) {
                int ia = hybr2iH[a][ha];
                int b = hybridsOfAtoms[a].friendAtomIndex[ha];
                tmp += Dnew.get(ia, ia)*Dnew.get(ia, ia);
                if (b != -1) {
                    int hb = hybridsOfAtoms[a].friendHybridIndex[ha];
                    int ib = hybr2iH[b][hb];
                    tmp += Dnew.get(ia, ib)*Dnew.get(ib, ia);
                }
            }
        }
        return tmp;        
    }
    //--------------------------------------------------------------------------
    /**
     * Iteratively optimizes both the hybrid connectivity (2C subset) and 
     * AO-to-hybrids transformation matrices
     */
    private void iterativeHybridOpt(int[][] hybr2iH, boolean allowAll, String iter_comment) {
        Matrix Dnew;
        //int[][] bonds_prev = null;
        double win_prev = 0.0;
        
        int iter = 0;        
        while (iter < nNAOs ) {

            
            Dnew = SDS_in_hybrid_basis(hybr2iH, null);
            //out.printf("PHI = %.12e tr(DNew**2) = %.12e %n", targ_func1(hybr2iH), -1.0);//Dnew.times(Dnew.transpose()).trace());
            //out.printf("tr(DNew**2) = %.12e %n", Dnew.times(Dnew.transpose()).trace());
            
            double win_new = reconnectHybrids(Dnew, hybr2iH, allowAll);    // changes hybridsOfAtoms[a] connectivity
            //out.printf("PHI = %.12e tr(DNew**2) = %.12e %n", targ_func1(hybr2iH), Dnew.times(Dnew.transpose()).trace());

            out.printf("%s iteration %2d: hybrids reconnected, target function = %12.7f%n", 
                    iter_comment, iter+1, win_new);

            if ((iter > 0) && ( win_new < win_prev - 1.0e-10 )) { // minus machiene precision!
                out.printf("WOW!!! new = %.12e < prev = %.12e ; new-prev = %.12e%n", win_new, win_prev, win_new-win_prev );
                out.printf("PHI = %.12e tr(DNew**2) = %.12e %n", targ_func1(hybr2iH), Dnew.times(Dnew.transpose()).trace());
                
            }
            if ((iter > 0) && (Math.abs( win_prev - win_new ) < target_function_conv_thresh )) { // at least one full iteration must be performed!
                out.printf("Done! ");
                break;
            }
            // else:
            win_prev = optimizeHybrids(); 
            //out.printf("PHI = %.12e tr(DNew**2) = %.12e %n", targ_func1(hybr2iH), Dnew.times(Dnew.transpose()).trace());
            out.printf("%s iteration %2d: hybrids optimized, target function = %12.7f%n", 
                    iter_comment, iter+1, win_prev);
            
/*            for (int a=0; a<nAtoms; a++)
                hybridsOfAtoms[a].backup_bonding();*/
            
            iter++;
        }
        out.printf("(in %d iterations)%n", iter);     
        out.println();
    }
    //--------------------------------------------------------------------------
    /**
     * @returns user-friendly names for LHOs/AHOs
     */
    private String[] create_Hybrid_Labels(String prefix, int[][] addr2iH) {
        String[] result = new String[nNAOs];
        for (int a=0; a<nAtoms; a++) {            
            for (int ha=0; ha<hybridsOfAtoms[a].NAO_indices.length; ha++) {
                int iH =  addr2iH[a][ ha ];
                result[ iH ] = String.format("%s%d:%s%d", Centers[a].Name, a+1, prefix, iH+1);
            }
        }        
        return result;
    }
    //--------------------------------------------------------------------------

    /**
     * Creates density-norm optimal localized orbitals.
     * Input:
     * @param SDS_NAO -- the density matrix in (orthonormal) NAO basis
     * @param NAO -- array indicating to which atom each NAO belongs
     * @para Atoms -- the total number of atoms
     */    
    
    private double target_function_conv_thresh = 1e-5;
    private int Opt_Max_Iter = 10000;
    
    public void createCLPOs(
            Matrix SDS_NAO, BasisFunction[] NAOs, AtomicCenter[] Centers,
            ono_options options) 
    {
        this.options = options;        
        target_function_conv_thresh = options.HybrOptOccConvThresh.get_double();
        Opt_Max_Iter = options.HybrOptMaxIter.get_int();
        
        nAtoms = Centers.length;
        this.Centers = Centers;
        
        this.SDS_NAO = SDS_NAO;
        this.nAtoms = nAtoms;
        this.nNAOs = NAOs.length;

        out.printf("Creating  LPOs (Localized Property-optimized Orbitals)%n%n%n");
        
        NAOs_at_center = get_NAOs_at_centers(NAOs);        

        hybridsOfAtoms = new AtomicHybrids[nAtoms]; // just alloc mem...
        for (int a=0; a<nAtoms; a++) { // ...and create objects for holding atomic hybrids
            hybridsOfAtoms[a] = new AtomicHybrids(NAOs_at_center[a]);
        }
        
        // create D[a][b<=a] submatrices
        DAB = diatomic_Submatrices();

        // prepare suitable hybrid guess and save it into hybridsOfAtoms[a].U matrices        
        out.println("Creating initial guess...");
        
        CS_Guess();        
        
        int[][] hybr2iH = getGlobalHybridIndices();//null); // must be called after CS_Guess();
        
        out.printf("Optimizing LPOs...%n%n");
        opt_Lewis_mode = false;
        iterativeHybridOpt(hybr2iH, true, "LPO");        

        /*this.LPO_to_AHO = new Matrix(nNAOs, nNAOs);
        this.NAO_to_AHO = new Matrix(nNAOs, nNAOs);
        this.LPO_labels = new String[nNAOs];
        createLOs("LPO ", new Matrix[]{ LPO_to_AHO  }, new Matrix[] { NAO_to_AHO }, LPO_labels);*/
        
        this.LPO_descript = new LOdescription(nNAOs);
        createLOs("LPO ",  LPO_descript);
        
        // Note: hybrid numbers/ordering could (?) have been modified/updated by iterativeHybridOpt
        // => re-create hybr2iH for correct labeling and further work
        // Subnote: createLOs() also calls getGlobalHybridIndices() internally
        int[][][] p_hybr2iH = new int[][][]{null};
        hybr2iH = getGlobalHybridIndices( );//p_hybr2iH ); 
        //atomic_AHOs_to_globalHybridIds = p_hybr2iH[0];
        
                
        out.println();
        //----------------------------------------------------------------------        
        printout.PrintStars();
        
        out.printf("Creating CLPOs (Chemist's Localized Property-optimized Orbitals)%n%n%n");
        out.println("Optimizing CLPOs...");
        opt_Lewis_mode = true;        
        iterativeHybridOpt(hybr2iH, false, "CLPO");

        // Note: hybrid numbers/ordering could (?) have been modified/updated by iterativeHybridOpt
        // => re-create hybr2iH for correct labeling and further work
        // Subnote: createLOs() also calls getGlobalHybridIndices() internally
        hybr2iH = getGlobalHybridIndices( );//p_hybr2iH ); // todo: refactor & remove
        //atomic_LHOs_to_globalHybridIds = p_hybr2iH[0];
        
        
        /*this.CLPO_to_LHO = new Matrix(nNAOs, nNAOs);
        this.NAO_to_LHO = new Matrix(nNAOs, nNAOs);
        this.CLPO_labels = new String[nNAOs];
        createLOs("CLPO", new Matrix[]{ CLPO_to_LHO  }, new Matrix[] { NAO_to_LHO }, CLPO_labels);*/
        this.CLPO_descript = new LOdescription(nNAOs);
        createLOs("CLPO", CLPO_descript);
                
        //LHO_labels = 
        CLPO_descript.Hybrid_labels = create_Hybrid_Labels("LHO", hybr2iH);
        
    }
    //--------------------------------------------------------------------------
    /**
     * Performs depth-first search to identify all molecules in the system
     * based on connectivity information provided by the localized orbitals @param LOs;
     * @param atomicCharges is used to compute the total charges of each framgent;
     * @returns the list of framgent ids (1-based) to which each atom belongs
     * Note: uses BDperAtomicPair[][] from LOs 
     * 
     */
    public int[] LOconnectivity(LOdescription LOs, double[] atomicCharges) {
        int[] fragmentOfAtom = new int[nAtoms]; // to which framgent (1-based) each atom belongs; value -1 means 'unknown'
        int[] valences = new int[nAtoms];
        int[][] neighs = new int[nAtoms][nAtoms]; // [a] = list of the atoms neighbouring to a; 
        int[] nNeighs = new int[nAtoms]; // [a] = number of valid elements in neighs[a] 

        // fill in the neighs[][], nNeighs[] and valences[] arrays based on adjacency matrix from LOs.BDperAtomicPair[][]
        for(int a=0; a<nAtoms; a++) {
            for(int b=a+1; b<nAtoms; b++) {
                if (LOs.BDperAtomicPair[a][b] > 0) {
                    valences[a]++;
                    valences[b]++;
                    neighs[a][nNeighs[a]++] = b; // Note that use adjacency matrix rather than raw bond list, so nNeighs is incremented not more than once per atomic pair!
                    neighs[b][nNeighs[b]++] = a;                    
                }
            }
        }
        
        // Run depth-first search to find all connected subgraphs

        int[] yetToVisit = new int[nAtoms];
        int lenYetToVisit = 0;
        for(int i=0; i<nAtoms; i++) 
            fragmentOfAtom[i] = -1; // 'not visited yet'
        int currentFragmentId = 0;  // the first fragment will be numbered as '1' (this is 
                                    // handy since currentFragmentId equals the total number of fragments when the search ends)

        for(int i=0; i<nAtoms; i++) {
            if (fragmentOfAtom[i] != -1)
                continue; //already visited            
            yetToVisit[lenYetToVisit++] = i;
            fragmentOfAtom[i] = ++currentFragmentId;
            
            while(lenYetToVisit > 0) {
                // pop the last element ...
                int c = yetToVisit[--lenYetToVisit];
                // ...and push all of its yet-not-visited neighbours instead (depth-first search)
                for(int j=0; j<nNeighs[c]; j++) {
                    int a2 = neighs[c][j];
                    if (fragmentOfAtom[a2] == -1) {
                        yetToVisit[lenYetToVisit++] = a2;
                        fragmentOfAtom[ a2 ] = currentFragmentId; // mark as visited
                    }
                }                
            }
        }
        
        // recompute NPA charges if atomicCharges == null -- TODO
        
        // print results

        out.printf("There are %d molecule(s) in the system%n", currentFragmentId);
        out.println("(the 'molecule' is defined as the set of atoms linked with BD orbitals)");
        
        double[] framgentCharges = new double[currentFragmentId];
        int[] nAtomsInMol = new int[currentFragmentId];
        for(int i=0; i<nAtoms; i++) {
            int f = fragmentOfAtom[i] - 1;
            framgentCharges[f] += atomicCharges[i];
            nAtomsInMol[f] ++;
            /*out.printf("%5s %3d %3d%n", 
                    String.format("%s%d", Centers[i].Name,i+1),
                    f+1, valences[i]);        */
        }
        
        out.println();
        //out.println("Molecular properties:");
        //out.println("MolId TotCharge  nAtoms  ListOfAtoms");
        out.println("Molecule   TotalNPA    NumberOf   ListOf ");
        out.println("id         charge      atoms      atoms...");
        for (int f=0; f<currentFragmentId; f++) {
            out.printf("%4d    %+11.5f  %7d     ", f+1, framgentCharges[f], nAtomsInMol[f]);
            for(int i=0; i<nAtoms; i++) {
                if ( fragmentOfAtom[i] == f+1 )
                    out.printf(" %s%d", Centers[i].Name, i+1);
                
            }
            out.println();
        }
        out.println();
        out.println("Note: The total NPA charge is computed as the sum of NPA charge of atoms belonging to each molecule");
        //out.println();
     
        return fragmentOfAtom;
        
    }
    
}
