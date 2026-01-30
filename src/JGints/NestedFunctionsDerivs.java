package JGints;

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
 */

import Polynom3D.*;

/**
 * An auxiliary class for calculating partial derivatives of any order of the
 * functions used by overlap integral evaluators
 *
 * Last.Rev.: 26.Oct.2014
 * 
 */
public class NestedFunctionsDerivs {
    final int max_deriv = 3;
    final int Cnk_max = 5; // should be not less than max_deriv

    public static double[][] Quick_CNK = new double[0][0]; // the values[n][k] of n!/k!/(n-k)!
    private double[] zetaDiffCi = new double[0]; // value[i] = (-1)^i/2^i * (2*i-1)!! / (2*i-1), valid for i>0
    
    public static double[][] Ank = new double[0][0]; // // A[n][mu] = n!/(n-mu)!
    //--------------------------------------------------------------------------
    /** Compute A[n][mu] = n!/(n-mu)!
     *
     */
    public static double[][] ensure_Ank_enough(int max_n) {
        if ((Ank.length-1) < max_n) {
            double[][] new_Ank = new double[max_n+1][];
            // just copy old values
            for (int i=0; i<Ank.length; ++i) {
                new_Ank[i] = new double[i+1];
                System.arraycopy(Ank[i], 0, new_Ank[i], 0, Ank[i].length);
            }
            // and compute new elements
            new_Ank[0] = new double[]{1.0}; // this is LEGAL in any case!
            int n_from = 1;
            if (Ank.length>0) n_from = Ank.length;
            for (int n=n_from; n<=max_n; ++n) {
                new_Ank[n] = new double[n+1];
                new_Ank[n][n] = new_Ank[n-1][n-1] * n;
                for (int mu=n-1; mu>=0; --mu)
                    new_Ank[n][mu] = new_Ank[n][mu+1]/(n-mu);
            }
            Ank = new_Ank; // update the array
        }
        return Ank;
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // A constructor
    public NestedFunctionsDerivs() {
        // fill in the array Quick_CNK
        ensure_Cnk_enough(Cnk_max);
        ensure_zetaDiffCi_enough(max_deriv);
    }
    //--------------------------------------------------------------------------
    // redirects to Polynom3D.ensure_Cnk_enough
    private void ensure_Cnk_enough(int Cnk_max) {
        Quick_CNK = Polynom3D.ensure_Cnk_enough(Cnk_max);
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether zetaDiffCi(n) with n up to @param max_n are available, and if not,
     * computes the missing values as
     *
     * zetaDiffCi[i] = (-1)^i/2^i * (2*i-1)!! / (2*i-1) = (-1)^i/2^i * (2*i-3)!! (for i>1)
     */
    private void ensure_zetaDiffCi_enough(int max_n) {
        if ((zetaDiffCi.length-1) < max_n) {
            double[] new_ci = new double[max_n+1];
            // just copy old values
            System.arraycopy(zetaDiffCi, 0, new_ci, 0, zetaDiffCi.length);
            
            // and compute new elements
            for (int i=zetaDiffCi.length; i<=max_n; ++i) {
                if (i==0)
                    new_ci[0] = -1.0;
                else
                    new_ci[i] = -new_ci[i-1]/2.0*(2*i-3);
            }

            zetaDiffCi = new_ci; // update the array
        }
    }


    //--------------------------------------------------------------------------
    /**
     * @returns derivatives
     *                  d^(i+j)  (                     )
     *  result[i][j] = --------- ( F(x,y) * G(x,y) * k )
     *                 dx^i dy^j (                     )
     * where
     *                    d^(i+j)
     * @param fd[i][j] = ---------- F(x,y)
     *                   dx^i dy^j
     *
     *                    d^(i+j)
     * @param fg[i][j] = ---------- G(x,y)
     *                   dx^i dy^j
     * 
     * and the result is valid for i<=@param maxI, j<=@param maxJ
     * 
     */
    public double[][] productDerivs(double[][] fd, double[][] gd, double k, int maxI, int maxJ) {
        // prepare an array of Quick_CNK
        if (maxI > maxJ) ensure_Cnk_enough(maxI); else ensure_Cnk_enough(maxJ);

        // do the job
        double[][] result = new double[maxI+1][maxJ+1];
        result[0][0] = fd[0][0] * gd[0][0]; // product value
        for(int i=0; i<=maxI; ++i)
            for (int j=0; j<=maxJ; ++j) {
                // now compute the (i,j)-the derivative of F*G
                result[i][j] = 0;
                for(int mu=0; mu<=i; ++mu)
                    for (int nu=0; nu<=j; ++nu)
                        result[i][j] += Quick_CNK[i][mu]*Quick_CNK[j][nu] * fd[mu][nu] * gd[i-mu][j-nu];
                result[i][j] *= k;
            }

        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * @returns the values of
     *
     *                       d^(i+j)            zeta1 * zab
     * result[i][j] = ------------------- ----------------------- = zab * zetaEff^(-1/2-i-j) * ((2*j-1)*zeta1 - 2*i*zeta2)* C[i+j]
     *                 dzeta1^i dzeta2^j   (zeta1 + zeta2)^(1/2)
     *
     * up to i<=@param n and j<=@param n (zetaEff = zeta1 + zeta2)
     *
     * where C[i] = (-1)^i/2^i * (2*i-1)!! / (2*i-1)
     * 
     * result[0][0] is filled with the function value
     *
     * if @param transposeResult is true, i and j indices of the result are swapped
     */    
    public double[][] abDerivs(double zeta1, double zeta2, double zab, int maxI, int maxJ, boolean transposeResult) {
        double zetaEff = zeta1 + zeta2;
        // prepare the coefficients in zetaDiffCi[]
        ensure_zetaDiffCi_enough(maxI+maxJ);
        // prepare the powers of zetaEff^(-1/2-i-j)
        double[] zefaEff_pwr = new double[maxI+maxJ+1]; // value[i] = zetaEff^(-1/2-i)
        zefaEff_pwr[0] = 1/Math.sqrt(zetaEff);
        for (int i=1; i<zefaEff_pwr.length; ++i)
            zefaEff_pwr[i] = zefaEff_pwr[i-1]/zetaEff;
        // create a result
        double[][] result ;
        if (transposeResult) {
            result = new double[maxJ+1][maxI+1];
            // 'transposed' fill version
            for (int i=0; i<=maxI; ++i)
                for (int j=0; j<=maxJ; ++j)
                    result[j][i] = zefaEff_pwr[i+j] * ((2*j-1)*zeta1 - 2*i*zeta2) * zab * zetaDiffCi[i+j];
        } else {
            // 'normal' fill version
            result = new double[maxI+1][maxJ+1];
            for (int i=0; i<=maxI; ++i)
                for (int j=0; j<=maxJ; ++j)
                    result[i][j] = zefaEff_pwr[i+j] * ((2*j-1)*zeta1 - 2*i*zeta2) * zab * zetaDiffCi[i+j];
        }

        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * @returns the values of
     *
     *                  d^(i+j)
     * result[i][j] = ----------- exp( F(x,y)  )
     *                 dx^i dy^j
     * up to i<=@param maxI and j<=@param maxJ
     * The only input is 
     * 
     *               d^(i+j)
     * @param fd = ----------- F(x,y)
     *              dx^i dy^j
     */
    public double[][] exponentialDerivs(double[][] fd, int maxI, int maxJ) /*throws Exception */{
        // allocate mem
        double[][] result = new double[maxI+1][maxJ+1];

        // calculate the function value - for result[0][0]-th element and further calculations
        double fValue = Math.exp(fd[0][0]); // the value of the function
        
        result[0][0] = fValue;

        // prepare an array of Quick_CNK
        if (maxI > maxJ) ensure_Cnk_enough(maxI); else ensure_Cnk_enough(maxJ);

        // generate result[0][0] -> result[1...maxI][0]:
        // calculate result[i][0] via result[0...i-1][0] by getting the i-th derivative  of f=exp(F) as
        // d^i f / d1^i = d^(i-1) (df/d1) / d1^(i-1) = d^(i-1) (exp(F)*dF/d1) / d1^(i-1)
        for(int i=1; i<=maxI; ++i) {
            result[i][0] = 0;
            for(int mu=0; mu<=i-1; ++mu) 
                    result[i][0] += Quick_CNK[i-1][mu]* result[mu][0] * fd[i-mu][0]; // generates result[i+1][0] from result[0...i][0]
        }
         
        // Now use the fact that we know all result[0...maxI][0] to calculate result[0...maxI][j], j=1...maxJ:
        // d^i d^(j+1)f/d2^j /d1^i = d^j (d^i(df/d2)/d1^i) /d2^j = d^j d^i(f * F'2 )/d1^i /d2^j, where j=0...maxJ-1, and i=0...maxI
        for(int i=0; i<=maxI; ++i) {
           // and now - fill in result[i][1...maxJ]===result[i][j+1] with j=0...maxJ-1
           for (int j=0; j<=maxJ-1; ++j) {
               // now compute the (i,j)-the derivative of A*B where A = exp(F) and B = F'2
               result[i][j+1] = 0;
               for(int mu=0; mu<=i; ++mu)
                   for (int nu=0; nu<=j; ++nu)
                       result[i][j+1] += Quick_CNK[i][mu]*Quick_CNK[j][nu] * result[mu][nu] * fd[i-mu][j-nu+1]; // generates result[i][j+1] from result[0...i]
           }
        }

         // compare results :)
         /*double diff = 0;
         for(int i=0; i<=maxI; ++i)
            for (int j=0; j<=maxJ; ++j) diff += (x[i][j] - result[i][j])*(x[i][j] - result[i][j]);
         diff = Math.sqrt(diff / maxI / maxJ);
         if (diff > 1E-15) {
             Jama.Matrix M = new  Jama.Matrix(x);
             M.print(15, 10);
             M = new  Jama.Matrix(result);
             M.print(15, 10);
        }
          * 
          */
         return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Computes
     *
     *                 d^(i+j)   (   zeta1 * zeta2        )
     * result[i][j] = ---------- ( - ------------- * DZ^2 ) = productDerivs(abDerivs(...), abDerivs(...)) _*
     *                dx^i dy^j  (   zeta1 + zeta2        )
     * 
     * for i<=@param maxI and j<=@param maxJ
     *
     *  _* indeed, the function to be differentiated is just a product of a(zeta1,zeta2) and -b(zeta1,zeta2)
     */
    public double[][] prefactorExponentDerivs(int maxI, int maxJ, double DZ, double zeta1, double zeta2) /*throws Exception */{

            //double[][] y = abDerivs(zeta2,zeta1,-DZ,maxJ,maxI, true);
            //double[][] r = productDerivs(abDerivs(zeta1,zeta2,DZ,maxI,maxJ,false), y, 1.0, maxI, maxJ);

            // new version
            double[] aux = new double[maxI+1 + maxJ+1];
            double minusZetaEffInv = -1/(zeta1+zeta2);
            aux[0] = 1.0;
            aux[1] = minusZetaEffInv; // aux[1] is ALWAYS valid!
            for (int i=2; i<aux.length; ++i)
                aux[i] = aux[i-1] * minusZetaEffInv * (i-1);

            double[][] result = new double[maxI+1][maxJ+1];
            double DZ2 = DZ*DZ;
            double zeta12 = zeta1*zeta1;
            double zeta22 = zeta2*zeta2;
            double tmp;

            // evaluate funtion value
            result[0][0] = zeta1*zeta2 * minusZetaEffInv * DZ2;

            // Evaluate derivatives as:
            //        (-1)^(i+j)                  (     i*zeta2^2 + j*zeta1^2          zeta1+zeta2  )
            //  -------------------- * (i+j-1)! * ( -  ----------------------- + i*j* ------------- )
            //  (zeta1+zeta2)^(i+j)               (         zeta1 + zeta2                i+j-1      )
            //                                                                   ~~~~~~~~~~~~~~~~~~~
            //                                                                  this term is zero, if i+j==1            

            for (int i=0; i<=maxI; ++i)
                for (int j=0; j<=maxJ; ++j)
                    if (i+j>0) { // do NOT overwrite result[0][0]
                        tmp = (i*zeta22 + j*zeta12)*minusZetaEffInv;
                        if (i>0 && j>0)
                            tmp += i*j*(zeta1+zeta2)/(i+j-1);
                        result[i][j] = aux[i+j] * tmp * (-DZ2);
                    }
            /*Jama.Matrix m1 = new Jama.Matrix(r);
            m1.print(13, 7);
            Jama.Matrix m2 = new Jama.Matrix(result);
            m2.print(13, 7);*/

            return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Computes coefficients for
     *
     *                 d^j  (     1   )                      1
     *                ----- ( ------- ) = result[n][j] * ---------
     *                dx^j  ( x^(n/2) )                   x^(n/2+j)
     *
     */
    public static double[][] xSqrtNInvDerivs(int maxN, int maxJ){
        double[][] result = new double[maxN+1][maxJ+1];
        for (int n=0; n<=maxN; ++n) {
            result[n][0] = 1.0;
            for (int j=1; j<=maxJ; ++j)
                result[n][j] = -result[n][j-1]*(n/2.0+j-1);
        }
        return result;
    }
    //--------------------------------------------------------------------------

    /**
     * And the most important method which computes
     *
     *                 d^(i+j)   (                     )
     * result[i][j] = ---------- ( W (z0(zeta1, zeta2) )
     *                dx^i dy^j  (                     )
     *
     * for i<=@param maxI and j<=@param maxJ
     *
     *  
     */
    /*
    public double[][] integral_Derivs(int maxI, int maxJ, double[] dW, double[][] dz) throws Exception {
        if (maxI > max_deriv || maxJ > max_deriv)
            throw new Exception(
                    String.format("integral_Derivs() is not designed to work with maxI = %d and maxJ = %d ",maxI,maxJ));

        double[][] result = new double[maxI+1][maxJ+1];
        result[0][0] = dW[0];

        if (maxI >= 0 && maxJ >= 1) result[0][1] = dW[1]*dz[0][1];
        if (maxI >= 0 && maxJ >= 2) result[0][2] = dW[2]*dz[0][1]*dz[0][1]+dW[1]*dz[0][2];
        if (maxI >= 0 && maxJ >= 3) result[0][3] = dW[1]*dz[0][3]+(3*dW[2]*dz[0][2]+dW[3]*dz[0][1]*dz[0][1])*dz[0][1];
        if (maxI >= 1 && maxJ >= 0) result[1][0] = dW[1]*dz[1][0];
        if (maxI >= 1 && maxJ >= 1) result[1][1] = dW[2]*dz[0][1]*dz[1][0]+dW[1]*dz[1][1];
        if (maxI >= 1 && maxJ >= 2) result[1][2] = dW[2]*dz[0][2]*dz[1][0]+dW[1]*dz[1][2]+(2*dW[2]*dz[1][1]+dW[3]*dz[1][0]*dz[0][1])*dz[0][1];
        if (maxI >= 1 && maxJ >= 3) result[1][3] = 3*dW[2]*dz[0][2]*dz[1][1]+dW[1]*dz[1][3]+dW[2]*dz[0][3]*dz[1][0]+(3*dW[3]*dz[1][0]*dz[0][2]+3*dW[2]*dz[1][2]+(3*dW[3]*dz[1][1]+dW[4]*dz[1][0]*dz[0][1])*dz[0][1])*dz[0][1];
        if (maxI >= 2 && maxJ >= 0) result[2][0] = dW[2]*dz[1][0]*dz[1][0]+dW[1]*dz[2][0];
        if (maxI >= 2 && maxJ >= 1) result[2][1] = 2*dW[2]*dz[1][0]*dz[1][1]+(dW[3]*dz[1][0]*dz[1][0]+dW[2]*dz[2][0])*dz[0][1]+dW[1]*dz[2][1];
        if (maxI >= 2 && maxJ >= 2) result[2][2] = 2*dW[2]*dz[1][0]*dz[1][2]+2*dW[2]*dz[1][1]*dz[1][1]+dW[1]*dz[2][2]+(dW[3]*dz[1][0]*dz[1][0]+dW[2]*dz[2][0])*dz[0][2]+(4*dW[3]*dz[1][0]*dz[1][1]+(dW[4]*dz[1][0]*dz[1][0]+dW[3]*dz[2][0])*dz[0][1])*dz[0][1]+2*dW[2]*dz[0][1]*dz[2][1];
        if (maxI >= 2 && maxJ >= 3) result[2][3] = 2*dW[2]*dz[1][0]*dz[1][3]+dW[1]*dz[2][3]+(dW[3]*dz[1][0]*dz[1][0]+dW[2]*dz[2][0])*dz[0][3]+6*dW[2]*dz[1][1]*dz[1][2]+6*dW[3]*dz[0][2]*dz[1][0]*dz[1][1]+(3*dW[2]*dz[2][2]+6*dW[3]*dz[1][0]*dz[1][2]+6*dW[3]*dz[1][1]*dz[1][1]+(3*dW[3]*dz[2][0]+3*dW[4]*dz[1][0]*dz[1][0])*dz[0][2]+(6*dW[4]*dz[1][0]*dz[1][1]+(dW[4]*dz[2][0]+dW[5]*dz[1][0]*dz[1][0])*dz[0][1])*dz[0][1])*dz[0][1]+(3*dW[3]*dz[0][1]*dz[0][1]+3*dW[2]*dz[0][2])*dz[2][1];
        if (maxI >= 3 && maxJ >= 0) result[3][0] = (3*dW[2]*dz[2][0]+dW[3]*dz[1][0]*dz[1][0])*dz[1][0]+dW[1]*dz[3][0];
        if (maxI >= 3 && maxJ >= 1) result[3][1] = 3*dW[3]*dz[1][0]*dz[1][0]*dz[1][1]+3*dW[2]*dz[1][1]*dz[2][0]+(3*dW[3]*dz[2][0]+dW[4]*dz[1][0]*dz[1][0])*dz[1][0]*dz[0][1]+3*dW[2]*dz[1][0]*dz[2][1]+dW[1]*dz[3][1]+dW[2]*dz[0][1]*dz[3][0];
        if (maxI >= 3 && maxJ >= 2) result[3][2] = dW[1]*dz[3][2]+(3*dW[2]*dz[2][2]+6*dW[3]*dz[1][1]*dz[1][1])*dz[1][0]+(3*dW[2]*dz[2][0]+3*dW[3]*dz[1][0]*dz[1][0])*dz[1][2]+(3*dW[3]*dz[2][0]+dW[4]*dz[1][0]*dz[1][0])*dz[1][0]*dz[0][2]+(6*dW[3]*dz[1][1]*dz[2][0]+6*dW[4]*dz[1][0]*dz[1][0]*dz[1][1]+(3*dW[4]*dz[2][0]+dW[5]*dz[1][0]*dz[1][0])*dz[1][0]*dz[0][1])*dz[0][1]+(6*dW[2]*dz[1][1]+6*dW[3]*dz[0][1]*dz[1][0])*dz[2][1]+2*dW[2]*dz[0][1]*dz[3][1]+(dW[3]*dz[0][1]*dz[0][1]+dW[2]*dz[0][2])*dz[3][0];
        if (maxI >= 3 && maxJ >= 3) result[3][3] = dW[1]*dz[3][3]+(9*dW[2]*dz[2][2]+6*dW[3]*dz[1][1]*dz[1][1])*dz[1][1]+3*dW[2]*dz[1][3]*dz[2][0]+(3*dW[2]*dz[2][3]+3*dW[3]*dz[1][3]*dz[1][0])*dz[1][0]+(3*dW[3]*dz[2][0]+dW[4]*dz[1][0]*dz[1][0])*dz[1][0]*dz[0][3]+18*dW[3]*dz[1][0]*dz[1][1]*dz[1][2]+(9*dW[3]*dz[1][1]*dz[2][0]+9*dW[4]*dz[1][0]*dz[1][0]*dz[1][1])*dz[0][2]+(3*dW[2]*dz[3][2]+(18*dW[4]*dz[1][1]*dz[1][1]+9*dW[3]*dz[2][2])*dz[1][0]+(9*dW[3]*dz[2][0]+9*dW[4]*dz[1][0]*dz[1][0])*dz[1][2]+(9*dW[4]*dz[2][0]+3*dW[5]*dz[1][0]*dz[1][0])*dz[1][0]*dz[0][2]+(9*dW[5]*dz[1][0]*dz[1][0]*dz[1][1]+9*dW[4]*dz[1][1]*dz[2][0]+(3*dW[5]*dz[2][0]+dW[6]*dz[1][0]*dz[1][0])*dz[1][0]*dz[0][1])*dz[0][1])*dz[0][1]+(9*dW[3]*dz[0][2]*dz[1][0]+9*dW[2]*dz[1][2]+(18*dW[3]*dz[1][1]+9*dW[4]*dz[1][0]*dz[0][1])*dz[0][1])*dz[2][1]+(3*dW[2]*dz[0][2]+3*dW[3]*dz[0][1]*dz[0][1])*dz[3][1]+(dW[2]*dz[0][3]+(3*dW[3]*dz[0][2]+dW[4]*dz[0][1]*dz[0][1])*dz[0][1])*dz[3][0];

        return result;
    }
     * 
     */
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

}
