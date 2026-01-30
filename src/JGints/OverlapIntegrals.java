package JGints;

/**
 * A file containing OverlapIntegrals class which is used for calculation
 * of overlap integrals between the GTO-type AO basis functions
 *
 * Version: 26.Oct.2013
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
 */

import Polynom3D.*;

//==============================================================================
/**
 * A class to calculate overlap integrals for exp(-zeta1*(r1-R01)^2)*Y_lm(r1-R01) * exp(-zeta2*(r2-R02)^2)*Y_lm(r2-R02)
 * @author timn
 */
public class OverlapIntegrals {
    private final boolean debugPrint = false;
    // gamma half-int (n) = Gamma (n+1/2)
    static double Gamma_hi(int n) {
        double result = Math.sqrt(Math.PI);
        for (int i=1; i<=n; i++)
            result *= (1/2.0+(i-1));
        return result;
    }
    //--------------------------------------------------------------------------
    // calculates
    //    / +inf
    //   | t^n · exp(-alpha·t^2)·dt
    //  / -inf
    //
    public static double primitive_int_1D(int n, double alpha) {
        if ((n%2)==1) return 0.0;
        // now for n == 2*k
        return Math.pow(alpha, -1/2.0 - n/2.0) * Gamma_hi(n/2);
    }
    //--------------------------------------------------------------------------
    // calculates
    //    / +inf
    //   | r^n · exp(-alpha·r^2)·dr
    //  / 0
    //
    public static double primitive_int_1D_Sphr(int n, double alpha) {
        if ((n%2)==0) return primitive_int_1D(n,alpha) / 2.0;
        // else:
        // reduce integral to:
        //    / +inf
        //   | r^(n-1) · exp(-alpha·r^2)·d(r^2)*0.5
        //  / 0
        // where n-1 = 2k, so that result will be
        //    / +inf
        //   | (alpha·t)^k · exp(-alpha·t)·d(alpha·t) * 0.5/alpha^(k+1)  = 0.5/alpha * k! / alpha^k
        //  / 0
        double result = 0.5/alpha;
        for (int i=1; i<=((n-1)/2); i++) {
            result *= i/alpha;
        }
        return result;
    }
    //--------------------------------------------------------------------------
    private final int L_MAX = 10; // should be not less than anywhere else in this project!
    private final int pwr_MAX = 2*L_MAX;
    // an array for quick access to C_nk values
    private double[][] Quick_CNK = new double[pwr_MAX+1][pwr_MAX+1];
    //--------------------------------------------------------------------------
    // A constructor
    public OverlapIntegrals() {
        // fill in the array Quick_CNK - it is used by BS_BS_Overlap
        for (int i=0; i<=pwr_MAX; i++)
            for (int j=0; j<=i; j++)
                Quick_CNK[i][j] = Polynom3D.C_nk(i, j);
    }
    //--------------------------------------------------------------------------
    /** Calculates an overlap between primitive basis functions
     * each of the form (Y_lm(r-R0) * exp(-zeta*(r-R0)^2), possibly with different R0s)
     */
    public double BS_BS_Overlap(double zeta1, Polynom3D YLM1, double[] R01, double zeta2, Polynom3D YLM2, double[] R02) {
        // use the gaussian product theorem
        double zeta_eff = zeta1 + zeta2;
        double[] R0eff = new double[3]; // new common center
        for (int mu=0; mu<3; mu++)
            R0eff[mu] = (zeta1*R01[mu] + zeta2*R02[mu]) / zeta_eff;
        double Prefactor = 0;
        for (int mu=0; mu<3; mu++)
            Prefactor += (R01[mu] - R02[mu])*(R01[mu] - R02[mu]);
        Prefactor = Math.exp(-zeta1*zeta2*Prefactor / zeta_eff);
        // now the angular part
        //YLM1.MultiplyBy(YLM2);
        double result = 0;
        // utilize the feature:
        //
        //    /                                                      / +inf
        //   |(x-X)^nx·(y-Y)^ny·(z-Z)^nz·exp(-alpha·(r-R)^2)·d3r  = | t^nx·exp(-alpha·t^2)·dt  · (...) · (...)
        //  /                                                      /
        //                                                          -inf
        // find maximum power of (r[mu]-Reff[mu])^n
        int nmax = 0;
        int ncurr;
        for (int i=0; i<YLM1.coefs.length; i++)
            for (int j=0; j<YLM2.coefs.length; j++)
                for (int mu=0; mu<3; mu++) {
                    ncurr = YLM1.powers[i][mu] + YLM2.powers[j][mu];
                    if (ncurr > nmax)  nmax = ncurr;
                }
        // precalculate 1D integrals: int( (t-t0)^n * exp(-zeta_eff*(t-t0)^2), t=-infinity..+infinity )
        // (we leaveall zero terms for simplicity)
        double[] Quick_Ints = new double[nmax+1];
        for (int i=0; i<=nmax; i++)
            Quick_Ints[i]= primitive_int_1D(i, zeta_eff);
        // prepare data for quick expansion of a polynom (t+t1)^n1 * (t+t2)^n2
        // where t1 = R0eff[mu]-R1[mu] and t2 = R0eff[mu]-R1[mu]
        // Note:  R0eff[mu]-R1[mu] = zeta2*(R02[mu]-R01[mu]))/zeta_eff
        //        R0eff[mu]-R2[mu] = zeta1*(R01[mu]-R02[mu]))/zeta_eff
        double[][] QuickExpans_1 = new double[nmax+1][3]; // powers of (R0eff[mu]-R1[mu])
        double[][] QuickExpans_2 = new double[nmax+1][3];
        double[][] cur_pow = new double[2][3]; // first index: t1 or t2; secod index: x,y,z
        for (int mu=0; mu<3; mu++) {
            cur_pow[0][mu] = 1;
            cur_pow[1][mu] = 1;
        }
        // now different powers of (R0eff[mu]-R1[mu]) and (R0eff[mu]-R2[mu])
        // 0-th power:
        for (int mu=0; mu<3; mu++) {
                QuickExpans_1[0][mu] = 1;
                QuickExpans_2[0][mu] = 1;
        }
        // higher powers
        for (int pwr = 1; pwr <= nmax; pwr++)  // powers of (R0eff[mu]-R1[mu]) (and (R0eff[mu]-R2[mu]))
            for (int mu=0; mu<3; mu++) {
                cur_pow[0][mu] *= (R0eff[mu]-R01[mu]);
                cur_pow[1][mu] *= (R0eff[mu]-R02[mu]);
                QuickExpans_1[pwr][mu] = cur_pow[0][mu];
                QuickExpans_2[pwr][mu] = cur_pow[1][mu];
            }
        // an auxiliary variable
        double Quick_C1;
        //
        if (debugPrint) {
            System.out.println("I'm about to integrate:");
            YLM1.Print();
            YLM2.Print();
        }
        // Now use them for integration
        double tmp, term;
        for (int i=0; i<YLM1.coefs.length; i++) {
            for (int j=0; j<YLM2.coefs.length; j++) {
                // evaluate one (i,j)-th term (X*Y*Z) of
                // INT [(r[mu]-R0eff[mu] + (R0eff[mu]-R1[mu])) ^ YLM1.powers[i][mu] *
                //      * (r[mu]-R0eff[mu] + (R0eff[mu]-R2[mu])) ^ YLM2.powers[j][mu] * exp(-zeta_eff*(r[mu]-R0eff[mu])^2) ]
                // term = (...)_x * (...)_y * (...)_z
                // INTEGRAL = sum( terms )
                term = YLM1.coefs[i]*YLM2.coefs[j];
                for (int mu=0; mu<3; mu++) {
                    // mu-th component of (i,j)-th term
                    tmp = 0;
                    for (int pwr1 = 0; pwr1 <= YLM1.powers[i][mu]; pwr1++) {  // power of "+t1"
                        Quick_C1 = Quick_CNK[YLM1.powers[i][mu]][pwr1] * QuickExpans_1[pwr1][mu];
                        //
                        for (int pwr2 = 0; pwr2 <= YLM2.powers[j][mu]; pwr2++)  // power of "+t2"
                            tmp += Quick_C1 *
                                   Quick_CNK[YLM2.powers[j][mu]][pwr2] * QuickExpans_2[pwr2][mu] *
                                   Quick_Ints[(YLM1.powers[i][mu]-pwr1)+(YLM2.powers[j][mu]-pwr2)];
                    }
                    term *= tmp;
                    if (tmp==0) break;
                }
                result += term;
            }
        }
        if (debugPrint) {
            System.out.println("_result = "+result);
            System.out.println("Prefactor = "+Prefactor);
        }
        return result * Prefactor;
    }
    // -------------------------------------------------------------------------

    // canculates overlap integral for radial parts
    //  F1 = SUM[C_i*exp(-zeta_i*(r-R01)^2 )] and F2 = SUM[C_j*exp(-zeta_j*(r-R02)^2 )] with R01 = R02
    // Returns INT( F1*F2, x,y,z=-infinity..+infinity ) = 4*Pi*INT( r^2*F1*F2, r=0..+infinity )
    public double SameCenterRadialOverlapIntegral(double[] coefs1, double[] exponents1, int L1, /* for F1 */
                                 double[] coefs2, double[] exponents2, int L2  /* for F2*/) {
        double result = 0;
        
        for (int c1=0; c1<coefs1.length; c1++)
            for (int c2=0; c2<coefs2.length; c2++)                
                result += coefs1[c1]*coefs2[c2]*primitive_int_1D_Sphr(L1+L2+2, exponents1[c1]+exponents2[c2]);
        return result * 4*Math.PI;
    }

}
//==============================================================================
