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
 * (c) Tymofii Nikolaienko, 2014
 */


/**
 * A class for evaluation of integrals
 *
 *     /(                                                                                   )|
 *    | ((x^2+y^2+(z-a)^2)^n1 * Y[L1][m](x,y,z-a) · (x^2+y^2+(z-b)^2)^n2 * Y[L2][m](x,y,z-b))|                     * df * sin(t)*dt  ===
 *   /  (                                                                                   )| (x=r*sin(t)*cos(f),
 *                                                                                           |  y=r*sin(t)*sin(f),
 *                                                                                           |  z=r*cos(t))
 * === Polynom[l1,n1,l2,n2,m](a, b, r)
 *
 *
 *  Class version: 18.Oct.2014
 *
 * @author timn
 */
public class AngularInts {

    public static double[][] coefs;

    private static double sqrt2 = Math.sqrt(2);
    private static double sqrt3 = Math.sqrt(3);
    private static double sqrt5 = Math.sqrt(5);
    private static double sqrt7 = Math.sqrt(7);
    private static double sqrt10 = Math.sqrt(10);
    private static double sqrtPi = Math.sqrt(Math.PI);
    private static double sqrtPi_3 = Math.pow(Math.PI, 3.0/2);
    //--------------------------------------------------------------------------
    // A 'workhorses' section:
    // fills the coefs[][] array with the coefficients at    __
    //   /                                                   \
    //  | Y[L1][m](x,y,z-a)*Y[L2][m](x,y,z-b)*exp(-r^2)*dr = / C[i][j] * a^i * b^j
    // /                                                     ~~
    //  where i = 0...L1 and j = 0...L2
    //    

     //===========     L1 = 0;  L2 = 0 //---   m = 0
     private static void eval_0_0_p0() {
        coefs[0][0] = sqrtPi_3;
     }
     //===========     L1 = 0;  L2 = 1 //---   m = 0
     private static void eval_0_1_p0() {
        coefs[0][1] = -sqrtPi_3*sqrt3;
     }
     //===========     L1 = 0;  L2 = 2 //---   m = 0
     private static void eval_0_2_p0() {
        coefs[0][2] = sqrtPi_3*sqrt5;
     }
     //===========     L1 = 0;  L2 = 3 //---   m = 0
     private static void eval_0_3_p0() {
        coefs[0][3] = -sqrtPi_3*sqrt7;
     }
     //===========     L1 = 0;  L2 = 4 //---   m = 0
     private static void eval_0_4_p0() {
        coefs[0][4] = 3*sqrtPi_3;
     }
     //===========     L1 = 1;  L2 = 1 //---   m = -1
     private static void eval_1_1_m1() {
        coefs[0][0] = 3.0/2*sqrtPi_3;
     }
     //---   m = 0
     private static void eval_1_1_p0() {
        coefs[1][1] = 3*sqrtPi_3;
        coefs[0][0] = 3.0/2*sqrtPi_3;
     }
     //---   m = 1
     private static void eval_1_1_p1() {
        coefs[0][0] = 3.0/2*sqrtPi_3;
     }
     //===========     L1 = 1;  L2 = 2 //---   m = -1
     private static void eval_1_2_m1() {
        coefs[0][1] = -3.0/2*sqrtPi_3*sqrt5;
     }
     //---   m = 0
     private static void eval_1_2_p0() {
        coefs[1][2] = -sqrtPi_3*sqrt3*sqrt5;
        coefs[0][1] = -sqrtPi_3*sqrt3*sqrt5;
     }
     //---   m = 1
     private static void eval_1_2_p1() {
        coefs[0][1] = -3.0/2*sqrtPi_3*sqrt5;
     }
     //===========     L1 = 1;  L2 = 3 //---   m = -1
     private static void eval_1_3_m1() {
        coefs[0][2] = 3.0/2*sqrtPi_3*sqrt2*sqrt7;
     }
     //---   m = 0
     private static void eval_1_3_p0() {
        coefs[1][3] = sqrtPi_3*sqrt3*sqrt7;
        coefs[0][2] = 3.0/2*sqrtPi_3*sqrt3*sqrt7;
     }
     //---   m = 1
     private static void eval_1_3_p1() {
        coefs[0][2] = 3.0/2*sqrtPi_3*sqrt2*sqrt7;
     }
     //===========     L1 = 1;  L2 = 4 //---   m = -1
     private static void eval_1_4_m1() {
        coefs[0][3] = -3.0/2*sqrtPi_3*sqrt3*sqrt10;
     }
     //---   m = 0
     private static void eval_1_4_p0() {
        coefs[1][4] = -3*sqrtPi_3*sqrt3;
        coefs[0][3] = -6*sqrtPi_3*sqrt3;
     }
     //---   m = 1
     private static void eval_1_4_p1() {
        coefs[0][3] = -3.0/2*sqrtPi_3*sqrt3*sqrt10;
     }
     //===========     L1 = 2;  L2 = 2 //---   m = -2
     private static void eval_2_2_m2() {
        coefs[0][0] = 15.0/4*sqrtPi_3;
     }
     //---   m = -1
     private static void eval_2_2_m1() {
        coefs[1][1] = 15.0/2*sqrtPi_3;
        coefs[0][0] = 15.0/4*sqrtPi_3;
     }
     //---   m = 0
     private static void eval_2_2_p0() {
        coefs[2][2] = 5*sqrtPi_3;
        coefs[1][1] = 10*sqrtPi_3;
        coefs[0][0] = 15.0/4*sqrtPi_3;
     }
     //---   m = 1
     private static void eval_2_2_p1() {
        coefs[1][1] = 15.0/2*sqrtPi_3;
        coefs[0][0] = 15.0/4*sqrtPi_3;
     }
     //---   m = 2
     private static void eval_2_2_p2() {
        coefs[0][0] = 15.0/4*sqrtPi_3;
     }
     //===========     L1 = 2;  L2 = 3 //---   m = -2
     private static void eval_2_3_m2() {
        coefs[0][1] = -15.0/4*sqrtPi_3*sqrt7;
     }
     //---   m = -1
     private static void eval_2_3_m1() {
        coefs[1][2] = -3.0/2*sqrtPi_3*sqrt5*sqrt2*sqrt7;
        coefs[0][1] = -3.0/2*sqrtPi_3*sqrt5*sqrt2*sqrt7;
     }
     //---   m = 0
     private static void eval_2_3_p0() {
        coefs[2][3] = -sqrtPi_3*sqrt5*sqrt7;
        coefs[1][2] = -3*sqrtPi_3*sqrt5*sqrt7;
        coefs[0][1] = -9.0/4*sqrtPi_3*sqrt5*sqrt7;
     }
     //---   m = 1
     private static void eval_2_3_p1() {
        coefs[1][2] = -3.0/2*sqrtPi_3*sqrt5*sqrt2*sqrt7;
        coefs[0][1] = -3.0/2*sqrtPi_3*sqrt5*sqrt2*sqrt7;
     }
     //---   m = 2
     private static void eval_2_3_p2() {
        coefs[0][1] = -15.0/4*sqrtPi_3*sqrt7;
     }
     //===========     L1 = 2;  L2 = 4 //---   m = -2
     private static void eval_2_4_m2() {
        coefs[0][2] = 45.0/4*sqrtPi_3*sqrt3;
     }
     //---   m = -1
     private static void eval_2_4_m1() {
        coefs[1][3] = 15.0/2*sqrtPi_3*sqrt3*sqrt2;
        coefs[0][2] = 45.0/4*sqrtPi_3*sqrt3*sqrt2;
     }
     //---   m = 0
     private static void eval_2_4_p0() {
        coefs[2][4] = 3*sqrtPi_3*sqrt5;
        coefs[1][3] = 12*sqrtPi_3*sqrt5;
        coefs[0][2] = 27.0/2*sqrtPi_3*sqrt5;
     }
     //---   m = 1
     private static void eval_2_4_p1() {
        coefs[1][3] = 15.0/2*sqrtPi_3*sqrt3*sqrt2;
        coefs[0][2] = 45.0/4*sqrtPi_3*sqrt3*sqrt2;
     }
     //---   m = 2
     private static void eval_2_4_p2() {
        coefs[0][2] = 45.0/4*sqrtPi_3*sqrt3;
     }
     //===========     L1 = 3;  L2 = 3 //---   m = -3
     private static void eval_3_3_m3() {
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //---   m = -2
     private static void eval_3_3_m2() {
        coefs[1][1] = 105.0/4*sqrtPi_3;
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //---   m = -1
     private static void eval_3_3_m1() {
        coefs[2][2] = 21*sqrtPi_3;
        coefs[1][1] = 42*sqrtPi_3;
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //---   m = 0
     private static void eval_3_3_p0() {
        coefs[3][3] = 7*sqrtPi_3;
        coefs[2][2] = 63.0/2*sqrtPi_3;
        coefs[1][1] = 189.0/4*sqrtPi_3;
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //---   m = 1
     private static void eval_3_3_p1() {
        coefs[2][2] = 21*sqrtPi_3;
        coefs[1][1] = 42*sqrtPi_3;
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //---   m = 2
     private static void eval_3_3_p2() {
        coefs[1][1] = 105.0/4*sqrtPi_3;
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //---   m = 3
     private static void eval_3_3_p3() {
        coefs[0][0] = 105.0/8*sqrtPi_3;
     }
     //===========     L1 = 3;  L2 = 4 //---   m = -3
     private static void eval_3_4_m3() {
        coefs[0][1] = -315.0/8*sqrtPi_3;
     }
     //---   m = -2
     private static void eval_3_4_m2() {
        coefs[1][2] = -45.0/4*sqrtPi_3*sqrt3*sqrt7;
        coefs[0][1] = -45.0/4*sqrtPi_3*sqrt3*sqrt7;
     }
     //---   m = -1
     private static void eval_3_4_m1() {
        coefs[2][3] = -3*sqrtPi_3*sqrt3*sqrt7*sqrt5;
        coefs[1][2] = -9*sqrtPi_3*sqrt3*sqrt7*sqrt5;
        coefs[0][1] = -45.0/8*sqrtPi_3*sqrt3*sqrt7*sqrt5;
     }
     //---   m = 0
     private static void eval_3_4_p0() {
        coefs[3][4] = -3*sqrtPi_3*sqrt7;
        coefs[2][3] = -18*sqrtPi_3*sqrt7;
        coefs[1][2] = -81.0/2*sqrtPi_3*sqrt7;
        coefs[0][1] = -45.0/2*sqrtPi_3*sqrt7;
     }
     //---   m = 1
     private static void eval_3_4_p1() {
        coefs[2][3] = -3*sqrtPi_3*sqrt3*sqrt7*sqrt5;
        coefs[1][2] = -9*sqrtPi_3*sqrt3*sqrt7*sqrt5;
        coefs[0][1] = -45.0/8*sqrtPi_3*sqrt3*sqrt7*sqrt5;
     }
     //---   m = 2
     private static void eval_3_4_p2() {
        coefs[1][2] = -45.0/4*sqrtPi_3*sqrt3*sqrt7;
        coefs[0][1] = -45.0/4*sqrtPi_3*sqrt3*sqrt7;
     }
     //---   m = 3
     private static void eval_3_4_p3() {
        coefs[0][1] = -315.0/8*sqrtPi_3;
     }
     //===========     L1 = 4;  L2 = 4 //---   m = -4
     private static void eval_4_4_m4() {
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = -3
     private static void eval_4_4_m3() {
        coefs[1][1] = 945.0/8*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = -2
     private static void eval_4_4_m2() {
        coefs[2][2] = 405.0/4*sqrtPi_3;
        coefs[1][1] = 405.0/2*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = -1
     private static void eval_4_4_m1() {
        coefs[3][3] = 45*sqrtPi_3;
        coefs[2][2] = 405.0/2*sqrtPi_3;
        coefs[1][1] = 2025.0/8*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = 0
     private static void eval_4_4_p0() {
        coefs[4][4] = 9*sqrtPi_3;
        coefs[3][3] = 72*sqrtPi_3;
        coefs[2][2] = 243*sqrtPi_3;
        coefs[1][1] = 270*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = 1
     private static void eval_4_4_p1() {
        coefs[3][3] = 45*sqrtPi_3;
        coefs[2][2] = 405.0/2*sqrtPi_3;
        coefs[1][1] = 2025.0/8*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = 2
     private static void eval_4_4_p2() {
        coefs[2][2] = 405.0/4*sqrtPi_3;
        coefs[1][1] = 405.0/2*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = 3
     private static void eval_4_4_p3() {
        coefs[1][1] = 945.0/8*sqrtPi_3;
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }
     //---   m = 4
     private static void eval_4_4_p4() {
        coefs[0][0] = 945.0/16*sqrtPi_3;
     }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /**
     * Redirects to one of the previous functions depending on arguments;
     * For proper operation l1 should be <= than l2
     */
    public static void eval(int l1, int l2, int m) throws Exception{
        // initialize radial[] array with zeros
        //for (int i=0; i<radial.length; ++i) radial[i] = 0;
        for (int i=0; i<coefs.length; ++i)
            for (int j=0; j<coefs[i].length; ++j) coefs[i][j] = 0;

        double val = 0;
         if ((l1 == 0) && (l2 == 0) && (m == 0)) eval_0_0_p0(); else
         if ((l1 == 0) && (l2 == 1) && (m == 0)) eval_0_1_p0(); else
         if ((l1 == 0) && (l2 == 2) && (m == 0)) eval_0_2_p0(); else
         if ((l1 == 0) && (l2 == 3) && (m == 0)) eval_0_3_p0(); else
         if ((l1 == 0) && (l2 == 4) && (m == 0)) eval_0_4_p0(); else
         if ((l1 == 1) && (l2 == 1) && (m == -1)) eval_1_1_m1(); else
         if ((l1 == 1) && (l2 == 1) && (m == 0)) eval_1_1_p0(); else
         if ((l1 == 1) && (l2 == 1) && (m == 1)) eval_1_1_p1(); else
         if ((l1 == 1) && (l2 == 2) && (m == -1)) eval_1_2_m1(); else
         if ((l1 == 1) && (l2 == 2) && (m == 0)) eval_1_2_p0(); else
         if ((l1 == 1) && (l2 == 2) && (m == 1)) eval_1_2_p1(); else
         if ((l1 == 1) && (l2 == 3) && (m == -1)) eval_1_3_m1(); else
         if ((l1 == 1) && (l2 == 3) && (m == 0)) eval_1_3_p0(); else
         if ((l1 == 1) && (l2 == 3) && (m == 1)) eval_1_3_p1(); else
         if ((l1 == 1) && (l2 == 4) && (m == -1)) eval_1_4_m1(); else
         if ((l1 == 1) && (l2 == 4) && (m == 0)) eval_1_4_p0(); else
         if ((l1 == 1) && (l2 == 4) && (m == 1)) eval_1_4_p1(); else
         if ((l1 == 2) && (l2 == 2) && (m == -2)) eval_2_2_m2(); else
         if ((l1 == 2) && (l2 == 2) && (m == -1)) eval_2_2_m1(); else
         if ((l1 == 2) && (l2 == 2) && (m == 0)) eval_2_2_p0(); else
         if ((l1 == 2) && (l2 == 2) && (m == 1)) eval_2_2_p1(); else
         if ((l1 == 2) && (l2 == 2) && (m == 2)) eval_2_2_p2(); else
         if ((l1 == 2) && (l2 == 3) && (m == -2)) eval_2_3_m2(); else
         if ((l1 == 2) && (l2 == 3) && (m == -1)) eval_2_3_m1(); else
         if ((l1 == 2) && (l2 == 3) && (m == 0)) eval_2_3_p0(); else
         if ((l1 == 2) && (l2 == 3) && (m == 1)) eval_2_3_p1(); else
         if ((l1 == 2) && (l2 == 3) && (m == 2)) eval_2_3_p2(); else
         if ((l1 == 2) && (l2 == 4) && (m == -2)) eval_2_4_m2(); else
         if ((l1 == 2) && (l2 == 4) && (m == -1)) eval_2_4_m1(); else
         if ((l1 == 2) && (l2 == 4) && (m == 0)) eval_2_4_p0(); else
         if ((l1 == 2) && (l2 == 4) && (m == 1)) eval_2_4_p1(); else
         if ((l1 == 2) && (l2 == 4) && (m == 2)) eval_2_4_p2(); else
         if ((l1 == 3) && (l2 == 3) && (m == -3)) eval_3_3_m3(); else
         if ((l1 == 3) && (l2 == 3) && (m == -2)) eval_3_3_m2(); else
         if ((l1 == 3) && (l2 == 3) && (m == -1)) eval_3_3_m1(); else
         if ((l1 == 3) && (l2 == 3) && (m == 0)) eval_3_3_p0(); else
         if ((l1 == 3) && (l2 == 3) && (m == 1)) eval_3_3_p1(); else
         if ((l1 == 3) && (l2 == 3) && (m == 2)) eval_3_3_p2(); else
         if ((l1 == 3) && (l2 == 3) && (m == 3)) eval_3_3_p3(); else
         if ((l1 == 3) && (l2 == 4) && (m == -3)) eval_3_4_m3(); else
         if ((l1 == 3) && (l2 == 4) && (m == -2)) eval_3_4_m2(); else
         if ((l1 == 3) && (l2 == 4) && (m == -1)) eval_3_4_m1(); else
         if ((l1 == 3) && (l2 == 4) && (m == 0)) eval_3_4_p0(); else
         if ((l1 == 3) && (l2 == 4) && (m == 1)) eval_3_4_p1(); else
         if ((l1 == 3) && (l2 == 4) && (m == 2)) eval_3_4_p2(); else
         if ((l1 == 3) && (l2 == 4) && (m == 3)) eval_3_4_p3(); else
         if ((l1 == 4) && (l2 == 4) && (m == -4)) eval_4_4_m4(); else
         if ((l1 == 4) && (l2 == 4) && (m == -3)) eval_4_4_m3(); else
         if ((l1 == 4) && (l2 == 4) && (m == -2)) eval_4_4_m2(); else
         if ((l1 == 4) && (l2 == 4) && (m == -1)) eval_4_4_m1(); else
         if ((l1 == 4) && (l2 == 4) && (m == 0)) eval_4_4_p0(); else
         if ((l1 == 4) && (l2 == 4) && (m == 1)) eval_4_4_p1(); else
         if ((l1 == 4) && (l2 == 4) && (m == 2)) eval_4_4_p2(); else
         if ((l1 == 4) && (l2 == 4) && (m == 3)) eval_4_4_p3(); else
         if ((l1 == 4) && (l2 == 4) && (m == 4)) eval_4_4_p4(); else
        throw new Exception(String.format("Unable to integrate YLMs for l1=%d, l2=%d, m=%d ! ", l1,l2,m));
        //return val;
    }
    //--------------------------------------------------------------------------
    // some pre-computed arrays needed by evaluate() method
    public static double[][] zetaEffCoefs = null; // coefs. at d^j( x^(i/2) )/dx^j = kappa[i][j] * 1/x^(i/2+j)
    public static double[][] Cnk; // C[n][k] = n!/k!/(n-k)!
    public static double[][] Ank; // A[n][k] = n!/(n-k)!
    //--------------------------------------------------------------------------
     /**
      * Evaluates the integral
      *                  1               /
      *  I = ------------------------ * |  Y[L1][m](x,y,z-a)*Y[L1][m](x,y,z-a)*exp(-r^2)*dr
      *     (zeta1+zeta2)^(L1+L2+3)/2   /
      *
      * where a = -DZ*zeta2/SQRT(zeta1+zeta2) and b = DZ*zeta1/SQRT(zeta1+zeta2)
      *
      * and its derivatives w.r.t. zeta1 and zeta2 up to (zeta1_derivs, zeta2_derivs)-th order
      *
      * @returns an array[i][j] = diff(I, zeta1$i, zeta2$j)
      *
      * Relies on a 'workhorse' this.eval() which has some limitations on maximum L1 and L2
      *
      * NOTE: proper number of the elements in zetaEffCoefs[][], Cnk[][] and Ank[][] arrays
      * SHOULD be available BEFORE calling this method
      *
      */
    public static double[][] evaluate(int l1, int l2, int m, int zeta1_derivs, int zeta2_derivs,
                                      double zeta1, double zeta2, double DZ) throws Exception{

        boolean transposed = false;
        // interchange {a,b} and {l1,l2} if l2<l1
        if (l2<l1) {
            // swap L
            int tmp = l1;
            l1 = l2;
            l2 = tmp;
            // swap zeta1, zeta2
            double TMP = zeta1;
            zeta1 = zeta2;
            zeta2 = TMP;
            // swap number of derivs
            tmp = zeta1_derivs;
            zeta1_derivs = zeta2_derivs;
            zeta2_derivs = tmp;
            // swap z coordinates
            DZ = -DZ;
            transposed = true;
        }

        // use a 'workhorse' - get the expansion of the intagral (w.o. (zeta1+zeta2)^(-...) ) over powers of a and b
        coefs = new double[l1+1][l2+1];
        eval(l1, l2, m);

        // Compute some auxiliary quantities:
        
        double zetaEff = zeta1+zeta2;
        double invZetaEff = 1/zetaEff;
        double sqrtInvZetaEffDZ = DZ * Math.sqrt(invZetaEff);

        // compute powers of zeta1 & zeta2
        int max_pwr_zetas = l1+l2;        
        double[] zeta1_pwrs = new double[max_pwr_zetas+1];
        zeta1_pwrs[0] = 1;
        double[] zeta2_pwrs = new double[max_pwr_zetas+1];
        zeta2_pwrs[0] = 1;
        for (int i=1; i<=max_pwr_zetas; ++i) {
            zeta1_pwrs[i] = zeta1_pwrs[i-1] * zeta1;
            zeta2_pwrs[i] = zeta2_pwrs[i-1] * zeta2;
        }

        // compute powers of (1/zetaEff)^k, k = 0...zeta1_derivs + zeta2_derivs 
        double[] zetaEffInvPwrs = new double[ zeta1_derivs + zeta2_derivs + 1 ];
        zetaEffInvPwrs[0] = 1.0;
        for (int i=1; i<zetaEffInvPwrs.length; ++i)
            zetaEffInvPwrs[i] = zetaEffInvPwrs[i-1] * invZetaEff;

        // compute d[k] = (1/zetaEff)^(L1+L2+3)/2 * (DZ)^k * (1/zeta_Eff)^(k/2), k = 0...(L1+L2)
        double[] zetaEffMDInvPwrs05 = new double[ l1 + l2+1 ];
        zetaEffMDInvPwrs05[0] = Math.pow(invZetaEff, (l1+l2+3)/2.0);
        for(int i=1; i<zetaEffMDInvPwrs05.length; ++i)
            zetaEffMDInvPwrs05[i] = zetaEffMDInvPwrs05[i-1] * sqrtInvZetaEffDZ;


        // allocate memory for the result
        double[][] result = new double[zeta1_derivs+1][zeta2_derivs+1];

        // compute each term and its derivatives
        for (int i=0; i<=l1; ++i)
            for (int j=0; j<=l2; ++j)
                if (coefs[i][j] != 0.0) {
                    // if there is a non-zero coefficient at (zeta1^i * zeta2^j) * 1/(zeta_Eff)^(i+j+...) == f * g,
                    // compute this term ...
                    if (i%2 == 0) // account for (-1)^i
                        result[0][0] += coefs[i][j] * zeta2_pwrs[i] * zeta1_pwrs[j] * zetaEffMDInvPwrs05[i+j];
                    else
                        result[0][0] -= coefs[i][j] * zeta2_pwrs[i] * zeta1_pwrs[j] * zetaEffMDInvPwrs05[i+j];
                    // ...and its derivatives
                    for (int n=0; n<=zeta1_derivs; ++n)
                        for (int k=0; k<=zeta2_derivs; ++k)
                            if (k+n>0) { // to prevent overwrining result[0][0]
                                double tmp = 0;
                                for(int mu=0; mu<=n && mu<=j; ++mu)
                                    for (int nu=0; nu<=k && nu<=i; ++nu)
                                        tmp += Cnk[n][mu]*Cnk[k][nu] *
                                                Ank[j][mu]*zeta1_pwrs[j-mu] * // d^mu (zeta1^i) / dzeta1^mu
                                                Ank[i][nu]*zeta2_pwrs[i-nu] * // d^nu (zeta2^j) / dzeta2^nu                                                
                                                zetaEffCoefs[i+j+l1+l2+3][n-mu+k-nu] * zetaEffInvPwrs[n+k-mu-nu]; // d^(n+k-nu-mu) ( zetaEff^(-i-j)/2 ) / d zeta1^(n-mu) / d zeta2^(k-nu) :
                                // a common coefficient for all terms in the sum_{mu,nu}
                                tmp *= zetaEffMDInvPwrs05[i+j] * coefs[i][j];
                                if (i%2 == 0)  // account for (-1)^i
                                    result[n][k] += tmp ;
                                else
                                    result[n][k] -= tmp ;
                            }
                }

        // check, whether we've swapped (L1, zeta1, A) <-> (L2 zeta2, B)  before computation
        if (transposed) {
            // transpose the result
            double[][] result2 = new double[zeta2_derivs+1][zeta1_derivs+1];
            for (int i=0; i<result.length; ++i)
                for (int j=0; j<result[i].length; ++j)
                    result2[j][i] = result[i][j];
            return result2;
        } else
            return result; // if there has been no swapping
    }

    //--------------------------------------------------------------------------
}
