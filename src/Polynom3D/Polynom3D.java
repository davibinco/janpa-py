package Polynom3D;

/**
 * A file containing Polynom3D class which is used for storing (and some
 * manipulation with) coefficients of polynoms like
 * SUM[ C_i * x^nx[i] * y^ny[i] * z^nz[i]  ]
 * and implements a very important C_nk() method
 *
 * Version: 08.Jan.2014
 * Created: 16.Nov.2013
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

import java.util.*;
//==============================================================================
/**
 *  A class representing P(x,y,z) = SUM[ C_i * x^nx[i] * y^ny[i] * z^nz[i]  ]
 *  (also provides a very important subroutine C_nk)
 */
public class Polynom3D {
    public double[] coefs;
    public int[][] powers;  // powrs[i][0] = nx[i], powrs[i][1] = ny[i], powrs[i][2] = nz[i]
    // The coefs and powers array should be as short as possible...
    public double[] coefs_temp;
    public int[][] powers_temp;
    // ...in contrast, coefs_temp and powers_temp array constitue a working buffer which should be rather
    // long if _Add_Term_To_Temp-like methods are to be used.
    public int temp_num_terms; // Number of valid terms in coefs_temp and powers_temp
    //--------------------------------------------------------------------------
    /** Create with empty values
     *
     */
    public Polynom3D(int NumTerms) {
        coefs = new double[NumTerms];
        powers = new int[NumTerms][3];
    }
    //--------------------------------------------------------------------------
    /** To be used for manual construction
     * 
     */
    public Polynom3D(double[] coefs_, int[][] powers_) {
        coefs = coefs_.clone();
        powers = new int[powers_.length][3];
        for (int i=0; i<powers_.length; i++)
            powers[i] = powers_[i].clone();
    }
    //--------------------------------------------------------------------------
    /** Another constructor which copies coefs and powers from src
     *
     */
    public Polynom3D(Polynom3D source) {
        coefs = source.coefs.clone();
        powers = new int[source.powers.length][3];
        for (int i=0; i<source.powers.length; i++)
            powers[i] = source.powers[i].clone();
    }
    //--------------------------------------------------------------------------
    /** Copies first NumTerms from existing arrays
     *
     */
    public void CreateFromArrays(int NumTerms, double[] coefs_, int[][] powers_) {
        // free memory
        coefs = null;
        powers = null;
        // allocate memory
        coefs = new double[NumTerms];
        powers = new int[NumTerms][3];
        System.arraycopy(coefs_, 0, coefs, 0, NumTerms);
        for (int i=0; i<NumTerms; i++)
            powers[i] = powers_[i].clone();
    }
    //--------------------------------------------------------------------------
    /** Sets length to 1, coef[0] to 0.0, powers[0] to {0,0,0} in MAIN arrays
     */
    public void SetZero() {
        coefs = new double[]{ 0.0 };
        powers = new int[][]{new int[]{0,0,0}};
    }
    //--------------------------------------------------------------------------
    /** Sets length to 1, coef[0] to 0.0, powers[0] to {0,0,0} in TEMP arrays
     */
    public void SetTempZero() {
            temp_num_terms = 1;
            coefs_temp = new double[]{ 0.0 };
            powers_temp = new int[][]{new int[]{0,0,0}};
    }
    //--------------------------------------------------------------------------
    /** Searches for the term with given nx,ny,nz and makes "+=COEF" with its coefficient;
     * If no term found - adds it with a given COEF
     */
    public void _Add_Term_To_Temp(int nx, int ny, int nz, double COEF) {
        boolean found = false;
        int i=0;
        while ((i < temp_num_terms) && (!found)) {
            found = (powers_temp[i][0]==nx) && (powers_temp[i][1]==ny) && (powers_temp[i][2]==nz);
            if (!found) i++;
        }
        if (found)
            coefs_temp[i] += COEF;
        else {
            powers_temp[temp_num_terms][0] = nx;
            powers_temp[temp_num_terms][1] = ny;
            powers_temp[temp_num_terms][2] = nz;
            coefs_temp[temp_num_terms] = COEF;
            temp_num_terms++;
        }
    }
    //--------------------------------------------------------------------------
    /** allocates memory for NumTerms in *_temp arrays
     * Note: we rely on Java's initialization of new arrays with zeros!
     */
    public void _Allocate_Temp(int NumTerms) {
        coefs_temp = new double[NumTerms];
        powers_temp = new int[NumTerms][3];
    }
    //--------------------------------------------------------------------------
    /** Removes all terms with abs(coef.) <= epsilon
     * ('<=' is important in order to make a call with epsilon = 0.0 have sense)
     * 
     */
    public void _Remove_small_terms_from_temp(double epsilon) {
        int i=0;
        while (i < temp_num_terms) {
            if ( Math.abs(coefs_temp[i]) <= epsilon) {
                // remove this term
                temp_num_terms--;
                for (int j=i; j<temp_num_terms; j++) {
                    coefs_temp[j] = coefs_temp[j+1];
                    powers_temp[j] = powers_temp[j+1];
                }
            } else i++;
        }
    }
    // copies temp_num_terms terms from *_temp to main arrays
    public void _Temp_To_Main() {
        CreateFromArrays(temp_num_terms, coefs_temp, powers_temp);
    }
    //--------------------------------------------------------------------------
    /** Copies to *_temp from main arrays
     *  NTemporaryTerms should be not less than powers.length (==coefs.length)
     */
    public void _Main_To_Temp(int NTemporaryTerms) {
        coefs_temp = new double[NTemporaryTerms];
        System.arraycopy(coefs, 0, coefs_temp, 0, coefs.length);
        powers_temp = new int[NTemporaryTerms][3];
        for (int i=0; i<powers.length; i++)
            powers_temp[i] = powers[i].clone();
        temp_num_terms = coefs.length;
    }
    //--------------------------------------------------------------------------
    /**
     * Multiplies this by Polynom3D B using temp as a buffer
     * self := self * B
     * Rev.: 03.May.2014
     */
    public void MultiplyBy(Polynom3D B) {
        temp_num_terms = 0; // global variable: number of acutal terms in the polynom product
        // allocate anough memory
        _Allocate_Temp( coefs.length  * B.coefs.length );
        // do optimized multiplication
        for (int i=0; i<coefs.length; i++) {
            for (int j=0; j<B.coefs.length; j++) {
                _Add_Term_To_Temp( powers[i][0] + B.powers[j][0],
                                   powers[i][1] + B.powers[j][1],
                                   powers[i][2] + B.powers[j][2],
                                   coefs[i] * B.coefs[j]);
            }
        }
        this._Remove_small_terms_from_temp(0.0);
        if (temp_num_terms == 0) SetTempZero(); // avoid zero lemgth of the polynom!
        // save result to Main
        _Temp_To_Main();
    }
    //--------------------------------------------------------------------------
    /**
     * Adds (Polynom3D B)*c to this using temp as a buffer
     * self := self + B*c
     * Rev.: 03.May.2014
     */
    public void AddPoly(Polynom3D B, double c) {
        temp_num_terms = 0; // global variable: number of acutal terms in the polynom product
        // allocate anough memory
        _Allocate_Temp( coefs.length + B.coefs.length );
        _Main_To_Temp(coefs.length + B.coefs.length);
        // do an addition
        for (int j=0; j<B.coefs.length; j++)
            _Add_Term_To_Temp( B.powers[j][0], B.powers[j][1], B.powers[j][2], c * B.coefs[j]);
        
        this._Remove_small_terms_from_temp(0.0);
        if (temp_num_terms == 0) SetTempZero(); // avoid zero lemgth of the polynom!
        // save result to Main
        _Temp_To_Main();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns laplacian of this polynom
     */
    public Polynom3D Laplacian() {
        Polynom3D result = new Polynom3D(this);
        result._Allocate_Temp(3 * this.coefs.length);
        temp_num_terms = 0;
        for (int i=0; i<this.coefs.length; ++i) {
            // (d/dx)^2
            if (this.powers[i][0] >= 2) 
                result._Add_Term_To_Temp(powers[i][0]-2, powers[i][1], powers[i][2], 
                        coefs[i]*powers[i][0]*(powers[i][0]-1));
            // +(d/dy)^2
            if (this.powers[i][1] >= 2) 
                result._Add_Term_To_Temp(powers[i][0], powers[i][1]-2, powers[i][2], 
                        coefs[i]*powers[i][1]*(powers[i][1]-1));
            // +(d/dz)^2
            if (this.powers[i][2] >= 2) 
                result._Add_Term_To_Temp(powers[i][0], powers[i][1], powers[i][2]-2, 
                        coefs[i]*powers[i][2]*(powers[i][2]-1));
        }
        result._Remove_small_terms_from_temp(0.0);
        if (result.temp_num_terms == 0) result.SetTempZero(); // avoid zero lemgth of the polynom!
        // save result to Main
        result._Temp_To_Main();
        
        return result;
    }

    //--------------------------------------------------------------------------
    /** returns n!/(k!*(n-k)!); k <= n
     * 
     */
    public static double C_nk(int n, int k) {
        if ((n-k)>k)
            return C_nk(n, n-k);    // since C_nk(n,n-k) = C_nk(n,k)
        else {
            double result = 1;
            // k<=n, so n!/k! = (k+1)*(k+2)*...*n  - has n-k multiplyers
            // while (n-k)!  = 1*2*...*(n-k) - has n-k multiplyers
            // so result *= (k+i)/i for i from 1 to (n-k)
            for (int i=1; i<=(n-k); i++)
                result *= ((1.0*k)/i + 1);
            return result;
        }
    }
    //--------------------------------------------------------------------------
    public static double[][] Quick_CNK = new double[0][0]; // the values[n][k] of n!/k!/(n-k)!
    //--------------------------------------------------------------------------
    /**
     * Checks whether C(n,k) with n,k up to @param max_n are available, and if not,
     * computes the missing values;
     * @returns a 'poiner' to the static array of C[n][k]
     *
     */
    public static double[][] ensure_Cnk_enough(int max_n) {
        if ((Quick_CNK.length-1) < max_n) {
            double[][] new_CNK = new double[max_n+1][];
            // just copy old values
            for (int i=0; i<Quick_CNK.length; ++i) {
                new_CNK[i] = new double[i+1];
                System.arraycopy(Quick_CNK[i], 0, new_CNK[i], 0, Quick_CNK[i].length);
            }
            // and compute new elements
            for (int i=Quick_CNK.length; i<=max_n; ++i) {
                new_CNK[i] = new double[i+1];
                for (int j=0; j<=i; j++)
                    new_CNK[i][j] = Polynom3D.C_nk(i, j);
            }

            Quick_CNK = new_CNK; // update the array
        }
        return Quick_CNK;
    }
    //--------------------------------------------------------------------------
    public String polynom_term_print_format = " + (%.4f)·x^%d·y^%d·z^%d ";
    /** prints a polynom
     * 
     */
    public void Print() {
        /*
        for (int i=0; i<num_terms; i++)
            System.out.printf(" + (%.3f)·x^%d·y^%d·z^%d ", Y_coefs[i], Y_powers[i][0],Y_powers[i][1],Y_powers[i][2]);
        System.out.println("\n=======================================");
         */
        for (int i=0; i<coefs.length; i++)
            System.out.printf(polynom_term_print_format, coefs[i], powers[i][0],
                    powers[i][1], powers[i][2]);
        System.out.println();
    }
    //--------------------------------------------------------------------------
    public void Print2(int L, int m) {
        String COEFS = new String();
        String PWRS = new String();
        Locale.setDefault(new Locale("en", "US"));
        for (int i=0; i<coefs.length; i++) {
            if (i>0) COEFS += ",";
            COEFS += String.format("%25.20f", coefs[i]);
            if (i>0) PWRS += ",";
            PWRS += String.format("{%2d,%2d,%2d}", powers[i][0], powers[i][1], powers[i][2]);
        }
        System.out.printf("   Quick_YLM[%2d][%2d] = new Polynom3D(new double[]{%s}, new int[][]{%s});\n", L,L+m,COEFS,PWRS);
        //Polynom3D A =
    }
    //--------------------------------------------------------------------------
    /**
     * Evaluates a polynom at a given point
     * Version: 08.Jan.2014
     */
    public double EvaluateAtPoint(double[] r) {
        double result = 0;
        int[] max_powers = new int[3]; // masimum powers of x,y and z used in this polynom
        for (int mu=0; mu<3; mu++) max_powers[mu] = 0; // this can be skipped ;)
        for (int i=0; i<powers.length; i++)
            for (int mu=0; mu<3; mu++)
                if (powers[i][mu] > max_powers[mu]) max_powers[mu] = powers[i][mu];

        // calculate all powers of x, y and z we need
        // allocate memory
        double[][] xyz_powers = new double[3][]; // 1-st index: x/y/z; 2-nd index: power of r[mu], xyz_powers[mu][0] = 1 for all mu
        for (int mu=0; mu<3; mu++) {
            xyz_powers[mu] = new double[max_powers[mu] + 1];
            xyz_powers[mu][0] = 1;
        }
        // do calculations
        for (int mu=0; mu<3; mu++) {
            for (int pwr = 1; pwr <= max_powers[mu]; pwr++)
                xyz_powers[mu][pwr] = xyz_powers[mu][pwr-1] * r[mu];
        }

        // calculate a polynom value
        for (int i=0; i<coefs.length; i++)            
            result += coefs[i] * xyz_powers[0][ powers[i][0] ] * xyz_powers[1][ powers[i][1] ] * xyz_powers[2][ powers[i][2] ];

        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * A simple routine for scaling the polynom coefs.
     * (Only 'main' polynom is affected)
     * 
     * Version: 08.Jan.2014
     */
    public void ScaleCoefsBy(double factor) {
        for (int cf=0; cf<coefs.length; cf++) coefs[cf] *= factor;
    }
    //--------------------------------------------------------------------------
    /** Differentiates polynom by r[mu] and returns NEW polynom as a result
     *  Rev.:03.May.2014
     */
    public Polynom3D Differentiate_mu(int mu /*0(x),1(y),2(z)*/) {
        Polynom3D result = new Polynom3D(this);
        for(int i=0; i<result.coefs.length; ++i) {
            if (result.powers[i][mu] == 0)
                result.coefs[i] = 0.0;
            else {
                result.coefs[i] *= result.powers[i][mu];
                --result.powers[i][mu];
            }
        }
        result._Main_To_Temp(result.coefs.length);
        result._Remove_small_terms_from_temp(0.0); // remove terms which are precisely zero
        // it may happen that no terms have been left (e.g., d( 10*y ) / dx)
        // In this case, leave a single term with coef = 0.0
        if (result.temp_num_terms == 0)
            result.SetTempZero();

        result._Temp_To_Main();

        // ready!
        return result;
    }
    //--------------------------------------------------------------------------

}
