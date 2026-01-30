package JGints;

/**
 * A class representing a basis function class suitable for storing
 * basis function data and calculating molecular integrals
 * (currently only overlap ones)
 *
 * Version: 08.Jan.2014
 * pre-Version: 29.Nov.2013
 * Created: 26.Oct.2013
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

import JGints.OverlapIntegrals;
import Polynom3D.*;

//==============================================================================
/**
 * A class representing basis function - a linear superposition
 *  chi(r) = Y_lm(r-R0)*( c1*exp(-z1*rho^2) + ... + cN*exp(-zN*rho^2) ),
 *  where rho = |r-R0|
 * Note: Quick_YLM field should be set before using OverlapWith() method
 * @author timn
 */
public class BasisFunction {
    public BasisFunction _next = null;
    public boolean IsSpherical = true;
    public int L,m; // L = 0...L_MAX, and
        // m = -L...L (if IsSpherical),
        // m = primitive gaussian type (if not IsSpherical): 0(S), 1..3(P), 4..9(D), 10..19(F), 20..34(G) (if !IsSpherical).
    public int additional_r_power = 0;  // for representing radial parts of the form r^additional_r_power·exp(-zeta·r^2)
    public int Center_ID = -1; // 1-based center ID 
    public double[] R0 = null; // center coordinates
    public double[] exponents;
    public double[] coefs;
    public int RadialPart_ID = -1; // ID of exponents/coefs set in some global database (may be used instead of exponents/coefs arrays)
    public double weight = 0.0; // used in NAO procedure only
    public boolean NRB = false; //  used in NAO procedure only
    public OverlapIntegrals OI;
    //--------------------------------------------------------------------------
    /** A constructor;
     * NumCoefs is a numebr of terms in the radial part
     * 
     */
    public BasisFunction(int L, int m, double[] Rcntr, /*double[] exps, double[] coefficients,*/ int NumCoefs) {
        this.L = L;
        this.m = m;
        if (Rcntr != null) this.R0 = Rcntr.clone();
        exponents = new double[NumCoefs];
        //System.arraycopy(exps, 0, exponents, 0, NumCoefs);     // copy first NumCoefs exponents
        coefs = new double[NumCoefs];
        //System.arraycopy(coefficients, 0, coefs, 0, NumCoefs); // copy first NumCoefs coefficients
        // OI to be used by OverlapWith()
        OI = new OverlapIntegrals();
    }
    //--------------------------------------------------------------------------
    /** Re-creates basis function by copying everything from source
     * May also be used to create 'uninitialized' basis function object
     *
     * Rev.: 15.Jan.2014
     */
    public BasisFunction(BasisFunction source) {
        if (source == null) return;
        //
        OI = source.OI;
        _next = source._next;
        IsSpherical = source.IsSpherical;
        L = source.L;
        m = source.m;
        additional_r_power = source.additional_r_power;
        Center_ID = source.Center_ID;
        if (source.R0 != null) this.R0 = source.R0.clone();
        if (source.exponents != null) exponents = source.exponents.clone();
        if (source.coefs != null)  coefs = source.coefs.clone();
        RadialPart_ID = source.RadialPart_ID;
        weight = source.weight;
        NRB = source.NRB;
        // clone source.Quick_YLM
        if (source.Quick_YLM != null) {
            Quick_YLM = new Polynom3D[source.Quick_YLM.length][];
            for (int i=0; i<Quick_YLM.length; i++) {
                Quick_YLM[i] = new Polynom3D[source.Quick_YLM[i].length];
                for (int j=0; j<Quick_YLM[i].length; j++)
                    Quick_YLM[i][j] = new Polynom3D(source.Quick_YLM[i][j]);
            }
        }
    }
    //--------------------------------------------------------------------------
    public Polynom3D[][] Quick_YLM = null; // must be set before calling OverlapWith
    //--------------------------------------------------------------------------
    /**
     * Modifies P as P:=P·(x^2+y^2+z^2)
     */
    public static void _Multiply_Pol3D_By_r2(Polynom3D P) {
        P._Main_To_Temp( P.coefs.length*3 ); // allocate memory for 'temp' storage
        P.temp_num_terms = 0; // remove all elements copied by _Main_To_Temp()
        // place 'main'*x^2 to 'temp
        for (int i=0; i<P.coefs.length; ++i)
            P._Add_Term_To_Temp(P.powers[i][0]+2, P.powers[i][1],   P.powers[i][2],   P.coefs[i]);
        // add 'main'*y^2 to 'temp'
        for (int i=0; i<P.coefs.length; ++i)
            P._Add_Term_To_Temp(P.powers[i][0],   P.powers[i][1]+2, P.powers[i][2],   P.coefs[i]);
        // add 'main'*z^2 to 'temp'
        for (int i=0; i<P.coefs.length; ++i)
            P._Add_Term_To_Temp(P.powers[i][0],   P.powers[i][1],   P.powers[i][2]+2, P.coefs[i]);
        // move 'temp' to 'main'
        P._Temp_To_Main();
    }
    //--------------------------------------------------------------------------
    /*
    private Polynom3D _Apply_Additional_r_power(Polynom3D src, int additional_r_power) throws Exception{
        Polynom3D result = new Polynom3D(src); // we MUST create a new object in order not to change src!
        if ((additional_r_power %2) != 0)
            throw new Exception("Can not calculate overlap if this.additional_r_power=="+additional_r_power+" is odd");
        // ok, additional_r_power = 2*k
        // use result := result · (r-R0)^2 operation for additional_r_power/2 times
        if (additional_r_power>0)
            for (int i=0; i<(additional_r_power/2); ++i)
                _Multiply_Pol3D_By_r2(result);
        return result;
    };
     */
    //--------------------------------------------------------------------------
    /**
     * @returns spherical function with L and m from @param f, multiplied by  (r-R0)^f.additional_r_power
     * If f.additional_r_power, result is a COPY of Quick_YLM[f.L][f.L+f.m]
     * If @param force_recreate, a copy of Quick_YLM[f.L][f.L+f.m] is created regardless of f.additional_r_power value
     *  Rev.:03.May.2014
     */
    public Polynom3D Get_YLM_with_Additional_r_power(BasisFunction f, boolean force_recreate) throws Exception{
        if (Quick_YLM == null)
            throw new Exception("Quick_YLM has not been set!");
        // simply get pointers to pre-calculated spherical functions
        Polynom3D result = Quick_YLM[f.L][f.L+f.m];
        if (force_recreate || (f.additional_r_power>0))
            result = new Polynom3D(result); // we MUST create a new object in order not to change src!
        // now account for (r-R0)^f.additional_r_power
        if (f.additional_r_power>0) {
            if ((f.additional_r_power %2) != 0)
                throw new Exception("Can not calculate overlap if this.additional_r_power=="+f.additional_r_power+" is odd");
            // ok, additional_r_power = 2*k            
            // use result := result · (r-R0)^2 operation for additional_r_power/2 times
            for (int i=0; i<(f.additional_r_power/2); ++i)
                _Multiply_Pol3D_By_r2(result);
            // ALTERNATIVE way is to use a table of pre-calculated YLM·r^(2k) (faster, but more memory)
        }
        return result;
    };
    //--------------------------------------------------------------------------
    /**
     * @returns INT[ YLM1 · SUM_i( c1[i]*exp(-zeta1[i]*(r-R1)^2) ) * YLM2 · SUM_j( c2[j]*exp(-zeta2[j]*(r-R2)^2) )]
     * This method is used by in OverlapWith() for calculating basis functions overlap
     * Note that in contrast to OverlapWith() NO MODIFICATIONS of YLM1/YLM2 are being made before calculation!
     *  Rev.:03.May.2014
     */
    public double _Term_By_Term_Product(BasisFunction a, Polynom3D YLM1, BasisFunction b, Polynom3D YLM2) {
        double result = 0;
        for (int i=0; i<a.coefs.length; i++)
            for (int j=0; j<b.coefs.length; j++) {
                    result += a.coefs[i] * b.coefs[j] *
                            OI.BS_BS_Overlap(a.exponents[i], YLM1, a.R0, b.exponents[j], YLM2, b.R0);
            }

        return result;
    }
    //--------------------------------------------------------------------------
    /** Does a simple term-by-term product
     *  Make sure that Quick_YLM and OI have been properly initialized!
     *  Rev.:03.May.2014
     */
    public double OverlapWith(BasisFunction f) throws Exception {
        // prepare spherical funtions, with ((r-R0)^additional_r_power terms accounted for!
        Polynom3D YLM1 = Get_YLM_with_Additional_r_power(this, false);
        Polynom3D YLM2 = Get_YLM_with_Additional_r_power(f, false);
        /*
        if (Quick_YLM == null) {
            throw new Exception("Quick_YLM has not been set!");
        } else {
            // simply get pointers to pre-calculated spherical functions
            YLM1 = Quick_YLM[L][L+m];
            YLM2 = Quick_YLM[f.L][f.L+f.m];
        }
        // account for ((r-R0)^additional_r_power
        YLM1 = _Apply_Additional_r_power(YLM1, additional_r_power);
        YLM2 = _Apply_Additional_r_power(YLM2, f.additional_r_power);
         */

        /*
        if ((additional_r_power %2) != 0)
            throw new Exception("Can not calculate overlap if this.additional_r_power=="+additional_r_power+" is odd");
        if ((f.additional_r_power %2) != 0)
            throw new Exception("Can not calculate overlap if f.additional_r_power=="+additional_r_power+" is odd");
        if (additional_r_power>0)
            YLM1 = new Polynom3D(YLM1); // prevent spoilong an element of Quick_YLM !!!
        for (int i=0; i<(additional_r_power/2); ++i)
            _Multiply_Pol3D_By_r2(YLM1);
        if (f.additional_r_power>0)
            YLM2 = new Polynom3D(YLM2); // prevent spoilong an element of Quick_YLM !!!
        for (int i=0; i<(f.additional_r_power/2); ++i)
            _Multiply_Pol3D_By_r2(YLM2);
        // ready!
         *
         */
        /*
        
        // perform integration - this code remains the same regardles of additional_r_power and f.additional_r_power
        double result = 0;
        for (int i=0; i<coefs.length; i++)
            for (int j=0; j<f.coefs.length; j++) {
                    result += coefs[i] * f.coefs[j] *
                            OI.BS_BS_Overlap(exponents[i], YLM1, R0, f.exponents[j], YLM2, f.R0);
            }
        return result;
         *
         */
        return _Term_By_Term_Product(this, YLM1, f, YLM2);
    }
    //--------------------------------------------------------------------------
    /** @returns INT( this · (-1/2·/\) · f , dr), where /\ is Laplacian operator
    public double Kinetic_OverlapWith(BasisFunction f) throws Exception {
        // prepare spherical funtions, with ((r-R0)^additional_r_power terms accounted for!
        Polynom3D YLM1 = _YLM_with_Additional_r_power(this, false);
        Polynom3D YLM2 = _YLM_with_Additional_r_power(f, true);
        // 
        
        // perform integration - this code remains the same regardles of additional_r_power and f.additional_r_power
        double result = 0;
        for (int i=0; i<coefs.length; i++)
            for (int j=0; j<f.coefs.length; j++) {
                    result += coefs[i] * f.coefs[j] *
                            OI.BS_BS_Overlap(exponents[i], YLM1, R0, f.exponents[j], YLM2, f.R0);
            }

        return result;
    }
     */
    //--------------------------------------------------------------------------
    /** Returns the value of this basis function at a point r-R0, or just r if R0 == null
     *  Quick_YLM should be set before using this method
     *  Version: 08.Jan.2014
     */
    public double EvaluateAtPoint(double[] r, double[] R0) {
        // dr is the point at which the function should in fact be evaluated
        double[] dr = r.clone();
        if (R0 != null)
            for (int mu=0; mu<3; mu++) dr[mu] -= R0[mu];

        // (dr)^2
        double r2 = 0;
            for (int mu=0; mu<3; mu++) r2 += dr[mu]*dr[mu];
        
        // Evaluate the radial part
        double RadialPart = 0;
        for (int i=0; i<exponents.length; i++)
            RadialPart += Math.exp(-exponents[i] * r2) * coefs[i];
        if (additional_r_power>0)
            RadialPart *= Math.pow(r2, additional_r_power/2.0); // (r2^(1/2))^additional_r_power

        // Evaluate the angular part
        double AngularPart = 0;
        if (!IsSpherical) {
            System.out.println("BasisFunction.EvaluateAtPoint still can not be used with cartesian functions!");
            return 0;
        }
        if (Quick_YLM == null) {
            System.out.println("BasisFunction.EvaluateAtPoint can not be used without setting Quick_YLM");
            return 0;
        }
        AngularPart = Quick_YLM[L][L+m].EvaluateAtPoint(dr);

        // return result :)
        return RadialPart * AngularPart;
    }
    //--------------------------------------------------------------------------
    /** A simple routine for convenient scaling of basis function contraction coefs.
     *  Version: 08.Jan.2014
     */
    public void ScaleCoefsBy(double factor) {
        for (int cf=0; cf<coefs.length; cf++) coefs[cf] *= factor;
    }
    //--------------------------------------------------------------------------
    /**
     *  A simple method for printing exponents and coefficients
     *  Version: 11.Jan.2014
     */
    public void Print() {
        System.out.printf(" %c %4d 1.00%n", "spdfgh".charAt(L), coefs.length);
        if (additional_r_power > 0)
            System.out.println(" and multiply radial part by  r^"+additional_r_power);
        for(int cf=0; cf<coefs.length; ++cf)
            //System.out.printf(" %17.10E %17.10E%n", exponents[cf], coefs[cf]);
            System.out.printf("%20.10f %20.10f %n", exponents[cf], coefs[cf]);
    }
    //--------------------------------------------------------------------------
    final double CF_SMALL_THRESHOLD = 1e-15;
    //--------------------------------------------------------------------------
    /**
     *  Removes primitive Gaussians with zero contraction coefs.
     *  This method is primarily used by the NwChem2Molden program
     *  Version: 25.Dec.2017
     */
    public void removeSmallContractionCoefs() {        
        int nonZeroCoefs = 0;
        for (int i=0; i<coefs.length; ++i) {
            if (Math.abs(coefs[i]) > CF_SMALL_THRESHOLD)
                nonZeroCoefs++;
        }
        if (nonZeroCoefs == coefs.length)
            return; // nothing to clean-up
        // else: create new arrays:
        double[] new_exponents = new double[nonZeroCoefs];
        double[] new_coefs = new double[nonZeroCoefs];
        int k = 0;
        for (int i=0; i<coefs.length; ++i) {
            if (Math.abs(coefs[i]) > CF_SMALL_THRESHOLD) {
                new_exponents[k] = exponents[i];
                new_coefs[k] = coefs[i];
                k++;
            }
        }
        // overwrite old exponents and coefs arrays with the newer ones
        exponents = new_exponents;
        coefs = new_coefs;
    }
    //--------------------------------------------------------------------------

}
