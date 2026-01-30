package JGints;

import java.io.*;
import Polynom3D.*;
import ProgramOptions.ProgTimer;
import java.util.Locale;
import moldenio.*;

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
 * (c) Tymofii Nikolaienko, 2014, 2018
 */


/**
 * Some molecular integrals over Gaussian-type basis functions implemented with JAVA
 * @author (c) TimN, 2014, 2018
 */
public class JGintsCyl {
    public BasisFunction[] basis;                       // basis functions to be integrated
    public RadialPartOfBasisFunction[] RadialParts;     // expansion of radial parts of basis functions indexed by basis[].RadialPartID over primitive gaussians
    public AtomicCenter[] centers;                      // coordinates of atoms referenced by basis[].Center_ID-1
    //public double[][] ylm_norms2;                       // ||norms||2 of spherical harmonics to be taken
    //public Polynom3D[][] ylm;                           // a set of spherical harmonics used by BasicOverlapIntegral() - TODO: eliminate!
                                                        // ylm_norms2 must match these expressions (user's responsibility!)!

    public double[][] OverlapMatrix;                    // overlap integrals matrix basis.length x basis.length
    public double[][][] DipoleMatrix;                   // same for dipole integrals: [bf1][bf2][0...2==x,y,z]
                                                        // Note that DipoleMatrix[bf1][bf2] and DipoleMatrix[bf2][bf1] do reference to the _same_ array

    PrintStream out = System.out;

    //--------------------------------------------------------------------------
    /**
     * Computes spherical coordinates of  RA->RB vector and returns:
     * 
     * result[0] = theta = t, result[1] = phi = i
     */
    private double[] ThetaPhi(double[] RA, double[] RB) {
        double[] dr = new double[3];
        for (int mu=0; mu<3; ++mu)
            dr[mu] = RB[mu] - RA[mu];
        double[] result = new double[2];
        result[0] = Math.atan2(Math.sqrt(dr[1]*dr[1]+dr[0]*dr[0]), dr[2]/*dz*/  );
        result[1] = Math.atan2(dr[1]/*dy*/, dr[0]/*dx*/);
        return result;
    }
    //--------------------------------------------------------------------------
    /* General algorithm:
     *   // 1-center normalization ints
     *   loop over all atoms(A)
     *      loop over all L@A
     *          fill in off-diagonal terms with zeros
     *          integrate all radial parts for this L
     *      end loop
     *   end loop
     *   // 2-center overlap ints
     *   loop over unique pairs of centers(A,B)
     *      compute A->B unit vector and sin,cos values
     *      generate all neccessary Winger matrices for ez'=A->B
     *      rotate (B-A) vector to new system
     *      loop over all radial parts R1@A, R2@B
     *          (LA,LB) = Ls of these radial parts
     *          compute 'diagonal' overlap integrals for m=-L_MIN...+L_MIN (note: +m/-m terms will be equal!)
     *          produce all overlap integrals with these L1,L2 in LAB system
     *      end loop
     *   end loop
     *   
     *
     *
     *
     */
    final int[][] P2moldenP = new int[][]{
        new int[]{0}, // S
        new int[]{1,-1,0}, // index = normal value of m(P), value = MOLDEN's value of m(P)
        new int[]{-2,-1,0,1,2},
        new int[]{-3,-2,-1,0,1,2,3},
        new int[]{-4,-3,-2,-1,0,1,2,3,4},
    };
    final int[][] moldenP2P = new int[][]{
        new int[]{0}, // S
        new int[]{0,1,-1}, // index = MOLDEN's value of m(P), value = normal value of m(P)
        new int[]{-2,-1,0,1,2},
        new int[]{-3,-2,-1,0,1,2,3},
        new int[]{-4,-3,-2,-1,0,1,2,3,4},
    };

    //--------------------------------------------------------------------------
    // returns Wigner matrix valid for indexing its elements with MOLDEN-style m's
    double[][] Wigner2MOLDEN(double[][] W, int L) {
        double[][] result = new double[2*L+1][2*L+1];
        for (int m=-L; m<=L; ++m) // MOLDEN's indexes
            for (int M=-L; M<=L; ++M)
                result[L+m][L+M] = W[L + moldenP2P[L][L+m]][L + moldenP2P[L][L+M]];
        return result;
    }
    //--------------------------------------------------------------------------
    double[][] WignerTranspose(double[][] W, int L) {
        double[][] result = new double[2*L+1][2*L+1];
        for (int m=-L; m<=L; ++m) 
            for (int M=-L; M<=L; ++M)
                result[L+m][L+M] = W[L+M][L+m];
        return result;
    }
    //--------------------------------------------------------------------------
    
    private double[] ra_new = new double[]{0,0,0};
    private double[] rb_new = new double[]{0,0,0};
    private Polynom3D YLM1, YLM2;
    private OverlapIntegrals OI = new OverlapIntegrals(); // is not used in the quickest version!
    private NestedFunctionsDerivs derivComputer = new NestedFunctionsDerivs();

    //--------------------------------------------------------------------------
    final private int PrimitiveGOverlapInt_resultsLEN = 4; // length of @param results array    
    /**
     * A central method which evaluates the integrals between two primitive Gaussian
     * functions, each multiplied by a spherical harmonics, centered at points
     * ZA and ZB of the Z axis.
     * (Uses cylindrical coordinate system!)
     * 
     * In particular, for the overlap integral it evaluates
     * 
     *            /
     * INT[m]  = | dxdydz · Y[L1][m](x,y,z-ZA)·(x^2+y^2+(z-ZA)^2)^(addit_r_pwr1/2)·
     *          /           Y[L2][m](x,y,z-ZB)·(x^2+y^2+(z-ZB)^2)^(addit_r_pwr2/2)·
     *                      exp(-zeta1·(x^2+y^2+(z-ZA)^2))·exp(-zeta2·(x^2+y^2+(z-ZB)^2))
     * 
     * for each m from -Lmin to +Lmin and places the results
     * into @param double[] "overlap_diag" (length = 2*Lmin+1)
     * 
     * It also computes non-zero elements of dipole integral matrix [m1][m2],
     * i.e.: 
     *  * dipoleZ_diag[m] == Z_{m,m}, where Z_{m1,m2}==def==< F1(L1,m1)|z|F2(L2,m2) > is
     *                                non-zero only for m1==m2
     *  * dipoleX_up[m] == X_{m,m+1} =/= X_{m+1,m}, 
     * where X_{m1,m2}==def==< F1(L1,m1)|x|F2(L2,m2) > is non-zero only for |m1-m2|==1
     * 
     * 
     * As the result,
     * @param results[0] is filled with (newly allocated) double[]   "overlap_diag"
     * @param results[1] is filled with (newly allocated) double[]   "dipoleZ_diag"
     * @param results[2] is filled with (newly allocated) double[L1+L2]   containing 0.5*<L1,m|rho|L2,m+1>, for m>=0 (!), m=-L1...L2-1
     * @param results[3] is filled with (newly allocated) double[L1+L2+2] containing 0.5*<L1,m|rho|L2,m-1>, for m>=0 (!), m=-L1...L2+1
     *                                                                            /!\~~~/!\
     *                                                                            ~~Attn!~~
     
     *
     *
    *  By using trig. formulas for cos(p +/- q) and sin(p +/- q) it can be shown that for __ m>=0 __:       valid for
        <L1,m | x | L2,m-1> = - <L1,m | y | L2,-m+1>  = 0.5 * Iminus[m] * Pi                                , m>=2 & m<=L1 & m<=L2+1
        <L1,1 | x | L2, 0 >                           = 1.0(!)*Iminus[1] * Pi = 0.2*2*Pi*Iminus[1]          , m==1 
      #                         <L1,1 | y | L2, 0  >  = 0                                                   , m==1 
    
        <L1,m | x | L2,m+1> =   <L1,m | y | L2,-m-1>  = 0.5 * Iplus[m]  * int( (Cos(m*phi))^2, phi=0..2*Pi ), m>=0 & m<=L1 & m<=L2-1
        <L1,-m| x |L2,-m-1> = - <L1,-m| y | L2, m+1>  = 0.5 * Iplus[m]  * int( (Sin(m*phi))^2, phi=0..2*Pi ), m>=0 & m<=L1 & m<=L2-1
    
        <L1,-m| x |L2,-m+1> =   <L1,-m| y | L2, m-1>  = 0.5 * Iminus[m] * Pi                                , m>=2 & m<=L1 & m<=L2+1
                                <L1,-1| y | L2, 0  >  = 1.0(!)*Iminus[1] * Pi = 0.5*2*Pi*Iminus[1]          , m==1
      # <L1,-1| x |L2, 0  >                           = 0                                                   , m==1
    Note that these '0.5*' are _already included_ into the returned arrays of Iplus, Iminus !
    
      where for the given L1, L2, m>=0:
        Iplus[m]  = int(rho*drho*dz * C_{L1,m}*P_{L1,m} * rho * C_{L2,m+1}*P_{L2,m+1}, rho=0..+infinity, z=-infinity..+infinity),
        Iminus[m] = int(rho*drho*dz * C_{L1,m}*P_{L1,m} * rho * C_{L2,m-1}*P_{L2,m-1}, rho=0..+infinity, z=-infinity..+infinity).
      Importantly, the derived equations employ the fact that C_{L,-m} = C_{L,m} and also P_{L,-m}(rho,z) = P_{L,m}(rho,z).
    */
    
    public static void PrimitiveGOverlapInt(double zeta1, int L1, double ZA, int addit_r_pwr1,
                           double zeta2, int L2, double ZB, int addit_r_pwr2,
                           int Lmin, // just min(L1,L2)
                           double[][] results,
                           boolean doDipoles) 
    {
        
        // use the gaussian product theorem
        double zeta_eff = zeta1 + zeta2;
        double z0 = (zeta1*ZA + zeta2*ZB) / zeta_eff;
        double Prefactor = Math.exp(-zeta1*zeta2*(ZB-ZA)*(ZB-ZA) / zeta_eff);

        int sz = L1 + L2 + addit_r_pwr1 + addit_r_pwr2 ; // max. power of variables (z or rho) which can be present in Y_{L1,m}*Y_{L2,m}
        // TODO: is '*2' needed here?
        int m_max = Lmin;

        // compute a series of int( exp(-zeta_eff*z^2) * z^n, z=-inf..+inf)
        double[] z_ints = new double[sz + 2]; // +1 is needed since sz is the max.power ([sz] must be valid), and another +1 is needed for z-dipole integral
        z_ints[0] = Math.sqrt( Math.PI / zeta_eff);
        for(int i=2; i<z_ints.length; i += 2) { // note: odd terms equal zero
            z_ints[i] = z_ints[i-2] * (i-1)/2.0 / zeta_eff; // - d int[-2] / d zeta_eff
        }

        // compute a series of int( exp(-zeta_eff*rho^2) * rho^n, rho=0..+inf),
        // partially using the available z^n integrals
        double[] rho_ints = new double[sz + 3]; // +1 is needed since sz is the max.power ([sz] must be valid), +1 for rho^1 in rho*drho, +1 for rho in x,y-dipoles
        rho_ints[1] = 1/2.0 / zeta_eff;
        for(int i=3; i<rho_ints.length; i+=2) {
            rho_ints[i] = rho_ints[i-2] * (i-1)/2.0 / zeta_eff;   // - d int[-2] / d zeta_eff
        }
/*        for(int i=0; i<z_ints.length; i+=2) {
            rho_ints[i] = z_ints[i] / 2;
        }*/
        // hmm, it's easier to re-compute even powers than to copy results from z_ints
        // taking care of different boundaries...
        rho_ints[0] = Math.sqrt( Math.PI / zeta_eff) * 0.5;
        for(int i=2; i<rho_ints.length; i += 2) { 
            rho_ints[i] = rho_ints[i-2] * (i-1)/2.0 / zeta_eff; 
        }
        
        // allocate memory to store results:
        results[0] = new double[2*Lmin + 1]; // 'm-diagonal' overlap <L1,m|L2,m>, -Lmin<=m<=Lmin

        if (doDipoles) {
            // z-dipole:
            results[1] = new double[Lmin + 1]; // z-dipole (m,m), <L1,m|z|L2,m>, 0<=m<=Mmin, note that Y_{L,-m}(rho,z) == Y_{L,m}(rho,z); need +1 to length to hold [Lmin]-th element

            // x,y-dipoles 'building blocks':
            results[2] = new double[Lmin+1];  // [L1+m] = 0.5*Iplus[m] = <L1,m|rho|L2,m+1> (without sin/cos factor), m=0...min(L1,L2-1), where
            // if L2>L1: min(L1,L2-1)=L1=Lmin;
            // if L2=L1: min(L1,L2-1)=L2-1=Lmin-1;
            // if L2<L1: min(L1,L2-1)=L2-1=Lmin-1;
            // => the largest possible last value of m is Lmin => need Lmin+1 length of array
            results[3] = new double[Lmin+2];// [L1+m] = 0.5*Iminus[m]= <L1,m|rho|L2,m-1> (without sin/cos factor), 
            // valid for m=1...min(L1,L2+1) (note that [0]-the element is not used!), where
            // if L2>L1: min(L1,L2+1)=L1=Lmin;
            // if L2=L1: min(L1,L2+1)=L1=Lmin;
            // if L2=L1-1: min(L1,L2+1)=L1=Lmin+1;
            // if L2<L1-1: min(L1,L2+1)=L2+1=Lmin+1;
            // => the largest possible last value of m is Lmin+1 => need Lmin+2 length of array
        }
        
        Polynom_rho_z M = new Polynom_rho_z(sz, sz); // will hold (rho,z)-dependent part of the (r1^2)^(n1)*<L1,m| and |L2,m'>*(r2^2)^n2 product
        
        // the terms of the (r1^2)^(n1)*<L1,m| and |L2,m'>*(r2^2)^n2 product, common to all further computations        
        M.mul_r2(addit_r_pwr1 / 2, ZA - z0); // max.power += addit_r_pwr1
        //M._printout();
        M.mul_r2(addit_r_pwr2 / 2, ZB - z0);
        //M._printout();        
        
        
        for (int m=0; m<=L1 && m-1<=L2; m++) {
            // Note that the associated Legendre polynomials used for MOLDEN spherical
            // functions are in fact identical for +|m| and -|m|  => make
            // computation for each |m| just once!

            M.pushToStack(); // save (r1^2)^n1 * (r2^2)^n2
            
            
            M.mul_AssocLeg(L1, m, ZA - z0); // max.power of z == max.power of rho == L1
            //M._printout();

            //whether to calc. <L1,m|L2,m> and <L1,m|z|L2,m> :
            boolean diag_available = (m<=L1) && (m<=L2);
            //whether to calc. <L1,m|L2,m+1> :
            boolean Iplus_available = (m+1 <= L2) && doDipoles; // other boundary: m=0 is always ok
            //whether to calc. <L1,m|L2,m-1> :
            boolean Iminus_available = (m >= 1) && (m-1<=L2) && doDipoles; //  //other boundary: m<=Lmin => m-1<Lmin<=L2 => m-1 is always below L2             
            // Each of these three blocks below leaves the M object unmodified in case
            // if any subsequent calc. blocks are to be executed

            // apply prefactor and convert associated Legendre polynomials we've used into the unity-normalized ones:
            double otherTerms1 = Prefactor * SphericalHarmonics.ALegNormalizer[L1][L1+m]; // will remain unchanged for x,y-dipole ints
            

            if (diag_available) {
                if (Iplus_available || Iminus_available)
                    M.pushToStack(); // save (r1^2)^n1 * (r2^2)^n2 * __ Y_{L1,m} __; not neccessary if we're not going to calc neither Iplus, nor Iminus
                
                M.mul_AssocLeg(L2, m, ZB - z0);   // not necessary if we're only going to evaluate Iminus
                //M._printout();

                double res_overlap = 0;
                double res_zDipole = 0;
                for (int i=0; i<=M.last_active_z_power; i++) {
                    //if (i+1 >= z_ints.length)  System.out.println();
                    double tmp = 0.0;
                    for(int j=0; j<=M.last_active_rho_power; j++) {
                        //if (j+1 >= rho_ints.length)  System.out.println();

                        tmp += M.cf_matrix[i][j] * rho_ints[j + 2*m + 1]; // +1 for rho in rho*drho
                    }
                    res_overlap += z_ints[i] * tmp;
                    res_zDipole += z_ints[i+1] * tmp; // +1 for additional power of z
                }

                double phi_int = Math.PI; // <(cos(m*j))^2> * 2*Pi = Pi, if m =/= 0
                if (m==0) {
                    phi_int *= 2;         // <(cos(m*j))^2> * 2*Pi = 2*Pi iff m==0
                }

                double otherTerms = otherTerms1 * phi_int * SphericalHarmonics.ALegNormalizer[L2][L2+m];        
                // Note that (-1)^m is never needed since MOLDEN uses such the spherical
                // functions in which Y_{L, -|m|} ~ sin(|m|*phi), Y_{L, +|m|} ~ cos(|m|*phi)

                // /!\: IMPORTANTLY, ALegNormalizer[L][L+|m|] == ALegNormalizer[L][L-|m|] 

                res_overlap *= otherTerms;
                res_zDipole *= otherTerms;

                // save overlap integral:
                results[0][ Lmin - m ] = res_overlap;
                results[0][ Lmin + m ] = res_overlap;

                // save z-dipole integral: we've just evaluated  res_zDipole == <fA_{L1,m}|(z-z0)|fB_{L2,m}>; 
                // it is the caller's responsibility to shift the origin _somewhere_ back
                if (doDipoles) {
                    results[1][ m ] = res_zDipole ;
                }
                
                
                if (Iplus_available || Iminus_available)
                    M.restoreFromStack(true); // discard mul. by |L2,m>, and set M = (r1^2)^n1 * (r2^2)^n2 * <L1,m|
            }
            
            
            // Now, do Iplus/Iminus, to be used for x,y-dipoles    
            
            // Compute <L1,m|rho|L2,m+1>, if possible
            if (Iplus_available) {
                if (Iminus_available)
                    M.pushToStack(); // save (r1^2)^n1 * (r2^2)^n2 * __ Y_{L1,m} __ ; not necessary, if we're not going to evaluate Iminus
                
                M.mul_AssocLeg(L2, m + 1, ZB - z0);
                //M._printout();

                double res_Iplus = 0.0;

                for (int i=0; i<=M.last_active_z_power; i++) {
                    double tmp = 0.0;
                    for(int j=0; j<=M.last_active_rho_power; j++) {
                        tmp += M.cf_matrix[i][j] * rho_ints[j + 2*m + 3]; // rho powers: +m for Y_{L1,m}, +(m+1) for Y_{L2,m+1}, +1 for rho*drho, +1 for rho*sin or rho*cos
                    }
                    res_Iplus += z_ints[i] * tmp;
                }

                // convert associated Legendre polynomials we've used into the unity-normalized ones
                // and apply prefactor
                res_Iplus *= otherTerms1;
                res_Iplus *= SphericalHarmonics.ALegNormalizer[L2][L2 + m+1]; 

                // save Iplus integral:
                results[2][ m ] = res_Iplus * 0.5; // 0.5: see the results[][] definition above

                if (Iminus_available)
                    M.restoreFromStack(true); // discard mul. by |L2,m+1>
            }

            // compute <L1,m|rho|L2,m-1>, if possible
            if (Iminus_available) {
                M.mul_AssocLeg(L2, m - 1, ZB - z0);
                //M._printout();

                double res_Iminus = 0.0;

                for (int i=0; i<=M.last_active_z_power; i++) {
                    double tmp = 0.0;
                    for(int j=0; j<=M.last_active_rho_power; j++) {
                        tmp += M.cf_matrix[i][j] * rho_ints[j + 2*m + 1]; // rho powers: +m for Y_{L1,m}, +(m-1) for Y_{L2,m-1}, +1 for rho*drho, +1 for rho*sin or rho*cos
                    }
                    res_Iminus += z_ints[i] * tmp;
                }

                // convert associated Legendre polynomials we've used into the unity-normalized ones
                // and apply prefactor
                res_Iminus *= otherTerms1;
                res_Iminus *= SphericalHarmonics.ALegNormalizer[L2][L2+m-1]; 

                // save Iminus integral:
                results[3][ m ] = res_Iminus * 0.5 ;  // 0.5: see the results[][] definition above

                /*if  ((nIterations == 2) && (itr == 0)) {
                    M.restoreFromStack(false); // restore the saved (r1^2)^n1 * (r2^2)^n2
                    m++;
                    M.mul_AssocLeg(L1, m, ZA - z0); // max.power of z == max.power of rho == L1

                    otherTerms1 = Prefactor * SphericalHarmonics.ALegNormalizer[L1][L1+m]; // will remain unchanged for x,y-dipole ints
                }*/
            }
            // Dipoles done!
            
            M.restoreFromStack(true); // restore the saved (r1^2)^n1 * (r2^2)^n2
        }
        
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /**
     * Properly multiplies @param W_LA . d[mu]_Z_frame . @param W_LB == d[mu]_Z_frame
     * where d[mu]Z_frame is (-LA...+LA) x (-LB...+LB) matrix taken 
     * from @param dipole_storage, for each of mu=0...2 (i.e., x,y,z)
     * ('_Z_frame' indicates that @param dipole_storage stores dipole 
     * integrals computed in axis, where A and B nuclei are on Z axis) ;
     * @param LA and @param LB are only used for convenience;
     * @returns array of computed d[mu]_originalAxis, each component being a _scalar_ .
     */
    private static double[] transformDipoleIntegrals(double[][] dipole_storage, 
            int LA, int LB, double[] W_LA, double[] W_LB) 
    {
        double[] result = new double[3];
        /*
        the code could look like this
            for (int ma = -LA; ma <= LA; ++ma) {
                for (int mb = -LB; mb <= LB; ++mb) {
                    result[mu] += W_LA[LA+ma] * (....)[mu] * W_LB[LB+mb]
                }
            }
        but due to nearly-diagonal structure of dipole_storage, the computation
        can be done in linear time
        */
        // elements involving m==0 on either side
        result[2] = W_LA[LA] * dipole_storage[LA][LB] * W_LB[LB]; // <L1,0|z|L2,0>
        if (LA>0) {
            result[1] = W_LA[LA-1] * dipole_storage[LA-1][LB] * W_LB[LB]; // <L1,-1|y|L2,0>
            result[0] = W_LA[LA+1] * dipole_storage[LA+1][LB] * W_LB[LB]; // <L1,1|x|L2,0>
        }
        if (LB>0) {
            result[1] += W_LA[LA] * dipole_storage[LA][LB-1] * W_LB[LB-1]; // <L1,0|y|L2,-1>
            result[0] += W_LA[LA] * dipole_storage[LA][LB+1] * W_LB[LB+1]; // <L1,0|x|L2,1>
        }
        int Lmin = LA;
        if (LB < LA) {
            Lmin = LB;
        }
        
        // do the rest of z-dipoles
        for (int m = 1; m <= Lmin; ++m) { /* equivalent to: m<=LA && m<= LB */
            result[2] += W_LA[LA+m] * dipole_storage[LA+m][LB+m] * W_LB[LB+m]; // <L1,m|z|L2,m>
            result[2] += W_LA[LA-m] * dipole_storage[LA-m][LB-m] * W_LB[LB-m]; // <L1,-m|z|L2,-m>
        }
        
        // Do  (-m-1, m) y-dipoles, for m=1.....   (-m-1 >= -LA => m<=LA-1)
        // and ( m+1,-m) y-dipoles, for m=1.....   (-m-1 >= -LA => m<=LA-1)
        // Do  ( m+1, m) x-dipoles for the same m-range
        // and (-m-1,-m) x-dipoles for the same m-range
        for (int m = 1; m <=LA-1 && m <= LB; ++m) {
            result[1] += W_LA[LA-m-1] * dipole_storage[LA-m-1][LB+m] * W_LB[LB+m]; // <L1,-m-1|y|L2, m>
            result[1] += W_LA[LA+m+1] * dipole_storage[LA+m+1][LB-m] * W_LB[LB-m]; // <L1, m+1|y|L2,-m>
            
            result[0] += W_LA[LA+m+1] * dipole_storage[LA+m+1][LB+m] * W_LB[LB+m]; // <L1, m+1|x|L2, m>
            result[0] += W_LA[LA-m-1] * dipole_storage[LA-m-1][LB-m] * W_LB[LB-m]; // <L1,-m-1|x|L2,-m>
        }
        
        // Do  (-m, m+1) y-dipoles, for m=1.....   (m+1 <=LB => m <= LB-1)
        // and ( m,-m-1) y-dipoles, for m=1.....   (-m-1 >= -LB => m<=LB-1)
        // Do  ( m, m+1) x-dipoles, for the same m-range
        // and (-m,-m-1) x-dipoles, for the same m-range
        for (int m = 1; m <=LA && m <= LB-1; ++m) {
            result[1] += W_LA[LA-m] * dipole_storage[LA-m][LB+m+1] * W_LB[LB+m+1]; // <L1,-m|y|L2, m+1>
            result[1] += W_LA[LA+m] * dipole_storage[LA+m][LB-m-1] * W_LB[LB-m-1]; // <L1, m|y|L2,-m-1>
            
            result[0] += W_LA[LA+m] * dipole_storage[LA+m][LB+m+1] * W_LB[LB+m+1]; // <L1, m|x|L2, m+1>
            result[0] += W_LA[LA-m] * dipole_storage[LA-m][LB-m-1] * W_LB[LB-m-1]; // <L1,-m|x|L2,-m-1>
        }

        return result;
    }
    //--------------------------------------------------------------------------    
    /**
     * MAIN METHOD FOR OVERLAP INTEGRAL MATRIX EVALUATION
     * Quick calculation of overlap integral matrix for basis functions with PURE spherical harmonics
     * Assumes that the coefficients in the basis functions correspond to the use of the
     * spherical harmonics normalized by 1.0 over the full-sphere
     * 
     */
    public void Build_Ovarlap_Matrix_PURE(/*double[][] ylm_norms2_over4Pi*/) {

        int nRadialParts = RadialParts.length;


        //----------------------------------------
        // Quick search arrays:
        //   * atomic center -> all of its radial parts
        //   * radial part -> all of basis functions where it is used
        // ***
        // and some statistics
        int max_L = 0, max_r2_pwr = 0;

        // index = radial part ID, values = IDs of basis functions having this radial part
        int[][] BFNsOf_RP = new int[nRadialParts][];
        int[] N_BFNsOf_RP = new int[nRadialParts];
        for (int b=0; b<basis.length; ++b) {
            ++N_BFNsOf_RP[basis[b].RadialPart_ID];
        }
        // RP->BF: allocate mem and zero N_BFNsOf_RP[] - we'll ues them as incices
        for (int r=0; r<nRadialParts; ++r) {
            BFNsOf_RP[r] = new int[N_BFNsOf_RP[r]];
            N_BFNsOf_RP[r] = 0;
        }
        // fill in BFNsOfAtom array
        { int r;
            for (int b=0; b<basis.length; ++b) {
                // write this number of function to its owning radial part
                r = basis[b].RadialPart_ID;
                BFNsOf_RP[r][ N_BFNsOf_RP[r] ] = b;
                ++N_BFNsOf_RP[r];
            }
        }
        // free auxiliary mem
        N_BFNsOf_RP = null;
        
        // build the array allowing for obtaining bf index with specified 'magnetic number' (m) value 
        // for the given radial part (which also determines the L number)
        int[][] RPart_m2bf = new int[nRadialParts][];
        for (int r=0; r<nRadialParts; ++r) {
            int L_rp = RadialParts[r].LUsedWith;
            RPart_m2bf[r] = new int[ 2*L_rp  + 1 ]; // corresponds to -L_rp...+L_rp
            for(int ib=0; ib<BFNsOf_RP[r].length; ib++) {
                int bf = BFNsOf_RP[r][ib];
                int true_m = this.moldenP2P[ L_rp ][ L_rp + basis[ bf ].m ];
                RPart_m2bf[r][ L_rp + true_m/*basis[ bf ].m*/ ] = bf;
            }
        }
        // Now: RPart_m2bf[ radialPartId ][ L + m] == bf having these: radial part & L & m .
        // Note that m is a 'true m', not 'molden m' (which, of course, differ for p functions only)
        
        
        // ***
        // index = 0-based Center ID, values = IDs of radial parts centered at this atom
        int[][] RPtsOfAtom = new int[centers.length][];
        int[] N_RPtsOfAtom = new int[centers.length];
        for (int rp=0; rp<nRadialParts; ++rp)
            ++N_RPtsOfAtom[ RadialParts[rp].CenterID-1 ];
        // allocate mem and zero N_RPtsOfAtom[] array (!)
        for (int c=0; c<centers.length; ++c) {
            RPtsOfAtom[c] = new int[ N_RPtsOfAtom[c] ];
            N_RPtsOfAtom[c] = 0;
        }
        int c;
        // fill RPtsOfAtom[] array in
        for (int rp=0; rp<nRadialParts; ++rp) {
            c = RadialParts[rp].CenterID-1;
            RPtsOfAtom[c][ N_RPtsOfAtom[c] ] = rp;
            ++N_RPtsOfAtom[c];

            // get statistics on max(L) and max(additional_r2_power):
            if (RadialParts[rp].Addit_r_power >  2*max_r2_pwr)
                max_r2_pwr = RadialParts[rp].Addit_r_power/2;
            if (RadialParts[rp].LUsedWith > max_L)
                max_L = RadialParts[rp].LUsedWith;
        }
        // free auxiliary mem
        N_RPtsOfAtom = null;
        //----------------------------------------
        // quick search arrays ready!
        
        //Polynom3D YLM1, YLM2;

        // 'global' variables, created once;
        int rp1, rp2, LA, LB, L_MIN, bf1, bf2;
        double zeta1, zeta2;
        double cf1, cf2;
        double result_overlap;
/*
        // Pre-compute zetaEffCoefs[][], Cnk[][] and Ank[][] arrays for AngularInts.(...):
        int max_deriv = max_r2_pwr + 1; // '+1' is for safe ))
        // prepare kappa[n][k]
        AngularInts.zetaEffCoefs = NestedFunctionsDerivs.xSqrtNInvDerivs(2*max_L+2*max_L+3 , max_deriv + max_deriv);
        // prepare C_nk
        AngularInts.Cnk = Polynom3D.ensure_Cnk_enough(max_deriv);
        // prepare A_nk
        AngularInts.Ank = NestedFunctionsDerivs.ensure_Ank_enough(max_L);
*/
        

        // The place to store result:
        OverlapMatrix = new double[basis.length][basis.length];
        DipoleMatrix = new double[basis.length][basis.length][3];
        for(int i=0; i<DipoleMatrix.length; ++i) { // alloc. the elements of this matrix. For now, make [i][j] and [j][i] point to the _same_ array
            DipoleMatrix[i][i] = new double[3];
            for(int j=i+1; j<DipoleMatrix.length; ++j) {
                DipoleMatrix[i][j] = new double[3];
                DipoleMatrix[j][i] = DipoleMatrix[i][j]; //   /!\
            }
        }

        // alloc the array holder:
        double[][] results_array = new double[PrimitiveGOverlapInt_resultsLEN][];
        boolean produce_dipoles = true; // in theory, setting it to false could speed-up calculations. TODO: this seem not to be completely the case...
        
        // **********************
        // 'Axial symmetry' case:
        // Compute 'off-diagonal' terms (overlap of basis functions at DIFFERENT centers)
        // **********************
        for (int A=0; A<centers.length; ++A)
            for (int B=(A+1); B<centers.length; ++B) {
                double[] RA = centers[A].R0.clone(); // RA is being changed below!!!
                double[] RB = centers[B].R0.clone();
                
                for (int mu=0; mu<3; ++mu) { RB[mu] -= RA[mu]; RA[mu] = 0; }
                
                // spherical coordinates of RA->RB w.r.t. LAB system
                double[] ti = this.ThetaPhi(RA, RB);
                double[][] rot_matrix = RealWignerAdapt.newAxisCoords(ti[0] /*theta*/, ti[1] /*phi*/);

                // rotate vectors
                double[] RA_new = new double[]{0,0,0}; // by definition!
                double d_ab = Math.sqrt(RB[0]*RB[0] + RB[1]*RB[1] + RB[2]*RB[2]);
                double[] RB_new = new double[]{0,0,d_ab}; // by definition RB_new has Z'-component only!               
                /* cartesian components of RA->RB w.r.t. LAB system
                double[] eZ = this.ThetaPhi_2_n(ti);
                this.rotate_vec_from_new(RB, eZ); -- this should give {0,0, |RB-RA|} */

                double[][][] sincos = RealWignerAdapt.create_sincos(ti[0], ti[1]);

                // get max L for this pair
                int L_MAX_A = 0;
                for (int i=0; i<RPtsOfAtom[A].length; ++i)
                    if (RadialParts[RPtsOfAtom[A][i]].LUsedWith > L_MAX_A)
                        L_MAX_A = RadialParts[RPtsOfAtom[A][i]].LUsedWith;
                int L_MAX_B = 0;
                for (int i=0; i<RPtsOfAtom[B].length; ++i)
                    if (RadialParts[RPtsOfAtom[B][i]].LUsedWith > L_MAX_B)
                        L_MAX_B = RadialParts[RPtsOfAtom[B][i]].LUsedWith;
                int L_MAX = L_MAX_A;
                if (L_MAX_B>L_MAX) L_MAX = L_MAX_B;

                // calculate Wigner matrices for this particular A->B direction
                double[][][] W = new double[L_MAX+1][][]; // 1-st index = L
                for (int L=0; L<=L_MAX; ++L) {
                    W[L] = RealWignerAdapt.Wigner(sincos, L);
                    W[L] = WignerTranspose(W[L], L); // its sad, so sa-a-d,...
                                                     // to transpose that pretty matrix...
                }


                // loop over radial parts centered at A
                for (int Irp1=0; Irp1<RPtsOfAtom[A].length; ++Irp1) {
                    rp1 = RPtsOfAtom[A][Irp1];
                    LA = RadialParts[rp1].LUsedWith;                    
                    for (int Irp2=0; Irp2<RPtsOfAtom[B].length; ++Irp2) {
                        // we loop over UNIQUE pairs of atoms => no double couting of radial parts
                        rp2 = RPtsOfAtom[B][Irp2];
                        LB = RadialParts[rp2].LUsedWith;
                        L_MIN = LA;
                        if (LB<L_MIN) L_MIN = LB;
                        
                        // alloc mem for 'diagonal' overlap integrals between the basis functions (not just primitive gaussians!)
                        double[] ovarlap_Diags = new double[2*L_MIN+1];
                        
                        double[][] dipole_storage = new double[2*LA+1][2*LB+1]; // to store dipole ints [LA+mA][LB+mB]
                        /* The use of the storage:  
                        1st   |
                        index:|
                         mA=  |      [-LB]        [-3]      [-2]      [-1]       [0]      [1]      [2]        [3]          [LB]   <-- 2-nd index == mB value
                        [-LA] [  _(-LA,-LB) ..........................................................................  _(-LA,LB)  ]
                              [     ...            ...         ...      ...       ...      ...      ...     ...             ...    ]
                         [-3] [  _(-3,-LB)  ... Z(-3,-3)  x(-3,-2)                               y(-3,2)            ...  _(-3,LB)  ]
                         [-2] [   (-2,-LB)  ... x(-2,-3)  Z(-2,-2)  x(-2,-1)            y(-2,1)            y(-2,3)  ...  _(-2,LB)  ]
                         [-1] [   (-1,-LB)  ...           x(-1,-2)  Z(-1,-1)   y(-1,0)           y(-1,2)            ...   (-1,LB)  ]
                         [ 0] [   _(0,-LB)  ...                      y(0,-1)   Z_<0,0>   x(0,1)                     ...   _(0,LB)  ]
                         [ 1] [   _(1,-LB)  ...            y(1,-2)              x(1,0)   Z(1,1)   x(1,2)            ...   _(1,LB)  ]
                         [ 2] [   _(2,-LB)  ...  y(2,-3)             y(2,-1)             x(2,1)   Z(2,2)    x(2,3)  ...   _(2,LB)  ]
                         [ 3] [   _(3,-LB)  ...            y(3,-2)                                x(3,2)    Z(3,3)  ...   _(3,LB)  ]
                              [     ...            ...         ...      ...       ...      ...      ...     ...             ...    ]
                         [LA] [   _(LA,-LB) ..........................................................................   _(LA,LB)  ]

                            i.e., Z component at main diagonal, X component at +/-1 from main diagonal, Y at +/-1 from the main 'anti-diagonal'
                        */
                        

                        // now loop over ALL pairs of exponents(i.e., primitive gaussians) in order to get 
                        // all integrals w.r.t. the 'zeta-middle' points
                        for (int expon1=0; expon1<RadialParts[rp1].Exponents.length; ++expon1) {
                            zeta1 = RadialParts[rp1].Exponents[expon1];
                            cf1 = RadialParts[rp1].Coefs[expon1];
                            //diags[expon1] = new double[molden.RadialParts[rp2].Exponents.length][];
                            for (int expon2=0; expon2<RadialParts[rp2].Exponents.length; ++expon2) {
                                zeta2 = RadialParts[rp2].Exponents[expon2];
                                
                                // evaluate (m1,m2)-matrix of the overlap and dipole integrals                                
                                PrimitiveGOverlapInt(zeta1, LA, RA_new[2], RadialParts[rp1].Addit_r_power,
                                                     zeta2, LB, RB_new[2], RadialParts[rp2].Addit_r_power,
                                                     L_MIN,
                                                     results_array, produce_dipoles);
                                // use results_array[0] to generate <sum(gsn)*Y_{L1,m} | sum(gsn)*Y_{L2,m)> overlap
                                // integrals between the whole basis functions ( gsn = sum{ cf_i*exp(-zeta_i*...) } ):

                                // combine overlap integrals with radial parts 'contraction coefs.' and store result                                
                                // 1) overlap 
                                double cf_prod12 = cf1 * RadialParts[rp2].Coefs[expon2];
                                for (int m=-L_MIN; m<=L_MIN; ++m) {
                                    ovarlap_Diags[L_MIN+m] += results_array[0][L_MIN+m] * cf_prod12;
                                }
                                
                                if (produce_dipoles) {
                                    // 2) z-dipole terms, TODO: check "Pi" part
                                    for (int m=0; m<=L_MIN; ++m) {
                                        // shift from zeta-middle point to RA={0,0,ZA}, and apply cf1*cf2:
                                        // results_array[1][m] = <L1,m|z-z0|L2,m> = <L1,m|z-zA|L2,m> + (zA-z0) * <L1,m|L2,m>,
                                        // so that dz w.r.t. RA is: <L1,m|z-zA|L2,m> = <L1,m|z-z0|L2,m> - (zA-z0) * <L1,m|L2,m>, where 
                                        // z0 = (zeta1*zA + zeta2*zB)/(zeta1+zeta2), zA-z0 = (zA-zB)*zeta2/(zeta1+zeta2), zA=0, zB=d_ab

                                        results_array[1][m] += d_ab * zeta2/(zeta1+zeta2) * results_array[0][L_MIN+m];

                                        dipole_storage[LA+m][LB+m] += results_array[1][m] * cf_prod12; // z-dipole
                                        dipole_storage[LA-m][LB-m] = dipole_storage[LA+m][LB+m]; // TODO: optimize -- do it later!
                                    }
                                    // Note that no shift is necessary for x,y-dipoles at this stage!
                                    // Note that at this stage we only shift dipoles to RA point (which will remain invariant after axis rotation)

                                    double[] Iplus  = results_array[2];
                                    double[] Iminus = results_array[3];

                                    // 2) x,y-dipole
                                    // L_MIN is never less than zero; 
                                    // considering m==0 separately avoids double-assignment for dipole_storage[LA+/-0][...]
                                    // and set non-zero matrix elements of x and y involving Y{L,0} (at either side)
                                    // <L1,1 | x | L2, 0 > = <L1,-1| y | L2, 0  > = 0.5 * Iminus[1] * 2?*Pi, m==1 
                                    // <L1,0 | x | L2, +1> = <L1,0 | y | L2,  -1>  = 0.5 * Iplus[0]  * 2*Pi, m == 0
                                    // Note that 0.5 is already included into Iplus/Iminus

                                    if (LB>0) {
                                        // y:
                                        dipole_storage[LA][LB-1] += cf_prod12 * Iplus[0] * 2*Math.PI; // <L1,0|y|L2,-1>; m==0 => 2*Pi
                                        // <L1,0|y|L2,+1> == 0
                                        // x:
                                        dipole_storage[LA][LB+1] = dipole_storage[LA][LB-1]; // <L1,0|x|L2,+1> == <L1,0|y|L2,-1>
                                        // <L1,0|x|L2,-1> == 0
                                    }
                                    if (LA>0) {
                                        // y:
                                        dipole_storage[LA-1][LB] += cf_prod12 * Iminus[1] *2* Math.PI; //<L1,-1|y|L2,0>; m==1 => 1*Pi
                                        // <L1,+1|y|L2,0>==0
                                        // x:
                                        dipole_storage[LA+1][LB] = dipole_storage[LA-1][LB]; // <L1,+1|x|L2,0> == <L1,-1|y|L2,0>
                                        // <L1,-1|x|L2,0> == 0
                                    }
                                    /* remaining non-zero elements are:
                                        <L1,m | x | L2,m-1> = - <L1,m | y | L2,-m+1>  = 0.5 * Iminus[m] * Pi  , m>=2 & m<=L1 & m<=L2+1
                                        <L1,-m| x |L2,-m+1> =   <L1,-m| y | L2, m-1>  = 0.5 * Iminus[m] * Pi  , m>=2 & m<=L1 & m<=L2+1

                                        <L1,m | x | L2,m+1> =   <L1,m | y | L2,-m-1>  = 0.5 * Iplus[m]  * Pi  , m>=1 & m<=L1 & m<=L2-1
                                        <L1,-m| x |L2,-m-1> = - <L1,-m| y | L2, m+1>  = 0.5 * Iplus[m]  * Pi  , m>=1 & m<=L1 & m<=L2-1
                                    */
                                    // the ones involving Iminus
                                    for (int m=2; (m <= LA/*L_MIN*/) && (m-1 <= LB); ++m) {
                                        double Iminus_term = Iminus[m] * cf_prod12 * Math.PI;
                                        dipole_storage[LA+m][LB+m-1] += Iminus_term;
                                        dipole_storage[LA+m][LB-m+1] -= Iminus_term;

                                        dipole_storage[LA-m][LB-m+1] += Iminus_term;
                                        dipole_storage[LA-m][LB+m-1] += Iminus_term;
                                    }
                                    // the ones involving Iplus
                                    for (int m=1; (m <= LA/*L_MIN*/) && (m+1 <= LB); ++m) {
                                        double Iplus_term = Iplus[m] * cf_prod12 * Math.PI;
                                        dipole_storage[LA+m][LB+m+1] += Iplus_term;
                                        dipole_storage[LA+m][LB-m-1] += Iplus_term;

                                        dipole_storage[LA-m][LB-m-1] += Iplus_term;
                                        dipole_storage[LA-m][LB+m+1] -= Iplus_term;
                                    }
                                } // if (produce_dipoles)
                                
                            } // for expon2
                        }
                        
                        // Ok, now ovarlap_Diags[] and dipole_storage[][] contains
                        // all non-zero (mA,mB)-conponents of overlap/dipole inegrals
                        // involving the whole basis functions (i.e., each already summed 
                        // over its constituting primitive gaussians)
                        
                        // Use these stored integrals to produce ALL integrals involving these radial parts
                        // by rotating spherical functions (and dipole vector) back to initial axis orientation

                        // loop over all basis functions which use these radial parts
                        for (int b1I=0; b1I<BFNsOf_RP[rp1].length; ++b1I) {
                            bf1 = BFNsOf_RP[rp1][b1I];
                            // номера функций, которую мы раскладываем, в НОРМАЛЬНОЙ m-нумерации
                            int fixed_mA = LA + moldenP2P[LA][LA+basis[bf1].m];
                            for (int b2I=0; b2I<BFNsOf_RP[rp2].length; ++b2I) {
                                bf2 = BFNsOf_RP[rp2][b2I];
                                // номера функций, которую мы раскладываем, в НОРМАЛЬНОЙ m-нумерации
                                int fixed_mB = LB + moldenP2P[LB][LB+basis[bf2].m];

                                result_overlap = 0;
                                for (int m=-L_MIN; m<=L_MIN; ++m) {
                                    result_overlap += ovarlap_Diags[L_MIN+m] * W[LA][fixed_mA][LA+m] * W[LB][fixed_mB][LB+m];                                    
                                }

                                OverlapMatrix[bf1][bf2] = result_overlap;
                                OverlapMatrix[bf2][bf1] = result_overlap;
                                
                                if (produce_dipoles) {
                                    // Do the same rotations of __spherical functions__ for the dipole moment                                
                                    double[] dipoles = transformDipoleIntegrals(dipole_storage, LA, LB, W[LA][fixed_mA], W[LB][fixed_mB]);

                                    // Transform the computed __vector__ back to original axis orientation

                                    // where to put the result:
                                    double[] dipoles_origAx = DipoleMatrix[bf1][bf2]; // points to the same (!already allocated!, as well) element as DipoleMatrix[bf2][bf1]


                                    for (int mu=0; mu<3; mu++) {
                                        dipoles_origAx[mu] = 0.0;
                                        for (int nu=0; nu<3; nu++) {
                                            dipoles_origAx[mu] += rot_matrix[mu][nu] * dipoles[nu];
                                        }
                                        // we've just computed <bf1| r_vec - RA |bf2> = <bf1| r_vec | bf2>  - RA <bf1 |bf2>
                                        // => need to shift this back to original origin by adding RA <bf1 |bf2> = RA * overlap[bf1][bf2]
                                        dipoles_origAx[mu] += centers[A].R0[mu] * result_overlap; // note than RA has been modified above! => need use original value!
                                    }

                                    // debug:
                                    // System.out.printf(" %d (L=%d) %d (L=%d) : [ %.7f  %.7f  %.7f] %n", bf1, LA, bf2, LB, dipoles_origAx[0], dipoles_origAx[1], dipoles_origAx[2]);
                                } // if produce_dipoles
                            }
                        } // for b1I
                        //ovarlap_Diags = null;
                        
                    } // for RadialPartIndex2
                } // for RadialPartIndex1

            } // for atom B
        // for atom A


        // **********************
        // 'Point symmetry' case
        // **********************
        // Calculate 'diagonal' terms of overlap matrix
        int addit_r_pwr1, addit_r_pwr2;
        for (int A=0; A<centers.length; ++A) {
            // loop over all radial parts used for this center
            for (int Irp1=0; Irp1<RPtsOfAtom[A].length; ++Irp1) {
                rp1 = RPtsOfAtom[A][Irp1];
                LA = RadialParts[rp1].LUsedWith;
                // Note that we CAN assume that the origin coincides with atom A
                double RP_product = 0;
                addit_r_pwr1 = RadialParts[rp1].Addit_r_power;

                // loop over all the other radial parts
                for (int Irp2=Irp1; Irp2<RPtsOfAtom[A].length; ++Irp2) {
                    rp2 = RPtsOfAtom[A][Irp2];
                    LB = RadialParts[rp2].LUsedWith;
                    boolean do_calc_overlap = (LB == LA);
                    boolean do_calc_dipole = (Math.abs(LB - LA) == 1) && produce_dipoles;
                    
                    /*if (RadialParts[rp2].LUsedWith != LA)
                        continue;*/                    
                    if ((!do_calc_overlap) && (!do_calc_dipole))
                        continue; //  // go to the next iteration, if: harmonics with different L will be orthogonal, AND dipole moment matrix element is zero

                    addit_r_pwr2 = RadialParts[rp2].Addit_r_power;

                    RP_product = 0;
                    double RP_product_LALB = 0; // used for dipoles
                                
                    for (int c1=0; c1<RadialParts[rp1].Coefs.length; ++c1) {
                        zeta1 = RadialParts[rp1].Exponents[c1];
                        cf1 = RadialParts[rp1].Coefs[c1];
                        for (int c2=0; c2<RadialParts[rp2].Coefs.length; ++c2) {
                            zeta2 = RadialParts[rp2].Exponents[c2];
                            cf2 = RadialParts[rp2].Coefs[c2];
                            // TODO: optimize by joining two integral routines into one, passing do_calc_dipole as a parameter
                            if (do_calc_overlap) {
                                RP_product += OverlapIntegrals.primitive_int_1D_Sphr(
                                        2/*·r^2*/+2*LA/*(r^l)^2*/+addit_r_pwr1+addit_r_pwr2/*(r^k1)*(r^k2)*/,
                                        zeta1 + zeta2) * cf1*cf2;
                            }
                            if (do_calc_dipole) {
                                RP_product_LALB += OverlapIntegrals.primitive_int_1D_Sphr(
                                        2/*·r^2*/+LA+LB+1/*(r^L1)*(r^L2)*r*/+addit_r_pwr1+addit_r_pwr2/*(r^k1)*(r^k2)*/,
                                        zeta1 + zeta2) * cf1*cf2;

                                /*PrimitiveGOverlapInt(zeta1, LA, 0.0, 0,
                                        zeta2, LB, 0.0, 0,
                                        (LA<LB)?LA:LB, r);
                                System.out.printf("%.7f %n", r[1][0] * cf1*cf2);
                                System.out.printf("%.7f %n", r[2][0] * cf1*cf2);
                                System.out.printf("%.7f %n", r[3][0] * cf1*cf2);*/                                
                            }
                        }
                    }
                    // use these results to produce same-m overlaps for this radial part of this center

                    // loop over ALL basis functions which use these radial parts
                    for (int bI1=0; bI1<BFNsOf_RP[rp1].length; ++bI1) {
                        bf1 = BFNsOf_RP[rp1][bI1];
                        int m1 = basis[bf1].m;
                        for (int bI2=0; bI2<BFNsOf_RP[rp2].length; ++bI2) {
                            // Ls are the same for rp1 and rp2 - we've already checked this! (see do_calc_overlap)
                            // are m's equal?
                            bf2 = BFNsOf_RP[rp2][bI2];
                            int m2 = basis[bf2].m;
                            if ((m2 == m1) && (do_calc_overlap))  {
                                result_overlap = RP_product ;  // note that ||Y_{lm}||^2 == 1  !
                                OverlapMatrix[bf1][bf2] = result_overlap;
                                OverlapMatrix[bf2][bf1] = result_overlap;
                            }
                        }
                    }
                    
                    if (do_calc_dipole) {
                        setAtomicDipolesNonDiag(RPart_m2bf, LA, LB, rp1, rp2, RP_product_LALB);                        
                    }
                    
                    // and all the other terms are zeros!                    
                    // Note, however, that these terms are zeros w.r.t. the nuclei.
                    // But since basis functions at the atom A in general need not be orthogonal,
                    // this can cause some elements of the dipole matrix (the diagonal ones, at least)
                    // be non-zero w.r.t. initial origin.
                    
                } // for 2-nd radial part of atom A                               
            } // for 1-st radial part of atom A
        } // for atom A

        
        // Now, when the full OverlapMatrix is available, we can safely perform 
        // a 'shift' for the 'diagonal' parts of the single-atomic dipole moments:
        if (produce_dipoles) {
            for (int A=0; A<centers.length; ++A) {
                for (int Irp1=0; Irp1<RPtsOfAtom[A].length; ++Irp1) {
                    rp1 = RPtsOfAtom[A][Irp1];
                    for (int Irp2=Irp1; Irp2<RPtsOfAtom[A].length; ++Irp2) {
                        rp2 = RPtsOfAtom[A][Irp2];

                        // We have < Y_lm1(r-R) * rp1(|r-R|) | x - Rx | Y_lm2(r-R) * rp2(|r-R|) > = < Y_lm1 * rp1 | x | Y_lm2 * rp2 > -Rx * < Y_lm1 * rp1 | Y_lm2 * rp2 >,
                        // so that < Y_lm1 * rp1 | x | Y_lm2 * rp2 > = < Y_lm1 * rp1 | x - Rx | Y_lm2 * rp2 >  + Rx * < Y_lm1 * rp1 | Y_lm2 * rp2 >
                        // Note that by now, < Y_lm1 * rp1 | x - Rx | Y_lm2 * rp2 > (possibly, equal to zero, possibly not!) is already saved in DipoleMatrix
                        // => we need only to += additional 'shift term'
                        // Thus, these remaining elements of dipole matrix reduces to adding the (already computed) overlap matrix

                        // In principle, this code works of when placed inside the previous 'global' rp1/rp2 loop,
                        // but it looks safer to have it here, when OverlapMatrix has been obtained

                        double[] R0 = centers[  A ].R0; // same as centers[  RadialParts[rp2].CenterID ].R0 since both rp1 and rp2 belong to the same atom by definition !

                        for (int bI1=0; bI1<BFNsOf_RP[rp1].length; ++bI1) {
                            bf1 = BFNsOf_RP[rp1][bI1];
                            for (int bI2=0; bI2<BFNsOf_RP[rp2].length; ++bI2) {                            
                                bf2 = BFNsOf_RP[rp2][bI2];
                                if (bf2 < bf1)
                                    continue; // we MUST do this to avoid double-counting due to the fact that DipoleMatrix[bf1][bf2] and DipoleMatrix[bf2][bf1] point to the same array 
                                for(int mu=0; mu<3; mu++) {
                                    DipoleMatrix[bf1][bf2] [mu] += R0[mu] * OverlapMatrix[bf1][bf2];
                                    //System.out.printf("%s: setting %d %d to %.5f %n", "xyz".charAt(mu), bf1, bf2, DipoleMatrix[bf1][bf2] [mu] );
                                }                            
                            }
                        }
                    }
                }
            }
        } // if produce_dipoles
        
        // READY!
    }
    //--------------------------------------------------------------------------
    /**
     * Computes matrix elements of the \vec{r} components between the basis
     * functions belonging to the same atom, and writes the computed values
     * into DipoleMatrix[][].
     * @param rp1, @param rp2 -- 'radial parts' numbers (which identify the corresponding sums of primitive gaussians),
     * @param LA, @param LB -- 'angular momentum' corresponding to rp1,rp2,
     * @param RPart_m2bf[][] -- gets [rp][L+m], where rp identifies the 'radial part' and m is the 'true' 
     *            (i.e., y,z,x for P functions) 'magnetic number', and yields the global index of 
     *            the basis function with the given radial part, L and m,
     * @param RP_product_LALB -- pre-computed integral < r^LA * rp1(r)| r | r^LB * rp2(r) >, where rp1(r) and rp2(r)
     *   denote the sums of primitive gaussians with (r^2)^additional_r2_power included.
     * 
     * NOTE: the full overlap matrix for the given radial parts rp1 and rp2 must be available before calling
     *       this method !
     * 
     * Rev. 22-Jul-2018
     */
    private void setAtomicDipolesNonDiag(int[][] RPart_m2bf, int LA, int LB, int rp1, int rp2, double RP_product_LALB) {
        int bf1;
        int bf2;
        
        // set the values of (only!) non-zero elements of x,y-dipoles:
        // the only non-zero elements are those having <LA, m | x or y | LB, m+/-1>
        // => loop only over possible m values which satisfy these ralations

        // z-dipoles:
        
        for (int m=0; m <= LA && m <= LB; ++m) {
            double YLM_factor;
            if (LA > LB) { 
                // <L+1==LA, m | z | L==LB, m>
                YLM_factor = (LB-m+1.0)*(LB+m+1.0)/(2*LB+1.0)/(2*LB+3.0);                                        
            } else { // note that |LA-LB|==1 due to do_calc_dipole==true
                // <L-1==LA, m | z | L==LB, m>
                YLM_factor = (LB-m)*(LB+m)/(2*LB-1.0)/(2*LB+1.0);
            }
            YLM_factor = Math.sqrt(YLM_factor);

            bf1 = RPart_m2bf[rp1][LA + m];
            bf2 = RPart_m2bf[rp2][LB + m];            
            DipoleMatrix[bf1][bf2][2] = YLM_factor * RP_product_LALB; // no need to change DipoleMatrix[bf2][bf1] 

            // debug:
            //System.out.printf("Z: %d (L=%d,m=%d) %d (L=%d,m=%d) : [ %.7f  %.7f  %.7f] %n", bf1, LA,m, bf2, LB,m, 0.0, 0.0, DipoleMatrix[bf1][bf2][2]);
            if (m>0) { // not needed  for m == 0
                bf1 = RPart_m2bf[rp1][LA - m];
                bf2 = RPart_m2bf[rp2][LB - m];
                DipoleMatrix[bf1][bf2][2] = YLM_factor * RP_product_LALB; // no need to change DipoleMatrix[bf2][bf1]             

                // debug:
                //System.out.printf("Z: %d (L=%d,m=%d) %d (L=%d,m=%d) : [ %.7f  %.7f  %.7f] %n", bf1, LA,m, bf2, LB,m, 0.0, 0.0, DipoleMatrix[bf1][bf2][2]);            
            }
        }

        // x,y-dipoles:
        

        double rhoPrefac = 0.5 * RP_product_LALB / Math.sqrt( (LA+LB+2.0) * (LA+LB) );  // |LA-LB|=1 => LA+LB =/= 0
        double rho_MatrElHalf;

        // In contrast do the definitions used for off-dialonal elements, here Iplus/Iminus 
        // are defined as the integrals of unity-normalized associated Legendre polynomials,
        // without '0.5' factor and without phi-dependent part and its normalization multipliers (e.g. 1/sqrt(Pi))
        // first, set the _non-zero_ elements with m == 0 on either side:
        /*                    
                <L1,1 | x | L2, 0 > = <L1,-1| y | L2, 0 > = Pi/sqrt(Pi * 2*Pi) * Iminus
                <L1,0 | x | L2,1> =   <L1,0 | y | L2,-1>  = Pi/sqrt(Pi * 2*Pi) * Iplus
        Note that /sqrt(...) factors arise here since here Iplus/Iminus include only assoc.Legendre polynomials
        normalization constants rather than the 'full' constants (including phi-part) from ALegNormalizer[][]
        Note also that 1/sqrt(2*Pi) arise only for m==0, whereas for all m!=0 those are 1/sqrt(Pi), and the
        phi-integral equals 0.5
        */
        if (LA>0) { // <L1,1 | x | L2, 0 > and <L1,-1| y | L2, 0 > are available
            int bf1P = RPart_m2bf[rp1][LA + 1];
            int bf1M = RPart_m2bf[rp1][LA - 1];
            bf2 = RPart_m2bf[rp2][LB + 0];
            if (LA>LB)
                rho_MatrElHalf =   Math.sqrt( 2* (LB+2)*(LB+1) ) * rhoPrefac; // <LB+1,1 | rho | LB,0>, '2*' works here as additional sqrt(2) factor
            else
                rho_MatrElHalf = - Math.sqrt( 2* LB*(LB-1) ) * rhoPrefac; // <LB-1,1 | rho | LB,0>, '2*' works here as additional sqrt(2) factor

            DipoleMatrix[bf1P][bf2] [0] = rho_MatrElHalf ; // <LA, 1|x|LB,0>
            DipoleMatrix[bf1M][bf2] [1] = rho_MatrElHalf ; // <LA,-1|y|LB,0>
            // debug:
            //System.out.printf("X: %d (L=%d,m=%d) %d (L=%d,m=%d) : [ %.7f  %.7f  %.7f] %n", bf1P, LA,1, bf2, LB,0, 0.0, 0.0, DipoleMatrix[bf1P][bf2][0]);
            //System.out.printf("Y: %d (L=%d,m=%d) %d (L=%d,m=%d) : [ %.7f  %.7f  %.7f] %n", bf1M, LA,-1, bf2, LB,0, 0.0, 0.0, DipoleMatrix[bf1M][bf2][1]);
            
        }
        if (LB>0) { // <L1,0 | x | L2, 1 > and <L1,0| y | L2,-1 > are available
            bf1 = RPart_m2bf[rp1][LA + 0];
            int bf2P = RPart_m2bf[rp2][LB + 1];
            int bf2M = RPart_m2bf[rp2][LB - 1];
            if (LA>LB)
                rho_MatrElHalf = - Math.sqrt( 2* LB*(LB+1) ) * rhoPrefac; // <LB+1,0 | rho | LB,1>, '2*' works here as additional sqrt(2) factor
            else
                rho_MatrElHalf =   Math.sqrt( 2* LB*(LB+1) ) * rhoPrefac; // <LB-1,0 | rho | LB,1>, '2*' works here as additional sqrt(2) factor

            DipoleMatrix[bf1][bf2P] [0] = rho_MatrElHalf ; // <LA, 1|x|LB,0>
            DipoleMatrix[bf1][bf2M] [1] = rho_MatrElHalf ; // <LA,-1|y|LB,0>
            // debug:
            //System.out.printf("X: %d (L=%d,m=%d) %d (L=%d,m=%d) : [ %.7f  %.7f  %.7f] %n", bf1, LA,0, bf2P, LB,1, 0.0, 0.0, DipoleMatrix[bf1][bf2P][0]);
            //System.out.printf("Y: %d (L=%d,m=%d) %d (L=%d,m=%d) : [ %.7f  %.7f  %.7f] %n", bf1, LA,0, bf2M, LB,-1, 0.0, 0.0, DipoleMatrix[bf1][bf2M][1]);
            
        }


        // Now, all elements with m!=0 on either side:                        
        //double rhoIntType1, rhoIntType2;

        // Do  ( m, m+1)  -dipoles for m=1... : m<=LA, m+1 <= LB
        // and (-m,-m-1)  -dipoles for m=1... : m<=LA, m+1 <= LB,
        // which are both proportional to Iplus: <LA,|m| | rho | LB,|m|+1> 
        // Note: rho matr.el. is _insensitive_ to signs of m-left & m-right due to the properties of assoc. Legendre polynomials !
        // Hence, first, reduce computation to positive m's on both sides, and only then consider if the m on the right is higher
        for (int m = 1; m <=LA && m <= LB-1; ++m) {
            /*rhoIntType1 =   Math.sqrt((LB-m+1)*(LB-m)) * rho_prefac; // <L+1,m|rho|L,m+1>
            rhoIntType2 = - Math.sqrt((LB+m+1)*(LB+m)) * rho_prefac; // <L+1,m|rho|L,m-1>*/
            if (LA>LB)
                rho_MatrElHalf = - Math.sqrt( (LB-m+1)*(LB-m) ) * rhoPrefac; // <LB+1,|m| | rho | LB,|m|+1>
            else
                rho_MatrElHalf =   Math.sqrt( (LB+m+1)*(LB+m) ) * rhoPrefac; // <LB-1,|m| | rho | LB,|m|+1>
            //P2moldenP[LA][LA+m] ,-m, ...
            int bf1P = RPart_m2bf[rp1][LA + m];     // <LA, m |
            int bf1M = RPart_m2bf[rp1][LA - m];     // <LA,-m |
            int bf2P = RPart_m2bf[rp2][LB + m + 1]; // |LB, m+1>
            int bf2M = RPart_m2bf[rp2][LB - m - 1]; // |LB,-m-1>

            // <LA,m|x|LB, m+1> ~ <LA,m|rho|LB,m+1>
            DipoleMatrix[bf1P][bf2P] [0] = rho_MatrElHalf ; // also writes to [bf2][bf1]
            // <LA,m|y|LB,-m-1> ~ <LA,m|rho|LB,m+1>
            DipoleMatrix[bf1P][bf2M] [1] = rho_MatrElHalf ; // also writes to [bf2][bf1]

            // <LA,-m|x|LB,-m-1> ~ <LA,m|rho|LB,m+1>
            DipoleMatrix[bf1M][bf2M] [0] = rho_MatrElHalf ; // also writes to [bf2][bf1]
            // <LA,-m|y|LB, m+1> ~ <LA,m|rho|LB,m+1>
            DipoleMatrix[bf1M][bf2P] [1] = - rho_MatrElHalf ; // also writes to [bf2][bf1]
        }

        // Do  ( m, m-1)  -dipoles for m=1... : m<=LA, m-1 <= LB, and m-1 must not be 0 => start from m=2
        // and (-m,-m+1)  -dipoles for m=1... : m<=LA, -(m-1) >= -LB i.e. m-1 <= LB
        // which are both proportional to Iplus: <LA,|m| | rho | LB,|m|-1>
        for (int m = 2; m <=LA && m <= LB+1; ++m) {
            if (LA>LB)
                rho_MatrElHalf =   Math.sqrt( (LB+m+1)*(LB+m) ) * rhoPrefac; // <LB+1,|m| | rho | LB,|m|-1>
            else
                rho_MatrElHalf = - Math.sqrt( (LB-m+1)*(LB-m) ) * rhoPrefac; // <LB-1,|m| | rho | LB,|m|-1>

            int bf1P = RPart_m2bf[rp1][LA + m];     // <LA, m |
            int bf1M = RPart_m2bf[rp1][LA - m];     // <LA,-m |
            int bf2P = RPart_m2bf[rp2][LB + m - 1]; // |LB, m-1>
            int bf2M = RPart_m2bf[rp2][LB - m + 1]; // |LB,-m+1>

            // <LA,m|x|LB, m-1> ~ <LA,m|rho|LB,m-1>
            DipoleMatrix[bf1P][bf2P] [0] = rho_MatrElHalf ; // also writes to [bf2][bf1]
            // <LA,m|y|LB,-m+1> ~ <LA,m|rho|LB,m-1>
            DipoleMatrix[bf1P][bf2M] [1] = - rho_MatrElHalf ; // also writes to [bf2][bf1]

            // <LA,-m|x|LB,-m+1> ~ <LA,m|rho|LB,m-1>
            DipoleMatrix[bf1M][bf2M] [0] = rho_MatrElHalf ; // also writes to [bf2][bf1]
            // <LA,-m|y|LB, m-1> ~ <LA,m|rho|LB,m-1>
            DipoleMatrix[bf1M][bf2P] [1] = rho_MatrElHalf ; // also writes to [bf2][bf1]
        }                
    }
    //--------------------------------------------------------------------------
    /**
     *  Gets a references to all basis-set-related information from @param molden object
     *  Note that this.ylm_norms2 should be set manually!
     */
    public void ImportBasisFromMolden(MOLDEN_IO molden) {
        basis = molden.Basis;

        centers = molden.Centers;
        RadialParts = molden.RadialParts;
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    private void _Naive(Polynom3D[][] ylm) {
        //Polynom3D[][] ylm;                           // a set of spherical harmonics used by BasicOverlapIntegral() - TODO: eliminate!
                                                        // ylm_norms2 must match these expressions (user's responsibility!)!
        
        //ylm = SphericalHarmonics.Get_Quick_YLM_Norm4PI();// .Get_Quick_YLM();
        //ylm = SphericalHarmonics.Get_Quick_YLM();
        
        ProgTimer t1 = new ProgTimer();
        Jama.Matrix S2 = new Jama.Matrix(basis.length, basis.length);
        for (int bf1=0; bf1<basis.length; ++bf1)
            for (int bf2=bf1; bf2<basis.length; ++bf2) {
                basis[bf1].Quick_YLM = ylm;
                try {
                    double o = basis[bf1].OverlapWith(basis[bf2]);
                    S2.set(bf1, bf2, o);
                    S2.set(bf2, bf1, o);
                } catch (Exception e) {};
            }
        t1.Stop();
        t1.Print();
        out.println("========================");
        //S2.print(13, 7);

        Jama.Matrix S = new Jama.Matrix(this.OverlapMatrix);
        //S.print(13, 7);
        Jama.Matrix DIFF = S2.minus(S);
        double mxdiff = 0;
        int b1=0,b2=0;
        for (int bf1=0; bf1<basis.length; ++bf1)
            for (int bf2=0; bf2<basis.length; ++bf2)
                if (Math.abs(DIFF.get(bf1, bf2)) > mxdiff) {
                    mxdiff = Math.abs(DIFF.get(bf1, bf2)) ;
                    b1 = bf1; b2 = bf2;
                    //out.printf("%.20E:  %.20f v.s. %.20f observed in: %d, %d%n",mxdiff, S2.get(b1, b2), S.get(b1, b2), b1, b2);
                }
        out.printf("%.20E:  %.20f v.s. %.20f observed in: %d, %d%n",mxdiff, S2.get(b1, b2), S.get(b1, b2), b1, b2);
                

    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // DEBUG PIECES GO BELOW
    //--------------------------------------------------------------------------
    {
        try {
           // _test();
        } catch (Exception e) {};
    }
    //--------------------------------------------------------------------------
    /** Quick v.s. Naive methods performance comparison
     */
    public void _test() throws Exception {
        moldenio.MOLDEN_IO molden = new moldenio.MOLDEN_IO();
        molden.Allow_additional_r_power = true;
        String fname = "M:\\abInitio\\JANPA\\tests_my\\conf_462_gRJ_vQz_CART_C2P.molden";
        molden.Load_From_MOLDEN(fname);
        molden.CoordsToAU();


        // convert to a convenient spherical harmonics
        //ylm = SphericalHarmonics.Get_Quick_YLM_Norm4PI();// .Get_Quick_YLM();
        //ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm);

        // ylm_norms2 for unity-normalized YLM's
        double[][] ylm_norms2  = new double[5][];        
        double[][] ylm_norms2_over4Pi  = new double[5][];        
        for (int L=0; L<ylm_norms2.length; L++) {
            ylm_norms2 [L] = new double[2*L+1];
            ylm_norms2_over4Pi[L] = new double[2*L+1];
            for (int m=0; m<2*L+1; m++) {
                ylm_norms2 [L][m] = 1;
                ylm_norms2_over4Pi[L][m] = 1/(4*Math.PI);
            }
        }        
        molden.ToUnnormalizedPrimitiveCoefs();

        //molden.UnNormalizePrimitives(ylm_norms2_over4Pi);
        /*basis = molden.Basis;
        centers = molden.Centers;
        RadialParts = molden.RadialParts;*/
        ImportBasisFromMolden(molden);

        // (spherical harmonic)·(r2^k) generator
        //BasisFunction tmp = new BasisFunction(0,0,null,0);
        //tmp.Quick_YLM = ylm;

        ProgramOptions.ProgTimer t1 = new ProgramOptions.ProgTimer();
        t1.Start();
        Build_Ovarlap_Matrix_PURE();//ylm_norms2_over4Pi);

        t1.Stop();
        t1.Print();

        // Naive method for comparison
        out.println("Naive method...");
        
        molden = new moldenio.MOLDEN_IO();
        molden.Allow_additional_r_power = true;
        molden.Load_From_MOLDEN(fname);        
        molden.CoordsToAU();
        Polynom3D[][] ylm = SphericalHarmonics.Get_Quick_YLM_Norm4PI();// .Get_Quick_YLM();
        ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm);
        molden.UnNormalizePrimitives(ylm_norms2);
        ImportBasisFromMolden(molden);
        _Naive(ylm);
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // result[0] = (nx,ny,0),
    // result[1] = (phi, sin(phi), cos(phi))
    private double[][] rot_vec(double[] Z_new) {
        double[][] result = new double[3][3];
        result[0][0] = Z_new[1]; // [(0,0,1) x nZ]_x = -nZ_Y
        result[0][1] = -Z_new[0];  // [(0,0,1) x nZ]_y = +nZ_X
        result[0][2] = 0;
        double n2 = Math.sqrt(Z_new[1]*Z_new[1] + Z_new[0]*Z_new[0]);
        result[1][1] = n2;
        result[1][2] = Z_new[2]; // cos(J) = // (0,0,1).nZ
        result[1][0] = Math.atan2(n2, Z_new[2]);
        result[0][0] /= n2;
        result[0][1] /= n2;
        return result;
    }
    //--------------------------------------------------------------------------
    private double[] rotate_vec_to_new(double[] r, double[] Z_new){
        double[][] n = rot_vec(Z_new);
        double cosJ = n[1][2];
        double sinJ = n[1][1];

        double[] result = new double[3];

        double nx = n[0][0];
        double ny = n[0][1];
        double an = r[0]*nx + r[1]*ny;

        result[0] = r[0]*cosJ+(1-cosJ)*nx*an-sinJ*(ny*r[2]-0*r[1]);
        result[1] = r[1]*cosJ+(1-cosJ)*ny*an-sinJ*(0*r[0]-nx*r[2]);
        result[2] = r[2]*cosJ+(1-cosJ)*0*an-sinJ*(nx*r[1]-ny*r[0]);
        return result;
    }
    //--------------------------------------------------------------------------
    private double[] rotate_vec_from_new(double[] r, double[] Z_new){
        double[][] n = rot_vec(Z_new);
        double cosJ = n[1][2];
        double sinJ = n[1][1];

        double[] result = new double[3];

        double nx = n[0][0];
        double ny = n[0][1];
        double an = r[0]*nx + r[1]*ny;

        result[0] = r[0]*cosJ+(1-cosJ)*nx*an+sinJ*(ny*r[2]-0*r[1]);
        result[1] = r[1]*cosJ+(1-cosJ)*ny*an+sinJ*(0*r[0]-nx*r[2]);
        result[2] = r[2]*cosJ+(1-cosJ)*0*an+sinJ*(nx*r[1]-ny*r[0]);
        return result;
    }
    //--------------------------------------------------------------------------
    private void print_vec(double[] r) {
        for (int mu=0; mu<3; ++mu) out.printf("%.5f\t", r[mu]);
        out.println();
    }
    //--------------------------------------------------------------------------

    private double[] nZ(double[] RA, double[] RB) {
        double[] result = new double[3];
        double nZ2 = 0;
        for (int mu=0; mu<3; ++mu) { result[mu] = RB[mu] - RA[mu]; nZ2 += result[mu]*result[mu]; }
        nZ2 = Math.sqrt(nZ2);
        for (int mu=0; mu<3; ++mu) result[mu] /= nZ2;
        return result;
    }
    //--------------------------------------------------------------------------
    private double[] ThetaPhi_2_n(double[] ti) {
        double[] n = new double[3];
        n[0] = Math.sin(ti[0])*Math.cos(ti[1]);
        n[1] = Math.sin(ti[0])*Math.sin(ti[1]);
        n[2] = Math.cos(ti[0]);
        return n;
    }
    //--------------------------------------------------------------------------
    /** A FUNCTION FOR DEBUG / TESTS
     *  NOT TO BE USED IN PRACTICE !!!!!!
     */
    private double Overlap(BasisFunction fA, BasisFunction fB) throws Exception{

        for (int mu=0; mu<3; ++mu) { fB.R0[mu] -= fA.R0[mu]; fA.R0[mu] = 0; }

        double result;
        // spherical coordinates of RA->RB w.r.t. LAB system
        double[] ti = this.ThetaPhi(fA.R0, fB.R0);
        // cartesian components of RA->RB w.r.t. LAB system
        double[] eZ = this.ThetaPhi_2_n(ti);

        double[][][] sincos = RealWignerAdapt.create_sincos(ti[0], ti[1]);
        double[][] WignerA = RealWignerAdapt.Wigner(sincos, fA.L);
        double[][] WignerB = RealWignerAdapt.Wigner(sincos, fB.L);

        WignerA = WignerTranspose(WignerA, fA.L);
        WignerB = WignerTranspose(WignerB, fB.L);

        double[] RA_new = this.rotate_vec_from_new(fA.R0, eZ); // TODO: simplify!
        double[] RB_new = this.rotate_vec_from_new(fB.R0, eZ); // TODO: simplify!

        BasisFunction tmpA = new BasisFunction(fA);
        tmpA.R0 = RA_new;
        tmpA.OI = new OverlapIntegrals();
        BasisFunction tmpB = new BasisFunction(fB);
        tmpB.R0 = RB_new;


        int L_MIN = fA.L;
        if (fB.L<L_MIN) L_MIN = fB.L;

        ///////////////
        // ONLY FOR DEBUG:
        /*
        result = 0;
        double[][] prods = new double[2*fA.L+1][2*fB.L+1];
        for (int m1=-fA.L; m1<=fA.L; ++m1) {
            for (int m2=-fB.L; m2<=fB.L; ++m2) {
                tmpA.m = P2moldenP[fA.L][fA.L+m1];
                tmpB.m = P2moldenP[fB.L][fB.L+m2];
                prods[fA.L+m1][fB.L+m2] = tmpA.OverlapWith(tmpB);
                out.printf("%.5f\t", prods[fA.L+m1][fB.L+m2]);
            }
            out.println();
        }
         *
         */


        double[] diags = new double[2*L_MIN+1];
        for (int m=-L_MIN; m<=L_MIN; ++m) {
                tmpA.m = P2moldenP[fA.L][fA.L+m];
                tmpB.m = P2moldenP[fB.L][fB.L+m];
                diags[L_MIN+m] = tmpA.OverlapWith(tmpB);
        }


        // номера функций, которую мы раскладываем, в НОРМАЛЬНОЙ нумерации
        int fixed_mA = fA.L + moldenP2P[fA.L][fA.L+fA.m];
        int fixed_mB = fB.L + moldenP2P[fB.L][fB.L+fB.m];

        result = 0;
        for (int m=-L_MIN; m<=L_MIN; ++m) {
            result += diags[L_MIN+m] * WignerA[fixed_mA][fA.L+m] * WignerB[fixed_mB][fB.L+m];
        }

        return result;
    }
    //--------------------------------------------------------------------------
    /* Some manipulations with coordinate axis change - FOR DEBUG PURPOSES ONLY!!!
     */
    private void ___aaa()throws Exception{
        moldenio.MOLDEN_IO A = new moldenio.MOLDEN_IO();
        //A.Load_From_MOLDEN("Z:\\DNA_software\\twfn\\OpenNAO\\CONVERTERS\\ORCA\\zver1_optimized_spher_noSymm_vdz2_m2m.molden");
        A.Load_From_MOLDEN("Z:\\DNA_software\\twfn\\OpenNAO\\CONVERTERS\\ORCA\\kr3_avtz_spher_ORCA_m2m.molden");
        //npa.dont_print_matrices = false;
        A.CoordsToAU();

        // convert to a convenient spherical harmonics
        Polynom3D[][] ylm = SphericalHarmonics.Get_Quick_YLM();
        double[][] ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm);
        A.UnNormalizePrimitives(ylm_norms2);


        for (int b =0; b<A.Basis.length; ++b)
            out.printf("%d  at center %d  L = %d %n", b, A.Basis[b].Center_ID  ,A.Basis[b].L);

        for (int b1=0; b1<A.Basis.length; ++b1)
            for (int b2=b1+1; b2<A.Basis.length; ++b2) {
                if ((A.Basis[b1].Center_ID != A.Basis[b2].Center_ID) && (A.Basis[b1].L>=1) && (A.Basis[b2].L>=1)) {
                    // perform a test
                    BasisFunction fA = A.Basis[ b1 ]; // 4
                    BasisFunction fB = A.Basis[ b2 ]; // 28

                    fA.Quick_YLM = ylm;
                    fB.Quick_YLM = ylm;
                    fA.OI = new OverlapIntegrals();
                    double exact = fA.OverlapWith(fB);
                    double quick = Overlap(fA, fB);
                    if ((Math.abs((exact-quick)/exact)>1E-12)&&(Math.abs(exact-quick)>1E-15)) {
                        out.printf("EXACT = %.20f   NEW = %.20f,  DIFF = %.3E%n", exact, quick, exact-quick);
                        //out.printf("NEW = %n%.20f%n", quick);
                    }
                }

            }
        BasisFunction fA = A.Basis[ 82 ]; // 4
        BasisFunction fB = A.Basis[ 28 ]; // 28

            fA.OI = new OverlapIntegrals();
            out.printf("EXACT = %n%.20f%n", fA.OverlapWith(fB));

            out.printf("NEW = %n%.20f%n", this.Overlap(fA, fB));


        fA.Quick_YLM = ylm;
        fB.Quick_YLM = ylm;

        double[] RA = fA.R0;
        double[] RB = fB.R0;

        double[] ti = ThetaPhi(RA,RB);
        double[] _nZ = ThetaPhi_2_n(ti);//nZ(RA,RB);
        //double[] _nZ2 = ThetaPhi_2_n(ThetaPhi(RB,RA));

        double[] r = new double[]{RA[0]+0.2*RB[0], RA[1]+0.3*RB[1], RA[2]+0.3*RB[2] };
        //double[] r = new double[]{RB[0]-RA[0], RB[1]-RA[1], RB[2]-RA[2] };

        //r = _nZ;
        /*out.print("ez'(theta,phi_) = "); print_vec(_nZ);
        out.print("TRUE ez' = "); print_vec(this.nZ(RA, RB));

        print_vec(r);*/
        /*print_vec(r_in_new);
        double[] x = rotate_vec_from_new(r_in_new, _nZ);
        print_vec(x);
         *
         */

        //Polynom3D[][] ylm = SphericalHarmonics.Get_Quick_YLM_Norm4PI();// .Get_Quick_YLM();
        // WIGNER CHECK
        int nPassed = 0;
        for (int N=0; N<10000; ++N) {

            r[0] = Math.random()*10; r[1] = Math.random()*10; r[2] = Math.random()*10;

            ti = new double[]{Math.random()*Math.PI, Math.random()*2*Math.PI} ;
            _nZ = ThetaPhi_2_n(ti);

            double[] r_in_new = rotate_vec_to_new(r, _nZ);

            final int _L_MAX = 2;
            int L = (int)(Math.random()*_L_MAX)/*0...3*/+1;
            int M = (int)(Math.random()*(L+1))/*0...L*/;
            if (Math.random() > 0.5) M = -M;
            M = M+L; // make it 0-based
            double ref = ylm[L][M/*L+P2moldenP[L][M]*/].EvaluateAtPoint(r);

            double[][] WD = RealWignerAdapt.Wigner(RealWignerAdapt.create_sincos(ti[0], ti[1]), L);
            WD = Wigner2MOLDEN(WD, L);

            double v = 0;
            for (int m=0; m<(2*L+1); ++m)
                v += WD[M][m]*ylm[L][m/*L+P2moldenP[L][m]*/].EvaluateAtPoint(r_in_new);

            if (Math.abs((ref-v)/ref) > 1E-12) {
                out.printf("L = %d; M = %d; r[0] = %.20f; r[1] = %.20f; r[2] = %.20f;%n // ref = %.20f v = %.20f%n",
                        L,M,r[0],r[1],r[2], ref, v);
            } else ++nPassed;
        }
        out.println("PASSED = "+nPassed);
            // rotate r





        //out.printf("%.20f%n", ylm[L][M].EvaluateAtPoint(r_in_new) );
        //out.printf("%.20f%n", ylm[L][M].EvaluateAtPoint(rotate_vec(r, _nZ2)) );


        out.println();

    }
    //--------------------------------------------------------------------------
    private void _test1() {    
        
        /*moldenio.MOLDEN_IO A = new moldenio.MOLDEN_IO();
        try {
        A.Load_From_MOLDEN("M:\\abInitio\\JANPA\\tests_my\\conf_462_gRJ_vQz_CART_C2P.molden");
        } catch(Exception e) {e.printStackTrace(); }
        A.CoordsToAU();*/
        // convert to a convenient spherical harmonics
        Polynom3D[][] ylm = SphericalHarmonics.Get_Quick_YLM();
        double[][] ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm);
        /*A.UnNormalizePrimitives(ylm_norms2);
        for (int b =0; b<A.Basis.length; ++b)
            out.printf("%d  at center %d  L = %d %n", b, A.Basis[b].Center_ID  ,A.Basis[b].L);*/

        int L1 = 2; int L2 = 4;
        int m1 = -1; int m2 = m1;
        double zeta1 = 1.0; double zeta2 = 1.0;
        double z1 = -0.3; 
        double z2 = 0.5;
        //double z2 = z1;
        
        // canonical: y,z,x -> molden: z,x,y
        int[] p_canonical2molden = new int[]{1,-1,0};
        int m1_molden = (L1==1)?p_canonical2molden[L1+m1]:m1;
        int m2_molden = (L2==1)?p_canonical2molden[L2+m2]:m2;
        BasisFunction fA = new BasisFunction(L1, m1_molden, new double[]{0,0,z1}, 1); //A.Basis[ b1 ]; // 4
        fA.coefs[0] = 1.0;
        fA.exponents[0] = zeta1;
        BasisFunction fB = new BasisFunction(L2, m2_molden, new double[]{0,0,z2}, 1);  //A.Basis[ b2 ]; // 28
        fB.coefs[0] = 1.0;
        fB.exponents[0] = zeta2;

        fA.Quick_YLM = ylm;
        fB.Quick_YLM = ylm;
        fA.OI = new OverlapIntegrals();
        fB.OI = new OverlapIntegrals();
        ylm[L1][L1+m1_molden].Print();
        double exact=0;
        double fA_nrm2 = 0; double fB_nrm2 = 0; 
        double r2Int_1 = OverlapIntegrals.primitive_int_1D_Sphr(2 + (L1+L1), (zeta1+zeta1) );
        double r2Int_2 = OverlapIntegrals.primitive_int_1D_Sphr(2 + (L2+L2), (zeta2+zeta2) );
        try {
          exact = fA.OverlapWith(fB);
          fB_nrm2 = fB.OverlapWith(fB) / r2Int_1; // leave only shp.harm. norm2
          fA_nrm2 = fA.OverlapWith(fA) / r2Int_2; // leave only shp.harm. norm2
        } catch(Exception e) {e.printStackTrace(); }
        //System.out.printf("exact = %.12f %n", exact );
        exact /=  Math.sqrt(fA_nrm2*fB_nrm2) ;
        System.out.printf("exact = %.12f %n", exact );
        
        //double nrm2_1 = 
        
        double newVal = 0;//BasicOverlapInt(zeta1, L1, z1, 0, zeta2, L2, z2, 0, m1);
        System.out.printf("  new = %.12f %n", newVal );
        
/*        
        int Lmax = 4;
        JGints rf = new JGints();
        for (int L1=0; L1<=Lmax; L1++) {
            for (int L2=0; L2<=Lmax; L2++) {
                double ref = JGints
            }
        }
        JGintsCyl a = new JGintsCyl();
        double zeta1 = 1.0;
        double zeta2 = 1.0;
        double z1 = 0.0, z2 = 0.0;
        int L = 4;
        for (int m=0; m<=L; m++) {
            double nrmMul = Polynom_rho_z.ALegNormalizer[L][L+m];
            
            // some pre-computed values of INT( (r^L * exp(-1.0*r^2))^2 * r^2, r=0..+infinity):
            double r2Int = OverlapIntegrals.primitive_int_1D_Sphr(2 + (L+L), (zeta1+zeta2) );
            
            System.out.printf("%.15f %n",                
                a.BasicOverlapInt(zeta1, L, z1, 0, zeta1, L, z2, 0, m) * nrmMul*nrmMul / r2Int
            );
        }
        */
    }
    
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        JGintsCyl a = new JGintsCyl();
        //JGints b = new JGints();
        a._test();
        
        if(true)
            return;
        
        //a._test1();
        
        double zeta1 = 1.0;
        double zeta2 = 1.0;
        double z1 = 0.0, z2 = 0.0;
        int L = 4;
        double[][] results = new double[a.PrimitiveGOverlapInt_resultsLEN][];
        a.PrimitiveGOverlapInt(zeta1, L, z1, 0, zeta1, L, z2, 0, L, results, true) ;
        
        for (int m=0; m<=L; m++) {
            double nrmMul = SphericalHarmonics.ALegNormalizer[L][L+m];
            
            // some pre-computed values of INT( (r^L * exp(-1.0*r^2))^2 * r^2, r=0..+infinity):
            double r2Int = OverlapIntegrals.primitive_int_1D_Sphr(2 + (L+L), (zeta1+zeta2) );
            //System.out.printf("%.5f%n", r2Int);            
            //double[] r2Int = new double[]{0.1566642671 /*L==0*/, 0.1174982003 /*L==1*/,
            //    0.1468727504 /*L==2*/, 0.2570273133 /*L==3*/, 0.5783114548 /*L==4*/, 1.590356501 /*L==5*/ };
            
            System.out.printf("%.15f %n",                
                results[0][L+m] / r2Int
            );
        }
        // 0,0 -> 1    * 4*Pi
        // 1,0 -> 1/4  * 4*Pi
        // 2,0 -> 1/12 * 4*Pi
        // 3,0 -> 3/80 * 4*Pi
        // 
    }

}
