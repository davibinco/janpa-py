package JGints;

import java.io.*;
import Polynom3D.*;
import ProgramOptions.ProgTimer;
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
 * (c) Tymofii Nikolaienko, 2014
 */


/**
 * Some molecular integrals over Gaussian-type basis functions implemented with JAVA
 * @author (c) TimN, 2014
 */
public class JGints {
    public BasisFunction[] basis;                       // basis functions to be integrated
    public RadialPartOfBasisFunction[] RadialParts;     // expansion of radial parts of basis functions indexed by basis[].RadialPartID over primitive gaussians
    public AtomicCenter[] centers;                      // coordinates of atoms referenced by basis[].Center_ID-1
    public double[][] ylm_norms2;                       // ||norms||2 of spherical harmonics to be taken
    public Polynom3D[][] ylm;                           // a set of spherical harmonics used by BasicOverlapIntegral() - TODO: eliminate!
                                                        // ylm_norms2 must match these expressions (user's responsibility!)!

    public double[][] OverlapMatrix;                    // overlap integrals matrix basis.length x basis.length

    PrintStream out = System.out;

    //--------------------------------------------------------------------------
    // result[0] = theta = t, result[1] = phi = i
    private double[] ThetaPhi(double[] RA, double[] RB) {
        double[] dr = new double[3];
        for (int mu=0; mu<3; ++mu)
            dr[mu] = RB[mu] - RA[mu];
        double[] result = new double[2];
        result[0] = Math.atan2(Math.sqrt(dr[1]*dr[1]+dr[0]*dr[0]), dr[2]/*dz*/  );
        result[1] = Math.atan2(dr[1]/*dy*/, dr[0]/*dx*/);
        return result;
    }
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
    private BasisFunction tmp = new BasisFunction(0,0,null,0);
    private Polynom3D YLM1, YLM2;
    private OverlapIntegrals OI = new OverlapIntegrals(); // is not used in the quickest version!
    private NestedFunctionsDerivs derivComputer = new NestedFunctionsDerivs();

    /**
     *                    /
     * Evaluates  INT  = | dxdydz · Y[L1][m](x,y,z-ZA)·(x^2+y^2+(z-ZA)^2)^(addit_r_pwr1/2)·
     *                  /           Y[L2][m](x,y,z-ZB)·(x^2+y^2+(z-ZB)^2)^(addit_r_pwr2/2)·
     *                              exp(-zeta1·(x^2+y^2+(z-ZA)^2))·exp(-zeta2·(x^2+y^2+(z-ZB)^2))
     *
     */
    public double BasicOverlapInt(double zeta1, int L1, double ZA, int addit_r_pwr1,
                           double zeta2, int L2, double ZB, int addit_r_pwr2,
                           int m) {
        double result = 0;
        double[] ra_new = new double[]{0,0,0};
        double[] rb_new = new double[]{0,0,0};

        ra_new[2] = ZA;
        rb_new[2] = ZB;
        try {
            boolean _dbg = false;
            if (_dbg) {
                // THE SIMPLEST WAY
                tmp.additional_r_power = addit_r_pwr1;
                tmp.L = L1;
                tmp.m = P2moldenP[L1][L1+m];
                YLM1 = tmp.Get_YLM_with_Additional_r_power(tmp, false);

                tmp.additional_r_power = addit_r_pwr2;
                tmp.L = L2;
                tmp.m = P2moldenP[L2][L2+m]; // note that m is the same !!!
                YLM2 = tmp.Get_YLM_with_Additional_r_power(tmp, false);

                // shift Y[L1,neccess_m]*r2^n and Y[L2,neccess_m]*r2^m to a 'middle point'
                // evaluate neccessary integrals and store them
            }
            if (_dbg) 
                result = OI.BS_BS_Overlap(zeta1, YLM1, ra_new, zeta2, YLM2, rb_new);
                            
            

            // MORE ADVANCED WAY!
            
            int max_r_pwr = L1 + L2 + 2*addit_r_pwr1 + 2*addit_r_pwr2; // dimensionalyty of the result

            // use the gaussian product theorem
            double zeta_eff = zeta1 + zeta2;
            double z0 = (zeta1*ZA + zeta2*ZB) / zeta_eff;
            double Prefactor = 0;
            Prefactor = Math.exp(-zeta1*zeta2*(ZB-ZA)*(ZB-ZA) / zeta_eff);


            // use those radial parts while computing angular overlaps
            double[][] AngularPart;
            //AngularInts.radial = new double[max_r_pwr+1];
            double sqrtZetaEff = Math.sqrt(zeta_eff);
            //int highest_Deriv =  + addit_r_pwr2/2; // w.r.t. zeta1 AND zeta2 (d^(highest_Deriv) / d zeta1^i / d zeta2 ^j, i+j=highest_Deriv)
            AngularPart = AngularInts.evaluate(L1, L2, m,
                    addit_r_pwr1/2, // N_zeta1_derivs
                    addit_r_pwr2/2, // N_zeta2_derivs
                    zeta1, zeta2, ZB-ZA);

            double resultPrelim = 0;

            resultPrelim = AngularPart[0][0]* Prefactor;// * Math.pow(zeta_eff, -L1/2.0 - L2/2.0 - 3/2.0);

            if (_dbg && (Math.abs(result/resultPrelim-1)>1E-10) && (addit_r_pwr1 + addit_r_pwr2==0))
                out.printf("TRUE = %.20E,  NEW = %.20E%n",result,resultPrelim);
                

            // Now Result[0][0]*Prefactor is a proper overlap integral in case if (addit_r_pwr1 == 0) && (addit_r_pwr2 == 0))

            // If addit_r_pwr1 > 0 or addit_r_pwr2 > 0, we can reduce to the previous case if we differentiate
            // that prev. case w.r.t. zeta1 and/or zeta2
            if (addit_r_pwr1 > 0 || addit_r_pwr2 > 0) {
                    // 1) compute all necessary derivatives of the Prefactor
                    //  1.1) compute all necessary derivatives of -zeta1*zeta2*DZ^2/(zeta1+zeta2)
                        double[][] expoDerivs = derivComputer.prefactorExponentDerivs(
                                                 addit_r_pwr1/2, addit_r_pwr2/2, ZB-ZA, zeta1, zeta2);
                    //  1.2) compute all necessary derivatives of exp(-zeta1*zeta2*DZ^2/(zeta1+zeta2))
                        double[][] prefactorDerivs = derivComputer.exponentialDerivs(expoDerivs, addit_r_pwr1/2, addit_r_pwr2/2);
                        /*
                    //  1.3) compute derivatives of zetaEff^(-L1/2-L2/2-3/2)
                        double[] zEff_drws = new double[2+addit_r_pwr1/2+addit_r_pwr2/2];
                        zEff_drws[0] = Math.pow(zeta_eff, -L1/2.0 - L2/2.0 - 3/2.0);
                        for (int i=1; i<zEff_drws.length; ++i)
                            zEff_drws[i] = zEff_drws[i-1]/zeta_eff * (-L1/2.0-L2/2.0-3.0/2-(i-1));
                        double[][] z12_drv = new double[1+addit_r_pwr1/2][1+addit_r_pwr2/2];
                        for (int i=0; i<z12_drv.length; ++i)
                            for (int j=0; j<z12_drv[i].length; ++j)
                                z12_drv[i][j] = zEff_drws[i+j];
                         * 
                         */
                    // 2) compute all necessary derivatives w.r.t. zeta1 and zeta2 of the spherical harmonic integral
                    //  Here we do like this:
                      //       /
                      //  I = | Y[L1][m](x,y,z-sqrt(zetaEff)*(ZA-z0))*Y[L2][m](x,y,z-sqrt(zetaEff)*(ZB-z0))*exp(-1.0*r^2)
                      //     /                 ~~~~~~~~~~~a~~~~~~~~~                 ~~~~~~~~~~~b~~~~~~~~~~
                    // already done!
                    // 3) combine results!
                        //double[][] prefactor1 = derivComputer.productDerivs(prefactorDerivs, z12_drv, 1.0, addit_r_pwr1/2, addit_r_pwr2/2);
                        double[][] finalDeriv = derivComputer.productDerivs(prefactorDerivs/*prefactor1*/, AngularPart, 1.0, addit_r_pwr1/2, addit_r_pwr2/2);

                        resultPrelim = finalDeriv[addit_r_pwr1/2][addit_r_pwr2/2];
                        // note: prefactorDerivs[0][0] already includes prefactor!!!
                        
                        if ((addit_r_pwr1/2+addit_r_pwr2/2)%2 == 1)
                            resultPrelim = -resultPrelim; // note (-1)^(n1+n2) which appears before derivatives
                        
                        //out.printf("%.20E v.s. %.20E", result,resultPrelim);
                    if (_dbg && (Math.abs(result/resultPrelim-1)>1E-10)  ) {
                        out.printf("TRUE = %.20E,  NEW = %.20E -- L1 = %d, L2 = %d ,m = %d ",result,resultPrelim,L1,L2,m);
                        out.printf("--- d1 = %f,  d2 = %f %n",addit_r_pwr1/2.0, addit_r_pwr2/2.0);
                    } else {
                        //out.printf("OK: L1 = %d, L2 = %d ,m = %d ",L1, L2, m);
                        //out.printf("--- d1 = %f,  d2 = %f %n",addit_r_pwr1/2.0, addit_r_pwr2/2.0);
                    }

                        

                    //resultPrelim = OI.BS_BS_Overlap(zeta1, YLM1, ra_new, zeta2, YLM2, rb_new);
            }


            result = resultPrelim;
            
        } catch (Exception e) {
            out.println("Exception in BasicOverlapInt: "+e.toString());
            e.printStackTrace();
        };
        return result;
    }
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    /**
     * MAIN METHOD FOR OVERLAP INTEGRAL MATRIX EVALUATION
     * Quick calculation of overlap integral matrix for basis functions with PURE spherical harmonics
     * 
     */
    public void Build_Ovarlap_Matrix_PURE() {

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
        double result;

        // Pre-compute zetaEffCoefs[][], Cnk[][] and Ank[][] arrays for AngularInts.(...):
        int max_deriv = max_r2_pwr + 1; // '+1' is for safe ))
        // prepare kappa[n][k]
        AngularInts.zetaEffCoefs = NestedFunctionsDerivs.xSqrtNInvDerivs(2*max_L+2*max_L+3 , max_deriv + max_deriv);
        // prepare C_nk
        AngularInts.Cnk = Polynom3D.ensure_Cnk_enough(max_deriv);
        // prepare A_nk
        AngularInts.Ank = NestedFunctionsDerivs.ensure_Ank_enough(max_L);

        
        tmp.Quick_YLM = ylm;   // needed by [OLD VERSION OF] BasicOverlapIntegral()

        // The place to store result:
        OverlapMatrix = new double[basis.length][basis.length];

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
                        // alloc mem
                        double[] Diags = new double[2*L_MIN+1];

                        // now loop over ALL pairs of exponents in order to get all 'middle points'
                        for (int expon1=0; expon1<RadialParts[rp1].Exponents.length; ++expon1) {
                            zeta1 = RadialParts[rp1].Exponents[expon1];
                            cf1 = RadialParts[rp1].Coefs[expon1];
                            //diags[expon1] = new double[molden.RadialParts[rp2].Exponents.length][];
                            for (int expon2=0; expon2<RadialParts[rp2].Exponents.length; ++expon2) {
                                zeta2 = RadialParts[rp2].Exponents[expon2];
                                //diags[expon1][expon2] = new double[2*L_MIN+1];
                                // calculate a middle point
                                // calculate prefactors and save them
                                // get the range of neccessary values of m (shared!)
                                for (int m=-L_MIN; m<=L_MIN; ++m) {
                                        // evaluate neccessary integrals and store them
                                        /*if (false) {
                                            // newer version:
                                            Diags[L_MIN+m] += JGintsCyl.BasicOverlapInt(zeta1, LA, RA_new[2], RadialParts[rp1].Addit_r_power,
                                                                              zeta2, LB, RB_new[2], RadialParts[rp2].Addit_r_power,
                                                                              m) *
                                                    cf1 * RadialParts[rp2].Coefs[expon2];
                                        } else {*/
                                            // default version                                        
                                            Diags[L_MIN+m] += BasicOverlapInt(zeta1, LA, RA_new[2], RadialParts[rp1].Addit_r_power,
                                                                              zeta2, LB, RB_new[2], RadialParts[rp2].Addit_r_power,
                                                                              m) *
                                                    cf1 * RadialParts[rp2].Coefs[expon2];
                                        //}
                                        
                                }
                            }
                        }
                        // use stored integrals to produce ALL integrals involving these radial parts

                        // loop over all basis functions which use these radial parts
                        for (int b1I=0; b1I<BFNsOf_RP[rp1].length; ++b1I) {
                            bf1 = BFNsOf_RP[rp1][b1I];
                            // номера функций, которую мы раскладываем, в НОРМАЛЬНОЙ m-нумерации
                            int fixed_mA = LA + moldenP2P[LA][LA+basis[bf1].m];
                            for (int b2I=0; b2I<BFNsOf_RP[rp2].length; ++b2I) {
                                bf2 = BFNsOf_RP[rp2][b2I];
                                // номера функций, которую мы раскладываем, в НОРМАЛЬНОЙ m-нумерации
                                int fixed_mB = LB + moldenP2P[LB][LB+basis[bf2].m];

                                result = 0;
                                for (int m=-L_MIN; m<=L_MIN; ++m) {
                                    result += Diags[L_MIN+m] * W[LA][fixed_mA][LA+m] * W[LB][fixed_mB][LB+m];
                                }

                                OverlapMatrix[bf1][bf2] = result;
                                OverlapMatrix[bf2][bf1] = result;
                            }
                        } // for b1I
                        Diags = null;
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
                    if (RadialParts[rp2].LUsedWith != LA)
                        continue; // go to the next iteration; harmonics with different L will be orthogonal!

                    addit_r_pwr2 = RadialParts[rp2].Addit_r_power;

                    RP_product = 0;
                    for (int c1=0; c1<RadialParts[rp1].Coefs.length; ++c1) {
                        zeta1 = RadialParts[rp1].Exponents[c1];
                        cf1 = RadialParts[rp1].Coefs[c1];
                        for (int c2=0; c2<RadialParts[rp2].Coefs.length; ++c2) {
                            cf2 = RadialParts[rp2].Coefs[c2];
                            RP_product += OverlapIntegrals.primitive_int_1D_Sphr(
                                    2/*·r^2*/+2*LA/*(r^l)^2*/+addit_r_pwr1+addit_r_pwr2/*(r^k1)*(r^k2)*/,
                                    zeta1 + RadialParts[rp2].Exponents[c2]) * cf1*cf2;
                        }
                    }
                    // use these results to produce same-m overlaps for this radial part of this center

                    // loop over ALL basis functions which use these radial parts
                    for (int bI1=0; bI1<BFNsOf_RP[rp1].length; ++bI1) {
                        bf1 = BFNsOf_RP[rp1][bI1];
                        int m1 = basis[bf1].m;
                        for (int bI2=0; bI2<BFNsOf_RP[rp2].length; ++bI2) {
                            // Ls are the same for rp1 and rp2 - we've already checked this!
                            // are m's equal?
                            bf2 = BFNsOf_RP[rp2][bI2];
                            if (basis[bf2].m == m1) {
                                result = RP_product * ylm_norms2[LA][LA + m1] * Math.PI*4;
                                OverlapMatrix[bf1][bf2] = result;
                                OverlapMatrix[bf2][bf1] = result;
                            }
                        }
                    }
                    // and all the other terms are zeros!
                    
                } // for 2-nd radial part of atom A
            } // for 1-st radial part of atom A
        } // for atom A

        // READY!
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
                }
        out.printf("%.20E:  %.20f v.s. %.20f%n",mxdiff, S2.get(b1, b2), S.get(b1, b2));
                

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
        molden.Load_From_MOLDEN("M:\\abInitio\\JANPA\\tests_my\\conf_462_gRJ_vQz_CART_C2P.molden");
        molden.CoordsToAU();


        // convert to a convenient spherical harmonics
        ylm = SphericalHarmonics.Get_Quick_YLM_Norm4PI();// .Get_Quick_YLM();
        ylm_norms2 = SphericalHarmonics.Get_Quick_YLM_Norm2(ylm);
        molden.UnNormalizePrimitives(ylm_norms2);
        /*basis = molden.Basis;
        centers = molden.Centers;
        RadialParts = molden.RadialParts;*/
        ImportBasisFromMolden(molden);

        // (spherical harmonic)·(r2^k) generator
        //BasisFunction tmp = new BasisFunction(0,0,null,0);
        tmp.Quick_YLM = ylm;

        ProgramOptions.ProgTimer t1 = new ProgramOptions.ProgTimer();
        t1.Start();
        Build_Ovarlap_Matrix_PURE();

        t1.Stop();
        t1.Print();

        // Naive method for comparison
        out.println("Naive method...");
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

}
