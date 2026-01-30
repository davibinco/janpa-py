package JGints;

import Polynom3D.*;

/**
 * A file implementing a SphericalHarmonics class which contains cartesian
 * representation of spherical harmonics suitable for using
 * with MOLDEN-compatible data.
 *
 * Version: 16.Nov.2013
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

/**
 *
 * @author (c) Tymofii Nikolaienko, 2014
 */
public class SphericalHarmonics {
    //--------------------------------------------------------------------------
    final public static int L_MAX = 4;

    // some constants for internal use
    private static double Sqrt3 = Math.sqrt(3);
    private static double Sqrt15 = Math.sqrt(15);
    private static double Sqrt10 = Math.sqrt(10);
    private static double Sqrt6 = Math.sqrt(6);
    private static double Sqrt7 = Math.sqrt(7);
    private static double Sqrt35 = Math.sqrt(35);
    private static double Sqrt350 = Math.sqrt(350);
    public final static double Sqrt2 = Math.sqrt(2);
    private static double Sqrt14 = Math.sqrt(14);    
    //--------------------------------------------------------------------------
    /** @return  Returns cartesian representation of 'pure' spherical functions
     * with L=0...L_MAX, which are convenient for using in conjunction with MOLDEN format
     * 
     */
    public static Polynom3D[][] Get_Quick_YLM() {
       /*
        * NOTE: The pure spherical harmonics defined below coincide with the ones
        * used by MOLDEN program up to a normalization constant.
        * Our harmonics have the following norms:
        *  || Y_lm ||^2 = 4*Pi     for L = 0 (S)
        *  || Y_lm ||^2 = 4*Pi/3   for L = 1 (P)
        *  || Y_lm ||^2 = 4*Pi/15  for L = 2 (D)
        *  || Y_lm ||^2 = 4*Pi/105 for L = 3 (F)
        *  || Y_lm ||^2 = 4*Pi/315 for L = 4 (G)
        * Such a normalization seems to be the most natural for being used with
        * molden files produced by ORCA package (as of ver. 3.0)
        * 
        * At the same time, different ab initio packages use slightly different sign conventions:
        * * PSI4, MOLDEN (and, probably, GAUSSIAN):
        *   The same as defined below
        * 
        * * NWCHEM:
        *   D(+1)_NWCHEM = - D(+1)_MOLDEN
        *   F(+1)_NWCHEM = - F(+1)_MOLDEN
        *   F(+3)_NWCHEM = - F(+3)_MOLDEN
        *   G(+1)_NWCHEM = - G(+1)_MOLDEN
        *   G(+3)_NWCHEM = - G(+3)_MOLDEN
        *
        * * ORCA:
        *   F(+3)_ORCA = - F(+3)_MOLDEN
        *   F(-3)_ORCA = - F(-3)_MOLDEN
        *   G(+3)_ORCA = - G(+3)_MOLDEN
        *   G(-3)_ORCA = - G(-3)_MOLDEN
        *   G(+4)_ORCA = - G(+4)_MOLDEN
        *   G(-4)_ORCA = - G(-4)_MOLDEN
        *
        */
       // Allocate memory
       Polynom3D[][] Quick_YLM = new Polynom3D[5][];
       // 1S
       Quick_YLM[ 0] = new Polynom3D[1];
       Quick_YLM[ 0][ 0] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 0, 0, 0}}); // S
       // 3P
       Quick_YLM[ 1] = new Polynom3D[3];
       Quick_YLM[ 1][ 0] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 0, 0, 1}}); // P(-1) = PZ
       Quick_YLM[ 1][ 1] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 1, 0, 0}}); // P( 0) = PX
       Quick_YLM[ 1][ 2] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 0, 1, 0}}); // P(+1) = PY
       // 5D
       Quick_YLM[ 2] = new Polynom3D[5];
       Quick_YLM[ 2][ 0] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 1, 1, 0}}); // D(-2) = XY
       Quick_YLM[ 2][ 1] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 0, 1, 1}}); // D(-1) = YZ
       Quick_YLM[ 2][ 2] = new Polynom3D(new double[]{   1.000/Sqrt3,  -0.500/Sqrt3,  -0.500/Sqrt3}, new int[][]{{ 0, 0, 2},{ 0, 2, 0},{ 2, 0, 0}}); // D( 0) = (Z^2 - 0.5*X^2 - 0.5*Y^2)/SQRT(3)
       Quick_YLM[ 2][ 3] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 1, 0, 1}}); // D(+1) = XZ
       Quick_YLM[ 2][ 4] = new Polynom3D(new double[]{   0.500,  -0.500}, new int[][]{{ 2, 0, 0},{ 0, 2, 0}}); // D(+2) = (X^2 - Y^2)/2
       // 7F
       Quick_YLM[ 3] = new Polynom3D[7];
       Quick_YLM[ 3][ 0] = new Polynom3D(new double[]{   1.500/Sqrt6,  -0.500/Sqrt6}, new int[][]{{ 2, 1, 0},{ 0, 3, 0}}); // F(-3) = y*(1.5*x^2-0.5*y^2)/SQRT(6)
       Quick_YLM[ 3][ 1] = new Polynom3D(new double[]{   1.000}, new int[][]{{ 1, 1, 1}}); // F(-2) = x*y*z
       Quick_YLM[ 3][ 2] = new Polynom3D(new double[]{   2.000/Sqrt10,  -0.500/Sqrt10,  -0.500/Sqrt10}, new int[][]{{ 0, 1, 2},{ 0, 3, 0},{ 2, 1, 0}}); // F(-1) = y*(2*z^2-0.5*y^2-0.5*x^2)/SQRT(10)
       Quick_YLM[ 3][ 3] = new Polynom3D(new double[]{   1.000/Sqrt15,  -1.500/Sqrt15,  -1.500/Sqrt15}, new int[][]{{ 0, 0, 3},{ 0, 2, 1},{ 2, 0, 1}}); // F( 0) = z*(z^2-1.5*y^2-1.5*x^2)/SQRT(15)
       Quick_YLM[ 3][ 4] = new Polynom3D(new double[]{   2.000/Sqrt10,  -0.500/Sqrt10,  -0.500/Sqrt10}, new int[][]{{ 1, 0, 2},{ 1, 2, 0},{ 3, 0, 0}}); // F(+1) = x*(2*z^2-0.5*y^2-0.5*x^2)/SQRT(10)
       Quick_YLM[ 3][ 5] = new Polynom3D(new double[]{   0.500,  -0.500}, new int[][]{{ 2, 0, 1},{ 0, 2, 1}}); // F(+2) = z*(0.5*x^2-0.5*y^2)
       Quick_YLM[ 3][ 6] = new Polynom3D(new double[]{   0.500/Sqrt6,  -1.500/Sqrt6}, new int[][]{{ 3, 0, 0},{ 1, 2, 0}}); // F(+3) = x*(0.5*x^2-1.5*y^2)/SQRT(6)
       // 9G
       Quick_YLM[ 4] = new Polynom3D[9];
       Quick_YLM[ 4][ 0] = new Polynom3D(new double[]{   0.500, -0.500}, new int[][]{{ 3, 1, 0},{ 1, 3, 0}}); // G(-4) = 0.5*x*y*(x^2-y^2) ~ S44
       Quick_YLM[ 4][ 1] = new Polynom3D(new double[]{   1.500/Sqrt2, -0.500/Sqrt2}, new int[][]{{ 2, 1, 1},{ 0, 3, 1}}); // G(-3) = 0.5*y*z*(3*x^2-y^2)/SQRT(2) ~ S43
       Quick_YLM[ 4][ 2] = new Polynom3D(new double[]{   3.000/Sqrt7,  -0.500/Sqrt7,  -0.500/Sqrt7}, new int[][]{{ 1, 1, 2},{ 1, 3, 0},{ 3, 1, 0}}); // G(-2) = x*y*(3*z^2-0.5*y^2-0.5*x^2)/SQRT(7) ~ S42
       Quick_YLM[ 4][ 3] = new Polynom3D(new double[]{  10.000/Sqrt350,  -7.500/Sqrt350,  -7.500/Sqrt350}, new int[][]{{ 0, 1, 3},{ 0, 3, 1},{ 2, 1, 1}}); // G(-1) = y*z*(10*z^2-7.5*y^2-7.5*x^2)/SQRT(350) ~ S41
       Quick_YLM[ 4][ 4] = new Polynom3D(new double[]{   1.000/Sqrt35,  -3.000/Sqrt35,  -3.000/Sqrt35,   0.3750/Sqrt35,   0.750/Sqrt35,   0.3750/Sqrt35}, new int[][]{{ 0, 0, 4},{ 0, 2, 2},{ 2, 0, 2},{ 0, 4, 0},{ 2, 2, 0},{ 4, 0, 0}}); // G( 0) = (z^4-3*y^2*z^2-3*x^2*z^2+3/8*y^4+3/4*x^2*y^2+3/8*x^4)/SQRT(35) ~ C40
       Quick_YLM[ 4][ 5] = new Polynom3D(new double[]{  10.000/Sqrt350,  -7.500/Sqrt350,  -7.500/Sqrt350}, new int[][]{{ 1, 0, 3},{ 1, 2, 1},{ 3, 0, 1}}); // G(+1) = x*z*(10*z^2-7.5*y^2-7.5*x^2)/SQRT(350) ~ C41
       Quick_YLM[ 4][ 6] = new Polynom3D(new double[]{   1.500/Sqrt7,  -1.500/Sqrt7, 0.250/Sqrt7,  -0.250/Sqrt7}, new int[][]{{ 2, 0, 2},{ 0, 2, 2},{ 0, 4, 0},{ 4, 0, 0}}); // G(+2) = (1.5*x^2*z^2-1.5*y^2*z^2 + 0.25*y^4-0.25*x^4)/SQRT(7) ~ C42
       Quick_YLM[ 4][ 7] = new Polynom3D(new double[]{   0.500/Sqrt2,  -1.500/Sqrt2}, new int[][]{{ 3, 0, 1},{ 1, 2, 1}}); // G(+3) = x*z*(0.5*x^2-1.5*y^2)/SQRT(2) ~ C43
       Quick_YLM[ 4][ 8] = new Polynom3D(new double[]{   0.125,  -0.750,   0.125}, new int[][]{{ 4, 0, 0},{ 2, 2, 0},{ 0, 4, 0}}); // G(+4) = 1/8*x^4-3/4*x^2*y^2+1/8*y^4 ~ C44
       
       /* Higher angular momentums can be added here,
        * but MOLDEN format explicitly defines only functions up to G
        */

       /*
       // for debug:
       final boolean OrcaSingConvention = false;
       if (OrcaSingConvention) {
           Quick_YLM[ 3][ 0].ScaleCoefsBy(-1); // F(-3)
           Quick_YLM[ 3][ 6].ScaleCoefsBy(-1); // F(+3)
           Quick_YLM[ 4][ 0].ScaleCoefsBy(-1); // G(-4)
           Quick_YLM[ 4][ 1].ScaleCoefsBy(-1); // G(-3)
           Quick_YLM[ 4][ 7].ScaleCoefsBy(-1); // G(+3)
           Quick_YLM[ 4][ 8].ScaleCoefsBy(-1); // G(+4) 
       }
        */

       return Quick_YLM;
    }
    //--------------------------------------------------------------------------
    /**
     * @return   Returns a set of pure spherical harmonics each having the squared norm equal to 4*Pi
     * Version: 11.Jan.2014
     *
     */
    public static Polynom3D[][] Get_Quick_YLM_Norm4PI() {
        Polynom3D[][] result = Get_Quick_YLM();
        double[][] norms2 = Get_Quick_YLM_Norm2(result);
        for (int L=0; L<result.length; L++)
            for (int lm=0; lm<result[L].length; lm++) result[L][lm].ScaleCoefsBy( 1.0 / Math.sqrt(norms2[L][lm]) );
        return result;
    }
    //--------------------------------------------------------------------------
    /*
     * Below there is a method for converting cartesian spherical functions defined as
     * 0:                                            1 (S, 1 function)
     * 1,2,3:                                        x,y,z (P, 3 functions)
     * 4,5,6,7,8,9:                                  xx, yy, zz, xy, xz, yz (D, 6 functions)
     * 10,11,12,13,14,15,16,17,18,19:                xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz (F, 10 functions)
     * 20,21,22,23,24,25,26,27,28,29,30,31,32,33,34: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
     *                                               xxyy xxzz yyzz xxyz yyxz zzxy (G, 15 functions)
     * to their spherical counterparts
     *
     */
    //--------------------------------------------------------------------------
    static int[][] molden_Cart_powers = new int[][] {
         /* S: 1 */ { 0,0,0 },

         /* P: x */ { 1,0,0 },
         /* P: y */ { 0,1,0 },
         /* P: z */ { 0,0,1 },

         /* D: xx */ { 2,0,0 },
         /* D: yy */ { 0,2,0 },
         /* D: zz */ { 0,0,2 },
         /* D: xy */ { 1,1,0 },
         /* D: xz */ { 1,0,1 },
         /* D: yz */ { 0,1,1 },

         /* F: xxx */ { 3,0,0 },
         /* F: yyy */ { 0,3,0 },
         /* F: zzz */ { 0,0,3 },
         /* F: xyy */ { 1,2,0 },
         /* F: xxy */ { 2,1,0 },
         /* F: xxz */ { 2,0,1 },
         /* F: xzz */ { 1,0,2 },
         /* F: yzz */ { 0,1,2 },
         /* F: yyz */ { 0,2,1 },
         /* F: xyz */ { 1,1,1 },

         /* G: xxxx */ { 4,0,0 },
         /* G: yyyy */ { 0,4,0 },
         /* G: zzzz */ { 0,0,4 },
         /* G: xxxy */ { 3,1,0 },
         /* G: xxxz */ { 3,0,1 },
         /* G: yyyx */ { 1,3,0 },
         /* G: yyyz */ { 0,3,1 },
         /* G: zzzx */ { 1,0,3 },
         /* G: zzzy */ { 0,1,3 },
         /* G: xxyy */ { 2,2,0 },
         /* G: xxzz */ { 2,0,2 },
         /* G: yyzz */ { 0,2,2 },
         /* G: xxyz */ { 2,1,1 },
         /* G: yyxz */ { 1,2,1 },
         /* G: zzxy */ { 1,1,2 },
    };
    //--------------------------------------------------------------------------
    /**
     * @return the norms of x^molden_Cart_powers[k][0]·y^molden_Cart_powers[k][1]·z^molden_Cart_powers[k][2]/r^(...)
     * divided over 4*Pi
     * Created: 15.Jan.2014
     */
    public static double[] molden_cart_norms2_over_4Pi() {
        double[] result = new double[molden_Cart_powers.length];
        for (int i=0; i<molden_Cart_powers.length; ++i) {
            result[i] =  OverlapIntegrals.primitive_int_1D(2*molden_Cart_powers[i][0], 2*1.0);
            result[i] *= OverlapIntegrals.primitive_int_1D(2*molden_Cart_powers[i][1], 2*1.0);
            result[i] *= OverlapIntegrals.primitive_int_1D(2*molden_Cart_powers[i][2], 2*1.0);
            /*
             * Now  result[i] = INT(  {x^nx·y^ny·z^nz·exp(-1.0·r^2)}^2, d3r ) =
             *  = INT(  {x^nx·y^ny·z^nz/r^(nx+ny+nz)}^2, dOmega ) · INT( {r^(nx+ny+nz)·exp(-1.0·r^2)}^2 * r^2, r=0..+infinity )
             * Hence,
             * norm2 = result[i] / INT( r^(2+2*(nx+ny+nz))·exp(-2*1.0·r^2), r=0..+infinity )
             * 
             */
            result[i] /= OverlapIntegrals.primitive_int_1D_Sphr(
                    2 + 2*(molden_Cart_powers[i][0]+molden_Cart_powers[i][1]+molden_Cart_powers[i][2]),
                    2*1.0);
            result[i] /= 4*Math.PI;
        }
        return result;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a set of expansions of primive cartesian harmonics (as defined above) over
     * 'pure' spherical functions
     *
     * Do not forget to multiply terms with L lower than YLM_Series.primary_L by 
     * (x^2+y^2+z^2)^((primary_L - L[i])/2) if neccessary
     *
     * @return An array of YLM_Series objects each containing an expansion over 'pure' YLMs
     * for one cartesian function (in the order defined above)
     */
    public static YLM_Series[] Cartesian_to_Pure() {
        YLM_Series[] result = new YLM_Series[35];

         result[ 0] = new YLM_Series( 0, new int[][]{ { 0, 0}  }, new double[]{1.0 } );  // x^0·y^0·z^0

         result[ 1] = new YLM_Series( 1, new int[][]{ { 1, 0}  }, new double[]{1.0 } );  // x^1·y^0·z^0
         result[ 2] = new YLM_Series( 1, new int[][]{ { 1, 1}  }, new double[]{1.0 } );  // x^0·y^1·z^0
         result[ 3] = new YLM_Series( 1, new int[][]{ { 1,-1}  }, new double[]{1.0 } );  // x^0·y^0·z^1

         result[ 4] = new YLM_Series( 2, new int[][]{ { 2, 0},  { 2, 2},  { 0, 0}  }, new double[]{-1/Sqrt3, 1.0, 1.0/3.0 } );  // x^2·y^0·z^0
         result[ 5] = new YLM_Series( 2, new int[][]{ { 2, 0},  { 2, 2},  { 0, 0}  }, new double[]{-1/Sqrt3, -1.0, 1.0/3.0 } );  // x^0·y^2·z^0
         result[ 6] = new YLM_Series( 2, new int[][]{ { 2, 0},  { 0, 0}  }, new double[]{2/Sqrt3, 1.0/3.0 } );  // x^0·y^0·z^2
         result[ 7] = new YLM_Series( 2, new int[][]{ { 2,-2}  }, new double[]{1.0 } );  // x^1·y^1·z^0
         result[ 8] = new YLM_Series( 2, new int[][]{ { 2, 1}  }, new double[]{1.0 } );  // x^1·y^0·z^1
         result[ 9] = new YLM_Series( 2, new int[][]{ { 2,-1}  }, new double[]{1.0 } );  // x^0·y^1·z^1

         result[10] = new YLM_Series( 3, new int[][]{ { 3, 1},  { 3, 3},  { 1, 0}  }, new double[]{-3/Sqrt10, Sqrt3/Sqrt2, 0.6 } );  // x^3·y^0·z^0
         result[11] = new YLM_Series( 3, new int[][]{ { 3,-3},  { 3,-1},  { 1, 1}  }, new double[]{-Sqrt3/Sqrt2, -3/Sqrt10, 0.6 } );  // x^0·y^3·z^0
         result[12] = new YLM_Series( 3, new int[][]{ { 3, 0},  { 1,-1}  }, new double[]{6/Sqrt15, 0.6 } );  // x^0·y^0·z^3
         result[13] = new YLM_Series( 3, new int[][]{ { 3, 1},  { 3, 3},  { 1, 0}  }, new double[]{-1/Sqrt10, -Sqrt3/Sqrt2, 0.2 } );  // x^1·y^2·z^0
         result[14] = new YLM_Series( 3, new int[][]{ { 3,-3},  { 3,-1},  { 1, 1}  }, new double[]{Sqrt3/Sqrt2, -1/Sqrt10, 0.2 } );  // x^2·y^1·z^0
         result[15] = new YLM_Series( 3, new int[][]{ { 3, 0},  { 3, 2},  { 1,-1}  }, new double[]{-3/Sqrt15, 1.0, 0.2 } );  // x^2·y^0·z^1
         result[16] = new YLM_Series( 3, new int[][]{ { 3, 1},  { 1, 0}  }, new double[]{4/Sqrt10, 0.2 } );  // x^1·y^0·z^2
         result[17] = new YLM_Series( 3, new int[][]{ { 3,-1},  { 1, 1}  }, new double[]{4/Sqrt10, 0.2 } );  // x^0·y^1·z^2
         result[18] = new YLM_Series( 3, new int[][]{ { 3, 0},  { 3, 2},  { 1,-1}  }, new double[]{-3/Sqrt15, -1.0, 0.2 } );  // x^0·y^2·z^1
         result[19] = new YLM_Series( 3, new int[][]{ { 3,-2}  }, new double[]{1.0 } );  // x^1·y^1·z^1
         
         result[20] = new YLM_Series( 4, new int[][]{ { 4, 0},  { 4, 2},  { 4, 4},  { 2, 0},  { 2, 2},  { 0, 0}  }, new double[]{3/Sqrt35, -2/Sqrt7, 1.0, -2*Sqrt3/7, 6.0/7.0, 0.2 } );  // x^4·y^0·z^0
         result[21] = new YLM_Series( 4, new int[][]{ { 4, 0},  { 4, 2},  { 4, 4},  { 2, 0},  { 2, 2},  { 0, 0}  }, new double[]{3/Sqrt35, 2/Sqrt7, 1.0, -2*Sqrt3/7, -6.0/7.0, 0.2 } );  // x^0·y^4·z^0
         result[22] = new YLM_Series( 4, new int[][]{ { 4, 0},  { 2, 0},  { 0, 0}  }, new double[]{8/Sqrt35, 4*Sqrt3/7,  0.2 } );  // x^0·y^0·z^4
         result[23] = new YLM_Series( 4, new int[][]{ { 4,-4},  { 4,-2},  { 2,-2}  }, new double[]{1.0, -1/Sqrt7, 3.0/7.0 } );  // x^3·y^1·z^0
         result[24] = new YLM_Series( 4, new int[][]{ { 4, 1},  { 4, 3},  { 2, 1}  }, new double[]{-3/Sqrt14, 1/Sqrt2, 3.0/7.0 } );  // x^3·y^0·z^1
         result[25] = new YLM_Series( 4, new int[][]{ { 4,-4},  { 4,-2},  { 2,-2}  }, new double[]{-1.0, -1/Sqrt7, 3.0/7.0 } );  // x^1·y^3·z^0
         result[26] = new YLM_Series( 4, new int[][]{ { 4,-3},  { 4,-1},  { 2,-1}  }, new double[]{-1/Sqrt2, -3/Sqrt14, 3.0/7.0 } );  // x^0·y^3·z^1
         result[27] = new YLM_Series( 4, new int[][]{ { 4, 1},  { 2, 1}  }, new double[]{4/Sqrt14, 3.0/7.0 } );  // x^1·y^0·z^3
         result[28] = new YLM_Series( 4, new int[][]{ { 4,-1},  { 2,-1}  }, new double[]{4/Sqrt14, 3.0/7.0 } );  // x^0·y^1·z^3
         result[29] = new YLM_Series( 4, new int[][]{ { 4, 0},  { 4, 4},  { 2, 0},  { 0, 0}  }, new double[]{1/Sqrt35, -1.0, -2*Sqrt3/21, 1.0/15.0 } );  // x^2·y^2·z^0
         result[30] = new YLM_Series( 4, new int[][]{ { 4, 0},  { 4, 2},  { 2, 0},  { 2, 2},  { 0, 0}  }, new double[]{-4/Sqrt35, 2/Sqrt7, Sqrt3/21, 1.0/7.0, 1.0/15.0 } );  // x^2·y^0·z^2
         result[31] = new YLM_Series( 4, new int[][]{ { 4, 0},  { 4, 2},  { 2, 0},  { 2, 2},  { 0, 0}  }, new double[]{-4/Sqrt35, -2/Sqrt7, Sqrt3/21, -1.0/7.0, 1.0/15.0 } );  // x^0·y^2·z^2
         result[32] = new YLM_Series( 4, new int[][]{ { 4,-3},  { 4,-1},  { 2,-1}  }, new double[]{1/Sqrt2, -1/Sqrt14, 1.0/7.0 } );  // x^2·y^1·z^1
         result[33] = new YLM_Series( 4, new int[][]{ { 4, 1},  { 4, 3},  { 2, 1}  }, new double[]{-1/Sqrt14, -1/Sqrt2, 1.0/7.0 } );  // x^1·y^2·z^1
         result[34] = new YLM_Series( 4, new int[][]{ { 4,-2},  { 2,-2}  }, new double[]{2/Sqrt7, 1.0/7.0 } );  // x^1·y^1·z^2
        
         return result;
    }
    //--------------------------------------------------------------------------

    /** Calculates the square of the actual norm (divided over 4*Pi) for the
     * spherical functions passed as a parameter(s)
     *
     * @param Quick_YLM An array of spherical functions
     * @return the square of the actual norm of Quick_YLM divided over 4*Pi
     * 
     */
    public static double[][] Get_Quick_YLM_Norm2(Polynom3D[][] Quick_YLM) {
        double[][] result;
        // allocate memory
        double RadialNorm2;
        result = new double[Quick_YLM.length][];
        // allocate memory & calculate norm2: NORM2 = INT(Y_LM(theta,phi))^2 · dOmega, over a sphere)
        for (int L=0; L<Quick_YLM.length; L++) {
            result[L] = new double[Quick_YLM[L].length];
            for (int m=0; m<Quick_YLM[L].length; m++) {
                // calculate norm of YLM
                BasisFunction BS = new BasisFunction(L, 0, new double[]{0,0,0}, 1);
                BS.coefs[0] = 1.0;
                BS.exponents[0] = 1.0; // some fake exponent
                BS.Quick_YLM = Quick_YLM;
                RadialNorm2 = OverlapIntegrals.primitive_int_1D_Sphr(2+2*BS.L, 2*1.0)*4*Math.PI;
                try {
                    result[L][m] = BS.OverlapWith(BS) / RadialNorm2; // OverlapWith has a protection against BS.Quick_YLM = null
                } catch (Exception E) {}
            }
        }        
        return result;
    }
    //--------------------------------------------------------------------------

    
    // Squared roots of associated Legendre polinomials, different from 0, +1 and -1
    // With their values available, it is quite easy to generate 'the most interesting' 
    // parts of associated Legendre polinomials as PRODUCT( ((z/r)^2 - el_k), k = 1...(L-|m|)/2 )
    // ( which is just a rewrite of PRODUCT( x-x_i, i=1..n) )
    // 
    // Each ALegRootsSquared[L][|m|] has (L-|m|)/2 elements
    public static final double[][][] ALegRootsSquared = new double[][][] {
        // L = 0: -- has no roots
        new double[][] {              
          new double[]{ }               // m = 0
        },                            

        // L = 1: -- has no roots     
        new double[][] {              
          new double[]{ },              // m = 0
          new double[]{ }               // m = 1
        },                            

        // L = 2:                     
        new double[][] {              
          new double[]{ 1.0/3.0 },      // m = 0
          new double[]{ },              // m = 1
          new double[]{ }               // m = 2
        },                            

        // L = 3:                     
        new double[][] {              
          new double[]{ 3.0/5.0 },      // m = 0
          new double[]{ 1.0/5.0 },      // m = 1
          new double[]{ },              // m = 2
          new double[]{ }               // m = 3
        },

        // L = 4:
        new double[][] {
          new double[]{ 0.11558710999704793517,  0.74155574714580920772 },  // m = 0, analytical form:  3.0/7.0 -/+ 2.0/7.0 * SQRT(6.0/5.0)
          new double[]{ 3.0/7.0 },      // m = 1
          new double[]{ 1.0/7.0 },      // m = 2
          new double[]{ },              // m = 3
          new double[]{ }               // m = 4
        },

        // L = 5: 
        new double[][] {
          new double[]{ 0.28994919792569030224,  0.82116191318542080892 },  // m = 0  5.0/9.0 -/+ 2.0/9.0 * SQRT(10.0/7.0)
          new double[]{ 0.08135701799384851519,  0.58530964867281815149 },  // m = 1  1.0/3.0 -/+ 2.0/3.0 * SQRT( 1.0/7.0)
          new double[]{ 1.0/3.0 },                                    // m = 2
          new double[]{ 1.0/9.0 },                                    // m = 3
          new double[]{ },                                            // m = 4
          new double[]{ }                                             // m = 5
        }/*,

        // L = 6: TODO: refine numbers
        new double[][] {
          new double[]{ 0.05693911596700739,  0.43719785275109463,  0.86949939491826234 },  // m = 0
          new double[]{ 0.21981919113894424,  0.68927171795196474 },  // m = 1
          new double[]{ 0.06278172029468174,  0.48267282515986437 },  // m = 2
          new double[]{ 0.27272727272727215 },   // 3.0/11.0  // m = 3
          new double[]{ 0.09090909090909115 },   // 1.0/11.0  // m = 4    
          new double[]{ },                                    // m = 5
          new double[]{ }                                     // m = 6
        },

        // L = 7: TODO: refine numbers
        new double[][] {
          new double[]{ 0.16471028689654230,  0.54986849921644299,  0.90080582927163033 },  // m = 0
          new double[]{ 0.04380616261458938,  0.35010910470801310,  0.75993088652355167 },  // m = 1
          new double[]{ 0.17716927308823283,  0.59206149614253545 },  // m = 2
          new double[]{ 0.05111562827041795,  0.41042283326804369 },  // m = 3
          new double[]{ 0.23076923076923037 },                // m = 4
          new double[]{ 0.07692307692307718 },                // m = 5
          new double[]{ },                                    // m = 6
          new double[]{ }                                     // m = 7
        }*/
    };

    
    
    /** Numerical values of a common multiplier of [L][L+m] associated Legendre polinomial
     * with (-1)^m included (in fact, this 'phase factor' is just not used by MOLDEN)
     * 
     **/ 
    public static double[][] ALegNormalizer = __ALegNormalizer_initializer(); // call initialization method defined below
    private static double[][] __ALegNormalizer_initializer() {
        double[][] ALegNormalizer ; // can be commented out
        /*
        The 'true' Y_{L,m} with  __ m >= 0 __ looks like this:

                     ( 2*L+1    (L-m)! )      1      (   rho^2  )^( |m|/2 )    d^(L+|m|)                    / cos(|m|*phi)*SQRT(2), m > 0
                 sqrt( ------ * ------ ) * ------  * ( ---------)           * ----------- ( (x^2-1)^L ) ) *(
                     (  4*Pi    (L+m)! )   2^L*L!    ( z^2+rho^2)             d x^(L+|m|)                   \ 1.0, if m == 0
                                         |                                                                |
                                         +----------- the true associated Legendre polynomial ------------+
                                          it has coef. at x^(L-|m|) equal to (2L)!/(L-|m|)! / (2^L*L!)
        Note that '*SQRT(2)' comes from the fact that for m>0 the real spherical harmonics are ~Re(Y_{L,m})~cos(m*phi),
        but in contrast to |Y_{L,m}|^2 ~ |e^(I*m*phi)|^2 = 1 which integrates to 2*Pi, in real case we have
        |Re(Y_{L,m})| ~ |cos(m*phi)|^2, which integrates to 2*Pi/2=Pi => need to mul. by SQRT(2) to preserve
        normalization (or, in other words, to modify the prefactor properly).
    
        Together with the first sqrt this leads to coef. equal to 
          (2L)! / (2^L*L!) / SQRT( (L-|m|)! * (L+|m|)! ) * SQRT((2L+1)/4/Pi)
        In contrast to that, 'our' associated Legendre polynomial has coef. at x^(L-|m|) equal to 1.0
    
    
        The 'true' Y_{L,m} with  __ m < 0 __ looks like this:
    
                     ( 2*L+1    (L- m )! )  (L-|m|)!       1      (   rho^2  )^( |m|/2 )    d^(L+|m|)                   
                 sqrt( ------ * -------- )  --------- * ------  * ( ---------)           * ----------- ( (x^2-1)^L ) ) * sin(|m|*phi)*SQRT(2)
                     (  4*Pi    (L+ m )! )  (L+|m|)!    2^L*L!    ( z^2+rho^2)             d x^(L+|m|)                  
                                         |                                                                |
                                         +----------- the true associated Legendre polynomial ------------+
                                          it has coef. at x^(L-|m|) equal to (2L)!/(L+|m|)! / (2^L*L!)    
        
        Notice that since m is negative, (L-m)! == (L+|m|)! and this _!!AGAIN!!_ leads to coef. equal to 
          (2L)! / (2^L*L!) / SQRT( (L-|m|)! * (L+|m|)! ) * SQRT((2L+1)/4/Pi)
        In contrast to that, 'our' associated Legendre polynomial has coef. at x^(L-|m|) equal to 1.0
    
        Note that the spherical functions used by MOLDEN do NOT use any kind of a prefactor (-1)^m.
        Instead, Y_{L, -|m|} ~ sin(|m|*phi) and Y_{L, +|m|} ~ cos(|m|*phi).
        
        */
        int Lsz = ALegRootsSquared.length;
        ALegNormalizer = new double[Lsz][];
        for (int L=0; L<Lsz; L++) {
            ALegNormalizer[L] = new double[2*L+1];            
            double tmp = 1.0; // compute (2*L)! / L! / L! / 2^L
            for (int i=1; i<=L; i++)
                tmp *= 2.0 - 1.0/i ; // i.e.,  (2.0*i)/i * (2.0*i-1.0)/i  / 2.0            
            ALegNormalizer[L][L] = tmp * Math.sqrt((2.0*L + 1.0)/4/Math.PI);
            for (int m=1; m<=L; m++) {
                ALegNormalizer[L][L + m] = ALegNormalizer[L][L + m-1] / Math.sqrt( (L+m)/(L-m+1.0) ); // /SQRT( (L+m)/(L-m+1) )
                if (m==1) {
                    // introduce SQRT(2), present at sin(|m|*phi) and cos(m*phi) for all m > 0
                    ALegNormalizer[L][L + m] *= Math.sqrt(2.0);
                }
                ALegNormalizer[L][L - m] = ALegNormalizer[L][L + m]; // attn.: NO (-1)^(...) is necessary !
            }
        }
        /*System.out.println(ALegNormalizer[3][3]);
        System.out.println(ALegNormalizer[3][3+1]);
        System.out.println(ALegNormalizer[3][3+2]);*/
        return ALegNormalizer;
    }
    
}
