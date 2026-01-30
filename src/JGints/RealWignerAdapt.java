package JGints;

/**
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
 * @author (c) Tymofii Nikolaienko, 2014
 */
public class RealWignerAdapt {
    /**
     * In the following we assume that the axis of the new coordinate system
     * are obtained by rotation around
     *   n || v = [ ez_new x ez_old ] = { ez_new_Y, -ez_new_X, 0 },
     * where ez_new = (RB-RA)/|RB-RA|, and ez_old = {0,0,1} (expressed in old axis).
     * Clearly, |v|==sin(alpha), where alpha is the angle by which _all_ the old
     * axes are to be rotated (around v/|v|) in order to get the new axis.
     * This procedure defines the coordinates of all new axis w.r.t. the old ones,
     * i.e., the whole transformation matrix.
     * Note that in the old axes the components of ez_new are defined as
     *  z_new = { sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) }.
     * Thus, the whole transfromation matrix is implicitly parameterized
     * by the values of theta and phi:
     * 
     */
    /**
     * @returns coordinates of the new axis w.r.t. old ones
     * @param Theta
     * @param Phi
     */
    public static double[][] newAxisCoords(double Theta, double Phi) {
        double sinT = Math.sin(Theta);
        double cosI = Math.cos(Phi);
        double sinI = Math.sin(Phi);
        double cosT = Math.cos(Theta);

        /* coordinates of ez_new w.r.t. old axis:
         Xp = sinT*cosI;
         Yp = sinT*sinI;
         Zp = cosT;
        
        Then the rotation vector N is:
         N = [ {Xp,Yp,Zp} x {0,0,1} ]  =>
         Nx = (Yp*1-Zp*0) = sinT*sinI
         Ny = (Zp*0-Xp*1) = -sinT*cosI
         Nz = (Xp*0-Yp*0) = 0
        */
        /* 
        // not used, in actual computations:        
        double  sinJ, nx, ny, nz;
        // normalize the rotation vector and compute sin and cos of rotation angle:
        // sinJ = sqrt(Nx^2+Ny^2+Nz^2) = |sinT|
        if (sinT > 0) {
                sinJ = sinT;   // sqrt(Nx^2+Ny^2+Nz^2);
                nx =   sinI;   // Nx/sinJ = sinT*sinI / |sinT|
                ny =  -cosI;   // Ny/sinJ = -sinT*cosI / |sinT|
                nz =   0;      // Nz/sinJ;
        } else {
            // never happens since 0<=theta<=Pi  => sinT >=0             
                sinJ = -sinT;  // sqrt(Nx^2+Ny^2+Nz^2);
                nx =   -sinI;  // Nx/sinJ = sinT*sinI / |sinT|
                ny =    cosI;  // Ny/sinJ = -sinT*cosI / |sinT|
                nz =    0;     // Nz/sinJ;
        }
        
        double cosJ = cosT; // < {Xp,Yp,Zp} . {0,0,1}> = Zp
        */
        
        /* Now the rotationa matrix is built as
           a_new = n * (n*a) + cosJ * (a - n*(n*a)) + sinJ * |a - n*(n*a)| * [(a - n*(n*a)) x n]/|[n x (a - n*(n*a))]| =
                 = n * (n*a) + cosJ * (a - n*(n*a)) + sinJ * |a - n*(n*a)| * [a x n]/|[n x a]|
         Note that |a - n*(n*a)|^2 = <(a - n*(n*a)).(a - n*(n*a))> = a^2 - (n*a)^2, 
          and that |[n x a]| = |n|^2*|a|^2 * (1 - cos(n^a)^2) = |n|^2*|a|^2 - (n*a)^2, where |n| = 1,
         so that |a - n*(n*a)| == |[n x a]|, and
           a_new = n * (n*a) + cosJ * (a - n*(n*a)) + sinJ * [a x n ] = 
                 = a + n*(n*a)*(1-cosJ) + sinJ * [a x n ].
                                          ~~~~~~~~~~~~~~~ == -sinJ * [n x a]
         In components:
            a_new_X = X*cosJ + (1-cosJ)*nx*an - sinJ*(ny*Z-nz*Y),
            a_new_Y = Y*cosJ + (1-cosJ)*ny*an - sinJ*(nz*X-nx*Z),
            a_new_Z = Z*cosJ + (1-cosJ)*nz*an - sinJ*(nx*Y-ny*X),
         wnere a = {X,Y,Z}, an = <a.n>.
        
         In particular, 
        a_old = ex_old = {1,0,0) => an = nx:
            a_new_X = cosJ + (1-cosJ)*nx*nx       = cosT*(nz^2+ny^2) + 1*nx*nx = cosT*cosI*cosI + sinI^2 = 1 - cosI^2*(1-Zp)
            a_new_Y = (1-cosJ)*ny*nx - sinJ*nz    = -(1-cosT)*sinI*cosI = -(1-Zp)*sinI*cosI
            a_new_Z = (1-cosJ)*nz*nx - sinJ*(-ny) = - sinT*cosI   = -Xp
        a_old = ey_old = {0,1,0) => an = ny:
            a_new_X = (1-cosJ)*nx*ny - sinJ*(-nz) = -(1-cosJ)*sinI*cosI
            a_new_Y = cosJ + (1-cosJ)*ny*ny       = 1*ny*ny + cosJ*(nx^2+nz^2) = 1*ny*ny + cosT*sinI*sinI = 1-sinI^2 + Zp*sinI^2 = 1 - sinI^2*(1-Zp)
            a_new_Z = (1-cosJ)*nz*ny - sinJ*nx    = - sinT*sinI   = -Yp
        a_old = ez_old = {0,0,1) => an = nz = 0:
            a_new_X = -sinJ*(ny)                  = sinT*cosI = Xp -- ok
            a_new_Y = -sinJ*(-nx)                 = sinT*sinI = Yp -- ok
            a_new_Z = cosJ                        = cosT      = Zp -- ok
        */
        double Xp = sinT*cosI;
        double Yp = sinT*sinI;
        double Zp = cosT;
    
        double csI = sinI*cosI; 
        double cI2 = cosI*cosI;
        double sI2 = sinI*sinI;
        
        return new double[][]{
            /*               ( ex_new )      (ey_new)     (ez_new)        */
            new double[]{ 1 - cI2*(1-Zp),     -csI*(1-Zp),   Xp  }, // x-component
            new double[]{    -csI*(1-Zp),  1 - sI2*(1-Zp),   Yp  }, // y-component
            new double[]{            -Xp,             -Yp,   Zp  }, // z-component
        };
    }
    //--------------------------------------------------------------------------
    //
    /**
     * evaluates c[0][n] = (cos(i))^(n+1), n=0,1,2,...
     *           c[1][n] = (cos(t))^(n+1), n=0,1,2,...
     *           s[0][n] = (sin(i))^(n+1), n=0,1,2,...
     *           s[1][n] = (sin(t))^(n+1), n=0,1,2,...
     * and @returns array(c,s)
     */
    public static double[][][] create_sincos(double t, double i) {
        int MAX_POW = 8; // up to 2*L_MAX
        double[][] c = new double[2][MAX_POW];
        double[][] s = new double[2][MAX_POW];
        c[0][0] = Math.cos(i); c[1][0] = Math.cos(t);
        s[0][0] = Math.sin(i); s[1][0] = Math.sin(t);
        for (int mu=0; mu<2; ++mu)
          for (int p=1; p<MAX_POW; ++p) {
              c[mu][p] = c[mu][p-1] * c[mu][0];
              s[mu][p] = s[mu][p-1] * s[mu][0];
          }

        return new double[][][]{c,s};
    }
    //--------------------------------------------------------------------------

    public final static double sqrt2 = Math.sqrt(2);
    public final static double sqrt3 = Math.sqrt(3);
    public final static double sqrt5 = Math.sqrt(5);
    public final static double sqrt7 = Math.sqrt(7);
    public final static double sqrt6 = Math.sqrt(6);
    public final static double sqrt10 = Math.sqrt(10);
    public final static double sqrt14 = Math.sqrt(14);
    public final static double sqrt15 = Math.sqrt(15);
    public final static double sqrt35 = Math.sqrt(35);
    public final static double sqrt70 = Math.sqrt(70);

    //
    //--------------------------------------------------------------------------
    // Real spherical harmonics transformation matrix for any L <= L_MAX
    public static double[][] Wigner(double[][][] sincos, int L){
        switch (L) {
            case 0: return Wigner_S(sincos);
            case 1: return Wigner_P(sincos);
            case 2: return Wigner_D(sincos);
            case 3: return Wigner_F(sincos);
            case 4: return Wigner_G(sincos);
        }
        return null;
    }

    //--------------------------------------------------------------------------
    // Real spherical harmonics transformation matrix for L = 0 (S)
    public static double[][] Wigner_S(double[][][] sincos){
        double[][] result = new double[1][1];
        result[0][0] = 1;
        return result;
    }

    //--------------------------------------------------------------------------
    // Real spherical harmonics transformation matrix for L = 1 (P)
    public static double[][] Wigner_P(double[][][] sincos){
        double[][] c = sincos[0];
        double[][] s = sincos[1];
        double[][] result = new double[3][3];

        // ORDERING: y,z,x  v.s. MOLDEN ordering: z,x,y
        result[0][0] = c[1][0]+(1-c[1][0])*c[0][1];
        result[0][1] = -s[1][0]*s[0][0];
        result[0][2] = c[0][0]*s[0][0]*(-1+c[1][0]);
        result[1][0] = -1.0 * result[0][1];
        result[1][1] = c[1][0];
        result[1][2] = s[1][0]*c[0][0];
        result[2][0] = 1.0 * result[0][2];
        result[2][1] = -1.0 * result[1][2];
        result[2][2] = 1+(-1+c[1][0])*c[0][1];

        return result;
    }
    //--------------------------------------------------------------------------
    // Real spherical harmonics transformation matrix for L = 2 (D)
    public static double[][] Wigner_D(double[][][] sincos){
        double[][] c = sincos[0];
        double[][] s = sincos[1];
        double[][] result = new double[5][5];

        //double t = Math.atan2(s[1][0], c[1][0]);
        //double i = Math.atan2(s[0][0], c[0][0]);

        result[0][0] = c[1][0]-4*c[1][0]*c[0][1]+4*c[1][0]*c[0][3]+2*c[0][1]-2*c[0][3]+2*c[1][1]*c[0][1]-2*c[0][3]*c[1][1];
        result[0][1] = s[1][0]*c[0][0]*(-2*c[1][0]+2*c[1][0]*c[0][1]+1-2*c[0][1]);
        result[0][2] = -(-1+c[1][1])*c[0][0]*s[0][0]*sqrt3;
        result[0][3] = -s[1][0]*s[0][0]*(1-2*c[0][1]+2*c[1][0]*c[0][1]);
        result[0][4] = s[0][0]*c[0][0]*(2*c[1][0]-4*c[1][0]*c[0][1]-1+2*c[0][1]+2*c[1][1]*c[0][1]-c[1][1]);
        result[1][0] = -1.0 * result[0][1];
        result[1][1] = -1+c[0][1]+2*c[1][1]-2*c[1][1]*c[0][1]+c[1][0]*c[0][1];
        result[1][2] = -c[1][0]*s[0][0]*sqrt3*s[1][0];
        result[1][3] = c[0][0]*s[0][0]*(-c[1][0]-1+2*c[1][1]);
        result[1][4] = s[1][0]*s[0][0]*(-c[1][0]+2*c[1][0]*c[0][1]-2*c[0][1]);
        result[2][0] = 1.0 * result[0][2];
        result[2][1] = -1.0 * result[1][2];
        result[2][2] = -1.0/2+3.0/2*c[1][1];
        result[2][3] = c[1][0]*s[1][0]*sqrt3*c[0][0];
        result[2][4] = -1.0/2*sqrt3+sqrt3*c[0][1]-sqrt3*c[1][1]*c[0][1]+1.0/2*sqrt3*c[1][1];
        result[3][0] = -1.0 * result[0][3];
        result[3][1] = 1.0 * result[1][3];
        result[3][2] = -1.0 * result[2][3];
        result[3][3] = 2*c[1][1]*c[0][1]+c[1][0]-c[1][0]*c[0][1]-c[0][1];
        result[3][4] = s[1][0]*c[0][0]*(2-2*c[0][1]+2*c[1][0]*c[0][1]-c[1][0]);
        result[4][0] = 1.0 * result[0][4];
        result[4][1] = -1.0 * result[1][4];
        result[4][2] = 1.0 * result[2][4];
        result[4][3] = -1.0 * result[3][4];
        result[4][4] = 2*c[0][3]*c[1][1]+1.0/2*c[1][1]-2*c[1][1]*c[0][1]+2*c[0][3]+1.0/2-2*c[0][1]+4*c[1][0]*c[0][1]-4*c[1][0]*c[0][3];

        return result;
    }
    //--------------------------------------------------------------------------
    // Real spherical harmonics transformation matrix for L = 3 (F)
    public static double[][] Wigner_F(double[][][] sincos){
        double[][] c = sincos[0];
        double[][] s = sincos[1];
        double[][] result = new double[7][7];

        result[0][0] = 4*c[0][5]+3.0/4*c[1][0]+9.0/4*c[0][1]-6*c[0][3]+1.0/4*c[1][2]+6*c[0][3]*c[1][2]-18*c[0][3]*c[1][1]+18*c[0][3]*c[1][0]+12*c[0][5]*c[1][1]-9.0/4*c[1][2]*c[0][1]-4*c[1][2]*c[0][5]-12*c[1][0]*c[0][5]-27.0/4*c[1][0]*c[0][1]+27.0/4*c[1][1]*c[0][1];
        result[0][1] = 1.0/2*s[1][0]*sqrt2*sqrt3*c[0][0]*(-5*c[0][1]+4*c[0][3]-8*c[0][3]*c[1][0]+1+c[1][1]-5*c[1][1]*c[0][1]+4*c[0][3]*c[1][1]-3*c[1][0]+10*c[1][0]*c[0][1]);
        result[0][2] = -1.0/4*sqrt5*sqrt3*c[1][0]-sqrt5*sqrt3*c[0][3]*c[1][1]+1.0/4*sqrt5*sqrt3*c[1][2]-3.0/4*sqrt5*sqrt3*c[0][1]-sqrt5*sqrt3*c[0][3]*c[1][0]+sqrt5*sqrt3*c[0][3]+sqrt5*sqrt3*c[0][3]*c[1][2]+3.0/4*sqrt5*sqrt3*c[0][1]*c[1][1]+5.0/4*sqrt5*sqrt3*c[0][1]*c[1][0]-5.0/4*sqrt5*sqrt3*c[0][1]*c[1][2];
        result[0][3] = 1.0/4*s[1][0]*sqrt10*s[0][0]*(1-4*c[0][1]-c[1][1]+4*c[1][1]*c[0][1]);
        result[0][4] = -1.0/4*(-3+4*c[0][1]+3*c[1][1]+c[1][0]-c[1][2]-4*c[1][1]*c[0][1]-4*c[1][0]*c[0][1]+4*c[1][2]*c[0][1])*s[0][0]*c[0][0]*sqrt5*sqrt3;
        result[0][5] = -1.0/4*sqrt3*sqrt2*s[1][0]*s[0][0]*(-16*c[0][3]*c[1][0]+12*c[1][0]*c[0][1]+1-6*c[0][1]+8*c[0][3]-6*c[1][1]*c[0][1]+8*c[0][3]*c[1][1]+c[1][1]);
        result[0][6] = 1.0/4*c[0][0]*s[0][0]*(16*c[0][1]-16*c[0][3]-3+3*c[1][2]-16*c[1][2]*c[0][1]+16*c[0][3]*c[1][2]-48*c[1][0]*c[0][1]+48*c[0][3]*c[1][0]-48*c[0][3]*c[1][1]-9*c[1][1]+48*c[1][1]*c[0][1]+9*c[1][0]);
        result[1][0] = -1.0 * result[0][1];
        result[1][1] = -1+2*c[0][3]*c[1][0]-6*c[0][3]*c[1][2]+8*c[0][3]*c[1][1]-2*c[1][0]*c[0][1]+4*c[0][1]-4*c[0][3]+2*c[1][1]-8*c[1][1]*c[0][1]+6*c[1][2]*c[0][1];
        result[1][2] = 1.0/2*(-c[0][1]+3*c[1][1]*c[0][1]-2*c[1][0]*c[0][1]-3*c[1][1]+1+c[1][0])*c[0][0]*s[1][0]*sqrt5*sqrt2;
        result[1][3] = -sqrt15*c[0][0]*s[0][0]*c[1][0]*(-1+c[1][1]);
        result[1][4] = -1.0/2*(-c[0][1]+3*c[1][1]*c[0][1]-2*c[1][0]*c[0][1]+c[1][0])*s[0][0]*s[1][0]*sqrt5*sqrt2;
        result[1][5] = c[0][0]*s[0][0]*(-2+c[1][0]+4*c[0][1]-3*c[1][2]+4*c[1][1]+6*c[1][2]*c[0][1]-8*c[1][1]*c[0][1]-2*c[1][0]*c[0][1]);
        result[1][6] = 1.0/2*s[1][0]*sqrt3*sqrt2*s[0][0]*(6*c[1][0]*c[0][1]-8*c[0][3]*c[1][0]-3*c[1][1]*c[0][1]+4*c[0][3]*c[1][1]-c[1][0]+4*c[0][3]-3*c[0][1]);
        result[2][0] = 1.0 * result[0][2];
        result[2][1] = -1.0 * result[1][2];
        result[2][2] = -11.0/4*c[1][0]+11.0/4*c[1][0]*c[0][1]+15.0/4*c[1][2]-15.0/4*c[1][2]*c[0][1]-1.0/4*c[0][1]+5.0/4*c[1][1]*c[0][1];
        result[2][3] = -1.0/4*(-1+5*c[1][1])*s[1][0]*sqrt6*s[0][0];
        result[2][4] = 1.0/4*(-11*c[1][0]-5*c[1][1]+1+15*c[1][2])*c[0][0]*s[0][0];
        result[2][5] = 1.0/4*s[1][0]*sqrt2*sqrt5*s[0][0]*(1-2*c[0][1]-4*c[1][0]*c[0][1]+6*c[1][1]*c[0][1]-3*c[1][1]);
        result[2][6] = -1.0/4*sqrt3*sqrt5*s[0][0]*c[0][0]*(-1-4*c[1][0]*c[0][1]-4*c[1][1]*c[0][1]+4*c[1][2]*c[0][1]-3*c[1][2]+c[1][1]+3*c[1][0]+4*c[0][1]);
        result[3][0] = -1.0 * result[0][3];
        result[3][1] = 1.0 * result[1][3];
        result[3][2] = -1.0 * result[2][3];
        result[3][3] = -3.0/2*c[1][0]+5.0/2*c[1][2];
        result[3][4] = 1.0/4*(-1+5*c[1][1])*c[0][0]*s[1][0]*sqrt6;
        result[3][5] = -1.0/2*sqrt15*c[1][0]+sqrt15*c[1][0]*c[0][1]-sqrt15*c[1][2]*c[0][1]+1.0/2*sqrt15*c[1][2];
        result[3][6] = -1.0/4*s[1][0]*sqrt10*c[0][0]*(3-4*c[0][1]-3*c[1][1]+4*c[1][1]*c[0][1]);
        result[4][0] = 1.0 * result[0][4];
        result[4][1] = -1.0 * result[1][4];
        result[4][2] = 1.0 * result[2][4];
        result[4][3] = -1.0 * result[3][4];
        result[4][4] = 15.0/4*c[1][2]*c[0][1]-11.0/4*c[1][0]*c[0][1]+1.0/4*c[0][1]-5.0/4*c[1][1]*c[0][1]-1.0/4+5.0/4*c[1][1];
        result[4][5] = 1.0/4*s[1][0]*sqrt2*sqrt5*c[0][0]*(1-2*c[0][1]-3*c[1][1]+6*c[1][1]*c[0][1]+4*c[1][0]-4*c[1][0]*c[0][1]);
        result[4][6] = -1.0/4*sqrt3*sqrt5-sqrt3*sqrt5*c[0][3]*c[1][2]+sqrt3*sqrt5*c[0][3]*c[1][0]+1.0/4*sqrt3*sqrt5*c[1][1]-sqrt3*sqrt5*c[0][3]-3.0/4*sqrt3*sqrt5*c[1][0]*c[0][1]-5.0/4*sqrt3*sqrt5*c[1][1]*c[0][1]+5.0/4*sqrt3*sqrt5*c[0][1]+3.0/4*sqrt3*sqrt5*c[1][2]*c[0][1]+sqrt3*sqrt5*c[0][3]*c[1][1];
        result[5][0] = -1.0 * result[0][5];
        result[5][1] = 1.0 * result[1][5];
        result[5][2] = -1.0 * result[2][5];
        result[5][3] = 1.0 * result[3][5];
        result[5][4] = -1.0 * result[4][5];
        result[5][5] = -1.0/2*c[1][0]+8*c[1][1]*c[0][1]-6*c[1][2]*c[0][1]+2*c[1][0]*c[0][1]-4*c[0][1]+4*c[0][3]+3.0/2*c[1][2]+6*c[0][3]*c[1][2]-2*c[0][3]*c[1][0]-8*c[0][3]*c[1][1];
        result[5][6] = 1.0/4*s[1][0]*sqrt3*sqrt2*c[0][0]*(8*c[0][3]*c[1][1]+3-10*c[0][1]+8*c[0][3]+20*c[1][0]*c[0][1]-16*c[0][3]*c[1][0]-4*c[1][0]+3*c[1][1]-10*c[1][1]*c[0][1]);
        result[6][0] = 1.0 * result[0][6];
        result[6][1] = -1.0 * result[1][6];
        result[6][2] = 1.0 * result[2][6];
        result[6][3] = -1.0 * result[3][6];
        result[6][4] = 1.0 * result[4][6];
        result[6][5] = -1.0 * result[5][6];
        result[6][6] = 4*c[1][2]*c[0][5]+12*c[1][0]*c[0][5]+1.0/4+27.0/4*c[1][0]*c[0][1]+9.0/4*c[1][2]*c[0][1]+6*c[0][3]+18*c[0][3]*c[1][1]-18*c[0][3]*c[1][0]-9.0/4*c[0][1]-4*c[0][5]+3.0/4*c[1][1]-6*c[0][3]*c[1][2]-12*c[0][5]*c[1][1]-27.0/4*c[1][1]*c[0][1];

        return result;
    }
    //--------------------------------------------------------------------------
    // Real spherical harmonics transformation matrix for L = 4 (G)
    public static double[][] Wigner_G(double[][][] sincos){
        double[][] c = sincos[0];
        double[][] s = sincos[1];
        double[][] result = new double[9][9];

        result[0][0] = 1.0/2*c[1][0]+16*c[1][3]*c[0][5]+32*c[1][2]*c[0][7]+32*c[0][7]*c[1][0]-10*c[0][3]-8*c[0][7]+16*c[0][5]+2*c[0][1]-64*c[0][5]*c[1][2]-64*c[0][5]*c[1][0]+40*c[0][3]*c[1][0]+40*c[1][2]*c[0][3]+1.0/2*c[1][2]-48*c[0][7]*c[1][1]-8*c[1][3]*c[0][7]-60*c[1][1]*c[0][3]-8*c[1][0]*c[0][1]-8*c[1][2]*c[0][1]+2*c[1][3]*c[0][1]-10*c[1][3]*c[0][3]+96*c[0][5]*c[1][1]+12*c[1][1]*c[0][1];
        result[0][1] = 1.0/4*s[1][0]*sqrt2*c[0][0]*(3+28*c[1][2]*c[0][1]+168*c[1][1]*c[0][3]-84*c[1][1]*c[0][1]+84*c[1][0]*c[0][1]-56*c[1][2]*c[0][3]-12*c[1][0]+9*c[1][1]+32*c[0][5]*c[1][2]-96*c[0][5]*c[1][1]-4*c[1][2]-32*c[0][5]-28*c[0][1]+56*c[0][3]+96*c[0][5]*c[1][0]-168*c[0][3]*c[1][0]);
        result[0][2] = -8*sqrt7*c[0][5]*c[1][2]+8*sqrt7*c[0][5]*c[1][0]+6*sqrt7*c[0][3]-4*sqrt7*c[0][5]-2*sqrt7*c[0][1]+4*sqrt7*c[1][3]*c[0][5]+1.0/2*sqrt7*c[1][2]-1.0/2*sqrt7*c[1][0]+2*sqrt7*c[1][3]*c[0][1]-5*sqrt7*c[1][2]*c[0][1]+5*sqrt7*c[1][0]*c[0][1]+12*sqrt7*c[1][2]*c[0][3]-12*sqrt7*c[0][3]*c[1][0]-6*sqrt7*c[1][3]*c[0][3];
        result[0][3] = -1.0/4*(1+4*c[1][2]-c[1][1]+12*c[1][0]*c[0][1]-12*c[1][2]*c[0][1]-8*c[0][3]*c[1][0]+8*c[1][2]*c[0][3]+8*c[1][1]*c[0][1]+8*c[0][3]-8*c[0][1]-4*c[1][0]-8*c[1][1]*c[0][3])*c[0][0]*s[1][0]*sqrt7*sqrt2;
        result[0][4] = 1.0/2*(-1+2*c[0][1]-c[1][3]+2*c[1][1]+2*c[1][3]*c[0][1]-4*c[1][1]*c[0][1])*c[0][0]*s[0][0]*sqrt35;
        result[0][5] = 1.0/4*(1+8*c[1][1]*c[0][1]+4*c[1][0]*c[0][1]-4*c[1][2]*c[0][1]-c[1][1]-8*c[1][1]*c[0][3]+8*c[1][2]*c[0][3]-8*c[0][3]*c[1][0]-8*c[0][1]+8*c[0][3])*s[0][0]*s[1][0]*sqrt7*sqrt2;
        result[0][6] = -(-1+c[1][0]+c[1][3]-4*c[0][3]-c[1][2]+4*c[0][1]-8*c[1][2]*c[0][3]+4*c[1][3]*c[0][3]+8*c[0][3]*c[1][0]+8*c[1][2]*c[0][1]-4*c[1][3]*c[0][1]-8*c[1][0]*c[0][1])*s[0][0]*c[0][0]*sqrt7;
        result[0][7] = -1.0/4*s[1][0]*sqrt2*s[0][0]*(1+3*c[1][1]-12*c[0][1]-32*c[0][5]+40*c[0][3]-96*c[0][5]*c[1][1]+12*c[1][2]*c[0][1]+36*c[1][0]*c[0][1]+32*c[0][5]*c[1][2]+96*c[0][5]*c[1][0]-120*c[0][3]*c[1][0]-40*c[1][2]*c[0][3]+120*c[1][1]*c[0][3]-36*c[1][1]*c[0][1]);
        result[0][8] = 1.0/2*s[0][0]*c[0][0]*(-1+16*c[1][3]*c[0][5]-64*c[0][5]*c[1][0]-c[1][3]-24*c[0][3]+4*c[1][2]-6*c[1][1]+16*c[0][5]+10*c[0][1]+4*c[1][0]+96*c[0][5]*c[1][1]-64*c[0][5]*c[1][2]+60*c[1][1]*c[0][1]-40*c[1][2]*c[0][1]+96*c[1][2]*c[0][3]-24*c[1][3]*c[0][3]-144*c[1][1]*c[0][3]-40*c[1][0]*c[0][1]+10*c[1][3]*c[0][1]+96*c[0][3]*c[1][0]);
        result[1][0] = -1.0 * result[0][1];
        result[1][1] = 3.0/4*c[1][1]+36*c[0][5]*c[1][2]-20*c[0][5]*c[1][0]+18*c[1][1]*c[0][3]+c[1][3]+27.0/4*c[0][1]-54*c[1][2]*c[0][3]+12*c[0][5]+30*c[0][3]*c[1][0]-18*c[0][3]+24*c[1][3]*c[0][3]-45.0/4*c[1][0]*c[0][1]-9*c[1][3]*c[0][1]+81.0/4*c[1][2]*c[0][1]-27.0/4*c[1][1]*c[0][1]-12*c[0][5]*c[1][1]-16*c[1][3]*c[0][5]-3.0/4;
        result[1][2] = 1.0/4*(3+30*c[1][1]*c[0][1]-24*c[1][1]*c[0][3]+8*c[0][3]-9*c[1][1]-10*c[0][1]+4*c[1][2]+16*c[1][2]*c[0][3]-20*c[1][2]*c[0][1])*c[0][0]*s[1][0]*sqrt14;
        result[1][3] = 1.0/4*sqrt7+sqrt7*c[1][3]-5.0/4*sqrt7*c[0][1]-5.0/4*sqrt7*c[1][1]+sqrt7*c[0][3]+4*sqrt7*c[1][3]*c[0][3]-3*sqrt7*c[1][2]*c[0][3]+3*sqrt7*c[0][3]*c[1][0]-5*sqrt7*c[1][1]*c[0][3]+9.0/4*sqrt7*c[1][2]*c[0][1]-9.0/4*sqrt7*c[1][0]*c[0][1]-5*sqrt7*c[1][3]*c[0][1]+25.0/4*sqrt7*c[1][1]*c[0][1];
        result[1][4] = 1.0/4*(-4*c[0][1]+4*c[1][1]*c[0][1]+1-c[1][1])*c[1][0]*s[0][0]*s[1][0]*sqrt70;
        result[1][5] = -1.0/4*(-1-20*c[1][1]*c[0][1]+12*c[1][0]*c[0][1]-12*c[1][2]*c[0][1]+5*c[1][1]+9*c[1][2]-4*c[1][3]+4*c[0][1]-9*c[1][0]+16*c[1][3]*c[0][1])*s[0][0]*c[0][0]*sqrt7;
        result[1][6] = -1.0/2*(-6*c[1][2]*c[0][1]+9*c[1][1]*c[0][1]+4*c[0][3]-3*c[0][1]+8*c[1][2]*c[0][3]-12*c[1][1]*c[0][3]+c[1][2])*s[0][0]*s[1][0]*sqrt14;
        result[1][7] = 1.0/4*s[0][0]*c[0][0]*(-9+15*c[1][0]-27*c[1][2]+12*c[1][3]-48*c[0][3]+48*c[0][1]+9*c[1][1]-48*c[1][1]*c[0][1]-64*c[1][3]*c[0][1]+144*c[1][2]*c[0][1]-80*c[1][0]*c[0][1]+64*c[1][3]*c[0][3]+48*c[1][1]*c[0][3]+80*c[0][3]*c[1][0]-144*c[1][2]*c[0][3]);
        result[1][8] = 1.0/4*s[1][0]*sqrt2*s[0][0]*(-36*c[1][1]*c[0][1]+120*c[1][1]*c[0][3]+96*c[0][5]*c[1][0]-96*c[0][5]*c[1][1]+32*c[0][5]*c[1][2]-120*c[0][3]*c[1][0]+36*c[1][0]*c[0][1]-40*c[1][2]*c[0][3]+12*c[1][2]*c[0][1]-32*c[0][5]-c[1][2]+40*c[0][3]-3*c[1][0]-12*c[0][1]);
        result[2][0] = 1.0 * result[0][2];
        result[2][1] = -1.0 * result[1][2];
        result[2][2] = -10*c[0][3]*c[1][0]+14*c[1][2]*c[0][3]-14*c[1][3]*c[0][3]+12*c[1][1]*c[0][3]+14*c[1][3]*c[0][1]+10*c[1][0]*c[0][1]-14*c[1][2]*c[0][1]-12*c[1][1]*c[0][1]-2*c[0][3]-5.0/2*c[1][0]+2*c[0][1]+7.0/2*c[1][2];
        result[2][3] = 1.0/4*(-1-14*c[1][1]*c[0][1]+16*c[1][0]+7*c[1][1]-28*c[1][2]+2*c[0][1]-16*c[1][0]*c[0][1]+28*c[1][2]*c[0][1])*c[0][0]*s[1][0]*sqrt2;
        result[2][4] = -1.0/2*(1+7*c[1][3]-8*c[1][1])*c[0][0]*s[0][0]*sqrt5;
        result[2][5] = -1.0/4*(-1+2*c[0][1]+7*c[1][1]-14*c[1][1]*c[0][1]-16*c[1][0]*c[0][1]+28*c[1][2]*c[0][1])*s[0][0]*s[1][0]*sqrt2;
        result[2][6] = (-1+7*c[1][2]-7*c[1][3]+6*c[1][1]+2*c[0][1]-5*c[1][0]+10*c[1][0]*c[0][1]+14*c[1][3]*c[0][1]-14*c[1][2]*c[0][1]-12*c[1][1]*c[0][1])*s[0][0]*c[0][0];
        result[2][7] = 1.0/4*s[1][0]*sqrt14*s[0][0]*(1+16*c[1][2]*c[0][3]-24*c[1][1]*c[0][3]-12*c[1][2]*c[0][1]+18*c[1][1]*c[0][1]+8*c[0][3]-3*c[1][1]-6*c[0][1]);
        result[2][8] = -1.0/2*sqrt7*c[0][0]*s[0][0]*(-1+16*c[0][3]*c[1][0]-16*c[1][2]*c[0][3]+8*c[1][3]*c[0][3]-8*c[1][3]*c[0][1]+16*c[1][2]*c[0][1]+4*c[1][0]-16*c[1][0]*c[0][1]-8*c[0][3]+8*c[0][1]-4*c[1][2]+c[1][3]);
        result[3][0] = -1.0 * result[0][3];
        result[3][1] = 1.0 * result[1][3];
        result[3][2] = -1.0 * result[2][3];
        result[3][3] = -7*c[1][3]*c[0][1]+27.0/4*c[1][1]*c[0][1]+3.0/4-27.0/4*c[1][1]-3.0/4*c[0][1]+7*c[1][3]+7.0/4*c[1][2]*c[0][1]-3.0/4*c[1][0]*c[0][1];
        result[3][4] = -1.0/4*(-3+7*c[1][1])*c[1][0]*s[0][0]*s[1][0]*sqrt10;
        result[3][5] = 1.0/4*(3-7*c[1][2]-27*c[1][1]+28*c[1][3]+3*c[1][0])*s[0][0]*c[0][0];
        result[3][6] = 1.0/2*(c[0][1]-7*c[1][1]*c[0][1]-8*c[1][0]*c[0][1]+14*c[1][2]*c[0][1]+4*c[1][0]-7*c[1][2])*s[0][0]*s[1][0]*sqrt2;
        result[3][7] = -1.0/4*sqrt7*s[0][0]*c[0][0]*(-3+3*c[1][2]+15*c[1][1]-3*c[1][0]+4*c[0][1]-12*c[1][3]-20*c[1][1]*c[0][1]+12*c[1][0]*c[0][1]-12*c[1][2]*c[0][1]+16*c[1][3]*c[0][1]);
        result[3][8] = -1.0/4*s[1][0]*sqrt2*sqrt7*s[0][0]*(-8*c[1][1]*c[0][3]+8*c[1][0]*c[0][1]+4*c[1][1]*c[0][1]-8*c[1][2]*c[0][1]-8*c[0][3]*c[1][0]+8*c[1][2]*c[0][3]+8*c[0][3]-4*c[0][1]+c[1][2]-c[1][0]);
        result[4][0] = 1.0 * result[0][4];
        result[4][1] = -1.0 * result[1][4];
        result[4][2] = 1.0 * result[2][4];
        result[4][3] = -1.0 * result[3][4];
        result[4][4] = 3.0/8+35.0/8*c[1][3]-15.0/4*c[1][1];
        result[4][5] = 1.0/4*(-3+7*c[1][1])*c[0][0]*c[1][0]*s[1][0]*sqrt10;
        result[4][6] = 1.0/4*sqrt5-1.0/2*sqrt5*c[0][1]+7.0/4*sqrt5*c[1][3]-2*sqrt5*c[1][1]-7.0/2*sqrt5*c[1][3]*c[0][1]+4*sqrt5*c[1][1]*c[0][1];
        result[4][7] = -1.0/4*c[1][0]*s[1][0]*sqrt70*c[0][0]*(4*c[1][1]*c[0][1]-4*c[0][1]+3-3*c[1][1]);
        result[4][8] = 1.0/8*sqrt35-sqrt35*c[1][3]*c[0][1]-2*sqrt35*c[1][1]*c[0][3]+sqrt35*c[1][3]*c[0][3]-sqrt35*c[0][1]+1.0/8*sqrt35*c[1][3]-1.0/4*sqrt35*c[1][1]+sqrt35*c[0][3]+2*sqrt35*c[1][1]*c[0][1];
        result[5][0] = -1.0 * result[0][5];
        result[5][1] = 1.0 * result[1][5];
        result[5][2] = -1.0 * result[2][5];
        result[5][3] = 1.0 * result[3][5];
        result[5][4] = -1.0 * result[4][5];
        result[5][5] = 3.0/4*c[0][1]+7.0/4*c[1][2]-3.0/4*c[1][0]+7*c[1][3]*c[0][1]-7.0/4*c[1][2]*c[0][1]-27.0/4*c[1][1]*c[0][1]+3.0/4*c[1][0]*c[0][1];
        result[5][6] = 1.0/2*(-1-7*c[1][1]*c[0][1]+7*c[1][1]-7*c[1][2]+c[0][1]+4*c[1][0]+14*c[1][2]*c[0][1]-8*c[1][0]*c[0][1])*c[0][0]*s[1][0]*sqrt2;
        result[5][7] = 5*sqrt7*c[1][1]*c[0][3]-4*sqrt7*c[1][3]*c[0][3]+3.0/4*sqrt7*c[1][2]-sqrt7*c[0][3]-3.0/4*sqrt7*c[1][0]+15.0/4*sqrt7*c[1][0]*c[0][1]+3.0/4*sqrt7*c[0][1]+3*sqrt7*c[1][3]*c[0][1]-15.0/4*sqrt7*c[1][2]*c[0][1]-15.0/4*sqrt7*c[1][1]*c[0][1]-3*sqrt7*c[1][0]*c[0][3]+3*sqrt7*c[1][2]*c[0][3];
        result[5][8] = -1.0/4*s[1][0]*sqrt2*sqrt7*c[0][0]*(4-c[1][0]-8*c[0][3]*c[1][0]+8*c[1][2]*c[0][3]-8*c[1][1]*c[0][3]+8*c[1][0]*c[0][1]-8*c[1][2]*c[0][1]+12*c[1][1]*c[0][1]-12*c[0][1]+8*c[0][3]+c[1][2]-4*c[1][1]);
        result[6][0] = 1.0 * result[0][6];
        result[6][1] = -1.0 * result[1][6];
        result[6][2] = 1.0 * result[2][6];
        result[6][3] = -1.0 * result[3][6];
        result[6][4] = 1.0 * result[4][6];
        result[6][5] = -1.0 * result[5][6];
        result[6][6] = 12*c[1][1]*c[0][1]-14*c[1][3]*c[0][1]-10*c[1][0]*c[0][1]+14*c[1][2]*c[0][1]+1.0/2-12*c[1][1]*c[0][3]-14*c[1][2]*c[0][3]+14*c[1][3]*c[0][3]+10*c[0][3]*c[1][0]+7.0/2*c[1][3]-2*c[0][1]-3*c[1][1]+2*c[0][3];
        result[6][7] = 1.0/2*s[1][0]*sqrt14*c[0][0]*(1-5*c[0][1]+4*c[0][3]-3*c[1][1]+3*c[1][2]+8*c[1][2]*c[0][3]-12*c[1][1]*c[0][3]-10*c[1][2]*c[0][1]+15*c[1][1]*c[0][1]);
        result[6][8] = -1.0/4*sqrt7-4*sqrt7*c[1][3]*c[0][5]+4*sqrt7*c[1][2]*c[0][1]-4*sqrt7*c[1][0]*c[0][1]-6*sqrt7*c[0][3]-5.0/2*sqrt7*c[1][3]*c[0][1]+8*sqrt7*c[0][5]*c[1][2]-8*sqrt7*c[0][5]*c[1][0]+4*sqrt7*c[0][5]+5.0/2*sqrt7*c[0][1]+1.0/4*sqrt7*c[1][3]+6*sqrt7*c[1][3]*c[0][3]-12*sqrt7*c[1][2]*c[0][3]+12*sqrt7*c[0][3]*c[1][0];
        result[7][0] = -1.0 * result[0][7];
        result[7][1] = 1.0 * result[1][7];
        result[7][2] = -1.0 * result[2][7];
        result[7][3] = 1.0 * result[3][7];
        result[7][4] = -1.0 * result[4][7];
        result[7][5] = 1.0 * result[5][7];
        result[7][6] = -1.0 * result[6][7];
        result[7][7] = -30*c[0][3]*c[1][0]+54*c[1][2]*c[0][3]-24*c[1][3]*c[0][3]-18*c[1][1]*c[0][3]+12*c[0][5]*c[1][1]+16*c[1][3]*c[0][5]+9*c[1][3]*c[0][1]+18*c[0][3]+9.0/4*c[1][2]-27.0/4*c[0][1]-12*c[0][5]-36*c[0][5]*c[1][2]+27.0/4*c[1][1]*c[0][1]-81.0/4*c[1][2]*c[0][1]-5.0/4*c[1][0]+20*c[0][5]*c[1][0]+45.0/4*c[1][0]*c[0][1];
        result[7][8] = 1.0/4*s[1][0]*sqrt2*c[0][0]*(4-84*c[1][1]*c[0][1]+28*c[1][2]*c[0][1]-168*c[0][3]*c[1][0]-96*c[0][5]*c[1][1]-9*c[1][0]+32*c[0][5]*c[1][2]-28*c[0][1]+56*c[0][3]+12*c[1][1]-32*c[0][5]-3*c[1][2]+96*c[0][5]*c[1][0]+168*c[1][1]*c[0][3]+84*c[1][0]*c[0][1]-56*c[1][2]*c[0][3]);
        result[8][0] = 1.0 * result[0][8];
        result[8][1] = -1.0 * result[1][8];
        result[8][2] = 1.0 * result[2][8];
        result[8][3] = -1.0 * result[3][8];
        result[8][4] = 1.0 * result[4][8];
        result[8][5] = -1.0 * result[5][8];
        result[8][6] = 1.0 * result[6][8];
        result[8][7] = -1.0 * result[7][8];
        result[8][8] = 1.0/8+48*c[0][7]*c[1][1]+8*c[1][3]*c[0][7]+8*c[1][0]*c[0][1]+8*c[1][2]*c[0][1]-40*c[0][3]*c[1][0]-12*c[1][1]*c[0][1]-2*c[1][3]*c[0][1]+10*c[1][3]*c[0][3]-32*c[0][7]*c[1][0]-32*c[1][2]*c[0][7]-40*c[1][2]*c[0][3]-16*c[0][5]+8*c[0][7]-2*c[0][1]+3.0/4*c[1][1]+10*c[0][3]-16*c[1][3]*c[0][5]+64*c[0][5]*c[1][2]-96*c[0][5]*c[1][1]+64*c[0][5]*c[1][0]+60*c[1][1]*c[0][3]+1.0/8*c[1][3];

        return result;
    }
    //--------------------------------------------------------------------------

}
