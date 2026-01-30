package Polynom3D;

/*
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
 * (c) Tymofii Nikolaienko, 2014
 */

// A simple class to store an expansion of a primitive cartesian x^n·y^m·z^k over
// pure spherical harmonics
public class YLM_Series {
    public int primary_L;
    public int[] Ls;
    public int[] Ms;
    public int[] indexes;      // an auxiliary array
    public double[] coefs;
    //--------------------------------------------------------------------------
    // Some constructors:
    //--------------------------------------------------------------------------
    public YLM_Series(int LMain, int[][] LMs /*1-st index: term number, 2-nd: 0(L), 1(M) */, double[] Coefs) {
        primary_L = LMain;
        Ls = new int[LMs.length];
        Ms = new int[LMs.length];
        indexes = new int[LMs.length];
        for (int k=0; k<LMs.length; k++) {
            Ls[k] = LMs[k][0];
            Ms[k] = LMs[k][1];
        }
        coefs = Coefs.clone();
    }
    //--------------------------------------------------------------------------
    public YLM_Series(YLM_Series source) {
        this.primary_L = source.primary_L;
        this.Ls = source.Ls.clone();
        this.Ms = source.Ms.clone();
        this.coefs = source.coefs.clone();
        this.indexes = source.indexes.clone();
    }
    //--------------------------------------------------------------------------
    public YLM_Series(int NTermsMax) {
        Ls = new int[NTermsMax];
        Ms = new int[NTermsMax];
        indexes = new int[NTermsMax];
        coefs = new double[NTermsMax];
    }
}
