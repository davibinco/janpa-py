package JGints;

import java.io.*;

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
 * @author (c) Tymofii Nikolaienko, 2014
 * @version: 17.05.2014
 */
public class RadialPartOfBasisFunction {

    public double[] Exponents = null;
    public double[] Coefs = null;
    public int CenterID = -1; // 1-based identifier of the nuclei !!!
    public int LUsedWith = -1;
    public int Addit_r_power = 0;
    //--------------------------------------------------------------------------
    public void alloc_mem(int nCoefs) {
        Exponents = new double[nCoefs];
        Coefs = new double[nCoefs];
    }
    //--------------------------------------------------------------------------
    public RadialPartOfBasisFunction CreateClone() {
        RadialPartOfBasisFunction result = new RadialPartOfBasisFunction();
        result.Exponents = Exponents.clone();
        result.Coefs = Coefs.clone();
        result.CenterID = CenterID;
        result.LUsedWith = LUsedWith;
        result.Addit_r_power = Addit_r_power;
        return result;
    }
    //--------------------------------------------------------------------------
    /** Save all data to @param out stream:
     * int CenterID | int LUsedWith | int Addit_r_power | int NExponents | double[] exponents | int NCoefs | double[] coefs
     */
    public void SaveToDataStream(DataOutputStream out) throws Exception {
        out.writeInt(CenterID);
        out.writeInt(LUsedWith);
        out.writeInt(Addit_r_power);
        out.writeInt(Exponents.length);
        for (int e=0; e<Exponents.length; ++e) out.writeDouble(Exponents[e]);
        out.writeInt(Coefs.length);
        for (int c=0; c<Coefs.length; ++c) out.writeDouble(Coefs[c]);
    }
    //--------------------------------------------------------------------------
    /** Reads all data from @param out stream:
     * int CenterID | int LUsedWith | int Addit_r_power | int NExponents | double[] exponents | int NCoefs | double[] coefs
     */
    public void ReadFromDataStream(DataInputStream in) throws Exception {
        CenterID = in.readInt();
        LUsedWith = in.readInt();
        Addit_r_power = in.readInt();
        Exponents = new double[ in.readInt() ];
        for (int e=0; e<Exponents.length; ++e) Exponents[e] = in.readDouble();
        Coefs = new double[  in.readInt() ];
        for (int c=0; c<Coefs.length; ++c) Coefs[c] = in.readDouble();
    }
    //--------------------------------------------------------------------------
}
