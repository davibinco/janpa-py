package JGints;

/**
 * A small file containing only a class declaration
 *
 * Version: 17.05.2014 / 16.Jan.2014 / 26.Oct.2013
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

import java.io.*;

public class MO {
    public double Energy;
    public double Occupancy;
    public int Spin; // spin*2: +1='alpha' of -1='beta'
    public double[] BS_Coefs;
    public int reserved = 0; // reserved, should be zero
    public MO _next = null;
    //--------------------------------------------------------------------------
    // a constructor that does nothing
    public MO() { }
    //--------------------------------------------------------------------------
    // a constructor that copies everything from source
    public MO(MO source) {
        BS_Coefs = source.BS_Coefs.clone();
        Energy = source.Energy;
        Occupancy = source.Occupancy;
        Spin = source.Spin;
        _next = source._next; // shouldn't it be here?
    }
    //--------------------------------------------------------------------------
    /** Save all data to @param out stream:
     * double E | double Occupancy | int Spin | int NCoefs | double[] coefs | int reserved
     */
    public void SaveToDataStream(DataOutputStream out) throws Exception {
        out.writeDouble(Energy);
        out.writeDouble(Occupancy);    
        out.writeInt(Spin);
        out.writeInt(BS_Coefs.length);
        for (int c=0; c<BS_Coefs.length; ++c) out.writeDouble(BS_Coefs[c]);
        out.writeInt(reserved);
    }
    //--------------------------------------------------------------------------
    /** Reads all data from @param out stream:
     * double E | double Occupancy | int Spin | int NCoefs | double[] coefs | int reserved
     */
    public void ReadFromDataStream(DataInputStream in) throws Exception {
        Energy = in.readDouble();
        Occupancy = in.readDouble();
        Spin = in.readInt();
        BS_Coefs = new double[in.readInt()];
        for (int c=0; c<BS_Coefs.length; ++c) BS_Coefs[c] = in.readDouble();
        reserved = in.readInt();
    }
    //--------------------------------------------------------------------------
}
//==============================================================================
