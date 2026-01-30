package JGints;

/**
 * A class representing an atomic center defined by its coordinates,
 * nucleus charge, unique ID and user-readable name
 *
 * Version: 17.05.2014 / 26.Oct.2013
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

public class AtomicCenter {
    public String Name;
    public int ID;
    public double Z;
    public double[] R0 = new double[3];
    public AtomicCenter _next = null;      // Let's be as simple as possible :)
    //--------------------------------------------------------------------------
    // an empty constructor
    public AtomicCenter() { }
    //--------------------------------------------------------------------------
    // a constructor which creates 'true' copy
    public AtomicCenter(AtomicCenter source ) {
        Name = new String(source.Name);
        ID = source.ID;
        Z = source.Z;
        R0 = source.R0.clone();
        _next = source._next;
    }
    //--------------------------------------------------------------------------
    /** Save all data to @param out stream:
     * int ID | (int) NameLength | (char[] name) | double Z | double X | double Y | double Z
     */
    public void SaveToDataStream(DataOutputStream out) throws Exception {
        out.writeInt(ID);
        out.writeInt(Name.length());
        out.writeChars(Name);
        out.writeDouble(Z);
        out.writeDouble(R0[0]); // X
        out.writeDouble(R0[1]); // Y
        out.writeDouble(R0[2]); // Z
    }
    //--------------------------------------------------------------------------
    /** Reads all data from @param out stream:
     * int ID | (int) NameLength | (char[] name) | double Z | double X | double Y | double Z
     */
    public void ReadFromDataStream(DataInputStream in) throws Exception {
        ID = in.readInt();
        char[] name = new char[ in.readInt() ];
        for(int i=0; i<name.length; ++i) name[i] = in.readChar();
        Name = String.valueOf(name);
        name = null;
        Z = in.readDouble();
        R0 = new double[3];
        for (int mu=0; mu<3; ++mu) R0[mu] = in.readDouble();
    }
    //--------------------------------------------------------------------------
    /**
     * @returns the distance from this atomic center to point r[]
     */
    public double distanceTo(double[] r) {
        double result = 0;
        for (int mu=0; mu<3; ++mu) result += (R0[mu]-r[mu])*(R0[mu]-r[mu]);
        return Math.sqrt( result );
    }
}
//==============================================================================
