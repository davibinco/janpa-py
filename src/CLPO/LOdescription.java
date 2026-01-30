package CLPO;
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
 */

import Jama.Matrix;

/**
 * Data structure containing a 'public' information about the localized 
 * orbitals and their constituting hybrids
 
 */
public class LOdescription { // 'LOs' stand for 'Localized Orbitals'

    /** Transformation matrix: 1-st index: NAO, 2-nd index: hybrid */
    public Matrix NAO_to_Hybrids;
    /** Transformation matrix: 1-st index: LPO, 2-nd index: hybrid */
    public Matrix  LO_to_Hybrids;
    /** Human-readable short description of LO */
    public String[] LO_labels;
    final static int LO_type_RY = 0;
    final static int LO_type_LP = 1;
    final static int LO_type_BD = 2;
    final static int LO_type_NB = 3;
    /** Type of LO assigned when it was being created: RY/LP/BD/NB; most useful
     for a posteriori distinguishing BDs from NBs; the type identification 
     constants are also declared in this class */
    public int[] LO_types;
    /** Human-readable short description of a constituent hybrid orbital */
    public String[] Hybrid_labels;
    /** [i] = 0-based index of the i-th hybrid 'parent' atom */
    public int[] hostAtomOfHybrid;
    /** [LOindex][0...numberOfHybridsNeededForThisLO-1] */
    public int[][] hybridsOfLO;
    /** number of BD orbitals per atomic pair (optional); symmetrical matrix structure is recommended */
    public int[][] BDperAtomicPair = null;
    
    //public int[][] atomicHybrids_to_globalHybrid; // [atomId][0...nHybr 'local' hybrid id]
    //public int[] partnerHybridGlobalIndex;
    //public int[] hostLOindex;               // [i] = global id of the LO which uses the i-th hybrid
    
    /**
     * Creates the object AND allocates all matrices and arrays
     */
    LOdescription(int nNAOs) {
        this.LO_to_Hybrids = new Matrix(nNAOs, nNAOs);
        this.NAO_to_Hybrids = new Matrix(nNAOs, nNAOs);
        this.LO_labels = new String[nNAOs];
        this.LO_types = new int[nNAOs];
        this.hostAtomOfHybrid = new int[nNAOs];
        //this.partnerHybridGlobalIndex = new int[nNAOs];
        //this.hostLOindex = new int[nNAOs];
        this.hybridsOfLO = new int[nNAOs][];
    }
    
    
}
