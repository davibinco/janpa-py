package CLPO;
/**
 * A class representing a set of hybrid orbitals for a given atom
 *
 * Created: 13.Jan.2018
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

import Jama.*;
import MatrixHelper.EigenEngine;
//import MatrixHelper.LDLtransform;
import MatrixHelper.Sorter;
        
/**
 * Represents a set of NAO->hybrids transformation on a single atom
 */
public class AtomicHybrids {
    public int[] NAO_indices;
    public Matrix U; // NAO->hybrids transformation: [local_NAO_index][local_hybrid_index]
    private Matrix U_saved;
    
    public int nValidHybrids = 0;
    public int[] friendAtomIndex;   // together with friendHybridIndex[] gives an    
    public int[] friendHybridIndex; // 'address' of the 'friend vector'
    private int[][] friendRefs_saved; // [0] = saved friendAtomIndex, [1] = saved friendHybridIndex
    //public int[] globalHybrIndices; // a 'global' address of the hybrid
    

    /**
     * Creates the U matrix
     * Must be called with the true array @param NAO_indices containing
     * the NAO indices belonging to the given atom
     */
    AtomicHybrids(int[] NAO_indices) {
        this.NAO_indices = NAO_indices;
        U = new Matrix( NAO_indices.length, NAO_indices.length );
        nValidHybrids = 0;
        friendAtomIndex = new int[NAO_indices.length];
        friendHybridIndex = new int[NAO_indices.length];        
        friendRefs_saved = new int[2][];        
    }
    //--------------------------------------------------------------------------
    /**
     * Copies this.U into this.U_saved by deep copy
     */
    public void backup_U() {
        this.U_saved = U.copy();
    }
    
    /**
     * restores the saved U by deep copy
     */
    public void restore_U() {
        U = this.U_saved.copy();
    }
    //--------------------------------------------------------------------------
    /**
     * backups friendAtomIndex and friendHybridIndex via deep copy
     */
    public void backup_bonding() {
        int n = friendAtomIndex.length;
        friendRefs_saved[0] = new int[n];
        System.arraycopy(friendAtomIndex, 0, friendRefs_saved[0], 0, n);
        
        n = friendHybridIndex.length;
        friendRefs_saved[1] = new int[n];
        System.arraycopy(friendHybridIndex, 0, friendRefs_saved[1], 0, n);
    }
    
    /**
     * restores saved data back into friendAtomIndex and friendHybridIndex via deep copy
     */    
    public void restore_bonding() {        
        int n = friendRefs_saved[0].length;
        friendAtomIndex = new int[n];
        System.arraycopy(friendRefs_saved[0], 0, friendAtomIndex, 0, n);
        
        n = friendRefs_saved[1].length;
        friendHybridIndex = new int[n];
        System.arraycopy(friendRefs_saved[1], 0, friendHybridIndex, 0, n);
        
    }    
    //--------------------------------------------------------------------------
    /**
     * @returns @param i-th column of the U matrix
     */
    public Matrix getHybrid(int i) {
        return U.getMatrix(0, NAO_indices.length-1, i, i);
    }
    //--------------------------------------------------------------------------
    
    /**
     * @returns scalar product of @param i-th hybrid and the 
     * column-vector @param v
     * The method avoids calling matrix/memory reorganization methods
     */
    public double hybridScalarMul(int i, Matrix v) {
        double result = 0.0;
        int nComponents = v.getRowDimension();
        for (int k=0; k<nComponents; k++)
            result += U.get(k, i) * v.get(k, 0);
        return result;
    }    
    //--------------------------------------------------------------------------
    
    /**
     * Inserts the given vector @param vec into the U matrix and increments 
     * nValidHybrids. 
     * @returns false in case of 'oversize'
     */
    public boolean appendHybrid(double[] vec) {
        if (nValidHybrids ==  NAO_indices.length)
            return false;
        for (int i=0; i<vec.length; i++)
            U.set(i, nValidHybrids, vec[i]);
        
        nValidHybrids++;
        return true;
    }
    //--------------------------------------------------------------------------
    
    /**
     * Compares the smallest eigenvalue of the hybrids overlap matrix
     * with the threshold value given in @param eigS_threshold
     * @returns true if min(eigenvals) > @param eigS_threshold
     * If @param eigenvalsS is not null, eigenvalsS[0] gets the array
     * of eigenvalues of the overlap matrix
     */
    public boolean testLinIndep(double eigS_threshold, double[][] eigenvalsS) {
        if ((nValidHybrids <= 1) && (eigenvalsS == null)) {
            return true;
        }
        Matrix Usmall = U.getMatrix(0, NAO_indices.length - 1, 0, nValidHybrids - 1);
        Matrix S = Usmall.transpose().times(Usmall);        
        
        if ((nValidHybrids <= 1) && (eigenvalsS != null)) {
            eigenvalsS[0] = new double[]{ S.get(0, 0)  }; // return a true scalar
            // product as an eigenvalue
            return true;
        }
        
        Object[] eigs = Sorter.sorted_eigenData(S.eig(), -1);
        double[] vals = (double[])eigs[0];
        if (eigenvalsS != null)
            eigenvalsS[0] = vals;
        return vals[vals.length - 1] > eigS_threshold;
    }
    //--------------------------------------------------------------------------
    
    /**
     * Adds @param vec into the U matrix, checks if the vectors in U are linearly 
     * independent (usind @param eigS_threshold) and @returns true if yes,
     * but otherwise returns false and removes the added vector from U
     * by decreasing nValidHybrids
     * In other words, U gets changed only if the returned value was true
     */
    public boolean appendHybridIfLinIndep(double eigS_threshold, double[] vec) {
        if (! appendHybrid(vec) )
            return false; // not that in this case nothing has changed
        // if we're here, the new vector has been inserted into U, and nValidHybrids
        // has been incremented
        if (testLinIndep(eigS_threshold, null))
            return true; // the vector vec is still in U, and nValidHybrids is valid
        else {
            // we must remove vec from U
            nValidHybrids--;
            return false;
        }
    }
    //--------------------------------------------------------------------------
    
    /**
     * Symmetrically orthogonalizes the hybrids saved in U
     */
    public void symmetrOrth() {
        if ( nValidHybrids <= 1 ) {
            return ; // nothing to do!
        }
        Matrix Usmall = U.getMatrix(0, NAO_indices.length - 1, 0, nValidHybrids - 1);
        // perform Symmetrical orthogonalization and save its result back into U
        U.setMatrix(0, NAO_indices.length - 1, 0, nValidHybrids - 1,
                EigenEngine.symmetrOrth(Usmall) );
    }
    //--------------------------------------------------------------------------
    /**
     * Finds and saves into the U matrix the vectors necessary to make the
     * number of hybrids equal to the number of NAOs for this atom.
     * The newly created hybrids are build in a way that they orthonormal
     * and orthogonal to the hybrids existing in U before this method is called.
     * What is more, the new hybrids are being built so that they were as
     * close to eigenvectors of @param D_AA as possible.
     * In case if @param D_AA == null, the new hybrids are just orthonormal
     * and orthogonal to the ones existing in U before the method is called.
     */
    /*public void complementHybridBasis(Matrix D_AA) {
        if ( nValidHybrids == 0 ) {
            // if there were no hybrids before => just create a unitary
            // matrix (if D_AA is null), or place all eigenvectors of D_AA
            // into U
            if (D_AA == null) {
                U = Matrix.identity(NAO_indices.length , NAO_indices.length);
            } else {
                U = D_AA.eig().getV();
            }
        } else {
            Matrix Usmall = U.getMatrix(0, NAO_indices.length - 1, 0, nValidHybrids - 1);
            // make a projection matrix
            Matrix P = Usmall.times(Usmall.transpose());
            // now the existing hybrids correspond to nValidHybrids largest
            // eigenvalues of P
            Object[] eig = Sorter.sorted_eigenData(P.eig(), -1);
            Matrix vecs = (Matrix)eig[1];
            Matrix supplBas = vecs.getMatrix(0, vecs.getRowDimension()-1, 
                    nValidHybrids, vecs.getColumnDimension()-1); // skip the first
            // nValidHybrids eigenvectors (from 0-th to nValidHybrids-1 - th)
            // corresponding to nValidHybrids largest eigenvalues
            Matrix D_newBas = 
                    LDLtransform.TransformSymmetricMatrixToNewBasis_LDLT(D_AA, 
                            supplBas.transpose()); // makes supplBas^T . D_AA . supplBas
            // now get eigenvectors of D_newBas and use them to prepare
            // (orthonormal) linear combinations of the supplBas vectors as
            // Vecs = supplBas.X, where X is eigenvectors of D_newBas
            
            eig = Sorter.sorted_eigenData(D_newBas.eig(), -1); // sort eigenvalues descendingly for nice print
            supplBas = supplBas.times( (Matrix)eig[1] );
            // now supplBas contains exactly the vectors we've wanted to build
            // Just insert them into U:
            U.setMatrix(0, NAO_indices.length-1, 
                    nValidHybrids, NAO_indices.length-1, 
                    supplBas);
        }
        // Finally, properly fill in the arrays with the hybrid properties
        while ( nValidHybrids < NAO_indices.length ) {
            friendAtomIndex[nValidHybrids] = -1;
            friendHybridIndex[nValidHybrids] = -1;
            nValidHybrids++;
        }
        // done!
    } */   
    //--------------------------------------------------------------------------
    
}
