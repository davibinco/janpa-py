package MatrixHelper;

/**
 * Some matrix manipulations which are useful for eigenproblem solving
 * 
 * Version: 24.Jan.2015
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
import java.io.*;
import ProgramOptions.*;


public class EigenEngine {

    public static PrintStream out = System.out;
    //--------------------------------------------------------------------------

    // Some options used _mainly_ by TransformMatrixToNewBasis() method:
    //public OptionParameter No_Enhanced_Matrix_Transform = new OptionParameter(null, "", false);
    //public OptionParameter MatrixTransform_Force_Symmetric = new OptionParameter(null, "", false);
    public boolean glbPrint = false;
    //--------------------------------------------------------------------------
    /**
     * @returns @param Q . @param A . (@param Q)^t where A is known to be symmetric
     * Note that Q can be rectangular (but A is square since it's symmetric!)
     * Ver. 11.Nov.2018
     * 
     */
    public static Matrix TransformSymmetricMatrixToNewBasis(Matrix A, Matrix Q) {
        // just use a 'semidirect' method: A = U + U^T => Anew = Q.U.Q^T + Q.U^T.Q^T = 
        //  = (Q.U).Q^T + Q.(Q.U)^T, where L is lower-triangular with Uii = Aii/2.0
        // Computing B = Q.U once requieres n*(n+(n-1)+...+1) ~ O(n^3/2) multiplications,
        // while computing B.Q^T takes O(n^3) multiplications. The final cost is thus O(1.5*n^3)
        // multiplications
        double[][] Aarr = A.getArray();
        int nDimA = Aarr.length;             // A: nDimA x nDimA ; U: nDimA x nDimA
        double[][] Qarr = Q.getArray();      // Q: nDimQ x nDimA
        int nDimQ = Qarr.length;                
        double[] bi = new double[nDimA];     // B = Q.U => B: nDimQ x nDimA        
        double[][] BQt = new double[nDimQ][nDimQ];  // BQt = B.Q^T => BQt: nDimQ x nDimQ

        for(int i=0; i<nDimQ; i++) {
            double[] ui = Qarr[i];
            for(int j=0; j<nDimA; j++) {
                bi[j] = 0; // !!!
                for (int k=0; k<j; k++) { // use only upper triangle of A, i.e., U[a][b]==0 @ a>b
                    bi[j] += ui[k] * Aarr[k][j]; // nonzero for k<j only
                }
                bi[j] += ui[j] * Aarr[j][j] / 2.0; // div2 should not degrade accuracy
            }
            
            // compute B . U^t
            double[] ri = BQt[i];
            for(int j=0; j<nDimQ; j++) {
                double[] qj = Qarr[j]; // good, since we'll have a sequential access
                for (int k=0; k<nDimA; k++) { 
                    ri[j] += bi[k] * qj[k];
                }
            }
        }
        
        // Use B.Q^t to create the result;
        // (we even don't need to allocate new memory!)
        for(int i=0; i<nDimQ; i++) {
            double[] but_i = BQt[i];
            for(int j=0; j<i; j++) {
                but_i[j] += BQt[j][i];
            }
            but_i[i] *= 2.0;
        }         
        for(int i=0; i<nDimQ; i++) {
            double[] but_i = BQt[i];
            for(int j=0; j<i; j++) {
                BQt[j][i] = but_i[j];
            }
        }         
        /*
        double[][] resArr = result.getArray();
        double tmp;
        for(int i=0; i<n; i++) {
            for(int j=0; j<i; j++) {
                tmp = BUt[i][j] + BUt[j][i];
                resArr[i][j] = tmp;
                resArr[j][i] = tmp;
            }
            resArr[i][i] = BUt[i][i] * 2.0 ;                
        }         
        */

        return new Matrix(BQt);        
    }
    //--------------------------------------------------------------------------
    
    /** Transforms matrix A to a new basis given by a transformation matrix U
     * @Returns U * A * U^T employing symmetric properties of A if isSymmetric flag
     * is set and
     *
     * Version: 11.Nov.2018 (LDLt replaced with a faster a more reliable TransformSymmetricMatrixToNewBasis() )
     * Created: 14.Oct.2014, minor changes: 24.Jan.2015
     * 
     */ 
    public Matrix TransformMatrixToNewBasis(Matrix A, Matrix U, boolean isSymmetric) {
        Matrix result;

        // Make it work even if input matrix has zero size!
        if (A.getArray().length == 0)
            return new Matrix(0,0);

        // now - a typical case, when the input matrix is not empty
        boolean do_direct = !isSymmetric; // whether to calculate U.A.U^T directly with 2 matrix multiplications
        // Do not use LDL-way if No_Enhanced_Matrix_Transform option is set to true:
        //if (No_Enhanced_Matrix_Transform.get_boolean()) do_direct = true; // 11.Nov.2018: removed // 07.Nov.2018: bugfix (was ' = false')

        if (do_direct) {
            result = U.times(A).times(U.transpose());           // requires 2*N^3 multiplications
            /*if (isSymmetric && MatrixTransform_Force_Symmetric.get_boolean())
                Matrix_Symmetrize_Check(result, "A matrix from TransformMatrixToNewBasis()");*/
        } else {
            // Special case of a symmetric matrix:
            result = TransformSymmetricMatrixToNewBasis(A, U);
            
            /*
            result = LDLtransform.TransformSymmetricMatrixToNewBasis_LDLT(A, U);
            if (result == null) {
                out.println();
                WarningManager.warning_printf("WARNING! LDL^T decomposition failed. Try using "+No_Enhanced_Matrix_Transform.Name+
                        " and "+MatrixTransform_Force_Symmetric.Name+" options in the command line.");
                out.println("Now trying a 'direct' method...");
                boolean old_option = MatrixTransform_Force_Symmetric.get_boolean();
                // force symmetric
                MatrixTransform_Force_Symmetric.set_data(true);
                // apply a 'direct' method*/
                //result = TransformMatrixToNewBasis(A,U, false /*!!!*/);
                /*//restore previous option value
                MatrixTransform_Force_Symmetric.set_data(old_option);
            }*/
        }
        return result;
    }
    //--------------------------------------------------------------------------
    
    final static public double machiene_eps = 1.0E-15; // upd.11.Nov.2018
    public static void Matrix_Symmetrize_Check(Matrix S, String caller_text) {
        for (int i=0; i<S.getArray().length; i++)
            for (int j=i; j<S.getArray().length; j++) {
                if (Math.abs(S.get(i, j) - S.get(j, i)) > machiene_eps)
                    WarningManager.warning_printf("WARNING: non-symmetric matrix found in %s ((%d,%d) = %.20E, (%d,%d) = %.20E), diff = %.5E%n",
                                    caller_text,  i+1,j+1,S.get(i, j)  ,j+1,i+1,S.get(j, i), S.get(i, j)-S.get(j, i));
                double sij = (S.get(i, j) + S.get(j, i))/2.0;
                S.set(j, i, sij);
                S.set(i, j, sij);
            }
    }
    
    //--------------------------------------------------------------------------
    /**  Calculates the square root of the matrix
     */
    public static Matrix Matrix_SQRT(Matrix S, boolean do_inverse) {
        // a trick for JAMA symmetry check (any overlap matrix is anyway intrinsicaly symmetrical!!!)
        Matrix_Symmetrize_Check(S, "Matrix_SQRT");
        EigenvalueDecomposition X = new EigenvalueDecomposition(S);
        Matrix Decomp05 = X.getD(); // Eigenvalues
        // replace eigenvalues with their square roots
        /*
        for (int i=0; i<S.getArray().length; i++) {
            Decomp05.set(i, i, Math.sqrt(Decomp05.get(i, i) ));
        }

        Matrix S05 = X.getV().times(Decomp05).times(X.getV().transpose()); // S^(1/2)   -- ineffective due to many multiplications with zeros
         *
         */
        // more efficient way:
        int N = S.getRowDimension();
        double[] sqrts = new double[N];
        if (do_inverse)
            for (int i=0; i<N; ++i) sqrts[i] = 1.0 / Math.sqrt(Decomp05.get(i, i) );    // generate S^(-1/2)
        else
            for (int i=0; i<N; ++i) sqrts[i] = Math.sqrt(Decomp05.get(i, i) );          //generate S^(+1/2)

        double[][] vecs = X.getV().getArray();

        Matrix S05 = new Matrix(N, N, 0.0);
        double tmp;
        for (int i=0; i<N; ++i)
            for (int j=0; j<N; ++j) {
                tmp = 0;
                for (int k=0; k<N; ++k) tmp += vecs[i][k] * sqrts[k] * vecs[j][k];
                S05.set(i, j, tmp);
            }

        return S05;
    }
    //--------------------------------------------------------------------------
    /**  Calculates the S^(-1/2)
     */
    public Matrix Matrix_minus05(Matrix S) {
        // check whether the input matrix is not empty:
        if (S.getArray().length == 0)
            return new Matrix(0,0);

        // Now - do an actual job
        // a trick for JAMA symmetry check (any overlap matrix is anyway intrinsicaly symmetrical!!!)
        Matrix_Symmetrize_Check(S, "Matrix_minus05");
        EigenvalueDecomposition X = new EigenvalueDecomposition(S);
        Matrix Decomp05 = X.getD(); // Eigenvalues
        // replace eigenvalues with their square roots
        /*
        for (int i=0; i<Decomp05.getArray().length; i++) {
            //out.println(" _eigenvalue: "+Decomp05.get(i, i));
            Decomp05.set(i, i, 1/Math.sqrt(Decomp05.get(i, i) ));
        }
        Matrix S05 = X.getV().times(Decomp05.times(X.getV().transpose())); // S^(1/2)
         *
         */
        return Matrix_SQRT(S, true);
    }
    //--------------------------------------------------------------------------
    static boolean EVD_DBG = false;
    /** solves A.x = lambda.S.x
     * @returns: result[0] = eigenVALUES, result[1] = eigenVECTORS
     */
    public Matrix[] Generalized_EVD_SymmMatr(Matrix A, Matrix S){
        Matrix[] result = new Matrix[2];
        // make S^(1/2)
        Matrix S05 = Matrix_SQRT(S, false);
        // now
        //    A . x = L . S . x   =>
        // => A . S^(-1/2) . S^1/2 . x = L . S^1/2T . S^1/2 . x
        // => S^1/2T^(-1) . A . S^(-1/2) . y = L . y,
        // where L is one of the eigenvalues, and y = S^1/2 . x, so that
        // original eigenvectors are: x = S^(-1/2) . y
        Matrix S05inv = S05.inverse();
        Matrix A_new = TransformMatrixToNewBasis(A, S05inv, true);
        //Matrix A_new = ((S05.transpose().inverse()).times(A)).times(S05inv);

        // solve eigenvalue problem for A_new
        EigenvalueDecomposition AEVD = A_new.eig();
        // get eigenvalues
        result[0] = AEVD.getD();
        if (EVD_DBG) {
            out.println("EigenVALUES: ");
            result[0].print(13, 7);
            out.println("Trace = "+result[0].trace());
        }
        // calculate eigenvectors
        // AEVD.getV() = y = S^(1/2).x => x = S^(-1/2).(AEVD.getV())
        result[1] = S05inv.times(AEVD.getV());
        // print eigenvectors
        if (EVD_DBG) {
            out.println("EigenVECTORS: ");
            result[1].print(13, 7);
        }
        if (glbPrint) {
            // check
            result[1].print(13, 7);
            (result[1].transpose()).times(S).times(result[1]).print(13, 7);
            //(result[1]).times(result[1].transpose()).print(13, 7);
            //(result[1].transpose()).times(result[1]).print(13, 7);
            S.print(13,7);
        }
        return result;
    }
    //--------------------------------------------------------------------------

    /**
     * @Returns the result of SVD decomposition of the matrix @param A
     * as A = result[0].result[1].result[2], where all matrices have appropriate
     * sizes and can be multiplied directly
     * The eigenvals argument can either be null, or have its 1-st dimension
     * equal to 1; in this case the array of eigenvalues is written into eigenvals[0][]
     */
    public static Matrix[] SVD(Matrix A, double[][] eigenvals) {
        SingularValueDecomposition svd;
        boolean isTransposed = false;
        
        // Jama works for rows >= columns        
        if (A.getRowDimension() >= A.getColumnDimension()) {
            svd = new SingularValueDecomposition(A);                            
        } else {
            svd = new SingularValueDecomposition(A.transpose());
            isTransposed = true;
        }
        // For an m-by-n matrix A (m >= n) Jama makes  A = U*S*V' with
        // m-by-n U, n-by-n diagonal S, n-by-n V 
        // and orders the 'eigenvalues' descendingly
        //System.out.println(isTransposed + " " + A.getRowDimension() + " " + A.getColumnDimension());
        
        
        //Matrix S = new Matrix(svd.getU().getColumnDimension(), svd.getV().getRowDimension());
        int nEigenval = svd.getSingularValues().length;
        /*for (int i = 0; i < nEigenval; i++)
            S.set(i, i, svd.getSingularValues()[i]);*/
        Matrix S = svd.getS();

        // 'export' eigenavlues as an array, if requested
        if (eigenvals != null) {
            eigenvals[0] = new double[ nEigenval ];
            System.arraycopy( svd.getSingularValues(), 0, eigenvals[0], 0, nEigenval);
        }

        // form and return result
        if (! isTransposed  ){
            return new Matrix[]{ svd.getU(), S, svd.getV().transpose()  };
        } else {
            // A^T = U*S*V' => A = (U*S*V')^T = V . S^T . U^T
            return new Matrix[]{ svd.getV(), S.transpose(), svd.getU().transpose()  };
        }

    }
    //--------------------------------------------------------------------------
    /**
     * @returns V.(V^T.V)^(-1/2), where V = @param Vecs
     * i.e., performs Symmetrical orthogonalization of the column-vectors
     * of @param Vecs by SVD decomposition
     * 
     */
    public static Matrix symmetrOrth(Matrix Vecs) {
        Matrix[] svd = EigenEngine.SVD(Vecs, null); // makes Vecs = (A==svd[0]).svd[1].(svd[2]==B)
        
        /*System.out.print("thsi is symmetrOrth");
        for (int i=0; i<svd[1].getColumnDimension(); i++) {
            System.out.printf("%.5e ", svd[1].get(i, i));
        } System.out.println();*/
        // Now, 
        // S = V^T.V = (A.L.B)^T.(A.L.B) = B^T . L . A^T . A . L . B = B^T . L^2 . B = 
        //   = B^T.L.B . B^T.L.B  => S^(1/2) = B^T.L.B (symmetrical square root,
        // provided that each element of L is non-negative)
        // Note that if V.shape == m x n, then V^T.V.shape == n x n,
        // and also the shapes of the matrices produced by SVD() are A = m x n
        // (not that A is rectangular!) and B = n x n (square)
        //  => S^(-1/2) = B^T.L^(-1).B
        // -- EigenEngine.Matrix_SQRT(Vecs.transpose().times(Vecs) , true).print(13, 7);
        //  => Vnew = V.S^(-1/2) = A.L.B . B^T.L^(-1).B = A.B   => Vnew = A.B
        /*int n = svd[1].getColumnDimension();
        Matrix d = new Matrix(n,n);
        for (int i=0; i<n; i++) {
            double L = svd[1].get(i,i);
            d.set(i, i, L / Math.sqrt(L*L + 1e-8)  );
        }*/
        //return svd[0].times( d  ).times(svd[2]);
                
        return svd[0].times(svd[2]);
    }
    //--------------------------------------------------------------------------

}
