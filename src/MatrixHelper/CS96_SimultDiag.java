package MatrixHelper;

import Jama.Matrix;

/**
 * Simlutanious (approximate) diagonalization of a set of matrices by
 * [ Cardoso, J. F., & Souloumiac, A. (1996). 
 * Jacobi angles for simultaneous diagonalization. 
 * SIAM journal on matrix analysis and applications, 17(1), 161-164.]
 * @author timAdmin
 */
public class CS96_SimultDiag {
    /**
     * Performs rotation A := U.A.U^T where (i,j)-subblock of U is 
     * ( c  s)  
     * (-s  c)  while other diagonal elements of U are ones
     * That is, U[i,i] = U[j,j] = c, U[i,j] = s, U[j,i] = s 
     * and in case if a!=i & a!=j & b!=i & b!=j: U[a,b] = 1 if a==b else 0,
     * and U[i,b] == 0 for all b except for i and j,
     * and U[a,j] == 0 for all a except for i and j.
     * Note that if A is symmetrical, it remains symmetrical after transformation
     */
    private static void _2x2_rotate(Matrix A, int i, int j, double c, double s) {
        /*Matrix U = Matrix.identity(A.getRowDimension(), A.getColumnDimension());
        U.set(i, i, c);
        U.set(j, j, c);
        U.set(i, j, s);
        U.set(j, i, -s);        
        Matrix Anew = U.times(A).times(U.transpose());*/
        // forst compute B = A.U^T and save it back into a
        double[][] a = A.getArray(); // the element array itself, *not* a copy!
        double aii_new = c*c*a[i][i] + s*s*a[j][j] + c*s*(a[i][j] + a[j][i]);
        double aij_new = c*c*a[i][j] - s*s*a[j][i] + c*s*(a[j][j] - a[i][i]);
        double ajj_new = c*c*a[j][j] + s*s*a[i][i] - c*s*(a[i][j] + a[j][i]);
        double aji_new = c*c*a[j][i] - s*s*a[i][j] + c*s*(a[j][j] - a[i][i]);
        // transform the elements with only one of the indices being equal to i or j
        for (int k=0; k<a[0].length; k++) {
            if ((k!=i) && (k!=j)) {
                double aik = a[i][k];
                //double ajk = a[j][k];
                a[i][k] =  c * aik + s * a[j][k];
                a[j][k] = -s * aik + c * a[j][k];
                //
                double aki = a[k][i];
                a[k][i] =  c * aki + s * a[k][j];
                a[k][j] = -s * aki + c * a[k][j];
            }
        }
        // transform the elements with both indices equal to either i or j
        a[i][i] = aii_new;
        a[j][j] = ajj_new;
        a[i][j] = aij_new;
        a[j][i] = aji_new;
        // no transformation is required for the rest of the elements
        //
        // check:
        //System.out.println(  Anew.minus( A ).normF() );
        /*Anew.print(13, 7);
        A.print(13, 7);
        Anew.minus( A ).print(13,7);*/
    }
    //--------------------------------------------------------------------------
    
    /**
     * Performes one iteration of the Jacobi rotations, 'rotates' the matrices
     * given in @param arr array in-place and @returns the created rotation matrix
     * @param n is the dimension of the (square) matrices in @param arr array
     * @param s_max[0] (*must* be allocated before calling the method!) 
     * gets the maximum value of the |sin| of rotation angle which caller can 
     * use to control the convergence of the method
     * @param U is initial value of the rotation matrix
     */
    private static Matrix iteration(Matrix [] arr, int n, double[] s_max, Matrix U) {
        //Matrix result = Matrix.identity(n, n);
        //double G00, G01, G11, 
        double h0, h1;
        double t_on, t_off;
        s_max[0] = 0.0;
        for (int i=0; i<n; i++) {
            for (int j=i+1; j<n; j++) {
                //G00 = 0.0; 
                //G01 = 0.0; 
                //G11 = 0.0;
                t_on = 0.0;  // sum of G[0,0] - sum of G[1,1] 
                t_off = 0.0; // sum of G[0,1]
                for (int k=0; k<arr.length; k++) {
                    h0 = arr[k].get(i, i) - arr[k].get(j, j);
                    h1 = arr[k].get(i, j) + arr[k].get(j, i);
                    //G00 += h0*h0;
                    //G01 += h0*h1;
                    //G11 += h1*h1;
                    t_on += h0*h0 - h1*h1; 
                    t_off += h0*h1;
                }
                t_off *= 2.0; // sum of G[0,1]+G[1,0]
                //double theta = 0.5 * Math.atan2( t_off, t_on + Math.sqrt( t_on*t_on + t_off * t_off) );
                double theta = Math.atan2(t_off, t_on) / 4.0;
                //System.out.printf("%d, %d: theta = %.5f%n", i,j,theta);
                double c = Math.cos(theta);
                double s = Math.sin(theta);                
                if (Math.abs(s) > s_max[0])
                    s_max[0] = Math.abs(s);
                // rotate all matrices
                for (int k=0; k<arr.length; k++) 
                    _2x2_rotate(arr[k], i, j, c, s);
                // update the total rotation matrix as result = U . result
                // r[a,b] = Sum(U[a,x] . r[x,b]) = Sum(U[a,x] . r[x,b], x!=i&x!=j) + 
                //    = U[a,i]*r[i,b] + U[a,j]*r[j,b] = (r[a,b] if a!=i&a!=j else 0) +
                //    + (c if a=i, -s if a=j, else 0)*r[i,b] + (c if a=j, s if a=i, else 0 * r[j,b]
                // => r[i,b] = c*r[i,b] + s*r[j,b] , r[j,b] = -s*r[i,b] + c*r[j,b], other - unchanged
                // Nothe that in these expressions b runs over all possible index values (including i and j !)
                for (int k=0; k<n; k++) {
                    double rik = U.get(i, k);
                    U.set(i, k,  c * rik + s * U.get(j, k));
                    U.set(j, k, -s * rik + c * U.get(j, k));
                }
            }
        }
        return U;
    }
    //--------------------------------------------------------------------------
    /**
     * The main method for performing *in-place* simultanious diag. of the 
     * matrices in @param arr array
     * @returns the found transformation matrix valid for Anew = U^T.A.U-type transformations
     * of the input matrices
     * @param iters[0] can contain a negative value or be null (in this case there is
     * no limitation for the number of iterations), or contain a maximum
     * number of iterations. In the latter case, iters[0] gets the number
     * of executed iterations when the method finishes.
     * @param u0 is initialization value for the desired rotation matrix (actually 
     * not used by the algorithm itself) (can be null which means an identity matrix)
     */
    public static Matrix simult_diag(Matrix[] arr, int[] iters, Matrix u0) {
        double thresh = 1e-7;
        int n = arr[0].getRowDimension();        
        Matrix result;
        if (u0 == null)
            result = Matrix.identity(n, n);
        else {
            result = u0.transpose();
            /*
            // transform the matrices in arr[] into this basis
            for (int i=0; i< arr.length; i++)
                arr[i] = MatrixHelper.LDLtransform.TransformSymmetricMatrixToNewBasis_LDLT(arr[i], result);
*/
        }
        
        double[] s_max = new double[]{ 0.0 };
        int nIter = 0;
        do {
            iteration(arr, n, s_max, result); // note that this works in convention when Anew = U.A.U^T
            nIter++;
        } while ((s_max[0] > thresh) && ( (iters == null) || (iters[0] < 0) || (nIter < iters[0]) ));
        if (iters != null)
            iters[0] = nIter;
        //for (int k=0; k<arr.length; k++) { arr[k].print(13, 7); }
        return result.transpose() ;
    }
    //--------------------------------------------------------------------------
}
