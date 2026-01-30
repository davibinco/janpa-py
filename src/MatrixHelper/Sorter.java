package MatrixHelper;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
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
 * (c) Tymofii Nikolaienko, 2016
 *
 */

/**
 * A collection of (static) sorting methods
 */
public class Sorter {
    
    //--------------------------------------------------------------------------
    // TODO: add boolean ascendingly parameter
    public static double[] merge_ascending(double[] a, double[] b) {
        double[] result = new double[a.length + b.length];
        int i = 0;
        int j = 0;
        int k = 0;
        while (k < result.length) {
            while ((i < a.length) && ((j == b.length) ||  (a[i] < b[j]))) {
                result[k++] = a[i];
                i++;
            }
            while ((j < b.length) && ((i == a.length) || (b[j] < a[i]))) {
                result[k++] = b[j];
                j++;
            }
        }
        return result;
    }
    //--------------------------------------------------------------------------
    private static void _test() {        
        
/*
        for (double  el: x) {
            System.out.printf("%.0f\t",el);
        }
        System.out.printf("");
        */
    }
    //--------------------------------------------------------------------------    
    /**
     * Sorts an array of doubles and @returns a remap such that
     * array[remap[i]] is sorted ascendingly (or descendingly, if !ascendingly)
     */
    // TODO: merge sort
    public static int[] array_sort(double[] array, boolean ascendingly) {
        
        int n = array.length;
        boolean done, found;
        int[] result = new int[n];
        for (int i=0; i<n; i++) result[i]=i;
        found = true;
        
        while (found) {
            found = false;
            for (int i=0; i<n-1; i++)
                if (( array[result[i]] < array[result[i+1]] ) ^ ascendingly) {
                    found = true;
                    // swap pnao-th and (pnao+1)-th functions
                    int tmp = result[i];
                    result[i] = result[i+1];
                    result[i+1] = tmp;
            }
        }
        return result;
    }
    //--------------------------------------------------------------------------
    public static void array_remap(double[] array, int[] remap) {
        double[] newArr = new double[array.length];
        for (int i=0; i<array.length; i++)
            newArr[i] = array[remap[i]];        
        for (int i=0; i<array.length; i++)
            array[i] = newArr[i];        
    }
    //--------------------------------------------------------------------------
    
    /** uses array_sort() to re-order the elements of @param array in-place
     Rev. 18.Feb.2018 **/
    public static int[] array_sort_inPlace(double[] array, boolean ascendingly) {
        //double[] newArr = new double[array.length];
        int[] remap = array_sort(array, ascendingly);
        array_remap(array, remap);
        /*for (int i=0; i<array.length; i++)
            newArr[i] = array[remap[i]];        
        for (int i=0; i<array.length; i++)
            array[i] = newArr[i];        */
        return remap;
    }
    //--------------------------------------------------------------------------
    
    
    // TODO: call array_sort
    
    /** @returns a remap array such that Eigenvals[ result[1] ] > Eigenvals[ result[2] ] > ...
     * 
     */
    public static int[] sort_eigenvalues(EigenvalueDecomposition eig, int NDim) {
        if (NDim == -1) {
            NDim = eig.getD().getColumnDimension();
        }
        int[] remap = new int[NDim]; // an array to return
        double[] eigenvals = new double[NDim]; // for readability of code only
        // initialize these arrays
        for (int k=0; k<NDim; k++) {
            remap[k] = k;
            eigenvals[k] = eig.getD().get(k, k);
        }
        // use bubble sort; we do not work with too long arrays anyway!
        boolean found = true;
        while (found) {
            found = false; // suppose that everything has already been sorted
            for (int k=0; k<(NDim-1); k++) {
                if (eigenvals[remap[k]] < eigenvals[remap[k+1]]) {
                    found = true;
                    // swap k-th and k+1-th elements
                    int tmp = remap[k];
                    remap[k] = remap[k+1];
                    remap[k+1] = tmp;
                }
            }
        }
        return remap;
    }
    //--------------------------------------------------------------------------
    /**
     * Sorts eigenvectors passed in @param eig according to their eigenvalues
     * @returns [ vals = double[] of eigenvalues sorted descendingly, 
     *            vecs = Matrix of eigenvectors corresponding to these eigenvalues (a-la .getV()) ]
     * Note: vecs[i-th component, of k-th eigenvector] (i.e., column-vectors are used)
     * 
     * Rev.: Aug.26, 2016
     * Rev.: Jan.14, 2018 -- added feature for nDim == -1
     * 
     */
    public static Object[] sorted_eigenData(EigenvalueDecomposition eig, int nDim) {
        if (nDim == -1) {
            nDim = eig.getD().getColumnDimension();
        }        
        int[] remap = sort_eigenvalues(eig, nDim);
        Object[] result = new Object[]{ new double[nDim], new Matrix(nDim, nDim) };
        Matrix V = eig.getV();
        for(int i=0; i<nDim; i++) {
            ((double[])result[0])[i] = eig.getRealEigenvalues()[ remap[i] ];
            ((Matrix)result[1]).setMatrix(0, nDim-1, remap[i], remap[i], 
                    V.getMatrix(0, nDim-1, i, i)
                    );
        }
        return result;
    }
    //--------------------------------------------------------------------------
        
    
}
