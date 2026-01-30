package JGints;

/**
 *
 * Represents a polynomial of two variables (named 'rho' and 'z') and performs
 * some manipulations with these polynomials (e.g., multiplication over
 * some pre-defined terms, saving/restoring onto/from the stack, etc.)
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
public class Polynom_rho_z {
    

    //--------------------------------------------------------------------------

    public double[][] cf_matrix = null;
    private double[][] tmp_matrix = null;
    
    /*private int[] first_nonzero_index = null;
    private int[] last_nonzero_index = null;*/
    public int last_active_rho_power = 0;
    public int last_active_z_power = 0;
    
    private Polynom_rho_z stackPointer = null;
    //--------------------------------------------------------------------------
    /**
     * Copies from @param source all the fields which can be properly 
     * copied by shallow (rather than the deep)-copying
     */
    private void fieldsCopy(Polynom_rho_z source, boolean doDeepCopy) {
        this.last_active_z_power = source.last_active_z_power;
        this.last_active_rho_power = source.last_active_rho_power;
        // note: we ignore stackPointer field !
        if (doDeepCopy) {
            cf_matrix = new double[source.cf_matrix.length][];
            for(int i=0; i<source.cf_matrix.length; i++) // todo: up to last_acrive_z_power only
                cf_matrix[i] = source.cf_matrix[i].clone();
        } else {
            cf_matrix = source.cf_matrix;
        }        
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a 'record' (rather than a fully functional object)
     * with a deep-copies of all 'meaningful' (without, e.g., tmp_matrix)
     * fields from a given @param parent 
     * !! Invoked ONLY by push/pop functionality !!
     */
    private Polynom_rho_z(Polynom_rho_z parent) {
        // note: we do not close tmp_matrix which is supposed to be a valid field
        // only within the object which 'hosts' the stack
        fieldsCopy(parent, true);
    }
    //--------------------------------------------------------------------------
    /**
     * Creates the polynomial initially equal to 1.0
     */
    public Polynom_rho_z(int max_z_pwr, int max_rho_pwr) {
        cf_matrix = new double[max_z_pwr+1][max_rho_pwr+1]; // +1 for 0-th power (i.e., const)
        cf_matrix[0][0] = 1.0;

        tmp_matrix = new double[max_z_pwr+1][max_rho_pwr+1];
        last_active_rho_power = 0;
        last_active_z_power = 0;
    }
    //--------------------------------------------------------------------------
    /**
     * Saves the current state of the object onto the stack
     */
    public void pushToStack() {
        // save the current state as a new object
        Polynom_rho_z backup = new Polynom_rho_z(this);
        // and put it on the top of the stack
        backup.stackPointer = this.stackPointer;        
        this.stackPointer = backup;
    }
    
    /**
     * Restores the last saved state
     */
    public void restoreFromStack(boolean deleteSourceFromStack) {
        // save the current state as a new object
        Polynom_rho_z backup = stackPointer;
        if (backup == null)
            return; // nothing to do!
        // else:        

        // 1) copy fields from backup properly
        fieldsCopy(backup, ! deleteSourceFromStack ); // we don't need a deep copy if the element on the stack is to be deleted
        
        // 2) reorganize stack if necessary
        if (deleteSourceFromStack) {
            stackPointer = backup.stackPointer; // remove from stack
        }
    }
    //--------------------------------------------------------------------------

    /** 
     * Swap 'pointers' cf_matrix <-> tmp_matrix, so that make the current
     * cf_matrix available as the tmp_matrix storage during the next call
     */
    private void _swap_old_new() {
        double[][] old_ptr = cf_matrix;
        cf_matrix = tmp_matrix;
        tmp_matrix = old_ptr;            
    }
    //--------------------------------------------------------------------------

    /**
     * multiplies the current polynom by (cfZ * (z-zA)^2 + cfRho * rho^2)
     * and converts the z argument into the form (z-z0)^2
     * Note that 
     * (z-zA)^2 = (z-z0-(zA-z0))^2 = (z-z0)^2 - 2*(z-z0)*(zA-z0) + (zA-z0)^2,
     * i.e., the final multiplication is by
     * (cfZ * (z-z0)^2 - (z-z0) * cfZ * 2*(zA-z0) + cfRho * rho^2 + cfZ*(zA-z0)^2)
     * 
     * hence, for reccurent impl.:
     * cf_matrix[i,j] := cfZ * cf_matrix[i-2,j] -(zA-z0)*2*cfZ * cf_matrix[i-1,j] + cfRho * cf_matrix[i,j-2] + cf_matrix[i,j] * cfZ*(zA-z0)^2
     * 
     */
    public void mul_sq(double cfZ, double cfRho, double zA_minus_z0) {
        //assert cf_matrix.length > 0; // TODO
        //TODO: eliminate the use of tmp_matrix by careful code refactoring
        
        //int nRhoSz = last_active_rho_power;//cf_matrix[0].length;
        // TODO: boundary check!
        //last_active_rho_power = cf_matrix[0].length-3;

        // mul. by: cfZ * (z-z0)^2
        for(int i=0; i<=last_active_z_power; ++i) {
            for (int j=0; j<=last_active_rho_power; ++j) {
                tmp_matrix[i+2][j] = cfZ * cf_matrix[i][j];
            }
        }
        // last_acrive_z_power += 2; is better done at the end, after adding a linear term;
        // this avoids unnecessary 'last_acrive_z_power-2' in the following 'for i' loops 
        
        // note that cf_matrix[0] and cf_matrix[1] were not modified in the prev. loop => zero them
        for (int j=0; j<=last_active_rho_power+2; ++j) { // +2 is important here since we're going to increase last_active_rho_power by 2 soon !
            tmp_matrix[0][j] = 0.0;
            tmp_matrix[1][j] = 0.0;
        }
        // fill with zeros also another block of tmp_matrix to which we're going to write to soon
        for (int i=0; i<=last_active_z_power+2; ++i) {
            tmp_matrix[i][last_active_rho_power+1] = 0.0;
            tmp_matrix[i][last_active_rho_power+2] = 0.0;
        }        
        
        // mul. by: -(z-z0) * cfZ * 2*(zA-z0)
        double tmp = cfZ * 2 * zA_minus_z0;
        for(int i=0; i<=last_active_z_power; ++i) {
            for (int j=0; j<=last_active_rho_power; ++j) {
                tmp_matrix[i+1][j] -= tmp * cf_matrix[i][j];
            }
        }
        
        // mul. by: cfRho * rho^2
        for(int i=0; i<=last_active_z_power; ++i) {
            // we're going to increase last_active_rho_power by 2 in a moment;
            // but now, although we expect that no access to those 'higher' (yet
            // unused!) powers has been occured yet, it seems better to zero
            // them to be on a safe side
            for (int j=0; j<=last_active_rho_power; ++j) {
                tmp_matrix[i][j+2] += cfRho * cf_matrix[i][j];
            }
        }                       
        // last_active_rho_power += 2; is better done at the end, after adding a linear term;
        // this avoids unnecessary 'last_active_rho_power-2' in the following 'for j' loops
        
        // mul. by: cfZ*(zA-z0)^2
        tmp =  cfZ * zA_minus_z0 * zA_minus_z0;
        for(int i=0; i<=last_active_z_power; ++i) {
            for (int j=0; j<=last_active_rho_power; ++j) {
                tmp_matrix[i][j] += tmp * cf_matrix[i][j];
            }
        }
        
        last_active_rho_power += 2;
        last_active_z_power += 2;

        _swap_old_new();
    }
    //--------------------------------------------------------------------------

    /**
     * multiplies the current polynom by (z-zA)
     * and converts the z argument into the form (z-z0)
     * Note that 
     * (z-zA) = (z-z0) - (zA-z0)
     * i.e., the final multiplication is by
     * ( (z-z0) - (zA-z0) )
     */
    public void mul_lin(double zA_minus_z0) {
        //assert cf_matrix.length > 0; // TODO

        // a faster (although more tricky) version which modifies cf_matrix in-place:
        // cf_matrix := Z * cf_matrix - Const * cf_matrix
        // in this case:
        // cf_matrix[X][j] := cf_matrix[X-1][j]                            , for X==last_active_z_power+1, and        
        // cf_matrix[i][j] := cf_matrix[i-1][j] - Const * cf_matrix[i][j]  , for 1 <= i <= last_active_z_power & i-- (!), and
        // cf_matrix[0  ][j] :=                 - Const * cf_matrix[0  ][j], for i == 0

        // mul. by: (z-z0), the src. terms ~ to z^(last_active_z_power)
        System.arraycopy(cf_matrix[last_active_z_power], 0, 
                         cf_matrix[last_active_z_power+1], 0,
                         last_active_rho_power+1); // +1 in order not to miss the [last_active_rho_power]-th element
        
        // mul. by: -(zA-z0) AND mul. by (z-z0) -- all 'intermediate' terms
        for(int i=last_active_z_power; i>=1; --i) {
            for (int j=0; j<=last_active_rho_power; ++j) {
                //cf_matrix[i][j] = cf_matrix[i-1][j] - zA_minus_z0 * cf_matrix[i][j];
                cf_matrix[i][j] *= -zA_minus_z0;  // seems to work faster than the prev. line !
                cf_matrix[i][j] += cf_matrix[i-1][j];
            }
        }
        
        // mul. by: -(zA-z0) -- the 0-th terms
        for (int j=0; j<=last_active_rho_power; ++j) {
            cf_matrix[0][j] *=  -zA_minus_z0;
        }
        

        last_active_z_power += 1;

        // note that we haven't used the temporary storage => no need to call  _swap_old_new();
        
        /*
        // mul. by: (z-z0)
        
        for(int i=0; i<=last_active_z_power; ++i) {
            // perform
            System.arraycopy(cf_matrix[i], 0, tmp_matrix[i+1], 0, last_active_rho_power+1); 
            // by interchanging the matrix rows:
            //double[] tmpPtr = tmp_matrix[i+1];
            //tmp_matrix[i+1] = cf_matrix[i];
            //cf_matrix[i] = tmpPtr;
        }
        
        // note that cf_matrix[0] was not modified => zero it
        for (int j=0; j<=last_active_rho_power; ++j) {
            tmp_matrix[0][j] = 0.0;
        }
        // mul. by: -(zA-z0)
        for(int i=0; i<=last_active_z_power; ++i) {
            for (int j=0; j<=last_active_rho_power; ++j) {
                tmp_matrix[i][j] -= zA_minus_z0 * cf_matrix[i][j];
            }
        }            

        last_active_z_power += 1;
        
        _swap_old_new();
        */
    }
    //--------------------------------------------------------------------------

    /**
     * Multiplies the current polynom by associated Legendre 
     * polynomial (without (1-z^2)^(m/2) and without normalization const.
     * Note that the sign of m is also ignored !
     */
    public void mul_AssocLeg(int L, int m, double zA_minus_z0) {
        double[] squared_roots = SphericalHarmonics.ALegRootsSquared[L][Math.abs(m)];
        int k = squared_roots.length; // (L - Math.abs(m)) / 2;
        for(int i=0; i<k; i++) {
            mul_sq(1 - squared_roots[i], - squared_roots[i], zA_minus_z0);
        }   
        if ( (L - Math.abs(m)) % 2 == 1 ) {
            // we must also multiply by (z-zA)
            mul_lin(zA_minus_z0);
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Multiplies the current polynom by ( (z-zA)^2 + rho^2 )^r2_pwr
     */        
    public void mul_r2(int r2_pwr, double zA_minus_z0) {
        for(int i=0; i<r2_pwr; i++) {
            mul_sq(1.0, 1.0, zA_minus_z0);
        }
    }
    //--------------------------------------------------------------------------
    public void printout() {
        for (int i=0; i<cf_matrix.length; i++) {
            boolean somethingPrinted = false;
            for (int j=0; j<cf_matrix[i].length; j++) {
                if (Math.abs(cf_matrix[i][j]) > 1e-10) {
                    System.out.printf(" %+.3f * z^%d * rho^%d ", cf_matrix[i][j], i, j);
                    somethingPrinted = true;
                }
            }            
            if (somethingPrinted)
                System.out.println();
        }
        System.out.println("---");
    }
    //--------------------------------------------------------------------------

    
}
