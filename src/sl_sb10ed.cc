/*

Copyright (C) 2009, 2010, 2011   Lukas F. Reichlin

This file is part of LTI Syncope.

LTI Syncope is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LTI Syncope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

H2 optimal state controller for a discrete-time system.
Uses SLICOT SB10ED by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.5

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10ed, SB10ED)
                 (int& N, int& M, int& NP,
                  int& NCON, int& NMEAS,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* AK, int& LDAK,
                  double* BK, int& LDBK,
                  double* CK, int& LDCK,
                  double* DK, int& LDDK,
                  double* RCOND,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  bool* BWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_sb10ed__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10ed__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10ED Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 6)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        
        int ncon = args(4).int_value ();
        int nmeas = args(5).int_value ();
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int np = c.rows ();     // np: number of outputs
        
        int lda = max (1, a.rows ());
        int ldb = max (1, b.rows ());
        int ldc = max (1, c.rows ());
        int ldd = max (1, d.rows ());
        
        int ldak = max (1, n);
        int ldbk = max (1, n);
        int ldck = max (1, ncon);
        int lddk = max (1, ncon);
        
        double tol = 0;
        
        // arguments out      
        Matrix ak (ldak, n);
        Matrix bk (ldbk, nmeas);
        Matrix ck (ldck, n);
        Matrix dk (lddk, nmeas);
        ColumnVector rcond (7);
        
        // workspace       
        int m2 = ncon;
        int m1 = m - m2;
        int np1 = np - nmeas;
        int np2 = nmeas;
        
        int q = max (m1, m2, np1, np2);
        int ldwork = 2*q*(3*q+2*n)+max(1,(n+q)*(n+q+6),q*(q+max(n,q,5)+1),
                      2*n*n+max(1,14*n*n+6*n+max(14*n+23,16*n),
                                q*(n+q+max(q,3))));
        int liwork = max (2*m2, 2*n, n*n, np2);
        
        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (bool, bwork, 2*n);
        
        // error indicator
        int info;


        // SLICOT routine SB10ED
        F77_XFCN (sb10ed, SB10ED,
                 (n, m, np,
                  ncon, nmeas,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  ak.fortran_vec (), ldak,
                  bk.fortran_vec (), ldbk,
                  ck.fortran_vec (), ldck,
                  dk.fortran_vec (), lddk,
                  rcond.fortran_vec (),
                  tol,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("h2syn: __sl_sb10ed__: exception in SLICOT subroutine SB10ED");

        static const char* err_msg[] = {
            "0: OK",
            "1: the matrix [A-exp(j*Theta)*I, B2; C1, D12] "
                "had not full column rank in respect to the tolerance EPS",
            "2: the matrix [A-exp(j*Theta)*I, B1; C2, D21] "
                "had not full row rank in respect to the tolerance EPS",
            "3: the matrix D12 had not full column rank in "
                "respect to the tolerance TOL",
            "4: the matrix D21 had not full row rank in respect "
                "to the tolerance TOL",
            "5: the singular value decomposition (SVD) algorithm "
                "did not converge (when computing the SVD of one of the matrices "
                "[A-I, B2; C1, D12], [A-I, B1; C2, D21], D12 or D21)",
            "6: the X-Riccati equation was not solved successfully",
            "7: the matrix Im2 + B2'*X2*B2 is not positive "
                "definite, or it is numerically singular (with "
                "respect to the tolerance TOL)",
            "8: the Y-Riccati equation was not solved successfully",
            "9: the matrix Ip2 + C2*Y2*C2' is not positive "
                "definite, or it is numerically singular (with "
                "respect to the tolerance TOL)",
            "10: the matrix Im2 + DKHAT*D22 is singular, or its "
                "estimated condition number is larger than or equal "
                "to 1/TOL"};

        error_msg ("h2syn", info, 10, err_msg);


        // return values
        retval(0) = ak;
        retval(1) = bk;
        retval(2) = ck;
        retval(3) = dk;
        retval(4) = rcond;
    }
    
    return retval;
}
