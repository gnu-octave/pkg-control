/*

Copyright (C) 2011   Lukas F. Reichlin

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

Positive feedback controller for a discrete-time system (D != 0).
Uses SLICOT SB10ZD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: August 2011
Version: 0.3

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10zd, SB10ZD)
                 (int& N, int& M, int& NP,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double& FACTOR,
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

// PKG_ADD: autoload ("__sl_sb10zd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10zd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10ZD Release 5.0\n\
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
        
        double factor = args(4).double_value ();
        double tol = args(5).double_value ();
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int np = c.rows ();     // np: number of outputs
        
        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, np);
        int ldd = max (1, np);
        
        int ldak = max (1, n);
        int ldbk = max (1, n);
        int ldck = max (1, m);
        int lddk = max (1, m);
        
        // arguments out
        Matrix ak (ldak, n);
        Matrix bk (ldbk, np);
        Matrix ck (ldck, n);
        Matrix dk (lddk, np);
        ColumnVector rcond (6);
        
        // workspace
        int liwork = 2 * max (n, m+np);
        int ldwork = 16*n*n + 5*m*m + 7*np*np + 6*m*n + 7*m*np +
                     7*n*np + 6*n + 2*(m + np) +
                     max (14*n+23, 16*n, 2*m-1, 2*np-1);

        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (bool, bwork, 2*n);
        
        // error indicator
        int info;


        // SLICOT routine SB10ZD
        F77_XFCN (sb10zd, SB10ZD,
                 (n, m, np,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  factor,
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
            error ("ncfsyn: __sl_sb10zd__: exception in SLICOT subroutine SB10ZD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the P-Riccati equation is not solved successfully",
            "2: the Q-Riccati equation is not solved successfully",
            "3: the iteration to compute eigenvalues or singular "
                "values failed to converge",
            "4: the matrix (gamma^2-1)*In - P*Q is singular",
            "5: the matrix Rx + Bx'*X*Bx is singular",
            "6: the matrix Ip + D*Dk is singular",
            "7: the matrix Im + Dk*D is singular",
            "8: the matrix Ip - D*Dk is singular",
            "9: the matrix Im - Dk*D is singular",
            "10: the closed-loop system is unstable"};

        error_msg ("ncfsyn", info, 10, err_msg);


        // resizing
        ak.resize (n, n);
        bk.resize (n, np);
        ck.resize (m, n);
        dk.resize (m, np);
        
        // return values
        retval(0) = ak;
        retval(1) = bk;
        retval(2) = ck;
        retval(3) = dk;
        retval(4) = rcond;
    }
    
    return retval;
}
