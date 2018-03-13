/*

Copyright (C) 2009-2016   Lukas F. Reichlin

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

H-infinity (sub)optimal state controller for a discrete-time system.
Uses SLICOT SB10DD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: December 2009
Version: 0.6

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10dd, SB10DD)
                 (F77_INT& N, F77_INT& M, F77_INT& NP,
                  F77_INT& NCON, F77_INT& NMEAS,
                  double& GAMMA,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  double* AK, F77_INT& LDAK,
                  double* BK, F77_INT& LDBK,
                  double* CK, F77_INT& LDCK,
                  double* DK, F77_INT& LDDK,
                  double* X, F77_INT& LDX,
                  double* Z, F77_INT& LDZ,
                  double* RCOND,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_LOGICAL* BWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb10dd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10dd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10DD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 7)
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
        
        F77_INT ncon = args(4).int_value ();
        F77_INT nmeas = args(5).int_value ();
        double gamma = args(6).double_value ();
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT np = TO_F77_INT (c.rows ());     // np: number of outputs
        
        F77_INT lda = max (1, TO_F77_INT (a.rows ()));
        F77_INT ldb = max (1, TO_F77_INT (b.rows ()));
        F77_INT ldc = max (1, TO_F77_INT (c.rows ()));
        F77_INT ldd = max (1, TO_F77_INT (d.rows ()));
        
        F77_INT ldak = max (1, n);
        F77_INT ldbk = max (1, n);
        F77_INT ldck = max (1, ncon);
        F77_INT lddk = max (1, ncon);
        
        F77_INT ldx = max (1, n);
        F77_INT ldz = max (1, n);
        
        double tol = 0;
        
        // arguments out
        Matrix ak (ldak, n);
        Matrix bk (ldbk, nmeas);
        Matrix ck (ldck, n);
        Matrix dk (lddk, nmeas);
        Matrix x (ldx, n);
        Matrix z (ldz, n);
        ColumnVector rcond (8);
        
        // workspace
        F77_INT m2 = ncon;
        F77_INT m1 = m - m2;
        F77_INT np1 = np - nmeas;
        F77_INT np2 = nmeas;
        
        F77_INT liwork = max (2*max (m2, n), m, m2+np2, n*n);
        F77_INT q = max (m1, m2, np1, np2);
        F77_INT ldwork = max ((n+q)*(n+q+6), 13*n*n + m*m + 2*q*q + n*(m+q) +
                     max (m*(m+7*n), 2*q*(8*n+m+2*q)) + 6*n +
                     max (14*n+23, 16*n, 2*n + max (m, 2*q), 3*max (m, 2*q)));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (F77_LOGICAL, bwork, 2*n);
        
        // error indicator
        F77_INT info;


        // SLICOT routine SB10DD
        F77_XFCN (sb10dd, SB10DD,
                 (n, m, np,
                  ncon, nmeas,
                  gamma,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  ak.fortran_vec (), ldak,
                  bk.fortran_vec (), ldbk,
                  ck.fortran_vec (), ldck,
                  dk.fortran_vec (), lddk,
                  x.fortran_vec (), ldx,
                  z.fortran_vec (), ldz,
                  rcond.fortran_vec (),
                  tol,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("hinfsyn: __sl_sb10dd__: exception in SLICOT subroutine SB10DD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the matrix [A-exp(j*Theta)*I, B2; C1, D12] "
                "had not full column rank",
            "2: the matrix | A-exp(j*Theta)*I, B1; C2, D21] "
                "had not full row rank",
            "3: the matrix D12 had not full column rank",
            "4: the matrix D21 had not full row rank",
            "5: the controller is not admissible "
                "(too small value of gamma)",
            "6: the X-Riccati equation was not solved "
                "successfully (the controller is not admissible or "
                "there are numerical difficulties)",
            "7: the Z-Riccati equation was not solved "
                "successfully (the controller is not admissible or "
                "there are numerical difficulties)",
            "8: the matrix Im2 + DKHAT*D22 is singular",
            "9: the singular value decomposition (SVD) algorithm "
                "did not converge (when computing the SVD of one of "
                "the matrices [A, B2; C1, D12], [A, B1; C2, D21], D12 or D21)"};

        error_msg ("hinfsyn", info, 9, err_msg);

        
        // return values
        retval(0) = ak;
        retval(1) = bk;
        retval(2) = ck;
        retval(3) = dk;
        retval(4) = rcond;
    }
    
    return retval;
}
