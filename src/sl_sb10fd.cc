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

H-infinity (sub)optimal state controller for a continuous-time system.
Uses SLICOT SB10FD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: December 2009
Version: 0.6

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10fd, SB10FD)
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
                  double* RCOND,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_LOGICAL* BWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb10fd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10fd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_sb10fd__ (@dots{})\n"
   "Wrapper for SLICOT function SB10FD.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
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
        
        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan () || d.any_element_is_inf_or_nan ())
          error ("__sl_sb10fd__: inputs must not contain NaN or Inf\n");

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
        
        double tol = 0;
        
        // arguments out
        Matrix ak (ldak, n);
        Matrix bk (ldbk, nmeas);
        Matrix ck (ldck, n);
        Matrix dk (lddk, nmeas);
        ColumnVector rcond (4);
        
        // workspace
        
        F77_INT m2 = ncon;
        F77_INT m1 = m - m2;
        F77_INT np1 = np - nmeas;
        F77_INT np2 = nmeas;
        
        F77_INT liwork = max (2*max (n, m-ncon, np-nmeas, ncon), n*n);
        F77_INT q = max (m1, m2, np1, np2);
        F77_INT ldwork = 2*q*(3*q+2*n) + max (1, (n+q)*(n+q+6), q*(q + max (n, q, 5) + 1),
                                          2*n*(n+2*q) + max (1, 4*q*q +
                                          max (2*q, 3*n*n + max (2*n*q, 10*n*n+12*n+5)),
                                          q*(3*n + 3*q + max (2*n, 4*q + max (n, q)))));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (F77_LOGICAL, bwork, 2*n);
        
        // error indicator
        F77_INT info;


        // SLICOT routine SB10FD
        F77_XFCN (sb10fd, SB10FD,
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
                  rcond.fortran_vec (),
                  tol,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("hinfsyn: __sl_sb10fd__: exception in SLICOT subroutine SB10FD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the matrix [A-j*omega*I, B2; C1, D12] had "
                "not full column rank in respect to the tolerance EPS",
            "2: the matrix [A-j*omega*I, B1; C2, D21] "
                "had not full row rank in respect to the tolerance EPS",
            "3: the matrix D12 had not full column rank in "
                "respect to the tolerance TOL",
            "4: the matrix D21 had not full row rank in respect "
                "to the tolerance TOL",
            "5: the singular value decomposition (SVD) algorithm "
                "did not converge (when computing the SVD of one of the matrices "
                "[A, B2; C1, D12], [A, B1; C2, D21], D12 or D21)",
            "6: the controller is not admissible (too small value "
                "of gamma)",
            "7: the X-Riccati equation was not solved "
                "successfully (the controller is not admissible or "
                "there are numerical difficulties)",
            "8: the Y-Riccati equation was not solved "
                "successfully (the controller is not admissible or "
                "there are numerical difficulties)",
            "9: the determinant of Im2 + Tu*D11HAT*Ty*D22 is "
                "zero [3]"};

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
