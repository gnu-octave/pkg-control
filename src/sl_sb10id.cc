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

Positive feedback controller for a continuous-time system.
Uses SLICOT SB10ID by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: July 2011
Version: 0.4

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10id, SB10ID)
                 (F77_INT& N, F77_INT& M, F77_INT& NP,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  double& FACTOR, F77_INT& NK,
                  double* AK, F77_INT& LDAK,
                  double* BK, F77_INT& LDBK,
                  double* CK, F77_INT& LDCK,
                  double* DK, F77_INT& LDDK,
                  double* RCOND,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_LOGICAL* BWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb10id__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10id__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_sb10id__ (@dots{})\n"
   "Wrapper for SLICOT function SB10ID.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 5)
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
          error ("__sl_sb10id__: inputs must not contain NaN or Inf\n");

        double factor = args(4).double_value ();
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT np = TO_F77_INT (c.rows ());     // np: number of outputs
        
        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, np);
        F77_INT ldd = max (1, np);
        
        F77_INT ldak = max (1, n);
        F77_INT ldbk = max (1, n);
        F77_INT ldck = max (1, m);
        F77_INT lddk = max (1, m);
        
        // arguments out
        F77_INT nk;
        Matrix ak (ldak, n);
        Matrix bk (ldbk, np);
        Matrix ck (ldck, n);
        Matrix dk (lddk, np);
        ColumnVector rcond (2);
        
        // workspace
        F77_INT liwork = max (2*n, n*n, m, np);
        F77_INT ldwork = 10*n*n + m*m + np*np + 2*m*n + 2*n*np + 4*n +
                     5 + max (1, 4*n*n + 8*n);

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (F77_LOGICAL, bwork, 2*n);
        
        // error indicator
        F77_INT info;


        // SLICOT routine SB10ID
        F77_XFCN (sb10id, SB10ID,
                 (n, m, np,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  factor, nk,
                  ak.fortran_vec (), ldak,
                  bk.fortran_vec (), ldbk,
                  ck.fortran_vec (), ldck,
                  dk.fortran_vec (), lddk,
                  rcond.fortran_vec (),
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("ncfsyn: __sl_sb10id__: exception in SLICOT subroutine SB10ID");

        static const char* err_msg[] = {
            "0: OK",
            "1: the X-Riccati equation is not solved successfully",
            "2: the Z-Riccati equation is not solved successfully",
            "3: the iteration to compute eigenvalues or singular "
                "values failed to converge",
            "4: the matrix Ip - D*Dk is singular",
            "5: the matrix Im - Dk*D is singular",
            "6: the closed-loop system is unstable"};

        error_msg ("ncfsyn", info, 6, err_msg);


        // resizing
        ak.resize (nk, nk);
        bk.resize (nk, np);
        ck.resize (m, nk);
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
