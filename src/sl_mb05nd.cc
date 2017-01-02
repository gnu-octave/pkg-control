/*

Copyright (C) 2014-2015   Thomas Vasileiou

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

Matrix exponential and integral for a real matrix.
Uses SLICOT MB05ND by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Thomas Vasileiou <thomas-v@wildmail.com>
Created: March 2014
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (mb05nd, MB05ND)
                 (F77_INT& N, double& DELTA,
                  double* A, F77_INT& LDA,
                  double* EX, F77_INT& LDEX,
                  double* EXINT, F77_INT& LDEXINT,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_mb05nd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_mb05nd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot MB05ND Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 3)
    {
        print_usage ();
    }
    else
    {
        // arguments in        
        Matrix a = args(0).matrix_value ();
        double delta = args(1).double_value ();
        double tol = args(2).double_value (); 
        
        F77_INT n = TO_F77_INT (a.rows ());
        F77_INT lda = max (1, n);
        F77_INT ldex = max (1, n);
        F77_INT ldexin = max (1, n);

        // arguments out
        Matrix ex (ldex, n);
        Matrix exint (ldexin, n);
        
        // workspace
        F77_INT ldwork = max (1, 2*n*n);                 // optimum performance
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT info = 0;


        // SLICOT routine MB05ND
        F77_XFCN (mb05nd, MB05ND,
                 (n, delta,
                  a.fortran_vec (), lda,
                  ex.fortran_vec (), ldex,
                  exint.fortran_vec (), ldexin,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_mb05nd__: exception in SLICOT subroutine MB05ND");
        
        if (info > 0)
        {
          if (info == n+1)
            info = 2;
          else  
            info = 1;
        }
        
            
        static const char* err_msg[] = {
            "0: OK",
            "1: some element of the denominator of the Pade "
                "approximation is zero, so the denominator "
                "is exactly singular.",
            "2: DELTA = (delta * frobenius norm of matrix A) is "
                "probably too large to permit meaningful computation.  "
                "That is, DELTA > SQRT(BIG), where BIG is a "
                "representable number near the overflow threshold of "
                "the machine."}; 
 
        error_msg ("__sl_mb05nd__", info, 2, err_msg);

        // return values
        retval(0) = ex;
        retval(1) = exint;
    }
    
    return retval;
}
