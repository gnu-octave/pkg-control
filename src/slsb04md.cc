/*

Copyright (C) 2010   Lukas F. Reichlin

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

Solution of continuous-time Sylvester equations.
Uses SLICOT SB04MD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: January 2010
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb04md, SB04MD)
                 (int& N, int& M,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* Z, int& LDZ,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_sb04md__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb04md__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB04MD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 3)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        
        int n = a.rows ();
        int m = b.rows ();
        
        int lda = max (1, n);
        int ldb = max (1, m);
        int ldc = max (1, n);
        int ldz = max (1, m);
        
        // arguments out
        Matrix z (ldz, m);
        
        // workspace
        int ldwork = max (1, 2*n*n + 8*n, 5*m, n + m);
        
        OCTAVE_LOCAL_BUFFER (int, iwork, 4*n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        int info;
        

        // SLICOT routine SB04MD
        F77_XFCN (sb04md, SB04MD,
                 (n, m,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  z.fortran_vec (), ldz,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("lyap: __sl_sb04md__: exception in SLICOT subroutine SB04MD");

        if (info != 0)
            error ("lyap: __sl_sb04md__: SB04MD returned info = %d", info);
        
        // return values
        retval(0) = c;
    }
    
    return retval;
}
