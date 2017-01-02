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

Convert descriptor state-space system into regular state-space form.
Uses SLICOT SB10JD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10jd, SB10JD)
                 (F77_INT& N, F77_INT& M, F77_INT& NP,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  double* E, F77_INT& LDE,
                  F77_INT& NSYS,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb10jd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10jd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10JD Release 5.0\n\
No argument checking.\n\
For internal use only.")
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
        Matrix e = args(4).matrix_value ();
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT np = TO_F77_INT (c.rows ());     // np: number of outputs
        
        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, np);
        F77_INT ldd = max (1, np);
        F77_INT lde = max (1, n);

        // arguments out
        F77_INT nsys;
        
        // workspace
        F77_INT ldwork = max (1, 2*n*n + 2*n + n*max (5, n + m + np));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        F77_INT info;


        // SLICOT routine SB10JD
        F77_XFCN (sb10jd, SB10JD,
                 (n, m, np,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  e.fortran_vec (), lde,
                  nsys,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_sb10jd__: exception in SLICOT subroutine SB10JD");

        if (info != 0)
            error ("__sl_sb10jd__: SB10JD returned info = %d", info);

        // resize
        a.resize (nsys, nsys);
        b.resize (nsys, m);
        c.resize (np, nsys);

        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = d;
    }
    
    return retval;
}
