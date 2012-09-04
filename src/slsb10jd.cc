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

Convert descriptor state-space system into regular state-space form.
Uses SLICOT SB10JD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2011
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10jd, SB10JD)
                 (int& N, int& M, int& NP,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* E, int& LDE,
                  int& NSYS,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_sb10jd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10jd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10JD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
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
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int np = c.rows ();     // np: number of outputs
        
        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, np);
        int ldd = max (1, np);
        int lde = max (1, n);

        // arguments out
        int nsys;
        
        // workspace
        int ldwork = max (1, 2*n*n + 2*n + n*max (5, n + m + np));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        int info;


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
