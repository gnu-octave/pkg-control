/*

Copyright (C) 2010, 2011   Lukas F. Reichlin

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

Minimal realization of state-space models.
Uses SLICOT TB01PD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.3

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb01pd, TB01PD)
                 (char& JOB, char& EQUIL,
                  int& N, int& M, int& P,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  int& NR,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_tb01pd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tb01pd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TB01PD Release 5.0\n\
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
        char job = 'M';
        char equil;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        double tol = args(3).double_value ();
        const int scaled = args(4).int_value ();
        
        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc;

        if (n == 0)
            ldc = 1;
        else
            ldc = max (1, m, p);

        b.resize (ldb, max (m, p));

        // arguments out
        int nr = 0;
        
        // workspace
        int liwork = n + max (m, p);
        int ldwork = max (1, n + max (n, 3*m, 3*p));

        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int info = 0;


        // SLICOT routine TB01PD
        F77_XFCN (tb01pd, TB01PD,
                 (job, equil,
                  n, m, p,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  nr,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("ss: minreal: __sl_tb01pd__: exception in SLICOT subroutine TB01PD");
            
        if (info != 0)
            error ("ss: minreal: __sl_tb01pd__: TB01PD returned info = %d", info);

        // resize
        a.resize (nr, nr);
        b.resize (nr, m);
        c.resize (p, nr);

        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = octave_value (nr);
    }
    
    return retval;
}
