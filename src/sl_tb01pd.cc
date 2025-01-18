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

Minimal realization of state-space models.
Uses SLICOT TB01PD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.5

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb01pd, TB01PD)
                 (char& JOB, char& EQUIL,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  F77_INT& NR,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_tb01pd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tb01pd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_tb01pd__ (@dots{})\n"
   "Wrapper for SLICOT function TB01PD.@*\n"
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
        char job = 'M';
        char equil;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        double tol = args(3).double_value ();
        const F77_INT scaled = args(4).int_value ();
        
        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan ())
        error ("__sl_tb01pd__: inputs must not contain NaN or Inf\n");

        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc;

        if (n == 0)
            ldc = 1;
        else
            ldc = max (1, m, p);

        b.resize (ldb, max (m, p));
        c.resize (ldc, n);

        // arguments out
        F77_INT nr = 0;
        
        // workspace
        F77_INT liwork = n + max (m, p);
        F77_INT ldwork = max (1, n + max (n, 3*m, 3*p));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT info = 0;


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
            error ("ss: minreal: __sl_tb01pd__: TB01PD returned info = %d", static_cast<int> (info));

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
