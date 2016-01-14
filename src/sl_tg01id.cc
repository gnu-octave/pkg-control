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

Staircase observability form for descriptor models.
Uses SLICOT TG01ID by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tg01id, TG01ID)
                 (char& JOBOBS,
                  char& COMPQ, char& COMPZ,
                  octave_idx_type& N, octave_idx_type& M, octave_idx_type& P,
                  double* A, octave_idx_type& LDA,
                  double* E, octave_idx_type& LDE,
                  double* B, octave_idx_type& LDB,
                  double* C, octave_idx_type& LDC,
                  double* Q, octave_idx_type& LDQ,
                  double* Z, octave_idx_type& LDZ,
                  octave_idx_type& NOBSV, octave_idx_type& NIUOBS,
                  octave_idx_type& NLBLCK,
                  octave_idx_type* CTAU,
                  double& TOL,
                  octave_idx_type* IWORK, double* DWORK,
                  octave_idx_type& INFO);
}

// PKG_ADD: autoload ("__sl_tg01id__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg01id__, args, nargout, "Slicot TG01ID Release 5.0")
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
        char jobobs = 'O';
        char compq = 'I';
        char compz = 'I';
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        double tol = args(4).double_value ();

        octave_idx_type n = a.rows ();      // n: number of states
        octave_idx_type m = b.columns ();   // m: number of inputs
        octave_idx_type p = c.rows ();      // p: number of outputs

        octave_idx_type lda = max (1, n);
        octave_idx_type lde = max (1, n);
        octave_idx_type ldb = max (1, n);
        octave_idx_type ldc = max (1, m, p);
        octave_idx_type ldq = max (1, n);
        octave_idx_type ldz = max (1, n);

        b.resize (ldb, max (m, p));

        // arguments out
        Matrix q (ldq, n);
        Matrix z (ldz, n);

        octave_idx_type nobsv;
        octave_idx_type niuobs;
        octave_idx_type nlblck;

        OCTAVE_LOCAL_BUFFER (octave_idx_type, ctau, n);
        
        // workspace
        octave_idx_type ldwork = max (n, 2*p);

        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, p);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        octave_idx_type info;


        // SLICOT routine TG01ID
        F77_XFCN (tg01id, TG01ID,
                 (jobobs,
                  compq, compz,
                  n, m, p,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  q.fortran_vec (), ldq,
                  z.fortran_vec (), ldz,
                  nobsv, niuobs,
                  nlblck,
                  ctau,
                  tol,
                  iwork, dwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_tg01id__: exception in SLICOT subroutine TG01ID");
            
        if (info != 0)
            error ("__sl_tg01id__: TG01ID returned info = %d", info);

        // resize
        a.resize (n, n);
        e.resize (n, n);
        b.resize (n, m);
        c.resize (p, n);
        q.resize (n, n);
        z.resize (n, n);

        // return values
        retval(0) = a;
        retval(1) = e;
        retval(2) = b;
        retval(3) = c;
        retval(4) = q;
        retval(5) = z;
        retval(6) = octave_value (nobsv);
    }
    
    return retval;
}
