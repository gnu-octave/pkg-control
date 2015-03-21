/*

Copyright (C) 2009-2015   Lukas F. Reichlin

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

Staircase controllability form for descriptor models.
Uses SLICOT TG01HD by courtesy of NICONET e.V.
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
    int F77_FUNC (tg01hd, TG01HD)
                 (char& JOBCON,
                  char& COMPQ, char& COMPZ,
                  octave_idx_type& N, octave_idx_type& M, octave_idx_type& P,
                  double* A, octave_idx_type& LDA,
                  double* E, octave_idx_type& LDE,
                  double* B, octave_idx_type& LDB,
                  double* C, octave_idx_type& LDC,
                  double* Q, octave_idx_type& LDQ,
                  double* Z, octave_idx_type& LDZ,
                  octave_idx_type& NCONT, octave_idx_type& NIUCON,
                  octave_idx_type& NRBLCK,
                  octave_idx_type* RTAU,
                  double& TOL,
                  octave_idx_type* IWORK, double* DWORK,
                  octave_idx_type& INFO);
}

// PKG_ADD: autoload ("__sl_tg01hd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg01hd__, args, nargout, "Slicot TG01HD Release 5.0")
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
        char jobcon = 'C';
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
        octave_idx_type ldc = max (1, p);
        octave_idx_type ldq = max (1, n);
        octave_idx_type ldz = max (1, n);

        // arguments out
        Matrix q (ldq, n);
        Matrix z (ldz, n);

        octave_idx_type ncont;
        octave_idx_type niucon;
        octave_idx_type nrblck;

        OCTAVE_LOCAL_BUFFER (octave_idx_type, rtau, n);
        
        // workspace
        octave_idx_type ldwork = max (n, 2*m);

        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        octave_idx_type info;


        // SLICOT routine TG01HD
        F77_XFCN (tg01hd, TG01HD,
                 (jobcon,
                  compq, compz,
                  n, m, p,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  q.fortran_vec (), ldq,
                  z.fortran_vec (), ldz,
                  ncont, niucon,
                  nrblck,
                  rtau,
                  tol,
                  iwork, dwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_tg01hd__: exception in SLICOT subroutine TG01HD");
            
        if (info != 0)
            error ("__sl_tg01hd__: TG01HD returned info = %d", info);

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
        retval(6) = octave_value (ncont);
    }
    
    return retval;
}
