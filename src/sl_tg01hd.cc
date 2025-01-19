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

Staircase controllability form for descriptor models.
Uses SLICOT TG01HD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tg01hd, TG01HD)
                 (char& JOBCON,
                  char& COMPQ, char& COMPZ,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double* A, F77_INT& LDA,
                  double* E, F77_INT& LDE,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* Q, F77_INT& LDQ,
                  double* Z, F77_INT& LDZ,
                  F77_INT& NCONT, F77_INT& NIUCON,
                  F77_INT& NRBLCK,
                  F77_INT* RTAU,
                  double& TOL,
                  F77_INT* IWORK, double* DWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_tg01hd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg01hd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_tg01hd__ (@dots{})\n"
   "Wrapper for SLICOT function TG01HD.@*\n"
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
        char jobcon = 'C';
        char compq = 'I';
        char compz = 'I';
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        double tol = args(4).double_value ();

        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan () || e.any_element_is_inf_or_nan ())
        error ("__sl_tg01hd__: inputs must not contain NaN or Inf\n");

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT lde = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldq = max (1, n);
        F77_INT ldz = max (1, n);

        // arguments out
        Matrix q (ldq, n);
        Matrix z (ldz, n);

        F77_INT ncont;
        F77_INT niucon;
        F77_INT nrblck;

        OCTAVE_LOCAL_BUFFER (F77_INT, rtau, n);
        
        // workspace
        F77_INT ldwork = max (n, 2*m);

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT info;


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
            error ("__sl_tg01hd__: TG01HD returned info = %d", static_cast<int> (info));

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
