/*

Copyright (C) 2025 Torsten Lilge

This function is part of the GNU Octave Control Package

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, see
<http://www.gnu.org/licenses/>.

Wrapper for the SLICOT function MB03RD for making a block diagonal
matrix from a schur matrix.

Author: Torsten Lilge <ttl-octave@mailbox.org>
Created: June 2025

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (mb03rd, MB03RD)
                 (char& JOB, char& SORT,
                  F77_INT& N,
                  double& PMAX,
                  double* A, F77_INT& LDA,
                  double* X, F77_INT& LDX,
                  F77_INT& NBLCKS,
                  F77_INT* BLSIZE,
                  double* WR,
                  double* WI,
                  double& TOL,
                  double* DWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_mb03rd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_mb03rd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_mb03rd__ (@dots{})\n"
   "Wrapper for SLICOT function MB03RD.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
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
        char job = 'U';
        char sort = 'N';
        double pmax = args(2).double_value ();
        Matrix a = args(1).matrix_value ();
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT lda = max (1, n);
        Matrix x = args(0).matrix_value ();
        F77_INT ldx = max (1, n);
        
        if (a.any_element_is_inf_or_nan () || x.any_element_is_inf_or_nan ())
          error ("__sl_mb03rd__: inputs must not contain NaN or Inf\n");

        F77_INT nblcks;
        OCTAVE_LOCAL_BUFFER (F77_INT, blsize, n);
        
        OCTAVE_LOCAL_BUFFER (double, wr, n);
        OCTAVE_LOCAL_BUFFER (double, wi, n);

        double tol = 0.;

        OCTAVE_LOCAL_BUFFER (double, dwork, n);
        
        // error indicator
        F77_INT info;

        // SLICOT routine MB03RD
        F77_XFCN (mb03rd, MB03RD,
                 (job, sort,
                  n,
                  pmax,
                  a.fortran_vec (), lda,
                  x.fortran_vec (), ldx,
                  nblcks,
                  blsize,
                  wr, wi,
                  tol,
                  dwork,
                  info));

        if (f77_exception_encountered)
          error ("bdschur: __sl_mb03rd__: exception in SLICOT subroutine MB03RD");

        if (info != 0)
          error ("__sl_tg01hd__: TG01HD returned info = %d", static_cast<int> (info));
     
        // return values
        RowVector blksz (nblcks);
        for (F77_INT i = 0; i < nblcks; i++)
          blksz.xelem (i) = blsize[i];
        
        // RowVector wreal (n);
        // RowVector wimag (n);
        // for (F77_INT i = 0; i < n; i++)
        //   {
        //     wreal.xelem (i) = wr[i];
        //     wimag.xelem (i) = wi[i];
        //   }

        retval(0) = x;
        retval(1) = a;
        retval(2) = blksz;
        // retval(3) = wreal;
        // retval(4) = wimag;
    }
    
    return retval;
}
