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

Gain of descriptor state-space models.  Based on SLICOT TB04BX.f.

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: March 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"
#include <complex>
#include <xpow.h>

extern "C"
{ 
    int F77_FUNC (tg04bx, TG04BX)
                 (F77_INT& IP, F77_INT& IZ,
                  double* A, F77_INT& LDA,
                  double* E,
                  double* B,
                  double* C,
                  double* D,
                  double* PR, double* PI,
                  double* ZR, double* ZI,
                  double& GAIN,
                  F77_INT* IWORK);
}

// PKG_ADD: autoload ("__sl_tg04bx__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg04bx__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_tg04bx__ (@dots{})\n"
   "Wrapper for SLICOT function TG04BX.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 9)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value (); 
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        Matrix d = args(4).matrix_value ();
        
        ColumnVector pr = args(5).column_vector_value ();
        ColumnVector pi = args(6).column_vector_value ();
        
        ColumnVector zr = args(7).column_vector_value ();
        ColumnVector zi = args(8).column_vector_value ();

        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan () || d.any_element_is_inf_or_nan () ||
            e.any_element_is_inf_or_nan ())
        error ("__sl_tg04bx__: inputs must not contain NaN or Inf\n");

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT ip = TO_F77_INT (pr.numel ());  // ip: number of finite poles
        F77_INT iz = TO_F77_INT (zr.numel ());  // iz: number of zeros
        
        // For ss, IP = n is always true.
        // However, dss models with poles at infinity
        // (filtered by pole.m) may have IP <= n
        
        // Take pr.length == pi.length == ip for granted,
        // and the same for iz, zr and zi.
        
        F77_INT lda = max (1, n);

        // arguments out
        double gain;

        // workspace
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, lda);

        
        F77_XFCN (tg04bx, TG04BX,
                 (ip, iz,
                  a.fortran_vec (), lda,
                  e.fortran_vec (),
                  b.fortran_vec (),
                  c.fortran_vec (),
                  d.fortran_vec (),
                  pr.fortran_vec (), pi.fortran_vec (),
                  zr.fortran_vec (), zi.fortran_vec (),
                  gain,
                  iwork));
                  
        if (f77_exception_encountered)
            error ("dss: zero: __sl_tg04bx__: exception in TG04BX");

        // return values
        retval(0) = octave_value (gain);
    }
    
    return retval;
}
