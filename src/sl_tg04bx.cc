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

Gain of descriptor state-space models.  Based on SLICOT TB04BX.f.

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: March 2011
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"
#include <complex>
#include <xpow.h>

extern "C"
{ 
    int F77_FUNC (tg04bx, TG04BX)
                 (int& IP, int& IZ,
                  double* A, int& LDA,
                  double* E,
                  double* B,
                  double* C,
                  double* D,
                  double* PR, double* PI,
                  double* ZR, double* ZI,
                  double& GAIN,
                  int* IWORK);
}

// PKG_ADD: autoload ("__sl_tg04bx__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg04bx__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TG04BX Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
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

        int n = a.rows ();      // n: number of states
        int ip = pr.length ();  // ip: number of finite poles
        int iz = zr.length ();  // iz: number of zeros
        
        // For ss, IP = n is always true.
        // However, dss models with poles at infinity
        // (filtered by pole.m) may have IP <= n
        
        // Take pr.length == pi.length == ip for granted,
        // and the same for iz, zr and zi.
        
        int lda = max (1, n);

        // arguments out
        double gain;

        // workspace
        OCTAVE_LOCAL_BUFFER (int, iwork, lda);

        
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
