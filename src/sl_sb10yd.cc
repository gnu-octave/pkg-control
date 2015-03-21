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

Fit FRD with SS model.
Uses SLICOT SB10YD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2011
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include <complex>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10yd, SB10YD)
                 (octave_idx_type& DISCFL, octave_idx_type& FLAG,
                  octave_idx_type& LENDAT,
                  double* RFRDAT, double* IFRDAT,
                  double* OMEGA,
                  octave_idx_type& N,
                  double* A, octave_idx_type& LDA,
                  double* B,
                  double* C,
                  double* D,
                  double& TOL,
                  octave_idx_type* IWORK, double* DWORK, octave_idx_type& LDWORK,
                  Complex* ZWORK, octave_idx_type& LZWORK,
                  octave_idx_type& INFO);
}

// PKG_ADD: autoload ("__sl_sb10yd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10yd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10YD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 6)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        Matrix rfrdat = args(0).matrix_value ();
        Matrix ifrdat = args(1).matrix_value ();
        Matrix omega = args(2).matrix_value ();
        octave_idx_type n = args(3).int_value ();
        octave_idx_type discfl = args(4).int_value ();
        octave_idx_type flag = args(5).int_value ();
        
        octave_idx_type lendat = omega.rows ();      // number of frequencies

        octave_idx_type lda = max (1, n);
        
        // arguments out
        Matrix a (lda, n);
        Matrix b (n, 1);
        Matrix c (1, n);
        Matrix d (1, 1);
        
        // workspace
        octave_idx_type liwork = max (2, 2*n + 1);
        octave_idx_type ldwork;
        octave_idx_type lzwork;
        octave_idx_type hnpts = 2048;
        
        octave_idx_type lw1 = 2*lendat + 4*hnpts;
        octave_idx_type lw2 =   lendat + 6*hnpts;
        octave_idx_type mn  = min (2*lendat, 2*n+1);
        octave_idx_type lw3;
        octave_idx_type lw4;
        
        if (n > 0)
        {
            lw3 = 2*lendat*(2*n+1) + max (2*lendat, 2*n+1)
                  + max (mn + 6*n + 4, 2*mn + 1);
            lzwork = lendat*(2*n+3);
        }
        else
        {
            lw3 = 4*lendat + 5;
            lzwork = lendat;
        }

        if (flag == 1)
            lw4 = max (n*n + 5*n, 6*n + 1 + min (1, n));
        else
            lw4 = 0;
       
        ldwork = max (2, lw1, lw2, lw3, lw4);
        
        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (Complex, zwork, lzwork);
        
        // tolerance
        double tol = 0;
        
        // error indicator
        octave_idx_type info;


        // SLICOT routine SB10YD
        F77_XFCN (sb10yd, SB10YD,
                 (discfl, flag,
                  lendat,
                  rfrdat.fortran_vec (), ifrdat.fortran_vec (),
                  omega.fortran_vec (),
                  n,
                  a.fortran_vec (), lda,
                  b.fortran_vec (),
                  c.fortran_vec (),
                  d.fortran_vec (),
                  tol,
                  iwork, dwork, ldwork,
                  zwork, lzwork,
                  info));

        if (f77_exception_encountered)
            error ("fitfrd: __sl_sb10yd__: exception in SLICOT subroutine SB10YD");
            
        if (info != 0)
            error ("fitfrd: __sl_sb10yd__: SB10YD returned info = %d", info);
        
        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = d;
        retval(4) = octave_value (n);
    }
    
    return retval;
}
