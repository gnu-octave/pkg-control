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

Solution of generalized Lyapunov equations.
Uses SLICOT SG03AD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: January 2010
Version: 0.3

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sg03ad, SG03AD)
                 (char& DICO, char& JOB,
                  char& FACT, char& TRANS,
                  char& UPLO,
                  F77_INT& N,
                  double* A, F77_INT& LDA,
                  double* E, F77_INT& LDE,
                  double* Q, F77_INT& LDQ,
                  double* Z, F77_INT& LDZ,
                  double* X, F77_INT& LDX,
                  double& SCALE,
                  double& SEP, double& FERR,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sg03ad__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sg03ad__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_sg03ad__ (@dots{})\n"
   "Wrapper for SLICOT function SG03AD.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 4)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char job = 'X';
        char fact = 'N';
        char trans = 'T';
        char uplo = 'U';    // ?!?
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix x = args(2).matrix_value ();
        F77_INT discrete = args(3).int_value ();
        
        if (a.any_element_is_inf_or_nan () || e.any_element_is_inf_or_nan () ||
            x.any_element_is_inf_or_nan ())
        error ("__sl_sg03ad__: inputs must not contain NaN or Inf\n");

        if (discrete == 0)
          dico = 'C';
        else
          dico = 'D';
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        
        F77_INT lda = max (1, n);
        F77_INT lde = max (1, n);
        F77_INT ldq = max (1, n);
        F77_INT ldz = max (1, n);
        F77_INT ldx = max (1, n);
                
        // arguments out
        double scale;
        double sep = 0;
        double ferr = 0;
        
        Matrix q (ldq, n);
        Matrix z (ldz, n);
        ColumnVector alphar (n);
        ColumnVector alphai (n);
        ColumnVector beta (n);
        
        // workspace
        F77_INT* iwork = 0;  // not referenced because job = X
        
        F77_INT ldwork = max (1, 8*n+16);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        F77_INT info;
        

        // SLICOT routine SG03AD
        F77_XFCN (sg03ad, SG03AD,
                 (dico, job,
                  fact, trans,
                  uplo,
                  n,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  q.fortran_vec (), ldq,
                  z.fortran_vec (), ldz,
                  x.fortran_vec (), ldx,
                  scale,
                  sep, ferr,
                  alphar.fortran_vec (), alphai.fortran_vec (),
                  beta.fortran_vec (),
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("lyap: __sl_sg03ad__: exception in SLICOT subroutine SG03AD");

        if (info != 0)
            error ("lyap: __sl_sg03ad__: SG03AD returned info = %d", static_cast<int> (info));

        // return values
        retval(0) = x;
        retval(1) = octave_value (scale);
    }
    
    return retval;
}
