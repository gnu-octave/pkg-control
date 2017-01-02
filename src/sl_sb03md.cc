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

Solution of Lyapunov equations.
Uses SLICOT SB03MD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: December 2009
Version: 0.4

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb03md, SB03MD)
                 (char& DICO, char& JOB,
                  char& FACT, char& TRANA,
                  F77_INT& N,
                  double* A, F77_INT& LDA,
                  double* U, F77_INT& LDU,
                  double* C, F77_INT& LDC,
                  double& SCALE,
                  double& SEP, double& FERR,
                  double* WR, double* WI,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb03md__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb03md__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB03MD Release 5.0\n\
No argument checking.\n\
For internal use only.")
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
        char dico;
        char job = 'X';
        char fact = 'N';
        char trana = 'T';
        
        Matrix a = args(0).matrix_value ();
        Matrix c = args(1).matrix_value ();
        F77_INT discrete = args(2).int_value ();
        
        if (discrete == 0)
          dico = 'C';
        else
          dico = 'D';
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        
        F77_INT lda = max (1, n);
        F77_INT ldu = max (1, n);
        F77_INT ldc = max (1, n);
        
        // arguments out
        double scale;
        double sep = 0;
        double ferr = 0;
        
        Matrix u (ldu, n);
        ColumnVector wr (n);
        ColumnVector wi (n);
        
        // workspace
        F77_INT* iwork = 0;  // not referenced because job = X
        
        F77_INT ldwork = max (n*n, 3*n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        F77_INT info;
        

        // SLICOT routine SB03MD
        F77_XFCN (sb03md, SB03MD,
                 (dico, job,
                  fact, trana,
                  n,
                  a.fortran_vec (), lda,
                  u.fortran_vec (), ldu,
                  c.fortran_vec (), ldc,
                  scale,
                  sep, ferr,
                  wr.fortran_vec (), wi.fortran_vec (),
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("lyap: __sl_sb03md__: exception in SLICOT subroutine SB03MD");

        if (info != 0)
            error ("lyap: __sl_sb03md__: SB03MD returned info = %d", info);
        
        // return values
        retval(0) = c;
        retval(1) = octave_value (scale);
    }
    
    return retval;
}
