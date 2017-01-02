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

Square-root solver for Lyapunov equations.
Uses SLICOT SB03OD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: January 2010
Version: 0.3

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb03od, SB03OD)
                 (char& DICO, char& FACT, char& TRANS,
                  F77_INT& N, F77_INT& M,
                  double* A, F77_INT& LDA,
                  double* Q, F77_INT& LDQ,
                  double* B, F77_INT& LDB,
                  double& SCALE,
                  double* WR, double* WI,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb03od__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb03od__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB03OD Release 5.0\n\
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
        char fact = 'N';
        char trans = 'N';
        // char trans = 'T';
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        F77_INT discrete = args(2).int_value ();
        
        if (discrete == 0)
          dico = 'C';
        else
          dico = 'D';
        
        F77_INT n = TO_F77_INT (a.rows ());
        F77_INT m = TO_F77_INT (b.rows ());
        // F77_INT m = TO_F77_INT (b.columns ());
        
        F77_INT lda = max (1, n);
        F77_INT ldq = max (1, n);
        F77_INT ldb = max (1, n, m);
        // F77_INT ldb = max (1, n);

        b.resize (ldb, n);
        // b.resize (ldb, max (m, n));
                
        // arguments out
        double scale;

        Matrix q (ldq, n);
        ColumnVector wr (n);
        ColumnVector wi (n);
        
        // workspace
        F77_INT ldwork = max (1, 4*n + min (m, n));
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        F77_INT info;
        

        // SLICOT routine SB03OD
        F77_XFCN (sb03od, SB03OD,
                 (dico, fact, trans,
                  n, m,
                  a.fortran_vec (), lda,
                  q.fortran_vec (), ldq,
                  b.fortran_vec (), ldb,
                  scale,
                  wr.fortran_vec (), wi.fortran_vec (),
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("lyapchol: __sl_sb03od__: exception in SLICOT subroutine SB03OD");

        if (info != 0)
            error ("lyapchol: __sl_sb03od__: SB03OD returned info = %d", info);

        // resize
        b.resize (n, n);
        
        // return values
        retval(0) = b;  // b has been overwritten by cholesky factor u
        retval(1) = octave_value (scale);
    }
    
    return retval;
}
