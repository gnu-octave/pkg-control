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

Balance state-space model.
Uses SLICOT TB01ID by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: May 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb01id, TB01ID)
                 (char& JOB,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double& MAXRED,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* SCALE,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_tb01id__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tb01id__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TB01ID Release 5.0\n\
No argument checking.\n\
For internal use only.")
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
        char job = 'A';
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        double maxred = args(3).double_value ();

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        

        // arguments out
        ColumnVector scale (n);

        
        // error indicators
        F77_INT info = 0;


        // SLICOT routine TB01ID
        F77_XFCN (tb01id, TB01ID,
                 (job,
                  n, m, p,
                  maxred,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  scale.fortran_vec (),
                  info));

        if (f77_exception_encountered)
            error ("ss: prescale: __sl_tb01id__: exception in SLICOT subroutine TB01ID");
            
        if (info != 0)
            error ("ss: prescale: __sl_tb01id__: TB01ID returned info = %d", info);


        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = octave_value (maxred);
        retval(4) = scale;
    }
    
    return retval;
}
