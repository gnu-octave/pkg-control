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

Hankel singular values.
Uses SLICOT AB13AD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: January 2010
Version: 0.4

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab13ad, AB13AD)
                 (char& DICO, char& EQUIL,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double& ALPHA,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  F77_INT& NS,
                  double* HSV,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ab13ad__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab13ad__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB13AD Release 5.0\n\
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
        char dico;
        char equil;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        F77_INT discrete = args(3).int_value ();
        double alpha = args(4).double_value ();
        const F77_INT scaled = args(5).int_value ();
        
        if (discrete == 0)
            dico = 'C';
        else
            dico = 'D';

        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        
        F77_INT lda = max (1, TO_F77_INT (a.rows ()));
        F77_INT ldb = max (1, TO_F77_INT (b.rows ()));
        F77_INT ldc = max (1, TO_F77_INT (c.rows ()));
        
        // arguments out
        F77_INT ns = 0;
        
        ColumnVector hsv (n);
        
        // workspace
        F77_INT ldwork = max (1, n*(max (n, m, p) + 5) + n*(n+1)/2);
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT info = 0;


        // SLICOT routine AB13AD
        F77_XFCN (ab13ad, AB13AD,
                 (dico, equil,
                  n, m, p,
                  alpha,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  ns,
                  hsv.fortran_vec (),
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("hsvd: __sl_ab13ad__: exception in SLICOT subroutine AB13AD");
            
        if (info != 0)
            error ("hsvd: __sl_ab13ad__: AB13AD returned info = %d", info);

        // resize
        hsv.resize (ns);
        
        // return values
        retval(0) = hsv;
        retval(1) = octave_value (ns);
    }
    
    return retval;
}
