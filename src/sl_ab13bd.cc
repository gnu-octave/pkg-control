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

H-2 norm of a SS model.
Uses SLICOT AB13BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.5

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    double F77_FUNC (ab13bd, AB13BD)
                    (char& DICO, char& JOBN,
                     F77_INT& N, F77_INT& M, F77_INT& P,
                     double* A, F77_INT& LDA,
                     double* B, F77_INT& LDB,
                     double* C, F77_INT& LDC,
                     double* D, F77_INT& LDD,
                     F77_INT& NQ,
                     double& TOL,
                     double* DWORK, F77_INT& LDWORK,
                     F77_INT& IWARN,
                     F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ab13bd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab13bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB13BD Release 5.\n\
No argument checking.\n\
For internal use only.")
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
        char dico;
        char jobn = 'H';
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        F77_INT discrete = args(4).int_value ();

        if (discrete == 0)
            dico = 'C';
        else
            dico = 'D';
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        
        F77_INT lda = max (1, TO_F77_INT (a.rows ()));
        F77_INT ldb = max (1, TO_F77_INT (b.rows ()));
        F77_INT ldc = max (1, TO_F77_INT (c.rows ()));
        F77_INT ldd = max (1, TO_F77_INT (d.rows ()));
        
        // arguments out
        double norm;
        F77_INT nq;
        
        // tolerance
        double tol = 0;
        
        // workspace
        F77_INT ldwork = max (1, m*(n+m) + max (n*(n+5), m*(m+2), 4*p ),
                             n*(max (n, p) + 4 ) + min (n, p));

        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        F77_INT iwarn;
        F77_INT info;


        // SLICOT routine AB13BD
        norm = F77_FUNC (ab13bd, AB13BD)
                        (dico, jobn,
                         n, m, p,
                         a.fortran_vec (), lda,
                         b.fortran_vec (), ldb,
                         c.fortran_vec (), ldc,
                         d.fortran_vec (), ldd,
                         nq,
                         tol,
                         dwork, ldwork,
                         iwarn,
                         info);

        if (f77_exception_encountered)
            error ("lti: norm: __sl_ab13bd__: exception in SLICOT subroutine AB13BD");
            
        if (info != 0)
            error ("lti: norm: __sl_ab13bd__: AB13BD returned info = %d", info);

        if (iwarn != 0)
            warning ("lti: norm: __sl_ab13bd__: AB13BD returned iwarn = %d", iwarn);
        
        // return value
        retval(0) = octave_value (norm);
    }
    
    return retval;
}
