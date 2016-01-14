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

Staircase controllability form.
Uses SLICOT AB01OD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: August 2010
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab01od, AB01OD)
                 (char& STAGES,
                  char& JOBU, char& JOBV,
                  octave_idx_type& N, octave_idx_type& M,
                  double* A, octave_idx_type& LDA,
                  double* B, octave_idx_type& LDB,
                  double* U, octave_idx_type& LDU,
                  double* V, octave_idx_type& LDV,
                  octave_idx_type& NCONT, octave_idx_type& INDCON,
                  octave_idx_type* KSTAIR,
                  double& TOL,
                  octave_idx_type* IWORK,
                  double* DWORK, octave_idx_type& LDWORK,
                  octave_idx_type& INFO);
}
 
// PKG_ADD: autoload ("__sl_ab01od__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab01od__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB01OD Release 5.0\n\
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
        char stages = 'F';
        char jobu = 'I';
        char jobv = 'N';    // not referenced because stages = 'F'
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        double tol = args(2).double_value ();

        octave_idx_type n = a.rows ();      // n: number of states
        octave_idx_type m = b.columns ();   // m: number of inputs

        octave_idx_type lda = max (1, n);
        octave_idx_type ldb = max (1, n);
        octave_idx_type ldu = max (1, n);
        octave_idx_type ldv = 1;

        // arguments out
        Matrix u (ldu, n);
        double* v = 0;          // not referenced because stages = 'F'

        octave_idx_type ncont;
        octave_idx_type indcon;

        OCTAVE_LOCAL_BUFFER (octave_idx_type, kstair, n);
        
        // workspace
        octave_idx_type ldwork = max (1, n + max (n, 3*m));

        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        octave_idx_type info;


        // SLICOT routine AB01OD
        F77_XFCN (ab01od, AB01OD,
                 (stages,
                  jobu, jobv,
                  n, m,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  u.fortran_vec (), ldu,
                  v, ldv,
                  ncont, indcon,
                  kstair,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_ab01od__: exception in SLICOT subroutine AB01OD");
            
        if (info != 0)
            error ("__sl_ab01od__: AB01OD returned info = %d", info);

        // resize
        a.resize (n, n);
        b.resize (n, m);
        u.resize (n, n);

        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = u;
        retval(3) = octave_value (ncont);
    }
    
    return retval;
}
