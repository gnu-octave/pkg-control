/*

Copyright (C) 2010   Lukas F. Reichlin

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
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab01od, AB01OD)
                 (char& STAGES,
                  char& JOBU, char& JOBV,
                  int& N, int& M,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* U, int& LDU,
                  double* V, int& LDV,
                  int& NCONT, int& INDCON,
                  int* KSTAIR,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}
 
// PKG_ADD: autoload ("__sl_ab01od__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab01od__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB01OD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
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

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldu = max (1, n);
        int ldv = 1;

        // arguments out
        Matrix u (ldu, n);
        double* v = 0;          // not referenced because stages = 'F'

        int ncont;
        int indcon;

        OCTAVE_LOCAL_BUFFER (int, kstair, n);
        
        // workspace
        int ldwork = max (1, n + max (n, 3*m));

        OCTAVE_LOCAL_BUFFER (int, iwork, m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int info;


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
