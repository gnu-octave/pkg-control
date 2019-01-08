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

Orthogonal canonical form.
Uses SLICOT TB01UD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb01ud, TB01UD)
                 (char& JOBZ,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  F77_INT& NCONT, F77_INT& INDCON,
                  F77_INT* NBLK,
                  double* Z, F77_INT& LDZ,
                  double* TAU,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_tb01ud__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tb01ud__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TB01UD Release 5.0\n\
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
        char jobz = 'I';
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        double tol = args(3).double_value ();

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldz = max (1, n);

        // arguments out
        Matrix z (ldz, n);

        F77_INT ncont;
        F77_INT indcon;

        OCTAVE_LOCAL_BUFFER (F77_INT, nblk, n);
        OCTAVE_LOCAL_BUFFER (double, tau, n);
        
        // workspace
        F77_INT ldwork = max (1, n, 3*m, p);

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT info;


        // SLICOT routine TB01UD
        F77_XFCN (tb01ud, TB01UD,
                 (jobz,
                  n, m, p,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  ncont, indcon,
                  nblk,
                  z.fortran_vec (), ldz,
                  tau,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_tb01ud__: exception in SLICOT subroutine TB01UD");
            
        if (info != 0)
            error ("__sl_tb01ud__: TB01UD returned info = %d", info);

        // resize
        a.resize (n, n);
        b.resize (n, m);
        c.resize (p, n);
        z.resize (n, n);

        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = z;
        retval(4) = octave_value (ncont);
    }
    
    return retval;
}
