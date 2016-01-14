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
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb01ud, TB01UD)
                 (char& JOBZ,
                  octave_idx_type& N, octave_idx_type& M, octave_idx_type& P,
                  double* A, octave_idx_type& LDA,
                  double* B, octave_idx_type& LDB,
                  double* C, octave_idx_type& LDC,
                  octave_idx_type& NCONT, octave_idx_type& INDCON,
                  octave_idx_type* NBLK,
                  double* Z, octave_idx_type& LDZ,
                  double* TAU,
                  double& TOL,
                  octave_idx_type* IWORK,
                  double* DWORK, octave_idx_type& LDWORK,
                  octave_idx_type& INFO);
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

        octave_idx_type n = a.rows ();      // n: number of states
        octave_idx_type m = b.columns ();   // m: number of inputs
        octave_idx_type p = c.rows ();      // p: number of outputs

        octave_idx_type lda = max (1, n);
        octave_idx_type ldb = max (1, n);
        octave_idx_type ldc = max (1, p);
        octave_idx_type ldz = max (1, n);

        // arguments out
        Matrix z (ldz, n);

        octave_idx_type ncont;
        octave_idx_type indcon;

        OCTAVE_LOCAL_BUFFER (octave_idx_type, nblk, n);
        OCTAVE_LOCAL_BUFFER (double, tau, n);
        
        // workspace
        octave_idx_type ldwork = max (1, n, 3*m, p);

        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        octave_idx_type info;


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
