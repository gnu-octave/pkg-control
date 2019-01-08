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

Minimal state-space representation (A,B,C,D) for a
proper transfer matrix T(s) given as either row or column
polynomial vectors over denominator polynomials, possibly with
uncancelled common terms.
Uses SLICOT TD04AD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: August 2011
Version: 0.3

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (td04ad, TD04AD)
                 (char& ROWCOL,
                  F77_INT& M, F77_INT& P,
                  F77_INT* INDEX,
                  double* DCOEFF, F77_INT& LDDCOE,
                  double* UCOEFF, F77_INT& LDUCO1, F77_INT& LDUCO2,
                  F77_INT& NR,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_td04ad__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_td04ad__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TD04AD Release 5.0\n\
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
        char rowcol = 'R';

        NDArray ucoeff = args(0).array_value ();
        Matrix dcoeff = args(1).matrix_value ();
        Matrix indexd = args(2).matrix_value ();
        double tol = args(3).double_value ();

        F77_INT p = TO_F77_INT (ucoeff.rows ());      // p: number of outputs
        F77_INT m = TO_F77_INT (ucoeff.columns ());   // m: number of inputs

        F77_INT lddcoe = max (1, p);     // TODO: handle case ucoeff.rows = 0
        F77_INT lduco1 = max (1, p);
        F77_INT lduco2 = max (1, m);
        
        F77_INT n = 0;
        OCTAVE_LOCAL_BUFFER (F77_INT, index, p);
        for (F77_INT i = 0; i < p; i++)
        {
            index[i] = TO_F77_INT (indexd.xelem (i));
            n += index[i];
        }

        // arguments out
        F77_INT nr = max (1, n);        // initialize to prevent crash if  info != 0
        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, m, p);
        F77_INT ldd = max (1, p);
        
        Matrix a (lda, n);
        Matrix b (ldb, max (m, p));
        Matrix c (ldc, n);
        Matrix d (ldd, m);

        // workspace
        F77_INT ldwork = max (1, n + max (n, 3*m, 3*p));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, n + max (m, p));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        F77_INT info;


        // SLICOT routine TD04AD
        F77_XFCN (td04ad, TD04AD,
                 (rowcol,
                  m, p,
                  index,
                  dcoeff.fortran_vec (), lddcoe,
                  ucoeff.fortran_vec (), lduco1, lduco2,
                  nr,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("tf2ss: __sl_td04ad__: exception in SLICOT subroutine TD04AD");

        if (info != 0)
            error ("tf2ss: __sl_td04ad__: TD04AD returned info = %d", info);

        // resize
        a.resize (nr, nr);
        b.resize (nr, m);
        c.resize (p, nr);
        d.resize (p, m);

        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = d;
    }

    return retval;
}
