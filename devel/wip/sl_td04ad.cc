/*

Copyright (C) 2009-2015   Lukas F. Reichlin

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
Version: 0.4

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (td04ad, TD04AD)
                 (char& ROWCOL,
                  octave_idx_type& M, octave_idx_type& P,
                  octave_idx_type* INDEX,
                  double* DCOEFF, octave_idx_type& LDDCOE,
                  double* UCOEFF, octave_idx_type& LDUCO1, octave_idx_type& LDUCO2,
                  octave_idx_type& NR,
                  double* A, octave_idx_type& LDA,
                  double* B, octave_idx_type& LDB,
                  double* C, octave_idx_type& LDC,
                  double* D, octave_idx_type& LDD,
                  double& TOL,
                  octave_idx_type* IWORK,
                  double* DWORK, octave_idx_type& LDWORK,
                  octave_idx_type& INFO);
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

    if (nargin != 5)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        NDArray ucoeff = args(0).array_value ();
        Matrix dcoeff = args(1).matrix_value ();
        Matrix indexd = args(2).matrix_value ();
        double tol = args(3).double_value ();
        octave_idx_type rowcol_sel = args(4).int_value ();
        
        bool isrow; 
        char rowcol;
        if (rowcol_sel == 1)
        {
            rowcol = 'R';
            isrow = true;
        }
        else if (rowcol_sel == 0)
        {
            rowcol = 'C';
            isrow = false;
        }
        else
            error ("__sl_td04ad__: fifth argument can only be 1 or 0");

        octave_idx_type p = ucoeff.rows ();      // p: number of outputs
        octave_idx_type m = ucoeff.columns ();   // m: number of inputs

        octave_idx_type porm;
        octave_idx_type lddcoe;                  // TODO: handle case ucoeff.rows = 0
        octave_idx_type lduco1;
        octave_idx_type lduco2;

        if (isrow)
        {
            porm = p;
            lddcoe = max (1, p);
            lduco1 = max (1, p);
            lduco2 = max (1, m);
        }
        else
        {
            porm = m;
            lddcoe = max (1, m);
            lduco1 = max (1, m, p);
            lduco2 = max (1, m, p);
        }

        octave_idx_type n = 0;
        OCTAVE_LOCAL_BUFFER (octave_idx_type, index, porm);
        for (octave_idx_type i = 0; i < porm; i++)
        {
            index[i] = indexd.xelem (i);
            n += index[i];
        }

        // arguments out
        octave_idx_type nr = max (1, n);        // initialize to prevent crash if  info != 0
        octave_idx_type lda = max (1, n);
        octave_idx_type ldb = max (1, n);
        octave_idx_type ldc = max (1, m, p);
        octave_idx_type ldd, dcol;
        if (isrow)
        {
            ldd = max (1, p);
            dcol = m;
        }
        else
        {
            // concat ucoeff to match the dimensions
            // use insert (is the same as concat)
            octave_idx_type kdcoef = ucoeff.dim3 ();
            if (m > p)
            {
                dim_vector dv (m-p, p, kdcoef);
                ucoeff.insert (NDArray (dv, 0.), p, 0);
            }
            else if (p > m)
            {
                dim_vector dv (p, p-m, kdcoef);
                ucoeff.insert (NDArray (dv, 0.), 0, m); 
            }
            ldd = max (1, m, p);
            dcol = max (m, p);
        }
        
        Matrix a (lda, n);
        Matrix b (ldb, max (m, p));
        Matrix c (ldc, n);
        Matrix d (ldd, dcol);

        // workspace
        octave_idx_type ldwork = max (1, n + max (n, 3*m, 3*p));

        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, n + max (m, p));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        octave_idx_type info;


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
