/*

Copyright (C) 2011   Lukas F. Reichlin

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
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"

extern "C"
{ 
    int F77_FUNC (td04ad, TD04AD)
                 (char& ROWCOL,
                  int& M, int& P,
                  int* INDEX,
                  double* DCOEFF, int& LDDCOE
                  double* UCOEFF, int& LDUCO1, int& LDUCO2,
                  int& NR,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

DEFUN_DLD (sltd04ad, args, nargout,
   "-*- texinfo -*-\n\
Slicot TD04AD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;

    if (nargin != 5)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char rowcol = 'R';

        NDArray ucoeff = args(0).array_value ();
        Matrix dcoeff = args(1).matrix_value ();

        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        const int scaled = args(4).int_value ();

        int m = ucoeff.size (2);   // m: number of inputs
        int p = ucoeff.size (1);   // p: number of outputs
        int md = n + 1;

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        int ldd = max (1, p);

        // arguments out
        int ldign = max (1, p);
        int ldigd = max (1, p);
        int lg = p * m * md;

        OCTAVE_LOCAL_BUFFER (int, ign, ldign*m);
        OCTAVE_LOCAL_BUFFER (int, igd, ldigd*m);
        
        Matrix ignm (ldign, m);
        Matrix igdm (ldigd, m);

        RowVector gn (lg);
        RowVector gd (lg);

        // tolerance
        double tol = 0;  // use default value

        // workspace
        int ldwork = max (1, n*(n + p) + max (n + max (n, p), n*(2*n + 5)));

        OCTAVE_LOCAL_BUFFER (int, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        int info;


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
            error ("are: sltd04ad: exception in SLICOT subroutine TD04AD");

        if (info != 0)
            error ("are: sltd04ad: TD04AD returned info = %d", info);

        // return values
        retval(0) = gn;
        retval(1) = gd;
        retval(2) = ignm;
        retval(3) = igdm;
        retval(4) = octave_value (md);
        retval(5) = octave_value (p);
        retval(6) = octave_value (m);
    }

    return retval;
}
