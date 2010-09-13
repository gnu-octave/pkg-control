/*

Copyright (C) 2010   Lukas F. Reichlin

This file is part of LTI Syncope.

LTI Syncope is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LTI Syncope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Solve algebraic Riccati equation.
Uses SLICOT SB02OD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: February 2010
Version: 0.2.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"
#include <complex>

extern "C"
{ 
    int F77_FUNC (sb02od, SB02OD)
                 (char& DICO, char& JOBB,
                  char& FACT, char& UPLO,
                  char& JOBL, char& SORT,
                  int& N, int& M, int& P,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* Q, int& LDQ,
                  double* R, int& LDR,
                  double* L, int& LDL,
                  double& RCOND,
                  double* X, int& LDX,
                  double* ALFAR, double* ALFAI,
                  double* BETA,
                  double* S, int& LDS,
                  double* T, int& LDT,
                  double* U, int& LDU,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  bool* BWORK,
                  int& INFO);
}

DEFUN_DLD (slsb02od, args, nargout, "Slicot SB02OD Release 5.0")
{
    int nargin = args.length ();
    octave_value_list retval;

    if (nargin != 7)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char jobb = 'B';
        char fact = 'N';
        char uplo = 'U';
        char jobl;
        char sort = 'S';

        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix q = args(2).matrix_value ();
        Matrix r = args(3).matrix_value ();
        Matrix l = args(4).matrix_value ();
        int discrete = args(5).int_value ();
        int ijobl = args(6).int_value ();

        if (discrete == 0)
            dico = 'C';
        else
            dico = 'D';

        if (ijobl == 0)
            jobl = 'Z';
        else
            jobl = 'N';

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = 0;              // p: number of outputs, not used because FACT = 'N'

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldq = max (1, n);
        int ldr = max (1, m);
        int ldl = max (1, n);

        // arguments out
        double rcond;

        int ldx = max (1, n);
        Matrix x (ldx, n);

        int nu = 2*n;
        ColumnVector alfar (nu);
        ColumnVector alfai (nu);
        ColumnVector beta (nu);

        int lds = max (1, 2*n + m);
        Matrix s (lds, lds);

        // unused output arguments
        int ldt = max (1, 2*n + m);
        OCTAVE_LOCAL_BUFFER (double, t, ldt * 2*n);

        int ldu = max (1, 2*n);
        OCTAVE_LOCAL_BUFFER (double, u, ldu * 2*n);

        // tolerance
        double tol = 0;  // use default value

        // workspace
        int liwork = max (1, m, 2*n);
        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);

        int ldwork = max (7*(2*n + 1) + 16, 16*n, 2*n + m, 3*m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        OCTAVE_LOCAL_BUFFER (bool, bwork, 2*n);

        // error indicator
        int info;


        // SLICOT routine SB02OD
        F77_XFCN (sb02od, SB02OD,
                 (dico, jobb,
                  fact, uplo,
                  jobl, sort,
                  n, m, p,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  q.fortran_vec (), ldq,
                  r.fortran_vec (), ldr,
                  l.fortran_vec (), ldl,
                  rcond,
                  x.fortran_vec (), ldx,
                  alfar.fortran_vec (), alfai.fortran_vec (),
                  beta.fortran_vec (),
                  s.fortran_vec (), lds,
                  t, ldt,
                  u, ldu,
                  tol,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("are: slsb02od: exception in SLICOT subroutine SB02OD");

        if (info != 0)
            error ("are: slsb02od: SB02OD returned info = %d", info);

        // assemble complex vector - adapted from DEFUN complex in data.cc
        alfar.resize (n);
        alfai.resize (n);
        beta.resize (n);

        ColumnVector zeror (n);
        ColumnVector zeroi (n);

        zeror = quotient (alfar, beta);
        zeroi = quotient (alfai, beta);
        
        ComplexColumnVector zero (n, Complex ());

        for (octave_idx_type i = 0; i < n; i++)
            zero.xelem (i) = Complex (zeror(i), zeroi(i));

        // return value
        retval(0) = x;
        retval(1) = zero;
    }

    return retval;
}
