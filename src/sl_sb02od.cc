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

Solve algebraic Riccati equation.
Uses SLICOT SB02OD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: February 2010
Version: 0.5

*/

#include <octave/oct.h>
#include "common.h"
#include <complex>

extern "C"
{ 
    int F77_FUNC (sb02od, SB02OD)
                 (char& DICO, char& JOBB,
                  char& FACT, char& UPLO,
                  char& JOBL, char& SORT,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* Q, F77_INT& LDQ,
                  double* R, F77_INT& LDR,
                  double* L, F77_INT& LDL,
                  double& RCOND,
                  double* X, F77_INT& LDX,
                  double* ALFAR, double* ALFAI,
                  double* BETA,
                  double* S, F77_INT& LDS,
                  double* T, F77_INT& LDT,
                  double* U, F77_INT& LDU,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_LOGICAL* BWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb02od__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb02od__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB02OD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
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
        F77_INT discrete = args(5).int_value ();
        F77_INT ijobl = args(6).int_value ();

        if (discrete == 0)
            dico = 'C';
        else
            dico = 'D';

        if (ijobl == 0)
            jobl = 'Z';
        else
            jobl = 'N';

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = 0;              // p: number of outputs, not used because FACT = 'N'

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldq = max (1, n);
        F77_INT ldr = max (1, m);
        F77_INT ldl = max (1, n);

        // arguments out
        double rcond;

        F77_INT ldx = max (1, n);
        Matrix x (ldx, n);

        F77_INT nu = 2*n;
        ColumnVector alfar (nu);
        ColumnVector alfai (nu);
        ColumnVector beta (nu);

        F77_INT lds = max (1, 2*n + m);
        Matrix s (lds, lds);

        // unused output arguments
        F77_INT ldt = max (1, 2*n + m);
        OCTAVE_LOCAL_BUFFER (double, t, ldt * 2*n);

        F77_INT ldu = max (1, 2*n);
        OCTAVE_LOCAL_BUFFER (double, u, ldu * 2*n);

        // tolerance
        double tol = 0;  // use default value

        // workspace
        F77_INT liwork = max (1, m, 2*n);
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);

        F77_INT ldwork = max (7*(2*n + 1) + 16, 16*n, 2*n + m, 3*m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        OCTAVE_LOCAL_BUFFER (F77_LOGICAL, bwork, 2*n);

        // error indicator
        F77_INT info;


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
            error ("are: __sl_sb02od__: exception in SLICOT subroutine SB02OD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the computed extended matrix pencil is singular, "
                "possibly due to rounding errors",
            "2: the QZ (or QR) algorithm failed",
            "3: reordering of the (generalized) eigenvalues "
                "failed",
            "4: after reordering, roundoff changed values of "
                "some complex eigenvalues so that leading eigenvalues "
                "in the (generalized) Schur form no longer satisfy "
                "the stability condition; this could also be caused "
                "due to scaling",
            "5: the computed dimension of the solution does not "
                "equal N",
            "6: a singular matrix was encountered during the "
                "computation of the solution matrix X"};

        error_msg ("are", info, 6, err_msg);


        // assemble complex vector - adapted from DEFUN complex in data.cc
        alfar.resize (n);
        alfai.resize (n);
        beta.resize (n);

        ColumnVector poler (n);
        ColumnVector polei (n);

        poler = quotient (alfar, beta);
        polei = quotient (alfai, beta);
        
        ComplexColumnVector pole (n, Complex ());

        for (F77_INT i = 0; i < n; i++)
            pole.xelem (i) = Complex (poler(i), polei(i));

        // return value
        retval(0) = x;
        retval(1) = pole;
    }

    return retval;
}
