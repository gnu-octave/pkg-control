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

Solution of algebraic Riccati equations for descriptor systems.
Uses SLICOT SG02AD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2010
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"
#include <complex>

extern "C"
{ 
    int F77_FUNC (sg02ad, SG02AD)
                 (char& DICO, char& JOBB,
                  char& FACT, char& UPLO,
                  char& JOBL, char& SCAL,
                  char& SORT, char& ACC,
                  int& N, int& M, int& P,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* B, int& LDB,
                  double* Q, int& LDQ,
                  double* R, int& LDR,
                  double* L, int& LDL,
                  double& RCONDU,
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
                  int& IWARN, int& INFO);
}

DEFUN_DLD (slsg02ad, args, nargout,
   "-*- texinfo -*-\n\
Slicot SG02AD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;

    if (nargin != 8)
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
        char scal = 'N';
        char sort = 'S';
        char acc = 'N';

        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix q = args(3).matrix_value ();
        Matrix r = args(4).matrix_value ();
        Matrix l = args(5).matrix_value ();
        int discrete = args(6).int_value ();
        int ijobl = args(7).int_value ();

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
        int lde = max (1, n);
        int ldb = max (1, n);
        int ldq = max (1, n);
        int ldr = max (1, m);
        int ldl = max (1, n);

        // arguments out
        double rcondu;

        int ldx = max (1, n);
        Matrix x (ldx, n);

        int nu = 2*n;
        ColumnVector alfar (nu);
        ColumnVector alfai (nu);
        ColumnVector beta (nu);

        // unused output arguments
        int lds = max (1, 2*n + m);
        OCTAVE_LOCAL_BUFFER (double, s, lds * lds);

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
        int iwarn;
        int info;


        // SLICOT routine SG02AD
        F77_XFCN (sg02ad, SG02AD,
                 (dico, jobb,
                  fact, uplo,
                  jobl, scal,
                  sort, acc,
                  n, m, p,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  q.fortran_vec (), ldq,
                  r.fortran_vec (), ldr,
                  l.fortran_vec (), ldl,
                  rcondu,
                  x.fortran_vec (), ldx,
                  alfar.fortran_vec (), alfai.fortran_vec (),
                  beta.fortran_vec (),
                  s, lds,
                  t, ldt,
                  u, ldu,
                  tol,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  iwarn, info));

        if (f77_exception_encountered)
            error ("are: slsg02ad: exception in SLICOT subroutine SG02AD");

        if (iwarn != 0)
            warning ("are: slsg02ad: solution may be inaccurate due to poor scaling or eigenvalues too close to the boundary of the stability domain");

        if (info != 0)
            error ("are: slsg02ad: SG02AD returned info = %d", info);

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
