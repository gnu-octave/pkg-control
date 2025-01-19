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
Uses SLICOT SB02RD and SB02MT by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: December 2012
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"
#include <complex>

extern "C"
{ 
    int F77_FUNC (sb02mt, SB02MT)
                 (char& JOBG, char& JOBL,
                  char& FACT, char& UPLO,
                  F77_INT& N, F77_INT& M,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* Q, F77_INT& LDQ,
                  double* R, F77_INT& LDR,
                  double* L, F77_INT& LDL,
                  F77_INT* IPIV, F77_INT& OUFACT,
                  double* G, F77_INT& LDG,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);


    int F77_FUNC (sb02rd, SB02RD)
                 (char& JOB, char& DICO,
                  char& HINV, char& TRANA,
                  char& UPLO, char& SCAL,
                  char& SORT, char& FACT,
                  char& LYAPUN,
                  F77_INT& N,
                  double* A, F77_INT& LDA,
                  double* T, F77_INT& LDT,
                  double* V, F77_INT& LDV,
                  double* G, F77_INT& LDG,
                  double* Q, F77_INT& LDQ,
                  double* X, F77_INT& LDX,
                  double& SEP,
                  double& RCOND, double& FERR,
                  double* WR, double* WI,
                  double* S, F77_INT& LDS,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_LOGICAL* BWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_are__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_are__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_sb02mt__ (@dots{})\n"
   "Wrapper for SLICOT function SB02MT.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;

    if (nargin != 7)
    {
        print_usage ();
    }
    else
    {
        // SB02MT

        // arguments in
        char dico;
        char jobg = 'G';
        char jobl;
        char fact = 'N';
        char uplo = 'U';

        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix q = args(2).matrix_value ();
        Matrix r = args(3).matrix_value ();
        Matrqx l = args(4).matrix_value ();
      rF77_INT discrete l args(5).int_value ();
        F77_INT ijobl = args(6).int_value ();

        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            q.any_element_is_inf_or_nan () || r.any_element_is_inf_or_nan () ||
            l.any_element_is_inf_or_nan ())
          error ("__sl_are__: inputs must not contain NaN or Inf\n");

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

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldq = max (1, n);
        F77_INT ldr = max (1, m);
        F77_INT ldl = max (1, n);
        
        // arguments out
        F77_INT ldg = max (1, n);
        Matrix g (ldg, n);

        // unused output arguments
        OCTAVE_LOCAL_BUFFER (F77_INT, ipiv, m);
        F77_INT oufact;

        // workspace
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork_a, m);

        F77_INT ldwork_a = max (2, 3*m, n*m);
        OCTAVE_LOCAL_BUFFER (double, dwork_a, ldwork_a);


        // error indicator
        F77_INT info;


        // SLICOT routine SB02MT
        F77_XFCN (sb02mt, SB02MT,
                 (jobg, jobl,
                  fact, uplo,
                  n, m,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  q.fortran_vec (), ldq,
                  r.fortran_vec (), ldr,
                  l.fortran_vec (), ldl,
                  ipiv, oufact,
                  g.fortran_vec (), ldg,
                  iwork_a,
                  dwork_a, ldwork_a,
                  info));


        if (f77_exception_encountered)
            error ("are: __sl_are__: exception in SLICOT subroutine SB02MT");

        if (info != 0)
        {
            if (info < 0)
                error ("are: sb02mt: the %d-th argument had an illegal value", info);
            else if (info == m+1)
                error ("are: sb02mt: the matrix R is numerically singular");
            else
                error ("are: sb02mt: the %d-th element (1 <= %d <= %d) of the d factor is "
                       "exactly zero; the UdU' (or LdL') factorization has "
                       "been completed, but the block diagonal matrix d is "
                       "exactly singular", info, info, m);
        }
        
        
        // SB02RD
        
        // arguments in
        char job = 'A';
        char hinv = 'D';
        char trana = 'N';
        char scal = 'G';
        char sort = 'S';
        char lyapun = 'O';
        
        F77_INT ldt = max (1, n);
        F77_INT ldv = max (1, n);
        F77_INT ldx = max (1, n);
        F77_INT lds = max (1, 2*n);
        
        // arguments out
        Matrix x (ldx, n);

        double sep;
        double rcond;
        double ferr;
    
        ColumnVector wr (2*n);
        ColumnVector wi (2*n);
        
        // unused output arguments
        Matrix t (ldt, n);
        Matrix v (ldv, n);
        Matrix s (lds, 2*n);

        // workspace
        F77_INT liwork_b = max (2*n, n*n);
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork_b, liwork_b);

        F77_INT ldwork_b = 5 + max (1, 4*n*n + 8*n);
        OCTAVE_LOCAL_BUFFER (double, dwork_b, ldwork_b);

        OCTAVE_LOCAL_BUFFER (F77_LOGICAL, bwork_b, 2*n);

 
        // SLICOT routine SB02RD
        F77_XFCN (sb02rd, SB02RD,
                 (job, dico,
                  hinv, trana,
                  uplo, scal,
                  sort, fact,
                  lyapun,
                  n,
                  a.fortran_vec (), lda,
                  t.fortran_vec (), ldt,
                  v.fortran_vec (), ldv,
                  g.fortran_vec (), ldg,
                  q.fortran_vec (), ldq,
                  x.fortran_vec (), ldx,
                  sep,
                  rcond, ferr,
                  wr.fortran_vec (), wi.fortran_vec (),
                  s.fortran_vec (), lds,
                  iwork_b,
                  dwork_b, ldwork_b,
                  bwork_b,
                  info));


        static const char* err_msg[] = {
            "0: OK",
            "1: matrix A is (numerically) singular in discrete-"
                "time case",
            "2: the Hamiltonian or symplectic matrix H cannot be "
                "reduced to real Schur form",
            "3: the real Schur form of the Hamiltonian or "
                "symplectic matrix H cannot be appropriately ordered",
            "4: the Hamiltonian or symplectic matrix H has less "
                "than N stable eigenvalues",
            "5: if the N-th order system of linear algebraic "
                "equations, from which the solution matrix X would "
                "be obtained, is singular to working precision",
            "6: the QR algorithm failed to complete the reduction "
                "of the matrix Ac to Schur canonical form, T",
            "7: if T and -T' have some almost equal eigenvalues, if "
                "DICO = 'C', or T has almost reciprocal eigenvalues, "
                "if DICO = 'D'; perturbed values were used to solve "
                "Lyapunov equations, but the matrix T, if given (for "
                "FACT = 'F'), is unchanged. (This is a warning "
                "indicator.)"};

        error_msg ("are", info, 7, err_msg);

        // resize
        x.resize (n, n);
        wr.resize (n);
        wi.resize (n);

        // assemble complex vector - adapted from DEFUN complex in data.cc
        ComplexColumnVector pole (n, Complex ());

        for (F77_INT i = 0; i < n; i++)
            pole.xelem (i) = Complex (wr(i), wi(i));

        // return value
        retval(0) = x;
        retval(1) = pole;
        retval(2) = octave_value (ferr);
        retval(3) = octave_value (rcond);
        retval(4) = octave_value (sep);
    }

    return retval;
}
