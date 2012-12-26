/*

Copyright (C) 2012   Lukas F. Reichlin

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
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"
#include <complex>

extern "C"
{ 
    int F77_FUNC (sb02mt, SB02MT)
                 (char& JOBG, char& JOBL,
                  char& FACT, char& UPLO,
                  int& N, int& M,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* Q, int& LDQ,
                  double* R, int& LDR,
                  double* L, int& LDL,
                  int* IPIV, int& OUFACT,
                  double* G, int& LDG,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);

                  

    int F77_FUNC (sb02rd, SB02RD)
                 (char& JOB, char& DICO,
                  char& HINV, char& TRANA,
                  char& UPLO, char& SCAL,
                  char& SORT, char& FACT,
                  char& LYAPUN,
                  int& N,
                  double* A, int& LDA,
                  double* T, int& LDT,
                  double* V, int& LDV,
                  double* G, int& LDG,
                  double* Q, int& LDQ,
                  double* X, int& LDX,
                  double& SEP,
                  double& RCOND, double& FERR
                  double* WR, double* WI,
                  double* S, int& LDS,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  bool* BWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_are__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_are__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB02RD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;

    if (nargin != 7)
    {
        print_usage ();
    }
    else
    {
        // SB02MT

        // arguments in
        char jobg = 'G';
        char jobl;
        char fact = 'N';
        char uplo = 'U';

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

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldq = max (1, n);
        int ldr = max (1, m);
        int ldl = max (1, n);
        
        // arguments out
        int ldg = max (1, n);
        Matrix g (ldg, n);

        // unused output arguments
        OCTAVE_LOCAL_BUFFER (int, ipiv, m);
        int oufact;

        // workspace
        OCTAVE_LOCAL_BUFFER (int, iwork_a, m);

        int ldwork_a = max (2, 3*m, n*m);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork_a);


        // error indicator
        int info;


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
                  iwork,
                  dwork, ldwork,
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
        char fact = 'N';
        char lyapun = 'O';
        
        int ldt = max (1, n);
        int ldv = max (1, n);
        
        // arguments out


        // unused output arguments
        Matrix t (ldt, n);
        Matrix v (ldv, n);

 
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
                  iwork,
                  dwork, ldwork,
                  bwork,
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


        // assemble complex vector - adapted from DEFUN complex in data.cc
        alfar.resize (n);
        alfai.resize (n);
        beta.resize (n);

        ColumnVector poler (n);
        ColumnVector polei (n);

        poler = quotient (alfar, beta);
        polei = quotient (alfai, beta);
        
        ComplexColumnVector pole (n, Complex ());

        for (octave_idx_type i = 0; i < n; i++)
            pole.xelem (i) = Complex (poler(i), polei(i));

        // return value
        retval(0) = x;
        retval(1) = pole;
    }

    return retval;
}
