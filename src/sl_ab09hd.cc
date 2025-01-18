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

Model reduction based on balanced stochastic truncation method.
Uses SLICOT AB09HD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab09hd, AB09HD)
                 (char& DICO, char& JOB, char& EQUIL, char& ORDSEL,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  F77_INT& NR,
                  double& ALPHA, double& BETA,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  F77_INT& NS,
                  double* HSV,
                  double& TOL1, double& TOL2,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_LOGICAL* BWORK,
                  F77_INT& IWARN, F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ab09hd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab09hd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_ab09hd__ (@dots{})\n"
   "Wrapper for SLICOT function AB09HD.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 13)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char job;
        char equil;
        char ordsel;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();

        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan () || d.any_element_is_inf_or_nan ())
          error ("__sl_ab09hd__: inputs must not contain NaN or Inf\n");
        
        const F77_INT idico = args(4).int_value ();
        const F77_INT iequil = args(5).int_value ();
        const F77_INT ijob = args(6).int_value ();
        
        F77_INT nr = args(7).int_value ();
        const F77_INT iordsel = args(8).int_value ();
        
        double alpha = args(9).double_value ();
        double beta = args(10).double_value ();
        double tol1 = args(11).double_value ();
        double tol2 = args(12).double_value ();

        switch (ijob)
        {
            case 0:
                job = 'B';
                break;
            case 1:
                job = 'F';
                break;
            case 2:
                job = 'S';
                break;
            case 3:
                job = 'P';
                break;
            default:
                error ("__sl_ab09hd__: argument job invalid");
        }

        if (idico == 0)
            dico = 'C';
        else
            dico = 'D';

        if (iequil == 0)
            equil = 'S';
        else
            equil = 'N';

        if (iordsel == 0)
            ordsel = 'F';
        else
            ordsel = 'A';

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldd = max (1, p);

        // arguments out
        F77_INT ns;
        ColumnVector hsv (n);

        // workspace
        F77_INT liwork = max (1, 2*n);
        F77_INT mb;
        
        if (beta == 0)
            mb = m;
        else
            mb = m + p;
        
        F77_INT ldwork = 2*n*n + mb*(n+p)
                     + max (2, n*(max (n,mb,p)+5), 2*n*p + max (p*(mb+2), 10*n*(n+1)));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (F77_LOGICAL, bwork, 2*n);
        
        // error indicators
        F77_INT iwarn = 0;
        F77_INT info = 0;


        // SLICOT routine AB09HD
        F77_XFCN (ab09hd, AB09HD,
                 (dico, job, equil, ordsel,
                  n, m, p,
                  nr,
                  alpha, beta,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  ns,
                  hsv.fortran_vec (),
                  tol1, tol2,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  iwarn, info));

        if (f77_exception_encountered)
            error ("bstmodred: exception in SLICOT subroutine AB09HD");


        static const char* err_msg[] = {
            "0: OK",
            "1: the computation of the ordered real Schur form of A "
                "failed",
            "2: the reduction of the Hamiltonian matrix to real "
                "Schur form failed",
            "3: the reordering of the real Schur form of the "
                "Hamiltonian matrix failed",
            "4: the Hamiltonian matrix has less than N stable "
                "eigenvalues",
            "5: the coefficient matrix U11 in the linear system "
                "X*U11 = U21 to determine X is singular to working "
                "precision",
            "6: BETA = 0 and D has not a maximal row rank",
            "7: the computation of Hankel singular values failed",
            "8: the separation of the ALPHA-stable/unstable diagonal "
                "blocks failed because of very close eigenvalues",
            "9: the resulting order of reduced stable part is less "
                "than the number of unstable zeros of the stable "
                "part"}; 

        static const char* warn_msg[] = {
            "0: OK",
            "1: with ORDSEL = 'F', the selected order NR is greater "
                "than NSMIN, the sum of the order of the "
                "ALPHA-unstable part and the order of a minimal "
                "realization of the ALPHA-stable part of the given "
                "system; in this case, the resulting NR is set equal "
                "to NSMIN.",
            "2: with ORDSEL = 'F', the selected order NR corresponds "
                "to repeated singular values for the ALPHA-stable "
                "part, which are neither all included nor all "
                "excluded from the reduced model; in this case, the "
                "resulting NR is automatically decreased to exclude "
                "all repeated singular values.",
            "3: with ORDSEL = 'F', the selected order NR is less "
                "than the order of the ALPHA-unstable part of the "
                "given system; in this case NR is set equal to the "
                "order of the ALPHA-unstable part."};
 
        error_msg ("bstmodred", info, 9, err_msg);
        warning_msg ("bstmodred", iwarn, 3, warn_msg);

        // resize
        a.resize (nr, nr);
        b.resize (nr, m);
        c.resize (p, nr);
        hsv.resize (ns);
        
        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = d;
        retval(4) = octave_value (nr);
        retval(5) = hsv;
        retval(6) = octave_value (ns);
    }
    
    return retval;
}
