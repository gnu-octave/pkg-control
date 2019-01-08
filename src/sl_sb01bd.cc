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

Pole assignment for a given matrix pair (A,B).
Uses SLICOT SB01BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.6

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb01bd, SB01BD)
                 (char& DICO,
                  F77_INT& N, F77_INT& M, F77_INT& NP,
                  double& ALPHA,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* WR, double* WI,
                  F77_INT& NFP, F77_INT& NAP, F77_INT& NUP,
                  double* F, F77_INT& LDF,
                  double* Z, F77_INT& LDZ,
                  double& TOL,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& IWARN, F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_sb01bd__", "__control_slicot_functions__.oct");     
DEFUN_DLD (__sl_sb01bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB01BD Release 5.0\n\
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
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        ColumnVector wr = args(2).column_vector_value ();
        ColumnVector wi = args(3).column_vector_value ();
        F77_INT discrete = args(4).int_value ();
        double alpha = args(5).double_value ();
        double tol = args(6).double_value ();
        
        if (discrete == 1)
            dico = 'D';
        else
            dico = 'C';
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT np = TO_F77_INT (wr.rows ());
        
        F77_INT lda = max (1, TO_F77_INT (a.rows ()));
        F77_INT ldb = max (1, TO_F77_INT (b.rows ()));
        F77_INT ldf = max (1, m);
        F77_INT ldz = max (1, n);
        
        // arguments out
        F77_INT nfp;
        F77_INT nap;
        F77_INT nup;
        
        Matrix f (ldf, n);
        Matrix z (ldz, n);
        
        // workspace
        F77_INT ldwork = max (1, 5*m, 5*n, 2*n+4*m);
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT iwarn;
        F77_INT info;


        // SLICOT routine SB01BD
        F77_XFCN (sb01bd, SB01BD,
                 (dico,
                  n, m, np,
                  alpha,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  wr.fortran_vec (), wi.fortran_vec (),
                  nfp, nap, nup,
                  f.fortran_vec (), ldf,
                  z.fortran_vec (), ldz,
                  tol,
                  dwork, ldwork,
                  iwarn, info));

        if (f77_exception_encountered)
            error ("place: __sl_sb01bd__: exception in SLICOT subroutine SB01BD");
            
        static const char* err_msg[] = {
            "0: OK",
            "1: the reduction of A to a real Schur form failed.",
            "2: a failure was detected during the ordering of the "
                "real Schur form of A, or in the iterative process "
                "for reordering the eigenvalues of Z'*(A + B*F)*Z "
                "along the diagonal.",
            "3: the number of eigenvalues to be assigned is less "
                "than the number of possibly assignable eigenvalues; "
                "NAP eigenvalues have been properly assigned, "
                "but some assignable eigenvalues remain unmodified.",
            "4: an attempt is made to place a complex conjugate "
                "pair on the location of a real eigenvalue. This "
                "situation can only appear when N-NFP is odd, "
                "NP > N-NFP-NUP is even, and for the last real "
                "eigenvalue to be modified there exists no available "
                "real eigenvalue to be assigned. However, NAP "
                "eigenvalues have been already properly assigned."};

        static const char* warn_msg[] = {
            "0: OK",
/* 0+%d: %d */  "violations of the numerical stability condition "
                "NORM(F) <= 100*NORM(A)/NORM(B) occurred during the "
                "assignment of eigenvalues."};

        error_msg ("place", info, 4, err_msg);
        warning_msg ("place", iwarn, 0, warn_msg, 0);


        // return values
        retval(0) = f;
        retval(1) = octave_value (nfp);
        retval(2) = octave_value (nap);
        retval(3) = octave_value (nup);
        retval(4) = z;
    }
    
    return retval;
}
