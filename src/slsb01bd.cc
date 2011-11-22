/*

Copyright (C) 2009, 2010   Lukas F. Reichlin

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
Version: 0.4

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"

extern "C"
{ 
    int F77_FUNC (sb01bd, SB01BD)
                 (char& DICO,
                  int& N, int& M, int& NP,
                  double& ALPHA,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* WR, double* WI,
                  int& NFP, int& NAP, int& NUP,
                  double* F, int& LDF,
                  double* Z, int& LDZ,
                  double& TOL,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}
     
DEFUN_DLD (slsb01bd, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB01BD Release 5.0\n\
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
        // arguments in
        char dico;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        ColumnVector wr = args(2).column_vector_value ();
        ColumnVector wi = args(3).column_vector_value ();
        int discrete = args(4).int_value ();
        double alpha = args(5).double_value ();
        double tol = args(6).double_value ();
        
        if (discrete == 1)
            dico = 'D';
        else
            dico = 'C';
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int np = wr.rows ();
        
        int lda = max (1, a.rows ());
        int ldb = max (1, b.rows ());
        int ldf = max (1, m);
        int ldz = max (1, n);
        
        // arguments out
        int nfp;
        int nap;
        int nup;
        
        Matrix f (ldf, n);
        Matrix z (ldz, n);
        
        // workspace
        int ldwork = max (1, 5*m, 5*n, 2*n+4*m);
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int iwarn;
        int info;


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
            error ("place: slsb01bd: exception in SLICOT subroutine SB01BD");
            
        if (info != 0)
            error ("place: slsb01bd: SB01BD returned info = %d", info);

        if (info != 0)
        {
            if (info < 0)
                error ("place: slsb01bd: the %d-th argument had an invalid value", info);
            else
            {
                switch (info)
                {
                    case 1:
                        error ("place: 1: the reduction of A to a real Schur form failed.");
                    case 2:
                        error ("place: 2: a failure was detected during the ordering of the "
                               "real Schur form of A, or in the iterative process "
                               "for reordering the eigenvalues of Z'*(A + B*F)*Z "
                               "along the diagonal.");
                    case 3:
                        error ("place: 3: the number of eigenvalues to be assigned is less "
                               "than the number of possibly assignable eigenvalues; "
                               "NAP eigenvalues have been properly assigned, "
                               "but some assignable eigenvalues remain unmodified.");
                    case 4:
                        error ("place: 4: an attempt is made to place a complex conjugate "
                               "pair on the location of a real eigenvalue. This "
                               "situation can only appear when N-NFP is odd, "
                               "NP > N-NFP-NUP is even, and for the last real "
                               "eigenvalue to be modified there exists no available "
                               "real eigenvalue to be assigned. However, NAP "
                               "eigenvalues have been already properly assigned.");
                    default:
                        error ("place: unknown error, info = %d", info);
                }
            }
        }

        if (iwarn != 0)
            warning ("place: slsb01bd: %d violations of the numerical stability condition "
                     "NORM(F) <= 100*NORM(A)/NORM(B) occured during the "
                     "assignment of eigenvalues.", iwarn);

        // return values
        retval(0) = f;
        retval(1) = octave_value (nfp);
        retval(2) = octave_value (nap);
        retval(3) = octave_value (nup);
        retval(4) = z;
    }
    
    return retval;
}
