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

Compute initial state vector x0
Uses IB01CD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: May 2012
Version: 0.1

*/

#include <octave/oct.h>
#include <octave/f77-fcn.h>
#include <octave/Cell.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ib01cd, IB01CD)
                 (char& JOBX0, char& COMUSE, char& JOB,
                  int& N, int& M, int& L,
                  int& NSMP,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* U, int& LDU,
                  double* Y, int& LDY,
                  double* X0,
                  double* V, int& LDV,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}

// PKG_ADD: autoload ("slib01cd", "devel_slicot_functions.oct");
DEFUN_DLD (slib01cd, args, nargout,
   "-*- texinfo -*-\n\
Slicot IB01CD Release 5.0\n\
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
        char jobx0 = 'X';
        char comuse = 'U';
        char jobbd = 'D';

        const Cell y_cell = args(0).cell_value ();
        const Cell u_cell = args(1).cell_value ();
        
        Matrix a = args(2).matrix_value ();
        Matrix b = args(3).matrix_value ();
        Matrix c = args(4).matrix_value ();
        Matrix d = args(5).matrix_value ();
        
        double rcond = args(6).double_value ();
        double tol_c = rcond;
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int l = c.rows ();      // l: number of outputs
        
        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, l);
        int ldd = max (1, l);

        // m and l are equal for all experiments, checked by iddata class
        int n_exp = y_cell.nelem ();            // number of experiments


        // arguments out
        Cell x0_cell (n_exp, 1);    // cell of initial state vectors x0

        // repeat for every experiment in the dataset
        // compute individual initial state vector x0 for every experiment        
        for (int i = 0; i < n_exp; i++)
        {
            Matrix y = y_cell.elem(i).matrix_value ();
            Matrix u = u_cell.elem(i).matrix_value ();
            
            int nsmp = y.rows ();   // nsmp: number of samples
            int ldv = max (1, n);
            
            int ldu;
        
            if (m == 0)
                ldu = 1;
            else                    // m > 0
                ldu = nsmp;

            int ldy = nsmp;

            // arguments out
            ColumnVector x0 (n);
            Matrix v (ldv, n);

            // workspace
            int liwork_c = n;     // if  JOBX0 = 'X'  and  COMUSE <> 'C'
            int ldwork_c;
            int t = nsmp;
   
            int ldw1_c = 2;
            int ldw2_c = t*l*(n + 1) + 2*n + max (2*n*n, 4*n);
            int ldw3_c = n*(n + 1) + 2*n + max (n*l*(n + 1) + 2*n*n + l*n, 4*n);

            ldwork_c = ldw1_c + n*( n + m + l ) + max (5*n, ldw1_c, min (ldw2_c, ldw3_c));

            OCTAVE_LOCAL_BUFFER (int, iwork_c, liwork_c);
            OCTAVE_LOCAL_BUFFER (double, dwork_c, ldwork_c);

            // error indicators
            int iwarn_c = 0;
            int info_c = 0;

            // SLICOT routine IB01CD
            F77_XFCN (ib01cd, IB01CD,
                     (jobx0, comuse, jobbd,
                      n, m, l,
                      nsmp,
                      a.fortran_vec (), lda,
                      b.fortran_vec (), ldb,
                      c.fortran_vec (), ldc,
                      d.fortran_vec (), ldd,
                      u.fortran_vec (), ldu,
                      y.fortran_vec (), ldy,
                      x0.fortran_vec (),
                      v.fortran_vec (), ldv,
                      tol_c,
                      iwork_c,
                      dwork_c, ldwork_c,
                      iwarn_c, info_c));


            if (f77_exception_encountered)
                error ("slib01cd: exception in SLICOT subroutine IB01CD");

            static const char* err_msg_c[] = {
                "0: OK",
                "1: the QR algorithm failed to compute all the "
                    "eigenvalues of the matrix A (see LAPACK Library "
                    "routine DGEES); the locations  DWORK(i),  for "
                    "i = g+1:g+N*N,  contain the partially converged "
                    "Schur form",
                "2: the singular value decomposition (SVD) algorithm did "
                    "not converge"};

            static const char* warn_msg_c[] = {
                "0: OK",
                "1: warning message not specified",
                "2: warning message not specified",
                "3: warning message not specified",
                "4: the least squares problem to be solved has a "
                    "rank-deficient coefficient matrix",
                "5: warning message not specified",
                "6: the matrix  A  is unstable;  the estimated  x(0) "
                    "and/or  B and D  could be inaccurate"};


            error_msg ("slib01cd", info_c, 2, err_msg_c);
            warning_msg ("slib01cd", iwarn_c, 6, warn_msg_c);
            
            x0_cell.elem(i) = x0;       // add x0 from the current experiment to cell of initial state vectors
        }
   
        
        // return values
        retval(0) = x0_cell;
    }
    
    return retval;
}
