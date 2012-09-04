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

Positive feedback controller for a discrete-time system (D == 0).
Uses SLICOT SB10KD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: July 2011
Version: 0.3

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb10kd, SB10KD)
                 (int& N, int& M, int& NP,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double& FACTOR,
                  double* AK, int& LDAK,
                  double* BK, int& LDBK,
                  double* CK, int& LDCK,
                  double* DK, int& LDDK,
                  double* RCOND,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  bool* BWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_sb10kd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb10kd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB10KD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 4)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        
        double factor = args(3).double_value ();
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int np = c.rows ();     // np: number of outputs
        
        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, np);
        
        int ldak = max (1, n);
        int ldbk = max (1, n);
        int ldck = max (1, m);
        int lddk = max (1, m);
        
        // arguments out
        Matrix ak (ldak, n);
        Matrix bk (ldbk, np);
        Matrix ck (ldck, n);
        Matrix dk (lddk, np);
        ColumnVector rcond (4);
        
        // workspace
        int liwork = 2 * max (n, np+m);
        int ldwork = 15*n*n + 6*n +
                     max (14*n+23, 16*n, 2*n+np+m, 3*(np+m)) +
                     max (n*n, 11*n*np + 2*m*m + 8*np*np + 8*m*n + 4*m*np + np);

        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (bool, bwork, 2*n);
        
        // error indicator
        int info;


        // SLICOT routine SB10KD
        F77_XFCN (sb10kd, SB10KD,
                 (n, m, np,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  factor,
                  ak.fortran_vec (), ldak,
                  bk.fortran_vec (), ldbk,
                  ck.fortran_vec (), ldck,
                  dk.fortran_vec (), lddk,
                  rcond.fortran_vec (),
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("ncfsyn: slsb10kd: exception in SLICOT subroutine SB10KD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the P-Riccati equation is not solved successfully",
            "2: the Q-Riccati equation is not solved successfully",
            "3: the X-Riccati equation is not solved successfully",
            "4: the iteration to compute eigenvalues failed to "
                "converge",
            "5: the matrix Rx + Bx'*X*Bx is singular",
            "6: the closed-loop system is unstable"};

        error_msg ("ncfsyn", info, 6, err_msg);


        // resizing
        ak.resize (n, n);
        bk.resize (n, np);
        ck.resize (m, n);
        dk.resize (m, np);
        
        // return values
        retval(0) = ak;
        retval(1) = bk;
        retval(2) = ck;
        retval(3) = dk;
        retval(4) = rcond;
    }
    
    return retval;
}
