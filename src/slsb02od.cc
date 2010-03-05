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
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>

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

int max (int a, int b)
{
    if (a > b)
        return a;
    else
        return b;
}

int max (int a, int b, int c)
{
    int d = max (a, b);
    
    return max (c, d);
}

int max (int a, int b, int c, int d)
{
    int e = max (a, b);
    int f = max (c, d);
    
    return max (e, f);
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
        
        NDArray a = args(0).array_value ();
        NDArray b = args(1).array_value ();
        NDArray q = args(2).array_value ();
        NDArray r = args(3).array_value ();
        NDArray l = args(4).array_value ();
        int digital = args(5).int_value ();
        int ijobl = args(6).int_value ();

        if (digital == 0)
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
        dim_vector dv (2);
        dv(0) = ldx;
        dv(1) = n;
        NDArray x (dv);
        
        // unused output arguments
        OCTAVE_LOCAL_BUFFER (double, alfar, 2*n);
        OCTAVE_LOCAL_BUFFER (double, alfai, 2*n);
        OCTAVE_LOCAL_BUFFER (double, beta, 2*n);
        
        int lds = max (1, 2*n + m);
        OCTAVE_LOCAL_BUFFER (double, s, lds*lds);
        
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
                  alfar, alfai,
                  beta,
                  s, lds,
                  t, ldt,
                  u, ldu,
                  tol,
                  iwork,
                  dwork, ldwork,
                  bwork,
                  info));

        if (f77_exception_encountered)
            error ("lti: norm: slsb02od: exception in SLICOT subroutine SB02OD");
            
        if (info != 0)
            error ("lti: norm: slsb02od: SB02OD returned info = %d", info);
        
        // return value
        retval(0) = x;
    }
    
    return retval;
}
