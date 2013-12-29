/*

Copyright (C) 2013   Thomas Vasileiou

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

Orthogonal reduction of a descriptor system to a SVD-like coordinate form.
Uses SLICOT TG01FD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Thomas Vasileiou <thomas-v@wildmail.com>
Created: September 2013
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tg01fd, TG01FD)
                 (char& COMPQ, char& COMPZ, char& JOBA,
                  int& L, int& N, int& M, int& P,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* Q, int& LDQ,
                  double* Z, int& LDZ,
                  int& RANKE, int& RNKA22,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_tg01fd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg01fd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TG01FD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 6)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char compq;
        char compz;
        char joba = 'T';
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        const int qz_flag = args(4).int_value ();
        double tol = args(5).double_value ();
        
        if (qz_flag == 0)
        {
            compq = 'N';
            compz = 'N';
        }
        else
        {
            compq = 'I';
            compz = 'I';
        }

        int l = a.rows ();
        int n = l;
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs

        int lda = max (1, l);
        int lde = max (1, l);
        int ldb = max (1, l);
        int ldc = max (1, p);
        int ldq = max (1, l);
        int ldz = max (1, n);

        // arguments out
        Matrix q(l, l, 0.);
        Matrix z(n, n, 0.);
        Matrix empty(0, 0);
        int ranke, rnka22;
        
        // workspace
        int ldwork = max (1, n+p, min (l,n) + max (3*n-1, m, l));
        OCTAVE_LOCAL_BUFFER (int, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int info = 0;


        // SLICOT routine TG01FD
        F77_XFCN (tg01fd, TG01FD,
                 (compq, compz, joba,
                  l, n, m, p,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  q.fortran_vec (), ldq,
                  z.fortran_vec (), ldz,
                  ranke, rnka22,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_tg01fd__: exception in SLICOT subroutine TG01FD");
            
        if (info != 0)
            error ("__sl_tg01fd__: TG01FD returned info = %d", info);

        // return values
        retval(0) = a;
        retval(1) = e;
        retval(2) = b;
        retval(3) = c;
        retval(4) = octave_value (ranke);
        retval(5) = octave_value (rnka22);
        if (qz_flag == 0)
        {
            retval(6) = empty;
            retval(7) = empty;
        }
        else
        {
            retval(6) = q;
            retval(7) = z;
        }
    }
    
    return retval;
}
