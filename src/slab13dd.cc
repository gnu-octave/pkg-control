/*

Copyright (C) 2009 - 2010   Lukas F. Reichlin

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

L-infinity norm of a SS model.
Uses SLICOT AB13DD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include <complex>

extern "C"
{ 
    int F77_FUNC (ab13dd, AB13DD)
                 (char& DICO, char& JOBE,
                  char& EQUIL, char& JOBD,
                  int& N, int& M, int& P,
                  double* FPEAK,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* GPEAK,
                  double& TOL,
                  int* IWORK, double* DWORK, int& LDWORK,
                  Complex* CWORK, int& LCWORK,
                  int& INFO);
}

int max (int a, int b)
{
    if (a > b)
        return a;
    else
        return b;
}

int min (int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}
     
DEFUN_DLD (slab13dd, args, nargout, "Slicot AB13DD Release 5.0")
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
        char dico;
        char jobe = 'I';
        char equil = 'N';
        char jobd = 'D';
        
        NDArray a = args(0).array_value ();
        NDArray b = args(1).array_value ();
        NDArray c = args(2).array_value ();
        NDArray d = args(3).array_value ();
        double* e = 0;
        int digital = args(4).int_value ();
        double tol = args(5).double_value ();
        
        if (digital == 0)
            dico = 'C';
        else
            dico = 'D';
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs
        
        int lda = max (1, a.rows ());
        int ldb = max (1, b.rows ());
        int ldc = max (1, c.rows ());
        int ldd = max (1, d.rows ());
        int lde = 1;
        
        dim_vector dv (1);
        dv(0) = 2;
        NDArray fpeak (dv);
        NDArray gpeak (dv);
        
        fpeak(0) = 0;
        fpeak(1) = 1;
        
        // workspace
        int ldwork = max (1, 15*n*n + p*p + m*m + (6*n+3)*(p+m) + 4*p*m +
                          n*m + 22*n + 7*min(p,m));
        int lcwork = max (1, (n+m)*(n+p) + 2*min(p,m) + max(p,m));
        
        OCTAVE_LOCAL_BUFFER (int, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (Complex, cwork, lcwork);
        
        // error indicator
        int info;


        // SLICOT routine AB13DD
        F77_XFCN (ab13dd, AB13DD,
                 (dico, jobe,
                  equil, jobd,
                  n, m, p,
                  fpeak.fortran_vec (),
                  a.fortran_vec (), lda,
                  e, lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  gpeak.fortran_vec (),
                  tol,
                  iwork, dwork, ldwork,
                  cwork, lcwork,
                  info));

        if (f77_exception_encountered)
            error ("lti: norm: slab13dd: exception in SLICOT subroutine AB13DD");
            
        if (info != 0)
            error ("lti: norm: slab13dd: AB13DD returned info = %d", info);
        
        // return values
        retval(0) = fpeak;
        retval(1) = gpeak;
    }
    
    return retval;
}
