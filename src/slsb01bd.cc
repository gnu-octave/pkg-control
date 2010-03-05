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

Pole assignment for a given matrix pair (A,B).
Uses SLICOT SB01BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.2.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>

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

int max (int a, int b)
{
    if (a > b)
        return a;
    else
        return b;
}

int max (int a, int b, int c, int d)
{
    int e = max (a, b);
    int f = max (c, d);
    
    return max (e, f);
}
     
DEFUN_DLD (slsb01bd, args, nargout, "Slicot SB01BD Release 5.0")
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
        
        NDArray a = args(0).array_value ();
        NDArray b = args(1).array_value ();
        NDArray wr = args(2).array_value ();
        NDArray wi = args(3).array_value ();
        int digital = args(4).int_value ();
        double alpha = args(5).double_value ();
        double tol = args(6).double_value ();
        
        if (digital == 1)
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
        
        dim_vector dv (2);
        dv(0) = ldf;
        dv(1) = n;
        NDArray f (dv);
        
        OCTAVE_LOCAL_BUFFER (double, z, ldz*n);
        
        // workspace
        int ldwork = max (1, 5*m, 5*n, 2*n+4*m);
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int iwarn;
        int info;


        // SLICOT routine AB13DD
        F77_XFCN (sb01bd, SB01BD,
                 (dico,
                  n, m, np,
                  alpha,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  wr.fortran_vec (), wi.fortran_vec (),
                  nfp, nap, nup,
                  f.fortran_vec (), ldf,
                  z, ldz,
                  tol,
                  dwork, ldwork,
                  iwarn, info));

        if (f77_exception_encountered)
            error ("place: slsb01bd: exception in SLICOT subroutine SB01BD");
            
        if (info != 0)
            error ("place: slsb01bd: SB01BD returned info = %d", info);
        
        // return values
        retval(0) = f;
        retval(1) = octave_value (iwarn);
        retval(2) = octave_value (nfp);
        retval(3) = octave_value (nap);
        retval(4) = octave_value (nup);
    }
    
    return retval;
}
