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

Hankel singular values.
Uses SLICOT AB13AD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: January 2010
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"

extern "C"
{ 
    int F77_FUNC (ab13ad, AB13AD)
                 (char& DICO, char& EQUIL,
                  int& N, int& M, int& P,
                  double& ALPHA,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  int& NS,
                  double* HSV,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}
     
DEFUN_DLD (slab13ad, args, nargout, "Slicot AB13AD Release 5.0")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 5)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char equil = 'N';
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        int digital = args(3).int_value ();
        double alpha = args(4).double_value ();
        
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
        
        // arguments out
        int ns;
        
        ColumnVector hsv (n);
        
        // workspace
        int ldwork = max (1, n*(max (n, m, p) + 5) + n*(n+1)/2);
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int info;


        // SLICOT routine AB13DD
        F77_XFCN (ab13ad, AB13AD,
                 (dico, equil,
                  n, m, p,
                  alpha,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  ns,
                  hsv.fortran_vec (),
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("place: slab13ad: exception in SLICOT subroutine AB13AD");
            
        if (info != 0)
            error ("place: slab13ad: AB13AD returned info = %d", info);
        
        // return values
        retval(0) = hsv;
        retval(1) = octave_value (ns);
    }
    
    return retval;
}
