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

Discrete-time <--> continuous-time systems conversion 
by a bilinear transformation.
Uses SLICOT AB04MD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2011
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab04md, AB04MD)
                 (char& TYPE,
                  int& N, int& M, int& P,
                  double& ALPHA, double& BETA,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_ab04md__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab04md__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB04MD Release 5.0\n\
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
        char type;

        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        
        double alpha = args(4).double_value ();
        double beta = args(5).double_value ();
        int discrete = args(6).int_value ();

        if (discrete == 0)
            type = 'C';
        else
            type = 'D';
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs
        
        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        int ldd = max (1, p);
        
        // workspace
        int ldwork = max (1, n);

        OCTAVE_LOCAL_BUFFER (int, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        int info;


        // SLICOT routine AB04MD
        F77_XFCN (ab04md, AB04MD,
                 (type,
                  n, m, p,
                  alpha, beta,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("__sl_ab04md__: exception in SLICOT subroutine AB04MD");

        if (info != 0)
            error ("__sl_ab04md__: AB04MD returned info = %d", info);
        
        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = d;
    }
    
    return retval;
}
