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

Balance descriptor state-space model.
Uses SLICOT TG01AD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: June 2011
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"

extern "C"
{ 
    int F77_FUNC (tg01ad, TG01AD)
                 (char& JOB,
                  int& L, int& N, int& M, int& P,
                  double& TRESH,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* LSCALE, double *RSCALE,
                  double* DWORK,
                  int& INFO);
}
     
DEFUN_DLD (sltg01ad, args, nargout,
   "-*- texinfo -*-\n\
Slicot TG01AD Release 5.0\n\
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
        char job = 'A';
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        double maxred = args(3).double_value ();

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        

        // arguments out
        ColumnVector scale (n);

        
        // error indicators
        int info = 0;


        // SLICOT routine TG01AD
        F77_XFCN (tg01ad, TG01AD,
                 (job,
                  n, m, p,
                  maxred,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  scale.fortran_vec (),
                  info));

        if (f77_exception_encountered)
            error ("ss: prescale: sltg01ad: exception in SLICOT subroutine TG01AD");
            
        if (info != 0)
            error ("ss: prescale: sltg01ad: TG01AD returned info = %d", info);


        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = octave_value (maxred);
        retval(4) = scale;
    }
    
    return retval;
}
