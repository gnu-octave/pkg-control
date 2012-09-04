/*

Copyright (C) 2010   Lukas F. Reichlin

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

Square-root solver for generalized Lyapunov equations.
Uses SLICOT SG03BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sg03bd, SG03BD)
                 (char& DICO, char& FACT, char& TRANS,
                  int& N, int& M,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* Q, int& LDQ,
                  double* Z, int& LDZ,
                  double* B, int& LDB,
                  double& SCALE,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_sg03bd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sg03bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SG03BD Release 5.0\n\
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
        char dico;
        char fact = 'N';
        char trans = 'N';
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        int discrete = args(3).int_value ();
        
        if (discrete == 0)
          dico = 'C';
        else
          dico = 'D';
        
        int n = a.rows ();
        int m = b.rows ();
        
        int lda = max (1, n);
        int lde = max (1, n);
        int ldq = max (1, n);
        int ldz = max (1, n);
        int ldb = max (1, m, n);

        int n1 = max (m, n);
        b.resize (ldb, n1);
                
        // arguments out
        double scale;

        Matrix q (ldq, n);
        Matrix z (ldz, n);
        ColumnVector alphar (n);
        ColumnVector alphai (n);
        ColumnVector beta (n);
        
        // workspace
        int ldwork = max (1, 4*n, 6*n-6);
        
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        int info;
        

        // SLICOT routine SG03BD
        F77_XFCN (sg03bd, SG03BD,
                 (dico, fact, trans,
                  n, m,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  q.fortran_vec (), ldq,
                  z.fortran_vec (), ldz,
                  b.fortran_vec (), ldb,
                  scale,
                  alphar.fortran_vec (), alphai.fortran_vec (),
                  beta.fortran_vec (),
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("lyap: __sl_sg03bd__: exception in SLICOT subroutine SG03BD");

        if (info != 0)
            error ("lyap: __sl_sg03bd__: SG03BD returned info = %d", info);

        // resize
        b.resize (n, n);
        
        // return values
        retval(0) = b;  // b has been overwritten by cholesky factor u
        retval(1) = octave_value (scale);
    }
    
    return retval;
}
