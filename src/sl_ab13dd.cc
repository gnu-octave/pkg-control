/*

Copyright (C) 2009-2016   Lukas F. Reichlin

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

L-infinity norm of a SS model.
Uses SLICOT AB13DD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.5

*/

#include <octave/oct.h>
#include <complex>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab13dd, AB13DD)
                 (char& DICO, char& JOBE,
                  char& EQUIL, char& JOBD,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double* FPEAK,
                  double* A, F77_INT& LDA,
                  double* E, F77_INT& LDE,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  double* GPEAK,
                  double& TOL,
                  F77_INT* IWORK, double* DWORK, F77_INT& LDWORK,
                  Complex* CWORK, F77_INT& LCWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ab13dd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab13dd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB13DD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 9)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char jobe;
        char equil;
        char jobd = 'D';
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        Matrix d = args(4).matrix_value ();
        F77_INT discrete = args(5).int_value ();
        F77_INT descriptor = args(6).int_value ();
        double tol = args(7).double_value ();
        const F77_INT scaled = args(8).int_value ();
        
        if (discrete == 0)
            dico = 'C';
        else
            dico = 'D';

        if (descriptor == 0)
            jobe = 'I';
        else
            jobe = 'G';

        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        
        F77_INT lda = max (1, n);
        F77_INT lde = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldd = max (1, p);
        
        ColumnVector fpeak (2);
        ColumnVector gpeak (2);
        
        fpeak(0) = 0;
        fpeak(1) = 1;
        
        // workspace
        F77_INT ldwork = max (1, 15*n*n + p*p + m*m + (6*n+3)*(p+m) + 4*p*m +
                          n*m + 22*n + 7*min(p,m));
        F77_INT lcwork = max (1, (n+m)*(n+p) + 2*min(p,m) + max(p,m));
        
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        OCTAVE_LOCAL_BUFFER (Complex, cwork, lcwork);
        
        // error indicator
        F77_INT info;


        // SLICOT routine AB13DD
        F77_XFCN (ab13dd, AB13DD,
                 (dico, jobe,
                  equil, jobd,
                  n, m, p,
                  fpeak.fortran_vec (),
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  gpeak.fortran_vec (),
                  tol,
                  iwork, dwork, ldwork,
                  cwork, lcwork,
                  info));

        if (f77_exception_encountered)
            error ("lti: norm: __sl_ab13dd__: exception in SLICOT subroutine AB13DD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the matrix E is (numerically) singular",
            "2: the (periodic) QR (or QZ) algorithm for computing "
                "eigenvalues did not converge",
            "3: the SVD algorithm for computing singular values did "
                "not converge",
            "4: the tolerance is too small and the algorithm did "
                "not converge"};

        error_msg ("__sl_ab13dd__", info, 4, err_msg);

        
        // return values
        retval(0) = fpeak;
        retval(1) = gpeak;
    }
    
    return retval;
}
