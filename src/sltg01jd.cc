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

Minimal realization of descriptor state-space models.
Uses SLICOT TG01JD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"

extern "C"
{ 
    int F77_FUNC (tg01jd, TG01JD)
                 (char& JOB, char& SYSTYP, char& EQUIL,
                  int& N, int& M, int& P,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  int& NR,
                  int* INFRED,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}
     
DEFUN_DLD (sltg01jd, args, nargout,
   "-*- texinfo -*-\n\
Slicot TG01JD Release 5.0\n\
No argument checking.\n\
For internal use only.")
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
        char job = 'I';
        char systyp = 'R';
        char equil = 'N';
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        double tol = args(4).double_value ();

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs

        int lda = max (1, n);
        int lde = max (1, n);
        int ldb = max (1, n);
        int ldc;

        if (n == 0)
            ldc = 1;
        else
            ldc = max (1, m, p);

        b.resize (ldb, max (m, p));

        // arguments out
        int nr;
        int infred[7];
        
        // workspace
        int liwork = n + max (m, p);
        int ldwork = max (n, 2*m, 2*p);
        // int ldwork = n * (2*n + m + p) + max (n, 2*m, 2*p);

        // FIXME: larger ldwork should give better results,
        //        but it breaks the test that Slicot provides.

        /*
        LDWORK  INTEGER
                The length of the array DWORK.
                LDWORK >= MAX(8*N,2*M,2*P), if EQUIL = 'S';
                LDWORK >= MAX(N,2*M,2*P),   if EQUIL = 'N'.
                If LDWORK >= MAX(2*N*N+N*M+N*P)+MAX(N,2*M,2*P) then more
                accurate results are to be expected by performing only
                those reductions phases (see METHOD), where effective
                order reduction occurs. This is achieved by saving the
                system matrices before each phase and restoring them if no
                order reduction took place.
        */

        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int info = 0;


        // SLICOT routine TG01JD
        F77_XFCN (tg01jd, TG01JD,
                 (job, systyp, equil,
                  n, m, p,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  nr,
                  infred,
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("dss: minreal: sltg01jd: exception in SLICOT subroutine TG01JD");
            
        if (info != 0)
            error ("dss: minreal: sltg01jd: TG01JD returned info = %d", info);

        // resize
        a.resize (nr, nr);
        e.resize (nr, nr);
        b.resize (nr, m);
        c.resize (p, nr);

        // return values
        retval(0) = a;
        retval(1) = e;
        retval(2) = b;
        retval(3) = c;
        retval(4) = octave_value (nr);
    }
    
    return retval;
}