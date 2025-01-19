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

Minimal realization of descriptor state-space models.
Uses SLICOT TG01JD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.5

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tg01jd, TG01JD)
                 (char& JOB, char& SYSTYP, char& EQUIL,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  double* A, F77_INT& LDA,
                  double* E, F77_INT& LDE,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  F77_INT& NR,
                  F77_INT* INFRED,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_tg01jd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tg01jd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_tg01jd__ (@dots{})\n"
   "Wrapper for SLICOT function TG01JD.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;

    if (nargin != 8)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char job;
        char systyp;
        char equil;
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value ();
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        double tol = args(4).double_value ();
        const F77_INT scaled = args(5).int_value ();
        const F77_INT ijob = args(6).int_value ();
        const F77_INT isystyp = args(7).int_value ();
        
        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan () || e.any_element_is_inf_or_nan ())
        error ("__sl_tg01jd__: inputs must not contain NaN or Inf\n");

        F77_INT c_factor = 1;

        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';
        
        switch (ijob)
        {
            case 0:
                job = 'I';
                c_factor = 2;
                break;
            case 1:
                job = 'C';
                break;
            case 2:
                job = 'O';
                break;
            default:
                error ("__sl_tg01jd__: argument job invalid");
        }

        switch (isystyp)
        {
            case 0:
                systyp = 'R';
                c_factor = 2;
                break;
            case 1:
                systyp = 'S';
                break;
            case 2:
                systyp = 'P';
                break;
            default:
                error ("__sl_tg01jd__: argument systyp invalid");
        }

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT lde = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc;

        if (n == 0)
            ldc = 1;
        else
            ldc = max (1, m, p);

        a.resize (lda, n);
        e.resize (lde, n);
        
        if (job == 'C')
            b.resize (ldb, m);
        else
            b.resize (ldb, max (m, p));

        c.resize (ldc, n);

        // arguments out
        F77_INT nr;
        F77_INT infred[7];
        
        // workspace
        F77_INT liwork = c_factor*n + max (m, p);
        F77_INT ldwork;
        // F77_INT ldwork = max (n, 2*m, 2*p);
        // F77_INT ldwork = n * (2*n + m + p) + max (n, 2*m, 2*p);
        
        if (equil == 'S')
            ldwork = max (8*n, 2*m, 2*p);
        else        // if EQUIL = 'N'
            ldwork = max (n, 2*m, 2*p);

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

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT info = 0;


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
            error ("dss: minreal: __sl_tg01jd__: exception in SLICOT subroutine TG01JD");
            
        if (info != 0)
            error ("dss: minreal: __sl_tg01jd__: TG01JD returned info = %d", static_cast<int> (info));

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
