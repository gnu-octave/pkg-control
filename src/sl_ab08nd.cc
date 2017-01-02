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

Invariant zeros of state-space models.
Uses SLICOT AB08ND by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.8

*/

#include <octave/oct.h>
#include "common.h"
#include <complex>
#include <xpow.h>

extern "C"
{ 
    int F77_FUNC (ab08nd, AB08ND)
                 (char& EQUIL,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  const double* A, F77_INT& LDA,
                  const double* B, F77_INT& LDB,
                  const double* C, F77_INT& LDC,
                  const double* D, F77_INT& LDD,
                  F77_INT& NU, F77_INT& RANK, F77_INT& DINFZ,
                  F77_INT& NKROR, F77_INT& NKROL, F77_INT* INFZ,
                  F77_INT* KRONR, F77_INT* KRONL,
                  double* AF, F77_INT& LDAF,
                  double* BF, F77_INT& LDBF,
                  double& TOL,
                  F77_INT* IWORK, double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
                                   
    int F77_FUNC (dggev, DGGEV)
                 (char& JOBVL, char& JOBVR,
                  F77_INT& N,
                  double* AF, F77_INT& LDAF,
                  double* BF, F77_INT& LDBF,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  double* VL, F77_INT& LDVL,
                  double* VR, F77_INT& LDVR,
                  double* WORK, F77_INT& LWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ab08nd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab08nd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB08ND Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 5)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char equil;
        
        const Matrix a = args(0).matrix_value ();
        const Matrix b = args(1).matrix_value ();
        const Matrix c = args(2).matrix_value ();
        const Matrix d = args(3).matrix_value ();
        const F77_INT scaled = args(4).int_value ();
        
        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';
        
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        
        F77_INT lda = max (1, TO_F77_INT (a.rows ()));
        F77_INT ldb = max (1, TO_F77_INT (b.rows ()));
        F77_INT ldc = max (1, TO_F77_INT (c.rows ()));
        F77_INT ldd = max (1, TO_F77_INT (d.rows ()));
        
        // arguments out
        F77_INT nu;
        F77_INT rank;
        F77_INT dinfz;
        F77_INT nkror;
        F77_INT nkrol;
        
        F77_INT ldaf = max (1, n + m);
        F77_INT ldbf = max (1, n + p);

        OCTAVE_LOCAL_BUFFER (F77_INT, infz, n);
        OCTAVE_LOCAL_BUFFER (F77_INT, kronr, 1 + max (n, m));
        OCTAVE_LOCAL_BUFFER (F77_INT, kronl, 1 + max (n, p));
        
        OCTAVE_LOCAL_BUFFER (double, af, ldaf * (n + min (p, m)));
        OCTAVE_LOCAL_BUFFER (double, bf, ldbf * (n + m));

        // workspace
        F77_INT s = max (m, p);
        F77_INT ldwork = max (1, max (s, n) + max (3*s-1, n+s));
        
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, s);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        F77_INT info;
        
        // tolerance
        double tol = 0;     // AB08ND uses DLAMCH for default tolerance

        // SLICOT routine AB08ND
        F77_XFCN (ab08nd, AB08ND,
                 (equil,
                  n, m, p,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  nu, rank, dinfz,
                  nkror, nkrol, infz,
                  kronr, kronl,
                  af, ldaf,
                  bf, ldbf,
                  tol,
                  iwork, dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("ss: zero: __sl_ab08nd__: exception in SLICOT subroutine AB08ND");
            
        if (info != 0)
            error ("ss: zero: __sl_ab08nd__: AB08ND returned info = %d", info);
        
        
        // DGGEV Part
        char jobvl = 'N';
        char jobvr = 'N';

        double* vl = 0;     // not referenced because jobvl = 'N'
        F77_INT ldvl = 1;
        double* vr = 0;     // not referenced because jobvr = 'N'
        F77_INT ldvr = 1;
        
        F77_INT lwork = max (1, 8*nu);
        OCTAVE_LOCAL_BUFFER (double, work, lwork);
        
        ColumnVector alphar (nu);
        ColumnVector alphai (nu);
        ColumnVector beta (nu);
        
        F77_INT info2;
        
        F77_XFCN (dggev, DGGEV,
                 (jobvl, jobvr,
                  nu,
                  af, ldaf,
                  bf, ldbf,
                  alphar.fortran_vec (), alphai.fortran_vec (),
                  beta.fortran_vec (),
                  vl, ldvl,
                  vr, ldvr,
                  work, lwork,
                  info2));
                                 
        if (f77_exception_encountered)
            error ("ss: zero: __sl_ab08nd__: exception in LAPACK subroutine DGGEV");
            
        if (info2 != 0)
            error ("ss: zero: __sl_ab08nd__: DGGEV returned info = %d", info2);

        // calculate gain
        octave_value gain = Matrix (0, 0);;

        if (m == 1 && p == 1)
        {
            if (nu < n)
                gain = c * xpow (a, double (n-1-nu)).matrix_value() * b;
            else
                gain = d;
        }

        // assemble complex vector - adapted from DEFUN complex in data.cc
        ColumnVector zeror (nu);
        ColumnVector zeroi (nu);

        zeror = quotient (alphar, beta);
        zeroi = quotient (alphai, beta);

        ComplexColumnVector zero (nu, Complex ());

        for (F77_INT i = 0; i < nu; i++)
            zero.xelem (i) = Complex (zeror(i), zeroi(i));

        // prepare additional outputs for info struct
        RowVector infzr (dinfz);
        RowVector kronrr (nkror);
        RowVector kronlr (nkrol);
        
        for (F77_INT i = 0; i < dinfz; i++)
            infzr.xelem (i) = infz[i];
        
        for (F77_INT i = 0; i < nkror; i++)
            kronrr.xelem (i) = kronr[i];
        
        for (F77_INT i = 0; i < nkrol; i++)
            kronlr.xelem (i) = kronl[i];

        // return values
        retval(0) = zero;
        retval(1) = gain;
        retval(2) = octave_value (rank);
        retval(3) = infzr;
        retval(4) = kronrr;
        retval(5) = kronlr;
    }
    
    return retval;
}
