/*

Copyright (C) 2009, 2010, 2011   Lukas F. Reichlin

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

Transmission zeros of state-space models.
Uses SLICOT AB08ND by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2009
Version: 0.5

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"
#include <complex>
#include <xpow.h>

extern "C"
{ 
    int F77_FUNC (ab08nd, AB08ND)
                 (char& EQUIL,
                  int& N, int& M, int& P,
                  const double* A, int& LDA,
                  const double* B, int& LDB,
                  const double* C, int& LDC,
                  const double* D, int& LDD,
                  int& NU, int& RANK, int& DINFZ,
                  int& NKROR, int& NKROL, int* INFZ,
                  int* KRONR, int* KRONL,
                  double* AF, int& LDAF,
                  double* BF, int& LDBF,
                  double& TOL,
                  int* IWORK, double* DWORK, int& LDWORK,
                  int& INFO);
                                   
    int F77_FUNC (dggev, DGGEV)
                 (char& JOBVL, char& JOBVR,
                  int& N,
                  double* AF, int& LDAF,
                  double* BF, int& LDBF,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  double* VL, int& LDVL,
                  double* VR, int& LDVR,
                  double* WORK, int& LWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_ab08nd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab08nd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB08ND Release 5.0\n\
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
        char equil;
        
        const Matrix a = args(0).matrix_value ();
        const Matrix b = args(1).matrix_value ();
        const Matrix c = args(2).matrix_value ();
        const Matrix d = args(3).matrix_value ();
        const int scaled = args(4).int_value ();
        
        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';
        
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs
        
        int lda = max (1, a.rows ());
        int ldb = max (1, b.rows ());
        int ldc = max (1, c.rows ());
        int ldd = max (1, d.rows ());
        
        // arguments out
        int nu;
        int rank;
        int dinfz;
        int nkror;
        int nkrol;
        
        int ldaf = max (1, n + m);
        int ldbf = max (1, n + p);

        OCTAVE_LOCAL_BUFFER (int, infz, n);
        OCTAVE_LOCAL_BUFFER (int, kronr, 1 + max (n, m));
        OCTAVE_LOCAL_BUFFER (int, kronl, 1 + max (n, p));
        
        OCTAVE_LOCAL_BUFFER (double, af, ldaf * (n + min (p, m)));
        OCTAVE_LOCAL_BUFFER (double, bf, ldbf * (n + m));

        // workspace
        int s = max (m, p);
        int ldwork = max (s, n) + max (3*s-1, n+s);
        
        OCTAVE_LOCAL_BUFFER (int, iwork, s);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        int info;
        
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
        int ldvl = 1;
        double* vr = 0;     // not referenced because jobvr = 'N'
        int ldvr = 1;
        
        int lwork = max (1, 8*nu);
        OCTAVE_LOCAL_BUFFER (double, work, lwork);
        
        ColumnVector alphar (nu);
        ColumnVector alphai (nu);
        ColumnVector beta (nu);
        
        int info2;
        
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
                gain = c * xpow (a, double (n-1-nu)) * b;
            else
                gain = d;
        }

        // assemble complex vector - adapted from DEFUN complex in data.cc
        ColumnVector zeror (nu);
        ColumnVector zeroi (nu);

        zeror = quotient (alphar, beta);
        zeroi = quotient (alphai, beta);

        ComplexColumnVector zero (nu, Complex ());

        for (octave_idx_type i = 0; i < nu; i++)
            zero.xelem (i) = Complex (zeror(i), zeroi(i));

        // return values
        retval(0) = zero;
        retval(1) = gain;
    }
    
    return retval;
}
