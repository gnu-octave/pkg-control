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

Finite Smith zeros of descriptor state-space models.
Uses SLICOT AG08BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.5

*/

#include <octave/oct.h>
#include "common.h"
#include <complex>
#include <xpow.h>

extern "C"
{ 
    int F77_FUNC (ag08bd, AG08BD)
                 (char& EQUIL,
                  F77_INT& L, F77_INT& N, F77_INT& M, F77_INT& P,
                  double* A, F77_INT& LDA,
                  double* E, F77_INT& LDE,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  F77_INT& NFZ, F77_INT& NRANK, F77_INT& NIZ, F77_INT& DINFZ,
                  F77_INT& NKROR, F77_INT& NINFE, F77_INT& NKROL,
                  F77_INT* INFZ,
                  F77_INT* KRONR, F77_INT* INFE, F77_INT* KRONL,
                  double& TOL,
                  F77_INT* IWORK, double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
                                   
    int F77_FUNC (dggev, DGGEV)
                 (char& JOBVL, char& JOBVR,
                  F77_INT& N,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  double* VL, F77_INT& LDVL,
                  double* VR, F77_INT& LDVR,
                  double* WORK, F77_INT& LWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ag08bd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ag08bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AG08BD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 6)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char equil;
        
        Matrix a = args(0).matrix_value ();
        Matrix e = args(1).matrix_value (); 
        Matrix b = args(2).matrix_value ();
        Matrix c = args(3).matrix_value ();
        Matrix d = args(4).matrix_value ();
        const F77_INT scaled = args(5).int_value ();
        
        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        F77_INT l = TO_F77_INT (a.rows ());      // l: number of states
        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        
        F77_INT lda = max (1, l);
        F77_INT lde = max (1, l);
        F77_INT ldb = max (1, l);

        if (m == 0)
            ldb = 1;

        F77_INT ldc = max (1, p);
        F77_INT ldd = max (1, p);
        
        // arguments out
        F77_INT nfz;
        F77_INT nrank;
        F77_INT niz;
        F77_INT dinfz;
        F77_INT nkror;
        F77_INT ninfe;
        F77_INT nkrol;

        OCTAVE_LOCAL_BUFFER (F77_INT, infz, n+1);
        OCTAVE_LOCAL_BUFFER (F77_INT, kronr, n+m+1);
        OCTAVE_LOCAL_BUFFER (F77_INT, infe, 1 + min (l+p, n+m));
        OCTAVE_LOCAL_BUFFER (F77_INT, kronl, l+p+1);

        // workspace
        F77_INT ldwork = max (l+p, m+n) * (m+n) + max (1, 5 * max (l+p, m+n));
        
        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, n + max (1, m));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        F77_INT info;
        
        // tolerance
        double tol = 0;     // AG08BD uses DLAMCH for default tolerance

        // SLICOT routine AG08BD
        F77_XFCN (ag08bd, AG08BD,
                 (equil,
                  l, n, m, p,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  nfz, nrank, niz, dinfz,
                  nkror, ninfe, nkrol,
                  infz,
                  kronr, infe, kronl,
                  tol,
                  iwork, dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("dss: zero: __sl_ag08bd__: exception in SLICOT subroutine AG08BD");
            
        if (info != 0)
            error ("dss: zero: __sl_ag08bd__: AG08BD returned info = %d", info);


        // DGGEV Part
        a.resize (nfz, nfz);  // Af
        e.resize (nfz, nfz);  // Ef

        lda = max (1, nfz);
        lde = max (1, nfz); 
        
        char jobvl = 'N';
        char jobvr = 'N';

        ColumnVector alphar (nfz);
        ColumnVector alphai (nfz);
        ColumnVector beta (nfz);

        double* vl = 0;     // not referenced because jobvl = 'N'
        F77_INT ldvl = 1;
        double* vr = 0;     // not referenced because jobvr = 'N'
        F77_INT ldvr = 1;
        
        F77_INT lwork = max (1, 8*nfz);
        OCTAVE_LOCAL_BUFFER (double, work, lwork);
        
        F77_INT info2;
        
        F77_XFCN (dggev, DGGEV,
                 (jobvl, jobvr,
                  nfz,
                  a.fortran_vec (), lda,
                  e.fortran_vec (), lde,
                  alphar.fortran_vec (), alphai.fortran_vec (),
                  beta.fortran_vec (),
                  vl, ldvl,
                  vr, ldvr,
                  work, lwork,
                  info2));
                                 
        if (f77_exception_encountered)
            error ("dss: zero: __sl_ag08bd__: exception in LAPACK subroutine DGGEV");
            
        if (info2 != 0)
            error ("dss: zero: __sl_ag08bd__: DGGEV returned info = %d", info2);

        // assemble complex vector - adapted from DEFUN complex in data.cc
        // LAPACK DGGEV.f says:
        //
        // Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
        // may easily over- or underflow, and BETA(j) may even be zero.
        // Thus, the user should avoid naively computing the ratio
        // alpha/beta.  However, ALPHAR and ALPHAI will be always less
        // than and usually comparable with norm(A) in magnitude, and
        // BETA always less than and usually comparable with norm(B).
        //
        // Since we need the zeros explicitly ...

        ColumnVector zeror (nfz);
        ColumnVector zeroi (nfz);

        zeror = quotient (alphar, beta);
        zeroi = quotient (alphai, beta);

        ComplexColumnVector zero (nfz, Complex ());

        for (F77_INT i = 0; i < nfz; i++)
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
        retval(1) = octave_value (nrank);
        retval(2) = infzr;
        retval(3) = kronrr;
        retval(4) = kronlr;
    }
    
    return retval;
}
