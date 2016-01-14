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
#include <f77-fcn.h>
#include "common.h"
#include <complex>
#include <xpow.h>

extern "C"
{ 
    int F77_FUNC (ag08bd, AG08BD)
                 (char& EQUIL,
                  octave_idx_type& L, octave_idx_type& N, octave_idx_type& M, octave_idx_type& P,
                  double* A, octave_idx_type& LDA,
                  double* E, octave_idx_type& LDE,
                  double* B, octave_idx_type& LDB,
                  double* C, octave_idx_type& LDC,
                  double* D, octave_idx_type& LDD,
                  octave_idx_type& NFZ, octave_idx_type& NRANK, octave_idx_type& NIZ, octave_idx_type& DINFZ,
                  octave_idx_type& NKROR, octave_idx_type& NINFE, octave_idx_type& NKROL,
                  octave_idx_type* INFZ,
                  octave_idx_type* KRONR, octave_idx_type* INFE, octave_idx_type* KRONL,
                  double& TOL,
                  octave_idx_type* IWORK, double* DWORK, octave_idx_type& LDWORK,
                  octave_idx_type& INFO);
                                   
    int F77_FUNC (dggev, DGGEV)
                 (char& JOBVL, char& JOBVR,
                  octave_idx_type& N,
                  double* A, octave_idx_type& LDA,
                  double* B, octave_idx_type& LDB,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  double* VL, octave_idx_type& LDVL,
                  double* VR, octave_idx_type& LDVR,
                  double* WORK, octave_idx_type& LWORK,
                  octave_idx_type& INFO);
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
        const octave_idx_type scaled = args(5).int_value ();
        
        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        octave_idx_type l = a.rows ();      // l: number of states
        octave_idx_type n = a.rows ();      // n: number of states
        octave_idx_type m = b.columns ();   // m: number of inputs
        octave_idx_type p = c.rows ();      // p: number of outputs
        
        octave_idx_type lda = max (1, l);
        octave_idx_type lde = max (1, l);
        octave_idx_type ldb = max (1, l);

        if (m == 0)
            ldb = 1;

        octave_idx_type ldc = max (1, p);
        octave_idx_type ldd = max (1, p);
        
        // arguments out
        octave_idx_type nfz;
        octave_idx_type nrank;
        octave_idx_type niz;
        octave_idx_type dinfz;
        octave_idx_type nkror;
        octave_idx_type ninfe;
        octave_idx_type nkrol;

        OCTAVE_LOCAL_BUFFER (octave_idx_type, infz, n+1);
        OCTAVE_LOCAL_BUFFER (octave_idx_type, kronr, n+m+1);
        OCTAVE_LOCAL_BUFFER (octave_idx_type, infe, 1 + min (l+p, n+m));
        OCTAVE_LOCAL_BUFFER (octave_idx_type, kronl, l+p+1);

        // workspace
        octave_idx_type ldwork = max (l+p, m+n) * (m+n) + max (1, 5 * max (l+p, m+n));
        
        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, n + max (1, m));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        octave_idx_type info;
        
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
        octave_idx_type ldvl = 1;
        double* vr = 0;     // not referenced because jobvr = 'N'
        octave_idx_type ldvr = 1;
        
        octave_idx_type lwork = max (1, 8*nfz);
        OCTAVE_LOCAL_BUFFER (double, work, lwork);
        
        octave_idx_type info2;
        
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

        for (octave_idx_type i = 0; i < nfz; i++)
            zero.xelem (i) = Complex (zeror(i), zeroi(i));
        
        // prepare additional outputs for info struct
        RowVector infzr (dinfz);
        RowVector kronrr (nkror);
        RowVector kronlr (nkrol);
        
        for (octave_idx_type i = 0; i < dinfz; i++)
            infzr.xelem (i) = infz[i];
        
        for (octave_idx_type i = 0; i < nkror; i++)
            kronrr.xelem (i) = kronr[i];
        
        for (octave_idx_type i = 0; i < nkrol; i++)
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
