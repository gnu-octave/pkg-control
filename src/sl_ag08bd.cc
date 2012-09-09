/*

Copyright (C) 2010, 2011   Lukas F. Reichlin

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

Transmission zeros of descriptor state-space models.
Uses SLICOT AG08BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: September 2010
Version: 0.3

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
                  int& L, int& N, int& M, int& P,
                  double* A, int& LDA,
                  double* E, int& LDE,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  int& NFZ, int& NRANK, int& NIZ, int& DINFZ,
                  int& NKROR, int& NINFE, int& NKROL,
                  int* INFZ,
                  int* KRONR, int* INFE, int* KRONL,
                  double& TOL,
                  int* IWORK, double* DWORK, int& LDWORK,
                  int& INFO);
                                   
    int F77_FUNC (dggev, DGGEV)
                 (char& JOBVL, char& JOBVR,
                  int& N,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* ALPHAR, double* ALPHAI,
                  double* BETA,
                  double* VL, int& LDVL,
                  double* VR, int& LDVR,
                  double* WORK, int& LWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_ag08bd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ag08bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AG08BD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
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
        const int scaled = args(5).int_value ();
        
        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        int l = a.rows ();      // l: number of states
        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs
        
        int lda = max (1, l);
        int lde = max (1, l);
        int ldb = max (1, l);

        if (m == 0)
            ldb = 1;

        int ldc = max (1, p);
        int ldd = max (1, p);
        
        // arguments out
        int nfz;
        int nrank;
        int niz;
        int dinfz;
        int nkror;
        int ninfe;
        int nkrol;

        OCTAVE_LOCAL_BUFFER (int, infz, n+1);
        OCTAVE_LOCAL_BUFFER (int, kronr, n+m+1);
        OCTAVE_LOCAL_BUFFER (int, infe, 1 + min (l+p, n+m));
        OCTAVE_LOCAL_BUFFER (int, kronl, l+p+1);

        // workspace
        int ldwork = max (l+p, m+n) * (m+n) + max (1, 5 * max (l+p, m+n));
        
        OCTAVE_LOCAL_BUFFER (int, iwork, n + max (1, m));
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicator
        int info;
        
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
        int ldvl = 1;
        double* vr = 0;     // not referenced because jobvr = 'N'
        int ldvr = 1;
        
        int lwork = max (1, 8*nfz);
        OCTAVE_LOCAL_BUFFER (double, work, lwork);
        
        int info2;
        
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

        // return values
        retval(0) = zero;
    }
    
    return retval;
}
