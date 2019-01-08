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

Model reduction based on Balance & Truncate (B&T) or
Singular Perturbation Approximation (SPA) method.
Uses SLICOT AB09ID by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab09id, AB09ID)
                 (char& DICO, char& JOBC, char& JOBO, char& JOB,
                  char& WEIGHT, char& EQUIL, char& ORDSEL,
                  F77_INT& N, F77_INT& M, F77_INT& P,
                  F77_INT& NV, F77_INT& PV, F77_INT& NW, F77_INT& MW,
                  F77_INT& NR,
                  double& ALPHA, double& ALPHAC, double& ALPHAO,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  double* AV, F77_INT& LDAV,
                  double* BV, F77_INT& LDBV,
                  double* CV, F77_INT& LDCV,
                  double* DV, F77_INT& LDDV,
                  double* AW, F77_INT& LDAW,
                  double* BW, F77_INT& LDBW,
                  double* CW, F77_INT& LDCW,
                  double* DW, F77_INT& LDDW,
                  F77_INT& NS,
                  double* HSV,
                  double& TOL1, double& TOL2,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& IWARN, F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_ab09id__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_ab09id__, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB09ID Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 25)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char jobc;
        char jobo;
        char job;
        char weight;
        char equil;
        char ordsel;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        
        const F77_INT idico = args(4).int_value ();
        const F77_INT iequil = args(5).int_value ();
        F77_INT nr = args(6).int_value ();
        const F77_INT iordsel = args(7).int_value ();
        double alpha = args(8).double_value ();
        const F77_INT ijob = args(9).int_value ();
                       
        Matrix av = args(10).matrix_value ();
        Matrix bv = args(11).matrix_value ();
        Matrix cv = args(12).matrix_value ();
        Matrix dv = args(13).matrix_value ();
      
        Matrix aw = args(14).matrix_value ();
        Matrix bw = args(15).matrix_value ();
        Matrix cw = args(16).matrix_value ();
        Matrix dw = args(17).matrix_value ();
        
        const F77_INT iweight = args(18).int_value ();
        const F77_INT ijobc = args(19).int_value ();
        double alphac = args(20).double_value ();
        const F77_INT ijobo = args(21).int_value ();
        double alphao = args(22).double_value ();

        double tol1 = args(23).double_value ();
        double tol2 = args(24).double_value ();

        if (idico == 0)
            dico = 'C';
        else
            dico = 'D';

        if (iequil == 0)
            equil = 'S';
        else
            equil = 'N';

        if (iordsel == 0)
            ordsel = 'F';
        else
            ordsel = 'A';

        if (ijobc == 0)
            jobc = 'S';
        else
            jobc = 'E';

        if (ijobo == 0)
            jobo = 'S';
        else
            jobo = 'E';

        switch (ijob)
        {
            case 0:
                job = 'B';
                break;
            case 1:
                job = 'F';
                break;
            case 2:
                job = 'S';
                break;
            case 3:
                job = 'P';
                break;
            default:
                error ("__sl_ab09id__: argument job invalid");
        }

        switch (iweight)
        {
            case 0:
                weight = 'N';
                break;
            case 1:
                weight = 'L';
                break;
            case 2:
                weight = 'R';
                break;
            case 3:
                weight = 'B';
                break;
            default:
                error ("__sl_ab09id__: argument weight invalid");
        }

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        
        F77_INT nv = TO_F77_INT (av.rows ());
        F77_INT pv = TO_F77_INT (cv.rows ());
        F77_INT nw = TO_F77_INT (aw.rows ());
        F77_INT mw = TO_F77_INT (bw.columns ());

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldd = max (1, p);

        F77_INT ldav = max (1, nv);
        F77_INT ldbv = max (1, nv);
        F77_INT ldcv = max (1, pv);
        F77_INT lddv = max (1, pv);

        F77_INT ldaw = max (1, nw);
        F77_INT ldbw = max (1, nw);
        F77_INT ldcw = max (1, m);
        F77_INT lddw = max (1, m);

        // arguments out
        F77_INT ns;
        ColumnVector hsv (n);

        // workspace
        F77_INT liwork;
        F77_INT liwrk1;
        F77_INT liwrk2;
        F77_INT liwrk3;

        switch (job)
        {
            case 'B':
                liwrk1 = 0;
                break;
            case 'F':
                liwrk1 = n;
                break;
            default:
                liwrk1 = 2*n;
        }

        if (nv == 0 || weight == 'R' || weight == 'N')
            liwrk2 = 0;
        else
            liwrk2 = nv + max (p, pv);

        if (nw == 0 || weight == 'L' || weight == 'N')
            liwrk3 = 0;
        else
            liwrk3 = nw + max (m, mw);

        liwork = max (3, liwrk1, liwrk2, liwrk3);

        F77_INT ldwork;
        F77_INT lminl;
        F77_INT lrcf;
        F77_INT lminr;
        F77_INT llcf;
        F77_INT lleft;
        F77_INT lright;

        if (nw == 0 || weight == 'L' || weight == 'N')
        {
            lrcf = 0;
            lminr = 0;
        }
        else
        {
            lrcf = mw*(nw+mw) + max (nw*(nw+5), mw*(mw+2), 4*mw, 4*m);
            if (m == mw)
                lminr = nw + max (nw, 3*m);
            else
                lminr = 2*nw*max (m, mw) + nw + max (nw, 3*m, 3*mw);
        }

        llcf = pv*(nv+pv) + pv*nv + max (nv*(nv+5), pv*(pv+2), 4*pv, 4*p);

        if (nv == 0 || weight == 'R' || weight == 'N')
            lminl = 0;
        else if (p == pv)
            lminl  = max (llcf, nv + max (nv, 3*p));
        else
            lminl  = max (p, pv) * (2*nv + max (p, pv)) + max (llcf, nv + max (nv, 3*p, 3*pv));


        if (pv == 0 || weight == 'R' || weight == 'N')
            lleft = n*(p+5);
        else
            lleft = (n+nv) * (n + nv + max (n+nv, pv) + 5);

        if (mw == 0 || weight == 'L' || weight == 'N')
            lright = n*(m+5);
        else
            lright = (n+nw) * (n + nw + max (n+nw, mw) + 5);

        ldwork =  max (lminl, lminr, lrcf,
                       2*n*n + max (1, lleft, lright, 2*n*n+5*n, n*max (m, p)));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT iwarn = 0;
        F77_INT info = 0;


        // SLICOT routine AB09ID
        F77_XFCN (ab09id, AB09ID,
                 (dico, jobc, jobo, job,
                  weight, equil, ordsel,
                  n, m, p,
                  nv, pv, nw, mw,
                  nr,
                  alpha, alphac, alphao,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  av.fortran_vec (), ldav,
                  bv.fortran_vec (), ldbv,
                  cv.fortran_vec (), ldcv,
                  dv.fortran_vec (), lddv,
                  aw.fortran_vec (), ldaw,
                  bw.fortran_vec (), ldbw,
                  cw.fortran_vec (), ldcw,
                  dw.fortran_vec (), lddw,
                  ns,
                  hsv.fortran_vec (),
                  tol1, tol2,
                  iwork,
                  dwork, ldwork,
                  iwarn, info));

        if (f77_exception_encountered)
            error ("modred: exception in SLICOT subroutine AB09ID");


        static const char* err_msg[] = {
            "0: OK",
            "1: the computation of the ordered real Schur form of A "
                "failed",
            "2: the separation of the ALPHA-stable/unstable "
                "diagonal blocks failed because of very close "
                "eigenvalues",
            "3: the reduction to a real Schur form of the state "
                "matrix of a minimal realization of V failed",
            "4: a failure was detected during the ordering of the "
                "real Schur form of the state matrix of a minimal "
                "realization of V or in the iterative process to "
                "compute a left coprime factorization with inner "
                "denominator",
            "5: if DICO = 'C' and the matrix AV has an observable "
                "eigenvalue on the imaginary axis, or DICO = 'D' and "
                "AV has an observable eigenvalue on the unit circle",
            "6: the reduction to a real Schur form of the state "
                "matrix of a minimal realization of W failed",
            "7: a failure was detected during the ordering of the "
                "real Schur form of the state matrix of a minimal "
                "realization of W or in the iterative process to "
                "compute a right coprime factorization with inner "
                "denominator",
            "8: if DICO = 'C' and the matrix AW has a controllable "
                "eigenvalue on the imaginary axis, or DICO = 'D' and "
                "AW has a controllable eigenvalue on the unit circle",
            "9: the computation of eigenvalues failed",
            "10: the computation of Hankel singular values failed"};

        static const char* warn_msg[] = {
            "0: OK",
            "1: with ORDSEL = 'F', the selected order NR is greater "
                "than NSMIN, the sum of the order of the "
                "ALPHA-unstable part and the order of a minimal "
                "realization of the ALPHA-stable part of the given "
                "system; in this case, the resulting NR is set equal "
                "to NSMIN.",
            "2: with ORDSEL = 'F', the selected order NR corresponds "
                "to repeated singular values for the ALPHA-stable "
                "part, which are neither all included nor all "
                "excluded from the reduced model; in this case, the "
                "resulting NR is automatically decreased to exclude "
                "all repeated singular values.",
            "3: with ORDSEL = 'F', the selected order NR is less "
                "than the order of the ALPHA-unstable part of the "
                "given system; in this case NR is set equal to the "
                "order of the ALPHA-unstable part.",
/* 10+%d: %d */ "violations of the numerical stability condition "
                "occurred during the assignment of eigenvalues in the "
                "SLICOT Library routines SB08CD and/or SB08DD."};


        error_msg ("modred", info, 10, err_msg);
        warning_msg ("modred", iwarn, 3, warn_msg, 10);

        // resize
        a.resize (nr, nr);
        b.resize (nr, m);
        c.resize (p, nr);
        hsv.resize (ns);
        
        // return values
        retval(0) = a;
        retval(1) = b;
        retval(2) = c;
        retval(3) = d;
        retval(4) = octave_value (nr);
        retval(5) = hsv;
        retval(6) = octave_value (ns);
    }
    
    return retval;
}
