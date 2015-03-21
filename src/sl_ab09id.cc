/*

Copyright (C) 2009-2015   Lukas F. Reichlin

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
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab09id, AB09ID)
                 (char& DICO, char& JOBC, char& JOBO, char& JOB,
                  char& WEIGHT, char& EQUIL, char& ORDSEL,
                  octave_idx_type& N, octave_idx_type& M, octave_idx_type& P,
                  octave_idx_type& NV, octave_idx_type& PV, octave_idx_type& NW, octave_idx_type& MW,
                  octave_idx_type& NR,
                  double& ALPHA, double& ALPHAC, double& ALPHAO,
                  double* A, octave_idx_type& LDA,
                  double* B, octave_idx_type& LDB,
                  double* C, octave_idx_type& LDC,
                  double* D, octave_idx_type& LDD,
                  double* AV, octave_idx_type& LDAV,
                  double* BV, octave_idx_type& LDBV,
                  double* CV, octave_idx_type& LDCV,
                  double* DV, octave_idx_type& LDDV,
                  double* AW, octave_idx_type& LDAW,
                  double* BW, octave_idx_type& LDBW,
                  double* CW, octave_idx_type& LDCW,
                  double* DW, octave_idx_type& LDDW,
                  octave_idx_type& NS,
                  double* HSV,
                  double& TOL1, double& TOL2,
                  octave_idx_type* IWORK,
                  double* DWORK, octave_idx_type& LDWORK,
                  octave_idx_type& IWARN, octave_idx_type& INFO);
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
        
        const octave_idx_type idico = args(4).int_value ();
        const octave_idx_type iequil = args(5).int_value ();
        octave_idx_type nr = args(6).int_value ();
        const octave_idx_type iordsel = args(7).int_value ();
        double alpha = args(8).double_value ();
        const octave_idx_type ijob = args(9).int_value ();
                       
        Matrix av = args(10).matrix_value ();
        Matrix bv = args(11).matrix_value ();
        Matrix cv = args(12).matrix_value ();
        Matrix dv = args(13).matrix_value ();
      
        Matrix aw = args(14).matrix_value ();
        Matrix bw = args(15).matrix_value ();
        Matrix cw = args(16).matrix_value ();
        Matrix dw = args(17).matrix_value ();
        
        const octave_idx_type iweight = args(18).int_value ();
        const octave_idx_type ijobc = args(19).int_value ();
        double alphac = args(20).double_value ();
        const octave_idx_type ijobo = args(21).int_value ();
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

        octave_idx_type n = a.rows ();      // n: number of states
        octave_idx_type m = b.columns ();   // m: number of inputs
        octave_idx_type p = c.rows ();      // p: number of outputs
        
        octave_idx_type nv = av.rows ();
        octave_idx_type pv = cv.rows ();
        octave_idx_type nw = aw.rows ();
        octave_idx_type mw = bw.columns ();

        octave_idx_type lda = max (1, n);
        octave_idx_type ldb = max (1, n);
        octave_idx_type ldc = max (1, p);
        octave_idx_type ldd = max (1, p);

        octave_idx_type ldav = max (1, nv);
        octave_idx_type ldbv = max (1, nv);
        octave_idx_type ldcv = max (1, pv);
        octave_idx_type lddv = max (1, pv);

        octave_idx_type ldaw = max (1, nw);
        octave_idx_type ldbw = max (1, nw);
        octave_idx_type ldcw = max (1, m);
        octave_idx_type lddw = max (1, m);

        // arguments out
        octave_idx_type ns;
        ColumnVector hsv (n);

        // workspace
        octave_idx_type liwork;
        octave_idx_type liwrk1;
        octave_idx_type liwrk2;
        octave_idx_type liwrk3;

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

        octave_idx_type ldwork;
        octave_idx_type lminl;
        octave_idx_type lrcf;
        octave_idx_type lminr;
        octave_idx_type llcf;
        octave_idx_type lleft;
        octave_idx_type lright;

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

        OCTAVE_LOCAL_BUFFER (octave_idx_type, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        octave_idx_type iwarn = 0;
        octave_idx_type info = 0;


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
                "occured during the assignment of eigenvalues in the "
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
