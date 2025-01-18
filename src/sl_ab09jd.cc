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

Model reduction based on Hankel-norm approximation method.
Uses SLICOT AB09JD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: July 2011
Version: 0.2

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ab09jd, AB09JD)
                 (char& JOBV, char& JOBW, char& JOBINV,
                  char& DICO, char& EQUIL, char& ORDSEL,
                  F77_INT& N, F77_INT& NV, F77_INT& NW, F77_INT& M, F77_INT& P,
                  F77_INT& NR,
                  double& ALPHA,
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

// PKG_ADD: autoload ("__sl_ab09jd__", "__control_slicot_functions__.oct");         
DEFUN_DLD (__sl_ab09jd__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __sl_ab09jd__ (@dots{})\n"
   "Wrapper for SLICOT function AB09JD.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_idx_type nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 22)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char jobv;
        char jobw;
        char jobinv;
        char dico;
        char equil;
        char ordsel;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();

        if (a.any_element_is_inf_or_nan () || b.any_element_is_inf_or_nan () ||
            c.any_element_is_inf_or_nan () || d.any_element_is_inf_or_nan ())
          error ("__sl_ab09jd__: inputs must not contain NaN or Inf\n");

        const F77_INT idico = args(4).int_value ();
        const F77_INT iequil = args(5).int_value ();
        F77_INT nr = args(6).int_value ();
        const F77_INT iordsel = args(7).int_value ();
        double alpha = args(8).double_value ();
                 
        const F77_INT ijobv = args(9).int_value ();       
        Matrix av = args(10).matrix_value ();
        Matrix bv = args(11).matrix_value ();
        Matrix cv = args(12).matrix_value ();
        Matrix dv = args(13).matrix_value ();

        if (av.any_element_is_inf_or_nan () || bv.any_element_is_inf_or_nan () ||
            cv.any_element_is_inf_or_nan () || dv.any_element_is_inf_or_nan ())
          error ("__sl_ab09jd__: inputs must not contain NaN or Inf\n");

        const F77_INT ijobw = args(14).int_value ();        
        Matrix aw = args(15).matrix_value ();
        Matrix bw = args(16).matrix_value ();
        Matrix cw = args(17).matrix_value ();
        Matrix dw = args(18).matrix_value ();

        if (aw.any_element_is_inf_or_nan () || bw.any_element_is_inf_or_nan () ||
            cw.any_element_is_inf_or_nan () || dw.any_element_is_inf_or_nan ())
          error ("__sl_ab09jd__: inputs must not contain NaN or Inf\n");

        const F77_INT ijobinv = args(19).int_value ();
        double tol1 = args(20).double_value ();
        double tol2 = args(21).double_value ();

        switch (ijobv)
        {
            case 0:
                jobv = 'N';
                break;
            case 1:
                jobv = 'V';
                break;
            case 2:
                jobv = 'I';
                break;
            case 3:
                jobv = 'C';
                break;
            case 4:
                jobv = 'R';
                break;
            default:
                error ("__sl_ab09jd__: argument jobv invalid");
        }

        switch (ijobw)
        {
            case 0:
                jobw = 'N';
                break;
            case 1:
                jobw = 'W';
                break;
            case 2:
                jobw = 'I';
                break;
            case 3:
                jobw = 'C';
                break;
            case 4:
                jobw = 'R';
                break;
            default:
                error ("__sl_ab09jd__: argument jobw invalid");
        }
            
        switch (ijobinv)
        {
            case 0:
                jobinv = 'N';
                break;
            case 1:
                jobinv = 'I';
                break;
            case 2:
                jobinv = 'A';
                break;
            default:
                error ("__sl_ab09jd__: argument jobinv invalid");
        }

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

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT nv = TO_F77_INT (av.rows ());
        F77_INT nw = TO_F77_INT (aw.rows ());
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldd = max (1, p);

        F77_INT ldav = max (1, nv);
        F77_INT ldbv = max (1, nv);
        F77_INT ldcv = max (1, p);
        F77_INT lddv = max (1, p);

        F77_INT ldaw = max (1, nw);
        F77_INT ldbw = max (1, nw);
        F77_INT ldcw = max (1, m);
        F77_INT lddw = max (1, m);

        // arguments out
        F77_INT ns;
        ColumnVector hsv (n);

        // workspace
        F77_INT liwork;
        F77_INT tmpc;
        F77_INT tmpd;

        if (jobv == 'N')
            tmpc = 0;
        else
            tmpc = max (2*p, nv+p+n+6, 2*nv+p+2);

        if (jobw == 'N')
            tmpd = 0;
        else
            tmpd = max (2*m, nw+m+n+6, 2*nw+m+2);
        
        if (dico == 'C')
            liwork = max (1, m, tmpc, tmpd);
        else
            liwork = max (1, n, m, tmpc, tmpd);

        F77_INT ldwork;
        F77_INT nvp = nv + p;
        F77_INT nwm = nw + m;
        F77_INT ldw1;
        F77_INT ldw2;
        F77_INT ldw3 = n*(2*n + max (n, m, p) + 5) + n*(n+1)/2;
        F77_INT ldw4 = n*(m+p+2) + 2*m*p + min (n, m) + max (3*m+1, min (n, m) + p);
        
        if (jobv == 'N')
        {
            ldw1 = 0;
        }
        else
        {
            ldw1 = 2*nvp*(nvp+p) + p*p + max (2*nvp*nvp + max (11*nvp+16, p*nvp),
                                              nvp*n + max (nvp*n+n*n, p*n, p*m));
        }

        if (jobw == 'N')
        {
            ldw2 = 0;
        }
        else
        {
            ldw2 = 2*nwm*(nwm+m) + m*m + max (2*nwm*nwm + max (11*nwm+16, m*nwm),
                                              nwm*n + max (nwm*n+n*n, m*n, p*m));
        }
        
        ldwork = max (ldw1, ldw2, ldw3, ldw4);

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        F77_INT iwarn = 0;
        F77_INT info = 0;


        // SLICOT routine AB09JD
        F77_XFCN (ab09jd, AB09JD,
                 (jobv, jobw, jobinv,
                  dico, equil, ordsel,
                  n, nv, nw, m, p,
                  nr,
                  alpha,
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
            error ("hnamodred: exception in SLICOT subroutine AB09JD");


        static const char* err_msg[] = {
            "0: OK",
            "1: the computation of the ordered real Schur form of A "
                "failed",

            "2: the separation of the ALPHA-stable/unstable "
                "diagonal blocks failed because of very close eigenvalues",

            "3: the reduction of AV to a real Schur form failed",

            "4: the reduction of AW to a real Schur form failed",

            "5: the reduction to generalized Schur form of the "
                "descriptor pair corresponding to the inverse of V "
                "failed",

            "6: the reduction to generalized Schur form of the "
                "descriptor pair corresponding to the inverse of W "
                "failed",

            "7: the computation of Hankel singular values failed",

            "8: the computation of stable projection in the "
                "Hankel-norm approximation algorithm failed",

            "9: the order of computed stable projection in the "
                "Hankel-norm approximation algorithm differs "
                "from the order of Hankel-norm approximation",

            "10: the reduction of AV-BV*inv(DV)*CV to a "
                "real Schur form failed",

            "11: the reduction of AW-BW*inv(DW)*CW to a "
                "real Schur form failed",

            "12: the solution of the Sylvester equation failed "
                "because the poles of V (if JOBV = 'V') or of "
                "conj(V) (if JOBV = 'C') are not distinct from "
                "the poles of G1 (see METHOD)",

            "13: the solution of the Sylvester equation failed "
                "because the poles of W (if JOBW = 'W') or of "
                "conj(W) (if JOBW = 'C') are not distinct from "
                "the poles of G1 (see METHOD)",

            "14: the solution of the Sylvester equation failed "
                "because the zeros of V (if JOBV = 'I') or of "
                "conj(V) (if JOBV = 'R') are not distinct from "
                "the poles of G1sr (see METHOD)",

            "15: the solution of the Sylvester equation failed "
                "because the zeros of W (if JOBW = 'I') or of "
                "conj(W) (if JOBW = 'R') are not distinct from "
                "the poles of G1sr (see METHOD)",

            "16: the solution of the generalized Sylvester system "
                "failed because the zeros of V (if JOBV = 'I') or "
                "of conj(V) (if JOBV = 'R') are not distinct from "
                "the poles of G1sr (see METHOD)",

            "17: the solution of the generalized Sylvester system "
                "failed because the zeros of W (if JOBW = 'I') or "
                "of conj(W) (if JOBW = 'R') are not distinct from "
                "the poles of G1sr (see METHOD)",

            "18: op(V) is not antistable",

            "19: op(W) is not antistable",

            "20: V is not invertible",

            "21: W is not invertible"};

        static const char* warn_msg[] = {
            "0: OK",
            "1: with ORDSEL = 'F', the selected order NR is greater "
                "than NSMIN, the sum of the order of the "
                "ALPHA-unstable part and the order of a minimal "
                "realization of the ALPHA-stable part of the given "
                "system. In this case, the resulting NR is set equal "
                "to NSMIN.",
            "2: with ORDSEL = 'F', the selected order NR is less "
                "than the order of the ALPHA-unstable part of the "
                "given system. In this case NR is set equal to the "
                "order of the ALPHA-unstable part."};

        error_msg ("hnamodred", info, 21, err_msg);
        warning_msg ("hnamodred", iwarn, 2, warn_msg);

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
