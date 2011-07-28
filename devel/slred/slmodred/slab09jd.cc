/*

Copyright (C) 2011   Lukas F. Reichlin

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
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.cc"

extern "C"
{ 
    int F77_FUNC (ab09jd, AB09JD)
                 (char& JOBV, char& JOBW, char& JOBINV,
                  char& DICO, char& EQUIL, char& ORDSEL,
                  int& N, int& NV, int& NW, int& M, int& P,
                  int& NR,
                  double& ALPHA,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* AV, int& LDAV,
                  double* BV, int& LDBV,
                  double* CV, int& LDCV,
                  double* DV, int& LDDV,
                  double* AW, int& LDAW,
                  double* BW, int& LDBW,
                  double* CW, int& LDCW,
                  double* DW, int& LDDW,
                  int& NS,
                  double* HSV,
                  double& TOL1, double& TOL2,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}
     
DEFUN_DLD (slab09jd, args, nargout,
   "-*- texinfo -*-\n\
Slicot AB09JD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
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
        
        Matrix av = args(4).matrix_value ();
        Matrix bv = args(5).matrix_value ();
        Matrix cv = args(6).matrix_value ();
        Matrix dv = args(7).matrix_value ();
        
        Matrix aw = args(8).matrix_value ();
        Matrix bw = args(9).matrix_value ();
        Matrix cw = args(10).matrix_value ();
        Matrix dw = args(11).matrix_value ();

        int nr = args(12).int_value ();
        double alpha = args(13).double_value ();
        
        const int ijobv = args(14).int_value ();
        const int ijobw = args(15).int_value ();
        const int ijobinv = args(16).int_value ();
        const int idico = args(17).int_value ();
        const int iequil = args(18).int_value ();
        const int iordsel = args(19).int_value ();
        
        double tol1 = args(20).double_value ();
        double tol2 = args(21).double_value ();
        
        if (idico == 0)
            dico = 'C';
        else
            dico = 'D';

        if (iequil == 0)
            equil = 'S';
        else
            equil = 'N';


        int n = a.rows ();      // n: number of states
        int nv = av.rows ();
        int nw = aw.rows ();
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        int ldd = max (1, p);

        int ldav = max (1, nv);
        int ldbv = max (1, nv);
        int ldcv = max (1, p);
        int lddv = max (1, p);

        int ldaw = max (1, nw);
        int ldbw = max (1, nw);
        int ldcw = max (1, m);
        int lddw = max (1, m);

        // arguments out
        int ns = 0;

        ColumnVector hsv (n);

        // workspace
        int liwork;
        int tmpc;
        int tmpd;

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

        int ldwork;
        int nvp = nv + p;
        int nwm = nw + m;
        int ldw1;
        int ldw2;
        int ldw3 = n*(2*n + max (n, m, p) + 5) + n*(n+1)/2;
        int ldw4 = n*(m+p+2) + 2*m*p + min (n, m) + max (3*m+1, min (n, m) + p);
        
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

        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int iwarn = 0;
        int info = 0;


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
            error ("hsvd: slab09jd: exception in SLICOT subroutine AB09JD");
            
        if (info != 0)
            error ("hsvd: slab09jd: AB09JD returned info = %d", info);

        // resize
        hsv.resize (ns);
        
        // return values
        retval(0) = hsv;
        retval(1) = octave_value (ns);
    }
    
    return retval;
}
