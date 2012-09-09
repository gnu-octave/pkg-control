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

TODO
Uses SLICOT SB16CD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: November 2011
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (sb16cd, SB16CD)
                 (char& DICO, char& JOBD, char& JOBMR, char& JOBCF,
                  char& ORDSEL,
                  int& N, int& M, int& P,
                  int& NCR,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* F, int& LDF,
                  double* G, int& LDG,
                  double* HSV,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}

// PKG_ADD: autoload ("__sl_sb16cd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb16cd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB16CD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 13)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char jobd;
        char jobmr;
        char jobcf;
        char ordsel;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        
        const int idico = args(4).int_value ();
        int ncr = args(5).int_value ();
        const int iordsel = args(6).int_value ();
        const int ijobd = args(7).int_value ();
        const int ijobmr = args(8).int_value ();
                       
        Matrix f = args(9).matrix_value ();
        Matrix g = args(10).matrix_value ();

        const int ijobcf = args(11).int_value ();
        double tol = args(12).double_value ();

        if (idico == 0)
            dico = 'C';
        else
            dico = 'D';

        if (iordsel == 0)
            ordsel = 'F';
        else
            ordsel = 'A';

        if (ijobd == 0)
            jobd = 'Z';
        else
            jobd = 'D';

        if (ijobcf == 0)
            jobcf = 'L';
        else
            jobcf = 'R';

        switch (ijobmr)
        {
            case 0:
                jobmr = 'B';
                break;
            case 1:
                jobmr = 'F';
                break;
            default:
                error ("__sl_sb16cd__: argument jobmr invalid");
        }


        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        int ldd;
        
        if (jobd == 'Z')
            ldd = 1;
        else
            ldd = max (1, p);

        int ldf = max (1, m);
        int ldg = max (1, n);

        // arguments out
        ColumnVector hsv (n);

        // workspace
        int liwork;

        if (jobmr == 'B')
            liwork = 0;
        else                // if JOBMR = 'F'
            liwork = n;

        int ldwork;
        int mp;

        if (jobcf == 'L')
            mp = m;
        else                // if JOBCF = 'R'
            mp = p;

        ldwork = 2*n*n + max (1, 2*n*n + 5*n, n*max(m,p),
                              n*(n + max(n,mp) + min(n,mp) + 6));


        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int iwarn = 0;
        int info = 0;


        // SLICOT routine SB16CD
        F77_XFCN (sb16cd, SB16CD,
                 (dico, jobd, jobmr, jobcf,
                  ordsel,
                  n, m, p,
                  ncr,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  f.fortran_vec (), ldf,
                  g.fortran_vec (), ldg,
                  hsv.fortran_vec (),
                  tol,
                  iwork,
                  dwork, ldwork,
                  iwarn, info));


        if (f77_exception_encountered)
            error ("fwcfconred: exception in SLICOT subroutine SB16CD");

        static const char* err_msg[] = {
            "0: OK",
            "1: eigenvalue computation failure",
            "2: the matrix A-L*C is not stable",
            "3: the matrix A-B*F is not stable",
            "4: the Lyapunov equation for computing the "
                "observability Grammian is (nearly) singular",
            "5: the Lyapunov equation for computing the "
                "controllability Grammian is (nearly) singular",
            "6: the computation of Hankel singular values failed"};

        static const char* warn_msg[] = {
            "0: OK",
            "1: with ORDSEL = 'F', the selected order NCR is "
                "greater than the order of a minimal realization "
                "of the controller.",
            "2: with ORDSEL = 'F', the selected order NCR "
                "corresponds to repeated singular values, which are "
                "neither all included nor all excluded from the "
                "reduced controller. In this case, the resulting NCR "
                "is set automatically to the largest value such that "
                "HSV(NCR) > HSV(NCR+1)."};

        error_msg ("fwcfconred", info, 6, err_msg);
        warning_msg ("fwcfconred", iwarn, 2, warn_msg);


        // resize
        a.resize (ncr, ncr);    // Ac
        g.resize (ncr, p);      // Bc
        f.resize (m, ncr);      // Cc
        // Dc = 0
        
        // return values
        retval(0) = a;
        retval(1) = g;
        retval(2) = f;
        retval(3) = octave_value (ncr);
        retval(4) = hsv;
    }
    
    return retval;
}
