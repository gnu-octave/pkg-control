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
Uses SLICOT SB16BD by courtesy of NICONET e.V.
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
    int F77_FUNC (sb16bd, SB16BD)
                 (char& DICO, char& JOBD, char& JOBMR, char& JOBCF,
                  char& EQUIL, char& ORDSEL,
                  int& N, int& M, int& P,
                  int& NCR,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* F, int& LDF,
                  double* G, int& LDG,
                  double* DC, int& LDDC,
                  double* HSV,
                  double& TOL1, double& TOL2,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}

// PKG_ADD: autoload ("__sl_sb16bd__", "__control_slicot_functions__.oct");         
DEFUN_DLD (__sl_sb16bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB16BD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 15)
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
        char equil;
        char ordsel;
        
        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        
        const int idico = args(4).int_value ();
        const int iequil = args(5).int_value ();
        int ncr = args(6).int_value ();
        const int iordsel = args(7).int_value ();
        const int ijobd = args(8).int_value ();
        const int ijobmr = args(9).int_value ();
                       
        Matrix f = args(10).matrix_value ();
        Matrix g = args(11).matrix_value ();

        const int ijobcf = args(12).int_value ();
        
        double tol1 = args(13).double_value ();
        double tol2 = args(14).double_value ();

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
            case 2:
                jobmr = 'S';
                break;
            case 3:
                jobmr = 'P';
                break;
            default:
                error ("__sl_sb16bd__: argument jobmr invalid");
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
        int lddc = max (1, m);

        // arguments out
        Matrix dc (lddc, p);
        ColumnVector hsv (n);

        // workspace
        int liwork;
        int pm;

        int ldwork;
        int lwr = max (1, n*(2*n+max(n,m+p)+5)+n*(n+1)/2);

        switch (jobmr)
        {
            case 'B':
                pm = 0;
                break;
            case 'F':
                pm = n;
                break;
            default:                // if JOBMR = 'S' or 'P'
                pm = max (1, 2*n);
        }

        if (ordsel == 'F' && ncr == n)
        {
            liwork = 0;
            ldwork = p*n;
        }
        else if (jobcf == 'L')
        {
            liwork = max (pm, m);
            ldwork = (n+m)*(m+p) + max (lwr, 4*m);
        }
        else                        // if JOBCF = 'R'
        {
            liwork = max (pm, p);
            ldwork = (n+p)*(m+p) + max (lwr, 4*p);
        }


        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int iwarn = 0;
        int info = 0;


        // SLICOT routine SB16BD
        F77_XFCN (sb16bd, SB16BD,
                 (dico, jobd, jobmr, jobcf,
                  equil, ordsel,
                  n, m, p,
                  ncr,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  f.fortran_vec (), ldf,
                  g.fortran_vec (), ldg,
                  dc.fortran_vec (), lddc,
                  hsv.fortran_vec (),
                  tol1, tol2,
                  iwork,
                  dwork, ldwork,
                  iwarn, info));


        if (f77_exception_encountered)
            error ("cfconred: exception in SLICOT subroutine SB16BD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the reduction of A-L*C to a real Schur form "
                "failed",
            "2: the matrix A-L*C is not stable (if DICO = 'C'), "
                "or not convergent (if DICO = 'D')",
            "3: the computation of Hankel singular values failed",
            "4: the reduction of A-B*F to a real Schur form "
                "failed",
            "5: the matrix A-B*F is not stable (if DICO = 'C'), "
                "or not convergent (if DICO = 'D')"};

        static const char* warn_msg[] = {
            "0: OK",
            "1: with ORDSEL = 'F', the selected order NCR is "
                "greater than the order of a minimal "
                "realization of the controller."};

        error_msg ("cfconred", info, 5, err_msg);
        warning_msg ("cfconred", iwarn, 1, warn_msg);


        // resize
        a.resize (ncr, ncr);    // Ac
        g.resize (ncr, p);      // Bc
        f.resize (m, ncr);      // Cc
        dc.resize (m, p);       // Dc
        
        // return values
        retval(0) = a;
        retval(1) = g;
        retval(2) = f;
        retval(3) = dc;
        retval(4) = octave_value (ncr);
        retval(5) = hsv;
    }
    
    return retval;
}
