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

Controller reduction based on Balance & Truncate (B&T) or
Singular Perturbation Approximation (SPA) method.
Uses SLICOT SB16AD by courtesy of NICONET e.V.
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
    int F77_FUNC (sb16ad, SB16AD)
                 (char& DICO, char& JOBC, char& JOBO, char& JOBMR,
                  char& WEIGHT, char& EQUIL, char& ORDSEL,
                  int& N, int& M, int& P,
                  int& NC, int& NCR,
                  double& ALPHA,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  double* AC, int& LDAC,
                  double* BC, int& LDBC,
                  double* CC, int& LDCC,
                  double* DC, int& LDDC,
                  int& NCS,
                  double* HSVC,
                  double& TOL1, double& TOL2,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}

// PKG_ADD: autoload ("__sl_sb16ad__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_sb16ad__, args, nargout,
   "-*- texinfo -*-\n\
Slicot SB16AD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 19)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char jobc;
        char jobo;
        char jobmr;
        char weight;
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
        double alpha = args(8).double_value ();
        const int ijobmr = args(9).int_value ();
                       
        Matrix ac = args(10).matrix_value ();
        Matrix bc = args(11).matrix_value ();
        Matrix cc = args(12).matrix_value ();
        Matrix dc = args(13).matrix_value ();
        
        const int iweight = args(14).int_value ();
        const int ijobc = args(15).int_value ();
        const int ijobo = args(16).int_value ();

        double tol1 = args(17).double_value ();
        double tol2 = args(18).double_value ();

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
                error ("__sl_sb16ad__: argument jobmr invalid");
        }

        switch (iweight)
        {
            case 0:
                weight = 'N';
                break;
            case 1:
                weight = 'O';
                break;
            case 2:
                weight = 'I';
                break;
            case 3:
                weight = 'P';
                break;
            default:
                error ("__sl_sb16ad__: argument weight invalid");
        }

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs
        
        int nc = ac.rows ();

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        int ldd = max (1, p);

        int ldac = max (1, nc);
        int ldbc = max (1, nc);
        int ldcc = max (1, m);
        int lddc = max (1, m);

        // arguments out
        int ncs;
        ColumnVector hsvc (n);

        // workspace
        int liwork;
        int liwrk1;
        int liwrk2;

        switch (jobmr)
        {
            case 'B':
                liwrk1 = 0;
                break;
            case 'F':
                liwrk1 = nc;
                break;
            default:
                liwrk1 = 2*nc;
        }

        if (weight == 'N')
            liwrk2 = 0;
        else
            liwrk2 = 2*(m+p);

        liwork = max (1, liwrk1, liwrk2);

        int ldwork;
        int lfreq;
        int lsqred;

        if (weight == 'N')
        {
            if (equil == 'N')           // if WEIGHT = 'N' and EQUIL = 'N'
                lfreq  = nc*(max (m, p) + 5);
            else                        // if WEIGHT = 'N' and EQUIL  = 'S'
                lfreq  = max (n, nc*(max (m, p) + 5));

        }
        else                            // if WEIGHT = 'I' or 'O' or 'P'
        {
            lfreq = (n+nc)*(n+nc+2*m+2*p) +
                     max ((n+nc)*(n+nc+max(n+nc,m,p)+7), (m+p)*(m+p+4));
        }

        lsqred = max (1, 2*nc*nc+5*nc);
        ldwork = 2*nc*nc + max (1, lfreq, lsqred);

        OCTAVE_LOCAL_BUFFER (int, iwork, liwork);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);
        
        // error indicators
        int iwarn = 0;
        int info = 0;


        // SLICOT routine SB16AD
        F77_XFCN (sb16ad, SB16AD,
                 (dico, jobc, jobo, jobmr,
                  weight, equil, ordsel,
                  n, m, p,
                  nc, ncr,
                  alpha,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  ac.fortran_vec (), ldac,
                  bc.fortran_vec (), ldbc,
                  cc.fortran_vec (), ldcc,
                  dc.fortran_vec (), lddc,
                  ncs,
                  hsvc.fortran_vec (),
                  tol1, tol2,
                  iwork,
                  dwork, ldwork,
                  iwarn, info));


        if (f77_exception_encountered)
            error ("conred: exception in SLICOT subroutine SB16AD");

        static const char* err_msg[] = {
            "0: OK",
            "1: the closed-loop system is not well-posed; "
                "its feedthrough matrix is (numerically) singular",
            "2: the computation of the real Schur form of the "
                "closed-loop state matrix failed",
            "3: the closed-loop state matrix is not stable",
            "4: the solution of a symmetric eigenproblem failed",
            "5: the computation of the ordered real Schur form "
                "of Ac failed",
            "6: the separation of the ALPHA-stable/unstable "
                "diagonal blocks failed because of very close eigenvalues",
            "7: the computation of Hankel singular values failed"};

        static const char* warn_msg[] = {
            "0: OK",
            "1: with ORDSEL = 'F', the selected order NCR is greater "
                "than NSMIN, the sum of the order of the "
                "ALPHA-unstable part and the order of a minimal "
                "realization of the ALPHA-stable part of the given "
                "controller; in this case, the resulting NCR is set "
                "equal to NSMIN.",
            "2: with ORDSEL = 'F', the selected order NCR "
                "corresponds to repeated singular values for the "
                "ALPHA-stable part of the controller, which are "
                "neither all included nor all excluded from the "
                "reduced model; in this case, the resulting NCR is "
                "automatically decreased to exclude all repeated "
                "singular values.",
            "3: with ORDSEL = 'F', the selected order NCR is less "
                "than the order of the ALPHA-unstable part of the "
                "given controller. In this case NCR is set equal to "
                "the order of the ALPHA-unstable part."};

        error_msg ("conred", info, 7, err_msg);
        warning_msg ("conred", iwarn, 3, warn_msg);


        // resize
        ac.resize (ncr, ncr);
        bc.resize (ncr, p);    // p: number of plant outputs
        cc.resize (m, ncr);    // m: number of plant inputs
        hsvc.resize (ncs);
        
        // return values
        retval(0) = ac;
        retval(1) = bc;
        retval(2) = cc;
        retval(3) = dc;
        retval(4) = octave_value (ncr);
        retval(5) = hsvc;
        retval(6) = octave_value (ncs);
    }
    
    return retval;
}
