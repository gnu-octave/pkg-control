/*

Copyright (C) 2012   Lukas F. Reichlin

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

SLICOT system identification
Uses SLICOT IB01AD, IB01BD and IB01CD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: March 2012
Version: 0.1

*/

#include <octave/oct.h>
#include <octave/f77-fcn.h>
#include <octave/Cell.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (ib01ad, IB01AD)
                 (char& METH, char& ALG, char& JOBD,
                  char& BATCH, char& CONCT, char& CTRL,
                  int& NOBR, int& M, int& L,
                  int& NSMP,
                  double* U, int& LDU,
                  double* Y, int& LDY,
                  int& N,
                  double* R, int& LDR,
                  double* SV,
                  double& RCOND, double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& IWARN, int& INFO);
}

// PKG_ADD: autoload ("slib01ad", "control_slicot_functions.oct");
DEFUN_DLD (slib01ad, args, nargout,
   "-*- texinfo -*-\n\
Slicot IB01AD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 10)
    {
        print_usage ();
    }
    else
    {
////////////////////////////////////////////////////////////////////////////////////
//      SLICOT IB01AD - preprocess the input-output data                          //
////////////////////////////////////////////////////////////////////////////////////

        // arguments in
        char meth_a;
        char meth_b;
        char alg;
        char jobd;
        char batch;
        char conct;
        char ctrl = 'N';
        
        const Cell y_cell = args(0).cell_value ();
        const Cell u_cell = args(1).cell_value ();
        int nobr = args(2).int_value ();
        int nuser = args(3).int_value ();
        
        const int imeth = args(4).int_value ();
        const int ialg = args(5).int_value ();
        const int iconct = args(6).int_value ();
        const int ictrl = args(7).int_value ();     // ignored
        
        double rcond = args(8).double_value ();
        double tol_a = args(9).double_value ();

        double tol_b = rcond;
        double tol_c = rcond;
        
            
        switch (imeth)
        {
            case 0:
                meth_a = 'M';
                meth_b = 'M';
                break;
            case 1:
                meth_a = 'N';
                meth_b = 'N';
                break;
            case 2:
                meth_a = 'N';    // no typo here
                meth_b = 'C';
                break;
            default:
                error ("slib01ad: argument 'meth' invalid");
        }

        switch (ialg)
        {
            case 0:
                alg = 'C';
                break;
            case 1:
                alg = 'F';
                break;
            case 2:
                alg = 'Q';
                break;
            default:
                error ("slib01ad: argument 'alg' invalid");
        }

        if (meth_a == 'M')
            jobd = 'M';
        else                    // meth_a == 'N'
            jobd = 'N';         // IB01AD.f says: This parameter is not relevant for METH = 'N'

        if (iconct == 0)
            conct = 'C';
        else
            conct = 'N';
/*
        if (ictrl == 0)
            ctrl = 'C';
        else
            ctrl = 'N';
*/
        // m and l are equal for all experiments, checked by iddata class
        int n_exp = y_cell.nelem ();            // number of experiments
        int m = u_cell.elem(0).columns ();      // m: number of inputs
        int l = y_cell.elem(0).columns ();      // l: number of outputs
        int nsmpl = 0;                          // total number of samples

        // arguments out
        int n;
        int ldr;
        
        if (meth_a == 'M' && jobd == 'M')
            ldr = max (2*(m+l)*nobr, 3*m*nobr);
        else if (meth_a == 'N' || (meth_a == 'M' && jobd == 'N'))
            ldr = 2*(m+l)*nobr;
        else
            error ("slib01ad: could not handle 'ldr' case");
        
        Matrix r (ldr, 2*(m+l)*nobr);
        ColumnVector sv (l*nobr);


        // repeat for every experiment in the dataset
        for (int i = 0; i < n_exp; i++)
        {
            if (n_exp == 1)
                batch = 'O';        // one block only
            else if (i == 0)
                batch = 'F';        // first block
            else if (i == n_exp-1)
                batch = 'L';        // last block
            else
                batch = 'I';        // intermediate block
      
            Matrix y = y_cell.elem(i).matrix_value ();
            Matrix u = u_cell.elem(i).matrix_value ();

            // y.rows == u.rows  is checked by iddata class
            // int m = u.columns ();   // m: number of inputs
            // int l = y.columns ();   // l: number of outputs
            int nsmp = y.rows ();   // nsmp: number of samples in the current experiment
            nsmpl += nsmp;          // nsmpl: total number of samples of all experiments

            // minimal nsmp size checked by __slicot_identification__.m
            if (batch == 'O')
            {
                if (nsmp < 2*(m+l+1)*nobr - 1)
                    error ("slident: require NSMP >= 2*(M+L+1)*NOBR - 1");
            }
            else
            {
                if (nsmp < 2*nobr)
                    error ("slident: require NSMP >= 2*NOBR");
            }
        
            int ldu;
        
            if (m == 0)
                ldu = 1;
            else                    // m > 0
                ldu = nsmp;

            int ldy = nsmp;

            // workspace
            int liwork_a;

            if (meth_a == 'N')            // if METH = 'N'
                liwork_a = (m+l)*nobr;
            else if (alg == 'F')        // if METH = 'M' and ALG = 'F'
                liwork_a = m+l;
            else                        // if METH = 'M' and ALG = 'C' or 'Q'
                liwork_a = 0;

            // TODO: Handle 'k' for DWORK

            int ldwork_a;
            int ns = nsmp - 2*nobr + 1;
        
            if (alg == 'C')
            {
                if (batch == 'F' || batch == 'I')
                {
                    if (conct == 'C')
                        ldwork_a = (4*nobr-2)*(m+l);
                    else    // (conct == 'N')
                        ldwork_a = 1;
                }
                else if (meth_a == 'M')   // && (batch == 'L' || batch == 'O')
                {
                    if (conct == 'C' && batch == 'L')
                        ldwork_a = max ((4*nobr-2)*(m+l), 5*l*nobr);
                    else if (jobd == 'M')
                        ldwork_a = max ((2*m-1)*nobr, (m+l)*nobr, 5*l*nobr);
                    else    // (jobd == 'N')
                        ldwork_a = 5*l*nobr;
                }
                else    // meth_b == 'N' && (batch == 'L' || batch == 'O')
                {
                    ldwork_a = 5*(m+l)*nobr + 1;
                }
            }
            else if (alg == 'F')
            {
                if (batch != 'O' && conct == 'C')
                    ldwork_a = (m+l)*2*nobr*(m+l+3);
                else if (batch == 'F' || batch == 'I')  // && conct == 'N'
                    ldwork_a = (m+l)*2*nobr*(m+l+1);
                else    // (batch == 'L' || '0' && conct == 'N')
                    ldwork_a = (m+l)*4*nobr*(m+l+1)+(m+l)*2*nobr;
            }
            else    // (alg == 'Q')
            {
                // int ns = nsmp - 2*nobr + 1;
                
                if (ldr >= ns && batch == 'F')
                {
                    ldwork_a = 4*(m+l)*nobr;
                }
                else if (ldr >= ns && batch == 'O')
                {
                    if (meth_a == 'M')
                        ldwork_a = max (4*(m+l)*nobr, 5*l*nobr);
                    else    // (meth == 'N')
                        ldwork_a = 5*(m+l)*nobr + 1;
                }
                else if (conct == 'C' && (batch == 'I' || batch == 'L'))
                {
                    ldwork_a = 4*(nobr+1)*(m+l)*nobr;
                }
                else    // if ALG = 'Q', (BATCH = 'F' or 'O', and LDR < NS), or (BATCH = 'I' or 'L' and CONCT = 'N')
                {
                    ldwork_a = 6*(m+l)*nobr;
                }
            }

            /*
            IB01AD.f Lines 438-445
            C     FURTHER COMMENTS
            C
            C     For ALG = 'Q', BATCH = 'O' and LDR < NS, or BATCH <> 'O', the
            C     calculations could be rather inefficient if only minimal workspace
            C     (see argument LDWORK) is provided. It is advisable to provide as
            C     much workspace as possible. Almost optimal efficiency can be
            C     obtained for  LDWORK = (NS+2)*(2*(M+L)*NOBR),  assuming that the
            C     cache size is large enough to accommodate R, U, Y, and DWORK.
            */

            ldwork_a = max (ldwork_a, (ns+2)*(2*(m+l)*nobr));

            /*
            IB01AD.f Lines 291-195:
            c             the workspace used for alg = 'q' is
            c                       ldrwrk*2*(m+l)*nobr + 4*(m+l)*nobr,
            c             where ldrwrk = ldwork/(2*(m+l)*nobr) - 2; recommended
            c             value ldrwrk = ns, assuming a large enough cache size.
            c             for good performance,  ldwork  should be larger.

            somehow ldrwrk and ldwork must have been mixed up here

            */


            OCTAVE_LOCAL_BUFFER (int, iwork_a, liwork_a);
            OCTAVE_LOCAL_BUFFER (double, dwork_a, ldwork_a);
        
            // error indicators
            int iwarn_a = 0;
            int info_a = 0;


            // SLICOT routine IB01AD
            F77_XFCN (ib01ad, IB01AD,
                     (meth_a, alg, jobd,
                      batch, conct, ctrl,
                      nobr, m, l,
                      nsmp,
                      u.fortran_vec (), ldu,
                      y.fortran_vec (), ldy,
                      n,
                      r.fortran_vec (), ldr,
                      sv.fortran_vec (),
                      rcond, tol_a,
                      iwork_a,
                      dwork_a, ldwork_a,
                      iwarn_a, info_a));


            if (f77_exception_encountered)
                error ("ident: exception in SLICOT subroutine IB01AD");

            static const char* err_msg[] = {
                "0: OK",
                "1: a fast algorithm was requested (ALG = 'C', or 'F') "
                    "in sequential data processing, but it failed; the "
                    "routine can be repeatedly called again using the "
                    "standard QR algorithm",
                "2: the singular value decomposition (SVD) algorithm did "
                    "not converge"};

            static const char* warn_msg[] = {
                "0: OK",
                "1: the number of 100 cycles in sequential data "
                    "processing has been exhausted without signaling "
                    "that the last block of data was get; the cycle "
                    "counter was reinitialized",
                "2: a fast algorithm was requested (ALG = 'C' or 'F'), "
                    "but it failed, and the QR algorithm was then used "
                    "(non-sequential data processing)",
                "3: all singular values were exactly zero, hence  N = 0 "
                    "(both input and output were identically zero)",
                "4: the least squares problems with coefficient matrix "
                    "U_f,  used for computing the weighted oblique "
                    "projection (for METH = 'N'), have a rank-deficient "
                    "coefficient matrix",
                "5: the least squares problem with coefficient matrix "
                    "r_1  [6], used for computing the weighted oblique "
                    "projection (for METH = 'N'), has a rank-deficient "
                    "coefficient matrix"};


            error_msg ("ident: IB01AD", info_a, 2, err_msg);
            warning_msg ("ident: IB01AD", iwarn_a, 5, warn_msg);
        }


        // resize
        int rs = 2*(m+l)*nobr;
        r.resize (rs, rs);

        
        // return values
        retval(0) = sv;
        retval(1) = octave_value (n);
    }
    
    return retval;
}
