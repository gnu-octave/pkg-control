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

Transfer matrix of a given state-space representation (A,B,C,D),
using the pole-zeros method.
Uses SLICOT TB04BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2010
Version: 0.1

*/

#include <octave/oct.h>
#include <f77-fcn.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb04bd, TB04BD)
                 (char& JOBD, char& ORDER, char& EQUIL,
                  int& N, int& M, int& P, int& MD,
                  double* A, int& LDA,
                  double* B, int& LDB,
                  double* C, int& LDC,
                  double* D, int& LDD,
                  int* IGN, int& LDIGN,
                  int* IGD, int& LDIGD,
                  double* GN, double* GD,
                  double& TOL,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

// PKG_ADD: autoload ("__sl_tb04bd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tb04bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TB04BD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    int nargin = args.length ();
    octave_value_list retval;

    if (nargin != 5)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char jobd = 'D';
        char order = 'D';
        char equil;

        Matrix a = args(0).matrix_value ();
        Matrix b = args(1).matrix_value ();
        Matrix c = args(2).matrix_value ();
        Matrix d = args(3).matrix_value ();
        const int scaled = args(4).int_value ();

        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        int n = a.rows ();      // n: number of states
        int m = b.columns ();   // m: number of inputs
        int p = c.rows ();      // p: number of outputs
        int md = n + 1;

        int lda = max (1, n);
        int ldb = max (1, n);
        int ldc = max (1, p);
        int ldd = max (1, p);

        // arguments out
        int ldign = max (1, p);
        int ldigd = max (1, p);
        int lg = p * m * md;

        OCTAVE_LOCAL_BUFFER (int, ign, ldign*m);
        OCTAVE_LOCAL_BUFFER (int, igd, ldigd*m);
        
        Matrix ignm (ldign, m);
        Matrix igdm (ldigd, m);

        RowVector gn (lg);
        RowVector gd (lg);

        // tolerance
        double tol = 0;  // use default value

        // workspace
        int ldwork = max (1, n*(n + p) + max (n + max (n, p), n*(2*n + 5)));

        OCTAVE_LOCAL_BUFFER (int, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        int info;


        // SLICOT routine TB04BD
        F77_XFCN (tb04bd, TB04BD,
                 (jobd, order, equil,
                  n, m, p, md,
                  a.fortran_vec (), lda,
                  b.fortran_vec (), ldb,
                  c.fortran_vec (), ldc,
                  d.fortran_vec (), ldd,
                  ign, ldign,
                  igd, ldigd,
                  gn.fortran_vec (), gd.fortran_vec (),
                  tol,
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("are: __sl_tb04bd__: exception in SLICOT subroutine TB04BD");

        if (info != 0)
            error ("are: __sl_tb04bd__: TB04BD returned info = %d", info);

        for (octave_idx_type i = 0; i < ldign*m; i++)
            ignm.xelem (i) = ign[i];
        
        for (octave_idx_type i = 0; i < ldigd*m; i++)
            igdm.xelem (i) = igd[i];

        // return values
        retval(0) = gn;
        retval(1) = gd;
        retval(2) = ignm;
        retval(3) = igdm;
        retval(4) = octave_value (md);
        retval(5) = octave_value (p);
        retval(6) = octave_value (m);
    }

    return retval;
}
