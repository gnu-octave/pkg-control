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

Transfer matrix of a given state-space representation (A,B,C,D),
using the pole-zeros method.
Uses SLICOT TB04BD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2010
Version: 0.3

*/

#include <octave/oct.h>
#include "common.h"

extern "C"
{ 
    int F77_FUNC (tb04bd, TB04BD)
                 (char& JOBD, char& ORDER, char& EQUIL,
                  F77_INT& N, F77_INT& M, F77_INT& P, F77_INT& MD,
                  double* A, F77_INT& LDA,
                  double* B, F77_INT& LDB,
                  double* C, F77_INT& LDC,
                  double* D, F77_INT& LDD,
                  F77_INT* IGN, F77_INT& LDIGN,
                  F77_INT* IGD, F77_INT& LDIGD,
                  double* GN, double* GD,
                  double& TOL,
                  F77_INT* IWORK,
                  double* DWORK, F77_INT& LDWORK,
                  F77_INT& INFO);
}

// PKG_ADD: autoload ("__sl_tb04bd__", "__control_slicot_functions__.oct");    
DEFUN_DLD (__sl_tb04bd__, args, nargout,
   "-*- texinfo -*-\n\
Slicot TB04BD Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_idx_type nargin = args.length ();
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
        const F77_INT scaled = args(4).int_value ();

        if (scaled == 0)
            equil = 'S';
        else
            equil = 'N';

        F77_INT n = TO_F77_INT (a.rows ());      // n: number of states
        F77_INT m = TO_F77_INT (b.columns ());   // m: number of inputs
        F77_INT p = TO_F77_INT (c.rows ());      // p: number of outputs
        F77_INT md = n + 1;

        F77_INT lda = max (1, n);
        F77_INT ldb = max (1, n);
        F77_INT ldc = max (1, p);
        F77_INT ldd = max (1, p);

        // arguments out
        F77_INT ldign = max (1, p);
        F77_INT ldigd = max (1, p);
        F77_INT lg = p * m * md;

        OCTAVE_LOCAL_BUFFER (F77_INT, ign, ldign*m);
        OCTAVE_LOCAL_BUFFER (F77_INT, igd, ldigd*m);

        RowVector gn (lg);
        RowVector gd (lg);

        // tolerance
        double tol = 0;  // use default value

        // workspace
        F77_INT ldwork = max (1, n*(n + p) + max (n + max (n, p), n*(2*n + 5)));

        OCTAVE_LOCAL_BUFFER (F77_INT, iwork, n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        F77_INT info;


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
            error ("ss2tf: __sl_tb04bd__: exception in SLICOT subroutine TB04BD");

        if (info != 0)
            error ("ss2tf: __sl_tb04bd__: TB04BD returned info = %d", info);


        // return values
        Cell num(p, m);
        Cell den(p, m);
        F77_INT ik, istr;

        for (ik = 0; ik < p*m; ik++)
        {
            istr = ik*md;
            num.xelem (ik) = gn.extract_n (istr, ign[ik]+1);
            den.xelem (ik) = gd.extract_n (istr, igd[ik]+1);
        }
       
        retval(0) = num;
        retval(1) = den; 
    }

    return retval;
}
