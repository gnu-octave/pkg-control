/*

Copyright (C) 2009 - 2010   Lukas F. Reichlin

This file is part of LTI Syncope.

LTI Syncope is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LTI Syncope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Solution of Lyapunov equations.
Uses SLICOT SB03MD by courtesy of NICONET e.V.
<http://www.slicot.org>

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: December 2009
Version: 0.2

*/

#include <octave/oct.h>
#include <f77-fcn.h>

extern "C"
{ 
    int F77_FUNC (sb03md, SB03MD)
                 (char& DICO, char& JOB,
                  char& FACT, char& TRANA,
                  int& N,
                  double* A, int& LDA,
                  double* U, int& LDU,
                  double* C, int& LDC,
                  double& SCALE,
                  double& SEP, double& FERR,
                  double* WR, double* WI,
                  int* IWORK,
                  double* DWORK, int& LDWORK,
                  int& INFO);
}

int max (int a, int b)
{
    if (a > b)
        return a;
    else
        return b;
}
     
DEFUN_DLD (slsb03md, args, nargout, "Slicot SB03MD Release 5.0")
{
    int nargin = args.length ();
    octave_value_list retval;
    
    if (nargin != 3)
    {
        print_usage ();
    }
    else
    {
        // arguments in
        char dico;
        char job = 'X';
        char fact = 'N';
        char trana = 'T';
        
        NDArray a = args(0).array_value ();
        NDArray c = args(1).array_value ();
        int dt = args(2).int_value ();
        
        if (dt == 0)
          dico = 'C';
        else
          dico = 'D';
        
        int n = a.rows ();      // n: number of states
        
        int lda = max (1, n);
        int ldu = max (1, n);
        int ldc = max (1, n);
        
        // arguments out
        double scale;
        double sep = 0;
        double ferr = 0;
        
        dim_vector dv_u (2);
        dv_u(0) = ldu;
        dv_u(1) = n;
        
        dim_vector dv (1);
        dv(0) = n;
        
        NDArray u (dv_u);
        NDArray wr (dv);
        NDArray wi (dv);
        
        // workspace
        int* iwork = 0;  // not referenced because job = X
        
        int ldwork = max (n*n, 3*n);
        OCTAVE_LOCAL_BUFFER (double, dwork, ldwork);

        // error indicator
        int info;
        

        // SLICOT routine SB03MD
        F77_XFCN (sb03md, SB03MD,
                 (dico, job,
                  fact, trana,
                  n,
                  a.fortran_vec (), lda,
                  u.fortran_vec (), ldu,
                  c.fortran_vec (), ldc,
                  scale,
                  sep, ferr,
                  wr.fortran_vec (), wi.fortran_vec (),
                  iwork,
                  dwork, ldwork,
                  info));

        if (f77_exception_encountered)
            error ("lyap: slsb03md: exception in SLICOT subroutine SB03MD");

        if (info != 0)
            error ("lyap: slsb03md: SB03MD returned info = %d", info);
        
        // return values
        retval(0) = c;
        retval(1) = octave_value (scale);
    }
    
    return retval;
}
