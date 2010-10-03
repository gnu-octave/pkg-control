## Copyright (C) 2009 - 2010   Lukas F. Reichlin
##
## This file is part of LTI Syncope.
##
## LTI Syncope is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LTI Syncope is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{f} =} place (@var{sys}, @var{p})
## @deftypefnx {Function File} {@var{f} =} place (@var{a}, @var{b}, @var{p})
## @deftypefnx {Function File} {[@var{f}, @var{nfp}, @var{nap}, @var{nup}] =} place (@var{sys}, @var{p}, @var{alpha})
## @deftypefnx {Function File} {[@var{f}, @var{nfp}, @var{nap}, @var{nup}] =} place (@var{a}, @var{b}, @var{p}, @var{alpha})
## Pole assignment for a given matrix pair (A,B) such that eig (A-B*F) = P.
## If parameter alpha is specified, poles with real parts (continuous time)
## or moduli (discrete time) below alpha are left untouched.
## Uses SLICOT SB01BD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system.
## @item p
## Desired eigenvalues of the closed-loop system state-matrix A-B*F.
## length (P) <= rows (A)
## @item alpha
## Specifies the maximum admissible value, either for real
## parts or for moduli, of the eigenvalues of A which will
## not be modified by the eigenvalue assignment algorithm.
## ALPHA >= 0 for discrete-time systems.
## @end table
##
## @strong{Outputs}
## @table @var
## @item f
## State feedback gain matrix.
## @item nfp
## The number of fixed poles, i.e. eigenvalues of A having
## real parts less than ALPHA, or moduli less than ALPHA.
## These eigenvalues are not modified by place.
## @item nap
## The number of assigned eigenvalues. NAP = N-NFP-NUP.
## @item nup
## The number of uncontrollable eigenvalues detected by the
## eigenvalue assignment algorithm.
## @end table
##
## @example
## @group
## Place is also suitable to design estimator gains:
## L = place (A.', C.', p).'
## L = place (sys.', p).'   % useful for discrete-time systems
## @end group
## @end example
## @end deftypefn

## Special thanks to Peter Benner from TU Chemnitz for his advice.
## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: December 2009
## Version: 0.3

function [f, nfp, nap, nup] = place (a, b, p = [], alpha = [], tol = [])

  if (nargin < 2 || nargin > 5)
    print_usage ();
  endif

  if (isa (a, "lti"))              # place (sys, p), place (sys, p, alpha), place (sys, p, alpha, tol)
    if (nargin > 4)                # nargin < 2 already tested
      print_usage ();
    else
      tol = alpha;
      alpha = p;
      p = b;
      sys = a;
      [a, b] = dssdata (sys, []);  # descriptor matrice e should have no influence
      discrete = ! isct (sys);     # treat tsam = -1 as continuous system
    endif
  else  # place (a, b, p), place (a, b, p, alpha), place (a, b, p, alpha, tol)
    if (nargin < 3)                # nargin > 5 already tested
      print_usage ();
    else
      if (! is_real_square_matrix (a) || ! is_real_matrix (b) || rows (a) != rows (b))
        error ("place: matrices a and b not conformal");
      endif
      discrete = 0;                # assume continuous system
    endif
  endif

  if (! isnumeric (p) || ! isvector (p) || isempty (p))  # p could be complex
    error ("place: p must be a vector");
  endif
  
  p = sort (reshape (p, [], 1));   # complex conjugate pairs must appear together
  wr = real (p);
  wi = imag (p);
  
  n = rows (a);                    # number of states
  np = length (p);                 # number of given eigenvalues
  
  if (np > n)
    error ("place: at most %d eigenvalues can be assigned for the given matrix a (%dx%d)",
            n, n, n);
  endif

  if (isempty (alpha))
    if (discrete)
      alpha = 0;
    else
      alpha = - norm (a, inf);
    endif
  endif
  
  if (isempty (tol))
    tol = 0;
  endif

  [f, iwarn, nfp, nap, nup] = slsb01bd (a, b, wr, wi, discrete, alpha, tol);
  f = -f;                          # A + B*F --> A - B*F
  
  if (iwarn)
    warning ("place: %d violations of the numerical stability condition NORM(F) <= 100*NORM(A)/NORM(B)",
              iwarn);
  endif

endfunction


## Test from "legacy" control package 1.0.*
%!shared A, B, C, P, Kexpected
%! A = [0, 1; 3, 2];
%! B = [0; 1];
%! C = [2, 1];  # C is needed for ss; it doesn't matter what the value of C is
%! P = [-1, -0.5];
%! Kexpected = [3.5, 3.5];
%!assert (place (ss (A, B, C), P), Kexpected, 2*eps);
%!assert (place (A, B, P), Kexpected, 2*eps);

## FIXME: Test from SLICOT example SB01BD fails with 4 eigenvalues in P
%!shared F, F_exp, ev_ol, ev_cl
%! A = [-6.8000   0.0000  -207.0000   0.0000
%!       1.0000   0.0000     0.0000   0.0000
%!      43.2000   0.0000     0.0000  -4.2000
%!       0.0000   0.0000     1.0000   0.0000];
%!
%! B = [ 5.6400   0.0000
%!       0.0000   0.0000
%!       0.0000   1.1800
%!       0.0000   0.0000];
%!
%! P = [-0.5000 + 0.1500i
%!      -0.5000 - 0.1500i];
#%!      -2.0000 + 0.0000i
#%!      -0.4000 + 0.0000i];
%!
%! ALPHA = -0.4;
%! TOL = 1e-8;
%!
%! F = place (A, B, P, ALPHA, TOL);
%!
%! F_exp = - [-0.0876  -4.2138   0.0837 -18.1412
%!            -0.0233  18.2483  -0.4259  -4.8120];
%!
%! ev_ol = sort (eig (A));
%! ev_cl = sort (eig (A - B*F));
%!
%!assert (F, F_exp, 1e-4);
%!assert (ev_ol(3:4), ev_cl(3:4), 1e-4);
