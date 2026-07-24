## Copyright (C) 2026 J. Román <jdaniel.roman2004@gmail.com>
##
## This file is part of the control package for GNU Octave.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{csys}, @var{T}] =} compreal (@var{sys})
## @deftypefnx {Function File} {[@var{csys}, @var{T}] =} compreal (@var{sys}, @var{realization})
##
## Returns the state-space companion controllable ors observable form of the input
## system as well as the transformation matrix.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @item realization
## The @var{realization} argument selects the companion form:
## @table @code
## @item 'c'
## Controllable companion form, the default if @var{realization} is omitted.
## @item 'o'
## Observable companion form.
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item csys
## State-space model in companian form.
## @item T
## Transformation matix.
## @end table
##
## Given the state vector
## @tex
## $$x(t)$$
## @end tex
## @ifnottex
## @example
## @group
##       x(t)
## @end group
## @end example
## @end ifnottex
##
## We look for a transformation matrix such that
##
## @tex
## $$ x(t)=Tz(t) $$
## @end tex
## @ifnottex
## @example
## @group
##       x(t)=Tz(t)
## @end group
## @end example
## @end ifnottex
##
## Where @math{T} is our transformation matrix and $z(t)$ is our new state vector.
## Our new state-space is given by the following relationships:
##
## @tex
## $$ A_z = T^{-1}AT $$
## $$ B_z = T^{-1}B $$
## $$ C_z = CT $$
## $$ D_z = D $$
## @end tex
## @ifnottex
## @example
## @group
## A_z = inv(T) A T
## B_z = inv(T) B
## C_z = C T
## D_z = D
## @end group
## @end example
## @end ifnottex
##
## Alorithm based on
## @enumerate
## @item
## T. Kailath, Linear systems. Prentice-Hall Englewood Cliffs, NJ, 1980.
## ISBM 0-13-536961-4. Eequation (15a) and (15b)
## @end enumerate
##
## @end deftypefn

## Author: J. Román <jdaniel.roman2004@gmail.com>
## Created: May 2026

function [csys, T] = compreal (sys, realization)

  ## Validation

  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if (issiso (sys) == 0)
    error ("compreal: system must be SISO");
  endif

  if (nargin == 1 || strcmp (realization, "c"))
    realization = "controllable";
  elseif (strcmp (realization, "o"))
    realization = "observable";
  else
    error ("compreal: realization must be 'c' or 'o'");
  endif



  ## Matrix T

  sys = ss (sys);

  switch (realization)

    case "controllable"

      ctrmatrix = ctrb (sys.A, sys.B);

      if (rank (ctrmatrix) < rows (ctrmatrix))
        error ("compreal: system is not controllable");
      endif

      T = ctrmatrix;
      csys = ss2ss (sys, inv(T));

    case "observable"

      obsmatrix = obsv (sys.A, sys.C);

      if (rank (obsmatrix) < rows (obsmatrix))
        error ("compreal: system is not observable");
      endif

      T = obsmatrix;
      csys = ss2ss (sys, T);

  endswitch

endfunction

%!shared sys, A, B, C, D
%! A = [0 1; -2 -3];
%! B = [0; 1];
%! C = [1 0];
%! D = 0;
%! sys = ss (A, B, C, D);

%!test
%! [csys, T] = compreal (sys, 'c');
%! assert (tfdata(tf(csys)), tfdata(tf(sys)), 1e-6);

%!test
%! [csys, T] = compreal (sys, 'c');
%! assert (T, ctrb (A, B), 1e-10);

%!test
%! [csys, T] = compreal (sys, 'c');
%! assert (csys.A, inv (T) * A * T, 1e-10);

%!test
%! [csys, T] = compreal (sys, 'o');
%! assert (tfdata (tf (csys)), tfdata (tf (sys)), 1e-6);

%!error <realization must be>
%! compreal (sys, 'x');

%!error <system must be SISO>
%! sysmimo = ss (A, B, eye(2), zeros(2,1));
%! compreal (sysmimo);
