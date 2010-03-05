## Copyright (C) 2009   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} mconnect (@var{sys}, @var{m})
## @deftypefnx {Function File} {@var{sys} =} mconnect (@var{sys}, @var{m}, @var{inputs}, @var{outputs})
## Arbitrary interconnections between the inputs and outputs of an LTI model.
## @example
## @group
## Solve the system equations of
## y(t) = G e(t)
## e(t) = u(t) + M y(t)
## in order to build
## y(t) = H u(t)
## The matrix M for a (p-by-m) system G has m rows and p columns.
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function sys = mconnect (sys, M, in_idx, out_idx)

  if (nargin != 2 && nargin != 4)
    print_usage ();
  endif

  [p, m] = size (sys);
  [mrows, mcols] = size (M);

  if (p != mcols || m != mrows)
    error ("mconnect: second argument must be a (%dx&d) matrix", m, p);
  endif

  if (! isreal (M))
    error ("mconnect: second argument must be a matrix with real coefficients");
  endif

  sys = __sysconnect__ (sys, M);

  if (nargin == 4)
    sys = __sysprune__ (sys, out_idx, in_idx);
  endif

endfunction