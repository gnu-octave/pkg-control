## Copyright (C) 2025 Torten Lilge
##
## This function is part of the GNU Octave Control Package
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} [@var{T}, @var{S}, @var{BLKSZ}] = bdschur (@var{A})
## @deftypefnx {Function File} [@var{T}, @var{S}, @var{BLKSZ}] = bdschur (@var{A}, @var{condmax})
##
## Compute a block-diagonal Schur decomposition.
##
## @strong{Inputs}
## @table @var
## @item A
## Matrix to be decomposed.
## @item condmax
## Maximum condition number (real scalar value greater or equal 1) of the
## resulting transformation matrix @var{T}. This parameter is optional.
## If not provided, a value of 1e4 is used.
## @end table
##
## @strong{Outputs}
## @table @var
## @item T
## Transformation matrix.
## @item S
## Block diagonalized matrix.
## @item BLKSZ
## Array with sizes of the blocks in @var{S}.
## @end table
##
## The resulting matrix @var{S} in block diagonal form is given by
## @example
## @group
##          -1
##     S = T  A T
## @end group
## @end example
##
## @end deftypefn

function [T, S, blksz] = bdschur (A, condmax = 1e4)

  if (nargin < 1) || (nargin > 2)
    print_usage ();
  endif

  if ! (is_real_matrix (A) && (size (A,1) == size (A,2)))
    error ("bdschur: first argument must be a real square matrix");
  endif

  if ! (is_real_scalar (condmax) && (condmax >= 1))
    error ("bdschur: second argument mus be a real scalar greater or equal 1");
  endif

  # Calculate Schur matrix first and from this the block diagonal schur matrix
  [T, S] = schur (A);
  [T, S, blksz] = __sl_mb03rd__ (T, S, condmax);

endfunction
