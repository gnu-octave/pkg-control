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
## @deftypefn {Function File} {@var{strvec} =} strseq (@var{str}, @var{idx})
## Return a cell vector of indexed strings by appending the indices @var{idx}
## to the string @var{str}.
##
## @example
## strseq ("x", 1:3) = @{"x1"; "x2"; "x3"@}
## strseq ("u", [1, 2, 5]) = @{"u1"; "u2"; "u5"@}
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function strvec = strseq (str, idx)

  if (nargin != 2 || ! ischar (str) || ! isnumeric (idx))
    print_usage ();
  endif

  n = numel (idx);
  strvec = cell (n, 1);

  for k = 1 : n
    strvec{k, 1} = sprintf ("%s%d", str, idx(k));
  endfor

endfunction