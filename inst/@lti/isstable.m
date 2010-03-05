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
## @deftypefn {Function File} {@var{bool} =} isstable (@var{sys})
## @deftypefnx {Function File} {@var{bool} =} isstable (@var{sys}, @var{tol})
## Determine whether LTI system is stable.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function bool = isstable (sys, tol = 0)

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif

  eigw = pole (sys);

  if (isct (sys))
    bool = all (real (eigw) < -tol*(1 + abs (eigw)));
  else
    bool = all (abs (eigw) < 1 - tol);
  endif

endfunction