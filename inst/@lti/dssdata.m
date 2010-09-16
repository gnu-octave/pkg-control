## Copyright (C) 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{a}, @var{b}, @var{c}, @var{d}, @var{e}, @var{tsam}] =} dssdata (@var{sys})
## Access descriptor state-space model data.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2010
## Version: 0.1

function [a, b, c, d, e, tsam] = dssdata (sys, flg = 0)

  ## NOTE: In case sys is not a dss model (matrice e empty),
  ##       dssdata (sys, []) returns e = [] whereas
  ##       dssdata (sys) returns e = eye (size (a))

  if (nargin > 2)
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    sys = ss (sys);
  endif

  [a, b, c, d, e] = __sys_data__ (sys);

  if (isempty (e) && ! isempty (flg))
    e = eye (size (a));  # return eye for ss models
  endif

  tsam = sys.tsam;

endfunction