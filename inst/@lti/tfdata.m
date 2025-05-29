# Copyright (C) 2009-2016   Lukas F. Reichlin
# Copyright (C) 2025        Torsten Lilge
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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{sys})
## @deftypefnx {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{sys}, @var{"vector"})
## @deftypefnx {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{sys}, @var{"tfpoly"})
##
## Access transfer function data.
## Argument @var{sys} is not limited to transfer function models.
## If @var{sys} is not a transfer function, it is converted automatically.
## @var{sys} can also be a real-valued gain matrix which is then
## interpreted as static gain system.
##
## @strong{Inputs}
## @table @var
## @item sys
## Any type of @acronym{LTI} model or a real-valued matrix which is
## interpreted as continous-time static gain.
## @item "v", "vector"
## For SISO models, return @var{num} and @var{den} directly as column vectors
## instead of cells containing a single column vector.
## @item "t", "tfpoly"
## Return the polynomials as tfpoly class.
## @end table
##
## @strong{Outputs}
## @table @var
## @item num
## Cell of numerator(s).  Each numerator is a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## num@{i,j@} contains the numerator polynomial from input j to output i.
## In the SISO case, a single vector is possible as well.
## @item den
## Cell of denominator(s).  Each denominator is a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## den@{i,j@} contains the denominator polynomial from input j to output i.
## In the SISO case, a single vector is possible as well.
## @item tsam
## Sampling time in seconds. If @var{sys} is a continuous-time model,
## a zero is returned.
## @end table
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>

function [num, den, tsam] = tfdata (sys, rtype = "cell")

  if (! isa (sys, 'lti'))
    if (! is_real_matrix (sys))
      error (["tfdata: has to be called with an @lti object ",...
                      "or with a real matrix (static gain)\n"]);
    else
      sys = tf (sys);
    endif
  endif

  if (! isa (sys, "tf"))
    sys = tf (sys);
  endif

  [num, den] = __sys_data__ (sys);

  tsam = sys.tsam;

  if (! strncmpi (rtype, "t", 1))                # != tfpoly
    [num, den] = __make_tf_polys_equally_long__ (tf (num,den));
    num = cellfun (@get, num, "uniformoutput", false);
    den = cellfun (@get, den, "uniformoutput", false);
  endif

  if (strncmpi (rtype, "v", 1) && issiso (sys))  # == vector
    num = num{1};
    den = den{1};
  endif

endfunction
