## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{num}, @var{den}, @var{tsam}] =} filtdata (@var{sys})
## @deftypefnx {Function File} {[@var{num}, @var{den}, @var{tsam}] =} filtdata (@var{sys}, @var{"vector"})
##
## Access discrete-time transfer function data in DSP format.
## Argument @var{sys} is not limited to transfer function models.
## If @var{sys} is not a transfer function, it is converted automatically.
## @var{sys} can also be a real-valued gain matrix which is then
## interpreted as static gain system with unspecified sampling time (-1).
##
## @strong{Inputs}
## @table @var
## @item sys
## Any type of discrete-time @acronym{LTI} model or a real-valued matrix which is
## interpreted as discrete-time static gain with unspecified sampling time -1.
## @item "v", "vector"
## For SISO models, return @var{num} and @var{den} directly as column vectors
## instead of cells containing a single column vector.
## @end table
##
## @strong{Outputs}
## @table @var
## @item num
## Cell of numerator(s).  Each numerator is a row vector
## containing the coefficients of the polynomial in ascending powers of z^-1.
## num@{i,j@} contains the numerator polynomial from input j to output i.
## In the SISO case, a single vector is possible as well.
## @item den
## Cell of denominator(s).  Each denominator is a row vector
## containing the coefficients of the polynomial in ascending powers of z^-1.
## den@{i,j@} contains the denominator polynomial from input j to output i.
## In the SISO case, a single vector is possible as well.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified, -1 is returned.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2012
## Version: 0.1

function [num, den, tsam] = filtdata (sys, rtype = "cell")

  if (nargin > 2)
    print_usage ();
  endif

  if (! isa (sys, 'lti'))
    if (! is_real_matrix (sys))
      error (["filtdata: has to be called with an @lti object ",...
                        "or with a real matrix (static gain)\n"]);
    else
      sys = tf (sys, [], -1);
    endif
  endif

  if (! isdt (sys))
    error ("lti: filtdata: require discrete-time system");
  endif

  [num, den, tsam] = tfdata (sys);  # order of tf is the desired one

  ## make numerator and denominator polynomials equally long
  ## by adding leading zeros
  lnum = cellfun (@length, num, "uniformoutput", false);
  lden = cellfun (@length, den, "uniformoutput", false);

  lmax = cellfun (@max, lnum, lden, "uniformoutput", false);

  num = cellfun (@prepad, num, lmax, "uniformoutput", false);
  den = cellfun (@prepad, den, lmax, "uniformoutput", false);

  if (strncmpi (rtype, "v", 1) && issiso (sys))
    num = num{1};
    den = den{1};
  endif

endfunction
