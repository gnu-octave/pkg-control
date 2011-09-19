## Copyright (C) 2009, 2010, 2011   Lukas F. Reichlin
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
## @deftypefnx {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{sys}, @var{"tfpoly"})
## Access transfer function data.
## Argument @var{sys} is not limited to transfer function models.
## If @var{sys} is not a transfer function, it is converted automatically.
##
## @strong{Inputs}
## @table @var
## @item sys
## Any type of LTI model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item num
## Cell of numerator(s).  Each numerator is a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## @item den
## Cell of denominator(s).  Each denominator is a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## @item tsam
## Sampling time in seconds.  If @var{sys} is a continuous-time model,
## a zero is returned.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function [num, den, tsam] = tfdata (sys, rtype = "vector")

  if (! isa (sys, "tf"))
    sys = tf (sys);
  endif

  [num, den] = __sys_data__ (sys); 

  tsam = sys.tsam;

  if (rtype(1) == "v")
    num = cellfun (@get, num, "uniformoutput", false);
    den = cellfun (@get, den, "uniformoutput", false);
  endif

endfunction
