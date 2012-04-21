## Copyright (C) 2012   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} filt (@var{num}, @var{den}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} filt (@var{num}, @var{den}, @var{tsam}, @dots{})
## Create discrete transfer function model from data in DSP format.
##
## @strong{Inputs}
## @table @var
## @item num
## Numerator or cell of numerators.  Each numerator must be a row vector
## containing the coefficients of the polynomial in ascending powers of z^-1.
## num@{i,j@} contains the numerator polynomial from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item den
## Denominator or cell of denominators.  Each denominator must be a row vector
## containing the coefficients of the polynomial in ascending powers of z^-1.
## den@{i,j@} contains the denominator polynomial from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified,
## default value -1 (unspecified) is taken.
## @item @dots{}
## Optional pairs of properties and values.
## Type @command{set (tf)} for more information.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Discrete-time transfer function model.
## @end table
##
## @strong{Example}
## @example
## @group
##                 3 z^-1     
## H(z^-1) = -------------------
##           1 + 4 z^-1 + 2 z^-2
##
## octave:1> H = filt ([0, 3], [1, 4, 2])
## 
## Transfer function 'H' from input 'u1' to output ...
## 
##            3 z     
##  y1:  -------------
##       z^2 + 4 z + 2
## 
## Sampling time: unspecified
## Discrete-time model.
## @end group
## @end example
##
## @seealso{tf}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2012
## Version: 0.1

function sys = filt (num = {}, den = {}, tsam = -1, varargin)

  switch (nargin)
    case 0              # filt ()
      sys = tf ();
      return;

    case 1              # filt (sys), filt (matrix)
      if (isa (num, "lti") || is_real_matrix (num))
        sys = tf (num);
        return;
      else
        print_usage ();
      endif

    otherwise           # filt (num, den, ...)
      if (! iscell (num))
        num = {num};
      endif
      if (! iscell (den))
        den = {den};
      endif

      ## convert from z^-1 to z
      ## expand each channel by z^x, where x is the largest exponent of z^-1 (z^-x)

      ## remove trailing zeros
      ## such that polynomials are as short as possible
      num = cellfun (@__remove_trailing_zeros__, num, "uniformoutput", false);
      den = cellfun (@__remove_trailing_zeros__, den, "uniformoutput", false);

      ## make numerator and denominator polynomials equally long
      ## by adding trailing zeros
      lnum = cellfun (@length, num, "uniformoutput", false);
      lden = cellfun (@length, den, "uniformoutput", false);

      lmax = cellfun (@max, lnum, lden, "uniformoutput", false);

      num = cellfun (@postpad, num, lmax, "uniformoutput", false);
      den = cellfun (@postpad, den, lmax, "uniformoutput", false);

      ## use standard tf constructor
      ## sys is stored and displayed in standard z form, not z^-1
      sys = tf (num, den, tsam, varargin{:});
  endswitch

endfunction


function p = __remove_trailing_zeros__ (p)

  idx = find (p != 0);
  
  if (isempty (idx))
    p = 0;
  else
    p = p(1 : idx(end));
  endif

endfunction