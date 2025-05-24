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
## @deftypefn {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{gain}, ...)
##
## Access transfer function data for a SISO static gain @var{GAIN}.@*
## @strong{For the full function for accessing transfer function of
## general LTI systems and with all optioinal inputs, @pxref{@@lti/tfdata}.}
##
## @strong{Inputs}
## @table @var
## @item gain
## Real valued scalar gain.
## @end table
##
## @strong{Outputs}
## @table @var
## @item num
## Numerator, (@var{gain} in this case).
## @item den
## Denominator (1 in this case).
## @item tsam
## Sampling time in seconds (0 in this case)
## @end table
##
## @seealso{@lti/tfdata}
## @end deftypefn

## Author: Torsten Lilge <ttl-octave@mailbox.org>

function [num, den, tsam] = tfdata (gain, varargin)

  if (! isscalar (gain)) || (! isreal (gain))
    error (["tfdata: has to be called with an @lti object ",...
                    "or with a real scalar (static gain)\n"]);
  endif

  ## Call @lti/tfdata after having made an transfer function from the gain
  if nargin > 1
    [num, den, tsam] = tfdata (tf(gain), varargin);
  else
    [num, den, tsam] = tfdata (tf(gain));
  endif

endfunction

