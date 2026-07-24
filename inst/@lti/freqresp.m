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
## @deftypefn{Function File} {@var{H} =} freqresp (@var{sys}, @var{w})
## Evaluate frequency response at given frequencies.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.
## @item w
## Vector of frequency values.
## @end table
##
## @strong{Outputs}
## @table @var
## @item H
## Array of frequency response.  For a system with m inputs and p outputs, the array @var{H}
## has dimensions [p, m, length (w)].
## The frequency response at the frequency w(k) is given by H(:,:,k).
## @end table
##
## @seealso{@@lti/dcgain}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function H = freqresp (sys, w)

  if (nargin != 2)           # case freqresp () not possible
    print_usage ();
  endif

  if (! is_real_vector (w))  # catches freqresp (sys, sys) and freqresp (w, sys) as well
    error ("freqresp: second argument 'w' must be a real-valued vector of frequencies");
  endif

  H = __freqresp__ (sys, w);
  H = __apply_freqresp_delay__ (sys, H, w, false);

endfunction

%!test  # no-delay system: unaffected (regression)
%! sys = tf (1, [1, 1]);
%! w = [0.1, 1, 5];
%! H = freqresp (sys, w);
%! rational_H = __freqresp__ (sys, w);
%! assert (H, rational_H);

%!test  # SISO continuous tf with InputDelay: freqresp applies the delay
%! sys = tf (1, [1, 1], "InputDelay", 0.3);
%! rational = tf (1, [1, 1]);
%! w = [0.1, 1, 5];
%! H = freqresp (sys, w);
%! expected = freqresp (rational, w) .* reshape (exp (-1i * w * 0.3), 1, 1, []);
%! assert (H, expected, 1e-10);

%!test  # SISO discrete tf with integer InputDelay: freqresp applies the delay
%! sys = tf (1, [1, -0.5], 0.1, "InputDelay", 2);
%! rational = tf (1, [1, -0.5], 0.1);
%! w = [0.1, 1, 5];
%! H = freqresp (sys, w);
%! expected = freqresp (rational, w) .* reshape (exp (-1i * w * 0.1 * 2), 1, 1, []);
%! assert (H, expected, 1e-10);

%!test  # MIMO tf with per-channel IODelay: freqresp applies each channel's own delay
%! sys = tf ({1, 1; 1, 1}, {[1, -0.5], [1, -0.6]; [1, -0.7], [1, -0.8]}, 0.1);
%! sys = set (sys, "IODelay", [1, 0; 0, 2]);
%! rational = tf ({1, 1; 1, 1}, {[1, -0.5], [1, -0.6]; [1, -0.7], [1, -0.8]}, 0.1);
%! total = [1, 0; 0, 2];
%! w = [0.1, 1, 5];
%! resp = freqresp (rational, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* exp (-1i * w * 0.1 * total(i,j));
%!   endfor
%! endfor
%! assert (freqresp (sys, w), expected, 1e-8);
