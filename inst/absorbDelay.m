## Copyright (C) 2026        Prateek Ganguli
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
## @deftypefn {Function File} {@var{sys} =} absorbDelay (@var{sys})
## Absorb the delays of a discrete-time @var{tf} or @var{zpk} model into
## the denominator or poles by adding integrators at the origin.
##
## For continuous-time delays or @var{ss} models, an error is raised.
##
## @strong{Inputs}
## @table @var
## @item sys
## Discrete-time @var{tf} or @var{zpk} system with delays.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Delay-free equivalent system (exact for discrete, integrators added).
## @end table
## @end deftypefn

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function sys = absorbDelay (sys)

  if (nargin != 1 || ! isa (sys, "lti"))
    print_usage ();
  endif

  if (! (isa (sys, "tf") || isa (sys, "zpk")))
    error ("absorbDelay: only tf and zpk models are supported (ss not yet implemented)");
  endif

  if (! hasdelay (sys))
    return;
  endif

  if (isct (sys))
    error ("absorbDelay: continuous-time delays cannot be absorbed exactly; use pade() to approximate them first");
  endif

  total = totaldelay (sys);
  origsys = sys;

  if (isa (sys, "zpk"))
    [z, p, k, tsam] = zpkdata (sys);
    [pr, pc] = size (p);

    for i = 1 : pr
      for j = 1 : pc
        n = total (i, j);
        if (n > 0)
          p{i,j} = [p{i,j}; zeros(n, 1)];
        endif
      endfor
    endfor

    sys = zpk (z, p, k, tsam);
  else
    [num, den, tsam] = tfdata (sys);
    [pr, pc] = size (den);

    for i = 1 : pr
      for j = 1 : pc
        n = total (i, j);
        if (n > 0)
          den{i,j} = [den{i,j}, zeros(1, n)];
        endif
      endfor
    endfor

    sys = tf (num, den, tsam);
  endif

  sys = set (sys, "lti", origsys);
  sys = set (sys, "InputDelay", 0, "OutputDelay", 0, "IODelay", 0);

endfunction


%!test  # discrete SISO tf with InputDelay: exact equivalence via freqresp
%! # freqresp() does not yet apply a system's stored delay (that is a
%! # later, separate phase), so the expected response is built by hand:
%! # delay-free rational response times the z^-k phase factor.
%! sys = tf (1, [1 -0.5], 0.1, "InputDelay", 2);
%! sys2 = absorbDelay (sys);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! expected = freqresp (tf (1, [1 -0.5], 0.1), w) .* reshape (exp (-1i * w * 0.1 * 2), 1, 1, []);
%! assert (freqresp (sys2, w), expected, 1e-10);

%!test  # discrete MIMO tf with differing per-channel IODelay
%! sys = tf ({1, 1; 1, 1}, {[1 -0.5], [1 -0.6]; [1 -0.7], [1 -0.8]}, 0.1);
%! sys = set (sys, "IODelay", [1, 0; 0, 2]);
%! sys2 = absorbDelay (sys);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! rational = tf ({1, 1; 1, 1}, {[1 -0.5], [1 -0.6]; [1 -0.7], [1 -0.8]}, 0.1);
%! total = [1, 0; 0, 2];
%! resp = freqresp (rational, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* exp (-1i * w * 0.1 * total(i,j));
%!   endfor
%! endfor
%! assert (freqresp (sys2, w), expected, 1e-8);

%!test  # zpk discrete with OutputDelay
%! sys = zpk ([], -0.5, 2, 0.1, "OutputDelay", 3);
%! sys2 = absorbDelay (sys);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! expected = freqresp (zpk ([], -0.5, 2, 0.1), w) .* reshape (exp (-1i * w * 0.1 * 3), 1, 1, []);
%! assert (freqresp (sys2, w), expected, 1e-10);

%!test  # no delay: no-op
%! sys = tf (1, [1 -0.5], 0.1);
%! sys2 = absorbDelay (sys);
%! [num1, den1] = tfdata (sys);
%! [num2, den2] = tfdata (sys2);
%! assert (num2, num1);
%! assert (den2, den1);
%! assert (hasdelay (sys2), false);

%!test  # continuous-time, no delay: no-op, no error
%! sys = tf (1, [1 1]);
%! sys2 = absorbDelay (sys);
%! assert (hasdelay (sys2), false);

%!error <continuous-time delays> absorbDelay (tf (1, [1 1], "InputDelay", 0.5))

%!error <ss not yet implemented> absorbDelay (ss (-1, 1, 1, 0, 0.1, "InputDelay", 1))
