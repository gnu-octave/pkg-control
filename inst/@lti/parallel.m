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
## @deftypefn{Function File} {@var{sys} =} parallel (@var{sys1}, @var{sys2})
## Parallel connection of two @acronym{LTI} systems.
##
## @strong{Block Diagram}
## @example
## @group
##     ..........................
##     :      +--------+        :
##     :  +-->|  sys1  |---+    :
##  u  :  |   +--------+   | +  :  y
## -------+                O--------->
##     :  |   +--------+   | +  :
##     :  +-->|  sys2  |---+    :
##     :      +--------+        :
##     :.........sys............:
##
## sys = parallel (sys1, sys2)
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function sys = parallel (sys1, sys2)

  if (nargin == 2)
    sys = sys1 + sys2;

    if ((isa (sys1, "lti") && hasdelay (sys1)) || (isa (sys2, "lti") && hasdelay (sys2)))
      d1 = totaldelay (sys1);
      d2 = totaldelay (sys2);

      if (! isequal (d1, d2))
        error ("parallel: mismatched delays require internal delay support (not yet implemented)");
      endif

      sys = set (sys, "InputDelay", 0, "OutputDelay", 0, "IODelay", d1);
    endif
  ## elseif (nargin == 6)

  ## TODO: implement "complicated" case sys = parallel (sys1, sys2, in1, in2, out1, out2)

  else
    print_usage ();
  endif

endfunction


%!test  # matching delays: exact composition
%! s1 = tf (1, [1 1], "InputDelay", 0.1, "OutputDelay", 0.2);
%! s2 = tf (1, [1 2], "InputDelay", 0.1, "OutputDelay", 0.2);
%! s = parallel (s1, s2);
%! assert (s.InputDelay, 0, 1e-10);
%! assert (s.OutputDelay, 0, 1e-10);
%! assert (s.IODelay, 0.3, 1e-10);

%!test  # no delay on either operand: result has no delay (regression)
%! s1 = tf (1, [1 1]);
%! s2 = tf (1, [1 2]);
%! s = parallel (s1, s2);
%! assert (hasdelay (s), false);

%!error <mismatched delays> parallel (tf (1, [1 1], "InputDelay", 0.1, "OutputDelay", 0.2), tf (1, [1 2], "InputDelay", 0.3, "OutputDelay", 0.4))

%!test  # one operand is a plain numeric gain (no delay): should not error (regression)
%! s1 = tf (1, [1 1]);
%! s = parallel (s1, 2);
%! assert (hasdelay (s), false);
