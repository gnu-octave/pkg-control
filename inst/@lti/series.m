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
## @deftypefn {Function File} {@var{sys} =} series (@var{sys1}, @var{sys2})
## @deftypefnx {Function File} {@var{sys} =} series (@var{sys1}, @var{sys2}, @var{outputs1}, @var{inputs2})
## Series connection of two @acronym{LTI} models.
##
## @strong{Block Diagram}
## @example
## @group
##     .....................................
##  u  :  +--------+ y1    u2  +--------+  :  y
## ------>|  sys1  |---------->|  sys2  |------->
##     :  +--------+           +--------+  :
##     :................sys.................
##
## sys = series (sys1, sys2)
## @end group
## @end example
## @example
## @group
##     .....................................
##     :                   v2  +--------+  :
##     :            ---------->|        |  :  y
##     :  +--------+ y1    u2  |  sys2  |------->
##  u  :  |        |---------->|        |  :
## ------>|  sys1  |       z1  +--------+  :
##     :  |        |---------->            :
##     :  +--------+                       :
##     :................sys.................
##
## outputs1 = [1]
## inputs2 = [2]
## sys = series (sys1, sys2, outputs1, inputs2)
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function sys = series (sys1, sys2, out1, in2)

  if (nargin == 2)
    sys = sys2 * sys1;

    [p1, m1] = size (sys1);
    [p2, m2] = size (sys2);

    if ((isa (sys1, "lti") && hasdelay (sys1)) || (isa (sys2, "lti") && hasdelay (sys2)))
      if (p1 != 1 || m2 != 1)
        error ("series: MIMO connections with delays are not yet supported (require internal delay support)");
      endif

      id1 = get (sys1, "inputdelay");
      od2 = get (sys2, "outputdelay");
      io_junction = get (sys1, "outputdelay") + get (sys2, "inputdelay") ...
                    + get (sys1, "iodelay") + get (sys2, "iodelay");

      sys = set (sys, "InputDelay", id1, "OutputDelay", od2, "IODelay", io_junction);
    endif
  elseif (nargin == 4)
    [p1, m1] = size (sys1);
    [p2, m2] = size (sys2);

    if ((isa (sys1, "lti") && hasdelay (sys1)) || (isa (sys2, "lti") && hasdelay (sys2)))
      error ("series: MIMO connections with delays are not yet supported (require internal delay support)");
    endif

    if (! is_real_vector (out1))
      error ("series: argument 3 (outputs1) invalid");
    endif

    if (! is_real_vector (in2))
      error ("series: argument 4 (inputs2) invalid");
    endif

    l_out1 = length (out1);
    l_in2 = length (in2);

    if (l_out1 > p1)
      error ("series: 'outputs1' has too many indices for 'sys1'");
    endif

    if (l_in2 > m2)
      error ("series: 'inputs2' has too many indices for 'sys2'");
    endif

    if (l_out1 != l_in2)
      error ("series: number of 'outputs1' and 'inputs2' indices must be equal");
    endif

    if (any (out1 > p1 | out1 < 1))
      error ("series: range of 'outputs1' indices exceeds dimensions of 'sys1'");
    endif

    if (any (in2 > m2 | in2 < 1))
      error ("series: range of 'inputs2' indices exceeds dimensions of 'sys2'");
    endif

    out_scl = full (sparse (1:l_out1, out1, 1, l_out1, p1));
    in_scl = full (sparse (in2, 1:l_in2, 1, m2, l_in2));
    
    ## NOTE: for-loop does NOT the same as
    ##       out_scl(1:l_out1, out1) = 1;
    ##       in_scl(in2, 1:l_out1) = 1;
    ##
    ## out_scl = zeros (l_out1, p1);
    ## in_scl = zeros (m2, l_in2);
    ##
    ## for k = 1 : l_out1
    ##   out_scl(k, out1(k)) = 1;
    ##   in_scl(in2(k), k) = 1;
    ## endfor

    scl = in_scl * out_scl;
    sys = sys2 * scl * sys1;
  else
    print_usage ();
  endif

endfunction


%!test  # single-channel cascade: exact delay composition
%! s1 = tf (1, [1 1], "InputDelay", 0.1, "OutputDelay", 0.2);
%! s2 = tf (1, [1 2], "InputDelay", 0.3, "OutputDelay", 0.4);
%! s = series (s1, s2);
%! assert (s.InputDelay, 0.1, 1e-10);
%! assert (s.OutputDelay, 0.4, 1e-10);
%! assert (s.IODelay, 0.5, 1e-10);

%!test  # no delay on either operand: result has no delay (regression)
%! s1 = tf (1, [1 1]);
%! s2 = tf (1, [1 2]);
%! s = series (s1, s2);
%! assert (hasdelay (s), false);

%!error <MIMO> series (tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}, "InputDelay", [0.1;0.2]), tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}))

%!test  # one operand is a plain numeric gain (no delay): should not error (regression)
%! s1 = tf (1, [1 1]);
%! s = series (s1, 2);
%! assert (hasdelay (s), false);

%!error <MIMO> series (tf (1, [1 1], "InputDelay", 0.1), tf (1, [1 2]), [1], [1])
