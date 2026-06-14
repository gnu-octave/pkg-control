## Copyright (C) 2026 Dmitri A. Sergatskov <dasergatskov@gmail.com>
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
## @deftypefn {Function File} {@var{K} =} dcgain (@var{sys})
##
## Compute the DC gain of @acronym{LTI} systems given as transfer function.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system created by tf(), ss(), dss(), etc.
## @end table
##
## @strong{Outputs}
## @table @var
## @item K
## DC gain matrix. For a system with m inputs and p outputs, the array @var{k}
## has dimensions [p, m].
## @end table
##
## @xref{@@lti/dcgain} for the method used for systems not represented by a transfer function.
##
## For a continuous-time system G(s) of the form
## @tex
## $$ G(s) = \frac{b_ns^n + \cdots + b_1s+ b_0}{a_ns^n + \cdots + a_1s+ a_0} $$
## @end tex
## @ifnottex
## @example
## @group
##           bn*s^n + ... + b1*s + b0
##   G(s) = --------------------------
##           an*s^n + ... + a1*s + a0
## @end group
## @end example
## @end ifnottex
## the gain is given by
## @tex
## $$ K = G(s=0) = \frac{b_0}{a_0} $$
## @end tex
## @ifnottex
## @example
## @group
##                 b0
##   K = G(s=0) = ----
##                 a0
## @end group
## @end example
## @end ifnottex
##
## In the discrete-time case the gain of
## @tex
## $$ G(z) = \frac{b_nz^n + \cdots + b_1z+ b_0}{a_nz^n + \cdots + a_1z+ a_0} $$
## @end tex
## @ifnottex
## @example
## @group
##           bn*z^n + ... + b1*z + b0
##   G(z) = --------------------------
##           an*z^n + ... + a1*z + a0
## @end group
## @end example
## @end ifnottex
## is given by
## @tex
## $$ K = G(z=1) = \frac{\sum_{i=0}^n b_i}{\sum_{i=0}^n a_i} $$
## @end tex
## @ifnottex
## @example
## @group
##                 bn + ... + b1 + b0
##   K = G(z=1) = --------------------
##                 an + ... + a1 + a0
## @end group
## @end example
## @end ifnottex
##
## @strong{Example}
## @example
## @group
## G = tf (@{[1 0],[1 1]@},@{[1 1],[1,0]@});
## K = dcgain (G)
## K =
##   0  Inf
## @end group
## @end example
##
## @seealso{@@lti/dcgain,@@lti/freqresp,@@tf/tf,@@ss/ss,dss}
## @end deftypefn

function gain = dcgain (sys)

  if (nargin != 1)
    print_usage ();
  endif

  [num, den] = tfdata (sys);

  if (isct (sys))
    gain = cellfun (@(n, d) n(end) ./ d(end), num, den);
  else
    gain = cellfun (@(n, d) polyval (n, 1) ./ polyval (d, 1), num, den);
  endif

endfunction

%!assert (dcgain (tf (1, [1, 1])), 1)
%!assert (dcgain (tf (2, [1, 1])), 2)
%!assert (dcgain (tf (1, [1, -0.5], 1)), 2)
%!assert (dcgain (tf (1, [1, -1], 1)), Inf)

%!test
%! num = [0, 0, 2.052877715426585e9, 4.416784109977014e13, 1.719831064785571e17];
%! den = [1, 8.380906856780281e4, 2.269817212624148e8, 4.344755797517798e12, 0];
%! gain = dcgain (tf (num, den));
%! assert (isinf (gain));
%! assert (gain > 0);
