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
## @deftypefn {Function File} {@var{est} =} estim (@var{sys}, @var{l})
## @deftypefnx {Function File} {@var{est} =} estim (@var{sys}, @var{l}, @var{sensors}, @var{known})
## Return state estimator for a given estimator gain.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @item l
## State feedback matrix.
## @item sensors
## Indices of measured output signals y from @var{sys}.  If omitted, all outputs are measured.
## @item known
## Indices of known input signals u (deterministic) to @var{sys}.  All other inputs to @var{sys}
## are assumed stochastic (w).  If argument @var{known} is omitted, no inputs u are known.
## @end table
##
## @strong{Outputs}
## @table @var
## @item est
## State-space model of estimator.
## @end table
##
## @strong{Block Diagram}
## @example
## @group
##                                  u  +-------+         ^
##       +---------------------------->|       |-------> y
##       |    +-------+     +       y  |  est  |         ^
## u ----+--->|       |----->(+)------>|       |-------> x
##            |  sys  |       ^ +      +-------+
## w -------->|       |       |
##            +-------+       | v
## @end group
## @end example
##
## @seealso{kalman, lqe, place}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.3

function est = estim (sys, l, sensors = [], known = [])

  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("estim: first argument must be an LTI system");
  endif

  sys = ss (sys);     # needed to get stname from tf models
  [a, b, c, d, e, tsam] = dssdata (sys, []);
  [inn, stn, outn, ing, outg] = get (sys, "inname", "stname", "outname", "ingroup", "outgroup");

  if (isempty (sensors))
    sensors = 1 : rows (c);
  endif
  
  if (ischar (sensors))
    sensors = {sensors};
  endif
  if (ischar (known))
    known = {known};
  endif

  if (iscell (sensors))
    tmp = cellfun (@(x) __str2idx__ (outg, outn, x, "out"), sensors(:), "uniformoutput", false);
    sensors = vertcat (tmp{:});
  endif
  if (iscell (known))
    tmp = cellfun (@(x) __str2idx__ (ing, inn, x, "in"), known(:), "uniformoutput", false);
    known = vertcat (tmp{:});
  endif

  m = length (known);
  n = rows (a);
  p = length (sensors);

  if (rows (l) != n)
    error ("estim: system '%s' has %d states, but the state estimator gain '%s' has %d rows", ...
            inputname (1), n, inputname (2), rows (l));
  endif

  if (columns (l) != p)
    error ("estim: estimator gain '%s' has %d columns, but argument 'known' contains %d indices", ...
            inputname (2), columns (l), p);
  endif

  b = b(:, known);
  c = c(sensors, :);
  d = d(sensors, known);

  stname = __labels__ (stn, "xhat");
  outname = vertcat (__labels__ (outn(sensors(:)), "yhat"), stname);
  inname = vertcat (__labels__ (inn(known(:)), "u"), __labels__ (outn(sensors(:)), "y"));

  f = a - l*c;
  g = [b - l*d, l];
  h = [c; eye(n)];
  j = [d, zeros(p, p); zeros(n, m), zeros(n, p)];
  ## k = e;

  est = dss (f, g, h, j, e, tsam);
  est = set (est, "inname", inname, "stname", stname, "outname", outname);

endfunction


%!shared m, m_exp
%! sys = ss (-2, 1, 1, 3);
%! est = estim (sys, 5);
%! [a, b, c, d] = ssdata (est);
%! m = [a, b; c, d];
%! m_exp = [-7, 5; 1, 0; 1, 0];
%!assert (m, m_exp, 1e-4);

%!shared m, m_exp
%! sys = ss (-1, 2, 3, 4);
%! est = estim (sys, 5);
%! [a, b, c, d] = ssdata (est);
%! m = [a, b; c, d];
%! m_exp = [-16, 5; 3, 0; 1, 0];
%!assert (m, m_exp, 1e-4);
