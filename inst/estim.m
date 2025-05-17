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
## @deftypefnx {Function File} {@var{est} =} estim (@var{sys}, @var{l}, @var{sensors}, @var{known}, @var{known})
## Return state estimator for a given estimator gain.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @item l
## State feedback matrix.
## @item sensors
## Indices of measured output signals y from @var{sys}. If omitted or empty,
## all outputs are measured.
## @item known
## Indices of known input signals u (deterministic) to @var{sys}.
## All other inputs to @var{sys} are assumed stochastic (w).
## If argument @var{known} is omitted or empty, no inputs u are known.
## @item type
## Type of the estimator for discrete-time systems. If set to 'delayed' the current
## estimation is based on y(k-1), if set to 'current' the current estimation is
## based on the lates mesaruement y(k). If omitted, the 'delayed' version is created.
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
## @strong{Remarks}
##
## The argument @var{type} is for discrete-time systems only. If set to 'current',
## the follwong prediction-correction scheme is used:
## @example
## @group
## ^           ^
## x*(k+1) = A x(k) + B u(k)
##    ^      ^        -1
##    x(k) = x*(k) + A   L (y(k) - C x*(k) - D u(k))
## @end group
## @end example
## The inverse fo the system matrix in the above equations is required
## for maintaining the desired observer error dynamics given by (A - LC).
##
## The advantage of this structure is that the current measurement y(k)
## is used for the current estiamted state and not for the next allowing
## the estimator to react to system disturbances faster. L is the
## observer feedback matrix for the common observer structure with
## the matrix (A - LC) being asymptotically stable, i.e. has
## eigenvalues strictly within the unit circle.
##
## @seealso{kalman, lqe, place}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.3

function est = estim (sys, l, sensors = [], known = [], type = 'delayed')

  if (nargin < 2 || nargin > 5)
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

  if strcmp (type, 'current')
    if isct (sys)
      warning ("kalman: ignoring 'type' parameter for continuous-time estimator\n");
      type = 'delayed';
    else
      if (cond (a) > 1e12)
        error ("estimd: system '%s' has noninvertibla system matrix", ...
               inputname (1));
      endif
    endif
  endif

  if strcmp (type, 'current')
    l = inv(a) * l;
    i_lc = eye(n,n) - l*c;
    f = a * i_lc;
    g = [ b-a*l*d, a*l ];
    h = [ c*i_lc
          i_lc ];
    j = [ -c*l*d, c*l;
            -l*d,   l ];
    ## k = e;
  else
    f = a - l*c;
    g = [b - l*d, l];
    h = [c; eye(n)];
    j = [d, zeros(p, p); zeros(n, m), zeros(n, p)];
    ## k = e;
  endif

  est = dss (f, g, h, j, e, tsam);
  est = set (est, "inname", inname, "stname", stname, "outname", outname);

endfunction

%!test
%!shared m, m_exp
%! sys = ss (-2, 1, 1, 3);
%! est = estim (sys, 5);
%! [a, b, c, d] = ssdata (est);
%! m = [a, b; c, d];
%! m_exp = [-7, 5; 1, 0; 1, 0];
%!assert (m, m_exp, 1e-4);

%!test
%!shared m, m_exp
%! sys = ss (-1, 2, 3, 4);
%! est = estim (sys, 5);
%! [a, b, c, d] = ssdata (est);
%! m = [a, b; c, d];
%! m_exp = [-16, 5; 3, 0; 1, 0];
%!assert (m, m_exp, 1e-4);

%!test
%! A = [ 0 1 0 ; 0 0 1 ; 0.5120 0.640 -0.800];
%! B = [ 0; 0; 1 ];
%! C = [ 0.1 0 0 ];
%! D = 1;
%! L = place (A',C',zeros(1,3))';   % Deadbeat
%! sysd = ss (A,B,C,D,1);
%! estd = estimd (sysd, L, [], 1);
%! x0  = [ .1 .1 .1 ];
%! xo0 = [ 0 0 0 ];
%! k = 0:1:25;
%! u = 0.1*sin(0.5*k).*cos(0.4*k).^2;
%! [y,t,x] = lsim (sysd, u, k, x0);
%! [yo,t,~] = lsim (estd, [u' y], k, xo0);
%! xo = yo(:,2:end);
%! e = xo - x;
%! assert (e(3,:), zeros(1,3), 1e-4);  % e already zero for k = 2

