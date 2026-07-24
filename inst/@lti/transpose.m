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
## @deftypefn{Overloaded Operator} {} transpose
## Transpose of @acronym{LTI} objects.  Used by Octave for "sys.'".
## Useful for dual problems, i.e. controllability and observability
## or designing estimator gains with @command{lqr} and @command{place}.
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: February 2010
## Version: 0.2

function sys = transpose (sys)

  if (nargin != 1)  # prevent sys = transpose (sys1, sys2, sys3, ...)
    error ("lti: transpose: this is an unary operator");
  endif

  [p, m] = size (sys);

  ## Capture before __transpose__ runs: every subclass hook (tf, ss, frd;
  ## zpk round-trips through tf) leaves these lti-level fields completely
  ## untouched, still shaped/valued for the pre-transpose orientation.
  indelay = sys.indelay;
  outdelay = sys.outdelay;
  iodelay = sys.iodelay;

  sys = __transpose__ (sys);

  sys.inname = repmat ({""}, p, 1);
  sys.outname = repmat ({""}, m, 1);
  sys.ingroup = struct ();
  sys.outgroup = struct ();

  ## Swap InputDelay/OutputDelay and transpose IODelay to match the new
  ## (m x p) orientation.  IODelay(i,j) is delay from input j to output i;
  ## transposing the system swaps which index is row vs. column.
  sys.indelay = outdelay;
  sys.outdelay = indelay;
  sys.iodelay = iodelay.';

endfunction


%!test  # no delay: unaffected (regression)
%! sys = tf ({1, 1; 1, 1}, {[1 1], [1 2]; [1 3], [1 4]});
%! sys_t = sys.';
%! assert (hasdelay (sys_t), false);
%! assert (size (sys_t), [2, 2]);

%!test  # tf, non-square MIMO: InputDelay/OutputDelay swap, IODelay transposes
%! sys = tf ({1, 1; 1, 1; 1, 1}, {[1 1], [1 2]; [1 3], [1 4]; [1 5], [1 6]});
%! sys = set (sys, "InputDelay", [0.1; 0.2], "OutputDelay", [0.3; 0.4; 0.5], ...
%!                 "IODelay", [1, 2; 3, 4; 5, 6]);
%! sys_t = sys.';
%! assert (size (sys_t), [2, 3]);
%! assert (get (sys_t, "inputdelay"), [0.3; 0.4; 0.5]);
%! assert (get (sys_t, "outputdelay"), [0.1; 0.2]);
%! assert (get (sys_t, "iodelay"), [1, 3, 5; 2, 4, 6]);

%!test  # ss, non-square MIMO: same swap/transpose, different subclass hook
%! sys = ss (-1, [1, 1], [1; 1; 1], zeros (3, 2));
%! sys = set (sys, "InputDelay", [0.1; 0.2], "OutputDelay", [0.3; 0.4; 0.5], ...
%!                 "IODelay", [1, 2; 3, 4; 5, 6]);
%! sys_t = sys.';
%! assert (isa (sys_t, "ss"));
%! assert (get (sys_t, "inputdelay"), [0.3; 0.4; 0.5]);
%! assert (get (sys_t, "outputdelay"), [0.1; 0.2]);
%! assert (get (sys_t, "iodelay"), [1, 3, 5; 2, 4, 6]);

%!test  # zpk, non-square MIMO: same swap/transpose via the tf round-trip hook
%! sys = zpk ({[-1], [-1]; [-1], [-1]}, {[-1], [-2]; [-3], [-4]}, [1, 1; 1, 1]);
%! sys = set (sys, "InputDelay", [0.1; 0.2], "OutputDelay", [0.3; 0.4], ...
%!                 "IODelay", [1, 2; 3, 4]);
%! sys_t = sys.';
%! assert (isa (sys_t, "zpk"));
%! assert (get (sys_t, "inputdelay"), [0.3; 0.4]);
%! assert (get (sys_t, "outputdelay"), [0.1; 0.2]);
%! assert (get (sys_t, "iodelay"), [1, 3; 2, 4]);

%!test  # square system with nonzero IODelay: catches value-staleness
%! # distinctly from the shape-staleness case above (a square system's
%! # IODelay is the same SIZE before and after transpose, so a bug that
%! # forgot the transpose entirely would not be caught by a shape check).
%! sys = tf ({1, 1; 1, 1}, {[1 1], [1 2]; [1 3], [1 4]});
%! sys = set (sys, "IODelay", [1, 2; 3, 4]);
%! sys_t = sys.';
%! assert (get (sys_t, "iodelay"), [1, 3; 2, 4]);
