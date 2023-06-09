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

  sys = __transpose__ (sys);

  sys.inname = repmat ({""}, p, 1);
  sys.outname = repmat ({""}, m, 1);
  sys.ingroup = struct ();
  sys.outgroup = struct ();

endfunction
