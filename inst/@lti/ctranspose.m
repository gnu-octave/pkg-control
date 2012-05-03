## Copyright (C) 2012   Lukas F. Reichlin
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
## Conjugate transpose of LTI objects.  Used by Octave for "sys'".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2012
## Version: 0.1

function sys = ctranspose (sys)

  if (nargin != 1)  # prevent sys = ctranspose (sys1, sys2, sys3, ...)
    error ("lti: ctranspose: this is an unary operator");
  endif

  [p, m] = size (sys);
  ct = isct (sys);

  sys = __ctranspose__ (sys, ct);

  sys.inname = repmat ({""}, p, 1);
  sys.outname = repmat ({""}, m, 1);

endfunction