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
## End indexing for @acronym{LTI} objects.
## Used by Octave for "sys(1:end, end-1)".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: December 2013
## Version: 0.1

function ret = end (sys, k, n)

  if (n != 2)
    error ("lti: end: require 2 indices in the expression");
  endif

  [p, m] = size (sys);

  switch (k)
    case 1
      ret = p;
    case 2
      ret = m;
    otherwise
      error ("lti: end: invalid expression index k = %d", k);
  endswitch

endfunction
