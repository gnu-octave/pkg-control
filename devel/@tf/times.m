## Copyright (C) 2009-2014   Lukas F. Reichlin
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
## Hadamard/Schur product of @acronym{TF} objects.
## Used by Octave for "sys1 .* sys2".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2014
## Version: 0.1

function sys = times (sys1, sys2)

  if (nargin != 2)                          # prevent sys = times (sys1, sys2, sys3, ...)
    error ("tf: times: this is a binary operator");
  endif

  if (! isa (sys1, "tf"))
    sys1 = tf (sys1);
  endif

  if (! isa (sys2, "tf"))
    sys2 = tf (sys2);
  endif

  [p1, m1] = size (sys1);
  [p2, m2] = size (sys2);
  
  if (p1 != p2 || m1 != m2)
    error ("tf: times: system dimensions incompatible: (%dx%d) .* (%dx%d)", ...
            p1, m1, p2, m2);
  endif
  
  num = cellfun (@mtimes, sys1.num, sys2.num, "uniformoutput", false)
  den = cellfun (@mtimes, sys1.den, sys2.den, "uniformoutput", false)

  ## FIXME: catch case  ss .* tf  or  tf .* ss, should error out
  ##        maybe this needs a dummy times.m in @ss and @frd

endfunction
