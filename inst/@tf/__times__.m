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
## Hadamard/Schur product of @acronym{TF} objects.
## Used by Octave for "sys1 .* sys2".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2014
## Version: 0.1

function sys = __times__ (sys1, sys2)

  if (! isa (sys1, "tf"))
    sys1 = tf (sys1);
  endif

  if (! isa (sys2, "tf"))
    sys2 = tf (sys2);
  endif

  sys = tf ();
  sys.lti = __lti_group__ (sys1.lti, sys2.lti, "times");
  
  sys.num = cellfun (@mtimes, sys1.num, sys2.num, "uniformoutput", false);
  sys.den = cellfun (@mtimes, sys1.den, sys2.den, "uniformoutput", false);

  if (sys1.tfvar == sys2.tfvar)
    sys.tfvar = sys1.tfvar;
  elseif (sys1.tfvar == "x")
    sys.tfvar = sys2.tfvar;
  else
    sys.tfvar = sys1.tfvar;
  endif

  if (sys1.inv || sys2.inv)
    sys.inv = true;
  endif

endfunction
