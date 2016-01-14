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
## Vertical concatenation of @acronym{TF} objects.
## Used by Octave for "[sys1; sys2]".
## Avoids conversion to state-space and back by overriding
## the general vertcat function for @acronym{LTI} objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2014
## Version: 0.1

function sys = vertcat (sys, varargin)

  sys = tf (sys);
  varargin = cellfun (@tf, varargin, "uniformoutput", false);

  for k = 1 : (nargin-1)
  
    sys1 = sys;
    sys2 = varargin{k};
    
    sys = tf ();
    sys.lti = __lti_group__ (sys1.lti, sys2.lti, "vertcat");
    
    [p1, m1] = size (sys1.num);
    [p2, m2] = size (sys2.num);
    
    if (m1 != m2)
      error ("tf: vertcat: number of system inputs incompatible: [(%dx%d); (%dx%d)]",
              p1, m1, p2, m2);
    endif
    
    sys.num = [sys1.num; sys2.num];
    sys.den = [sys1.den; sys2.den];
    
    if (strcmp (sys1.tfvar, sys2.tfvar))
      sys.tfvar = sys1.tfvar;
    elseif (strcmp (sys1.tfvar, "x"))
      sys.tfvar = sys2.tfvar;
    else
      sys.tfvar = sys1.tfvar;
    endif

    if (sys1.inv || sys2.inv)
      sys.inv = true;
    endif

  endfor

endfunction
