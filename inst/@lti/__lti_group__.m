## Copyright (C) 2009   Lukas F. Reichlin
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
## Block diagonal concatenation of two LTI models.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function retlti = __lti_group__ (lti1, lti2)

  retlti = lti ();

  retlti.inname = [lti1.inname;
                   lti2.inname];

  retlti.outname = [lti1.outname;
                    lti2.outname];

  if (nfields (lti1.ingroup) || nfields (lti2.ingroup))
    m1 = numel (lti1.inname);
    lti2.ingroup = structfun (@(x) x + m1, lti2.ingroup, "uniformoutput", false);
    retlti.ingroup = __merge_struct__ (lti1.ingroup, lti2.ingroup, "in");
  endif
  
  if (nfields (lti1.outgroup) || nfields (lti2.outgroup))
    p1 = numel (lti1.outname);
    lti2.outgroup = structfun (@(x) x + p1, lti2.outgroup, "uniformoutput", false);
    retlti.outgroup = __merge_struct__ (lti1.outgroup, lti2.outgroup, "out");
  endif

  if (lti1.tsam == lti2.tsam)
    retlti.tsam = lti1.tsam;
  elseif (lti1.tsam == -2)
    retlti.tsam = lti2.tsam;
  elseif (lti2.tsam == -2)
    retlti.tsam = lti1.tsam;
  elseif (lti1.tsam == -1 && lti2.tsam > 0)
    retlti.tsam = lti2.tsam;
  elseif (lti2.tsam == -1 && lti1.tsam > 0)
    retlti.tsam = lti1.tsam;
  else
    error ("lti_group: systems must have identical sampling times");
  endif

endfunction


function ret = __merge_struct__ (a, b, iostr)

  ## FIXME: this is too complicated;
  ##        isn't there a simple function for this task?

  ## bug #40224: orderfields (struct ()) errors out in Octave 3.6.4
  ## therefore use nfields to check for empty structs
  if (nfields (a))
    a = orderfields (a);
  endif
  if (nfields (b))
    b = orderfields (b);
  endif

  fa = fieldnames (a);
  fb = fieldnames (b);
  [fi, ia, ib] = intersect (fa, fb);
  ca = struct2cell (a);
  cb = struct2cell (b);

  for k = numel (fi) : -1 : 1
    ca{ia(k)} = vertcat (ca{ia(k)}(:), cb{ib(k)}(:));
    fb(ib(k)) = [];
    cb(ib(k)) = [];
  endfor
  
  ret = cell2struct ([ca; cb], [fa; fb]);

endfunction
