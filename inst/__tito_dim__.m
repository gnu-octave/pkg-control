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

## TODO

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2014
## Version: 0.1

function [nmeas, ncon] = __tito_dim__ (P)

  outgroup = P.outgroup;
  ingroup = P.ingroup;
  
  ## FIXME: add checks whether indices of V and U are in ascending order
  ##        and at the end of the outputs/inputs
  
  if (isfield (outgroup, "V"))
    nmeas = numel (outgroup.V);
  else
    error ("tito_dim: missing outgroup 'V'");
  endif
  
  if (isfield (ingroup, "U"))
    ncon = numel (ingroup.U);
  else
    error ("tito_dim: missing ingroup 'U'");
  endif

endfunction
