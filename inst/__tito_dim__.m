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

## Extract nmeas and ncon from plant P which has been partitioned by mktito.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2014
## Version: 0.1

function [nmeas, ncon] = __tito_dim__ (P)

  [p, m] = size (P);
  outgroup = P.outgroup;
  ingroup = P.ingroup;
  
  if (! isfield (outgroup, "V"))
    error ("tito_dim: missing outgroup 'V'");
  endif
  
  if (! isfield (ingroup, "U"))
    error ("tito_dim: missing ingroup 'U'");
  endif
  
  nmeas = numel (outgroup.V);
  ncon = numel (ingroup.U);
  
  ## check whether indices of V and U are in ascending order
  ## and at the end of the outputs/inputs
  
  if (! isequal (outgroup.V(:), (p-nmeas+1:p)(:)))
    error ("tito_dim: outgroup 'V' invalid");
  endif
  
  if (! isequal (ingroup.U(:), (m-ncon+1:m)(:)))
    error ("tito_dim: ingroup 'U' invalid");
  endif

endfunction
