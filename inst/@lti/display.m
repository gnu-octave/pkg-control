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
## Display routine for LTI objects.  Called by its child classes.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function display (sys)

  if (nfields2 (sys.ingroup) > 0)
    __disp_group__ (sys.ingroup, "Input");
  endif
  
  if (nfields2 (sys.outgroup) > 0)
    __disp_group__ (sys.outgroup, "Output");
  endif

  if (! isempty (sys.name))
    disp (["Name: ", sys.name]);
  endif

  if (sys.tsam > 0)
    disp (sprintf ("Sampling time: %g s", sys.tsam));
  elseif (sys.tsam == -1)
    disp ("Sampling time: unspecified");
  endif

endfunction


function __disp_group__ (group, io)

  name = fieldnames (group);
  idx = struct2cell (group);

  cellfun (@(name, idx) printf ("%s group '%s' = %s\n", io, name, mat2str (idx(:).')), ...
                                name, idx, "uniformoutput", false);

endfunction
