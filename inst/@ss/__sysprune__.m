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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Submodel extraction and reordering for SS objects.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function sys = __sysprune__ (sys, out_idx, in_idx, st_idx = ":")

  sys.lti = __ltiprune__ (sys.lti, out_idx, in_idx);

  sys.a = sys.a(st_idx, st_idx);
  sys.b = sys.b(st_idx, in_idx);
  sys.c = sys.c(out_idx, st_idx);
  sys.d = sys.d(out_idx, in_idx);

  sys.stname = sys.stname(st_idx);

endfunction