## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Select a subset of inputs/outputs from a ZPK model (lossless indexing).

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function sys = __sys_prune__ (sys, out_idx, in_idx)

  [sys.lti, out_idx, in_idx] = __lti_prune__ (sys.lti, out_idx, in_idx);

  sys.z = sys.z(out_idx, in_idx);
  sys.p = sys.p(out_idx, in_idx);
  sys.k = sys.k(out_idx, in_idx);

endfunction
