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
## Compute frequency response of ZPK model by delegating to TF.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function H = __freqresp__ (sys, w, cellflag = false)
  H = __freqresp__ (tf (sys), w, cellflag);
endfunction
