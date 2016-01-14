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
## Hadamard/Schur product of @acronym{SS} objects.
## Used by Octave for "sys1 .* sys2".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2014
## Version: 0.1

function sys = __times__ (sys1, sys2)

  error ("ss: times: Hadamard/Schur product is only supported for transfer functions");

endfunction
