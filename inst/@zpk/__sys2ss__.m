## Copyright (C) 2026  Mitchell Thompkins
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
## Convert ZPK model to state-space representation.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function [retsys, retlti] = __sys2ss__ (sys)

  ## explicit representation change. goes via the (lossy) polynomial form,
  ## same numerics as the former zpk-returns-tf wrapper
  [retsys, retlti] = __sys2ss__ (tf (sys));

endfunction
