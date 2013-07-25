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
## @deftypefn {Function File} {@var{p} =} pole (@var{sys})
## Compute poles of @acronym{LTI} system.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item p
## Poles of @var{sys}.
## @end table
##
## @strong{Algorithm}@*
## For (descriptor) state-space models, @command{pole}
## relies on Octave's @command{eig}.
## For @acronym{SISO} transfer functions, @command{pole}
## uses Octave's @command{roots}.
## @acronym{MIMO} transfer functions are converted to
## a minimal state-space representation for the
## computation of the poles.
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function pol = pole (sys)

  if (nargin > 1)
    print_usage ();
  endif

  pol = __pole__ (sys);

endfunction
