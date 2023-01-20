## Copyright (C) 2022  Torsten Lilge
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
## @deftypefn {Function File} {@var{bool} =} isstaticgain (@var{sys})
## Determine whether @acronym{LTI} model is a static gain.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool = 0
## @var{sys} is a dynamical system
## @item bool = 1
## @var{sys} is a static gain
## @end table
## @end deftypefn

## Author: Torsten Lilge <ttl-octave@mailbox.org>
## Created: October 2022
## Version: 0.1

function static_gain = isstaticgain (ltisys)

  if (nargin == 0)
    print_usage ();
  endif

  static_gain = isempty (ltisys.a);

endfunction
