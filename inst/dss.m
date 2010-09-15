## Copyright (C) 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} dss (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} dss (@var{d})
## @deftypefnx {Function File} {@var{sys} =} dss (@var{a}, @var{b}, @var{c}, @var{d}, @var{e})
## @deftypefnx {Function File} {@var{sys} =} dss (@var{a}, @var{b}, @var{c}, @var{d}, @var{e}, @var{tsam})
## Create or convert to state-space model.
##
## @strong{Inputs}
## @table @var
## @item a
## State transition matrix.
## @item b
## Input matrix.
## @item c
## Measurement matrix.
## @item d
## Feedthrough matrix.
## @item e
## Descriptor matrix.
## @item tsam
## Sampling time. If @var{tsam} is not specified, a continuous-time
## model is assumed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Descriptor state-space model.
## @end table
##
## @seealso{ss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2010
## Version: 0.1

function sys = dss (varargin)

  switch (nargin)
    case {0, 2, 3, 4}
      print_usage ();

    case 1  # static gain
      sys = ss (varargin{1});

    otherwise  # general case
      sys = ss (varargin{[1:4, 6:end]}, "e", varargin{5});
  endswitch

endfunction