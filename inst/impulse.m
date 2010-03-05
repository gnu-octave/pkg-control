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
## @deftypefn{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys}, @var{t})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys}, @var{tfinal})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys}, @var{tfinal}, @var{dt})
## Impulse response of LTI system.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function [y_r, t_r, x_r] = impulse (sys, tfinal = [], dt = [])

  if (nargin == 0 || nargin > 3)
    print_usage ();
  endif

  [y, t, x] = __timeresp__ (sys, "impulse", ! nargout, tfinal, dt);

  if (nargout)
    y_r = y;
    t_r = t;
    x_r = x;
  endif

endfunction
