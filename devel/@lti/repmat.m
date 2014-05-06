## Copyright (C) 2009-2014   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{rsys} =} repmat (@var{sys}, @var{m}, @var{n})
## @deftypefnx {Function File} {@var{rsys} =} repmat (@var{sys}, [@var{m}, @var{n}])
## @deftypefnx {Function File} {@var{rsys} =} repmat (@var{sys}, @var{m})
## TODO
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2014
## Version: 0.1

function sys = repmat (sys, x = 1, y = x)

  if (! isa (sys, "lti"))                               # nargin always >= 1
    error ("lti: repmat: first argument must be an LTI system");
  endif

  switch (nargin)
    case 2
      if (is_real_scalar (x))                           # repmat (sys, m)
        ## nothing to do here, y = x
      elseif (is_real_vector (x) && length (x) == 2)    # repmat (sys, [m, n])
        x = x(1);
        y = y(2);
      else
        error ("lti: repmat: second argument invalid");
      endif
    case 3                                              # repmat (sys, m, n)
      if (! is_real_scalar (x, y))
        error ("lti: repmat: dimensions 'm' and 'n' must be real integers");
      endif
    otherwise
      print_usage ();
  endswitch

  tmp = cell (x, 1);
  tmp(1:x) = sys;
  sys = horzcat (tmp{:});
  
  tmp = cell (y, 1);
  tmp(1:y) = sys;
  sys = vertcat (tmp{:});

endfunction
