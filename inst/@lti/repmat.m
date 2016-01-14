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
## @deftypefn {Function File} {@var{rsys} =} repmat (@var{sys}, @var{m}, @var{n})
## @deftypefnx {Function File} {@var{rsys} =} repmat (@var{sys}, [@var{m}, @var{n}])
## @deftypefnx {Function File} {@var{rsys} =} repmat (@var{sys}, @var{m})
## Form a block transfer matrix of @var{sys} with @var{m} copies vertically
## and @var{n} copies horizontally.  If @var{n} is not specified, it is set to @var{m}.
## @code{repmat (sys, 2, 3)} is equivalent to @code{[sys, sys, sys; sys, sys, sys]}.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2014
## Version: 0.1

function sys = repmat (sys, x, y)

  switch (nargin)
    case 2
      if (is_real_scalar (x))                           # repmat (sys, m)
        y = x;
      elseif (is_real_vector (x) && length (x) == 2)    # repmat (sys, [m, n])
        y = x(2);
        x = x(1);
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

  [p, m] = size (sys);
  
  out_idx = repmat (1:p, 1, x);
  in_idx = repmat (1:m, 1, y);
  
  sys = __sys_prune__ (sys, out_idx, in_idx);
  
endfunction
