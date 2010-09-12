## Copyright (C) 1996, 2000, 2002, 2003, 2004, 2005, 2007
##               Auburn University.  All rights reserved.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} issample (@var{ts})
## @deftypefnx {Function File} {} issample (@var{ts}, @var{flg})
## Return true if @var{ts} is a valid sampling time
## (real, scalar, > 0). If a second argument != 0
## is passed, -1 and 0 become valid sample times as well.
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: July 1995

## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: September 2009
## Version: 0.2

function bool = issample (tsam, flg = 0)

  if (nargin < 1 || nargin > 2)
    print_usage (); 
  endif

  if (flg == 0)  # refuse -1 and 0
    bool = (isreal (tsam) && isscalar (tsam) && (tsam > 0));       
  else  # allow -1 and 0
    bool = (isreal (tsam) && isscalar (tsam) && (tsam >= 0 || tsam == -1));
  endif

endfunction


## flg == 0
%!assert (issample (1), true)
%!assert (issample (pi), true)
%!assert (issample (0), false)
%!assert (issample (-1), false)
%!assert (issample (-1, 0), false)
%!assert (issample ("a"), false)
%!assert (issample (eye (2)), false)
%!assert (issample (2+2i), false)

## flg != 0
%!assert (issample (-1, 1), true)
%!assert (issample (-1, -1), true)
%!assert (issample (pi, 1), true)
%!assert (issample (0, 1), true)
%!assert (issample ("b", 1), false)
%!assert (issample (-1, "ab"), true)
%!assert (issample (rand (3,2), 1), false)
%!assert (issample (2+2i, 1), false)
