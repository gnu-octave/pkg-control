## Copyright (C) 2008 Luca Favatella <slackydeb@gmail.com>
##
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} isct (@var{sys})
## Return true if the LTI system @var{sys} is continuous-time, false otherwise.
##
## @seealso{isdt}
## @end deftypefn

## Author: Luca Favatella <slackydeb@gmail.com>
## Version: 0.2.1

function retval = isct (sys)
  if (nargin != 1)
    print_usage ();
  else
    retval = (! is_digital (sys));
  endif
endfunction


%!test
%! A = [-1 0 0; 1/2 -1 0; 1/2 0 -1];
%! B = [1 0; 0 -1; 0 1];
%! C = [0 0 1; 1 1 0];
%! assert (isct (ss (A, B, C)));

%!test
%! A = [-1 0 0; 1/2 -1 0; 1/2 0 -1];
%! B = [1 0; 0 -1; 0 1];
%! C = [0 0 1; 1 1 0];
%! D = zeros (rows (C), columns (B));
%! Ts = 0.1;
%! assert (! isct (ss (A, B, C, D, Ts)));