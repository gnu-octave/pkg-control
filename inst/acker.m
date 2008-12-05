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
## @deftypefn {Function File} {@var{K} =} acker (@var{a}, @var{b}, @var{p})
## A wrapper for place (@var{a}, @var{b}, @var{p}).
##
## @seealso{place}
## @end deftypefn

## Author: Luca Favatella <slackydeb@gmail.com>
## Version: 0.1

                                # TODO: modify this function to use
                                # Ackermann's formula instead of being
                                # only a wrapper for place(a,b,p)

function K = acker (a, b, p)
  if (nargin == 3)
    K = place (a, b, p);
  else
    print_usage ();
  endif
endfunction


%!test
%! A = [0 1; 3 2];
%! B = [0; 1];
%! P = [-1 -0.5];
%! Kexpected = [3.5 3.5];
%! assert (acker (A, B, P), Kexpected);