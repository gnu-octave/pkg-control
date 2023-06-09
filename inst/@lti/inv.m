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
## @deftypefn{Function File} {@var{SYSI} = } inv (@var{SYS})
## Inversion of @acronym{LTI} objects.
##
## @strong{Inputs}
## @table @var
## @item SYS
## System to be inverted.
## @end table
##
## @strong{Outputs}
## @table @var
## @item SYSI
## Inverteted system of @var{SYS}.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function retsys = inv (sys)

  if (nargin != 1)  # prevent sys = inv (sys1, sys2, sys3, ...)
    error ("lti: inv: this is an unary operator");
  endif

  [p, m] = size (sys);

  if (p != m)
    error ("lti: inv: system must be square");
  endif

  retsys = __sys_inverse__ (sys);

  ## handle i/o names
  retsys.inname = sys.outname;
  retsys.outname = sys.inname;
  retsys.ingroup = sys.outgroup;
  retsys.outgroup = sys.ingroup; 

endfunction


## inverse of state-space models
## test from SLICOT AB07ND
## result differs intentionally from slicot
## to prevent states x_inv = -x
%!shared M, Me
%! A = [ 1.0   2.0   0.0
%!       4.0  -1.0   0.0
%!       0.0   0.0   1.0 ];
%!
%! B = [ 1.0   0.0
%!       0.0   1.0
%!       1.0   0.0 ];
%!
%! C = [ 0.0   1.0  -1.0
%!       0.0   0.0   1.0 ];
%!
%! D = [ 4.0   0.0
%!       0.0   1.0 ];
%!
%! sys = ss (A, B, C, D);
%! sysinv = inv (sys);
%! [Ai, Bi, Ci, Di] = ssdata (sysinv);
%! M = [Ai, Bi; Ci, Di];
%!
%! Ae = [ 1.0000   1.7500   0.2500
%!        4.0000  -1.0000  -1.0000
%!        0.0000  -0.2500   1.2500 ];
%!
%! Be = [-0.2500   0.0000
%!        0.0000  -1.0000
%!       -0.2500   0.0000 ];
%!
%! Ce = [ 0.0000   0.2500  -0.2500
%!        0.0000   0.0000   1.0000 ];
%!
%! De = [ 0.2500   0.0000
%!        0.0000   1.0000 ];
%!
%! Me = [Ae, -Be; -Ce, De];
%!
%!assert (M, Me, 1e-4);
