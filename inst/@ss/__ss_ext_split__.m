## Copyright (C) 2026        Prateek Ganguli
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
## Split a discretized/continuized extended system (built by
## __ss_ext_build__) back into the ordinary dynamics (a/b/c/d/e) plus the
## delay-port fields (b2/c2/d12/d21/d22), reattaching them to a copy of the
## original InternalDelay-carrying ss model.
##
## @var{nu} is the number of columns of the original (non-extended) B; the
## first @var{nu} columns of @code{ext_sys.b} are the new B, the remainder
## the new B2.  @var{ny} is the number of rows of the original C; analogous
## row split for C/C2.  D is split into its four quadrants accordingly.
##
## For internal use only, called from @lti/c2d.m and @lti/d2c.m.

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function newsys = __ss_ext_split__ (origsys, ext_sys, nu, ny)

  newsys = origsys;
  newsys.a = ext_sys.a;
  newsys.e = ext_sys.e;
  newsys.b = ext_sys.b(:, 1:nu);
  newsys.b2 = ext_sys.b(:, nu+1:end);
  newsys.c = ext_sys.c(1:ny, :);
  newsys.c2 = ext_sys.c(ny+1:end, :);
  newsys.d = ext_sys.d(1:ny, 1:nu);
  newsys.d12 = ext_sys.d(1:ny, nu+1:end);
  newsys.d21 = ext_sys.d(ny+1:end, 1:nu);
  newsys.d22 = ext_sys.d(ny+1:end, nu+1:end);

endfunction
