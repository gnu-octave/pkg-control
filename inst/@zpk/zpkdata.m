## Copyright (C) 2009-2016   Lukas F. Reichlin
## Copyright (C) 2026        Mitchell Thompkins
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
## @deftypefn {Function File} {[@var{z}, @var{p}, @var{k}, @var{tsam}] =} zpkdata (@var{sys})
## @deftypefnx {Function File} {[@var{z}, @var{p}, @var{k}, @var{tsam}] =} zpkdata (@var{sys}, @var{"v"})
## Access zero-pole-gain data of a ZPK model directly, without polynomial
## conversion.
## @end deftypefn

function [z, p, k, tsam] = zpkdata (sys, rtype = "cell")

  [z, p, k] = __sys_data__ (sys);
  tsam = get (sys, "tsam");

  if (strncmpi (rtype, "v", 1) && issiso (sys))
    z = z{1};
    p = p{1};
  endif

endfunction


%!test
%! ze = {[1]; [-2; 0]};
%! pe = {[-1; 0]; [-4; -3; -1]};
%! ke = [5; 10];
%! sys = zpk (ze, pe, ke);
%! [zo, po, ko] = zpkdata (sys);
%! assert (zo, ze);
%! assert (po, pe);
%! assert (ko, ke);

%!test
%! sys = zpk ([], [-1; -2], 3);
%! [z, p, k, tsam] = zpkdata (sys, 'v');
%! assert (isempty (z));
%! assert (sort (p), [-2; -1]);
%! assert (k, 3);
%! assert (tsam, 0);
