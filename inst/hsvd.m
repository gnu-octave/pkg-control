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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {@var{hsv} =} hsvd (@var{sys})
## @deftypefnx{Function File} {@var{hsv} =} hsvd (@var{sys}, @var{"offset"}, @var{alpha})
## Hankel singular values of the stable part of an LTI model.  If no output arguments are
## given, the Hankel singular values are displayed in a plot.
## Uses SLICOT AB13AD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2010
## Version: 0.2

function hsv_r = hsvd (sys, prop = "offset", val = 1e-8)

  if (! nargin)
    print_usage ();
  endif
  
  if (! isa (sys, "lti"))
    error ("hsvd: first argument must be an LTI system");
  endif
  
  if (! strcmp (tolower (prop(1)), "o"))
    error ("hsvd: second argument invalid");
  endif
  
  [a, b, c] = ssdata (sys);

  discrete = ! isct (sys);
  
  if (discrete)
    alpha = 1 - val;
  else
    alpha = - val;
  endif
  
  [hsv, ns] = slab13ad (a, b, c, discrete, alpha);
  
  if (nargout)
    hsv_r = hsv;
  else
    bar ((1:ns) + (rows (a) - ns), hsv);
    title (["Hankel Singular Values of Stable Part of ", inputname(1)]);
    xlabel ("State");
    ylabel ("State Energy");
    grid ("on");
  endif

endfunction


%!shared hsv, hsv_exp
%! a = [ -0.04165  0.0000  4.9200  -4.9200  0.0000  0.0000  0.0000
%!       -5.2100  -12.500  0.0000   0.0000  0.0000  0.0000  0.0000
%!        0.0000   3.3300 -3.3300   0.0000  0.0000  0.0000  0.0000
%!        0.5450   0.0000  0.0000   0.0000 -0.5450  0.0000  0.0000
%!        0.0000   0.0000  0.0000   4.9200 -0.04165 0.0000  4.9200
%!        0.0000   0.0000  0.0000   0.0000 -5.2100 -12.500  0.0000
%!        0.0000   0.0000  0.0000   0.0000  0.0000  3.3300 -3.3300];
%!
%! b = [  0.0000   0.0000
%!       12.5000   0.0000
%!        0.0000   0.0000
%!        0.0000   0.0000
%!        0.0000   0.0000
%!        0.0000   12.500
%!        0.0000   0.0000];
%!
%! c = [  1.0000   0.0000  0.0000   0.0000  0.0000  0.0000  0.0000
%!        0.0000   0.0000  0.0000   1.0000  0.0000  0.0000  0.0000
%!        0.0000   0.0000  0.0000   0.0000  1.0000  0.0000  0.0000];
%!
%! sys = ss (a, b, c);
%! hsv = hsvd (sys);
%!
%! hsv_exp = [2.5139; 2.0846; 1.9178; 0.7666; 0.5473; 0.0253; 0.0246];
%!
%!assert (hsv, hsv_exp, 1e-4);
