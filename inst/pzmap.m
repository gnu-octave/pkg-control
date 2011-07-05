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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} pzmap (@var{sys})
## @deftypefnx {Function File} {[@var{p}, @var{z}] =} pzmap (@var{sys})
## Plot the poles and zeros of an LTI system in the complex plane.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.1

function [pol_r, zer_r] = pzmap (sys)

  ## TODO: multiplot feature:   pzmap (sys1, "b", sys2, "r", ...)

  if (nargin != 1)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("pzmap: argument must be an LTI system");
  endif

  pol = pole (sys);
  zer = zero (sys);

  if (! nargout)
    pol_re = real (pol);
    pol_im = imag (pol);
    zer_re = real (zer);
    zer_im = imag (zer);

    plot (pol_re, pol_im, "xb", zer_re, zer_im, "or")
    grid ("on")  
    title (["Pole-Zero Map of ", inputname(1)])
    xlabel ("Real Axis")
    ylabel ("Imaginary Axis")
  else
    pol_r = pol;
    zer_r = zer;
  endif
  
endfunction