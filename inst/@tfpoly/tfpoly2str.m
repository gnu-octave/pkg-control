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
## @deftypefn {Function File} {@var{str} =} tfpoly2str (@var{p})
## @deftypefnx {Function File} {@var{str} =} tfpoly2str (@var{p}, @var{tfvar})
## Return the string of a polynomial with string @var{tfvar} as variable.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function str = tfpoly2str (p, tfvar = "x")

  str = "";

  lp = numel (p.poly);

  if (lp > 0)
    ## first element (highest order)
    a = p.poly(1);

    [a,cs,p1,p2] = __fix_sign__ (a, true);

    if (lp == 1)
      str = [cs, p1, num2str(a, 4), p2];
    else
      if (a == 1)
        str = [cs, __variable__(tfvar, lp-1)];
      else
        str = [cs, p1, num2str(a,4), p2, " ", __variable__(tfvar, lp-1)];
      endif
    endif

    if (lp > 1)
      ## elements in the middle
      for k = 2 : lp-1
        a = p.poly(k);

        if (a != 0)

          [a,cs,p1,p2] = __fix_sign__ (a, false);

          if (a == 1)
            str = [str, cs, __variable__(tfvar, lp-k)];
          else
            str = [str, cs, p1, num2str(a,4), p2, " ", __variable__(tfvar, lp-k)];
          endif
        endif
      endfor

      ## last element (lowest order)
      a = p.poly(lp);

      if (a != 0)
        [a,cs,p1,p2] = __fix_sign__ (a, false);
        str = [str, cs, p1, num2str(a,4), p2];
      endif

    endif
  endif

endfunction


function [a,cs,p1,p2] = __fix_sign__ (a,first)

  if (real (a) < 0)
    if (first)
      cs = "-";
    else
      cs = " - ";
    endif
    a = -a;
  else
    if (first)
      cs = "";
    else
      cs = " + ";
    endif
  endif

  if (imag(a) != 0)
    p1 = "("; p2 = ")";
  else
    p1 = ""; p2 = "";
  endif

endfunction


function str = __variable__ (tfvar, n)

  if (n == 1)
    str = tfvar;
  else
    str = [tfvar, "^", num2str(n)];
  endif

endfunction
