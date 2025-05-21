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
## Return the string of polynomial vector @var{p} with string @var{tfvar^-1}
## as variable.  Note that there is an almost identical function for the
## @command{tfpoly} class which returns a string with @var{tfvar}
## (not @var{tfvar^-1}) as variable.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2012
## Version: 0.1

function str = tfpoly2str (p, tfvar = "x")

  ## TODO: simplify this ugly code

  str = "";

  lp = numel (p);

  if (lp > 0)               # first element (lowest order)
    idx = find (p);         # first non-zero element
    if (isempty (idx))
      str = "0";
      return;
    else
      idx = idx(1);
    endif
    a = p(idx);

    [a,cs,p1,p2] = __fix_sign__ (a, true);

    if (idx == 1)
      str = [cs, p1, num2str(a, 4), p2];
    else
      if (a == 1)
        str = [cs, __variable__(tfvar, idx-1)];
      else
        str = [cs, p1, num2str(a,4), p2, " ", __variable__(tfvar, idx-1)];
      endif
    endif

    if (lp > idx)           # remaining elements of higher order
      for k = idx+1 : lp
        a = p(k);

        if (a != 0)
          [a,cs,p1,p2] = __fix_sign__ (a,false);
          if (a == 1)
            str = [str, cs, __variable__(tfvar, k-1)];
          else
            str = [str, cs, p1, num2str(a,4), p2, " ", __variable__(tfvar, k-1)];
          endif
        endif
      endfor

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

  str = [tfvar, "^-", num2str(n)];

endfunction
