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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Minimal realization of SS models. The physical meaning of states is lost.

## Algorithm based on sysmin by A. Scottedward Hodel
## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function retsys = __minreal__ (sys, tol)

  A = sys.a;
  B = sys.b;
  C = sys.c;

  if (! isempty (A))
    if (tol == "def")
      [cflg, Uc] = isctrb (A, B);
    else
      [cflg, Uc] = isctrb (A, B, tol);
    endif
    
    if (! cflg)
      if (! isempty (Uc))
        A = Uc.' * A * Uc;
        B = Uc.' * B;
        C = C * Uc;
      else
        A = B = C = [];
      endif
    endif
  endif

  if (! isempty (A))
    if (tol == "def")
      [oflg, Uo] = isobsv (A, C);
    else
      [oflg, Uo] = isobsv (A, C, tol);
    endif

    if (! oflg)
      if (! isempty (Uo))
        A = Uo.' * A * Uo;
        B = Uo.' * B;
        C = C * Uo;
      else
        A = B = C = [];
      endif
    endif
  endif

  retsys = ss (A, B, C, sys.d);
  retsys.lti = sys.lti;  # retain i/o names and tsam

endfunction
