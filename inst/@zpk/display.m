## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Display routine for ZPK objects.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function display (sys)

  sysname = inputname (1);
  [inname, outname, tsam] = __lti_data__ (sys.lti);

  [inname, m] = __labels__ (inname, "u");
  [outname, p] = __labels__ (outname, "y");

  disp ("");
  disp (["Zero/pole/gain model '", sysname, "':"]);

  for nu = 1 : m
    for ny = 1 : p
      z_ch = sys.z{ny, nu};
      p_ch = sys.p{ny, nu};
      k_ch = sys.k(ny, nu);

      disp ("");
      disp ([outname{ny}, " <- ", inname{nu}, ":"]);
      disp (["gain:  ", num2str(k_ch)]);

      if (isempty (z_ch))
        disp ("zeros: (none)");
      else
        disp ("zeros:");
        for i = 1 : length (z_ch)
          disp (["  ", num2str(z_ch(i), 6)]);
        endfor
      endif

      if (isempty (p_ch))
        disp ("poles: (none)");
      else
        disp ("poles:");
        for i = 1 : length (p_ch)
          disp (["  ", num2str(p_ch(i), 6)]);
        endfor
      endif
    endfor
  endfor

  disp ("");
  display (sys.lti);

  if (tsam == -1)
    disp ("Static gain.");
  elseif (tsam == 0)
    disp ("Continuous-time model.");
  else
    disp ("Discrete-time model.");
  endif

endfunction
