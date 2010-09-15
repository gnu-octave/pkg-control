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
## @deftypefn {Function File} {@var{sys} =} feedback (@var{sys1})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{"+"})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{"+"})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout}, @var{"+"})
## Feedback connection of two LTI models.
## @example
## @group
##  u    +         +--------+             y
## ------>(+)----->|  sys1  |-------+------->
##         ^ -     +--------+       |
##         |                        |
##         |       +--------+       |
##         +-------|  sys2  |<------+
##                 +--------+
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function sys = feedback (sys1, sys2_or_fbsign, feedin_or_fbsign, feedout, fbsign)

  [p1, m1] = size (sys1);

  switch (nargin)
    case 1  # sys = feedback (sys)
      if (p1 != m1)
        error ("feedback: argument must be a square system");
      endif

      fbsign = -1;
      sys2 = eye (p1);
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 2
      if (ischar (sys2_or_fbsign))  # sys = feedback (sys, "+")
        if (p1 != m1)
          error ("feedback: argument must be a square system");
        endif

        fbsign = checkfbsign (sys2_or_fbsign);
        sys2 = eye (p1);
        feedin = 1 : m1;
        feedout = 1 : p1;
      else  # sys = feedback (sys1, sys2)
        fbsign = -1;
        sys2 = sys2_or_fbsign;
        feedin = 1 : m1;
        feedout = 1 : p1;
      endif

    case 3  # sys = feedback (sys1, sys2, "+")
      fbsign = checkfbsign (feedin_or_fbsign);
      sys2 = sys2_or_fbsign;
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 4  # sys = feedback (sys1, sys2, feedin, feedout)
      fbsign = -1;
      sys2 = sys2_or_fbsign;
      feedin = feedin_or_fbsign;

    case 5  # sys = feedback (sys1, sys2, feedin, feedout, "+")
      fbsign = checkfbsign (fbsign);
      sys2 = sys2_or_fbsign;
      feedin = feedin_or_fbsign;

    otherwise
      print_usage ();

  endswitch

  [p2, m2] = size (sys2);

  l_feedin = length (feedin);
  l_feedout = length (feedout);

  if (l_feedin != p2)
    error ("feedback: feedin indices: %d, outputs sys2: %d", l_feedin, p2);
  endif

  if (l_feedout != m2)
    error ("feedback: feedout indices: %d, inputs sys2: %d", l_feedout, m2);
  endif

  if (any (feedin > m1 | feedin < 1))
    error ("feedback: range of feedin indices exceeds dimensions of sys1");
  endif

  if (any (feedin > p1 | feedin < 1))
    error ("feedback: range of feedout indices exceeds dimensions of sys1");
  endif

  M11 = zeros (m1, p1);
  M12 = zeros (m1, p2);
  M21 = zeros (m2, p1);
  M22 = zeros (m2, p2);

  for k = 1 : l_feedin
    M12(feedin(k), k) = fbsign;
  endfor

  for k = 1 : l_feedout
    M21(k, feedout(k)) = 1;
  endfor

  M = [M11, M12;
       M21, M22];

  in_idx = 1 : m1;
  out_idx = 1 : p1;

  sys = __sys_group__ (sys1, sys2);
  sys = __sys_connect__ (sys, M);
  sys = __sys_prune__ (sys, out_idx, in_idx);

endfunction


function fbsign = checkfbsign (fbsign)

  if (isnumeric (fbsign))
    fbsign = sign (fbsign);
  elseif (ischar (fbsign))
    if (strcmp (fbsign, "+"))
      fbsign = +1;
    elseif (strcmp (fbsign, "-"))
      fbsign = -1;
    else
      error ("feedback: invalid feedback sign string");
    endif
  else
    error ("feedback: invalid feedback sign type");
  endif

endfunction