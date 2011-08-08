## Copyright (C) 2009, 2010, 2011   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} feedback (@var{sys1})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{"+"})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{"+"})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout}, @var{"+"})
## Feedback connection of two LTI models.
##
## @strong{Inputs}
## @table @var
## @item sys1
## LTI model of forward transmission.  @code{[p1, m1] = size (sys1)}.
## @item sys2
## LTI model of backward transmission.
## If not specified, an identity matrix of appropriate size is taken.
## @item feedin
## Vector containing indices of inputs to @var{sys1} which are involved in the feedback loop.
## The number of @var{feedin} indices and outputs of @var{sys2} must be equal.
## If not specified, @code{1:m1} is taken.
## @item feedout
## Vector containing indices of outputs from @var{sys1} which are to be connected to @var{sys2}.
## The number of @var{feedout} indices and inputs of @var{sys2} must be equal.
## If not specified, @code{1:p1} is taken.
## @item "+"
## Positive feedback sign.  If not specified, @var{"-"} for a negative feedback interconnection
## is assumed.  @var{+1} and @var{-1} are possible as well, but only from the third argument
## onward due to ambiguity.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Resulting LTI model.
## @end table
##
## @strong{Block Diagram}
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
## Version: 0.4

function sys = feedback (sys1, sys2, feedin, feedout, fbsign = -1)

  [p1, m1] = size (sys1);

  switch (nargin)
    case 1                          # sys = feedback (sys)
      if (p1 != m1)
        error ("feedback: argument must be a square system");
      endif
      sys2 = eye (p1);
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 2
      if (ischar (sys2))            # sys = feedback (sys, "+")
        if (p1 != m1)
          error ("feedback: first argument must be a square system");
        endif
        fbsign = __check_fbsign__ (sys2);
        sys2 = eye (p1);
      endif                         # sys = feedback (sys1, sys2)
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 3                          # sys = feedback (sys1, sys2, "+")
      fbsign = __check_fbsign__ (feedin);
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 4                          # sys = feedback (sys1, sys2, feedin, feedout)
      ## nothing needs to be done here
      ## case 4 required to prevent "otherwise"

    case 5                          # sys = feedback (sys1, sys2, feedin, feedout, "+")
      fbsign = __check_fbsign__ (fbsign);

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
  M22 = zeros (m2, p2);

  M12 = full (sparse (feedin, 1:l_feedin, fbsign, m1, p2));
  M21 = full (sparse (1:l_feedout, feedout, 1, m2, p1));
  
  ## NOTE: for-loops do NOT the same as
  ##       M12(feedin, 1:l_feedin) = fbsign;
  ##       M21(1:l_feedout, feedout) = 1;
  ##
  ## M12 = zeros (m1, p2);
  ## M21 = zeros (m2, p1);
  ##
  ## for k = 1 : l_feedin
  ##   M12(feedin(k), k) = fbsign;
  ## endfor
  ##
  ## for k = 1 : l_feedout
  ##   M21(k, feedout(k)) = 1;
  ## endfor

  M = [M11, M12;
       M21, M22];

  in_idx = 1 : m1;
  out_idx = 1 : p1;

  sys = __sys_group__ (sys1, sys2);
  sys = __sys_connect__ (sys, M);
  sys = __sys_prune__ (sys, out_idx, in_idx);

endfunction


function fbsign = __check_fbsign__ (fbsign)

  if (is_real_scalar (fbsign))
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
