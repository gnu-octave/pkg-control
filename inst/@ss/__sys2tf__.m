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
## SS to TF conversion.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.7

function [retsys, retlti] = __sys2tf__ (sys)

  sg_flag = false;                              # static gain flag

  try
    [a, b, c, d, tsam, scaled] = ssdata (sys);  # system could be a descriptor model

    if (tsam == -2 || isempty (a))              # static gain
      sg_flag = true;
    else
      [num, den] = __sl_tb04bd__ (a, b, c, d, scaled);
    endif
  catch
    if (strcmp (lasterror.identifier, "dss:improper"))
      ## sys is non-realizable, therefore ssdata failed.
      if (issiso (sys))
        [num, den] = __siso_ss2tf__ (sys);
      else
        [p, m] = size (sys);
        num = cell (p, m);
        den = cell (p, m);
        for i = 1 : p
          for j = 1 : m
            idx = substruct ("()", {i, j});
            tmp = subsref (sys, idx);             # extract siso model
            tmp = sminreal (tmp);                 # sminreal is more suitable than minreal here
            [n, d] = __siso_ss2tf__ (tmp);
            num(i, j) = n;
            den(i, j) = d;
          endfor
        endfor
      endif
      tsam = get (sys, "tsam");
    else
      ## something else happened
      rethrow (lasterror);
    endif
  end_try_catch

  if (sg_flag)
    retsys = tf (d);
  else
    retsys = tf (num, den, tsam);               # tsam needed to set appropriate tfvar
  endif
  retlti = sys.lti;                             # preserve lti properties

endfunction


function [num, den] = __siso_ss2tf__ (sys)

  if (isempty (sys.a))                          # static gain
    num = sys.d;
    den = 1;
  else                                          # default case
    [zer, gain] = zero (sys);
    pol = pole (sys);
    num = gain * real (poly (zer));
    den = real (poly (pol));
  endif

endfunction
