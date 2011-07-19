## Copyright (C) 2010, 2011   Lukas F. Reichlin
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
## Display routine for FRD objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: February 2010
## Version: 0.1

function display (sys)

  sysname = inputname (1);
  [inname, outname, tsam] = __lti_data__ (sys.lti);

  [inname, m] = __labels__ (inname, "u");
  [outname, p] = __labels__ (outname, "y");
  
  w = __freq2str__ (sys.w);

  disp ("");

  for k = 1 : m
    disp (["Frequency response \"", sysname, "\" from input \"", inname{k}, "\" to output ..."]);
    disp ("");
    __disp_resp__ (sys.H(:,k,:), w, outname);
  endfor

  display (sys.lti);  # display sampling time

  if (tsam == -2)
    disp ("Static gain.");
  elseif (tsam == 0)
    disp ("Continuous-time frequency response.");
  else
    disp ("Discrete-time frequency response.");
  endif

endfunction


function __disp_resp__ (H, w, outname)

  p = rows (H);       # number of outputs
  len = size (H, 3);  # number of frequencies

  H = mat2cell (H, ones (1, p), 1, len)(:);
  H = cellfun (@__resp2str__, H, outname, "uniformoutput", false);
  
  tsize = terminal_size ();
  
  ## TODO: Show as many outputs on one line as the terminal width allows
  
  for k = 1 : p
    disp ([w, H{k}]);
    disp ("");
  endfor

endfunction


function str = __freq2str__ (w, title = "w [rad/s]")

  len = rows (w);
  str = __vec2str__ (w);
  line = repmat ("-", 1, max (columns (str), columns (title)));
  str = strvcat (title, line, str);
  space = repmat ("  ", len+2, 1);
  str = [space, str];
  
endfunction


function str = __resp2str__ (H, outname)

  len = length (H);
  real_str = __vec2str__ (real (H), "  ");
  imag_str = __vec2str__ (imag (H), "i");
  str = [real_str, imag_str];
  line = repmat ("-", 1, max (columns (str), columns (outname)));
  str = strvcat (outname, line, str);
  space = repmat ("   ", len+2, 1);
  str = [space, str];

endfunction


function col = __vec2str__ (vec, str)

  vec = num2cell (vec(:));
  col = cellfun (@(x) sprintf ("%.4f", x), vec, "uniformoutput", false);
  col = strjust (char (col), "right");
  if (nargin > 1)
    col = [col, repmat(str, length (vec), 1)];
  endif

endfunction
