## Copyright (C) 2009 - 2010   Lukas F. Reichlin
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
## Display routine for TF objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function display (sys)

  [inname, outname] = __lti_data__ (sys.lti);

  m = numel (inname);
  p = numel (outname);

  if (m == 0 || isequal ("", inname{:}))
    inname = strseq ("u", 1:m);
  else
    inname = __mark_empty_names__ (inname);
  endif

  if (p == 0 || isequal ("", outname{:}))
    outname = strseq ("y", 1:p);
  else
    outname = __mark_empty_names__ (outname);
  endif

  disp ("");

  for nu = 1 : m
    disp (["Transfer function \"", inputname(1), "\" from input \"", inname{nu}, "\" to output ..."]);
    for ny = 1 : p
      __disp_frac__ (sys.num{ny, nu}, sys.den{ny, nu}, sys.tfvar, outname{ny});
    endfor
  endfor

  display (sys.lti);

endfunction


function __disp_frac__ (num, den, tfvar, name)

  MAX_LEN = 12;  # max length of output name

  if (num == 0)
    str = [" ", name, ":  0"];
  elseif (den == 1)
    str = [" ", name, ":  "];
    numstr = tfpoly2str (num, tfvar);
    str = [str, numstr];
  ##elseif (length (den) == 1)  # de-comment for non-development use
  ##  str = [" ", name, ":  "];
  ##  num = num * (1/get (den));
  ##  numstr = tfpoly2str (num, tfvar);
  ##  str = [str, numstr];
  else
    numstr = tfpoly2str (num, tfvar);
    denstr = tfpoly2str (den, tfvar);
    fracstr = repmat ("-", 1, max (length (numstr), length (denstr)));

    str = strjust (strvcat (numstr, fracstr, denstr), "center");

    namestr = [name, ":  "];
    namestr = strjust (strvcat (" ", namestr, " "), "left");
    namestr = namestr(:, 1 : min (MAX_LEN, end));
    namestr = horzcat (repmat (" ", 3, 1), namestr);

    str = [namestr, str];
  endif

  disp (str);
  disp ("");

endfunction

