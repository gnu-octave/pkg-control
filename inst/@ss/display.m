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
## Display routine for SS objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function display (sys)

  [inname, outname, tsam] = __lti_data__ (sys.lti);
  stname = sys.stname;

  m = numel (inname);
  p = numel (outname);
  n = numel (stname);

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

  if (n == 0 || isequal ("", stname{:}))
    stname = strseq ("x", 1:n);
  else
    stname = __mark_empty_names__ (stname);
  endif

  disp ("");

  if (! isempty (sys.a))
    __disp_mat__ (sys.a, [inputname(1), ".a"], stname, stname);
    __disp_mat__ (sys.b, [inputname(1), ".b"], stname, inname);
    __disp_mat__ (sys.c, [inputname(1), ".c"], outname, stname);
  endif

  __disp_mat__ (sys.d, [inputname(1), ".d"], outname, inname);

  display (sys.lti);  # display sampling time

  if (tsam == -1)
    disp ("Static gain.");
  elseif (tsam == 0)
    disp ("Continuous-time model.");
  else
    disp ("Discrete-time model.");
  endif

endfunction


function __disp_mat__ (m, mname, rname, cname)

  MAX_LEN = 12;  # max length of row name and column name
  [mrows, mcols] = size (m);

  row_name = strjust (strvcat (" ", rname), "left");
  row_name = row_name(:, 1 : min (MAX_LEN, end));
  row_name = horzcat (repmat (" ", mrows+1, 3), row_name);

  mat = cell (1, mcols);
  
  for k = 1 : mcols
    cname{k} = cname{k}(:, 1 : min (MAX_LEN, end));
    acol = vertcat (cname(k), cellstr (deblank (num2str (m(:, k), 4))));
    mat{k} = strjust (strvcat (acol{:}), "right");
  endfor

  lcols = cellfun ("size", mat, 2);
  lcols_max = 2 + max (horzcat (lcols, 1));

  for k = 1 : mcols
    mat{k} = horzcat (repmat (" ", mrows+1, lcols_max-lcols(k)), mat{k});
  endfor

  tsize = terminal_size ();
  dispcols = max (1, floor ((tsize(2) - columns (row_name)) / lcols_max));
  disprows = max (1, ceil (mcols / dispcols));

  disp ([mname, " ="]);

  for k = 1 : disprows
    disp (horzcat (row_name, mat{1+(k-1)*dispcols : min (mcols, k*dispcols)}));
    disp ("");
  endfor

endfunction