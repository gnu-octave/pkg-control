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
## Return frequency response H and frequency vector w.
## If w is empty, it will be calculated by __frequency_vector__.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.7

function [H, w, sty, leg] = __frequency_response__ (caller, args, nout = 0) 

  ## CALLER         | MIMO  | RANGE | CELL  |
  ## ---------------+-------+-------+-------+
  ## bode           | false | std   | false |
  ## bodemag        | false | std   | false |
  ## margin         | false | std   | false |
  ## nichols        | false | ext   | false |
  ## nyquist        | false | ext   | false |
  ## sensitivity    | false | ext   | false |
  ## sigma          | true  | std   | true  |

  mimoflag = false;
  cellflag = false;
  wbounds = "std";
  
  if (strcmp (caller, {"sigma"}))
    mimoflag = true;
    cellflag = true;
  endif
  
  if (any (strcmp (caller, {"nyquist", "nichols", "sensitivity"})))
    wbounds = "ext";
  endif

  sys_idx = cellfun (@isa, args, {"lti"});          # look for LTI models
  frd_idx = cellfun (@isa, args, {"frd"});          # look for FRD models
  w_idx = cellfun (@is_real_vector, args);          # look for frequency vectors
  r_idx = cellfun (@iscell, args);                  # look for frequency ranges {wmin, wmax}
  s_idx = cellfun (@ischar, args);                  # look for strings (style arguments)

  inv_idx = ! (sys_idx | w_idx | r_idx | s_idx);    # look for invalid arguments

  if (any (inv_idx))
    warning ("%s: argument(s) number %s are invalid and are being ignored", ...
             caller, mat2str (find (inv_idx)(:).'));
  endif

  if (nnz (sys_idx) == 0)
    error ("%s: require at least one LTI model", caller);
  endif

  if (nout > 0 && (nnz (sys_idx) > 1 || any (s_idx)))
    evalin ("caller", "print_usage ()");
  endif

  if (! mimoflag && ! all (cellfun (@issiso, args(sys_idx))))
    error ("%s: require SISO systems", caller);
  endif

  if (any (find (s_idx) < find (sys_idx)(1)))
    warning ("%s: strings in front of first LTI model are being ignored", caller);
  endif

  ## determine frequency vectors
  if (any (r_idx))                                  # if there are frequency ranges
    if (nnz (r_idx) > 1)
      warning ("%s: several frequency ranges specified, taking the last one", caller);
    endif
    r = args(r_idx){end};
    if (numel (r) == 2 && issample (r{1}) && issample (r{2}) && r{1} < r{2})
      w = __frequency_vector__ (args(sys_idx), wbounds, r{1}, r{2});
    else
      error ("%s: the cell defining the desired frequency range is invalid", caller);
    endif
  elseif (any (w_idx))                              # are there any frequency vectors?
    if (nnz (r_idx) > 1)
      warning ("%s: several frequency vectors specified, taking the last one", caller);
    endif
    w = args(w_idx){end};
    w = repmat ({w}, 1, nnz (sys_idx));
  else                                              # there are neither frequency ranges nor vectors    
    w = __frequency_vector__ (args(sys_idx), wbounds);
  endif

  ## temporarily save frequency vectors of FRD models
  w_frd = cellfun (@(x) get (x, "w"), args(frd_idx), "uniformoutput", false);
  w(frd_idx) = {[]};                                # freqresp returns all frequencies of FRD models for w=[]

  ## compute frequency response H for all LTI models
  H = cellfun (@__freqresp__, args(sys_idx), w, {cellflag}, "uniformoutput", false);

  ## restore frequency vectors of FRD models in w
  w(frd_idx) = w_frd;

  ## extract plotting styles
  tmp = cumsum (sys_idx);
  tmp(sys_idx | ! s_idx) = 0;
  n = nnz (sys_idx);
  sty = arrayfun (@(x) args(tmp == x), 1:n, "uniformoutput", false);

  ## get system names for legend
  ## "''" needed for  bode (lticell{:})
  leg = cell (1, n);
  idx = find (sys_idx);
  for k = 1:n
    leg{k} = evalin ("caller", sprintf ("inputname(%d)", idx(k)), "''");
  endfor

endfunction
