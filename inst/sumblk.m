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
## @deftypefn{Function File} {@var{S} =} sumblk (@var{formula})
## @deftypefnx{Function File} {@var{S} =} sumblk (@var{formula}, @var{n})
## Create summing junction @var{S} from string @var{formula}
## for name-based interconnections.
##
## @strong{Inputs}
## @table @var
## @item formula
## String containing the formula of the summing junction,
## e.g. @code{e = r - y + d} 
## @item n
## Signal size.  Default value is 1.
## @end table
##
## @strong{Outputs}
## @table @var
## @item S
## State-space model of the summing junction.
## @end table
##
## @strong{Example}
## @example
## @group
## octave:1> S = sumblk ('e = r - y + d')
## 
## S.d =
##        r   y   d
##    e   1  -1   1
## 
## Static gain.
## octave:2> S = sumblk ('e = r - y + d', 2)
## 
## S.d =
##        r1  r2  y1  y2  d1  d2
##    e1   1   0  -1   0   1   0
##    e2   0   1   0  -1   0   1
## 
## Static gain.
## @end group
## @end example
## @seealso{connect}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2013
## Version: 0.2

function s = sumblk (formula, n = 1)

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif
  
  if (! ischar (formula))
    error ("sumblk: require string as first argument");
  endif
  
  if (! is_real_scalar (n) || n < 1)
    error ("sumblk: require integer as second argument");
  endif
  
  if (length (strfind (formula, "=")) != 1)
    error ("sumblk: formula requires exactly one '='");
  endif

  ## if the first signal has no sign, add a plus
  if (isempty (regexp (formula, "=\\s*[+-]", "once")))
    formula = regexprep (formula, "=", "=+", "once");
  endif
  
  ## extract operators, remove "=" from formula
  operator = regexp (formula, "[=+-]", "match");
  if (! strcmp (operator{1}, "="))
    error ("sumblk: formula has misplaced '='");
  endif
  operator = operator(2:end);
  formula = formula(formula != "=");
  
  ## extract signal names
  signal = regexp (formula, "\\s*[+-]\\s*", "split");
  if (any (cellfun (@isempty, signal)))
    error ("sumblk: formula is missing some input/output names");
  endif
  outname = signal(1);
  inname = signal(2:end);
  signs = ones (1, numel (inname));
  signs(strcmp (operator, "-")) = -1;
    
  d = kron (signs, eye (n));
  s = ss (d);
  
  ## NOTE: the dark side returns a tf, but i prefer an ss model
  ##       because in general, transfer functions and mimo
  ##       interconnections don't mix well
  
  if (n > 1)
    outgroup = struct (outname{1}, 1:n);
    outname = strseq (outname{1}, 1:n);
    idx = 1 : n*numel (inname);
    idx = reshape (idx, n, []);
    idx = num2cell (idx, 1);
    ingroup = cell2struct (idx, inname, 2);
    tmp = cellfun (@strseq, inname, {1:n}, "uniformoutput", false);
    inname = vertcat (tmp{:});
    s = set (s, "outgroup", outgroup, "ingroup", ingroup);
  endif

  s = set (s, "inname", inname, "outname", outname);

endfunction
