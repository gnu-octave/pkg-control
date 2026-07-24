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
## @deftypefn {Function File} {@var{sys} =} append (@var{sys1}, @var{sys2}, @dots{}, @var{sysN})
## Group @acronym{LTI} models by appending their inputs and outputs.
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function sys = append (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  sys = blkdiag (varargin{:});

  any_delay = any (cellfun (@(s) isa (s, "lti") && hasdelay (s), varargin));

  if (any_delay)
    indelay = cell (nargin, 1);
    outdelay = cell (nargin, 1);
    iodelay = cell (nargin, 1);

    for k = 1 : nargin
      if (isa (varargin{k}, "lti"))
        [indelay{k}, outdelay{k}, iodelay{k}] = get (varargin{k}, "inputdelay", "outputdelay", "iodelay");
      else
        [p, m] = size (varargin{k});
        indelay{k} = zeros (m, 1);
        outdelay{k} = zeros (p, 1);
        iodelay{k} = zeros (p, m);
      endif
    endfor

    sys = set (sys, "InputDelay", vertcat (indelay{:}), ...
                    "OutputDelay", vertcat (outdelay{:}), ...
                    "IODelay", blkdiag (iodelay{:}));
  endif

endfunction


%!test  # 3-operand append: delays carry through per-block, no cross terms
%! s1 = tf (1, [1 1], "InputDelay", 0.1, "OutputDelay", 0.2);
%! s2 = tf (1, [1 2], "InputDelay", 0.3, "OutputDelay", 0.4);
%! s3 = tf (1, [1 3]);
%! s = append (s1, s2, s3);
%! assert (s.InputDelay, [0.1; 0.3; 0], 1e-10);
%! assert (s.OutputDelay, [0.2; 0.4; 0], 1e-10);
%! assert (s.IODelay, zeros (3, 3), 1e-10);

%!test  # no delay on any operand: result has no delay (regression)
%! s = append (tf (1, [1 1]), tf (1, [1 2]));
%! assert (hasdelay (s), false);

%!test  # mix of dynamic (with delay) and static-gain (numeric) operands
%! s1 = tf (1, [1 1], "InputDelay", 0.1, "OutputDelay", 0.2);
%! s = append (s1, eye (2));
%! assert (s.InputDelay, [0.1; 0; 0], 1e-10);
%! assert (s.OutputDelay, [0.2; 0; 0], 1e-10);
%! assert (s.IODelay, zeros (3, 3), 1e-10);
