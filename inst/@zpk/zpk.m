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
## @deftypefn {Function File} {@var{s} =} zpk (@var{'s'})
## @deftypefnx {Function File} {@var{z} =} zpk (@var{'z'}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{k}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{z}, @var{p}, @var{k}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{z}, @var{p}, @var{k}, @var{tsam}, @dots{})
## Create zero-pole-gain model.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model to be converted to zero-pole-gain form.
## @item z
## Cell of vectors containing the zeros for each channel.
## z@{i,j@} contains the zeros from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item p
## Cell of vectors containing the poles for each channel.
## p@{i,j@} contains the poles from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item k
## Matrix of gains.  k(i,j) contains the gain from input j to output i.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified, a
## continuous-time model is assumed.
## @end table
## @strong{Outputs}
## @table @var
## @item sys
## Zero-pole-gain model.
## @end table
## @strong{Option Keys and Values}
## @table @var
## @item 'InputName'
## String or cell array of strings specifying the names of the inputs.
## @item 'OutputName'
## String or cell array of strings specifying the names of the outputs.
## @end table
## Type @code{set (zpk)} for further options.
##
## @seealso{@@tf/tf, @@ss/ss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.2

function sys = zpk (varargin)

  ## model precedence: zpk > tf > double
  superiorto ("tf", "double");

  if (nargin >= 1 && ischar (varargin{1}))   # zpk ('s'),  zpk ('z', tsam)
    ## transfer-function variable shorthand. still returns a TF model
    sys = tf (varargin{:});
    return;
  endif

  if (nargin == 1 && isa (varargin{1}, "zpk"))     # zpk (zpksys)
    sys = varargin{1};
    return;
  endif

  if (nargin == 1 && isa (varargin{1}, "lti"))     # zpk (ltisys)

   ## conversion from tf/ss/frd
   [num, den, tsam] = tfdata (varargin{1});
   num = cellfun (@__remove_leading_zeros__, num, 'uniformoutput', false);
   den = cellfun (@__remove_leading_zeros__, den, 'uniformoutput', false);

   z = cellfun (@roots, num, "uniformoutput", false);
   p = cellfun (@roots, den, "uniformoutput", false);
   k = cellfun (@(n,d)  n(1)/d(1), num, den);

  varargin = {};

  else

    z = {};  p = {};  k = [];                # default values
    tsam = 0;                                # default sampling time

    [mat_idx, opt_idx] = __lti_input_idx__ (varargin);

    switch (numel (mat_idx))
      case 0
        tsam = -1;                           # empty model
      case 1                                 # zpk (k)  ->  static gain
        k = varargin{mat_idx};
        tsam = -1;
      case 3                                 # zpk (z, p, k)
        [z, p, k] = varargin{mat_idx};
      case 4                                 # zpk (z, p, k, tsam)
        [z, p, k, tsam] = varargin{mat_idx};
        if (isempty (tsam) && is_real_matrix (tsam))
          tsam = -1;
        elseif (! issample (tsam, -10))
          error ("zpk: invalid sampling time");
        endif
      otherwise
        print_usage ();
    endswitch

    varargin = varargin(opt_idx);

    if (! iscell (z))
      z = {z};
    endif

    if (! iscell (p))
      p = {p};
    endif

    ## static gain: one empty zero/pole vector per channel of k
    if (isempty (z) && isempty (p) && ! isempty (k))
      z = repmat ({[]}, size (k));
      p = repmat ({[]}, size (k));
    endif

    if (! size_equal (z, p, k))
      error ("zpk: arguments 'z', 'p' and 'k' must have equal dimensions");
    endif

    if (! is_zp_vector (z{:}, 1))  # last argument 1 needed if z is empty cell
      error ("zpk: first argument 'z' must be a vector or a cell of vectors");
    endif

    if (! is_zp_vector (p{:}, 1))
      error ("zpk: second argument 'p' must be a vector or a cell of vectors");
    endif

    if (! is_real_matrix (k))
      error ("zpk: third argument 'k' must be a real-valued gain matrix");
    endif

  endif

  ## ensure column vectors; store data verbatim, no polynomial conversion
  z = cellfun (@(v) v(:), z, "uniformoutput", false);
  p = cellfun (@(v) v(:), p, "uniformoutput", false);

  [p_out, m_in] = size (k);

  zdata = struct ("z", {z}, "p", {p}, "k", k);
  ltisys = lti (p_out, m_in, tsam);
  sys = class (zdata, "zpk", ltisys);

  if (numel (varargin) > 0)
    sys = set (sys, varargin{:});
  endif

endfunction


%!test
%! sys = zpk ([], [-1], 1);
%! assert (isa (sys, 'zpk'));
%! assert (get (sys, 'tsam'), 0);

%!test
%! sys = zpk ([], [-0.5], 2, 0.1);
%! assert (isa (sys, 'zpk'));
%! assert (get (sys, 'tsam'), 0.1);

%!test
%! ze = {[1]; [-2; 0]};
%! pe = {[-1; 0]; [-4; -3; -1]};
%! ke = [5; 10];
%! sys = zpk (ze, pe, ke);
%! assert (isa (sys, 'zpk'));
%! assert (size (sys), [2, 1]);

%!test
%! ## static gain
%! sys = zpk (5);
%! assert (isa (sys, 'zpk'));
%! assert (isstaticgain (sys));

%!test
%! ## conversion from tf keeps class
%! sys = zpk (tf (1, [1 3 2]));
%! assert (isa (sys, 'zpk'));
%! [~, p] = zpkdata (sys, 'v');
%! assert (sort (real (p)), [-2; -1], 1e-10);
