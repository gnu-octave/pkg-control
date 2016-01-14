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
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{z}, @var{p}, @var{k}, @var{tsam}, @dots{})
## Create transfer function model from zero-pole-gain data.
## This is just a stop-gap compatibility wrapper since zpk
## models are not yet implemented.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model to be converted to transfer function.
## @item z
## Cell of vectors containing the zeros for each channel.
## z@{i,j@} contains the zeros from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item p
## Cell of vectors containing the poles for each channel.
## p@{i,j@} contains the poles from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item k
## Matrix containing the gains for each channel.
## k(i,j) contains the gain from input j to output i.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified,
## a continuous-time model is assumed.
## @item @dots{}
## Optional pairs of properties and values.
## Type @command{set (tf)} for more information.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Transfer function model.
## @end table
##
## @seealso{tf, ss, dss, frd}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.2

function sys = zpk (varargin)

  if (nargin <= 1)                  # zpk (),  zpk (sys),  zpk (k),  zpk ('s')
    sys = tf (varargin{:});
    return;
  elseif (nargin == 2 ...
          && ischar (varargin{1}))  # zpk ('z', tsam)
    sys = tf (varargin{:});
    return;
  endif

  z = {}; p = {}; k = [];           # default values
  tsam = 0;                         # default sampling time

  [mat_idx, opt_idx] = __lti_input_idx__ (varargin);

  switch (numel (mat_idx))
    case 1
      k = varargin{mat_idx};
    case 3
      [z, p, k] = varargin{mat_idx};
    case 4
      [z, p, k, tsam] = varargin{mat_idx};
      if (isempty (tsam) && is_real_matrix (tsam))
        tsam = -1;
      elseif (! issample (tsam, -10))
        error ("zpk: invalid sampling time");
      endif
    case 0
      ## nothing to do here, just prevent case 'otherwise'
    otherwise
      print_usage ();
  endswitch

  varargin = varargin(opt_idx);

  if (isempty (z) && isempty (p) && is_real_matrix (k))
    sys = tf (k, varargin{:});
    return;
  endif

  if (! iscell (z))
    z = {z};
  endif

  if (! iscell (p))
    p = {p};
  endif

  if (! size_equal (z, p, k))
    error ("zpk: arguments 'z', 'p' and 'k' must have equal dimensions");
  endif

  ## NOTE: accept [], scalars and vectors but not matrices as 'z' and 'p'
  ##       because  poly (matrix)  returns the characteristic polynomial
  ##       if the matrix is square!

  if (! is_zp_vector (z{:}, 1))  # last argument 1 needed if z is empty cell
    error ("zpk: first argument 'z' must be a vector or a cell of vectors");
  endif

  if (! is_zp_vector (p{:}, 1))
    error ("zpk: second argument 'p' must be a vector or a cell of vectors")
  endif

  if (! is_real_matrix (k))
    error ("zpk: third argument 'k' must be a real-valued gain matrix");
  endif

  num = cellfun (@(zer, gain) real (gain * poly (zer)), z, num2cell (k), "uniformoutput", false);
  den = cellfun (@(pol) real (poly (pol)), p, "uniformoutput", false);

  sys = tf (num, den, tsam, varargin{:});

endfunction
