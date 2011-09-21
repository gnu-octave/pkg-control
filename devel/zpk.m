## Copyright (C) 2011   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{s} =} zpk (@var{"s"})
## @deftypefnx {Function File} {@var{z} =} zpk (@var{"z"}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{k})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{z}, @var{p}, @var{k}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} zpk @var{z}, @var{p}, @var{k}, @var{tsam}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} zpk (@var{z}, @var{p}, @var{k}, @var{tsam}, @dots{})
## Create transfer function model from zero-pole-gain data.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model to be converted to transfer function.
## @item num
## Numerator or cell of numerators.  Each numerator must be a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## @item den
## Denominator or cell of denominators.  Each denominator must be a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified, a continuous-time
## model is assumed.
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
## @strong{Example}
## @example
## @group
## @end group
## @end example
##
## @seealso{tf, ss, dss, frd}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.1

function sys = zpk (z = {}, p = {}, k = [], varargin)

  switch (nargin)
    case 0
      sys = tf ();
      return;

    case 1
      if (isa (z, "lti") || is_real_matrix (z) || ischar (z))
        sys = tf (z);
        return;
      else
        print_usage ();
      endif

    case 2
      if (ischar (z) && issample (p, -1))
        sys = tf (z, p);
      else
        print_usage ();
      endif

    otherwise
      if (! iscell (z))
        z = {z};
      endif
      if (! iscell (p))
        p = {p};
      endif
      num = cellfun (@(zer, gain) real (gain * poly (zer)), z, num2cell (k), "uniformoutput", false);
      den = cellfun (@(pol) real (poly (pol)), p, "uniformoutput", false);
      sys = tf (num, den, varargin{:});
    endswitch

endfunction

