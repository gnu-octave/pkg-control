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
## @deftypefn {Function File} {@var{s} =} tf (@var{'s'})
## @deftypefnx {Function File} {@var{z} =} tf (@var{'z'}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{mat}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{num}, @var{den}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{num}, @var{den}, @var{tsam}, @dots{})
## Create or convert to transfer function model.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model to be converted to transfer function.
## @item mat
## Gain matrix to be converted to static transfer function.
## @item num
## Numerator or cell of numerators.  Each numerator must be a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## num@{i,j@} contains the numerator polynomial from input j to output i.
## In the SISO case, a single vector is accepted as well.
## @item den
## Denominator or cell of denominators.  Each denominator must be a row vector
## containing the coefficients of the polynomial in descending powers of
## the transfer function variable.
## den@{i,j@} contains the denominator polynomial from input j to output i.
## In the SISO case, a single vector is accepted as well.
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
## @strong{Option Keys and Values}
## @table @var
## @item 'num'
## Numerator.  See 'Inputs' for details.
##
## @item 'den'
## Denominator.  See 'Inputs' for details.
##
## @item 'tfvar'
## String containing the transfer function variable.
##
## @item 'inv'
## Logical.  True for negative powers of the transfer function variable.
##
## @item 'tsam'
## Sampling time.  See 'Inputs' for details.
##
## @item 'inname'
## The name of the input channels in @var{sys}.
## Cell vector of length m containing strings.
## Default names are @code{@{'u1', 'u2', ...@}}
##
## @item 'outname'
## The name of the output channels in @var{sys}.
## Cell vector of length p containing strings.
## Default names are @code{@{'y1', 'y2', ...@}}
##
## @item 'ingroup'
## Struct with input group names as field names and
## vectors of input indices as field values.
## Default is an empty struct.
##
## @item 'outgroup'
## Struct with output group names as field names and
## vectors of output indices as field values.
## Default is an empty struct.
##
## @item 'name'
## String containing the name of the model.
##
## @item 'notes'
## String or cell of string containing comments.
##
## @item 'userdata'
## Any data type.
## @end table
##
## @strong{Example}
## @example
## @group
## octave:1> s = tf ('s');
## octave:2> G = 1/(s+1)
##
## Transfer function 'G' from input 'u1' to output ...
## 
##         1  
##  y1:  -----
##       s + 1
## 
## Continuous-time model.
## @end group
## @end example
## @example
## @group
## octave:3> z = tf ('z', 0.2);
## octave:4> H = 0.095/(z-0.9)
## 
## Transfer function 'H' from input 'u1' to output ...
## 
##        0.095 
##  y1:  -------
##       z - 0.9
## 
## Sampling time: 0.2 s
## Discrete-time model.
## @end group
## @end example
## @example
## @group
## octave:5> num = @{[1, 5, 7], [1]; [1, 7], [1, 5, 5]@};
## octave:6> den = @{[1, 5, 6], [1, 2]; [1, 8, 6], [1, 3, 2]@};
## octave:7> sys = tf (num, den)
## @end group
## @end example
## 
## @example
## @group
## Transfer function 'sys' from input 'u1' to output ...
## 
##       s^2 + 5 s + 7
##  y1:  -------------
##       s^2 + 5 s + 6
## 
##           s + 7    
##  y2:  -------------
##       s^2 + 8 s + 6
## @end group
## @end example
## 
## @example
## @group
## Transfer function 'sys' from input 'u2' to output ...
## 
##         1  
##  y1:  -----
##       s + 2
## 
##       s^2 + 5 s + 5
##  y2:  -------------
##       s^2 + 3 s + 2
## 
## Continuous-time model.
## octave:8> 
## @end group
## @end example
##
## @seealso{filt, ss, dss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.4

function sys = tf (varargin)

  ## model precedence: frd > ss > zpk > tf > double
  ## inferiorto ("frd", "ss", "zpk");           # error if de-commented. bug in octave?
  superiorto ("double");

  if (nargin == 1)
    if (isa (varargin{1}, "tf"))                # tf (tfsys)
      sys = varargin{1};
      return;
    elseif (isa (varargin{1}, "lti"))           # tf (ltisys)
      [sys, lti] = __sys2tf__ (varargin{1});
      sys.lti = lti;
      return;
    elseif (ischar (varargin{1}))               # s = tf ('s')
      sys = tf ([1, 0], 1, "tfvar", varargin{:});
      return;
    endif
  elseif (nargin == 2 ...
          && ischar (varargin{1}) ...
          && is_zp_vector (varargin{2}) ...
          && length (varargin{2}) <= 1)         # z = tf ('z', tsam)
    sys = tf ([1, 0], 1, varargin{2}, "tfvar", varargin{[1,3:end]});
    return;
  endif
  
  num = {}; den = {};                           # default transfer matrix
  tsam = -2;                                    # default sampling time

  [mat_idx, opt_idx, obj_flg] = __lti_input_idx__ (varargin);

  switch (numel (mat_idx))
    case 1
      num = varargin{mat_idx};
    case 2
      [num, den] = varargin{mat_idx};
      tsam = 0;
    case 3
      [num, den, tsam] = varargin{mat_idx};
      if (isempty (tsam) && is_real_matrix (tsam))
        tsam = -1;
      elseif (! issample (tsam, -10))
        error ("tf: invalid sampling time");
      endif
    case 0
      ## nothing to do here, just prevent case 'otherwise'
    otherwise
      print_usage ();
  endswitch
  
  varargin = varargin(opt_idx);
  if (obj_flg)
    varargin = horzcat ({"lti"}, varargin);
  endif
  
  [num, den, tsam, tfvar] = __adjust_tf_data__ (num, den, tsam);
  [p, m] = __tf_dim__ (num, den);              # determine number of outputs and inputs

  tfdata = struct ("num", {num},
                   "den", {den},
                   "tfvar", tfvar,
                   "inv", false);              # struct for tf-specific data

  ltisys = lti (p, m, tsam);                   # parent class for general lti data

  sys = class (tfdata, "tf", ltisys);          # create tf object

  if (numel (varargin) > 0)                    # if there are any properties and values, ...
    sys = set (sys, varargin{:});              # use the general set function
  endif

endfunction
