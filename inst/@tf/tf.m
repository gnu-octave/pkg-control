## Copyright (C) 2009, 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{s} =} tf (@var{"s"})
## @deftypefnx {Function File} {@var{z} =} tf (@var{"z"}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{num}, @var{den}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{num}, @var{den}, @var{tsam}, @dots{})
## Create or convert to transfer function model.
##
## @strong{Inputs}
## @table @var
## @item num
## Numerator or cell of numerators. Row vector containing the exponents
## of the polynomial in descending order.
## @item den
## Denominator or cell of denominators. Row vector containing the exponents
## of the polynomial in descending order.
## @item tsam
## Sampling time. If @var{tsam} is not specified, a continuous-time
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
## octave:1> s = tf ("s");
## octave:2> G = 1/(s+1)
##
## Transfer function "G" from input "u1" to output ...
##         1  
##  y1:  -----
##       s + 1
##
## octave:3> z = tf ("z", 0.2);
## octave:4> H = 0.095/(z-0.9)
## 
## Transfer function "H" from input "u1" to output ...
##        0.095 
##  y1:  -------
##       z - 0.9
## 
## Sampling time: 0.2 s
## octave:5> 
## @end group
## @end example
##
## @seealso{ss, dss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function sys = tf (num = {}, den = {}, varargin)

  ## model precedence: frd > ss > zpk > tf > double
  ## inferiorto ("frd", "ss", "zpk");      # error if de-commented. bug in octave?
  superiorto ("double");

  argc = 0;

  switch (nargin)
    case 0
      tsam = -1;
      tfvar = "x";                         # undefined

    case 1
      if (isa (num, "tf"))                 # already in tf form
        sys = num;
        return;
      elseif (isa (num, "lti"))            # another lti object
        [sys, numlti] = __sys2tf__ (num);
        sys.lti = numlti;                  # preserve lti properties
        return;
      elseif (is_real_matrix (num))        # static gain
        num = num2cell (num);
        num = __vec2tfpoly__ (num);
        [p, m] = size (num);
        den = tfpolyones (p, m);
        tsam = -1;
        tfvar = "x";                       # undefined
      elseif (ischar (num))                # s = tf ("s")
        tfvar = num;
        num = __vec2tfpoly__ ([1, 0]);
        den = __vec2tfpoly__ ([1]);
        tsam = 0;
      else
        print_usage ();
      endif

    case 2
      if (ischar (num) && issample (den))  # z = tf ("z", 0.3)
        tfvar = num;
        tsam = den;
        num = __vec2tfpoly__ ([1, 0]);
        den = __vec2tfpoly__ ([1]);
      else                                 # sys = tf (num, den)
        num = __vec2tfpoly__ (num);
        den = __vec2tfpoly__ (den);
        tfvar = "s";
        tsam = 0;
      endif

    otherwise                              # default case
      num = __vec2tfpoly__ (num);
      den = __vec2tfpoly__ (den);
      argc = numel (varargin);

      if (issample (varargin{1}, 0))       # sys = tf (num, den, tsam, "prop1, "val1", ...)
        tsam = varargin{1};
        argc--;
        
        if (varargin{1} == 0)
          tfvar = "s";
        elseif (varargin{1} == -1)
          tfvar = "x";
        else
          tfvar = "z";
        endif
        
        if (argc > 0)
          varargin = varargin(2:end);
        endif
      else                                 # sys = tf (num, den, "prop1, "val1", ...)
        tsam = 0;
        tfvar = "s";
      endif

  endswitch


  [p, m] = __tf_dim__ (num, den);

  tfdata = struct ("num", {num},
                   "den", {den},
                   "tfvar", tfvar);

  ltisys = lti (p, m, tsam);

  sys = class (tfdata, "tf", ltisys);

  if (argc > 0)
    sys = set (sys, varargin{:});
  endif

endfunction

