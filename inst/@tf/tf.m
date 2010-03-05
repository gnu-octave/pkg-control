## Copyright (C) 2009   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{s} =} tf (@var{"s"})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{num}, @var{den})
## @deftypefnx {Function File} {@var{sys} =} tf (@var{num}, @var{den}, @var{tsam})
## Create or convert to transfer function model.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function sys = tf (num, den, varargin)

  ## model precedence: frd > ss > zpk > tf > double
  %inferiorto ("frd", "ss", "zpk");  # error if de-commented. bug in octave?
  superiorto ("double");

  argc = 0;

  switch (nargin)
    case 0
      num = {};
      den = {};
      tsam = -1;
      tfvar = "x";  # undefined

    case 1
      if (isa (num, "tf"))  # already in tf form
        sys = num;
        return;
      elseif (isa (num, "lti"))  # another lti object
        [sys, numlti] = __sys2tf__ (num);
        sys.lti = numlti;  # preserve lti properties
        return;
      elseif (isnumeric (num))  # static gain
        num = num2cell (num);
        num = __conv2tfpolycell__ (num);
        [p, m] = size (num);
        den = tfpolyones (p, m);
        tsam = -1;
        tfvar = "x";  # undefined
      elseif (ischar (num))  # s = tf ("s")
        tfvar = num;
        num = __conv2tfpolycell__ ([1, 0]);
        den = __conv2tfpolycell__ ([1]);
        tsam = 0;
      else
        print_usage ();
      endif

    case 2
      if (ischar (num) && issample (den))  # z = tf ("z", 0.3)
        tfvar = num;
        tsam = den;
        num = __conv2tfpolycell__ ([1, 0]);
        den = __conv2tfpolycell__ ([1]);
      else  # sys = tf (num, den)
        num = __conv2tfpolycell__ (num);
        den = __conv2tfpolycell__ (den);
        tfvar = "s";
        tsam = 0;
      endif

    otherwise  # default case
      num = __conv2tfpolycell__ (num);
      den = __conv2tfpolycell__ (den);
      argc = numel (varargin);

      if (issample (varargin{1}, 1))  # sys = tf (num, den, tsam, "prop1, "val1", ...)
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
      else  # sys = tf (num, den, "prop1, "val1", ...)
        tsam = 0;
        tfvar = "s";
      endif

  endswitch


  [p, m] = __tfnddim__ (num, den);

  tfdata = struct ("num", {num},
                   "den", {den},
                   "tfvar", tfvar);

  ltisys = lti (p, m, tsam);

  sys = class (tfdata, "tf", ltisys);

  if (argc > 0)
    sys = set (sys, varargin{:});
  endif

endfunction

