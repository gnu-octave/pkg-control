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
## @deftypefn {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{sys})
## @deftypefnx {Function File} {[@var{num}, @var{den}, @var{tsam}] =} tfdata (@var{sys}, @var{"tfpoly"})
## Access transfer function data.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function [num, den, tsam] = tfdata (sys, rtype = "vector")

  if (! isa (sys, "tf"))
    sys = tf (sys);
  endif

  sysdata = __sys_data__ (sys); 

  num = sysdata.num;
  den = sysdata.den;
  tsam = sys.tsam;

  if (rtype == "vector")
    num = cellfun ("@tfpoly/get", num, "uniformoutput", false);
    den = cellfun ("@tfpoly/get", den, "uniformoutput", false);
  endif

endfunction