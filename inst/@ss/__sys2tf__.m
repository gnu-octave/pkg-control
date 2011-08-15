## Copyright (C) 2009, 2011   Lukas F. Reichlin
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
## SS to TF conversion.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function [retsys, retlti] = __sys2tf__ (sys)

  [a, b, c, d] = ssdata (sys);                # system could be a descriptor model

  [num, den, ign, igd, md, p, m] = sltb04bd (a, b, c, d);

  num = reshape (num, md, p, m);
  den = reshape (den, md, p, m);

  num = mat2cell (num, md, ones(1,p), ones(1,m));
  den = mat2cell (den, md, ones(1,p), ones(1,m));

  num = squeeze (num);
  den = squeeze (den);

  ign = mat2cell (ign, ones(1,p), ones(1,m));
  igd = mat2cell (igd, ones(1,p), ones(1,m));

  num = cellfun (@(x, y) x(1:y+1), num, ign, "uniformoutput", false);
  den = cellfun (@(x, y) x(1:y+1), den, igd, "uniformoutput", false);

  retsys = tf (num, den, get (sys, "tsam"));  # tsam needed to set appropriate tfvar
  retlti = sys.lti;                           # preserve lti properties
  
  ## FIXME: sys = tf (ss (5))

endfunction
