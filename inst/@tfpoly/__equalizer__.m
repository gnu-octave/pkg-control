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
## Make two polynomials equally long by adding leading zeros
## to the shorter one. For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function [peq1, peq2] = __equalizer__ (p1, p2)

  lp1 = length (p1.poly);
  lp2 = length (p2.poly);
  lmax = max (lp1, lp2);

  leadzer1 = zeros (1, lmax - lp1);
  leadzer2 = zeros (1, lmax - lp2);

  peq1 = p1;
  peq2 = p2;

  peq1.poly = [leadzer1, p1.poly];
  peq2.poly = [leadzer2, p2.poly];

endfunction