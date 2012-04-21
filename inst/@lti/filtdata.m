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
## @deftypefn {Function File} {[@var{z}, @var{p}, @var{k}, @var{tsam}] =} zpkdata (@var{sys})
## @deftypefnx {Function File} {[@var{z}, @var{p}, @var{k}, @var{tsam}] =} zpkdata (@var{sys}, @var{"v"})
## Access zero-pole-gain data.
##
## @strong{Inputs}
## @table @var
## @item sys
## Any type of LTI model.
## @item "v", "vector"
## For SISO models, return @var{z} and @var{p} directly as column vectors
## instead of cells containing a single column vector.
## @end table
##
## @strong{Outputs}
## @table @var
## @item z
## Cell of column vectors containing the zeros for each channel.
## z@{i,j@} contains the zeros from input j to output i.
## @item p
## Cell of column vectors containing the poles for each channel.
## p@{i,j@} contains the poles from input j to output i.
## @item k
## Matrix containing the gains for each channel.
## k(i,j) contains the gain from input j to output i.
## @item tsam
## Sampling time in seconds.  If @var{sys} is a continuous-time model,
## a zero is returned.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.1

function [num, den, tsam] = filtdata (sys, rtype = "cell")

  [num, den, tsam] = tfdata (sys);
  
  ## make numerator and denominator polynomials equally long
  ## by adding leading zeros
  lnum = cellfun (@length, num, "uniformoutput", false);
  lden = cellfun (@length, den, "uniformoutput", false);

  lmax = cellfun (@max, lnum, lden, "uniformoutput", false);

  num = cellfun (@prepad, num, lmax, "uniformoutput", false);
  den = cellfun (@prepad, den, lmax, "uniformoutput", false);
      
  ## remove trailing zeros
  ## such that polynomials are as short as possible
  num = cellfun (@__remove_trailing_zeros__, num, "uniformoutput", false);
  den = cellfun (@__remove_trailing_zeros__, den, "uniformoutput", false);

  if (lower (rtype(1)) == "v" && issiso (sys))
    num = num{1};
    den = den{1};
  endif

endfunction
