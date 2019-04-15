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
## Conjugate transpose of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2012
## Version: 0.1

function sys = __ctranspose__ (sys, ct)

  num = sys.num;
  den = sys.den;
  
  if (ct)   # continuous-time
    num = cellfun (@conj_ct, num, "uniformoutput", false);
    den = cellfun (@conj_ct, den, "uniformoutput", false);
  else      # discrete-time
    ## Both num and den must be the same length 
    ## for the flip command,
    ## so add leading zeros to the shortest one.
    ng=get(num{1});
    dg=get(den{1});
    nng=numel(ng);
    ndg=numel(dg);
    if (nng>ndg)
      dg=  [zeros(1,nng-ndg) dg];
    else
      ng= [zeros(1,ndg-nng) ng] ;
    endif
    ## the flip effectively replaces every z with 1/z 
    num={tfpoly(flip(ng))};
    den={tfpoly(flip(dg))}; 
  endif  

  sys.num = num.';
  sys.den = den.';

endfunction
