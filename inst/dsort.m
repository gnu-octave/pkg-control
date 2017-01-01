## Copyright (C) 2016  Mark Bronsfeld
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{s} =} dsort(@var{p})
## @deftypefnx {Function File} {[@var{s}, @var{ndx}] =} dsort(@var{p})
## Sort discrete-time poles by magnitude (in decreasing order).
##
## @strong{Inputs}
## @table @var
## @item p
## Input vector containing discrete-time poles.
## @end table
##
## @strong{Outputs}
## @table @var
## @item s
## Vector with sorted poles.
## @item ndx
## Vector containing the indices used in the sort.
## @end table
##
## @seealso{eig, esort, pole, pzmap, sort, zero}
## @end deftypefn

## Author: Mark Bronsfeld <m.brnsfld@googlemail.com>
## Created: December 2016
## Version: 0.2

function [s, ndx] = dsort(p)
        if(nargin == 1)
                if(!isvector(p))
                error("dsort: argument must be a vector");
                endif
        else
                print_usage();
        endif

        [s, ndx] = sort(p, 'descend');
endfunction

%!shared s_exp, ndx_exp, s_obs, ndx_obs
%! p = [-0.2410+0.5573i;
%!      -0.2410-0.5573i;
%!      0.1503;
%!      -0.0972;
%!      -0.2590];
%! s_exp = [    -0.2410+0.5573i;
%!              -0.2410-0.5573i;
%!              -0.2590;
%!              0.1503;
%!              -0.0972];
%! ndx_exp = [  1;
%!              2;
%!              5;
%!              3;
%!              4];
%! [s_obs, ndx_obs] = dsort(p);
%!assert(s_obs, s_exp, 0);
%!assert(ndx_obs, ndx_exp, 0);
