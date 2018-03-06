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
## @deftypefn {Function File} {@var{p} =} pole (@var{sys})
## Compute poles of @acronym{LTI} system.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item p
## Poles of @var{sys}.
## @end table
##
## @strong{Algorithm}@*
## For (descriptor) state-space models and system/state matrices, @command{pole}
## relies on Octave's @command{eig}.
## For @acronym{SISO} transfer functions, @command{pole}
## uses Octave's @command{roots}.
## @acronym{MIMO} transfer functions are converted to
## a @emph{minimal} state-space representation for the
## computation of the poles.
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Contributor: Mark Bronsfeld <m.brnsfld@googlemail.com>
## Created: October 2009
## Version: 0.2

function pol = pole (sys)

   if(nargin == 1) # pole(sys)
     if(!(isa(sys, "lti")) && issquare(sys))
       pol = eig(sys);
     elseif(isa(sys, "lti"))
       pol = __pole__(sys);
     else
       error("pole: argument must be an LTI system");
     endif
   else
     print_usage();
   endif

endfunction

 
%!shared pol_exp, pol_obs
%! A = [-1, 0,  0; 
%!          0,  -2, 0; 
%!          0,  0,  -3];
%! pol_exp = [-3; 
%!                   -2; 
%!                   -1];
%! pol_obs = pole(A);
%!assert(pol_obs, pol_exp, 0);

## Poles of descriptor state-space model
%!shared pol, pol_exp, infp, kronr, kronl, infp_exp, kronr_exp, kronl_exp
%! A = [  1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     1     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     1     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0
%!        0     0     0     0     0     0     0     0     1 ];
%!
%! E = [  0     0     0     0     0     0     0     0     0
%!        1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0 ];
%!
%! B = [ -1     0     0
%!        0     0     0
%!        0     0     0
%!        0    -1     0
%!        0     0     0
%!        0     0     0
%!        0     0    -1
%!        0     0     0
%!        0     0     0 ];
%!
%! C = [  0     1     1     0     3     4     0     0     2
%!        0     1     0     0     4     0     0     2     0
%!        0     0     1     0    -1     4     0    -2     2 ];
%!
%! D = [  1     2    -2
%!        0    -1    -2
%!        0     0     0 ];
%!
%! sys = dss (A, B, C, D, E, "scaled", true);
%! [pol, ~, infp, kronr, kronl] = __sl_ag08bd__ (A, E, [], [], [], true);
%!
%! pol_exp = zeros (0,1);
%!
%! infp_exp = [0, 3];
%! kronr_exp = zeros (1,0);
%! kronl_exp = zeros (1,0);
%!
%!assert (pol, pol_exp, 1e-4);
%!assert (infp, infp_exp);
%!assert (kronr, kronr_exp);
%!assert (kronl, kronl_exp);
