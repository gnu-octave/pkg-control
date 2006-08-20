## Copyright (C) 2000 Ben Sapp.  All rights reserved.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## This is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
## for more details.

## feedback(sys1,sys2)
##
## Filter the output of sys1 through sys2 and subtract it from the input.
##
##              _____________
##     +        |            |
## u --->0----->|    sys1    |------->
##       |-     |____________|   |
##       |                       |
##       |      _____________    |
##       |      |            |   |
##       -------|    sys2    |----
##              |____________|
##
## This only works for SISO systems.
 
## Author: Ben Sapp <bsapp@lanl.gov>
 
function out = feedback(sys1,sys2)
  if (nargin != 2)
    error("only 2 arguements accepted");
  endif
  if(!is_siso(sys1) || !is_siso(sys2))
    error("only single input single output systems supported");
  endif
  T = sysgroup(sys1,sys2);
  T = sysdup(T,2,[]);
  T = sysscale(T,diag([1,1,-1]),[]);
  T = sysconnect(T,3,1);
  T = sysconnect(T,1,2);
  out = sysprune(T,1,1);
endfunction
