## Copyright (C) 2007 ZYED ben Mohamed EL HIDRI.  All rights reserved.
## Copyright (C) 2009 Luca Favatella <slackydeb@gmail.com>
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

##unitfeedback(sys1)
##
## Creates the System Sys2(s)  from the system Sys1(s)
## when we have a negative feedback. 
##               ____________
##     +   e    |            |         
## u --->0----->|   Sys1(s)  |------------> y
##       ^-     |____________|      |
##       |                          |
##       |                          |
##       |                          |
##       --------------<-------------
##
## and Returns :
##            Y(s)             Sys1(s)
## Sys2(s)=  ------  =  -------------------
##            U(s)         1   +    Sys1(s)
##
##
## the feedback is negative.
##               ____________
##              |            |         
## u ---------->|   Sys2(s)  |------------> y
##              |____________|      
##
## This only works for SISO systems.
 
## Author: Zyed El Hidri <zyedm79@yahoo.com>

function out = unitfeedback (sys)
  
  if (nargin != 1)
    error ("only 1 argument accepted"); 	
  endif

  if (! is_siso (sys)) 	
    error("only single input single output systems supported"); 	
  endif
  
  out = sysfeedback (sys);
endfunction


%!test
%!
%! ## open loop system: 1 / s
%! [num den] = sys2tf (unitfeedback (tf (1, [1 0])));
%!
%! ## closed loop system: 1 / (s + 1)
%! num_exp = 1;
%! den_exp = [1 1];
%!
%! assert ([num den], [num_exp den_exp]);
