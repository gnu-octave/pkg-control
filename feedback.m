## feedback(sys1,sys2)
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
 
## Author: Ben Sapp <mailto:bsapp@lanl.gov>
 
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
