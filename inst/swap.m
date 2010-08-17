## This sophisticated function is needed for Octave 3.2.x
## because krylov.m has been moved from control-1.0.x
## to Octave core whereas swap.m remained in control-1.0.x.
## krylov.m is currently in use by sminreal and minreal
## for state-space models.

function [b, a] = swap (a, b)

endfunction
