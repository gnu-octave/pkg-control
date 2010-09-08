## ==============================================================================
## The tests below usually segfault after consecutive calls
## ==============================================================================

test test_minreal

## ss: minreal
%!shared Ar, Br, Cr
%!
%! A = [ -2  -2.4
%!        0  -4.4 ];
%!
%! B = [ 0.6
%!       0.6 ];
%!
%! C = [ 4  -4 ];
%!
%! [Ar, Br, Cr, nr] = sltb01pd (A, B, C, 0);
%!
%!assert (isempty (Ar), true);
%!assert (isempty (Br), true);
%!assert (isempty (Cr), true);


## ss: minreal (SLICOT TB01PD)
%!shared M, Me
%! A = [ 1.0   2.0   0.0
%!       4.0  -1.0   0.0
%!       0.0   0.0   1.0 ];
%!
%! B = [ 1.0
%!       0.0
%!       1.0 ];
%!
%! C = [ 0.0   1.0  -1.0
%!       0.0   0.0   1.0 ];
%!
%! D = zeros (2, 1);
%!
%! [Ar, Br, Cr, nr] = sltb01pd (A, B, C, 0.0);
%! M = [Ar, Br; Cr, D];
%!
%! Ae = [ 1.0000  -1.4142   1.4142
%!       -2.8284   0.0000   1.0000
%!        2.8284   1.0000   0.0000 ];
%!
%! Be = [-1.0000
%!        0.7071
%!        0.7071 ];
%!
%! Ce = [ 0.0000   0.0000  -1.4142
%!        0.0000   0.7071   0.7071 ];
%!
%! De = zeros (2, 1);
%!
%! Me = [Ae, Be; Ce, De];
%!
%!assert (M, Me, 1e-4);