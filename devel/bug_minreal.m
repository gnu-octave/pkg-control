P = ss (-1, 1, 1, 0)
Pi = inv (P)
minreal (Pi)

%{
 ** On entry to TG01JD parameter number 20 had an illegal value
error: __sl_tg01jd__: exception encountered in Fortran subroutine tg01jd_
error: called from:
error:   /Users/lukas/control/inst/@ss/__minreal__.m at line 39, column 19
error:   /Users/lukas/control/inst/@lti/minreal.m at line 38, column 7
%}


ss (inv (tf (Boeing707 )))

%{
octave(211,0x7fff70ddbcc0) malloc: *** error for object 0x108d0e2d8: incorrect checksum for freed object - object was probably modified after being freed.
*** set a breakpoint in malloc_error_break to debug
panic: Abort trap -- stopping myself...
%}