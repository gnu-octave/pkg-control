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

%{
octave:1> ss (inv (tf (Boeing707 )))
i =  1
j =  1
ans =

    2    3    4    5    6    7    8    9   10   11   12   13   14

i =  1
j =  2
ans =

    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14

error: __proper_tf2ss__: subscript indices must be either positive integers or logicals
error: called from:
error:   /Users/lukas/control/inst/@tf/__sys2ss__.m at line 122, column 47
error:   /Users/lukas/control/inst/@tf/__sys2ss__.m at line 54, column 15
error:   /Users/lukas/control/inst/@ss/ss.m at line 165, column 14
octave:1> 
%}