num = {[1, 5, 6], [1, 2], [1, 2, -1, 2, 2]; [1, 8, 6], [1, 3, 2], [3, 4, 0]};
den = {[1, 5, 7], [1], [1, -1, 1]; [1, 7], [1, 5, 5], [1, 1]};

%num = {[1, 5, 6], [1, 2], [1, 2, -1, 2, 2]; [1, 8, 6], [1, 3, 2], [3, 4, 0]}.';
%den = {[1, 5, 7], [1], [1, -1, 1]; [1, 7], [1, 5, 5], [1, 1]}.';

%num = {[1, 5, 6], [1, 2]; [1, 8, 6], [1, 3, 2]};
%den = {[1, 5, 7], [1]; [1, 7], [1, 5, 5]};


sys = tf (num, den)

P = ss (sys)


ss (inv (tf (Boeing707 )))

%{
octave:2> test_tf2dss 

Transfer function "sys" from input "u1" to output ...

      s^2 + 5 s + 6
 y1:  -------------
      s^2 + 5 s + 7

      s^2 + 8 s + 6
 y2:  -------------
          s + 7    

Transfer function "sys" from input "u2" to output ...

 y1:  s + 2

      s^2 + 3 s + 2
 y2:  -------------
      s^2 + 5 s + 5

Transfer function "sys" from input "u3" to output ...

      s^4 + 2 s^3 - s^2 + 2 s + 2
 y1:  ---------------------------
              s^2 - s + 1        

      3 s^2 + 4 s
 y2:  -----------
         s + 1   

Continuous-time model.
octave(211,0x7fff70ddbcc0) malloc: *** error for object 0x108d0e2d8: incorrect checksum for freed object - object was probably modified after being freed.
*** set a breakpoint in malloc_error_break to debug
panic: Abort trap -- stopping myself...
%}