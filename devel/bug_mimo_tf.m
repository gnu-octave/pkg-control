f = tf (Boeing707)

sys = f/f

% this makes more practical sense
g = f \ f

%{
Transfer function 'g' from input 'thrust' to output ...

          1 s^5 + 1.963 s^4 + 1.827 s^3 + 0.6224 s^2 + 0.06773 s + 0.0146
 thrust:  ---------------------------------------------------------------
           s^5 + 1.963 s^4 + 1.827 s^3 + 0.6224 s^2 + 0.06773 s + 0.0146 

                                   -4.464e-13                          
 rudder:  -------------------------------------------------------------
          s^5 + 1.963 s^4 + 1.827 s^3 + 0.6224 s^2 + 0.06773 s + 0.0146

Transfer function 'g' from input 'rudder' to output ...

             -3.286e-14 s^3 - 2.685e-14 s^2 - 6.53e-14 s - 9.663e-14   
 thrust:  -------------------------------------------------------------
          s^5 + 1.963 s^4 + 1.827 s^3 + 0.6224 s^2 + 0.06773 s + 0.0146

          1 s^5 + 1.963 s^4 + 1.827 s^3 + 0.6223 s^2 + 0.06773 s + 0.0146
 rudder:  ---------------------------------------------------------------
           s^5 + 1.963 s^4 + 1.827 s^3 + 0.6224 s^2 + 0.06773 s + 0.0146 

Continuous-time model.
%}

% this result is pretty close to

%{
octave:26> tf (eye (2))

Transfer function 'ans' from input 'u1' to output ...

 y1:  1

 y2:  0

Transfer function 'ans' from input 'u2' to output ...

 y1:  0

 y2:  1

Static gain.
%}

% according to my expectations
