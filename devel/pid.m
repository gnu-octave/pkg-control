function C = pid (Kp = 1, Ki = 0, Kd = 0, Tf = 0)

  if (! is_real_scalar (Kp, Ki, Kd, Tf) || nargin > 4)
    print_usage ();
  endif

  if (nargin < 2)
    C = tf (Kp);
  else
    C = tf ([Kp*Tf+Kd, Kp+Ki*Tf, Ki], [Tf, 1, 0]);
  endif

endfunction
