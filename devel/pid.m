## -*- texinfo -*-
## @deftypefn{Function File} {@var{C} =} pid (@var{Kp})
## @deftypefnx{Function File} {@var{C} =} pid (@var{Kp}, @var{Ki})
## @deftypefnx{Function File} {@var{C} =} pid (@var{Kp}, @var{Ki}, @var{Kd})
## @deftypefnx{Function File} {@var{C} =} pid (@var{Kp}, @var{Ki}, @var{Kd}, @var{Tf})
## Return the transfer function @var{C} of the @acronym{PID} controller
## in parallel form with first-order roll-off.
##
## @example
## @group
##              Ki      Kd s
## C(s) = Kp + ---- + --------
##              s     Tf s + 1
## @end group
## @end example
## @end deftypefn

function C = pid (Kp = 1, Ki = 0, Kd = 0, Tf = 0)

  if (! is_real_scalar (Kp, Ki, Kd, Tf) || nargin > 4)
    print_usage ();
  endif

  if (Kd == 0)    # catch cases like  pid (2, 0, 0, 3)
    Tf = 0;
  endif

  if (Ki == 0)    # minimal realization if  num(3) == 0  and  den(3) == 0
    C = tf ([Kp*Tf+Kd, Kp], [Tf, 1]);
  else
    C = tf ([Kp*Tf+Kd, Kp+Ki*Tf, Ki], [Tf, 1, 0]);
  endif

endfunction
