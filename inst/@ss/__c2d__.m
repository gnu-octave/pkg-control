## Copyright (C) 2009-2016   Lukas F. Reichlin
## Copyright (C)             Torsten Lilge
##
## This file is part of LTI Syncope.
##
## LTI Syncope is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LTI Syncope is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Convert the continuous SS model into its discrete-time equivalent.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009

function sys = __c2d__ (sys, tsam, method = "zoh", w0 = 0)

  switch (method(1))
    case {"z", "s"}                # {"zoh", "std"}
      [sys.a, sys.b, sys.c, sys.d, sys.e] = __dss2ss__ (sys.a, sys.b, sys.c, sys.d, sys.e);
      [sys.a, tmp] = __sl_mb05nd__ (sys.a, tsam, eps); 
      sys.b = tmp * sys.b;         # G

    case {"f"}            # {"foh"}
                          # http://people.duke.edu/~hpgavin/cee541/LTI.pdf page 33
      [sys.a, sys.b, sys.c, sys.d, sys.e] = __dss2ss__ (sys.a, sys.b, sys.c, sys.d, sys.e);
      n = size (sys.a,1);
      m = size (sys.b,2);
      M = [ sys.a      sys.b      zeros(n,m);...
            zeros(m,n) zeros(m,m) eye(m,m);...
            zeros(m,n) zeros(m,m) zeros(m,m) ];
      expM = __sl_mb05nd__ (M, tsam, eps);
      sys.a = expM(1:n,1:n);
      Bd    = expM(1:n,n+1:n+m);
      Bdd   = expM(1:n,n+m+1:end);
      Bd0   = Bd - Bdd/tsam;      # input matrix for u(k)
      Bd1   = Bdd/tsam;           # input matrix for u(k+1)
      sys.b = Bd0 + sys.a*Bd1;    # input matrix for ss with new state z(k) = x(k) - Bd1*u(k)
      sys.d = sys.d + sys.c*Bd1;  # throughput for ss with new state z
      sys = ss (sys.a, sys.b, sys.c, sys.d, tsam, 'userdata', Bd1); # store Bd1 in userdata

    case {"t", "b", "p"}           # {"tustin", "bilin", "prewarp"}
      if (method(1) == "p")        # prewarping
        beta = w0 / tan (w0*tsam/2);
      else
        beta = 2/tsam;
      endif
      if (isempty (sys.e))
        [sys.a, sys.b, sys.c, sys.d] = __sl_ab04md__ (sys.a, sys.b, sys.c, sys.d, 1, beta, false);
      else
        [sys.a, sys.b, sys.c, sys.d, sys.e] = __dss_bilin__ (sys.a, sys.b, sys.c, sys.d, sys.e, beta, false);
      endif

    case "i"                       # "impulse"

      if (any (sys.d(:)))
        error ("c2d: impuls invariant discrete-time models only supported for systems without direct feedthrough\n");
      endif

      ## cont-time: u(t) = delta(t)
      ##            x(0) = Phi(0)*B = B  (Phi(0) = 0)
      ##            x(t) = Phi(t)*B
      ##            y(0) = C*B
      ##            y(t) = C*x(t)
      ## disc-time  x(k+1) = Ad*x(k) + Bd*u(k)
      ##            y(k)   = Cd*x(k) + Dd*u(k)
      ##        =>  Ad = Phi(T),  Bd = Phi(T)*B*T, Cd = C, Dd = C*B*T
      [sys.a, sys.b, sys.c, sys.d, sys.e] = __dss2ss__ (sys.a, sys.b, sys.c, sys.d, sys.e);
      [sys.a, tmp] = __sl_mb05nd__ (sys.a, tsam, eps); 
      sys.d = sys.c * sys.b * tsam;
      sys.b = sys.a * sys.b * tsam;
      
    case "m"                       # "matched"
      tmp = ss (c2d (zpk (sys), tsam, method));
      sys.e = tmp.e;
      sys.a = tmp.a;
      sys.b = tmp.b;
      sys.c = tmp.c;
      sys.d = tmp.d;

    otherwise
      error ("ss: c2d: '%s' is an invalid or missing method", method);
  endswitch

endfunction
