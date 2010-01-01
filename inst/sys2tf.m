## Copyright (C) 1996, 2000, 2002, 2004, 2005, 2006, 2007
##               Auburn University.  All rights reserved.
##
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{num}, @var{den}, @var{tsam}, @var{inname}, @var{outname}] =} sys2tf (@var{sys})
## Extract transfer function data from a system data structure.
##
## See @command{tf} for parameter descriptions.
##
## @strong{Example}
## @example
## octave:1> sys=ss([1 -2; -1.1,-2.1],[0;1],[1 1]);
## octave:2> [num,den] = sys2tf(sys)
## num = 1.0000  -3.0000
## den = 1.0000   1.1000  -4.3000
## @end example
## @end deftypefn

## Author: R. Bruce Tenison <btenison@eng.auburn.edu>
## Created: June 24, 1994
## modified to make sys2tf by A. S. Hodel Aug 1995
## modified again for updated system format by John Ingram July 1996

function [num, den, tsam, inname, outname] = sys2tf (Asys)

  if (nargin != 1)
    print_usage ();
  endif

  if (! isstruct (Asys))
    error ("Asys must be a system data structure (see ss, tf, zp)");
  elseif (! is_siso (Asys))
    [n, nz, m, p] = sysdimensions (Asys);
    error ("system is not SISO: %d inputs, %d outputs", m, p);
  endif

  Asys = sysupdate (Asys, "tf");          # just in case

  num = Asys.num;
  den = Asys.den;

  tsam = sysgettsam (Asys);
  inname = sysgetsignals (Asys, "in");
  outname = sysgetsignals (Asys, "out");

endfunction


## This test case was suggested by Christian Riesch.
## An unexpected result is shown. It is not clear if the bug is in tf,
## sysmult, feedback, sys2tf.
##
## Reference:
## http://sourceforge.net/mailarchive/message.php?msg_name=95DC1AA8EC908B48939B72CF375AA5E30E2AED15%40alice.at.omicron.at
##
%!shared sysUNIT, numA, denA, sysB, numCL_exp, denCL_exp
%! sysUNIT = tf (1, 1, 1);
%! numA = [1 -.7];
%! denA = [1 -2 1];
%! sysB = sysmult (tf ([8e-9 -5.6e-9], [1 -1], 1),
%!                 tf ([125e6], [1 -1], 1));
%! numCL_exp = [1 -2 1];
%! denCL_exp = [1 -1 0.3];
%!test
%! [num_obs den_obs] = sys2tf (sysB);
%! assert ([num_obs den_obs], [numA denA], 1e-8);
%!test
%! [num_obs den_obs] = sys2tf (feedback (sysUNIT, tf (numA, denA, 1)));
%! assert ([num_obs den_obs], [numCL_exp denCL_exp], 1e-12);
%!xtest
%! ## If this code produces no error, probably the bug was fixed forgetting
%! ## to remove this test case.
%! [num_obs den_obs] = sys2tf (feedback (sysUNIT, sysB));
%! assert ([num_obs den_obs], [numCL_exp denCL_exp]);