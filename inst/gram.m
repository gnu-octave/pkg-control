## Copyright (C) 1996, 2000, 2003, 2004, 2005, 2007
##               Auburn University. All rights reserved.
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
## @deftypefn {Function File} {@var{W} =} gram (@var{sys}, @var{mode})
## @deftypefnx {Function File} {@var{Wc} =} gram (@var{a}, @var{b})
## @code{gram (@var{sys}, 'c')} returns the controllability gramian of
## the (continuous- or discrete-time) system @var{sys}.
## @code{gram (@var{sys}, 'o')} returns the observability gramian of the
## (continuous- or discrete-time) system @var{sys}.
## @code{gram (@var{a}, @var{b})} returns the controllability gramian
## @var{Wc} of the continuous-time system @math{dx/dt = a x + b u};
## i.e., @var{Wc} satisfies @math{a Wc + m Wc' + b b' = 0}.
##
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>

                                # TODO: substitute is_stable with isstable

function W = gram (argin1, argin2)

  if (nargin != 2)
    print_usage ();
  else

    if (! ischar (argin2)) ## the function was called as "gram (a, b)"
      a = argin1;
      b = argin2;
      if (! is_stable (a))
        error ("the continuous-time system must have a stable 'a' matrix");
      else

        ## let lyap do the error checking about dimensions
        W = lyap (a', b*b');
      endif

    else ## the function was called as "gram (sys, mode)"
      if (! (strcmp (argin2, 'c') || strcmp (argin2, 'o')))
        print_usage ();
      else
        if (strcmp (argin2, 'c'))
          a = argin1.a;
          b = argin1.b;
        elseif (strcmp (argin2, 'o'))
          a = argin1.a';
          b = argin1.c';
        endif

        if (isct (argin1))
          if (! is_stable (argin1))
            error ("the continuous-time system must be stable");
          else

            ## let lyap do the error checking about dimensions
            W = lyap (a', b*b');
          endif
        elseif (isdt (argin1))
          if (! is_stable (argin1))
            error ("the discrete-time system must be stable");
          else

            ## let dlyap do the error checking about dimensions
            W = dlyap (a, b*b');
          endif
        else
          error ("strange behaviour in isct/isdt: if you can reproduce \
              this, please submit a bug report");
        endif
      endif
    endif

  endif

endfunction


%!test
%! a = [-1 0 0; 1/2 -1 0; 1/2 0 -1];
%! b = [1 0; 0 -1; 0 1];
%! c = [0 0 1; 1 1 0]; ## it doesn't matter what the value of c is
%! Wc = gram (ss (a, b, c), 'c');
%! assert (a * Wc + Wc * a' + b * b', zeros (size (a)))

%!test
%! a = [-1 0 0; 1/2 -1 0; 1/2 0 -1];
%! b = [1 0; 0 -1; 0 1]; ## it doesn't matter what the value of b is
%! c = [0 0 1; 1 1 0];
%! Wo = gram (ss (a, b, c), 'o');
%! assert (a' * Wo + Wo * a + c' * c, zeros (size (a)))

%!test
%! a = [-1 0 0; 1/2 -1 0; 1/2 0 -1];
%! b = [1 0; 0 -1; 0 1];
%! Wc = gram (a, b);
%! assert (a * Wc + Wc * a' + b * b', zeros (size (a)))

%!test
%! a = [-1 0 0; 1/2 1 0; 1/2 0 -1] / 2;
%! b = [1 0; 0 -1; 0 1];
%! c = [0 0 1; 1 1 0]; ## it doesn't matter what the value of c is
%! d = zeros (rows (c), columns (b)); ## it doesn't matter what the value of d is
%! Ts = 0.1; ## Ts != 0
%! Wc = gram (ss (a, b, c, d, Ts), 'c');
%! assert (a * Wc * a' - Wc + b * b', zeros (size (a)), 1e-12)

%!test
%! a = [-1 0 0; 1/2 1 0; 1/2 0 -1] / 2;
%! b = [1 0; 0 -1; 0 1]; ## it doesn't matter what the value of b is
%! c = [0 0 1; 1 1 0];
%! d = zeros (rows (c), columns (b)); ## it doesn't matter what the value of d is
%! Ts = 0.1; ## Ts != 0
%! Wo = gram (ss (a, b, c, d, Ts), 'o');
%! assert (a' * Wo * a - Wo + c' * c, zeros (size (a)), 1e-12)