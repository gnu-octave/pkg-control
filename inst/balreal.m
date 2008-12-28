## Copyright (C) 2008 Luca Favatella <slackydeb@gmail.com>
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
## @deftypefn {Function File} {[@var{sysb}, @var{g}] =} balreal (@var{sys})
## @deftypefnx {Function File} {[@var{sysb}, @var{g}, @var{T}, @var{Ti}] =} balreal (@var{sys})
## Balanced realization of the continuous-time LTI system @var{sys}.
##
## @strong{Input}
## @table @var
## @item sys
## Stable, controllable and observable continuous-time LTI system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sysb
## Balanced realization of @var{sys}.
## @item g
## Diagonal of the balanced gramians.
## @item T
## State transformation to convert @var{sys} to @var{sysb}.
## @item Ti
## Inverse of T.
## @end table
##
## @seealso{gram}
## @end deftypefn

## Author: Luca Favatella <slackydeb@gmail.com>
## Version: 0.2.4

                                # TODO
                                # improve continuous-time system test
                                # test against discrete-time systems
                                # clean asserts
                                # clean names Wc2/Wc and Wo2/Wc
                                # substitute is_stable with isstable

function [sysb, g, T, Ti] = balreal (sys)

  if (nargin != 1)
    print_usage ();
  elseif (! is_stable (sys))
    error ("sys must be stable");
  elseif (! is_controllable (sys))
    error ("sys must be controllable");
  elseif (! is_observable (sys))
    error ("sys must be observable");
  else
    ASSERT_TOL = 1e-12; ## DEBUG

    ## step 1: compute T1
    Wc2 = gram (sys, 'c');
    [Vc, Sc2, trash] = svd (Wc2);
    assert (trash, Vc, ASSERT_TOL); ## DEBUG
    T1 = Vc * sqrt (Sc2);

    ## step 2: apply T1
    Atilde = inv (T1) * sys.a * T1;
    Btilde = inv (T1) * sys.b;
    Ctilde = sys.c * T1;
    Wc2tilde = gram (Atilde, Btilde); ## DEBUG
    assert (Wc2tilde, eye (size (Wc2tilde)), ASSERT_TOL); ## DEBUG

    ## step 3: compute T2
    Wo2tilde = gram (Atilde', Ctilde');
    [Vo, Sigma2, trash] = svd (Wo2tilde);
    assert (trash, Vo, ASSERT_TOL); ## DEBUG
    T2 = Vo * (Sigma2 ^ (- 1/4));

    ## step 4: apply T2
    Abal = inv (T2) * Atilde * T2;
    Bbal = inv (T2) * Btilde;
    Cbal = Ctilde * T2;


    ## prepare return values
    ## sysb
    sysb = ss (Abal, Bbal, Cbal, sys.d);

    ## g
    Wc2bal = gram (Abal, Bbal);
    g = diag (Wc2bal);
    Wo2bal = gram (Abal', Cbal'); ## DEBUG
    assert (Wc2bal, Wo2bal, ASSERT_TOL); ## DEBUG
    Wo2 = gram (sys, 'o'); ## DEBUG
    Sigma = diag (sqrt (sort (eig (Wc2 * Wo2), 'descend'))); ## DEBUG
    assert (Wc2bal, Sigma, ASSERT_TOL); ## DEBUG
    assert (Wo2bal, Sigma, ASSERT_TOL); ## DEBUG

    ## T and Ti
    T = inv (T1 * T2);
    Ti = inv (T);
  endif
endfunction


%!test
%! a = [-1 0 0; 1/2 -1 0; 1/2 0 -1];
%! b = [1 0; 0 -1; 0 1];
%! c = [0 0 1; 1 1 0];
%! balreal (ss (a, b, c));