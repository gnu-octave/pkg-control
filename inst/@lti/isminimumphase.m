## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{bool} =} isminimumphase (@var{sys})
## @deftypefnx {Function File} {@var{bool} =} isminimumphase (@var{sys}, @var{tol})
## Determine whether @acronym{LTI} system has asymptotically stable zero dynamics.
## According to the definition of Byrnes/Isidori [1], the zeros
## of a minimum-phase system must be strictly
## inside the left complex half-plane (continuous-time case)
## or inside the unit circle (discrete-time case).
## Note that the poles are not tested.
##
## M. Zeitz [2] discusses the inconsistent definitions of the minimum-phase property
## in a German paper.  The abstract in English states the following [2]:
## 
## Originally, the minimum phase property has been defined by H. W. Bode [3]
## in order to characterize the unique relationship between gain and phase of
## the frequency response.
## With regard to the design of digital filters, another definition of minimum
## phase is used and a filter is said to be minimum phase if both the filter
## and its inverse are asymptotically stable.
## Finally, systems with asymptotically stable zero dynamics are named as
## minimum phase by C. I. Byrnes and A. Isidori [1].
## Due to the inconsistent definitions, avoiding the minimum phase property
## for control purposes is advocated and the well-established criteria of
## Hurwitz or Ljapunow to describe the stability of filters and zero dynamics
## are recommended.
## 
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.
## @item tol
## Optional tolerance.
## @var{tol} must be a real-valued, non-negative scalar.
## Default value is 0.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool
## True if the system is minimum-phase and false otherwise.
## @end table
##
## @example
## @group
##   real (z) < -tol*(1 + abs (z))    continuous-time
##   abs (z) < 1 - tol                discrete-time
## @end group
## @end example
##
## @strong{References}@*
## [1] Byrnes, C.I. and Isidori, A.
## @cite{A Frequency Domain Philosophy for Nonlinear Systems}.
## IEEE Conf. Dec. Contr. 23, pp. 1569–1573, 1984.
##
## [2] Zeitz, M.
## @cite{Minimum phase – no relevant property of automatic control!}.
## at – Automatisierungstechnik. Volume 62, Issue 1, pp. 3–10, 2014.
##
## [3] Bode, H.W.
## @cite{Network Analysis and Feedback Amplifier Design}.
## D. Van Nostrand Company, pp. 312-318, 1945.
## pp. 341-351, 1992.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2011
## Version: 0.3

function bool = isminimumphase (sys, tol = 0)

  if (nargin > 2)
    print_usage ();
  endif

  z = zero (sys);
  ct = isct (sys);

  bool = __is_stable__ (z, ct, tol);

endfunction


%!assert (isminimumphase (tf (1, [1, 0])), true);
