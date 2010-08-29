## Copyright (C) 1998, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
##               Kai P. Mueller.
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
## @deftypefn {Function File} {@var{retval} =} isstabilizable (@var{sys}, @var{tol})
## @deftypefnx {Function File} {@var{retval} =} isstabilizable (@var{a}, @var{b}, @var{tol}, @var{dflg})
## Logical check for system stabilizability (i.e., all unstable modes are controllable). 
## Returns 1 if the system is stabilizable, 0 if the system is not stabilizable, -1 
## if the system has non stabilizable modes at the imaginary axis (unit circle for 
## discrete-time systems.
##
## Test for stabilizability is performed via Hautus Lemma. If 
## @iftex
## @tex
## @var{dflg}$\neq$0
## @end tex
## @end iftex
## @ifinfo 
## @var{dflg}!=0
## @end ifinfo
## assume that discrete-time matrices (a,b) are supplied.
## @seealso{size, rows, columns, length, ismatrix, isscalar, isvector, is_observable, is_stabilizable, is_detectable}
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: August 1993
## Updated by A. S. Hodel (scotte@eng.auburn.edu) Aubust, 1995 to use krylovb
## Updated by John Ingram (ingraje@eng.auburn.edu) July, 1996 to accept systems

## FIXME: where has the version which uses krylovb gone?
## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: October 2009
## Version: 0.1

function retval = isstabilizable (a, b = [], tol = [], dflg = 0)

  if (nargin < 1 || nargin > 4)
    print_usage ();
  elseif (isa (a, "lti"))  # system passed
    if (nargin > 2)
      print_usage ();
    endif
    tol = b;
    dflg = isdt (a);
    [a, b] = ssdata (a);
  elseif (nargin == 1)  # a,b arguments sent directly
    print_usage ();
  endif

  if (isempty (tol))
    tol = 200 * eps;
  endif

  ## Checking dimensions
  n = rows (a);
  [nr, m] = size (b);

  if (! issquare (a))
    error ("isstabilizable: a must be square");
  endif
  
  if (nr != n)
    error ("isstabilizable:  (a,b) not conformal");
  endif
  
  ## Computing the eigenvalue of A
  L = eig (a);
  retval = 1;
  specflag = 0;

  for k = 1 : n
    if (dflg)
      ## Discrete time case
      rL = abs (L(k));
      if (rL >= 1)
        H = [eye(n)*L(k)-a, b];
        f = (rank (H, tol) == n);
        if (f == 0)
          retval = 0;
          if (rL == 1)
            specflag = 1;
          endif
        endif
      endif
    else
      ## Continuous time case
      rL = real (L(k));
      if (rL >= 0)
        H = [eye(n)*L(k)-a, b];
        f = (rank (H, tol) == n);
        if (f == 0)
          retval = 0;
          if (rL == 0)
	        specflag = 1;
          endif
        endif
      endif
    endif
  endfor

  if (specflag == 1)
    ## This means that the system has uncontrollable modes at the imaginary axis 
    ## (or at the unit circle for discrete time systems)
    retval = -1;
  endif

endfunction
