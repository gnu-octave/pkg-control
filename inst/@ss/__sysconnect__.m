## Copyright (C) 2009   Lukas F. Reichlin
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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{retsys} =} __sysconnect__ (@var{sys}, @var{M})
## This function is part of the Model Abstraction Layer. No argument checking.
## For internal use only.
## @example
## @group
## Problem: Solve the system equations of
## .
## x(t) = A x(t) + B e(t)
##
## y(t) = C x(t) + D e(t)
##
## e(t) = u(t) + M y(t)
##
## in order to build
## .
## x(t) = F x(t) + G u(t)
##
## y(t) = H x(t) + J u(t)
##
## Solution: Laplace Transformation
## s X(s) = A X(s) + B U(s) + B M Y(s)                       [1]
##
## Y(s) = C X(s) + D U(s) + D M Y(s)                         [2]
##
## solve [2] for Y(s)
## Y(s) = [I - D M]^(-1) C X(s)  +  [I - D M]^(-1) D U(s)
##
## substitute Z = [I - D M]^(-1)
## Y(s) = Z C X(s) + Z D U(s)                                [3]
##
## insert [3] in [1], solve for X(s)
## X(s) = [s I - (A + B M Z C)]^(-1) (B + B M Z D) U(s)      [4]
##
## inserting [4] in [3] finally yields
## Y(s) = Z C [s I - (A + B M Z C)]^(-1) (B + B M Z D) U(s)  +  Z D U(s)
##        \ /        \_____ _____/       \_____ _____/          \ /
##         H               F                   G                 J
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function sys = __sysconnect__ (sys, M)

  [p, m] = size (sys);

  A = sys.a;
  B = sys.b;
  C = sys.c;
  D = sys.d;

  I = eye (p);
  Z = I - D*M;

  if (rcond (Z) < eps)  # check for singularity
    error ("ss: sysconnect: (I - D*M) not invertible because of algebraic loop");
  endif

  Z = inv (Z);

  sys.a = A + B*M*Z*C;  # F
  sys.b = B + B*M*Z*D;  # G
  sys.c = Z*C;          # H
  sys.d = Z*D;          # J

endfunction