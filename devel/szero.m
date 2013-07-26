## Copyright (C) 2011   Ferdinand Svaricek, UniBw Munich.
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
## @deftypefn {Function File} {@var{z} =} szero (@var{sys})
## @deftypefnx {Function File} {[@var{z}, @var{rnk}] =} szero (@var{sys})
## Compute @emph{system} zeros of @acronym{LTI} model.
## According to [1], system zeros are defined as the set
## of zeros which includes both the transmission and
## decoupling zero sets for a given system.
## In case you are not sure which zeros you need and you're
## just looking for something like the @emph{vanilla} zeros,
## use function @command{zero} instead.
##
## @strong{Inputs}
## @table @var
## @item sys
## State-space model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item z
## System zeros of @var{sys} as defined in [1] and [2].
## The system zeros include in all cases
## (square, non-square, degenerate or non-degenerate system) 
## all transmission and decoupling zeros.
## @item rnk
## The normal rank of the transfer function matrix.
## @end table
##
## @strong{References}@*
## [1] MacFarlane, A. and Karcanias, N.
## @cite{Poles and zeros of linear multivariable systems:
## a survey of the algebraic, geometric and complex-variable
## theory}.  Int. J. Control, vol. 24, pp. 33-74, 1976.@*
## [2] Rosenbrock, H.H.
## @cite{Correction to 'The zeros of a system'}.
## Int. J. Control, vol. 20, no. 3, pp. 525-527, 1974.@*
##
## @seealso{zero, tzero}
## @end deftypefn

## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: July 2013
## Version: 0.1

function [z, Rank] = szero (sys)

  if (nargin != 1)
    print_usage ();
  endif

  [a, b, c, d] = ssdata (sys);
  [pp, mm] = size (sys);
  nn = rows (a);

  ## Tolerance for intersection of zeros
  Zeps = 10 * sqrt ((nn+pp)*(nn+mm)) * eps * norm (a,'fro');

  [z, ~, Rank] = zero (sys);

  ## System is not degenerated and square
  if (Rank == 0 || (Rank == min(pp,mm) && mm == pp))
    z = sort (z(:));
    return;
  endif

  ## System (A,B,C,D) is degenerated and/or non-square
  z = [];

  ## Computation of the greatest common divisor of all minors of the 
  ## Rosenbrock system matrix that have the following form
  ##
  ##    1, 2, ..., n, n+i_1, n+i_2, ..., n+i_k
  ##   P
  ##    1, 2, ..., n, n+j_1, n+j_2, ..., n+j_k
  ## 
  ## with k = Rank.

  NKP = nchoosek (1:pp, Rank);
  [IP, JP] = size (NKP);
  NKM = nchoosek (1:mm, Rank);
  [IM, JM] = size (NKM);

  for i = 1:IP
    for j = 1:JP
      k = NKP(i,j);
      C1(j,:) = c(k,:);         # Build C of dimension (Rank x n)
    endfor
    for ii = 1:IM
      for jj = 1:JM
        k = NKM(ii,jj);
        B1(:,jj) = b(:,k);      # Build B of dimension (n x Rank)
      endfor
      [z1, ~, rank1] = zero (ss (a, B1, C1, zeros (Rank, Rank)));
      if (rank1 == Rank)
        if (isempty (z1))
          z = z1;               # Subsystem has no zeros -> system has no system zeros
          return;
        else
          if (isempty (z))
            z = z1;             # Zeros of the first subsystem
          else                  # Compute intersection of z and z1 with tolerance Zeps 
            z2 = [];
            for ii=1:length(z)
              for jj=1:length(z1)
                if (abs (z(ii)-z1(jj)) < Zeps)
                  z2(length(z2)+1) = z(ii);
                  z1(jj) = [];
                  break;
                endif
              endfor
            endfor
            z = z2;             # System zeros are the common zeros of all subsystems
          endif
        endif
      endif
    endfor
  endfor

  z = sort (z(:));

endfunction


## Example taken from Paper [1]
%!shared zo, ze
%! A = diag ([1, 1, 3, -4, -1, 3]);
%! 
%! B = [  0,  -1
%!       -1,   0
%!        1,  -1
%!        0,   0
%!        0,   1
%!       -1,  -1  ];
%!        
%! C = [  1,  0,  0,  1,  0,  0
%!        0,  1,  0,  1,  0,  1
%!        0,  0,  1,  0,  0,  1  ];
%!         
%! D = zeros (3, 2);
%! 
%! SYS = ss (A, B, C, D);
%! zo = szero (SYS);
%! 
%! ze = [-4; -1; 2];
%! 
%!assert (zo, ze, 1e-4); 
