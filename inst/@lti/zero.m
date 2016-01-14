## Copyright (C) 2009-2016   Lukas F. Reichlin
## Copyright (C) 2011        Ferdinand Svaricek, UniBw Munich.
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
## @deftypefn {Function File} {@var{z} =} zero (@var{sys})
## @deftypefnx {Function File} {@var{z} =} zero (@var{sys}, @var{type})
## @deftypefnx {Function File} {[@var{z}, @var{k}, @var{info}] =} zero (@var{sys})
## Compute zeros and gain of @acronym{LTI} model.
## By default, @command{zero} computes the invariant zeros,
## also known as Smith zeros.  Alternatively, when called with
## a second input argument, @command{zero} can also compute
## the system zeros, transmission zeros, input decoupling zeros
## and output decoupling zeros.  See paper [1] for an explanation
## of the various zero flavors as well as for further details.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @item type
## String specifying the type of zeros:
## @table @var
## @item 'system', 's'
## Compute the system zeros.
## The system zeros include in all cases
## (square, non-square, degenerate or non-degenerate system) 
## all transmission and decoupling zeros.
## @item 'invariant', 'inv'
## Compute invariant zeros.  Default selection.
## @item 'transmission', 't'
## Compute transmission zeros.  Transmission zeros
## are a subset of the invariant zeros.  
## The transmission zeros are the zeros of the
## Smith-McMillan form of the transfer function matrix.
## @item 'input', 'inp', 'id'
## Compute input decoupling zeros.  The input decoupling zeros are
## also known as the uncontrollable eigenvalues of the pair (A,B).
## @item 'output', 'o', 'od'
## Compute output decoupling zeros.  The output decoupling zeros are
## also known as the unobservable eigenvalues of the pair (A,C).
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item z
## Depending on argument @var{type}, @var{z} contains the
## invariant (default), system, transmission, input decoupling
## or output decoupling zeros of @var{sys} as defined in [1].
## @item k
## Gain of @acronym{SISO} system @var{sys}.  For @acronym{MIMO}
## systems, an empty matrix @code{[]} is returned.
## @item info
## Struct containing additional information.  For details,
## see the documentation of @acronym{SLICOT} routines
## @acronym{AB08ND} and @acronym{AG08BD}.
## @item info.rank
## The normal rank of the transfer function matrix (regular state-space models)
## or of the system pencil (descriptor state-space models).
## @item info.infz
## Contains information on the infinite elementary divisors as follows:
## the system has info.infz(i) infinite elementary divisors of degree i,
## where i=1,2,...,length(info.infz).
## @item info.kronr
## Right Kronecker (column) indices.
## @item info.kronl
## Left Kronecker (row) indices.
## @end table
##
## @strong{Examples}
## @example
## @group
## [z, k, info] = zero (sys)        # invariant zeros
## z = zero (sys, 'system')         # system zeros
## z = zero (sys, 'invariant')      # invariant zeros
## z = zero (sys, 'transmission')   # transmission zeros
## z = zero (sys, 'output')         # output decoupling zeros
## z = zero (sys, 'input')          # input decoupling zeros
## @end group
## @end example
##
## @strong{Algorithm}@*
## For (descriptor) state-space models, @command{zero}
## relies on SLICOT AB08ND and AG08BD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}
## For @acronym{SISO} transfer functions, @command{zero}
## uses Octave's @command{roots}.
## @acronym{MIMO} transfer functions are converted to
## a @emph{minimal} state-space representation for the
## computation of the zeros.
##
## @strong{References}@*
## [1] MacFarlane, A. and Karcanias, N.
## @cite{Poles and zeros of linear multivariable systems:
## a survey of the algebraic, geometric and complex-variable
## theory}.  Int. J. Control, vol. 24, pp. 33-74, 1976.@*
## [2] Rosenbrock, H.H.
## @cite{Correction to 'The zeros of a system'}.
## Int. J. Control, vol. 20, no. 3, pp. 525-527, 1974.@*
## [3] Svaricek, F.
## @cite{Computation of the structural invariants of linear
## multivariable systems with an extended version of the
## program ZEROS}.
## Systems & Control Letters, vol. 6, pp. 261-266, 1985.@*
## [4] Emami-Naeini, A. and Van Dooren, P.
## @cite{Computation of zeros of linear multivariable systems}.
## Automatica, vol. 26, pp. 415-430, 1982.@*
##
## @end deftypefn

## TODO: write a short summary about the characteristics of the
##       various zero flavors and add it to the docstring.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function [zer, gain, info] = zero (sys, type = "invariant")

  if (nargin > 2)
    print_usage ();
  endif

  if (strncmpi (type, "invariant", 3))                              # invariant zeros, default
    [zer, gain, info] = __zero__ (sys, nargout);
  elseif (strncmpi (type, "transmission", 1))                       # transmission zeros
    [zer, gain, info] = zero (minreal (sys));
  elseif (strncmpi (type, "input", 3) || strncmpi (type, "id", 2))  # input decoupling zeros
    [a, b, c, d, e, tsam] = dssdata (sys, []);
    tmp = dss (a, b, zeros (0, columns (a)), zeros (0, columns (b)), e, tsam);
    [zer, gain, info] = zero (tmp);
  elseif (strncmpi (type, "output", 1))                             # output decoupling zeros
    [a, b, c, d, e, tsam] = dssdata (sys, []);
    tmp = dss (a, zeros (rows (a), 0), c, zeros (rows (c), 0), e, tsam);
    [zer, gain, info] = zero (tmp);
  elseif (strncmpi (type, "system", 1))                             # system zeros
    [zer, gain, info] = __szero__ (sys);
  else
    error ("zero: type '%s' invalid", type);
  endif

endfunction


## Function for computing the system zeros.
## Adapted from Ferdinand Svaricek's szero.m
function [z, gain, info] = __szero__ (sys)

  ## TODO: support descriptor state-space models
  ##       with singular 'E' matrices

  [a, b, c, d] = ssdata (sys);
  [pp, mm] = size (sys);
  nn = rows (a);

  ## Tolerance for intersection of zeros
  Zeps = 10 * sqrt ((nn+pp)*(nn+mm)) * eps * norm (a,'fro');

  [z, gain, info] = zero (ss (a, b, c, d));    # zero (sys) lets descriptor test fail
  Rank = info.rank;

  ## System is not degenerated and square
  if (Rank == 0 || (Rank == min(pp,mm) && mm == pp))
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
      [z1, ~, info1] = zero (ss (a, B1, C1, zeros (Rank, Rank)));
      rank1 = info1.rank;
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
                  z2(end+1) = z(ii);
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

endfunction


## Invariant zeros of state-space models
##
## Results from the "Dark Side" 7.5 and 7.8
##
##  -13.2759
##   12.5774
##  -0.0155
##
## Results from Scilab 5.2.0b1 (trzeros)
##
##  - 13.275931  
##    12.577369  
##  - 0.0155265
##
%!shared z, z_exp
%! A = [   -0.7   -0.0458     -12.2        0
%!            0    -0.014   -0.2904   -0.562
%!            1   -0.0057      -1.4        0
%!            1         0         0        0 ];
%!
%! B = [  -19.1      -3.1
%!      -0.0119   -0.0096
%!        -0.14     -0.72
%!            0         0 ];
%!
%! C = [      0         0        -1        1
%!            0         0     0.733        0 ];
%!
%! D = [      0         0
%!       0.0768    0.1134 ];
%!
%! sys = ss (A, B, C, D, "scaled", true);
%! z = sort (zero (sys));
%!
%! z_exp = sort ([-13.2759; 12.5774; -0.0155]);
%!
%!assert (z, z_exp, 1e-4);


## Invariant zeros of regular state-space models
%!shared z, z_exp, info, rank_exp, infz_exp, kronr_exp, kronl_exp
%! A = [  1.0   0.0   0.0   0.0   0.0   0.0
%!        0.0   1.0   0.0   0.0   0.0   0.0
%!        0.0   0.0   3.0   0.0   0.0   0.0
%!        0.0   0.0   0.0  -4.0   0.0   0.0
%!        0.0   0.0   0.0   0.0  -1.0   0.0
%!        0.0   0.0   0.0   0.0   0.0   3.0 ];
%!
%! B = [  0.0  -1.0
%!       -1.0   0.0
%!        1.0  -1.0
%!        0.0   0.0
%!        0.0   1.0
%!       -1.0  -1.0 ];
%!
%! C = [  1.0   0.0   0.0   1.0   0.0   0.0
%!        0.0   1.0   0.0   1.0   0.0   1.0
%!        0.0   0.0   1.0   0.0   0.0   1.0 ];
%!
%! D = [  0.0   0.0
%!        0.0   0.0
%!        0.0   0.0 ];
%!
%! sys = ss (A, B, C, D, "scaled", true);
%! [z, ~, info] = zero (sys);
%!
%! z_exp = [ 2.0000
%!          -1.0000 ];
%!
%! rank_exp = 2;
%! infz_exp = 2;
%! kronr_exp = zeros (1, 0);
%! kronl_exp = 2;
%!
%!assert (z, z_exp, 1e-4);
%!assert (info.rank, rank_exp);
%!assert (info.infz, infz_exp);
%!assert (info.kronr, kronr_exp);
%!assert (info.kronl, kronl_exp);


## Invariant zeros of descriptor state-space models
%!shared z, z_exp, info, rank_exp, infz_exp, kronr_exp, kronl_exp
%! A = [  1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     1     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     1     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0
%!        0     0     0     0     0     0     0     0     1 ];
%!
%! E = [  0     0     0     0     0     0     0     0     0
%!        1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0 ];
%!
%! B = [ -1     0     0
%!        0     0     0
%!        0     0     0
%!        0    -1     0
%!        0     0     0
%!        0     0     0
%!        0     0    -1
%!        0     0     0
%!        0     0     0 ];
%!
%! C = [  0     1     1     0     3     4     0     0     2
%!        0     1     0     0     4     0     0     2     0
%!        0     0     1     0    -1     4     0    -2     2 ];
%!
%! D = [  1     2    -2
%!        0    -1    -2
%!        0     0     0 ];
%!
%! sys = dss (A, B, C, D, E, "scaled", true);
%! [z, ~, info] = zero (sys);
%!
%! z_exp = 1;
%!
%! rank_exp = 11;
%! infz_exp = [0, 1];
%! kronr_exp = 2;
%! kronl_exp = 1;
%!
%!assert (z, z_exp, 1e-4);
%!assert (info.rank, rank_exp);
%!assert (info.infz, infz_exp);
%!assert (info.kronr, kronr_exp);
%!assert (info.kronl, kronl_exp);


## Gain of descriptor state-space models
%!shared p, pi, z, zi, k, ki, p_tf, pi_tf, z_tf, zi_tf, k_tf, ki_tf
%! P = ss (-2, 3, 4, 5);
%! Pi = inv (P);
%!
%! p = pole (P);
%! [z, k] = zero (P);
%!
%! pi = pole (Pi);
%! [zi, ki] = zero (Pi);
%!
%! P_tf = tf (P);
%! Pi_tf = tf (Pi);
%!
%! p_tf = pole (P_tf);
%! [z_tf, k_tf] = zero (P_tf);
%!
%! pi_tf = pole (Pi_tf);
%! [zi_tf, ki_tf] = zero (Pi_tf);
%!
%!assert (p, zi, 1e-4);
%!assert (z, pi, 1e-4);
%!assert (k, inv (ki), 1e-4);
%!assert (p_tf, zi_tf, 1e-4);
%!assert (z_tf, pi_tf, 1e-4);
%!assert (k_tf, inv (ki_tf), 1e-4);


## Example taken from Paper [1]
## Regular state-space system
%!shared z_inv, z_tra, z_inp, z_out, z_sys, z_inv_e, z_tra_e, z_inp_e, z_out_e, z_sys_e
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
%!
%! z_inv = zero (SYS);
%! z_tra = zero (SYS, "transmission");
%! z_inp = zero (SYS, "input decoupling");
%! z_out = zero (SYS, "output decoupling");
%! z_sys = zero (SYS, "system");
%!
%! z_inv_e = [2; -1];
%! z_tra_e = [2];
%! z_inp_e = [-4];
%! z_out_e = [-1];
%! z_sys_e = [-4, -1, 2];
%! 
%!assert (z_inv, z_inv_e, 1e-4); 
%!assert (z_tra, z_tra_e, 1e-4); 
%!assert (z_inp, z_inp_e, 1e-4); 
%!assert (z_out, z_out_e, 1e-4); 
%!assert (z_sys, z_sys_e, 1e-4);


## Example taken from Paper [1]
## Well, this is not exactly a descriptor state-space model,
## but it is the best thing I have right now and it is better
## than no test at all.  The routine for the system zeros works
## only for descriptor state-space models with regular 'E' matrices.
%!shared z_inv, z_tra, z_inp, z_out, z_sys, z_inv_e, z_tra_e, z_inp_e, z_out_e, z_sys_e
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
%! E = eye (6);
%! 
%! SYS = dss (A, B, C, D, E);
%!
%! z_inv = zero (SYS);
%! z_tra = zero (SYS, "transmission");
%! z_inp = zero (SYS, "input decoupling");
%! z_out = zero (SYS, "output decoupling");
%! z_sys = zero (SYS, "system");
%!
%! z_inv_e = [2; -1];
%! z_tra_e = [2];
%! z_inp_e = [-4];
%! z_out_e = [-1];
%! z_sys_e = [-4, -1, 2];
%! 
%!assert (z_inv, z_inv_e, 1e-4); 
%!assert (z_tra, z_tra_e, 1e-4); 
%!assert (z_inp, z_inp_e, 1e-4); 
%!assert (z_out, z_out_e, 1e-4); 
%!assert (z_sys, z_sys_e, 1e-4);
