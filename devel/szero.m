%SZERO  System zeros of LTI systems.
%% 
%%    Z = SZERO(A,B,C,D) returns the system zeros of the LTI 
%%    system (A,B,C,D). The system zeros include in all cases
%%    (square, non-square, degenerate or non-degenerate system) 
%%    all transmission and decoupling zeros. 
%%
%%    [Z,RANK] = SZERO(A,B,C,D) also returns the normal rank of 
%%    the transfer function matrix.
%% 
%%
%% References: H.H. Rosenbrock:
%%             Correction to 'The zeros of a system',
%%             INT.J. CONTROL, 1974, VOL. 20, No.3, 525-527.
%%
%%    See also TZERO.
%%
%%   Copyright 2011 Ferdinand Svaricek, UniBw Munich.
%%
%%   $Revision: 1.1 $  $Date: 2013/07/24 14:40:00 $
%%

function [z, Rank] = szero (sys)

  [a, b, c, d] = ssdata (sys);

  nn = size(a,1);
  pp = size(c,1);
  mm = size(b,2);

  % Tolerance for intersection of zeros
  Zeps = 10*sqrt((nn+pp)*(nn+mm))*eps*norm(a,'fro');

  [z, ~, Rank] = zero (ss (a,b,c,d));
  z = sort (z(:));

  if Rank == 0
    return;
  endif

  if Rank==min(pp,mm) & mm==pp
  % System (A,B,C,D) is not degenerated and square
    return
  else % System (A,B,C,D) is degenerated and/or non-square
    z = [];
    %
    % Computation of the greatest common divisor of all minors of the 
    % Rosenbrock system matrix that have the following form
    %
    %    1, 2, ..., n, n+i_1, n+i_2, ..., n+i_k
    %   P
    %    1, 2, ..., n, n+j_1, n+j_2, ..., n+j_k
    % 
    % with k = Rank.
    %
    NKP = nchoosek(1:pp,Rank);
    [IP,JP] = size(NKP);
    NKM = nchoosek(1:mm,Rank);
    [IM,JM] = size(NKM);
    for i = 1:IP
      for j = 1:JP
        k = NKP(i,j);
        C1(j,:) = c(k,:); % Build C of dimension (Rank x n)
      endfor
      for ii = 1:IM
        for jj = 1:JM
          k = NKM(ii,jj);
          B1(:,jj) = b(:,k); % Build B of dimension (n x Rank)
        endfor
        [z1,~,rank1] = zero (ss (a, B1, C1, zeros (Rank, Rank)));
        if rank1 == Rank
          if isempty(z1)
            z = sort (z1(:));    % Subsystem has no zeros -> system has no system zeros
            return;
          else
            if isempty(z)
              z = z1;  % Zeros of the first subsystem
            else    % Compute intersection of z and z1 with tolerance Zeps 
              z2=[];
              for ii=1:length(z)
                for jj=1:length(z1)
                  if abs(z(ii)-z1(jj)) < Zeps
                    z2(length(z2)+1)=z(ii);
                    z1(jj)=[];
                    break
                  endif
                endfor
              endfor
              z = z2; % System zeros are the common zeros of all subsystems
            endif
          endif
        endif
      endfor
    endfor
  endif
  
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
