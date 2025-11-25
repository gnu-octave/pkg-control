## Copyright (C) 2009-2015   Lukas F. Reichlin
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
## TF to SS conversion.
## Reference:
## Varga, A.: Computation of irreducible generalized state-space realizations. 
## Kybernetika, 26:89-106, 1990

## Special thanks to Vasile Sima and Andras Varga for their advice.
## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.9

function [retsys, retlti] = __sys2ss__ (sys)

  ## TODO: determine appropriate tolerance from number of inputs
  ##       (since we multiply all denominators in a row), index, ...
  ##       default tolerance from TB01UD is TOLDEF = N*N*EPS 

  ## SECRET WISH: a routine which accepts individual denominators for
  ##              each channel and which supports descriptor systems

  [p, m] = size (sys);
  [num, den] = tfdata (sys);

  len_num = cellfun (@length, num);
  len_den = cellfun (@length, den);

  ## check for properness  
  ## tfpoly ensures that there are no leading zeros
  tmp = len_num > len_den;
  if (any (tmp(:)))       # non-proper transfer function
    ## separation into strictly proper and polynomial part
    [numq, numr] = cellfun (@deconv, num, den, "uniformoutput", false);
    numr = cellfun (@(num_v, len_d) num_v(max(1,end-len_d+2):end), ...
                    numr, num2cell (len_den), "uniformoutput", false);
    numq = cellfun (@__remove_leading_zeros__, numq, "uniformoutput", false);
    numr = cellfun (@__remove_leading_zeros__, numr, "uniformoutput", false);

    ## minimal state-space realization for the proper part
    [a1, b1, c1] = __proper_tf2ss__ (numr, den, p, m);
    e1 = eye (size (a1));

    ## minimal realization for the polynomial part   
    [e2, a2, b2, c2] = __polynomial_tf2ss__ (numq, p, m);

    ## assemble irreducible descriptor realization
    e = blkdiag (e1, e2);
    a = blkdiag (a1, a2);
    b = vertcat (b1, b2);
    c = horzcat (c1, c2);
    retsys = dss (a, b, c, [], e);
  else                    # proper transfer function
    [a, b, c, d] = __proper_tf2ss__ (num, den, p, m);
    retsys = ss (a, b, c, d);
  endif

  retlti = sys.lti;       # preserve lti properties such as tsam

endfunction


## transfer function to state-space conversion for proper models
function [a, b, c, d] = __proper_tf2ss__ (num, den, p, m)
  
  ## allocate variables - check also for zero numerators
  rowmat = zeros (p, m);
  colmat = zeros (p, m);
  zeronum = cellfun (@any, num);

  ## check for equal denominators per rows
  ## estimate maximum ss order
  tot_roword = 0;
  for ik = 1 : p
    id_cnt = 1;
    for il = 1 : m
      if ((rowmat(ik,il) == 0) && zeronum(ik,il))
        eqden = cellfun (@isequal, den(ik,:), {den{ik,il}});
        eqden(!zeronum(ik,:)) = 0;
        if (any (eqden(:)))
          rowmat(ik,:) = rowmat(ik,:) + id_cnt*eqden;
          id_cnt = id_cnt + 1;
          tot_roword = tot_roword + length (den{ik,il}) - 1;
        endif
      endif
    endfor
  endfor

  ## check for equal denominators per columns
  ## estimate maximum ss order
  tot_colord = 0;
  for ik = 1 : m
    id_cnt = 1;
    for il = 1 : p
      if (colmat(il,ik) == 0 && zeronum(il,ik))
        eqden = cellfun (@isequal, den(:,ik), {den{il,ik}});
        eqden(!zeronum(:,ik)) = 0;
        if (any (eqden(:)))
          colmat(:,ik) = colmat(:,ik) + id_cnt*eqden;
          id_cnt = id_cnt + 1;
          tot_colord = tot_colord + length (den{il,ik}) - 1;
        end
      endif
    endfor
  endfor
  
  ## zero system, fast return
  if (tot_roword == 0)
    a = [];
    b = zeros (0,m);
    c = zeros (p,0)
    d = zeros (p,m);
    return;
  endif
  
  ## decide if we proceed row or column-wise:
  ## use the one with the least estimated order
  ## or if the orders are the same, use the one
  ## with the least calls to TD04AD
  nsys_row = max (rowmat(:));
  nsys_col = max (colmat(:));
  if (tot_colord == tot_roword)
    if (nsys_col < nsys_row)
      tot_roword = tot_roword + 1;
    else
      tot_colord = tot_colord + 1;
    endif
  endif
  
  if (tot_colord < tot_roword)
    isrow = false;
    nsys = nsys_col;
    pp = m;
    idxmat = colmat.';
    ord_est = tot_colord;
    den = den.';
  else
    isrow = true;
    nsys = nsys_row;
    pp = p;
    idxmat = rowmat;
    ord_est = tot_roword;
  endif
  
  ## calculate matrices
  a = zeros (ord_est);
  b = zeros (ord_est, m);
  c = zeros (p, ord_est);
  d = zeros (p, m);
  nadd = 1;
  len_num = cellfun (@length, num);
  len_den = cellfun (@length, den);
  
  for ik = 1 : nsys
    
    ## find denominators
    len_den_tmp = ones (pp, 1);
    den_idx = zeros (pp, 1);
    for il = 1 : pp
      fidx = find (idxmat(il,:) == ik, 1);
      if ! (isempty (fidx))
        den_idx(il) = fidx;
        len_den_tmp(il) = len_den(il, fidx);
      endif
    endfor
    
    ## create arrays and fill in the data
    ## in a way that Slicot TD04AD can use
    max_len_den = max (len_den_tmp(:));
    ucoeff = zeros (p, m, max_len_den);
    dcoeff = zeros (pp, max_len_den);
    index = len_den_tmp-1;
    
    for il = 1 : pp
      if (den_idx(il) != 0)
        len = len_den_tmp(il);
        dcoeff(il, 1:len) = den{il,den_idx(il)};
        if (isrow)
          for ij = 1 : m
            if (idxmat(il,ij) == ik)
              ucoeff(il, ij, len-len_num(il,ij)+1 : len) = num{il,ij};
            endif
          endfor
        else
          for ij = 1 : p
            if (idxmat(il,ij) == ik)
              ucoeff(ij, il, len-len_num(ij,il)+1 : len) = num{ij,il};
            endif
          endfor
        endif
      else
        dcoeff(il,1) = 1;
      endif
    endfor
    
    [as, bs, cs, ds] = __sl_td04ad__ (ucoeff, dcoeff, index, 0, isrow);
    
    ## concat arrays
    ncurr = size (as, 1);
    nfin = nadd + ncurr - 1;
    a(nadd:nfin,nadd:nfin) = as;
    b(nadd:nfin,:) = bs;
    c(:,nadd:nfin) = cs;
    d = d + ds;
    
    nadd =  nadd + ncurr;
  endfor
  
  ## return subarrays
  ## (maybe there was zero-pole cancelation)
  a = a(1:nfin,1:nfin);
  b = b(1:nfin,:);
  c = c(:,1:nfin);

endfunction


## realization of the polynomial part according to Andras' paper
function [e2, a2, b2, c2] = __polynomial_tf2ss__ (numq, p, m)

  len_numq = cellfun (@length, numq);
  max_len_numq = max (len_numq(:));
  numq = cellfun (@(x) prepad (x, max_len_numq, 0, 2), numq, "uniformoutput", false);
  f = @(y) cellfun (@(x) x(y), numq);
  s = 1 : max_len_numq;
  D = arrayfun (f, s, "uniformoutput", false);

  e2 = diag (ones (p*(max_len_numq-1), 1), -p);
  a2 = eye (p*max_len_numq);
  b2 = vertcat (D{:});
  c2 = horzcat (zeros (p, p*(max_len_numq-1)), -eye (p));

  ## remove uncontrollable part
  [a2, e2, b2, c2] = __sl_tg01jd__ (a2, e2, b2, c2, 0.0, true, 1, 2);

endfunction


## remove leading zeros from polynomial vector
function p = __remove_leading_zeros__ (p)

  idx = find (p != 0, 1);

  if (isempty (idx))
    p = 0;
  else
    p = p(idx : end);
  endif

endfunction
