## Copyright (C) 2009-2015   Lukas F. Reichlin
## Copyright (C) 2016 Douglas A. Stewart
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
## along with LTI Syncope. If not, see <http://www.gnu.org/licenses/>.


## -*- texinfo -*-
## @deftypefn  {Function File} {[@var{b_out}, @var{a_out}] =} imp_invar (@var{b}, @var{a}, @var{fs}, @var{tol})
## @deftypefnx {Function File} {[@var{b_out}, @var{a_out}] =} imp_invar (@var{b}, @var{a}, @var{fs})
## @deftypefnx {Function File} {[@var{b_out}, @var{a_out}] =} imp_invar (@var{b}, @var{a})
## @deftypefnx {Function File} {[@var{sys_out}] =} imp_invar (@var{b}, @var{a}, @var{fs}. @var{tol})
## @deftypefnx {Function File} {[@var{sys_out}] =} imp_invar (@var{sys_in}, @var{fs}, @var{tol})
## Converts analog filter with coefficients @var{b} and @var{a} and/or @var{sys_in} to digital,
## conserving impulse response.
##
## If @var{fs} is not specified, or is an empty vector, it defaults to 1Hz.
##
## If @var{tol} is not specified, it defaults to 0.0001 (0.1%)
## This function does the inverse of invimpinvar so that the following example should
## restore the original values of @var{a} and @var{b}.
##
## @command{invimpinvar} implements the reverse of this function.
## @example
## [b, a] = impinvar (b, a);
## [b, a] = invimpinvar (b, a);
## @end example
##
## Reference: Thomas J. Cavicchi (1996) ``Impulse invariance and multiple-order
## poles''. IEEE transactions on signal processing, Vol 44 (9): 2344--2347
##
## @seealso{zoh, bilinear, invimpinvar}
## @end deftypefn


function [bz az] = imp_invar (b , a , fs , tol = 1e-4)
  ## This funtion will accept both a
  ## sys variable as input and/or
  ## numerator, denominator as input.

  if (nargin < 1)
    print_usage;
  endif
  
  if (isa (b, "tf") == 1)
  ## the input is an LTI object  
  ## therefore inputs are (sys,fs,tol)
  ## so b is sys
  ## and a is fs
  ## and fs is tol
  
    if (exist("fs","var") != 0)
      tol = fs;
    else
      tol=0.0001;
    endif 
    
    if (exist ("a") == 1)
      fs=a;
    else
      fs=1;
    endif
    
    [b a] = tfdata (b , "v");
  else
     ## the input is vectors
    if (exist ("fs") == 0)
      fs = 1;
    endif
  endif
    
  if (isempty (fs))
    fs = 1;
  endif
   
  if (isempty (tol))
    tol = 1e-4;
  endif
    
  T = 1/fs;
  [r, p, k, e] = residue (b, a); # partial fraction expansion

  if (length (k) > 0) # Greater than zero means we cannot do impulse invariance
    error("Order numerator >= order denominator");
  endif

  
  p = -p;
  p = -exp(-p.*T);
  L = length (r);
  if (L == 1)
    retn = [r(1) 0];
    retd = [1 p(1)];
  else
    rr = nump (r(1), e(1), T, [1 p(1)]);
    [retn retd] = acpfwm (rr, [1 p(1)], r(2), [1 p(2)], e(2), T);
    for k=3:L
      [retn retd] = acpfwm (retn, retd, r(k), [1 p(k)], e(k), T);
    endfor
  endif

  rn = real (retn);
  rn (abs (rn) < 1e-14) = 0; # Get rid of any imaginary parts from round off errors
  rd = real (retd);
  sys11 = tf (rn, rd, T);
  sys2 = minreal (sys11, tol); # Use this to remove the common roots.
  
  if (nargout() < 2)
    bz = sys2;
  else
    [bz, az] = tfdata (sys2, "v");
  endif
  
endfunction


function ret = addvecs2 (v11, v22)
  ## Add 2 vectors when they have different lengths.
  l1 = numel (v11);
  l2 = numel (v22);
  if (l1 != l2)
    if (l1 > l2)
      v22 = prepad (v22,l1);
    else
      v11 = prepad (v11, l2);
    endif
  endif
  ret = v11 + v22;
endfunction


function ret = addvecs1 (v11, v22)
  ml = max (numel (v11), numel (v22));
  ret = prepad1 (v11, ml) .+ prepad1 (v22, ml);
endfunction



function ret = addvecs (v11, v22)
  ## Add 2 vectors when they have different lengths.
  l1 = numel (v11);
  l2 = numel (v22);
  if (l1 != l2)
    if (l1 > l2)
     v22=[zeros(1,l1-l2) v22];
    else
     v11=[zeros(1,l2-l1) v11];
    endif
  endif
  ret = v11 + v22;
endfunction
   


function v1 = convm (v1, v2, m)
  ## If m=0 then return v1
  ## This function mutiplies 2 poynomials together
  ## with the second poly raised to the power M
  ## If v1 == v2 then when m=1 you get v1^2
  ## If v1 != v2 then you get v1*v2^m

  while (m > 0)
    v1 = conv (v1, v2);
    m--;
  endwhile
endfunction




function ret = nump (c1,m,T,p)
  ## create the numerator plynomial
  ## for the Z transform.
  ## c1 is the gain for this term
  ## m is the multiplicity of this term
  ## T is the time period 1/fs
  ## p is the denominator (pole) of this term

  p=-p(2); # get the root value
  switch m  # multiplicity of the term.
    case{1}
      ret = [c1 0];
    case{2}
      ret = [c1*T*p 0];
    case{3};
      ret = c1*([T^2 *p, T^2 *p^2, 0]) / 2;
    case{4}
      ret = c1*[T^3 *p, 4*T^3 *p^2, T^3 *p^3] / 6;
    case(5)
      ret = c1*[T^4 *p, 11*T^4 *p^2,  11*T^4 *p^3, T^4 *p^4] / 24;
    case(6)
      ret = c1*[T^5 *p, 26*T^5 *p^2, 66*T^5 *p^3, 26*T^5 *p^4, T^5 *p^5] / 120;
    case(7)
      ret = c1*[T^6 *p, 57*T^6 *p^2, 302*T^6 *p^3, 302*T^6 *p^4, 57*T^6 *p^5, T^6 *p^6] / factorial(6);
    otherwise
      error ("Too many repeated roots");
  endswitch
endfunction



function [retn retd] = acpfwm (c1,d1,c2,d2,m,T)
  ## add_complex_poly_fracion_with_multiplicity(a,b,m)
  ## Partial fraction expansion creates
  ## a set of polynomials  Cn/(s+Rn)^m
  ## Take 2 of these and add them together.
  ## c1/(s+d1) is assumed to have a multiplicity of 1
  ## and c2/(s+d2)^m has a multiplicity of m

  ## first do (s+d2)^m
  d2m = convm (d2, d2, m-1);
  ## Now calculate the numerator.
  c2m = nump (c2, m, T, d2);
  ## Cross multiply.
  n1 = conv (c1, d2m);
  n2 = conv (c2m, d1);
  ## And put it all together.
  retn = addvecs (n1, n2);
  retd = conv (d1, d2m);
endfunction

