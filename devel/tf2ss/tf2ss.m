function retsys = tf2ss ()

num = {[1, 5, 7], [1]; [1, 7], [1, 5, 5]};
den = {[1, 5, 6], [1, 2]; [1, 8, 6], [1, 3, 2]};
sys = tf (num, den);

%{
sys = tf (1, [1, 0])
sys = tf (1, [1, 1])

sys = tf (1, conv ([1, 1, 1], [1, 4, 6, 4, 1]))

sys = tf ()
sys = tf ("s")
%}

sys = tf (WestlandLynx);
tol = sqrt (eps)

  [p, m] = size (sys);
  [num, den] = tfdata (sys);
  
  numc = cell (p, m);
  denc = cell (p, 1);
  
  ## multiply all denominators in a row and
  ## update each numerator accordingly 
  for i = 1 : p
    denc(i) = __conv__ (den{i,:});
    for j = 1 : m
      idx = setdiff (1:m, j);
      numc(i,j) = __conv__ (num{i,j}, den{i,idx});
    endfor
  endfor

  len_numc = cellfun (@length, numc);
  len_denc = cellfun (@length, denc);
  
  ## tfpoly ensures that there are no leading zeros
  tmp = len_numc > repmat (len_denc, 1, m);
  if (any (tmp(:)))
    error ("tf: tf2ss: system must be proper");
  endif

  max_len_denc = max (len_denc(:));
  ucoeff = zeros (p, m, max_len_denc);
  dcoeff = zeros (p, max_len_denc);
  index = len_denc-1;

  for i = 1 : p
    len = len_denc(i);
    dcoeff(i, 1:len) = denc{i};
    for j = 1 : m
      ucoeff(i, j, len-len_numc(i,j)+1 : len) = numc{i,j};
    endfor
  endfor
index, prod (index), eps*prod (index)
  [a, b, c, d] = sltd04ad (ucoeff, dcoeff, index, tol);

  retsys = ss (a, b, c, d);

endfunction


function vec = __conv__ (vec, varargin)

  if (nargin == 1)
    return;
  else
    for k = 1 : nargin-1
      vec = conv (vec, varargin{k});
    endfor
  endif

endfunction
