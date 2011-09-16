function retsys = tf2ss ()

num = {[1, 5, 7], [1]; [1, 7], [1, 5, 5]};
den = {[1, 5, 6], [1, 2]; [1, 8, 6], [1, 3, 2]};
sys = tf (num, den)

% __conv__ (num{:})

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
  
  max_len_numc = max (len_numc(:));
  max_len_denc = max (len_denc(:));

  ## tfpoly ensures that there are no leading zeros
  if (length (max_len_numc) > length (max_len_denc))
    error ("tf: tf2ss: system must be proper");
  endif

  ucoeff = zeros (p, m, max_len_denc);
  dcoeff = zeros (p, max_len_denc);
  index = (max_len_denc-1) * ones (p, 1);

  for i = 1 : p
    for j = 1 : m
      len = len_numc(i,j);
      ucoeff(i, j, max_len_denc-len+1 : max_len_denc) = numc{i,j};
    endfor
    len = len_denc(i);
    dcoeff(i, max_len_denc-len+1 : max_len_denc) = denc{i};
  endfor
ucoeff, dcoeff, index
  [a, b, c, d] = sltd04ad (ucoeff, dcoeff, index);

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
