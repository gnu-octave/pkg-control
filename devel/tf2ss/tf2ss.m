num = {[1, 5, 7], [1]; [1, 7], [1, 5, 5]};
den = {[1, 5, 6], [1, 2]; [1, 8, 6], [1, 3, 2]};
sys = tf (num, den)

  [p, m] = size (sys);
  [num, den] = tfdata (sys);

  ## row denominators

  len_num = cellfun (@length, num);
  len_den = cellfun (@length, den);
  
  max_len_num = max (len_num(:));
  max_len_den = max (len_den(:));

  ## tfpoly ensures that there are no leading zeros
  if (length (max_len_num) > length (max_len_den))
    error ("tf: tf2ss: system must be proper");
  endif

  ucoeff = zeros (p, m, max_len_den);
  dcoeff = zeros (p, max_len_den);

  for i = 1 : p
    for j = 1 : m
      len = len_num(i,j);
      ucoeff(i, j, max_len_den-len+1 : max_len_den) = num{i,j};
    endfor
    len = len_row_den(i);
    dcoeff(i, max_len_den-len+1 : max_len_den) = row_den{i};
  endfor

  [a, b, c, d] = sltd04ad (ucoeff, dcoeff, index);

  retsys = ss (a, b, c, d);

