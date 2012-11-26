function db = mag2db (mag)

  db = 20 .* log10 (mag);
  db(mag < 0) = NaN;

endfunction
