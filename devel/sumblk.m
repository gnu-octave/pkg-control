function s = sumblk (formula, n = 1)

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif
  
  if (! ischar (formula))
    error ("sumblk: require string as first argument");
  endif
  
  if (! is_real_scalar (n))
    error ("sumblk: require integer as second argument");
  endif
  
  if (length (strfind (formula, "=")) != 1)
    error ("sumblk: formula requires exactly one '='");
  endif

  ## if the first signal has no sign, add a plus
  if (isempty (regexp (formula, "=\\s*[+-]", "once")))
    formula = regexprep (formula, "=", "=+", "once");
  endif
  
  ## remove "="
  token = regexp (formula, "[=+-]", "match");
  token = token(2:end);
  formula = formula(formula != "=");
  
  ## extract signal names
  signal = regexp (formula, "\\s*[+-]\\s*", "split");
  if (any (cellfun (@isempty, signal)))
    error ("sumblk: formula is missing some input/output names");
  endif
  outname = signal(1);
  inname = signal(2:end);
  signs = ones (1, numel (inname));
  signs(strcmp (token, "-")) = -1;
  
  if (n != 1)
    outname = strseq (outname{1}, 1:n);
    tmp = cellfun (@strseq, inname, {1:n}, "uniformoutput", false);
    inname = vertcat (tmp{:});
  endif
  
  d = kron (signs, eye (n));
  s = ss (d);
  s = set (s, "inname", inname, "outname", outname);

endfunction
