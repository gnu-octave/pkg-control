function key = __match_key__ (str, props, caller = "match_key")

  if (! ischar (str))
    error ("%s: property name must be a string", caller);
  endif

  ## exact matching - needed for e.g. iddata properties "u" and "userdata"
  idx = strcmpi (str, props);
  n = sum (idx);
  
  if (n == 1)           # 1 exact match
    key = lower (str);
    return;
  elseif (n > 1)        # props are not unique, this would be a bug in the control package
    error ("%s: property name '%s' is ambiguous", caller, str);
  endif
  
  ## partial matching - n was zero
  idx = strncmpi (str, props, length (str));
  n = sum (idx);
  
  if (n == 1)
    key = lower (props{idx});
    return;
  elseif (n > 1)
    error ("%s: property name '%s' is ambiguous", caller, str);
  endif
  
  error ("%s: property name '%s' is unknown", caller, str);
  
endfunction
