## FIXME: avoid duplicated code from __lti_group__

function retlti = __lti_vertcat__ (lti1, lti2)

  retlti = lti ();
  
  retlti.inname = lti1.inname                   # FIXME: return empty inname
  retlti.outname = [lti1.outname; lti2.outname];
  
  if (nfields (lti1.outgroup) || nfields (lti2.outgroup))
    p1 = numel (lti1.outname);
    lti2_outgroup = structfun (@(x) x + p1, lti2.outgroup, "uniformoutput", false);
    retlti.outgroup = __merge_struct__ (lti1.outgroup, lti2_outgroup, "out");
  endif

  if (lti1.tsam == lti2.tsam)
    retlti.tsam = lti1.tsam;
  elseif (lti1.tsam == -2)
    retlti.tsam = lti2.tsam;
  elseif (lti2.tsam == -2)
    retlti.tsam = lti1.tsam;
  elseif (lti1.tsam == -1 && lti2.tsam > 0)
    retlti.tsam = lti2.tsam;
  elseif (lti2.tsam == -1 && lti1.tsam > 0)
    retlti.tsam = lti1.tsam;
  else
    error ("lti_group: systems must have identical sampling times");
  endif

endfunction


function ret = __merge_struct__ (a, b, iostr)

  ## FIXME: this is too complicated;
  ##        isn't there a simple function for this task?

  ## bug #40224: orderfields (struct ()) errors out in Octave 3.6.4
  ## therefore use nfields to check for empty structs
  if (nfields (a))
    a = orderfields (a);
  endif
  if (nfields (b))
    b = orderfields (b);
  endif

  fa = fieldnames (a);
  fb = fieldnames (b);
  [fi, ia, ib] = intersect (fa, fb);
  ca = struct2cell (a);
  cb = struct2cell (b);

  for k = numel (fi) : -1 : 1
    ca{ia(k)} = vertcat (ca{ia(k)}(:), cb{ib(k)}(:));
    fb(ib(k)) = [];
    cb(ib(k)) = [];
  endfor
  
  ret = cell2struct ([ca; cb], [fa; fb]);

endfunction
