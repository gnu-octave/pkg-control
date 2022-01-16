## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## Subscripted reference for @acronym{LTI} objects.
## Used by Octave for "sys = sys(2:4, :)" or "val = sys.prop".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.5

function a = subsref (a, s)

  if (numel (s) == 0)
    return;
  endif

  switch (s(1).type)
    case "()"
      idx = s(1).subs;
      switch (numel (idx))
        case 2
          a = __sys_prune__ (a, idx{1}, idx{2});
        case 1
          a = freqresp (a, idx{1});
        case 0
          ## NOTE: for example, Octave returns  a = a()  for  a = rand (2, 2)
          ##       so we do the same for LTI models.
          ## a = a;  # superfluous
          ## error (["lti: subsref: invalid subscripted reference of type '()'.  ", ...
          ##         "need one or two arguments inside the brackets."]);
        otherwise
          ## bug #45314  <https://savannah.gnu.org/bugs/?45314>
          error (["lti: subsref: arrays of linear models are not supported (yet).  ", ...
                  "as a workaround, try to use cells of LTI models instead."]);
      endswitch
    case "."
      fld = s(1).subs;
      a = get (a, fld);
      ## warning ("lti: subsref: do not use subsref for development");
    otherwise
      error ("lti: subsref: invalid subscript type '%s'", s(1).type);
  endswitch
  
  a = subsref (a, s(2:end));

endfunction


## lti: subsref
%!shared a
%! s = tf ("s");
%! G = (s+1)*s*5/(s+1)/(s^2+s+1);
%! a = G(1,1).num{1,1}(1);
%!assert (a, 5, 1e-4);
