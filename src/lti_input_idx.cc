/*

Copyright (C) 2009-2016   Lukas F. Reichlin

This file is part of LTI Syncope.

LTI Syncope is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LTI Syncope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

First, this function implemented the following Octave function in C++

function [mat_idx, opt_idx] = __lti_input_idx__ (varargin)

  str_idx = find (cellfun (@ischar, varargin));

  if (isempty (str_idx))
    mat_idx = 1 : nargin;
    opt_idx = [];
  else
    mat_idx = 1 : str_idx(1)-1;
    opt_idx = str_idx(1) : nargin;
  endif

endfunction

Later on, the C++ function was extended such that
it recognizes classes in some cases.  See comment
block in the code for details.  I know this looks
like a horrible definition for a function, but it
is exactly the behavior I need.

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2015
Version: 0.1

*/


#include <octave/oct.h>

// PKG_ADD: autoload ("__lti_input_idx__", "__control_helper_functions__.oct");    
DEFUN_DLD (__lti_input_idx__, args, ,
       "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{mat_idx}, @var{opt_idx}, @var{obj_flg}] =} __lti_input_idx__ (@var{args})\n\
Return some indices for cell @var{args}.  For internal use only.\n\
Read the source code in @code{lti_input_idx.cc} for details.\n\
@end deftypefn")
{
  octave_value_list retval;
  octave_idx_type nargin = args.length ();

  // first, check whether a cell is passed
  if (nargin == 1 && args(0).is_defined () && args(0).iscell ())
  {
    octave_idx_type len = args(0).cell_value().numel();
    octave_idx_type idx = len;
    octave_idx_type offset = 0;

    // if the cell is not empty, look for the first string
    for (octave_idx_type i = 0; i < len; i++)
    {
      if (args(0).cell_value().elem(i).is_string ())
      {
        idx = i;
        break;
      }
    }

    // ** If the element before the first string is an object,
    //    then this object belongs to the option index.
    //      ss (d, ltisys, 'key', val)
    // ** If there is no string at all in cell args(0), check
    //    whether the last element in args(0) is an object.
    //    If so, this object also belongs to the option index.
    //      tf (num, den, ltisys)
    // ** All other objects (before built-in types (except chars)
    //    and after strings) are not recognized as objects.
    //      ss (a, b, ltisys, c, d, 'key', val, 'lti', ltisys)
    if (len > 0 && idx > 0
        && args(0).cell_value().elem(idx-1).isobject ())
    {
      offset = 1;
    }

    Range mat_idx (1, idx-offset);
    Range opt_idx (idx+1-offset, len);

    retval(2) = offset;     // abused as logical in the LTI constructors
    retval(1) = opt_idx;
    retval(0) = mat_idx;
  }
  else
  {
    print_usage ();
  }

  return retval;
}
