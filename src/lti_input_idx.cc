/*
Implements the following Octave function in C++

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
*/


#include <octave/oct.h>

// PKG_ADD: autoload ("__lti_input_idx__", "__control_helper_functions__.oct");    
DEFUN_DLD (__lti_input_idx__, args, ,
       "-*- texinfo -*-\n\
@deftypefn {Built-in Function} {} ischar (@var{x})\n\
Return true if @var{x} is a character array.\n\
@seealso{isfloat, isinteger, islogical, isnumeric, iscellstr, isa}\n\
@end deftypefn")
{
  octave_value_list retval;
  octave_idx_type nargin = args.length ();

  // first, check whether a cell is passed
  if (nargin == 1 && args(0).is_defined () && args(0).is_cell ())
  {
    octave_idx_type len = args(0).cell_value().nelem();
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
    // ** If there is no string at all in cell args(0), check
    //    whether the last element in args(0) is an object.
    //    If so, this object also belongs to the option index.
    // ** All other objects (before built-in types (except chars)
    //    and after strings) are not recognized as objects.
    if (len > 0 && idx > 0
        && args(0).cell_value().elem(idx-1).is_object ())
    {
      offset = 1;
    }

    Range mat_idx (1, idx-offset);
    Range opt_idx (idx+1-offset, len);

    retval(2) = offset;
    retval(1) = opt_idx;
    retval(0) = mat_idx;
  }
  else
  {
    print_usage ();
  }

  return retval;
}
