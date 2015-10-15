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

  if (nargin == 1 && args(0).is_defined () && args(0).is_cell ())
  {
    octave_idx_type len = args(0).cell_value().nelem();
    octave_idx_type idx = len;
    
    for (octave_idx_type i = 0; i < len; i++)
    {
      if (args(0).cell_value().elem(i).is_string ())
      {
        idx = i;
        break;
      }
    }
    
    Range mat_idx (1, idx);
    Range opt_idx (idx+1, len);

    retval(1) = opt_idx;
    retval(0) = mat_idx;
  }
  else
  {
    print_usage ();
  }

  return retval;
}
