#include <octave/oct.h>
#include <octave/ov-struct.h>

DEFUN_DLD (nfields2, args, ,
       "-*- texinfo -*-\n\
@deftypefn {Built-in Function} {} nfields (@var{s})\n\
Return the number of fields of the structure @var{s}.\n\
@end deftypefn")
{
  octave_value retval;

  int nargin = args.length ();

  if (nargin == 1 && args(0).is_map ())
    {
      retval = static_cast<double> (args(0).nfields ());
    }
  else
    print_usage ();

  return retval;
}
