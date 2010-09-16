#include <octave/oct.h>

DEFUN_DLD (is_real_scalar, args, nargout,
   "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} is_real_scalar (@var{a}, @dots{})\n\
Return true if arguent is a real scalar.\n\
@seealso{isreal, ndims}\n\
@end deftypefn")
{
    octave_value retval = true;
    int nargin = args.length ();

    if (nargin == 0)
    {
        print_usage ();
    }
    else
    {
        for (int i = 0; i < nargin; i++)
        {
            if (args(i).ndims () != 2 || ! args(i).is_scalar_type ()
                || ! args(i).is_numeric_type () || ! args(i).is_real_type ())
            {
                retval = false;
                break;
            }
        }
    }

    return retval;
}
