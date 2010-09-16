#include <octave/oct.h>

DEFUN_DLD (is_ss_mat, args, nargout,
   "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} is_ss_mat (@var{a}, @dots{})\n\
Return true if arguents look like a state-space matrix.\n\
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
            if (args(i).ndims () != 2 || ! args(i).is_numeric_type ()
                || ! args(i).is_real_type ())
            {
                retval = false;
                break;
            }
        }
    }

    return retval;
}
