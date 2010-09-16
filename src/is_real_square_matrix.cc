#include <octave/oct.h>

DEFUN_DLD (is_real_square_matrix, args, nargout,
   "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} is_real_square_matrix (@var{a}, @dots{})\n\
Return true if arguments look like a square state-space matrix.\n\
Avoid nasty stuff like @code{true = isreal (""a"")}\n\
@seealso{is_real_matrix, is_real_vector, is_real_scalar}\n\
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
            if (args(i).ndims () != 2 || args(i).rows () != args(i).columns ()
                || ! args(i).is_numeric_type () || ! args(i).is_real_type ())
            {
                retval = false;
                break;
            }
        }
    }

    return retval;
}
