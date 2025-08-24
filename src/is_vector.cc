/*

Copyright (C) 2024  Torsten Lilge

This file is part of the GNU Octave control package

GNU Octave is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This Software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GNU Octave control pakcage.
If not, see <http://www.gnu.org/licenses/>.

Return true if all arguments are real-valued vectors and false otherwise.

Author: Torsten Lilge <ttl-octave@mailbox.org>
Created: November 2024

*/

#include <octave/oct.h>
#include "config.h"

// PKG_ADD: autoload ("is_vector", "__control_helper_functions__.oct");    
DEFUN_DLD (is_vector, args, nargout,
   "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} is_vector (@var{a}, @dots{})\n\
Return true if all arguments are vectors and false otherwise.\n\
@var{[]} is not a valid vector.\n\
Avoid nasty stuff like @code{true = isreal (\"a\")}\n\
@seealso{is_matrix, is_real_vector}\n\
@end deftypefn")
{
    octave_value retval = true;
    octave_idx_type nargin = args.length ();

    if (nargin == 0)
    {
        print_usage ();
    }
    else
    {
        for (octave_idx_type i = 0; i < nargin; i++)
        {
            if (args(i).ndims () != 2 || ! (args(i).rows () == 1 || args(i).columns () == 1)
                || ! args(i).isnumeric ())
            {
                retval = false;
                break;
            }
        }
    }

    return retval;
}
