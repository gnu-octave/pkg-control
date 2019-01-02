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

Return true if all arguments zero-pole vectors and false otherwise.

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: October 2015
Version: 0.1

*/

#include <octave/oct.h>

// PKG_ADD: autoload ("is_zp_vector", "__control_helper_functions__.oct");    
DEFUN_DLD (is_zp_vector, args, nargout,
   "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} is_zp_vector (@var{a}, @dots{})\n\
Return true if all arguments are zero-pole vectors and false otherwise.\n\
@var{[]} is a valid zero-pole vector.\n\
Avoid nasty stuff like @code{true = isreal (\"a\")}\n\
@seealso{is_real_matrix, is_real_square_matrix, is_real_vector, is_real_scalar}\n\
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
            if (args(i).ndims () != 2
                || (args(i).rows () > 1 && args(i).columns () > 1)
                || ! args(i).isnumeric ()
                || ! (args(i).iscomplex () || args(i).isreal ()))
            {
                retval = false;
                break;
            }
        }
    }

    return retval;
}
