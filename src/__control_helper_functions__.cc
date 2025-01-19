#include "is_real_scalar.cc"
#include "is_real_vector.cc"
#include "is_real_matrix.cc"
#include "is_real_square_matrix.cc"
#include "is_matrix.cc"
#include "is_vector.cc"
#include "is_zp_vector.cc"
#include "lti_input_idx.cc"


// stub function to avoid gen_doc_cache warning upon package installation
DEFUN_DLD (__control_helper_functions__, args, nargout,
   "-*- texinfo -*-\n"
   "@deftypefn {} __control_helper_functions__ (@dots{})\n"
   "Helper functions for the control package.@*\n"
   "For internal use only.\n"
   "@end deftypefn")
{
    octave_value_list retval;
    error ("__control_helper_functions__: for internal use only");
    return retval;
}

