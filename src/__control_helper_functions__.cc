#include "is_real_scalar.cc"
#include "is_real_vector.cc"
#include "is_real_matrix.cc"
#include "is_real_square_matrix.cc"
#include "is_matrix.cc"
#include "is_zp_vector.cc"
#include "lti_input_idx.cc"

// #include "nfields2.cc"  // delete this if support for Octave 3.8 gets dropped



// stub function to avoid gen_doc_cache warning upon package installation
DEFUN_DLD (__control_helper_functions__, args, nargout,
   "-*- texinfo -*-\n\
Helper functions for the control package.\n\
For internal use only.")
{
    octave_value_list retval;
    error ("__control_helper_functions__: for internal use only");
    return retval;
}

