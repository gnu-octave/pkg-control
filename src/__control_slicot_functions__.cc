#include "sl_ab08nd.cc"  // transmission zeros of state-space models
#include "sl_ab13dd.cc"  // L-infinity norm
#include "sl_sb10hd.cc"  // H-2 controller synthesis - continuous-time
#include "sl_sb10ed.cc"  // H-2 controller synthesis - discrete-time
#include "sl_ab13bd.cc"  // H-2 norm
#include "sl_sb01bd.cc"  // Pole assignment
#include "sl_sb10fd.cc"  // H-infinity controller synthesis - continuous-time
#include "sl_sb10dd.cc"  // H-infinity controller synthesis - discrete-time
#include "sl_sb03md.cc"  // Lyapunov equations
#include "sl_sb04md.cc"  // Sylvester equations - continuous-time
#include "sl_sb04qd.cc"  // Sylvester equations - discrete-time
#include "sl_sg03ad.cc"  // generalized Lyapunov equations
#include "sl_sb02od.cc"  // algebraic Riccati equations
#include "sl_ab13ad.cc"  // Hankel singular values
#include "sl_ab01od.cc"  // staircase form using orthogonal transformations
#include "sl_tb01pd.cc"  // minimal realization of state-space models
#include "sl_sb03od.cc"  // Cholesky factor of Lyapunov equations
#include "sl_sg03bd.cc"  // Cholesky factor of generalized Lyapunov equations
#include "sl_ag08bd.cc"  // transmission zeros of descriptor state-space models
#include "sl_tg01jd.cc"  // minimal realization of descriptor state-space models
#include "sl_tg01hd.cc"  // controllability staircase form of descriptor state-space models
#include "sl_tg01id.cc"  // observability staircase form of descriptor state-space models
#include "sl_sg02ad.cc"  // solution of algebraic Riccati equations for descriptor systems
#include "sl_tg04bx.cc"  // gain of descriptor state-space models
#include "sl_tb01id.cc"  // scaling of state-space models
#include "sl_tg01ad.cc"  // scaling of descriptor state-space models
#include "sl_sb10id.cc"  // H-infinity loop shaping - continuous-time
#include "sl_sb10kd.cc"  // H-infinity loop shaping - discrete-time - strictly proper case
#include "sl_sb10zd.cc"  // H-infinity loop shaping - discrete-time - proper case
#include "sl_tb04bd.cc"  // State-space to transfer function conversion
#include "sl_ab04md.cc"  // bilinear transformation
#include "sl_sb10jd.cc"  // descriptor to regular state-space conversion
#include "sl_td04ad.cc"  // transfer function to state-space conversion
#include "sl_tb01ud.cc"  // controllable block Hessenberg realization
#include "sl_ab09hd.cc"  // balanced stochastic truncation model reduction
#include "sl_ab09id.cc"  // balanced truncation & singular perturbation approximation model reduction
#include "sl_ab09jd.cc"  // hankel-norm approximation model reduction
#include "sl_sb16ad.cc"  // balanced truncation & singular perturbation approximation controller reduction
#include "sl_sb16bd.cc"  // coprime factorization state-feedback controller reduction
#include "sl_sb16cd.cc"  // frequency-weighted coprime factorization state-feedback controller reduction
#include "sl_sb10yd.cc"  // fit state-space model to frequency response data
#include "sl_ident.cc"   // system identification
#include "sl_ib01cd.cc"  // compute initial state vector
#include "sl_ib01ad.cc"  // compute singular values
#include "sl_are.cc"     // solve ARE with Schur vector approach and scaling


// stub function to avoid gen_doc_cache warning upon package installation
DEFUN_DLD (__control_slicot_functions__, args, nargout,
   "-*- texinfo -*-\n\
Slicot Release 5.0\n\
No argument checking.\n\
For internal use only.")
{
    octave_value_list retval;
    error ("__control_slicot_functions__: for internal use only");
    return retval;
}

