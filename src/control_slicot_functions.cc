/* TODO: * create common.h
 *       * add common.cc to makefile (control_slicot_functions.oct)
 */

#include <octave/oct.h>
#include "common.cc"    // common code for all functions of the oct-file

#include "slab08nd.cc"  // transmission zeros of state-space models
#include "slab13dd.cc"  // L-infinity norm
#include "slsb10hd.cc"  // H-2 controller synthesis - continuous-time
#include "slsb10ed.cc"  // H-2 controller synthesis - discrete-time
#include "slab13bd.cc"  // H-2 norm
#include "slsb01bd.cc"  // Pole assignment
#include "slsb10fd.cc"  // H-infinity controller synthesis - continuous-time
#include "slsb10dd.cc"  // H-infinity controller synthesis - discrete-time
#include "slsb03md.cc"  // Lyapunov equations
#include "slsb04md.cc"  // Sylvester equations - continuous-time
#include "slsb04qd.cc"  // Sylvester equations - discrete-time
#include "slsg03ad.cc"  // generalized Lyapunov equations
#include "slsb02od.cc"  // algebraic Riccati equations
#include "slab13ad.cc"  // Hankel singular values
#include "slab01od.cc"  // staircase form using orthogonal transformations
#include "sltb01pd.cc"  // minimal realization of state-space models
#include "slsb03od.cc"  // Cholesky factor of Lyapunov equations
#include "slsg03bd.cc"  // Cholesky factor of generalized Lyapunov equations
#include "slag08bd.cc"  // transmission zeros of descriptor state-space models
#include "sltg01jd.cc"  // minimal realization of descriptor state-space models
#include "sltg01hd.cc"  // controllability staircase form of descriptor state-space models
#include "sltg01id.cc"  // observability staircase form of descriptor state-space models
#include "slsg02ad.cc"  // solution of algebraic Riccati equations for descriptor systems
#include "sltg04bx.cc"  // gain of descriptor state-space models
#include "sltb01id.cc"  // scaling of state-space models
#include "sltg01ad.cc"  // scaling of descriptor state-space models
#include "slsb10id.cc"  // H-infinity loop shaping - continuous-time
#include "slsb10kd.cc"  // H-infinity loop shaping - discrete-time - strictly proper case
#include "slsb10zd.cc"  // H-infinity loop shaping - discrete-time - proper case
#include "sltb04bd.cc"  // State-space to transfer function conversion
#include "slab04md.cc"  // bilinear transformation
#include "slsb10jd.cc"  // descriptor to regular state-space conversion
#include "sltd04ad.cc"  // transfer function to state-space conversion
#include "sltb01ud.cc"  // controllable block Hessenberg realization
#include "slab09hd.cc"  // balanced stochastic truncation model reduction
#include "slab09id.cc"  // balanced truncation & singular perturbation approximation model reduction
#include "slab09jd.cc"  // hankel-norm approximation model reduction
#include "slsb16ad.cc"  // balanced truncation & singular perturbation approximation controller reduction
#include "slsb16bd.cc"  // coprime factorization state-feedback controller reduction
#include "slsb16cd.cc"  // frequency-weighted coprime factorization state-feedback controller reduction
#include "slsb10yd.cc"  // fit state-space model to frequency response data
