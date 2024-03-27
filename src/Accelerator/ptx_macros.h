// -*-c++-*-
#ifndef STORMM_PTX_MACROS_H
#define STORMM_PTX_MACROS_H

#include "copyright.h"

#ifdef STORMM_USE_CUDA
#  define SHFL_DOWN(a, b)    __shfl_down_sync(0xffffffff, a, b, 32)
#  define SHFL_XOR(a, b)     __shfl_xor_sync(0xffffffff, a, b, 32)
#  define SHFL_UP(a, b)      __shfl_up_sync(0xffffffff, a, b, 32)
#  define SHFL(a, b)         __shfl_sync(0xffffffff, a, b, 32)
#  define SYNCWARP           syncWarp()
#  define BALLOT(predicate)  __ballot_sync(0xffffffff, predicate)

// Reduce all elements of the warp such that the sum appears in lane 0
#  define WARP_REDUCE_DOWN(var) { \
  var += __shfl_down_sync(0xffffffff, var, 16, 32); \
  var += __shfl_down_sync(0xffffffff, var,  8, 32); \
  var += __shfl_down_sync(0xffffffff, var,  4, 32); \
  var += __shfl_down_sync(0xffffffff, var,  2, 32); \
  var += __shfl_down_sync(0xffffffff, var,  1, 32); \
}

// Compute an inclusive prefix sum over all elements of the warp.  An example of this operation:
// [  5  3  9  1  0  5 ] -> [  5  8 17 18 18 23 ]
#  define INCLUSIVE_WARP_PREFIXSUM(var, tgx) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up_sync(0xffffffff, var, 1, 32); \
  var += ((tgx &  3) ==  3) * __shfl_up_sync(0xffffffff, var, 2, 32); \
  var += ((tgx &  7) ==  7) * __shfl_up_sync(0xffffffff, var, 4, 32); \
  var += ((tgx & 15) == 15) * __shfl_up_sync(0xffffffff, var, 8, 32); \
  var += (tgx == 31) * __shfl_up_sync(0xffffffff, var, 16, 32); \
  var += ((tgx & 15) == 7 && tgx > 16) * __shfl_up_sync(0xffffffff, var, 8, 32); \
  var += ((tgx &  7) == 3 && tgx > 8)  * __shfl_up_sync(0xffffffff, var, 4, 32); \
  var += ((tgx &  3) == 1 && tgx > 4)  * __shfl_up_sync(0xffffffff, var, 2, 32); \
  var += ((tgx &  1) == 0 && tgx >= 2) * __shfl_up_sync(0xffffffff, var, 1, 32); \
}

// Compute an exclusive prefix sum over all elements of the warp.  If the total of the prefix sum
// is required, i.e. the upper limit of the final bin, use EXLCUSIVE_WARP_PREFIXSUM_SAVETOTAL.  An
// example of this operation: [  5  3  9  1  0  5 ] -> [  0  5  8 17 18 18 23 ]
#  define EXCLUSIVE_WARP_PREFIXSUM(var, tgx) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up_sync(0xffffffff, var, 1, 32); \
  var += ((tgx &  3) ==  3) * __shfl_up_sync(0xffffffff, var, 2, 32); \
  var += ((tgx &  7) ==  7) * __shfl_up_sync(0xffffffff, var, 4, 32); \
  var += ((tgx & 15) == 15) * __shfl_up_sync(0xffffffff, var, 8, 32); \
  var += (tgx == 31) * __shfl_up_sync(0xffffffff, var, 16, 32); \
  var += ((tgx & 15) == 7 && tgx > 16) * __shfl_up_sync(0xffffffff, var, 8, 32); \
  var += ((tgx &  7) == 3 && tgx > 8)  * __shfl_up_sync(0xffffffff, var, 4, 32); \
  var += ((tgx &  3) == 1 && tgx > 4)  * __shfl_up_sync(0xffffffff, var, 2, 32); \
  var += ((tgx &  1) == 0 && tgx >= 2) * __shfl_up_sync(0xffffffff, var, 1, 32); \
  var = __shfl_up_sync(0xffffffff, var, 1, 32); \
  if (tgx == 0) { \
    var = 0; \
  } \
}

// Compute an exclusive prefix sum over all elements of the warp and retain the total in an
// auxiliary variable which is then broadcast to all threads.
#  define EXCLUSIVE_WARP_PREFIXSUM_SAVETOTAL(var, tgx, result) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up_sync(0xffffffff, var, 1, 32); \
  var += ((tgx &  3) ==  3) * __shfl_up_sync(0xffffffff, var, 2, 32); \
  var += ((tgx &  7) ==  7) * __shfl_up_sync(0xffffffff, var, 4, 32); \
  var += ((tgx & 15) == 15) * __shfl_up_sync(0xffffffff, var, 8, 32); \
  var += (tgx == 31) * __shfl_up_sync(0xffffffff, var, 16, 32); \
  var += ((tgx & 15) == 7 && tgx > 16) * __shfl_up_sync(0xffffffff, var, 8, 32); \
  var += ((tgx &  7) == 3 && tgx > 8)  * __shfl_up_sync(0xffffffff, var, 4, 32); \
  var += ((tgx &  3) == 1 && tgx > 4)  * __shfl_up_sync(0xffffffff, var, 2, 32); \
  var += ((tgx &  1) == 0 && tgx >= 2) * __shfl_up_sync(0xffffffff, var, 1, 32); \
  result = __shfl_sync(0xffffffff, var, 31, 32); \
  var = __shfl_up_sync(0xffffffff, var,  1, 32); \
  if (tgx == 0) { \
    var = 0; \
  } \
}
#endif // STORMM_USE_CUDA

#ifdef STORMM_USE_HIP
#  define SHFL_DOWN(a, b)    __shfl_down(0xffffffffffffffff, a, b, 64)
#  define SHFL_XOR(a, b)     __shfl_xor(0xffffffffffffffff, a, b, 64)
#  define SHFL_UP(a, b)      __shfl_up(0xffffffffffffffff, a, b, 64)
#  define SHFL(a, b)         __shfl(0xffffffffffffffff, a, b, 64)
#  define SYNCWARP
#  define BALLOT(predicate)  __ballot(0xffffffffffffffff, predicate)

// Reduce all elements of the warp such that the sum appears in lane 0
#  define WARP_REDUCE_DOWN(var) { \
  var += __shfl_down(var, 32); \
  var += __shfl_down(var, 16); \
  var += __shfl_down(var,  8); \
  var += __shfl_down(var,  4); \
  var += __shfl_down(var,  2); \
  var += __shfl_down(var,  1); \
}

// Compute an inclusive prefix sum over all elements of the warp.  An example of this operation:
// [  5  3  9  1  0  5 ] -> [  5  8 17 18 18 23 ]
#  define INCLUSIVE_WARP_PREFIXSUM(var, tgx) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up(var,  1); \
  var += ((tgx &  3) ==  3) * __shfl_up(var,  2); \
  var += ((tgx &  7) ==  7) * __shfl_up(var,  4); \
  var += ((tgx & 15) == 15) * __shfl_up(var,  8); \
  var += ((tgx & 31) == 31) * __shfl_up(var, 16); \
  var += (tgx == 63) * __shfl_up(var, 32); \
  var += ((tgx & 31) == 15 && tgx > 32) * __shfl_up(var, 16); \
  var += ((tgx & 15) ==  7 && tgx > 16) * __shfl_up(var,  8); \
  var += ((tgx &  7) ==  3 && tgx >  8) * __shfl_up(var,  4); \
  var += ((tgx &  3) ==  1 && tgx >  4) * __shfl_up(var,  2); \
  var += ((tgx &  1) ==  0 && tgx >= 2) * __shfl_up(var,  1); \
}

// Compute an exclusive prefix sum over all elements of the warp.  If the total of the prefix sum
// is required, i.e. the upper limit of the final bin, use EXLCUSIVE_WARP_PREFIXSUM_SAVETOTAL.  An
// example of this operation: [  5  3  9  1  0  5 ] -> [  0  5  8 17 18 18 23 ]
#  define EXCLUSIVE_WARP_PREFIXSUM(var, tgx) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up(var,  1); \
  var += ((tgx &  3) ==  3) * __shfl_up(var,  2); \
  var += ((tgx &  7) ==  7) * __shfl_up(var,  4); \
  var += ((tgx & 15) == 15) * __shfl_up(var,  8); \
  var += ((tgx & 31) == 31) * __shfl_up(var, 16); \
  var += (tgx == 63) * __shfl_up(var, 32); \
  var += ((tgx & 31) == 15 && tgx > 32) * __shfl_up(var, 16); \
  var += ((tgx & 15) ==  7 && tgx > 16) * __shfl_up(var,  8); \
  var += ((tgx &  7) ==  3 && tgx >  8) * __shfl_up(var,  4); \
  var += ((tgx &  3) ==  1 && tgx >  4) * __shfl_up(var,  2); \
  var += ((tgx &  1) ==  0 && tgx >= 2) * __shfl_up(var,  1); \
  var = __shfl_up(var, 1); \
  if (tgx == 0) { \
    var = 0; \
  } \
}

// Compute an exclusive prefix sum over all elements of the warp and retain the total in an
// auxiliary variable which is then broadcast to all threads.
#  define EXCLUSIVE_WARP_PREFIXSUM_SAVETOTAL(var, tgx, result) \
{ \
  var += ((tgx &  1) ==  1) * __shfl_up(var,  1); \
  var += ((tgx &  3) ==  3) * __shfl_up(var,  2); \
  var += ((tgx &  7) ==  7) * __shfl_up(var,  4); \
  var += ((tgx & 15) == 15) * __shfl_up(var,  8); \
  var += ((tgx & 31) == 31) * __shfl_up(var, 16); \
  var += (tgx == 63) * __shfl_up(var, 32); \
  var += ((tgx & 31) == 15 && tgx > 32) * __shfl_up(var, 16); \
  var += ((tgx & 15) ==  7 && tgx > 16) * __shfl_up(var,  8); \
  var += ((tgx &  7) ==  3 && tgx >  8) * __shfl_up(var,  4); \
  var += ((tgx &  3) ==  1 && tgx >  4) * __shfl_up(var,  2); \
  var += ((tgx &  1) ==  0 && tgx >= 2) * __shfl_up(var,  1); \
  result = __shfl(var, 63); \
  var = __shfl_up(var, 1); \
  if (tgx == 0) { \
    var = 0; \
  } \
}
#endif // STORMM_USE_HIP

#endif
