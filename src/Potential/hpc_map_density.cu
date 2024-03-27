// -*-c++-*-
#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Numerics/split_fixed_precision.h"
#include "map_density.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

#include "Accelerator/syncwarp.cui"
#include "Math/bspline.cui"
#include "Numerics/accumulation.cui"
#include "cellgrid_imaging.cui"

// Compile the __shared__ memory density accumulation kernels.  The format of each name is
// "kSA" + {l,s} + {i, r} + {d,f} + [4, 6] + {d,s} + "MapDensity".  The {l,s} branch indicates
// whether the coordinates in the cell grid have a short (32-bit) or long (64-bit) representation.
// The {i,r} branch indicates whether the coordinate rpresentation is real [r] or fixed-precision
// integer [i].  The final {d,f} branch indicates whether the calculations are to be performed in
// double- or single-precision, and is followed by the interpolation order (each order must get its
// own kernel in the interest of register conservation).  The final {d,s} branch indicates whether
// to carry out accumulations in double (95-bit) or single (63-bit) accumulations, necessary to
// reduce register pressure hat would otherwise be incurred by combining branches of the innermost
// loop and to take all accumulation buffers into __shared__.  Begin with the double-precision
// kernels.
#define ACC_MODE_DOUBLE
#define CALC_MODE_DOUBLE
#define TCALC double
#define TCALC2 double2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsid4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsid5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsid6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlid4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlid5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlid6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrd4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrd5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrd6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlrd4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrd5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrd6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef CALC_MODE_DOUBLE
#undef ACC_MODE_DOUBLE
  
#define CALC_MODE_DOUBLE
#define TCALC double
#define TCALC2 double2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsid4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsid5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsid6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlid4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlid5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlid6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrd4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrd5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrd6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlrd4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrd5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrd6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#      undef DENSITY_SPREADING_THREADS
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef CALC_MODE_DOUBLE
  
// Define the single-precision __shared__ memory density accumulation kernels.
#define ACC_MODE_DOUBLE
#define TCALC float
#define TCALC2 float2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsif4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsif5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsif6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlif4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlif5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlif6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrf4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrf5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrf6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlrf4dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrf5dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrf6dMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef ACC_MODE_DOUBLE
  
#define TCALC float
#define TCALC2 float2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsif4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsif5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsif6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlif4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlif5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlif6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrf4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrf5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrf6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlrf4sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrf5sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrf6sMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC

// Define additional __shared__ accumulation kernels for use when the accumulation can be
// expected to stay within one 32- or 64-bit accumulator, not requiring the overflow bits.
#define SHORT_FORMAT_ACCUMULATION
#define ACC_MODE_DOUBLE
#define CALC_MODE_DOUBLE
#define TCALC double
#define TCALC2 double2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsid4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsid5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsid6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlid4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlid5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlid6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrd4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrd5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrd6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlrd4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrd5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrd6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef CALC_MODE_DOUBLE
#undef ACC_MODE_DOUBLE
  
#define CALC_MODE_DOUBLE
#define TCALC double
#define TCALC2 double2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsid4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsid5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsid6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlid4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlid5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlid6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrd4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrd5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrd6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG 
#    define ORDER 4
#      define KERNEL_NAME kSAlrd4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrd5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrd6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG 
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef CALC_MODE_DOUBLE
  
// Define the single-precision __shared__ memory density accumulation kernels.
#define ACC_MODE_DOUBLE
#define TCALC float
#define TCALC2 float2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsif4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsif5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsif6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlif4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlif5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlif6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrf4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrf5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrf6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlrf4dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrf5dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrf6dsfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef ACC_MODE_DOUBLE
  
#define TCALC float
#define TCALC2 float2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME kSAsif4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsif5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsif6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define T4 llint4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlif4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlif5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlif6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME kSAsrf4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAsrf5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAsrf6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define T4 double4
#  define TMAT_IS_LONG
#    define ORDER 4
#      define KERNEL_NAME kSAlrf4ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME kSAlrf5ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME kSAlrf6ssfMapDensity
#        include "map_density_shracc.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef TMAT_IS_LONG
#  undef T4
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef SHORT_FORMAT_ACCUMULATION

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryShrAccQMapKernelRequirements(const PrecisionModel calc_prec,
                                                            const PrecisionModel acc_prec,
                                                            const bool overflow_needed,
                                                            const size_t cg_tmat,
                                                            const int order) {
  cudaFuncAttributes result;
  switch (calc_prec) {
  case PrecisionModel::DOUBLE:
    switch (acc_prec) {
    case PrecisionModel::DOUBLE:
      if (overflow_needed) {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsid4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsid5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsid6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlid4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlid5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlid6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrd4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrd5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrd6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrd4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrd5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrd6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      else {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsid4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsid5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsid6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlid4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlid5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlid6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrd4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrd5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrd6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrd4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrd5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrd6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      break;
    case PrecisionModel::SINGLE:
      if (overflow_needed) {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsid4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsid5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsid6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlid4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlid5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlid6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrd4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrd5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrd6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrd4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrd5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrd6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      else {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsid4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsid5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsid6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsid6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlid4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlid5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlid6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlid6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrd4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrd5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrd6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrd6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrd4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrd5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrd6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrd6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (acc_prec) {
    case PrecisionModel::DOUBLE:
      if (overflow_needed) {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsif4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsif5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsif6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlif4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlif5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlif6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrf4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrf5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrf6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrf4dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf4dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrf5dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf5dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrf6dMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf6dMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      else {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsif4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsif5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsif6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlif4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlif5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlif6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrf4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrf5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrf6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrf4dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf4dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrf5dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf5dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrf6dsfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf6dsfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      break;
    case PrecisionModel::SINGLE:
      if (overflow_needed) {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsif4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsif5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsif6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlif4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlif5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlif6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrf4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrf5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrf6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrf4sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf4sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrf5sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf5sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrf6sMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf6sMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      else {
        if (cg_tmat == int_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsif4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsif5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsif6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsif6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == llint_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlif4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlif5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlif6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlif6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == float_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAsrf4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAsrf5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAsrf6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAsrf6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
        else if (cg_tmat == double_type_index) {
          if (order == 4 && cudaFuncGetAttributes(&result, kSAlrf4ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf4ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 5 &&
                   cudaFuncGetAttributes(&result, kSAlrf5ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf5ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
          else if (order == 6 &&
                   cudaFuncGetAttributes(&result, kSAlrf6ssfMapDensity) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kSAlrf6ssfMapDensity.",
                  "queryShrAccQMapKernelRequirements");
          }
        }
      }
      break;
    }
    break;
  }
  return result;
}

// Compile the double-precision naive density mapping kernels.  The format of each name is
// "k" + {l,s} + {i, r} + {d,f} + [4, 6] + "MapDensity".  The {l,s} branch indicates whether the
// coordinates in the cell grid have a short (32-bit) or long (64-bit) representation.  The {i,r}
// branch indicates whether the coordinate rpresentation is real [r] or fixed-precision integer
// [i].  The final {d,f} branch indicates whether the calculations are to be performed in single-
// or double-precision.  The letter codes are followed by the interpolation order (each order must
// get its own interpolation order in the interest of register conservation).
#define TCALC double
#define TCALC2 double2
#define TCALC_IS_DOUBLE
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME ksid4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksid5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksid6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define TMAT_IS_LONG
#  define T4 llint4
#    define ORDER 4
#      define KERNEL_NAME klid4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klid5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klid6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME ksrd4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksrd5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksrd6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define TMAT_IS_LONG
#  define T4 double4
#    define ORDER 4
#      define KERNEL_NAME klrd4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klrd5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klrd6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC_IS_DOUBLE
#undef TCALC2
#undef TCALC

// Compile the single-precision naive density mapping kernels
#define TCALC float
#define TCALC2 float2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME ksif4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksif5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksif6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define TMAT_IS_LONG
#  define T4 llint4
#    define ORDER 4
#      define KERNEL_NAME klif4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klif5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klif6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME ksrf4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksrf5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksrf6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define TMAT_IS_LONG
#  define T4 double4
#    define ORDER 4
#      define KERNEL_NAME klrf4MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klrf5MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klrf6MapDensity
#      include "map_density.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryGeneralQMapKernelRequirements(const PrecisionModel prec,
                                                             const size_t cg_tmat,
                                                             const int order) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (cg_tmat == int_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksid4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksid4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksid5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksid5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksid6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksid6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    else if (cg_tmat == llint_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klid4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klid4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klid5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klid5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klid6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klid6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    else if (cg_tmat == float_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksrd4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrd4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksrd5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrd5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksrd6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrd6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    else if (cg_tmat == double_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klrd4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrd4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klrd5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrd5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klrd6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrd6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    break;
  case PrecisionModel::SINGLE:
    if (cg_tmat == int_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksif4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksif4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksif5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksif5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksif6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksif6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    else if (cg_tmat == llint_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klif4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klif4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klif5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klif5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klif6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klif6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    else if (cg_tmat == float_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksrf4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrf4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksrf5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrf5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksrf6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrf6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    else if (cg_tmat == double_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klrf4MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrf4MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klrf5MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrf5MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klrf6MapDensity) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrf6MapDensity.",
              "queryGeneralQMapKernelRequirements");
      }
    }
    break;
  }
  return result;
}  

//-------------------------------------------------------------------------------------------------
extern void launchShrAccDensityKernel(PMIGridWriter *pm_wrt, const bool overflow,
                                      MMControlKit<double> *ctrl,
                                      const CellGridReader<void, void, void, void> &v_cgr,
                                      const size_t cg_tmat,
                                      const SyNonbondedKit<double, double2> &synbk,
                                      const int2 lp) {
  matchThemes(pm_wrt->theme, v_cgr.theme);
  switch (pm_wrt->mode) {
  case PrecisionModel::DOUBLE:
    if (cg_tmat == int_type_index) {
      const CellGridReader<int, void, double, int4> cgr = restoreType<int, void,
                                                                      double, int4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsid4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsid4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsid5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsid5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsid6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsid6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == llint_type_index) {
      const CellGridReader<llint, void, double, llint4> cgr = restoreType<llint, void,
                                                                          double, llint4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlid4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlid4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlid5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlid5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlid6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlid6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == float_type_index) {
      const CellGridReader<float, void, double, float4> cgr = restoreType<float, void,
                                                                          double, float4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsrd4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrd4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsrd5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrd5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsrd6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrd6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == double_type_index) {
      const CellGridReader<double, void,
                           double, double4> cgr = restoreType<double, void,
                                                              double, double4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlrd4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrd4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlrd5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrd5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlrd6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrd6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    if (cg_tmat == int_type_index) {
      const CellGridReader<int, void, double, int4> cgr = restoreType<int, void,
                                                                      double, int4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsid4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsid4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsid5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsid5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsid6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsid6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == llint_type_index) {
      const CellGridReader<llint, void, double, llint4> cgr = restoreType<llint, void,
                                                                          double, llint4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlid4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlid4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlid5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlid5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlid6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlid6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == float_type_index) {
      const CellGridReader<float, void, double, float4> cgr = restoreType<float, void,
                                                                          double, float4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsrd4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrd4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsrd5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrd5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsrd6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrd6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == double_type_index) {
      const CellGridReader<double, void,
                           double, double4> cgr = restoreType<double, void,
                                                              double, double4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlrd4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrd4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlrd5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrd5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlrd6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrd6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchShrAccDensityKernel(PMIGridWriter *pm_wrt, const bool overflow,
                                      MMControlKit<float> *ctrl,
                                      const CellGridReader<void, void, void, void> &v_cgr,
                                      const size_t cg_tmat,
                                      const SyNonbondedKit<float, float2> &synbk, const int2 lp) {
  matchThemes(pm_wrt->theme, v_cgr.theme);
  switch (pm_wrt->mode) {
  case PrecisionModel::DOUBLE:
    if (cg_tmat == int_type_index) {
      const CellGridReader<int, void, float, int4> cgr = restoreType<int, void,
                                                                     float, int4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsif4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsif4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsif5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsif5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsif6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsif6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == llint_type_index) {
      const CellGridReader<llint, void, float, llint4> cgr = restoreType<llint, void,
                                                                         float, llint4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlif4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlif4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlif5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlif5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlif6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlif6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == float_type_index) {
      const CellGridReader<float, void, float, float4> cgr = restoreType<float, void,
                                                                         float, float4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsrf4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrf4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsrf5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrf5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsrf6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrf6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == double_type_index) {
      const CellGridReader<double, void, float, double4> cgr = restoreType<double, void,
                                                                           float, double4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlrf4dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrf4dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlrf5dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrf5dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlrf6dMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrf6dsfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    if (cg_tmat == int_type_index) {
      const CellGridReader<int, void, float, int4> cgr = restoreType<int, void,
                                                                     float, int4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsif4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsif4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsif5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsif5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsif6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsif6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == llint_type_index) {
      const CellGridReader<llint, void, float, llint4> cgr = restoreType<llint, void,
                                                                         float, llint4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlif4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlif4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlif5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlif5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlif6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlif6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == float_type_index) {
      const CellGridReader<float, void, float, float4> cgr = restoreType<float, void,
                                                                         float, float4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAsrf4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrf4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAsrf5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrf5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAsrf6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAsrf6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    else if (cg_tmat == double_type_index) {
      const CellGridReader<double, void, float, double4> cgr = restoreType<double, void,
                                                                           float, double4>(v_cgr);
      switch (pm_wrt->order) {
      case 4:
        if (overflow) kSAlrf4sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrf4ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 5:
        if (overflow) kSAlrf5sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrf5ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      case 6:
        if (overflow) kSAlrf6sMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        else kSAlrf6ssfMapDensity<<<lp.x, lp.y>>>(*pm_wrt, *ctrl, cgr, synbk);
        break;
      default:
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchGenPrpDensityKernel(PMIGridAccumulator *pm_acc,
                                      const CellGridReader<void, void, void, void> &v_cgr,
                                      const size_t cg_tmat,
                                      const SyNonbondedKit<double, double2> &synbk,
                                      const int2 lp) {
  matchThemes(pm_acc->theme, v_cgr.theme);
  if (cg_tmat == int_type_index) {
    const CellGridReader<int, void, double, int4> cgr = restoreType<int, void,
                                                                    double, int4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      ksid4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      ksid5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      ksid6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == llint_type_index) {
    const CellGridReader<llint, void, double, llint4> cgr = restoreType<llint, void,
                                                                        double, llint4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      klid4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      klid5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      klid6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == float_type_index) {
    const CellGridReader<float, void, double, float4> cgr = restoreType<float, void,
                                                                        double, float4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      ksrd4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      ksrd5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      ksrd6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == double_type_index) {
    const CellGridReader<double, void, double, double4> cgr = restoreType<double, void,
                                                                          double, double4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      klrd4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      klrd5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      klrd6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchGenPrpDensityKernel(PMIGridAccumulator *pm_acc,
                                      const CellGridReader<void, void, void, void> &v_cgr,
                                      const size_t cg_tmat,
                                      const SyNonbondedKit<float, float2> &synbk, const int2 lp) {
  matchThemes(pm_acc->theme, v_cgr.theme);
  if (cg_tmat == int_type_index) {
    const CellGridReader<int, void, float, int4> cgr = restoreType<int, void,
                                                                   float, int4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      ksif4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      ksif5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      ksif6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == llint_type_index) {
    const CellGridReader<llint, void, float, llint4> cgr = restoreType<llint, void,
                                                                       float, llint4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      klif4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      klif5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      klif6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == float_type_index) {
    const CellGridReader<float, void, float, float4> cgr = restoreType<float, void,
                                                                       float, float4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      ksrf4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      ksrf5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      ksrf6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == double_type_index) {
    const CellGridReader<double, void, float, double4> cgr = restoreType<double, void,
                                                                         float, double4>(v_cgr);
    switch (pm_acc->order) {
    case 4:
      klrf4MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 5:
      klrf5MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    case 6:
      klrf6MapDensity<<<lp.x, lp.y>>>(*pm_acc, cgr, synbk);
      break;
    default:
      break;
    }
  }
}

} // namespace energy
} // namespace stormm
