// -*-c++-*-
#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Numerics/split_fixed_precision.h"
#include "cellgrid.h"
#include "gather_forces.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

#include "Accelerator/syncwarp.cui"
#include "Math/bspline.cui"
#include "Numerics/accumulation.cui"
#include "cellgrid_imaging.cui"

// Compile double-precision naive force interpolation kernels.  The format of each name is
// "k" + {l, s} + {i, r} + {d, f} + [4, 6] + "GatherForces".  As in density mapping, the {l, s}
// branch indicates whether the coordinates in the cell grid have a short (32-bit) or long (64-bit)
// representation, while the following {i, r} letters denote fixed-precision integral [i] or real
// [r] representations in those local particle positions.  The final {d, f} branch indicates
// whether calculations are to be performed in single-precision [f] or double-precision [d].  Each
// series of letter codes is followed by a specific interpolation order (which are enumerated in
// the interest of register conservation).
#define INTERPOLATE_FORCES
#define TACC llint
#define TCALC double
#define TCALC2 double2
#define TCALC_IS_DOUBLE
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME ksid4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksid5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksid6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define TMAT_IS_LONG
#  define T4 llint4
#    define ORDER 4
#      define KERNEL_NAME klid4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klid5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klid6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME ksrd4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksrd5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksrd6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define TMAT_IS_LONG
#  define T4 double4_16a
#    define ORDER 4
#      define KERNEL_NAME klrd4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klrd5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klrd6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC_IS_DOUBLE
#undef TCALC2
#undef TCALC
#undef TACC
  
// Compile the single-precision force interpolation kernels.  See above for definitions of the
// alphanumeric codes in each kernel name.
#define TACC int
#define TCALC float
#define TCALC2 float2
#  define TMAT int
#  define T4 int4
#    define ORDER 4
#      define KERNEL_NAME ksif4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksif5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksif6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT llint
#  define TMAT_IS_LONG
#  define T4 llint4
#    define ORDER 4
#      define KERNEL_NAME klif4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klif5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klif6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  define TMAT_IS_REAL
#  define TMAT float
#  define T4 float4
#    define ORDER 4
#      define KERNEL_NAME ksrf4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME ksrf5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME ksrf6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT
#  define TMAT double
#  define TMAT_IS_LONG
#  define T4 double4_16a
#    define ORDER 4
#      define KERNEL_NAME klrf4GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 5
#      define KERNEL_NAME klrf5GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#    define ORDER 6
#      define KERNEL_NAME klrf6GatherForces
#      include "map_pmi_basic.cui"
#      undef KERNEL_NAME
#    undef ORDER
#  undef T4
#  undef TMAT_IS_LONG
#  undef TMAT
#  undef TMAT_IS_REAL
#undef TCALC2
#undef TCALC
#undef TACC
#undef INTERPOLATE_FORCES
  
//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryGeneralForceGatheringKernelRequirements(const PrecisionModel prec,
                                                                       const size_t cg_tmat,
                                                                       const int order) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (cg_tmat == int_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksid4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksid4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksid5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksid5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksid6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksid6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    else if (cg_tmat == llint_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klid4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klid4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klid5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klid5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klid6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klid6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    else if (cg_tmat == float_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksrd4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrd4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksrd5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrd5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksrd6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrd6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    else if (cg_tmat == double_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klrd4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrd4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klrd5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrd5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klrd6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrd6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    break;
  case PrecisionModel::SINGLE:
    if (cg_tmat == int_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksif4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksif4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksif5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksif5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksif6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksif6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    else if (cg_tmat == llint_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klif4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klif4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klif5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klif5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klif6GatherForces) != cudaSuccess) {
	rtErr("Error obtaining attributes for kernel klif6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    else if (cg_tmat == float_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, ksrf4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrf4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, ksrf5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrf5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, ksrf6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ksrf6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    else if (cg_tmat == double_type_index) {
      if (order == 4 && cudaFuncGetAttributes(&result, klrf4GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrf4GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 5 && cudaFuncGetAttributes(&result, klrf5GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrf5GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
      else if (order == 6 && cudaFuncGetAttributes(&result, klrf6GatherForces) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel klrf6GatherForces.",
              "queryGeneralForceGatheringKernelRequirements");
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchGenForceGatheringKernel(CellGridWriter<void, void, void, void> *v_cgw,
                                          const PMIGridReader &pm_rdr, const size_t cg_tmat,
                                          const PsSynthesisBorders &pssb,
                                          const SyNonbondedKit<double, double2> &synbk,
                                          const int2 lp) {
  matchThemes(pm_rdr.theme, v_cgw->theme);
  if (cg_tmat == int_type_index) {
    CellGridWriter<int, llint, double, int4> cgw = restoreType<int, llint, double, int4>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      ksid4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      ksid5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      ksid6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == llint_type_index) {
    CellGridWriter<llint, llint, double, llint4> cgw = restoreType<llint, llint,
                                                                   double, llint4>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      klid4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      klid5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      klid6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == float_type_index) {
    CellGridWriter<float, llint, double, float4> cgw = restoreType<float, llint,
                                                                   double, float4>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      ksrd4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      ksrd5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      ksrd6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == double_type_index) {
    CellGridWriter<double, llint, double, double4_16a> cgw = restoreType<double, llint,
                                                                     double, double4_16a>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      klrd4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      klrd5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      klrd6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchGenForceGatheringKernel(CellGridWriter<void, void, void, void> *v_cgw,
                                          const PMIGridReader &pm_rdr, const size_t cg_tmat,
                                          const PsSynthesisBorders &pssb,
                                          const SyNonbondedKit<float, float2> &synbk,
                                          const int2 lp) {
  matchThemes(pm_rdr.theme, v_cgw->theme);
  if (cg_tmat == int_type_index) {
    CellGridWriter<int, int, float, int4> cgw = restoreType<int, int, float, int4>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      ksif4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      ksif5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      ksif6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == llint_type_index) {
    CellGridWriter<llint, int, float, llint4> cgw = restoreType<llint, int, float, llint4>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      klif4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      klif5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      klif6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == float_type_index) {
    CellGridWriter<float, int, float, float4> cgw = restoreType<float, int, float, float4>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      ksrf4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      ksrf5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      ksrf6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
  else if (cg_tmat == double_type_index) {
    CellGridWriter<double, int, float, double4_16a> cgw = restoreType<double, int,
                                                                  float, double4_16a>(v_cgw);
    switch (pm_rdr.order) {
    case 4:
      klrf4GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 5:
      klrf5GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    case 6:
      klrf6GatherForces<<<lp.x, lp.y>>>(cgw, pm_rdr, pssb, synbk);
      break;
    default:
      break;
    }
  }
}

} // namespace energy
} // namespace stormm
