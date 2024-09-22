// -*-c++-*-
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "cellgrid.h"

namespace stormm {
namespace energy {

#include "Numerics/accumulation.cui"
  
/// \brief Templated, simple kernel to initialize forces within a CellGrid object.  The process in
///        this kernel will be replicated in other, fused kernels for highly optimized processes,
///        but this standalone capability facilitates prototyping and helps illustrate the manner
///        in which the CellGrid operates.  This kernel will initialize forces for the current
///        cell layout, based on the production of the CellGrid abstract.  Thus, an abstract that
///        has been produced based on the original CellGrid object's current time cycle setting is
///        what should be submitted, once the contents of image_cell_limits or
///        image_cell_limits_alt have been settled.  This routine will not initialize warp-specific
///        allocations for local work to be used during the processing of each cell's neighbor
///        list, which must occur repeatedly over the course of a force calculation, each time a
///        warp moves to work on a new cell.
///
/// \param cgw  Writeable abstract of the cell grid, containing pointers to forces to initialize
template <typename T, typename Tacc, typename Tcalc, typename T4>
__global__ void __launch_bounds__(large_block_size, 1)
kInitializeForces(CellGridWriter<T, Tacc, Tcalc, T4> cgw) {

  // Each warp will proceed to initialize one cell.  No asynchronicity is needed for this kernel as
  // optimization, which in this case would be keeping the memory bus saturated throughout the life
  // of the kernel, is not critical here.
  const int warps_per_kernel = ((blockDim.x * gridDim.x) >> warp_bits);
  int warp_pos = (((blockIdx.x * blockDim.x) + threadIdx.x) >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  const Tacc value_zero = 0;
  while (warp_pos < cgw.total_cell_count) {
    const uint2 cell_lims = cgw.cell_limits[warp_pos];
    const uint natom = (cell_lims.y >> 16);
    const uint hlim = cell_lims.x + natom;
    for (uint image_idx = cell_lims.x + lane_idx; image_idx < hlim; image_idx += warp_size_int) {
      cgw.xfrc[image_idx] = value_zero;
      cgw.yfrc[image_idx] = value_zero;
      cgw.zfrc[image_idx] = value_zero;
      cgw.xfrc_ovrf[image_idx] = 0;
      cgw.yfrc_ovrf[image_idx] = 0;
      cgw.zfrc_ovrf[image_idx] = 0;
    }
    warp_pos += warps_per_kernel;
  }
}

/// \brief Templated, simple kernel for transferring forces from a CellGrid object to all
///        corresponing atoms in the associated coordiante synthesis.
///
/// \param cgr       Read-only abstract for the neighbor list object
/// \param poly_psw  Writeable abstract of the coordinate synthesis
template <typename T, typename Tacc, typename Tcalc, typename T4>
__global__ void __launch_bounds__(large_block_size, 1)
kContributeForces(const CellGridReader<T, Tacc, Tcalc, T4> cgr, const size_t ct_acc,
                  const size_t int_type_code, PsSynthesisWriter poly_psw) {
  int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
  const int last_sys = poly_psw.system_count - 1;
  const int atom_limit = poly_psw.atom_starts[last_sys] + poly_psw.atom_counts[last_sys];
  const size_t image_limit = cgr.total_cell_count * cgr.cell_base_capacity;
  while (pos < atom_limit) {
    const size_t img_idx = __ldcv(&cgr.img_atom_idx[pos]);
    if (img_idx < image_limit) {
      if (ct_acc == int_type_code) {
        const llint ixfrc = int63ToLongLong(__ldg(&cgr.xfrc[img_idx]),
                                            __ldg(&cgr.xfrc_ovrf[img_idx]));
        const llint iyfrc = int63ToLongLong(__ldg(&cgr.yfrc[img_idx]),
                                            __ldg(&cgr.yfrc_ovrf[img_idx]));
        const llint izfrc = int63ToLongLong(__ldg(&cgr.zfrc[img_idx]),
                                            __ldg(&cgr.zfrc_ovrf[img_idx]));
        if (poly_psw.frc_bits <= force_scale_nonoverflow_bits) {
          __stwt(&poly_psw.xfrc[pos], __ldcv(&poly_psw.xfrc[pos]) + ixfrc);
          __stwt(&poly_psw.yfrc[pos], __ldcv(&poly_psw.yfrc[pos]) + iyfrc);
          __stwt(&poly_psw.zfrc[pos], __ldcv(&poly_psw.zfrc[pos]) + izfrc);
        }
        else {
          const int95_t n_ixfrc = int95Sum(__ldcv(&poly_psw.xfrc[pos]),
                                           __ldcv(&poly_psw.xfrc_ovrf[pos]), ixfrc, 0);
          const int95_t n_iyfrc = int95Sum(__ldcv(&poly_psw.yfrc[pos]),
                                           __ldcv(&poly_psw.yfrc_ovrf[pos]), iyfrc, 0);
          const int95_t n_izfrc = int95Sum(__ldcv(&poly_psw.zfrc[pos]),
                                           __ldcv(&poly_psw.zfrc_ovrf[pos]), izfrc, 0);
          __stwt(&poly_psw.xfrc[pos], n_ixfrc.x);
          __stwt(&poly_psw.yfrc[pos], n_iyfrc.x);
          __stwt(&poly_psw.zfrc[pos], n_izfrc.x);
          __stwt(&poly_psw.xfrc_ovrf[pos], n_ixfrc.y);
          __stwt(&poly_psw.yfrc_ovrf[pos], n_iyfrc.y);
          __stwt(&poly_psw.zfrc_ovrf[pos], n_izfrc.y);
        }
      }
      else {
        if (poly_psw.frc_bits <= force_scale_nonoverflow_bits) {
          __stwt(&poly_psw.xfrc[pos],
                 __ldcv(&poly_psw.xfrc[pos]) + __ldg(&cgr.xfrc[img_idx]));
          __stwt(&poly_psw.yfrc[pos],
                 __ldcv(&poly_psw.yfrc[pos]) + __ldg(&cgr.yfrc[img_idx]));
          __stwt(&poly_psw.zfrc[pos],
                 __ldcv(&poly_psw.zfrc[pos]) + __ldg(&cgr.zfrc[img_idx]));
        }
        else {
          const int95_t n_ixfrc = int95Sum(__ldcv(&poly_psw.xfrc[pos]),
                                           __ldcv(&poly_psw.xfrc_ovrf[pos]),
                                           __ldg(&cgr.xfrc[img_idx]),
                                           __ldg(&cgr.xfrc_ovrf[img_idx]));
          const int95_t n_iyfrc = int95Sum(__ldcv(&poly_psw.yfrc[pos]),
                                           __ldcv(&poly_psw.yfrc_ovrf[pos]),
                                           __ldg(&cgr.yfrc[img_idx]),
                                           __ldg(&cgr.yfrc_ovrf[img_idx]));
          const int95_t n_izfrc = int95Sum(__ldcv(&poly_psw.zfrc[pos]),
                                           __ldcv(&poly_psw.zfrc_ovrf[pos]),
                                           __ldg(&cgr.zfrc[img_idx]),
                                           __ldg(&cgr.zfrc_ovrf[img_idx]));
          __stwt(&poly_psw.xfrc[pos], n_ixfrc.x);
          __stwt(&poly_psw.yfrc[pos], n_iyfrc.x);
          __stwt(&poly_psw.zfrc[pos], n_izfrc.x);
          __stwt(&poly_psw.xfrc_ovrf[pos], n_ixfrc.y);
          __stwt(&poly_psw.yfrc_ovrf[pos], n_iyfrc.y);
          __stwt(&poly_psw.zfrc_ovrf[pos], n_izfrc.y);
        }
      }
    }
    pos += blockDim.x * gridDim.x;
  }
}

/// \brief Evaluate a switch over possible actions related to the CellGrid object.
///
/// Overloaded:
///   - Provide only the CellGrid abstract for processes involving only its organization
///   - Provide a read-only abstract for the coordinate synthesis to perform selected processes
///   - Provide a writeable abstract for the coordinate synthesis to perform selected processes
///
/// \param cgw       The cell grid abstract, containing accumulator arrays and all critical bounds
/// \param cgr       The cell grid abstract, containing accumulator arrays and all critical bounds
///                  but lacking the prior image and working arrays (as the struct is read-only)
/// \param poly_psw  Writeable abstract for the coordiante synthesis (it is trusted that the bounds
///                  of systems in this object correspond to those used to construct the cell grid)
/// \param poly_psr  Read-only abstract for the coordiante synthesis (it is trusted that the bounds
///                  of systems in this object correspond to those used to construct the cell grid)
/// \param gpu       Details of the GPU that will perform the initializations
/// \param process   The action to perform (certain actions in overloads with improper inputs, e.g.
///                  lacking a writeable coordinate synthesis, will raise runtime errors)
/// \{
template <typename T, typename Tacc, typename Tcalc, typename T4>
void executeCellGridAction(CellGridWriter<T, Tacc, Tcalc, T4> *cgw, const GpuDetails &gpu,
                           const CellGridAction process) {
  bool problem = false;
  switch (process) {
  case CellGridAction::INIT_FORCES:
    kInitializeForces<T, Tacc, Tcalc, T4><<<gpu.getSMPCount(), large_block_size>>>(*cgw);
    break;
  case CellGridAction::XFER_FORCES:
    problem = true;
    break;
  }
  if (problem) {
    rtErr("A coordinate synthesis abstract must be provided in order to execute " +
          getEnumerationName(process) + " with a CellGrid object.  Const-ness of either object "
          "may also be important.", "executCellGridAction");
  }
}

template <typename T, typename Tacc, typename Tcalc, typename T4>
void executeCellGridAction(CellGridWriter<T, Tacc, Tcalc, T4> *cgw,
                           const PsSynthesisReader &poly_psr, const GpuDetails &gpu,
                           const CellGridAction process) {
  bool problem = false;
  switch (process) {
  case CellGridAction::INIT_FORCES:
    kInitializeForces<T, Tacc, Tcalc, T4><<<gpu.getSMPCount(), large_block_size>>>(*cgw);
    break;
  case CellGridAction::XFER_FORCES:
    problem = true;
    break;
  }
  if (problem) {
    rtErr("A writeable coordinate synthesis abstract must be provided in order to execute " +
          getEnumerationName(process) + ".", "executeCellGridAction");
  }
}

template <typename T, typename Tacc, typename Tcalc, typename T4>
void executeCellGridAction(const CellGridReader<T, Tacc, Tcalc, T4> &cgr,
                           PsSynthesisWriter *poly_psw, const GpuDetails &gpu,
                           const CellGridAction process) {
  bool problem = false;
  switch (process) {
  case CellGridAction::INIT_FORCES:
    problem = true;
    break;
  case CellGridAction::XFER_FORCES:
    {
      const size_t ct_acc = std::type_index(typeid(Tacc)).hash_code();
      kContributeForces<T, Tacc,
                        Tcalc, T4><<<gpu.getSMPCount(), large_block_size>>>(cgr, ct_acc,
                                                                            int_type_index,
                                                                            *poly_psw);
    }
    break;
  }
  if (problem) {
    rtErr("A writeable cell grid stract must be provided in order to execute " +
          getEnumerationName(process) + ".", "executeCellGridAction");
  }
}
/// \}

/// \brief Unroll the matrix data and coordinate data types for the CellGridWriter abstract to
///        eventually launch the appropriate templated initialization kernel.  The matrix and
///        coordinate representations are coupled to one another, meaning that providing one will
///        imply the identity of the other for the purposes of type restoration.  Overloading and
///        descriptions of parameters follow executeCellGridAction() above.
/// \{
template <typename T, typename Tcalc, typename T4>
void unrollLaunchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_acc,
                                const GpuDetails &gpu, const CellGridAction process) {
  if (tc_acc == int_type_index) {
    CellGridWriter<T, int, Tcalc, T4> rcgw = restoreType<T, int, Tcalc, T4>(cgw);
    executeCellGridAction(&rcgw, gpu, process);
  }
  else if (tc_acc == llint_type_index) {
    CellGridWriter<T, llint, Tcalc, T4> rcgw = restoreType<T, llint, Tcalc, T4>(cgw);
    executeCellGridAction(&rcgw, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "unrollLaunchCellGridAction");
  }
}

template <typename T, typename Tcalc, typename T4>
void unrollLaunchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_acc,
                                const PsSynthesisReader &poly_psr, const GpuDetails &gpu,
                                const CellGridAction process) {
  if (tc_acc == int_type_index) {
    CellGridWriter<T, int, Tcalc, T4> rcgw = restoreType<T, int, Tcalc, T4>(cgw);
    executeCellGridAction(&rcgw, poly_psr, gpu, process);
  }
  else if (tc_acc == llint_type_index) {
    CellGridWriter<T, llint, Tcalc, T4> rcgw = restoreType<T, llint, Tcalc, T4>(cgw);
    executeCellGridAction(&rcgw, poly_psr, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "unrollLaunchCellGridAction");
  }
}

template <typename T, typename Tcalc, typename T4>
void unrollLaunchCellGridAction(const CellGridReader<void, void, void, void> &cgr,
                                const size_t tc_acc, PsSynthesisWriter *poly_psw,
                                const GpuDetails &gpu, const CellGridAction process) {
  if (tc_acc == int_type_index) {
    const CellGridReader<T, int, Tcalc , T4> rcgr = restoreType<T, int, Tcalc, T4>(cgr);
    executeCellGridAction(rcgr, poly_psw, gpu, process);
  }
  else if (tc_acc == llint_type_index) {
    const CellGridReader<T, llint, Tcalc, T4> rcgr = restoreType<T, llint, Tcalc, T4>(cgr);
    executeCellGridAction(rcgr, poly_psw, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "unrollLaunchCellGridAction");
  }
}
/// \}

} // namespace energy
} // namespace stormm
