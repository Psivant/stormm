// -*-c++-*-
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "cellgrid.h"

namespace stormm {
namespace energy {

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
      cgw.yfrc[image_idx] = value_zero;
      cgw.xfrc_ovrf[image_idx] = 0;
      cgw.yfrc_ovrf[image_idx] = 0;
      cgw.zfrc_ovrf[image_idx] = 0;
    }
    warp_pos += warps_per_kernel;
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
template <typename T, typename Tacc, typename Tcalc, typename Tcrd>
void executeCellGridAction(CellGridWriter<T, Tacc, Tcalc, Tcrd> *cgw, const GpuDetails &gpu,
                           const CellGridAction process) {
  bool problem = false;
  switch (process) {
  case CellGridAction::INIT_FORCES:
    kInitializeForces<T, Tacc, Tcalc, Tcrd><<<gpu.getSMPCount(), large_block_size>>>(*cgw);
    break;
  case CellGridAction::XFER_FORCES:
    problem = true;
    break;
  case CellGridAction::UPDATE_IMG_COORD:
    break;
  case CellGridAction::UPDATE_IMG_CELLS:
    break;
  }
  if (problem) {
    rtErr("A coordinate synthesis abstract must be provided in order to execute " +
          getEnumerationName(process) + " with a CellGrid object.  Const-ness of either object "
          "may also be important.", "executCellGridAction");
  }
}

template <typename T, typename Tacc, typename Tcalc, typename Tcrd>
void executeCellGridAction(CellGridWriter<T, Tacc, Tcalc, Tcrd> *cgw,
                           const PsSynthesisReader &poly_psr, const GpuDetails &gpu,
                           const CellGridAction process) {
  bool problem = false;
  switch (process) {
  case CellGridAction::INIT_FORCES:
    kInitializeForces<T, Tacc, Tcalc, Tcrd><<<gpu.getSMPCount(), large_block_size>>>(*cgw);
    break;
  case CellGridAction::XFER_FORCES:
    problem = true;
    break;
  case CellGridAction::UPDATE_IMG_COORD:
    break;
  case CellGridAction::UPDATE_IMG_CELLS:
    break;
  }
  if (problem) {
    rtErr("A writeable coordinate synthesis abstract must be provided in order to execute " +
          getEnumerationName(process) + ".", "executeCellGridAction");
  }
}

template <typename T, typename Tacc, typename Tcalc, typename Tcrd>
void executeCellGridAction(const CellGridReader<T, Tacc, Tcalc, Tcrd> &cgw,
                           PsSynthesisWriter *poly_psr, const GpuDetails &gpu,
                           const CellGridAction process) {
  bool problem = false;
  switch (process) {
  case CellGridAction::INIT_FORCES:
    problem = true;
    break;
  case CellGridAction::XFER_FORCES:
    break;
  case CellGridAction::UPDATE_IMG_COORD:
    break;
  case CellGridAction::UPDATE_IMG_CELLS:
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
template <typename T, typename Tcalc, typename Tcrd>
void unrollLaunchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_acc,
                                const GpuDetails &gpu, const CellGridAction process) {
  if (tc_acc == int_type_index) {
    CellGridWriter<T, int, Tcalc, Tcrd> rcgw = restoreType<T, int, Tcalc, Tcrd>(cgw);
    executeCellGridAction(&rcgw, gpu, process);
  }
  else if (tc_acc == llint_type_index) {
    CellGridWriter<T, llint, Tcalc, Tcrd> rcgw = restoreType<T, llint, Tcalc, Tcrd>(cgw);
    executeCellGridAction(&rcgw, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "unrollLaunchCellGridAction");
  }
}

template <typename T, typename Tcalc, typename Tcrd>
void unrollLaunchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_acc,
                                const PsSynthesisReader &poly_psr, const GpuDetails &gpu,
                                const CellGridAction process) {
  if (tc_acc == int_type_index) {
    CellGridWriter<T, int, Tcalc, Tcrd> rcgw = restoreType<T, int, Tcalc, Tcrd>(cgw);
    executeCellGridAction(&rcgw, poly_psr, gpu, process);
  }
  else if (tc_acc == llint_type_index) {
    CellGridWriter<T, llint, Tcalc, Tcrd> rcgw = restoreType<T, llint, Tcalc, Tcrd>(cgw);
    executeCellGridAction(&rcgw, poly_psr, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "unrollLaunchCellGridAction");
  }
}

template <typename T, typename Tcalc, typename Tcrd>
void unrollLaunchCellGridAction(const CellGridReader<void, void, void, void> &cgr,
                                const size_t tc_acc, PsSynthesisWriter *poly_psw,
                                const GpuDetails &gpu, const CellGridAction process) {
  if (tc_acc == int_type_index) {
    const CellGridReader<T, int, Tcalc , Tcrd> rcgr = restoreType<T, int, Tcalc, Tcrd>(cgr);
    executeCellGridAction(rcgr, poly_psw, gpu, process);
  }
  else if (tc_acc == llint_type_index) {
    const CellGridReader<T, llint, Tcalc, Tcrd> rcgr = restoreType<T, llint, Tcalc, Tcrd>(cgr);
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
