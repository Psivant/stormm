// -*-c++-*-
#include "copyright.h"
#include "hpc_cellgrid.cuh"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
  void launchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_mat,
                          const size_t tc_acc, const GpuDetails &gpu,
                          const CellGridAction process) {

  // Restoring the type of the CellGrid object carries many possibilities, although only the type
  // of the accumulators is essential to their initialization.
  if (tc_mat == double_type_index) {
    unrollLaunchCellGridAction<double, double, double4>(cgw, tc_acc, gpu, process);
  }
  else if (tc_mat == float_type_index) {
    unrollLaunchCellGridAction<float, float, float4>(cgw, tc_acc, gpu, process);
  }
  else if (tc_mat == int_type_index) {
    unrollLaunchCellGridAction<int, float, int4>(cgw, tc_acc, gpu, process);
  }
  else if (tc_mat == llint_type_index) {
    unrollLaunchCellGridAction<llint, double, llint4>(cgw, tc_acc, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "launchCellGridAction");
  }
}

//-------------------------------------------------------------------------------------------------
void launchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_mat,
                          const size_t tc_acc, const PsSynthesisReader &poly_psr,
                          const GpuDetails &gpu, const CellGridAction process) {

  // Restoring the type of the CellGrid object carries many possibilities, although only the type
  // of the accumulators is essential to their initialization.
  if (tc_mat == double_type_index) {
    unrollLaunchCellGridAction<double, double, double4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else if (tc_mat == float_type_index) {
    unrollLaunchCellGridAction<float, float, float4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else if (tc_mat == int_type_index) {
    unrollLaunchCellGridAction<int, float, int4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else if (tc_mat == llint_type_index) {
    unrollLaunchCellGridAction<llint, double, llint4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "launchCellGridAction");
  }
}

//-------------------------------------------------------------------------------------------------
void launchCellGridAction(const CellGridReader<void, void, void, void> &cgr, const size_t tc_mat,
                          const size_t tc_acc, PsSynthesisWriter *poly_psw, const GpuDetails &gpu,
                          const CellGridAction process) {

  // Restoring the type of the CellGrid object carries many possibilities, although only the type
  // of the accumulators is essential to their initialization.
  if (tc_mat == double_type_index) {
    unrollLaunchCellGridAction<double, double, double4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else if (tc_mat == float_type_index) {
    unrollLaunchCellGridAction<float, float, float4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else if (tc_mat == int_type_index) {
    unrollLaunchCellGridAction<int, float, int4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else if (tc_mat == llint_type_index) {
    unrollLaunchCellGridAction<llint, double, llint4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "launchCellGridAction");
  }
}

} // namespace energy
} // namespace stormm
