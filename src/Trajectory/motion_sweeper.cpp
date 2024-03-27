#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "motion_sweeper.h"

namespace stormm {
namespace trajectory {

using stmath::invertSquareMatrix;
using stmath::roundUp;
using synthesis::PsSynthesisReader;
  
//-------------------------------------------------------------------------------------------------
MotionSweepWriter::MotionSweepWriter(const int nwu_in, const int4* work_units_in,
                                     const double* total_mass_in, const double com_scale_in,
                                     const double mv_scale_in, const double inrt_scale_in,
                                     llint* xcom_in, llint* ycom_in, llint* zcom_in,
                                     llint* xcom_nxt_in, llint* ycom_nxt_in, llint* zcom_nxt_in,
                                     int* xcom_ovrf_in, int* ycom_ovrf_in, int* zcom_ovrf_in,
                                     int* xcom_nxt_ovrf_in, int* ycom_nxt_ovrf_in,
                                     int* zcom_nxt_ovrf_in, llint* xmv_in, llint* ymv_in,
                                     llint* zmv_in, llint* xmv_nxt_in, llint* ymv_nxt_in,
                                     llint* zmv_nxt_in, int* xmv_ovrf_in, int* ymv_ovrf_in,
                                     int* zmv_ovrf_in, int* xmv_nxt_ovrf_in, int* ymv_nxt_ovrf_in,
                                     int* zmv_nxt_ovrf_in, llint* rxmv_in, llint* rymv_in,
                                     llint* rzmv_in, llint* rxmv_nxt_in, llint* rymv_nxt_in,
                                     llint* rzmv_nxt_in, int* rxmv_ovrf_in, int* rymv_ovrf_in,
                                     int* rzmv_ovrf_in, int* rxmv_nxt_ovrf_in,
                                     int* rymv_nxt_ovrf_in, int* rzmv_nxt_ovrf_in, llint* inrt_in,
                                     llint* inrt_nxt_in, int* inrt_ovrf_in,
                                     int* inrt_nxt_ovrf_in) :
    nwu{nwu_in}, work_units{work_units_in}, total_mass{total_mass_in}, com_scale{com_scale_in},
    mv_scale{mv_scale_in}, inrt_scale{inrt_scale_in}, xcom{xcom_in}, ycom{ycom_in}, zcom{zcom_in},
    xcom_nxt{xcom_nxt_in}, ycom_nxt{ycom_nxt_in}, zcom_nxt{zcom_nxt_in}, xcom_ovrf{xcom_ovrf_in},
    ycom_ovrf{ycom_ovrf_in}, zcom_ovrf{zcom_ovrf_in}, xcom_nxt_ovrf{xcom_nxt_ovrf_in},
    ycom_nxt_ovrf{ycom_nxt_ovrf_in}, zcom_nxt_ovrf{zcom_nxt_ovrf_in}, xmv{xmv_in}, ymv{ymv_in},
    zmv{zmv_in}, xmv_nxt{xmv_nxt_in}, ymv_nxt{ymv_nxt_in}, zmv_nxt{zmv_nxt_in},
    xmv_ovrf{xmv_ovrf_in}, ymv_ovrf{ymv_ovrf_in}, zmv_ovrf{zmv_ovrf_in},
    xmv_nxt_ovrf{xmv_nxt_ovrf_in}, ymv_nxt_ovrf{ymv_nxt_ovrf_in}, zmv_nxt_ovrf{zmv_nxt_ovrf_in},
    rxmv{rxmv_in}, rymv{rymv_in}, rzmv{rzmv_in}, rxmv_nxt{rxmv_nxt_in}, rymv_nxt{rymv_nxt_in},
    rzmv_nxt{rzmv_nxt_in}, rxmv_ovrf{rxmv_ovrf_in}, rymv_ovrf{rymv_ovrf_in},
    rzmv_ovrf{rzmv_ovrf_in}, rxmv_nxt_ovrf{rxmv_nxt_ovrf_in}, rymv_nxt_ovrf{rymv_nxt_ovrf_in},
    rzmv_nxt_ovrf{rzmv_nxt_ovrf_in}, inrt{inrt_in}, inrt_nxt{inrt_nxt_in}, inrt_ovrf{inrt_ovrf_in},
    inrt_nxt_ovrf{inrt_nxt_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
MotionSweepReader::MotionSweepReader(const int nwu_in, const int4* work_units_in,
                                     const double* total_mass_in, const double com_scale_in,
                                     const double mv_scale_in, const double inrt_scale_in,
                                     const llint* xcom_in, const llint* ycom_in,
                                     const llint* zcom_in, const llint* xcom_nxt_in,
                                     const llint* ycom_nxt_in, const llint* zcom_nxt_in,
                                     const int* xcom_ovrf_in, const int* ycom_ovrf_in,
                                     const int* zcom_ovrf_in, const int* xcom_nxt_ovrf_in,
                                     const int* ycom_nxt_ovrf_in, const int* zcom_nxt_ovrf_in,
                                     const llint* xmv_in, const llint* ymv_in, const llint* zmv_in,
                                     const llint* xmv_nxt_in, const llint* ymv_nxt_in,
                                     const llint* zmv_nxt_in, const int* xmv_ovrf_in,
                                     const int* ymv_ovrf_in, const int* zmv_ovrf_in,
                                     const int* xmv_nxt_ovrf_in, const int* ymv_nxt_ovrf_in,
                                     const int* zmv_nxt_ovrf_in, const llint* rxmv_in,
                                     const llint* rymv_in, const llint* rzmv_in,
                                     const llint* rxmv_nxt_in, const llint* rymv_nxt_in,
                                     const llint* rzmv_nxt_in, const int* rxmv_ovrf_in,
                                     const int* rymv_ovrf_in, const int* rzmv_ovrf_in,
                                     const int* rxmv_nxt_ovrf_in, const int* rymv_nxt_ovrf_in,
                                     const int* rzmv_nxt_ovrf_in, const llint* inrt_in,
                                     const llint* inrt_nxt_in, const int* inrt_ovrf_in,
                                     const int* inrt_nxt_ovrf_in) :
    nwu{nwu_in}, work_units{work_units_in}, total_mass{total_mass_in}, com_scale{com_scale_in},
    mv_scale{mv_scale_in}, inrt_scale{inrt_scale_in}, xcom{xcom_in}, ycom{ycom_in}, zcom{zcom_in},
    xcom_nxt{xcom_nxt_in}, ycom_nxt{ycom_nxt_in}, zcom_nxt{zcom_nxt_in}, xcom_ovrf{xcom_ovrf_in},
    ycom_ovrf{ycom_ovrf_in}, zcom_ovrf{zcom_ovrf_in}, xcom_nxt_ovrf{xcom_nxt_ovrf_in},
    ycom_nxt_ovrf{ycom_nxt_ovrf_in}, zcom_nxt_ovrf{zcom_nxt_ovrf_in}, xmv{xmv_in}, ymv{ymv_in},
    zmv{zmv_in}, xmv_nxt{xmv_nxt_in}, ymv_nxt{ymv_nxt_in}, zmv_nxt{zmv_nxt_in},
    xmv_ovrf{xmv_ovrf_in}, ymv_ovrf{ymv_ovrf_in}, zmv_ovrf{zmv_ovrf_in},
    xmv_nxt_ovrf{xmv_nxt_ovrf_in}, ymv_nxt_ovrf{ymv_nxt_ovrf_in}, zmv_nxt_ovrf{zmv_nxt_ovrf_in},
    rxmv{rxmv_in}, rymv{rymv_in}, rzmv{rzmv_in}, rxmv_nxt{rxmv_nxt_in}, rymv_nxt{rymv_nxt_in},
    rzmv_nxt{rzmv_nxt_in}, rxmv_ovrf{rxmv_ovrf_in}, rymv_ovrf{rymv_ovrf_in},
    rzmv_ovrf{rzmv_ovrf_in}, rxmv_nxt_ovrf{rxmv_nxt_ovrf_in}, rymv_nxt_ovrf{rymv_nxt_ovrf_in},
    rzmv_nxt_ovrf{rzmv_nxt_ovrf_in}, inrt{inrt_in}, inrt_nxt{inrt_nxt_in}, inrt_ovrf{inrt_ovrf_in},
    inrt_nxt_ovrf{inrt_nxt_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
MotionSweepReader::MotionSweepReader(const MotionSweepWriter *w) :
    nwu{w->nwu}, work_units{w->work_units}, com_scale{w->com_scale}, mv_scale{w->mv_scale},
    inrt_scale{w->inrt_scale}, xcom{w->xcom}, ycom{w->ycom}, zcom{w->zcom}, xcom_nxt{w->xcom_nxt},
    ycom_nxt{w->ycom_nxt}, zcom_nxt{w->zcom_nxt}, xcom_ovrf{w->xcom_ovrf}, ycom_ovrf{w->ycom_ovrf},
    zcom_ovrf{w->zcom_ovrf}, xcom_nxt_ovrf{w->xcom_nxt_ovrf}, ycom_nxt_ovrf{w->ycom_nxt_ovrf},
    zcom_nxt_ovrf{w->zcom_nxt_ovrf}, total_mass{w->total_mass}, xmv{w->xmv}, ymv{w->ymv},
    zmv{w->zmv}, xmv_nxt{w->xmv_nxt}, ymv_nxt{w->ymv_nxt}, zmv_nxt{w->zmv_nxt},
    xmv_ovrf{w->xmv_ovrf}, ymv_ovrf{w->ymv_ovrf}, zmv_ovrf{w->zmv_ovrf},
    xmv_nxt_ovrf{w->xmv_nxt_ovrf}, ymv_nxt_ovrf{w->ymv_nxt_ovrf}, zmv_nxt_ovrf{w->zmv_nxt_ovrf},
    rxmv{w->rxmv}, rymv{w->rymv}, rzmv{w->rzmv}, rxmv_nxt{w->rxmv_nxt}, rymv_nxt{w->rymv_nxt},
    rzmv_nxt{w->rzmv_nxt}, rxmv_ovrf{w->rxmv_ovrf}, rymv_ovrf{w->rymv_ovrf},
    rzmv_ovrf{w->rzmv_ovrf}, rxmv_nxt_ovrf{w->rxmv_nxt_ovrf}, rymv_nxt_ovrf{w->rymv_nxt_ovrf},
    rzmv_nxt_ovrf{w->rzmv_nxt_ovrf}, inrt{w->inrt}, inrt_nxt{w->inrt_nxt}, inrt_ovrf{w->inrt_ovrf},
    inrt_nxt_ovrf{w->inrt_nxt_ovrf}
{}


//-------------------------------------------------------------------------------------------------
MotionSweepReader::MotionSweepReader(const MotionSweepWriter &w) :
    nwu{w.nwu}, work_units{w.work_units}, com_scale{w.com_scale}, mv_scale{w.mv_scale},
    inrt_scale{w.inrt_scale}, xcom{w.xcom}, ycom{w.ycom}, zcom{w.zcom}, xcom_nxt{w.xcom_nxt},
    ycom_nxt{w.ycom_nxt}, zcom_nxt{w.zcom_nxt}, xcom_ovrf{w.xcom_ovrf}, ycom_ovrf{w.ycom_ovrf},
    zcom_ovrf{w.zcom_ovrf}, xcom_nxt_ovrf{w.xcom_nxt_ovrf}, ycom_nxt_ovrf{w.ycom_nxt_ovrf},
    zcom_nxt_ovrf{w.zcom_nxt_ovrf}, total_mass{w.total_mass}, xmv{w.xmv}, ymv{w.ymv},
    zmv{w.zmv}, xmv_nxt{w.xmv_nxt}, ymv_nxt{w.ymv_nxt}, zmv_nxt{w.zmv_nxt}, xmv_ovrf{w.xmv_ovrf},
    ymv_ovrf{w.ymv_ovrf}, zmv_ovrf{w.zmv_ovrf}, xmv_nxt_ovrf{w.xmv_nxt_ovrf},
    ymv_nxt_ovrf{w.ymv_nxt_ovrf}, zmv_nxt_ovrf{w.zmv_nxt_ovrf}, rxmv{w.rxmv}, rymv{w.rymv},
    rzmv{w.rzmv}, rxmv_nxt{w.rxmv_nxt}, rymv_nxt{w.rymv_nxt}, rzmv_nxt{w.rzmv_nxt},
    rxmv_ovrf{w.rxmv_ovrf}, rymv_ovrf{w.rymv_ovrf}, rzmv_ovrf{w.rzmv_ovrf},
    rxmv_nxt_ovrf{w.rxmv_nxt_ovrf}, rymv_nxt_ovrf{w.rymv_nxt_ovrf}, rzmv_nxt_ovrf{w.rzmv_nxt_ovrf},
    inrt{w.inrt}, inrt_nxt{w.inrt_nxt}, inrt_ovrf{w.inrt_ovrf}, inrt_nxt_ovrf{w.inrt_nxt_ovrf}
{}

//-------------------------------------------------------------------------------------------------
MotionSweeper::MotionSweeper(const PhaseSpaceSynthesis *poly_ps,
                             const int momentum_bit_count_in,
                             const int center_of_mass_bit_count_in,
                             const int inertia_bit_count_in) :
    cycle_position{CoordinateCycle::WHITE},
    center_of_mass_bit_count{center_of_mass_bit_count_in},
    momentum_bit_count{momentum_bit_count_in},
    inertia_bit_count{inertia_bit_count_in},
    xcom_white{HybridKind::POINTER, "mosw_xcom_wh"},
    ycom_white{HybridKind::POINTER, "mosw_ycom_wh"},
    zcom_white{HybridKind::POINTER, "mosw_zcom_wh"},
    xcom_black{HybridKind::POINTER, "mosw_xcom_bk"},
    ycom_black{HybridKind::POINTER, "mosw_ycom_bk"},
    zcom_black{HybridKind::POINTER, "mosw_zcom_bk"},
    xcom_overflow_white{HybridKind::POINTER, "mosw_ovrf_xcom_wh"},
    ycom_overflow_white{HybridKind::POINTER, "mosw_ovrf_ycom_wh"},
    zcom_overflow_white{HybridKind::POINTER, "mosw_ovrf_zcom_wh"},
    xcom_overflow_black{HybridKind::POINTER, "mosw_ovrf_xcom_bk"},
    ycom_overflow_black{HybridKind::POINTER, "mosw_ovrf_ycom_bk"},
    zcom_overflow_black{HybridKind::POINTER, "mosw_ovrf_zcom_bk"},
    total_mass{HybridKind::ARRAY},
    xmv_white{HybridKind::POINTER, "mosw_xmv_wh"},
    ymv_white{HybridKind::POINTER, "mosw_ymv_wh"},
    zmv_white{HybridKind::POINTER, "mosw_zmv_wh"},
    xmv_black{HybridKind::POINTER, "mosw_xmv_bk"},
    ymv_black{HybridKind::POINTER, "mosw_ymv_bk"},
    zmv_black{HybridKind::POINTER, "mosw_zmv_bk"},
    xmv_overflow_white{HybridKind::POINTER, "mosw_xmv_ovrf_wh"},
    ymv_overflow_white{HybridKind::POINTER, "mosw_ymv_ovrf_wh"},
    zmv_overflow_white{HybridKind::POINTER, "mosw_zmv_ovrf_wh"},
    xmv_overflow_black{HybridKind::POINTER, "mosw_xmv_ovrf_bk"},
    ymv_overflow_black{HybridKind::POINTER, "mosw_ymv_ovrf_bk"},
    zmv_overflow_black{HybridKind::POINTER, "mosw_zmv_ovrf_bk"},
    rxmv_white{HybridKind::POINTER, "mosw_rxmv_wh"},
    rymv_white{HybridKind::POINTER, "mosw_rymv_wh"},
    rzmv_white{HybridKind::POINTER, "mosw_rzmv_wh"},
    rxmv_black{HybridKind::POINTER, "mosw_rxmv_bk"},
    rymv_black{HybridKind::POINTER, "mosw_rymv_bk"},
    rzmv_black{HybridKind::POINTER, "mosw_rzmv_bk"},
    rxmv_overflow_white{HybridKind::POINTER, "mosw_rxmv_ovrf_wh"},
    rymv_overflow_white{HybridKind::POINTER, "mosw_rymv_ovrf_wh"},
    rzmv_overflow_white{HybridKind::POINTER, "mosw_rzmv_ovrf_wh"},
    rxmv_overflow_black{HybridKind::POINTER, "mosw_rxmv_ovrf_bk"},
    rymv_overflow_black{HybridKind::POINTER, "mosw_rymv_ovrf_bk"},
    rzmv_overflow_black{HybridKind::POINTER, "mosw_rzmv_ovrf_bk"},
    inertial_tensor_white{HybridKind::POINTER, "mosw_inrt_wh"},
    inertial_tensor_black{HybridKind::POINTER, "mosw_inrt_bk"},
    inertial_tensor_overflow_white{HybridKind::POINTER, "mosw_inrt_ovrf_wh"},
    inertial_tensor_overflow_black{HybridKind::POINTER, "mosw_inrt_ovrf_bk"},
    llint_data{HybridKind::ARRAY, "mosw_llint_data"},
    int_data{HybridKind::ARRAY, "mosw_int_data"},
    work_unit_count{0},
    work_units{HybridKind::ARRAY, "mosw_work_units"},
    poly_ps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps)}
{
  allocate();
  
  // Load the total masses for all systems
  const PsSynthesisReader poly_psr = poly_ps_ptr->data();
  const std::vector<AtomGraph*>& unique_ag = poly_ps_ptr->getUniqueTopologies();
  const int n_unique_ag = poly_ps_ptr->getUniqueTopologyCount();
  std::vector<double> unique_masses(n_unique_ag);
  for (int i = 0; i < n_unique_ag; i++) {
    unique_masses[i] = unique_ag[i]->getTotalMass();   
  }
  total_mass.resize(poly_psr.system_count);
  double* tmass_ptr = total_mass.data();
  for (int i = 0; i < poly_psr.system_count; i++) {
    tmass_ptr[i] = unique_masses[poly_ps_ptr->getUniqueTopologyIndex(i)];
  }
  
  // Determine the work units
  int nbatch = 0;
  for (int i = 0; i < poly_psr.system_count; i++) {
    nbatch += (poly_psr.atom_counts[i] + small_block_size - 1) / small_block_size;
  }
  work_units.resize(nbatch);
  work_unit_count = nbatch;
  int wu_idx = 0;
  int4* wu_ptr = work_units.data();
  for (int i = 0; i < poly_psr.system_count; i++) {
    int system_first_wu = 1;
    for (int j = 0; j < poly_psr.atom_counts[i]; j += small_block_size) {
      const int natom_wu = (j + small_block_size < poly_psr.atom_counts[i]) ?
                           small_block_size : poly_psr.atom_counts[i] - j;
      const int start_idx = poly_psr.atom_starts[i] + j;
      const int4 twu = { start_idx, start_idx + natom_wu, i, system_first_wu };
      wu_ptr[wu_idx] = twu;
      system_first_wu = 0;
      wu_idx++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
MotionSweeper::MotionSweeper(const MotionSweeper &original) :
    cycle_position{original.cycle_position},
    center_of_mass_bit_count{original.center_of_mass_bit_count},
    momentum_bit_count{original.momentum_bit_count},
    inertia_bit_count{original.inertia_bit_count},
    xcom_white{original.xcom_white},
    ycom_white{original.ycom_white},
    zcom_white{original.zcom_white},
    xcom_black{original.xcom_black},
    ycom_black{original.ycom_black},
    zcom_black{original.zcom_black},
    xcom_overflow_white{original.xcom_overflow_white},
    ycom_overflow_white{original.ycom_overflow_white},
    zcom_overflow_white{original.zcom_overflow_white},
    xcom_overflow_black{original.xcom_overflow_black},
    ycom_overflow_black{original.ycom_overflow_black},
    zcom_overflow_black{original.zcom_overflow_black},
    total_mass{original.total_mass},
    xmv_white{original.xmv_white},
    ymv_white{original.ymv_white},
    zmv_white{original.zmv_white},
    xmv_black{original.xmv_black},
    ymv_black{original.ymv_black},
    zmv_black{original.zmv_black},
    xmv_overflow_white{original.xmv_overflow_white},
    ymv_overflow_white{original.ymv_overflow_white},
    zmv_overflow_white{original.zmv_overflow_white},
    xmv_overflow_black{original.xmv_overflow_black},
    ymv_overflow_black{original.ymv_overflow_black},
    zmv_overflow_black{original.zmv_overflow_black},
    rxmv_white{original.rxmv_white},
    rymv_white{original.rymv_white},
    rzmv_white{original.rzmv_white},
    rxmv_black{original.rxmv_black},
    rymv_black{original.rymv_black},
    rzmv_black{original.rzmv_black},
    rxmv_overflow_white{original.rxmv_overflow_white},
    rymv_overflow_white{original.rymv_overflow_white},
    rzmv_overflow_white{original.rzmv_overflow_white},
    rxmv_overflow_black{original.rxmv_overflow_black},
    rymv_overflow_black{original.rymv_overflow_black},
    rzmv_overflow_black{original.rzmv_overflow_black},
    inertial_tensor_white{original.inertial_tensor_white},
    inertial_tensor_black{original.inertial_tensor_black},
    inertial_tensor_overflow_white{original.inertial_tensor_overflow_white},
    inertial_tensor_overflow_black{original.inertial_tensor_overflow_black},
    llint_data{original.llint_data},
    int_data{original.int_data},
    work_unit_count{original.work_unit_count},
    work_units{original.work_units},
    poly_ps_ptr{original.poly_ps_ptr}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
MotionSweeper::MotionSweeper(MotionSweeper &&original) :
    cycle_position{original.cycle_position},
    center_of_mass_bit_count{original.center_of_mass_bit_count},
    momentum_bit_count{original.momentum_bit_count},
    inertia_bit_count{original.inertia_bit_count},
    xcom_white{std::move(original.xcom_white)},
    ycom_white{std::move(original.ycom_white)},
    zcom_white{std::move(original.zcom_white)},
    xcom_black{std::move(original.xcom_black)},
    ycom_black{std::move(original.ycom_black)},
    zcom_black{std::move(original.zcom_black)},
    xcom_overflow_white{std::move(original.xcom_overflow_white)},
    ycom_overflow_white{std::move(original.ycom_overflow_white)},
    zcom_overflow_white{std::move(original.zcom_overflow_white)},
    xcom_overflow_black{std::move(original.xcom_overflow_black)},
    ycom_overflow_black{std::move(original.ycom_overflow_black)},
    zcom_overflow_black{std::move(original.zcom_overflow_black)},
    total_mass{std::move(original.total_mass)},
    xmv_white{std::move(original.xmv_white)},
    ymv_white{std::move(original.ymv_white)},
    zmv_white{std::move(original.zmv_white)},
    xmv_black{std::move(original.xmv_black)},
    ymv_black{std::move(original.ymv_black)},
    zmv_black{std::move(original.zmv_black)},
    xmv_overflow_white{std::move(original.xmv_overflow_white)},
    ymv_overflow_white{std::move(original.ymv_overflow_white)},
    zmv_overflow_white{std::move(original.zmv_overflow_white)},
    xmv_overflow_black{std::move(original.xmv_overflow_black)},
    ymv_overflow_black{std::move(original.ymv_overflow_black)},
    zmv_overflow_black{std::move(original.zmv_overflow_black)},
    rxmv_white{std::move(original.rxmv_white)},
    rymv_white{std::move(original.rymv_white)},
    rzmv_white{std::move(original.rzmv_white)},
    rxmv_black{std::move(original.rxmv_black)},
    rymv_black{std::move(original.rymv_black)},
    rzmv_black{std::move(original.rzmv_black)},
    rxmv_overflow_white{std::move(original.rxmv_overflow_white)},
    rymv_overflow_white{std::move(original.rymv_overflow_white)},
    rzmv_overflow_white{std::move(original.rzmv_overflow_white)},
    rxmv_overflow_black{std::move(original.rxmv_overflow_black)},
    rymv_overflow_black{std::move(original.rymv_overflow_black)},
    rzmv_overflow_black{std::move(original.rzmv_overflow_black)},
    inertial_tensor_white{std::move(original.inertial_tensor_white)},
    inertial_tensor_black{std::move(original.inertial_tensor_black)},
    inertial_tensor_overflow_white{std::move(original.inertial_tensor_overflow_white)},
    inertial_tensor_overflow_black{std::move(original.inertial_tensor_overflow_black)},
    llint_data{std::move(original.llint_data)},
    int_data{std::move(original.int_data)},
    work_unit_count{original.work_unit_count},
    work_units{std::move(original.work_units)},
    poly_ps_ptr{original.poly_ps_ptr}
{}

//-------------------------------------------------------------------------------------------------
MotionSweeper& MotionSweeper::operator=(const MotionSweeper &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy the POINTER-kind Hybrid objects, temporarily leaving them targeted to the other object.
  cycle_position = other.cycle_position;
  center_of_mass_bit_count = other.center_of_mass_bit_count;
  momentum_bit_count = other.momentum_bit_count;
  inertia_bit_count = other.inertia_bit_count;
  xcom_white = other.xcom_white;
  ycom_white = other.ycom_white;
  zcom_white = other.zcom_white;
  xcom_black = other.xcom_black;
  ycom_black = other.ycom_black;
  zcom_black = other.zcom_black;
  xcom_overflow_white = other.xcom_overflow_white;
  ycom_overflow_white = other.ycom_overflow_white;
  zcom_overflow_white = other.zcom_overflow_white;
  xcom_overflow_black = other.xcom_overflow_black;
  ycom_overflow_black = other.ycom_overflow_black;
  zcom_overflow_black = other.zcom_overflow_black;
  total_mass = other.total_mass;
  xmv_white = other.xmv_white;
  ymv_white = other.ymv_white;
  zmv_white = other.zmv_white;
  xmv_black = other.xmv_black;
  ymv_black = other.ymv_black;
  zmv_black = other.zmv_black;
  xmv_overflow_white = other.xmv_overflow_white;
  ymv_overflow_white = other.ymv_overflow_white;
  zmv_overflow_white = other.zmv_overflow_white;
  xmv_overflow_black = other.xmv_overflow_black;
  ymv_overflow_black = other.ymv_overflow_black;
  zmv_overflow_black = other.zmv_overflow_black;
  rxmv_white = other.rxmv_white;
  rymv_white = other.rymv_white;
  rzmv_white = other.rzmv_white;
  rxmv_black = other.rxmv_black;
  rymv_black = other.rymv_black;
  rzmv_black = other.rzmv_black;
  rxmv_overflow_white = other.rxmv_overflow_white;
  rymv_overflow_white = other.rymv_overflow_white;
  rzmv_overflow_white = other.rzmv_overflow_white;
  rxmv_overflow_black = other.rxmv_overflow_black;
  rymv_overflow_black = other.rymv_overflow_black;
  rzmv_overflow_black = other.rzmv_overflow_black;
  inertial_tensor_white = other.inertial_tensor_white;
  inertial_tensor_black = other.inertial_tensor_black;
  inertial_tensor_overflow_white = other.inertial_tensor_overflow_white;
  inertial_tensor_overflow_black = other.inertial_tensor_overflow_black;
  llint_data = other.llint_data;
  int_data = other.int_data;
  work_unit_count = other.work_unit_count;
  work_units = other.work_units;
  poly_ps_ptr = other.poly_ps_ptr;

  // Repair the pointers with the allocation function (the ARRAY-kind memory will not be resized
  // as all such members already have the correct sizes).
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
MotionSweeper& MotionSweeper::operator=(MotionSweeper &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  cycle_position = other.cycle_position;
  center_of_mass_bit_count = other.center_of_mass_bit_count;
  momentum_bit_count = other.momentum_bit_count;
  inertia_bit_count = other.inertia_bit_count;
  xcom_white = std::move(other.xcom_white);
  ycom_white = std::move(other.ycom_white);
  zcom_white = std::move(other.zcom_white);
  xcom_black = std::move(other.xcom_black);
  ycom_black = std::move(other.ycom_black);
  zcom_black = std::move(other.zcom_black);
  xcom_overflow_white = std::move(other.xcom_overflow_white);
  ycom_overflow_white = std::move(other.ycom_overflow_white);
  zcom_overflow_white = std::move(other.zcom_overflow_white);
  xcom_overflow_black = std::move(other.xcom_overflow_black);
  ycom_overflow_black = std::move(other.ycom_overflow_black);
  zcom_overflow_black = std::move(other.zcom_overflow_black);
  total_mass = std::move(other.total_mass);
  xmv_white = std::move(other.xmv_white);
  ymv_white = std::move(other.ymv_white);
  zmv_white = std::move(other.zmv_white);
  xmv_black = std::move(other.xmv_black);
  ymv_black = std::move(other.ymv_black);
  zmv_black = std::move(other.zmv_black);
  xmv_overflow_white = std::move(other.xmv_overflow_white);
  ymv_overflow_white = std::move(other.ymv_overflow_white);
  zmv_overflow_white = std::move(other.zmv_overflow_white);
  xmv_overflow_black = std::move(other.xmv_overflow_black);
  ymv_overflow_black = std::move(other.ymv_overflow_black);
  zmv_overflow_black = std::move(other.zmv_overflow_black);
  rxmv_white = std::move(other.rxmv_white);
  rymv_white = std::move(other.rymv_white);
  rzmv_white = std::move(other.rzmv_white);
  rxmv_black = std::move(other.rxmv_black);
  rymv_black = std::move(other.rymv_black);
  rzmv_black = std::move(other.rzmv_black);
  rxmv_overflow_white = std::move(other.rxmv_overflow_white);
  rymv_overflow_white = std::move(other.rymv_overflow_white);
  rzmv_overflow_white = std::move(other.rzmv_overflow_white);
  rxmv_overflow_black = std::move(other.rxmv_overflow_black);
  rymv_overflow_black = std::move(other.rymv_overflow_black);
  rzmv_overflow_black = std::move(other.rzmv_overflow_black);
  inertial_tensor_white = std::move(other.inertial_tensor_white);
  inertial_tensor_black = std::move(other.inertial_tensor_black);
  inertial_tensor_overflow_white = std::move(other.inertial_tensor_overflow_white);
  inertial_tensor_overflow_black = std::move(other.inertial_tensor_overflow_black);
  llint_data = std::move(other.llint_data);
  int_data = std::move(other.int_data);
  work_unit_count = other.work_unit_count;
  work_units = std::move(other.work_units);
  poly_ps_ptr = other.poly_ps_ptr;
  return *this;
}
  
//-------------------------------------------------------------------------------------------------
MotionSweeper::MotionSweeper(const PhaseSpaceSynthesis &poly_ps,
                             const int momentum_bit_count_in,
                             const int center_of_mass_bit_count_in,
                             const int inertia_bit_count_in) :
  MotionSweeper(poly_ps.getSelfPointer(), momentum_bit_count_in, center_of_mass_bit_count_in,
                inertia_bit_count_in)
{}

//-------------------------------------------------------------------------------------------------
int MotionSweeper::getSystemCount() const {
  return poly_ps_ptr->getSystemCount();
}

//-------------------------------------------------------------------------------------------------
int MotionSweeper::getWorkUnitCount() const {
  return work_unit_count;
}

//-------------------------------------------------------------------------------------------------
CoordinateCycle MotionSweeper::getCyclePosition() const {
  return cycle_position;
}

//-------------------------------------------------------------------------------------------------
int MotionSweeper::getMomentumBitCount() const {
  return momentum_bit_count;
}

//-------------------------------------------------------------------------------------------------
int MotionSweeper::getCenterOfMassBitCount() const {
  return center_of_mass_bit_count;
}

//-------------------------------------------------------------------------------------------------
int MotionSweeper::getInertialTensorBitCount() const {
  return inertia_bit_count;
}

//-------------------------------------------------------------------------------------------------
double3 MotionSweeper::getNetVelocity(const int idx, const HybridTargetLevel tier) const {
  double3 result;
  const double tmass = total_mass.readHost(idx) * pow(2.0, momentum_bit_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      result.x = hostInt95ToDouble(xmv_white.readHost(idx), xmv_overflow_white.readHost(idx));
      result.y = hostInt95ToDouble(ymv_white.readHost(idx), ymv_overflow_white.readHost(idx));
      result.z = hostInt95ToDouble(zmv_white.readHost(idx), zmv_overflow_white.readHost(idx));
      break;
    case CoordinateCycle::BLACK:
      result.x = hostInt95ToDouble(xmv_black.readHost(idx), xmv_overflow_black.readHost(idx));
      result.y = hostInt95ToDouble(ymv_black.readHost(idx), ymv_overflow_black.readHost(idx));
      result.z = hostInt95ToDouble(zmv_black.readHost(idx), zmv_overflow_black.readHost(idx));
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      result.x = hostInt95ToDouble(xmv_white.readDevice(idx), xmv_overflow_white.readDevice(idx));
      result.y = hostInt95ToDouble(ymv_white.readDevice(idx), ymv_overflow_white.readDevice(idx));
      result.z = hostInt95ToDouble(zmv_white.readDevice(idx), zmv_overflow_white.readDevice(idx));
      break;
    case CoordinateCycle::BLACK:
      result.x = hostInt95ToDouble(xmv_black.readDevice(idx), xmv_overflow_black.readDevice(idx));
      result.y = hostInt95ToDouble(ymv_black.readDevice(idx), ymv_overflow_black.readDevice(idx));
      result.z = hostInt95ToDouble(zmv_black.readDevice(idx), zmv_overflow_black.readDevice(idx));
      break;
    }
    break;
#endif
  }
  result.x /= tmass;
  result.y /= tmass;
  result.z /= tmass;
  return result;
}

//-------------------------------------------------------------------------------------------------
double3 MotionSweeper::getCenterOfMass(const int idx, const HybridTargetLevel tier) const {
  double3 result;
  const double tmass = total_mass.readHost(idx) * pow(2.0, center_of_mass_bit_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      result.x = hostInt95ToDouble(xcom_white.readHost(idx), xcom_overflow_white.readHost(idx));
      result.y = hostInt95ToDouble(ycom_white.readHost(idx), ycom_overflow_white.readHost(idx));
      result.z = hostInt95ToDouble(zcom_white.readHost(idx), zcom_overflow_white.readHost(idx));
      break;
    case CoordinateCycle::BLACK:
      result.x = hostInt95ToDouble(xcom_black.readHost(idx), xcom_overflow_black.readHost(idx));
      result.y = hostInt95ToDouble(ycom_black.readHost(idx), ycom_overflow_black.readHost(idx));
      result.z = hostInt95ToDouble(zcom_black.readHost(idx), zcom_overflow_black.readHost(idx));
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      result.x = hostInt95ToDouble(xcom_white.readDevice(idx),
                                   xcom_overflow_white.readDevice(idx));
      result.y = hostInt95ToDouble(ycom_white.readDevice(idx),
                                   ycom_overflow_white.readDevice(idx));
      result.z = hostInt95ToDouble(zcom_white.readDevice(idx),
                                   zcom_overflow_white.readDevice(idx));
      break;
    case CoordinateCycle::BLACK:
      result.x = hostInt95ToDouble(xcom_black.readDevice(idx),
                                   xcom_overflow_black.readDevice(idx));
      result.y = hostInt95ToDouble(ycom_black.readDevice(idx),
                                   ycom_overflow_black.readDevice(idx));
      result.z = hostInt95ToDouble(zcom_black.readDevice(idx),
                                   zcom_overflow_black.readDevice(idx));
      break;
    }
    break;
#endif
  }
  result.x /= tmass;
  result.y /= tmass;
  result.z /= tmass;
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> MotionSweeper::getInertialTensor(const int idx,
                                                     const HybridTargetLevel tier) const {
  std::vector<double> result(9), cmps(6);
  std::vector<llint> prim;
  std::vector<int> ovrf;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      prim = inertial_tensor_white.readHost(6 * idx, 6);
      ovrf = inertial_tensor_overflow_white.readHost(6 * idx, 6);
      break;
    case CoordinateCycle::BLACK:
      prim = inertial_tensor_black.readHost(6 * idx, 6);
      ovrf = inertial_tensor_overflow_black.readHost(6 * idx, 6);
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      prim = inertial_tensor_white.readDevice(6 * idx, 6);
      ovrf = inertial_tensor_overflow_white.readDevice(6 * idx, 6);
      break;
    case CoordinateCycle::BLACK:
      prim = inertial_tensor_black.readDevice(6 * idx, 6);
      ovrf = inertial_tensor_overflow_black.readDevice(6 * idx, 6);
      break;
    }
    break;
#endif
  }
  const double inv_inrt_scale = pow(2.0, -inertia_bit_count);
  for (int i = 0; i < 6; i++) {
    cmps[i] = hostInt95ToDouble(prim[i], ovrf[i]) * inv_inrt_scale;
  }
  result[0] =  cmps[3] + cmps[5];
  result[1] = -cmps[1];
  result[2] = -cmps[2];
  result[3] = -cmps[1];
  result[4] =  cmps[0] + cmps[5];
  result[5] = -cmps[4];
  result[6] = -cmps[2];
  result[7] = -cmps[4];
  result[8] =  cmps[0] + cmps[3];
  return result;
}

//-------------------------------------------------------------------------------------------------
double3 MotionSweeper::getAngularVelocity(const int idx, const HybridTargetLevel tier) const {
  const std::vector<double> itns = getInertialTensor(idx, tier);
  std::vector<double> invi(9);
  invertSquareMatrix(itns, &invi);
  double3 angm;
  const double rmv_fac = pow(2.0, -momentum_bit_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      angm = { rxmv_white.readHost(idx) * rmv_fac, rymv_white.readHost(idx) * rmv_fac,
               rzmv_white.readHost(idx) * rmv_fac };
      break;
    case CoordinateCycle::BLACK:
      angm = { rxmv_black.readHost(idx) * rmv_fac, rymv_black.readHost(idx) * rmv_fac,
               rzmv_black.readHost(idx) * rmv_fac };
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (cycle_position) {
    case CoordinateCycle::WHITE:
      angm = { rxmv_white.readDevice(idx) * rmv_fac, rymv_white.readDevice(idx) * rmv_fac,
               rzmv_white.readDevice(idx) * rmv_fac };
      break;
    case CoordinateCycle::BLACK:
      angm = { rxmv_black.readDevice(idx) * rmv_fac, rymv_black.readDevice(idx) * rmv_fac,
               rzmv_black.readDevice(idx) * rmv_fac };
      break;
    }
    break;
#endif
  }
  return { (invi[0] * angm.x) + (invi[3] * angm.y) + (invi[6] * angm.z),
           (invi[1] * angm.x) + (invi[4] * angm.y) + (invi[7] * angm.z),
           (invi[2] * angm.x) + (invi[5] * angm.y) + (invi[8] * angm.z) };
}

//-------------------------------------------------------------------------------------------------
void MotionSweeper::updateCyclePosition() {
  switch (cycle_position) {
  case CoordinateCycle::WHITE:
    cycle_position = CoordinateCycle::BLACK;
    break;
  case CoordinateCycle::BLACK:
    cycle_position = CoordinateCycle::WHITE;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
MotionSweepWriter MotionSweeper::data(const CoordinateCycle orientation,
                                      const HybridTargetLevel tier) {
  const double com_fac = pow(2.0, center_of_mass_bit_count);
  const double mv_fac = pow(2.0, momentum_bit_count);
  const double inrt_fac = pow(2.0, inertia_bit_count);
  switch (orientation) {
  case CoordinateCycle::WHITE:
    return MotionSweepWriter(work_unit_count, work_units.data(tier), total_mass.data(tier),
                             com_fac, mv_fac, inrt_fac, xcom_white.data(tier),
                             ycom_white.data(tier), zcom_white.data(tier), xcom_black.data(tier),
                             ycom_black.data(tier), zcom_black.data(tier),
                             xcom_overflow_white.data(tier), ycom_overflow_white.data(tier),
                             zcom_overflow_white.data(tier), xcom_overflow_black.data(tier),
                             ycom_overflow_black.data(tier), zcom_overflow_black.data(tier),
                             xmv_white.data(tier), ymv_white.data(tier), zmv_white.data(tier),
                             xmv_black.data(tier), ymv_black.data(tier), zmv_black.data(tier),
                             xmv_overflow_white.data(tier), ymv_overflow_white.data(tier),
                             zmv_overflow_white.data(tier), xmv_overflow_black.data(tier),
                             ymv_overflow_black.data(tier), zmv_overflow_black.data(tier),
                             rxmv_white.data(tier), rymv_white.data(tier), rzmv_white.data(tier),
                             rxmv_black.data(tier), rymv_black.data(tier), rzmv_black.data(tier),
                             rxmv_overflow_white.data(tier), rymv_overflow_white.data(tier),
                             rzmv_overflow_white.data(tier), rxmv_overflow_black.data(tier),
                             rymv_overflow_black.data(tier), rzmv_overflow_black.data(tier),
                             inertial_tensor_white.data(tier), inertial_tensor_black.data(tier),
                             inertial_tensor_overflow_white.data(tier),
                             inertial_tensor_overflow_black.data(tier));
  case CoordinateCycle::BLACK:
    return MotionSweepWriter(work_unit_count, work_units.data(tier), total_mass.data(tier),
                             com_fac, mv_fac, inrt_fac, xcom_black.data(tier),
                             ycom_black.data(tier), zcom_black.data(tier), xcom_white.data(tier),
                             ycom_white.data(tier), zcom_white.data(tier),
                             xcom_overflow_black.data(tier), ycom_overflow_black.data(tier),
                             zcom_overflow_black.data(tier), xcom_overflow_white.data(tier),
                             ycom_overflow_white.data(tier), zcom_overflow_white.data(tier),
                             xmv_black.data(tier), ymv_black.data(tier), zmv_black.data(tier),
                             xmv_white.data(tier), ymv_white.data(tier), zmv_white.data(tier),
                             xmv_overflow_black.data(tier), ymv_overflow_black.data(tier),
                             zmv_overflow_black.data(tier), xmv_overflow_white.data(tier),
                             ymv_overflow_white.data(tier), zmv_overflow_white.data(tier),
                             rxmv_black.data(tier), rymv_black.data(tier), rzmv_black.data(tier),
                             rxmv_white.data(tier), rymv_white.data(tier), rzmv_white.data(tier),
                             rxmv_overflow_black.data(tier), rymv_overflow_black.data(tier),
                             rzmv_overflow_black.data(tier), rxmv_overflow_white.data(tier),
                             rymv_overflow_white.data(tier), rzmv_overflow_white.data(tier),
                             inertial_tensor_black.data(tier), inertial_tensor_white.data(tier),
                             inertial_tensor_overflow_black.data(tier),
                             inertial_tensor_overflow_white.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MotionSweepWriter MotionSweeper::data(const HybridTargetLevel tier) {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const MotionSweepReader MotionSweeper::data(const CoordinateCycle orientation,
                                            const HybridTargetLevel tier) const {
  const double com_fac = pow(2.0, center_of_mass_bit_count);
  const double mv_fac = pow(2.0, momentum_bit_count);
  const double inrt_fac = pow(2.0, inertia_bit_count);
  switch (orientation) {
  case CoordinateCycle::WHITE:
    return MotionSweepReader(work_unit_count, work_units.data(tier), total_mass.data(tier),
                             com_fac, mv_fac, inrt_fac, xcom_white.data(tier),
                             ycom_white.data(tier), zcom_white.data(tier), xcom_black.data(tier),
                             ycom_black.data(tier), zcom_black.data(tier),
                             xcom_overflow_white.data(tier), ycom_overflow_white.data(tier),
                             zcom_overflow_white.data(tier), xcom_overflow_black.data(tier),
                             ycom_overflow_black.data(tier), zcom_overflow_black.data(tier),
                             xmv_white.data(tier), ymv_white.data(tier), zmv_white.data(tier),
                             xmv_black.data(tier), ymv_black.data(tier), zmv_black.data(tier),
                             xmv_overflow_white.data(tier), ymv_overflow_white.data(tier),
                             zmv_overflow_white.data(tier), xmv_overflow_black.data(tier),
                             ymv_overflow_black.data(tier), zmv_overflow_black.data(tier),
                             rxmv_white.data(tier), rymv_white.data(tier), rzmv_white.data(tier),
                             rxmv_black.data(tier), rymv_black.data(tier), rzmv_black.data(tier),
                             rxmv_overflow_white.data(tier), rymv_overflow_white.data(tier),
                             rzmv_overflow_white.data(tier), rxmv_overflow_black.data(tier),
                             rymv_overflow_black.data(tier), rzmv_overflow_black.data(tier),
                             inertial_tensor_white.data(tier), inertial_tensor_black.data(tier),
                             inertial_tensor_overflow_white.data(tier),
                             inertial_tensor_overflow_black.data(tier));
  case CoordinateCycle::BLACK:
    return MotionSweepReader(work_unit_count, work_units.data(tier), total_mass.data(tier),
                             com_fac, mv_fac, inrt_fac, xcom_black.data(tier),
                             ycom_black.data(tier), zcom_black.data(tier), xcom_white.data(tier),
                             ycom_white.data(tier), zcom_white.data(tier),
                             xcom_overflow_black.data(tier), ycom_overflow_black.data(tier),
                             zcom_overflow_black.data(tier), xcom_overflow_white.data(tier),
                             ycom_overflow_white.data(tier), zcom_overflow_white.data(tier),
                             xmv_black.data(tier), ymv_black.data(tier), zmv_black.data(tier),
                             xmv_white.data(tier), ymv_white.data(tier), zmv_white.data(tier),
                             xmv_overflow_black.data(tier), ymv_overflow_black.data(tier),
                             zmv_overflow_black.data(tier), xmv_overflow_white.data(tier),
                             ymv_overflow_white.data(tier), zmv_overflow_white.data(tier),
                             rxmv_black.data(tier), rymv_black.data(tier), rzmv_black.data(tier),
                             rxmv_white.data(tier), rymv_white.data(tier), rzmv_white.data(tier),
                             rxmv_overflow_black.data(tier), rymv_overflow_black.data(tier),
                             rzmv_overflow_black.data(tier), rxmv_overflow_white.data(tier),
                             rymv_overflow_white.data(tier), rzmv_overflow_white.data(tier),
                             inertial_tensor_black.data(tier), inertial_tensor_white.data(tier),
                             inertial_tensor_overflow_black.data(tier),
                             inertial_tensor_overflow_white.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const MotionSweepReader MotionSweeper::data(const HybridTargetLevel tier) const {
  return data(cycle_position, tier);
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void MotionSweeper::uploadWorkUnits() {
  work_units.upload();
}

//-------------------------------------------------------------------------------------------------
void MotionSweeper::downloadWorkUnits() {
  work_units.download();
}

//-------------------------------------------------------------------------------------------------
void MotionSweeper::uploadAll() {
  llint_data.upload();
  int_data.upload();
  total_mass.upload();
  work_units.upload();
}

//-------------------------------------------------------------------------------------------------
void MotionSweeper::downloadAll() {
  llint_data.download();
  int_data.download();
  total_mass.download();
  work_units.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void MotionSweeper::allocate() {
  const int nsys = poly_ps_ptr->getSystemCount();
  const int padded_nsys = roundUp(nsys, warp_size_int);
  const int padded_inrt = roundUp(nsys * 6, warp_size_int);
  llint_data.resize((padded_nsys * 18) + (2 * padded_inrt));
  int_data.resize((padded_nsys * 18) + (2 * padded_inrt));
  xcom_white.setPointer(&llint_data,                       0, nsys);
  ycom_white.setPointer(&llint_data,             padded_nsys, nsys);
  zcom_white.setPointer(&llint_data,         2 * padded_nsys, nsys);
  xcom_black.setPointer(&llint_data,         3 * padded_nsys, nsys);
  ycom_black.setPointer(&llint_data,         4 * padded_nsys, nsys);
  zcom_black.setPointer(&llint_data,         5 * padded_nsys, nsys);
  xmv_white.setPointer(&llint_data,          6 * padded_nsys, nsys);
  ymv_white.setPointer(&llint_data,          7 * padded_nsys, nsys);
  zmv_white.setPointer(&llint_data,          8 * padded_nsys, nsys);
  xmv_black.setPointer(&llint_data,          9 * padded_nsys, nsys);
  ymv_black.setPointer(&llint_data,         10 * padded_nsys, nsys);
  zmv_black.setPointer(&llint_data,         11 * padded_nsys, nsys);
  rxmv_white.setPointer(&llint_data,        12 * padded_nsys, nsys);
  rymv_white.setPointer(&llint_data,        13 * padded_nsys, nsys);
  rzmv_white.setPointer(&llint_data,        14 * padded_nsys, nsys);
  rxmv_black.setPointer(&llint_data,        15 * padded_nsys, nsys);
  rymv_black.setPointer(&llint_data,        16 * padded_nsys, nsys);
  rzmv_black.setPointer(&llint_data,        17 * padded_nsys, nsys);
  inertial_tensor_white.setPointer(&llint_data,  18 * padded_nsys               , 6 * nsys);
  inertial_tensor_black.setPointer(&llint_data, (18 * padded_nsys) + padded_inrt, 6 * nsys);
  xcom_overflow_white.setPointer(&int_data,                       0, nsys);
  ycom_overflow_white.setPointer(&int_data,             padded_nsys, nsys);
  zcom_overflow_white.setPointer(&int_data,         2 * padded_nsys, nsys);
  xcom_overflow_black.setPointer(&int_data,         3 * padded_nsys, nsys);
  ycom_overflow_black.setPointer(&int_data,         4 * padded_nsys, nsys);
  zcom_overflow_black.setPointer(&int_data,         5 * padded_nsys, nsys);
  xmv_overflow_white.setPointer(&int_data,          6 * padded_nsys, nsys);
  ymv_overflow_white.setPointer(&int_data,          7 * padded_nsys, nsys);
  zmv_overflow_white.setPointer(&int_data,          8 * padded_nsys, nsys);
  xmv_overflow_black.setPointer(&int_data,          9 * padded_nsys, nsys);
  ymv_overflow_black.setPointer(&int_data,         10 * padded_nsys, nsys);
  zmv_overflow_black.setPointer(&int_data,         11 * padded_nsys, nsys);
  rxmv_overflow_white.setPointer(&int_data,        12 * padded_nsys, nsys);
  rymv_overflow_white.setPointer(&int_data,        13 * padded_nsys, nsys);
  rzmv_overflow_white.setPointer(&int_data,        14 * padded_nsys, nsys);
  rxmv_overflow_black.setPointer(&int_data,        15 * padded_nsys, nsys);
  rymv_overflow_black.setPointer(&int_data,        16 * padded_nsys, nsys);
  rzmv_overflow_black.setPointer(&int_data,        17 * padded_nsys, nsys);
  inertial_tensor_overflow_white.setPointer(&int_data,  18 * padded_nsys               , 6 * nsys);
  inertial_tensor_overflow_black.setPointer(&int_data, (18 * padded_nsys) + padded_inrt, 6 * nsys);
}

} // namespace trajectory
} // namespace stormm
