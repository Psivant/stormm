#include "copyright.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "convolution_manager.h"
#include "pme_util.h"

namespace stormm {
namespace energy {

using card::HybridKind;
using synthesis::SyNonbondedKit;
using topology::UnitCellType;

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(const PMIGrid *pmig_in, const double ewald_coefficient_in) :
  system_count{pmig_in->getSystemCount()},
  ewald_coefficient{ewald_coefficient_in},
  system_offsets{static_cast<size_t>(system_count), "cvol_offsets"},
  self_ecorr{HybridKind::POINTER, "cvol_ecorr"},
  b_prefactor_a{HybridKind::POINTER, "cvol_bpref_a"},
  b_prefactor_b{HybridKind::POINTER, "cvol_bpref_b"},
  b_prefactor_c{HybridKind::POINTER, "cvol_bpref_c"},
  m_values_a{HybridKind::POINTER, "cvol_mval_a"},
  m_values_b{HybridKind::POINTER, "cvol_mval_b"},
  m_values_c{HybridKind::POINTER, "cvol_mval_c"},
  mshift_values_a{HybridKind::POINTER, "cvol_mshift_a"},
  mshift_values_b{HybridKind::POINTER, "cvol_mshift_b"},
  mshift_values_c{HybridKind::POINTER, "cvol_mshift_c"},
  c_prefactor_a{HybridKind::POINTER, "cvol_cpref_a"},
  c_prefactor_b{HybridKind::POINTER, "cvol_cpref_b"},
  c_prefactor_c{HybridKind::POINTER, "cvol_cpref_c"},
  double_data{HybridKind::ARRAY, "cvol_dbl_data"},
  sp_self_ecorr{HybridKind::POINTER, "sp_cvol_ecorr"},
  sp_b_prefactor_a{HybridKind::POINTER, "sp_cvol_bpref_a"},
  sp_b_prefactor_b{HybridKind::POINTER, "sp_cvol_bpref_b"},
  sp_b_prefactor_c{HybridKind::POINTER, "sp_cvol_bpref_c"},
  sp_m_values_a{HybridKind::POINTER, "sp_cvol_mval_a"},
  sp_m_values_b{HybridKind::POINTER, "sp_cvol_mval_b"},
  sp_m_values_c{HybridKind::POINTER, "sp_cvol_mval_c"},
  sp_mshift_values_a{HybridKind::POINTER, "sp_cvol_mshift_a"},
  sp_mshift_values_b{HybridKind::POINTER, "sp_cvol_mshift_b"},
  sp_mshift_values_c{HybridKind::POINTER, "sp_cvol_mshift_c"},
  sp_c_prefactor_a{HybridKind::POINTER, "sp_cvol_cpref_a"},
  sp_c_prefactor_b{HybridKind::POINTER, "sp_cvol_cpref_b"},
  sp_c_prefactor_c{HybridKind::POINTER, "sp_cvol_cpref_c"},
  float_data{HybridKind::ARRAY, "cvol_flt_data"},
  pmig_ptr{const_cast<PMIGrid*>(pmig_in)},
  poly_ag_ptr{const_cast<AtomGraphSynthesis*>(pmig_in->getTopologySynthesisPointer())}
{
  allocate();
  const int ordr = pmig_ptr->getInterpolationOrder();
  const UnitCellType uc = poly_ag_ptr->getUnitCellType();
  for (int i = 0; i < system_count; i++) {
    const uint4 i_gdims = pmig_ptr->getGridDimensions(i);
    const std::vector<double> ba = pmeLoadBPrefactor(ordr, i_gdims.x);
    const std::vector<double> bb = pmeLoadBPrefactor(ordr, i_gdims.y);
    const std::vector<double> bc = pmeLoadBPrefactor(ordr, i_gdims.z);
    const int sysi_offset = system_offsets.readHost(i);
    b_prefactor_a.putHost(ba, sysi_offset, i_gdims.x);
    b_prefactor_b.putHost(bb, sysi_offset, i_gdims.y);
    b_prefactor_c.putHost(bc, sysi_offset, i_gdims.z);
    const std::vector<float> fba(ba.begin(), ba.end());    
    const std::vector<float> fbb(bb.begin(), bb.end());    
    const std::vector<float> fbc(bc.begin(), bc.end());
    sp_b_prefactor_a.putHost(fba, sysi_offset, i_gdims.x);
    sp_b_prefactor_b.putHost(fbb, sysi_offset, i_gdims.y);
    sp_b_prefactor_c.putHost(fbc, sysi_offset, i_gdims.z);
    const std::vector<double> mval_a = pmeLoadMVec(i_gdims.x);
    const std::vector<double> mval_b = pmeLoadMVec(i_gdims.y);
    const std::vector<double> mval_c = pmeLoadMVec(i_gdims.z);
    m_values_a.putHost(mval_a, sysi_offset, i_gdims.x);
    m_values_b.putHost(mval_b, sysi_offset, i_gdims.y);
    m_values_c.putHost(mval_c, sysi_offset, i_gdims.z);
    const std::vector<float> sp_mval_a(mval_a.begin(), mval_a.end());
    const std::vector<float> sp_mval_b(mval_b.begin(), mval_b.end());
    const std::vector<float> sp_mval_c(mval_c.begin(), mval_c.end());
    sp_m_values_a.putHost(sp_mval_a, sysi_offset, i_gdims.x);
    sp_m_values_b.putHost(sp_mval_b, sysi_offset, i_gdims.y);
    sp_m_values_c.putHost(sp_mval_c, sysi_offset, i_gdims.z);
    const std::vector<double> mvs_a = pmeLoadMVecShift(i_gdims.x);
    const std::vector<double> mvs_b = pmeLoadMVecShift(i_gdims.y);
    const std::vector<double> mvs_c = pmeLoadMVecShift(i_gdims.z);
    mshift_values_a.putHost(mvs_a, sysi_offset, i_gdims.x);
    mshift_values_b.putHost(mvs_a, sysi_offset, i_gdims.y);
    mshift_values_c.putHost(mvs_a, sysi_offset, i_gdims.z);
    const std::vector<float> sp_mvs_a(mvs_a.begin(), mvs_a.end());
    const std::vector<float> sp_mvs_b(mvs_b.begin(), mvs_b.end());
    const std::vector<float> sp_mvs_c(mvs_c.begin(), mvs_c.end());
    sp_mshift_values_a.putHost(sp_mvs_a, sysi_offset, i_gdims.x);
    sp_mshift_values_b.putHost(sp_mvs_b, sysi_offset, i_gdims.y);
    sp_mshift_values_c.putHost(sp_mvs_c, sysi_offset, i_gdims.z);
    switch (uc) {
    case UnitCellType::NONE:
      rtErr("Convolutions in reiprocal space can only be performed with periodic systems.  A " +
            getEnumerationName(UnitCellType::ORTHORHOMBIC) + " or " +
            getEnumerationName(UnitCellType::TRICLINIC) + " unit cell type must be specified for "
            "system index " + std::to_string(i) + ".", "ConvolutionManager");
    case UnitCellType::ORTHORHOMBIC:

      // The C mesh prefactors can be loaded as a further optimization for rectilinear boxes.
      break;
    case UnitCellType::TRICLINIC:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(const PMIGrid &pmig_in, const double ewald_coefficient_in) :
  ConvolutionManager(pmig_in.getSelfPointer(), ewald_coefficient_in)
{}

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(const ConvolutionManager &original) :
    system_count{original.system_count},
    ewald_coefficient{original.ewald_coefficient},
    system_offsets{original.system_offsets},
    self_ecorr{original.self_ecorr},
    b_prefactor_a{original.b_prefactor_a},
    b_prefactor_b{original.b_prefactor_b},
    b_prefactor_c{original.b_prefactor_c},
    m_values_a{original.m_values_a},
    m_values_b{original.m_values_b},
    m_values_c{original.m_values_c},
    mshift_values_a{original.mshift_values_a},
    mshift_values_b{original.mshift_values_b},
    mshift_values_c{original.mshift_values_c},
    c_prefactor_a{original.c_prefactor_a},
    c_prefactor_b{original.c_prefactor_b},
    c_prefactor_c{original.c_prefactor_c},
    double_data{original.double_data},
    sp_self_ecorr{original.sp_self_ecorr},
    sp_b_prefactor_a{original.sp_b_prefactor_a},
    sp_b_prefactor_b{original.sp_b_prefactor_b},
    sp_b_prefactor_c{original.sp_b_prefactor_c},
    sp_m_values_a{original.sp_m_values_a},
    sp_m_values_b{original.sp_m_values_b},
    sp_m_values_c{original.sp_m_values_c},
    sp_mshift_values_a{original.sp_mshift_values_a},
    sp_mshift_values_b{original.sp_mshift_values_b},
    sp_mshift_values_c{original.sp_mshift_values_c},
    sp_c_prefactor_a{original.sp_c_prefactor_a},
    sp_c_prefactor_b{original.sp_c_prefactor_b},
    sp_c_prefactor_c{original.sp_c_prefactor_c},
    float_data{original.float_data},
    pmig_ptr{const_cast<PMIGrid*>(original.pmig_ptr)},
    poly_ag_ptr{const_cast<AtomGraphSynthesis*>(original.poly_ag_ptr)}
{}

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(ConvolutionManager &&original) :
    system_count{original.system_count},
    ewald_coefficient{original.ewald_coefficient},
    system_offsets{std::move(original.system_offsets)},
    self_ecorr{std::move(original.self_ecorr)},
    b_prefactor_a{std::move(original.b_prefactor_a)},
    b_prefactor_b{std::move(original.b_prefactor_b)},
    b_prefactor_c{std::move(original.b_prefactor_c)},
    m_values_a{std::move(original.m_values_a)},
    m_values_b{std::move(original.m_values_b)},
    m_values_c{std::move(original.m_values_c)},
    mshift_values_a{std::move(original.mshift_values_a)},
    mshift_values_b{std::move(original.mshift_values_b)},
    mshift_values_c{std::move(original.mshift_values_c)},
    c_prefactor_a{std::move(original.c_prefactor_a)},
    c_prefactor_b{std::move(original.c_prefactor_b)},
    c_prefactor_c{std::move(original.c_prefactor_c)},
    double_data{std::move(original.double_data)},
    sp_self_ecorr{std::move(original.sp_self_ecorr)},
    sp_b_prefactor_a{std::move(original.sp_b_prefactor_a)},
    sp_b_prefactor_b{std::move(original.sp_b_prefactor_b)},
    sp_b_prefactor_c{std::move(original.sp_b_prefactor_c)},
    sp_m_values_a{std::move(original.sp_m_values_a)},
    sp_m_values_b{std::move(original.sp_m_values_b)},
    sp_m_values_c{std::move(original.sp_m_values_c)},
    sp_mshift_values_a{std::move(original.sp_mshift_values_a)},
    sp_mshift_values_b{std::move(original.sp_mshift_values_b)},
    sp_mshift_values_c{std::move(original.sp_mshift_values_c)},
    sp_c_prefactor_a{std::move(original.sp_c_prefactor_a)},
    sp_c_prefactor_b{std::move(original.sp_c_prefactor_b)},
    sp_c_prefactor_c{std::move(original.sp_c_prefactor_c)},
    float_data{std::move(original.float_data)},
    pmig_ptr{original.pmig_ptr},
    poly_ag_ptr{original.poly_ag_ptr}
{}

//-------------------------------------------------------------------------------------------------
ConvolutionManager& ConvolutionManager::operator=(const ConvolutionManager &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  system_count = other.system_count;
  ewald_coefficient = other.ewald_coefficient;
  system_offsets = other.system_offsets;
  self_ecorr = other.self_ecorr;
  b_prefactor_a = other.b_prefactor_a;
  b_prefactor_b = other.b_prefactor_b;
  b_prefactor_c = other.b_prefactor_c;
  m_values_a = other.m_values_a;
  m_values_b = other.m_values_b;
  m_values_c = other.m_values_c;
  mshift_values_a = other.mshift_values_a;
  mshift_values_b = other.mshift_values_b;
  mshift_values_c = other.mshift_values_c;
  c_prefactor_a = other.c_prefactor_a;
  c_prefactor_b = other.c_prefactor_b;
  c_prefactor_c = other.c_prefactor_c;
  double_data = other.double_data;
  sp_self_ecorr = other.sp_self_ecorr;
  sp_b_prefactor_a = other.sp_b_prefactor_a;
  sp_b_prefactor_b = other.sp_b_prefactor_b;
  sp_b_prefactor_c = other.sp_b_prefactor_c;
  sp_m_values_a = other.sp_m_values_a;
  sp_m_values_b = other.sp_m_values_b;
  sp_m_values_c = other.sp_m_values_c;
  sp_mshift_values_a = other.sp_mshift_values_a;
  sp_mshift_values_b = other.sp_mshift_values_b;
  sp_mshift_values_c = other.sp_mshift_values_c;
  sp_c_prefactor_a = other.sp_c_prefactor_a;
  sp_c_prefactor_b = other.sp_c_prefactor_b;
  sp_c_prefactor_c = other.sp_c_prefactor_c;
  float_data = other.float_data;
  pmig_ptr = other.pmig_ptr;
  poly_ag_ptr = other.poly_ag_ptr;

  // Use the allocator to repair pointers
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
ConvolutionManager& ConvolutionManager::operator=(ConvolutionManager &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  system_count = other.system_count;
  ewald_coefficient = other.ewald_coefficient;
  system_offsets = other.system_offsets;
  self_ecorr = other.self_ecorr;
  b_prefactor_a = other.b_prefactor_a;
  b_prefactor_b = other.b_prefactor_b;
  b_prefactor_c = other.b_prefactor_c;
  m_values_a = other.m_values_a;
  m_values_b = other.m_values_b;
  m_values_c = other.m_values_c;
  mshift_values_a = other.mshift_values_a;
  mshift_values_b = other.mshift_values_b;
  mshift_values_c = other.mshift_values_c;
  c_prefactor_a = other.c_prefactor_a;
  c_prefactor_b = other.c_prefactor_b;
  c_prefactor_c = other.c_prefactor_c;
  double_data = other.double_data;
  sp_self_ecorr = other.sp_self_ecorr;
  sp_b_prefactor_a = other.sp_b_prefactor_a;
  sp_b_prefactor_b = other.sp_b_prefactor_b;
  sp_b_prefactor_c = other.sp_b_prefactor_c;
  sp_m_values_a = other.sp_m_values_a;
  sp_m_values_b = other.sp_m_values_b;
  sp_m_values_c = other.sp_m_values_c;
  sp_mshift_values_a = other.sp_mshift_values_a;
  sp_mshift_values_b = other.sp_mshift_values_b;
  sp_mshift_values_c = other.sp_mshift_values_c;
  sp_c_prefactor_a = other.sp_c_prefactor_a;
  sp_c_prefactor_b = other.sp_c_prefactor_b;
  sp_c_prefactor_c = other.sp_c_prefactor_c;
  float_data = other.float_data;
  pmig_ptr = other.pmig_ptr;
  poly_ag_ptr = other.poly_ag_ptr;

  // No pointer repair is needed in most move assignments and move operations
  return *this;
}

//-------------------------------------------------------------------------------------------------
double ConvolutionManager::getEwaldCoefficient() const {
  return ewald_coefficient;
}

//-------------------------------------------------------------------------------------------------
double ConvolutionManager::getCoulombConstant() const {
  return poly_ag_ptr->getCoulombConstant();
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ConvolutionManager::upload() {
  system_offsets.upload();
  double_data.upload();
  float_data.upload();
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::download() {
  system_offsets.download();
  double_data.download();
  float_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::allocate() {
  system_offsets.resize(system_count);
  const PMIGridReader pmir = pmig_ptr->data();
  int nreal = roundUp(system_count, warp_size_int);
  for (int i = 0; i < system_count; i++) {
    const uint4 i_dims = pmir.dims[i];
    const int maxdim = std::max(std::max(i_dims.x, i_dims.y), i_dims.z);
    system_offsets.putHost(nreal, i);
    nreal += roundUp(maxdim, warp_size_int);
  }
  double_data.resize(13 * nreal);
  self_ecorr.setPointer(&double_data,               0);
  b_prefactor_a.setPointer(&double_data,        nreal);
  b_prefactor_b.setPointer(&double_data,    2 * nreal);
  b_prefactor_c.setPointer(&double_data,    3 * nreal);
  m_values_a.setPointer(&double_data,       4 * nreal);
  m_values_b.setPointer(&double_data,       5 * nreal);
  m_values_c.setPointer(&double_data,       6 * nreal);
  mshift_values_a.setPointer(&double_data,  7 * nreal);
  mshift_values_b.setPointer(&double_data,  8 * nreal);
  mshift_values_c.setPointer(&double_data,  9 * nreal);
  c_prefactor_a.setPointer(&double_data,   10 * nreal);
  c_prefactor_b.setPointer(&double_data,   11 * nreal);
  c_prefactor_c.setPointer(&double_data,   12 * nreal);
  float_data.resize(13 * nreal);
  sp_self_ecorr.setPointer(&float_data,              0);
  sp_b_prefactor_a.setPointer(&float_data,        nreal);
  sp_b_prefactor_b.setPointer(&float_data,    2 * nreal);
  sp_b_prefactor_c.setPointer(&float_data,    3 * nreal);
  sp_m_values_a.setPointer(&float_data,       4 * nreal);
  sp_m_values_b.setPointer(&float_data,       5 * nreal);
  sp_m_values_c.setPointer(&float_data,       6 * nreal);
  sp_mshift_values_a.setPointer(&float_data,  7 * nreal);
  sp_mshift_values_b.setPointer(&float_data,  8 * nreal);
  sp_mshift_values_c.setPointer(&float_data,  9 * nreal);
  sp_c_prefactor_a.setPointer(&float_data,   10 * nreal);
  sp_c_prefactor_b.setPointer(&float_data,   11 * nreal);
  sp_c_prefactor_c.setPointer(&float_data,   12 * nreal);
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::computeSystemSelfEnergies() {
  const SyNonbondedKit<double, double2> poly_nbk = poly_ag_ptr->getDoublePrecisionNonbondedKit();
  for (int i = 0; i < poly_nbk.nsys; i++) {
    double tmp_ec = 0.0;
    const size_t jlim = poly_nbk.atom_offsets[i] + poly_nbk.atom_counts[i];
    for (size_t j = poly_nbk.atom_offsets[i]; j < jlim; j++) {
      tmp_ec += poly_nbk.charge[j] * poly_nbk.charge[j];
    }
    tmp_ec *= -poly_nbk.coulomb * ewald_coefficient * sqrt(inverse_pi);
  }
}

} // namespace energy
} // namespace stormm
