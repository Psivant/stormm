#include "copyright.h"
#include "Accelerator/cuda_wrappers.h"
#include "Math/math_enumerators.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "convolution_manager.h"
#include "pme_util.h"

namespace stormm {
namespace energy {

#ifdef STORMM_USE_HPC
using card::getHpcErrorString;
#endif
using card::HybridFormat;
using card::HybridKind;
using stmath::Normalization;
using synthesis::SyNonbondedKit;
using topology::UnitCellType;

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(const PMIGrid *pmig_in, const double ewald_coefficient_in,
                                       const GpuDetails &gpu) :
    system_count{pmig_in->getSystemCount()},
    hpc_work_unit_count{0},
    ewald_coefficient{ewald_coefficient_in},
    prefactor_offsets{static_cast<size_t>(system_count), "cvol_offsets"},
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
#ifdef STORMM_USE_HPC
    operating_tier{(gpu == null_gpu) ? HybridTargetLevel::HOST : HybridTargetLevel::DEVICE},
#else
    operating_tier{HybridTargetLevel::HOST},
#endif
    fft_group_count{0}, fft_groups{}, fft_group_bounds{},
#ifdef STORMM_USE_HPC
    dp_frequency_data{HybridKind::ARRAY, "dp_freq_stash",
                      (gpu == null_gpu) ? HybridFormat::HOST_MOUNTED : HybridFormat::DEVICE_ONLY},
    sp_frequency_data{HybridKind::ARRAY, "sp_freq_stash",
                      (gpu == null_gpu) ? HybridFormat::HOST_MOUNTED : HybridFormat::DEVICE_ONLY},
#else
    dp_frequency_data{HybridKind::ARRAY, "dp_freq_stash", HybridFormat::HOST_ONLY},
    sp_frequency_data{HybridKind::ARRAY, "sp_freq_stash", HybridFormat::HOST_ONLY},
#endif
    frequency_offsets{HybridKind::ARRAY, "cvol_freq_offsets"},
    fft_operations{},
    work_units{HybridKind::ARRAY, "cvol_work_units"},
    pmig_ptr{const_cast<PMIGrid*>(pmig_in)},
    poly_ag_ptr{const_cast<AtomGraphSynthesis*>(pmig_in->getTopologySynthesisPointer())}
{
  allocate();
  const int ordr = pmig_ptr->getInterpolationOrder();
  const UnitCellType uc = poly_ag_ptr->getUnitCellType();
  for (int i = 0; i < system_count; i++) {
    const uint4 i_gdims = pmig_ptr->getGridDimensions(i);
    const AtomGraph *iag_ptr = poly_ag_ptr->getSystemTopologyPointer(i);
    const NonbondedKit<double> inbk = iag_ptr->getDoublePrecisionNonbondedKit();
    double qself_i = 0.0;
    for (int j = 0; j < inbk.natom; j++) {
      qself_i += inbk.charge[j] * inbk.charge[j];
    }
    qself_i *= -poly_ag_ptr->getCoulombConstant() * ewald_coefficient / sqrt(pi);
    self_ecorr.putHost(qself_i, i);
    std::vector<double> ba = pmeLoadBPrefactor(ordr, i_gdims.x);
    const std::vector<double> bb = pmeLoadBPrefactor(ordr, i_gdims.y);
    const std::vector<double> bc = pmeLoadBPrefactor(ordr, i_gdims.z);
    const int sysi_offset = prefactor_offsets.readHost(i);
    const double coulomb_c = pmig_in->getTopologySynthesisPointer()->getCoulombConstant();
    for (int j = 0; j < i_gdims.x; j++) {
      ba[j] *= coulomb_c;
    }
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
    mshift_values_b.putHost(mvs_b, sysi_offset, i_gdims.y);
    mshift_values_c.putHost(mvs_c, sysi_offset, i_gdims.z);
    const std::vector<float> sp_mvs_a(mvs_a.begin(), mvs_a.end());
    const std::vector<float> sp_mvs_b(mvs_b.begin(), mvs_b.end());
    const std::vector<float> sp_mvs_c(mvs_c.begin(), mvs_c.end());
    sp_mshift_values_a.putHost(sp_mvs_a, sysi_offset, i_gdims.x);
    sp_mshift_values_b.putHost(sp_mvs_b, sysi_offset, i_gdims.y);
    sp_mshift_values_c.putHost(sp_mvs_c, sysi_offset, i_gdims.z);
    switch (uc) {
    case UnitCellType::NONE:
      rtErr("Convolutions in reciprocal space can only be performed with periodic systems.  A " +
            getEnumerationName(UnitCellType::ORTHORHOMBIC) + " or " +
            getEnumerationName(UnitCellType::TRICLINIC) + " unit cell type must be specified for "
            "system index " + std::to_string(i) + ".", "ConvolutionManager");
    case UnitCellType::ORTHORHOMBIC:
      {
        // The C mesh prefactors can be loaded as a further optimization for rectilinear boxes.
        // These prefactors require knowledge of the exact lengths of each system's rectilinear
        // box, and will not be accurate
        const std::vector<double> ca = pmeLoadOrthoCPrefactor(i_gdims, ewald_coefficient,
                                                              m_values_a.readHost(),
                                                              UnitCellAxis::A);
        const std::vector<double> cb = pmeLoadOrthoCPrefactor(i_gdims, ewald_coefficient,
                                                              m_values_b.readHost(),
                                                              UnitCellAxis::B);
        const std::vector<double> cc = pmeLoadOrthoCPrefactor(i_gdims, ewald_coefficient,
                                                              m_values_c.readHost(),
                                                              UnitCellAxis::C);
        c_prefactor_a.putHost(ca, sysi_offset, i_gdims.x);
        c_prefactor_b.putHost(cb, sysi_offset, i_gdims.y);
        c_prefactor_c.putHost(cc, sysi_offset, i_gdims.z);
        const std::vector<float> sp_ca(ca.begin(), ca.end());
        const std::vector<float> sp_cb(cb.begin(), cb.end());
        const std::vector<float> sp_cc(cc.begin(), cc.end());
        sp_c_prefactor_a.putHost(sp_ca, sysi_offset, i_gdims.x);
        sp_c_prefactor_b.putHost(sp_cb, sysi_offset, i_gdims.y);
        sp_c_prefactor_c.putHost(sp_cc, sysi_offset, i_gdims.z);
      }
      break;
    case UnitCellType::TRICLINIC:
      break;
    }
  }
  
  // Prepare for FFT operations
  const PMIGridReader pmigr = PMIGridReader(pmig_ptr->data());
  frequency_offsets.resize(system_count);
  switch (operating_tier) {
  case HybridTargetLevel::HOST:
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    this->upload();
    break;
#endif
  }
  switch (pmig_ptr->getFFTStaging()) {
  case FFTMode::IN_PLACE:

    // Set the complex type starting marks for each frequency-space grid.
    for (int i = 0; i < system_count; i++) {
      frequency_offsets.putHost(pmigr.dims[i].w / 2, i);
    }
    break;
  case FFTMode::OUT_OF_PLACE:
    {
      // Create padded offsets for each system's frequency-space grid, then allocate memory for
      // the out-of-place transform results.
      int ijcon = 0;
      uint vol_prog = 0;
      for (int i = 0; i < system_count; i++) {
        uint4 ig_len = pmigr.dims[i];
        uint ig_vol = roundUp<uint>(((ig_len.x / 2) + 1) * ig_len.y * ig_len.z, 16);
        frequency_offsets.putHost(vol_prog, i);
        vol_prog += ig_vol;
      }
      switch (pmigr.mode) {
      case PrecisionModel::DOUBLE:
        dp_frequency_data.resize(vol_prog);
        break;
      case PrecisionModel::SINGLE: 
        sp_frequency_data.resize(vol_prog);
        break;
      }
    }
    break;
  }
  makeFFTGroups(gpu);

  // Create work units to organize GPU operations
  makeWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(const PMIGrid &pmig_in, const double ewald_coefficient_in,
                                       const GpuDetails &gpu) :
    ConvolutionManager(pmig_in.getSelfPointer(), ewald_coefficient_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(const ConvolutionManager &original) :
    system_count{original.system_count},
    hpc_work_unit_count{original.hpc_work_unit_count},
    ewald_coefficient{original.ewald_coefficient},
    prefactor_offsets{original.prefactor_offsets},
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
    operating_tier{original.operating_tier},
    fft_group_count{original.fft_group_count},
    fft_groups{original.fft_groups},
    fft_group_bounds{original.fft_group_bounds},
    dp_frequency_data{original.dp_frequency_data},
    sp_frequency_data{original.sp_frequency_data},
    frequency_offsets{original.frequency_offsets},
    fft_operations{original.fft_operations},
    work_units{original.work_units},
    pmig_ptr{const_cast<PMIGrid*>(original.pmig_ptr)},
    poly_ag_ptr{const_cast<AtomGraphSynthesis*>(original.poly_ag_ptr)}
{}

//-------------------------------------------------------------------------------------------------
ConvolutionManager::ConvolutionManager(ConvolutionManager &&original) :
    system_count{original.system_count},
    hpc_work_unit_count{original.hpc_work_unit_count},
    ewald_coefficient{original.ewald_coefficient},
    prefactor_offsets{std::move(original.prefactor_offsets)},
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
    operating_tier{original.operating_tier},
    fft_group_count{original.fft_group_count},
    fft_groups{std::move(original.fft_groups)},
    fft_group_bounds{std::move(original.fft_group_bounds)},
    dp_frequency_data{std::move(original.dp_frequency_data)},
    sp_frequency_data{std::move(original.sp_frequency_data)},
    frequency_offsets{std::move(original.frequency_offsets)},
    fft_operations{std::move(original.fft_operations)},
    work_units{std::move(original.work_units)},
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
  hpc_work_unit_count = other.hpc_work_unit_count;
  ewald_coefficient = other.ewald_coefficient;
  prefactor_offsets = other.prefactor_offsets;
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
  operating_tier = other.operating_tier;
  fft_group_count = other.fft_group_count;
  fft_groups = other.fft_groups;
  fft_group_bounds = other.fft_group_bounds;
  dp_frequency_data = other.dp_frequency_data;
  sp_frequency_data = other.sp_frequency_data;
  frequency_offsets = other.frequency_offsets;
  fft_operations = other.fft_operations;
  work_units = other.work_units;
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
  hpc_work_unit_count = other.hpc_work_unit_count;
  ewald_coefficient = other.ewald_coefficient;
  prefactor_offsets = std::move(other.prefactor_offsets);
  self_ecorr = std::move(other.self_ecorr);
  b_prefactor_a = std::move(other.b_prefactor_a);
  b_prefactor_b = std::move(other.b_prefactor_b);
  b_prefactor_c = std::move(other.b_prefactor_c);
  m_values_a = std::move(other.m_values_a);
  m_values_b = std::move(other.m_values_b);
  m_values_c = std::move(other.m_values_c);
  mshift_values_a = std::move(other.mshift_values_a);
  mshift_values_b = std::move(other.mshift_values_b);
  mshift_values_c = std::move(other.mshift_values_c);
  c_prefactor_a = std::move(other.c_prefactor_a);
  c_prefactor_b = std::move(other.c_prefactor_b);
  c_prefactor_c = std::move(other.c_prefactor_c);
  double_data = std::move(other.double_data);
  sp_self_ecorr = std::move(other.sp_self_ecorr);
  sp_b_prefactor_a = std::move(other.sp_b_prefactor_a);
  sp_b_prefactor_b = std::move(other.sp_b_prefactor_b);
  sp_b_prefactor_c = std::move(other.sp_b_prefactor_c);
  sp_m_values_a = std::move(other.sp_m_values_a);
  sp_m_values_b = std::move(other.sp_m_values_b);
  sp_m_values_c = std::move(other.sp_m_values_c);
  sp_mshift_values_a = std::move(other.sp_mshift_values_a);
  sp_mshift_values_b = std::move(other.sp_mshift_values_b);
  sp_mshift_values_c = std::move(other.sp_mshift_values_c);
  sp_c_prefactor_a = std::move(other.sp_c_prefactor_a);
  sp_c_prefactor_b = std::move(other.sp_c_prefactor_b);
  sp_c_prefactor_c = std::move(other.sp_c_prefactor_c);
  float_data = std::move(other.float_data);
  operating_tier = other.operating_tier;
  fft_group_count = other.fft_group_count;
  fft_groups = std::move(other.fft_groups);
  fft_group_bounds = std::move(other.fft_group_bounds);
  dp_frequency_data = std::move(other.dp_frequency_data);
  sp_frequency_data = std::move(other.sp_frequency_data);
  frequency_offsets = std::move(other.frequency_offsets);
  fft_operations = std::move(other.fft_operations);
  work_units = std::move(other.work_units);
  pmig_ptr = other.pmig_ptr;
  poly_ag_ptr = other.poly_ag_ptr;

  // No pointer repair is needed in most move assignments or move operations performed in STORMM
  // classes.  This object is typical.
  return *this;
}

//-------------------------------------------------------------------------------------------------
int ConvolutionManager::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int ConvolutionManager::getHPCWorkUnitCount() const {
  return hpc_work_unit_count;
}

//-------------------------------------------------------------------------------------------------
double ConvolutionManager::getEwaldCoefficient() const {
  return ewald_coefficient;
}

//-------------------------------------------------------------------------------------------------
double ConvolutionManager::getCoulombConstant() const {
  return poly_ag_ptr->getCoulombConstant();
}

//-------------------------------------------------------------------------------------------------
const PMIGrid* ConvolutionManager::getPMIGridPointer() const {
  return pmig_ptr;
}

//-------------------------------------------------------------------------------------------------
const AtomGraphSynthesis* ConvolutionManager::getTopologySynthesisPointer() const {
  return poly_ag_ptr;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ConvolutionManager::getSelfEcorr(const PrecisionModel prec) const {
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return self_ecorr.readHost();
  case PrecisionModel::SINGLE:
    {
      const size_t n = sp_self_ecorr.size();
      std::vector<double> result(n);
      for (size_t i = 0; i < n; i++) {
        result[i] = sp_self_ecorr.readHost(i);
      }
      return result;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int ConvolutionManager::getFFTGroupCount() const {
  return fft_group_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ConvolutionManager::getFFTGroupSystemList(const int group_index) const {
  if (group_index < 0 || group_index >= fft_group_count) {
    rtErr("Group index " + std::to_string(group_index) + " is invalid for an object with " +
          std::to_string(fft_group_count) + " FFT groups.", "ConvolutionManager",
          "getFFTGroupSystemList");
  }
  const int n_problems = fft_group_bounds[group_index + 1] - fft_group_bounds[group_index];
  std::vector<int> result(n_problems);
  for (int i = 0; i < n_problems; i++) {
    result[i] = fft_groups[fft_group_bounds[group_index] + i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
HybridTargetLevel ConvolutionManager::getOperatingTier() const {
  return operating_tier;
}

//-------------------------------------------------------------------------------------------------
const ConvolutionReader<double, double2> ConvolutionManager::dpData() const {
  const HybridTargetLevel op_t = operating_tier;
  return ConvolutionReader<double,
                           double2>(system_count, ewald_coefficient, prefactor_offsets.data(op_t),
                                    self_ecorr.data(op_t), b_prefactor_a.data(op_t),
                                    b_prefactor_b.data(op_t), b_prefactor_c.data(op_t),
                                    m_values_a.data(op_t), m_values_b.data(op_t),
                                    m_values_c.data(op_t), mshift_values_a.data(op_t),
                                    mshift_values_b.data(op_t), mshift_values_c.data(op_t),
                                    c_prefactor_a.data(op_t), c_prefactor_b.data(op_t),
                                    c_prefactor_c.data(op_t), dp_frequency_data.data(op_t),
                                    frequency_offsets.data(op_t), &fft_operations);
}

//-------------------------------------------------------------------------------------------------
ConvolutionWriter<double, double2> ConvolutionManager::dpData() {
  const HybridTargetLevel op_t = operating_tier;
  return ConvolutionWriter<double,
                           double2>(system_count, hpc_work_unit_count, ewald_coefficient,
                                    prefactor_offsets.data(op_t), self_ecorr.data(op_t),
                                    b_prefactor_a.data(op_t), b_prefactor_b.data(op_t),
                                    b_prefactor_c.data(op_t), m_values_a.data(op_t),
                                    m_values_b.data(op_t), m_values_c.data(op_t),
                                    mshift_values_a.data(op_t), mshift_values_b.data(op_t),
                                    mshift_values_c.data(op_t), c_prefactor_a.data(op_t),
                                    c_prefactor_b.data(op_t), c_prefactor_c.data(op_t),
                                    dp_frequency_data.data(op_t), frequency_offsets.data(op_t),
                                    work_units.data(op_t), &fft_operations);
}

//-------------------------------------------------------------------------------------------------
const ConvolutionReader<float, float2> ConvolutionManager::spData() const {
  const HybridTargetLevel op_t = operating_tier;
  return ConvolutionReader<float,
                           float2>(system_count, ewald_coefficient, prefactor_offsets.data(op_t),
                                   sp_self_ecorr.data(op_t), sp_b_prefactor_a.data(op_t),
                                   sp_b_prefactor_b.data(op_t), sp_b_prefactor_c.data(op_t),
                                   sp_m_values_a.data(op_t), sp_m_values_b.data(op_t),
                                   sp_m_values_c.data(op_t), sp_mshift_values_a.data(op_t),
                                   sp_mshift_values_b.data(op_t), sp_mshift_values_c.data(op_t),
                                   sp_c_prefactor_a.data(op_t), sp_c_prefactor_b.data(op_t),
                                   sp_c_prefactor_c.data(op_t), sp_frequency_data.data(op_t),
                                   frequency_offsets.data(op_t), &fft_operations);
}

//-------------------------------------------------------------------------------------------------
ConvolutionWriter<float, float2> ConvolutionManager::spData() {
  const HybridTargetLevel op_t = operating_tier;
  return ConvolutionWriter<float,
                           float2>(system_count, hpc_work_unit_count, ewald_coefficient,
                                   prefactor_offsets.data(op_t), sp_self_ecorr.data(op_t),
                                   sp_b_prefactor_a.data(op_t), sp_b_prefactor_b.data(op_t),
                                   sp_b_prefactor_c.data(op_t), sp_m_values_a.data(op_t),
                                   sp_m_values_b.data(op_t), sp_m_values_c.data(op_t),
                                   sp_mshift_values_a.data(op_t), sp_mshift_values_b.data(op_t),
                                   sp_mshift_values_c.data(op_t), sp_c_prefactor_a.data(op_t),
                                   sp_c_prefactor_b.data(op_t), sp_c_prefactor_c.data(op_t),
                                   sp_frequency_data.data(op_t), frequency_offsets.data(op_t),
                                   work_units.data(op_t), &fft_operations);
}

//-------------------------------------------------------------------------------------------------
const ConvolutionManager* ConvolutionManager::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::forwardFFT() {
  for (int i = 0; i < fft_group_count; i++) {
    fft_operations[i].forwardFFT();
  }
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::applyGreensFunction(const GpuDetails &gpu, ScoreCard *sc) {
  const PhaseSpaceSynthesis *poly_ps_ptr = pmig_ptr->getCoordinateSynthesisPointer();
  const PMIGridWriter pmigw = pmig_ptr->data(operating_tier);
  const PsSynthesisBorders pssb = poly_ps_ptr->borders(operating_tier);
  switch (poly_ag_ptr->getUnitCellType()) {
  case UnitCellType::NONE:
    rtErr("A PME treatment of the energy is not possible without periodic boundary conditions.",
          "ConvolutionManager", "applyGreensFunction");
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:

    // Perform the convolution on the resource for which the FFTs are prepared, and at the same
    // precision as the FFTs.  There would be little benefit in terms of precision to perform the
    // element-wise multiplication, or composition of the B and C meshes, in double-precision
    // while the FFTs before and after such an operation are performed in single-precision.
    switch (operating_tier) {
    case HybridTargetLevel::HOST:
      switch (pmig_ptr->getMode()) {
      case PrecisionModel::DOUBLE:
        {
          ConvolutionWriter<double, double2> cvolw = this->dpData();
          if (sc == nullptr) {
            pmeGreensFunction<double, double2>(&cvolw, pssb, pmigw);
          }
          else {
            ScoreCardWriter scw = sc->data();
            pmeGreensFunction<double, double2>(&cvolw, pssb, pmigw, &scw);
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          ConvolutionWriter<float, float2> cvolw = this->spData();
          if (sc == nullptr) {
            pmeGreensFunction<float, float2>(&cvolw, pssb, pmigw);
          }
          else {
            ScoreCardWriter scw = sc->data();
            pmeGreensFunction<float, float2>(&cvolw, pssb, pmigw, &scw);
          }
        }
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        const PhaseSpaceSynthesis *poly_ps = pmig_ptr->getCoordinateSynthesisPointer();
        switch (pmig_ptr->getMode()) {
        case PrecisionModel::DOUBLE:
          {
            ConvolutionWriter<double, double2> cvolw = this->dpData();
            if (sc == nullptr) {
              pmeGreensFunction(&cvolw, pssb, pmigw, gpu);
            }
            else {
              ScoreCardWriter scw = sc->data(HybridTargetLevel::DEVICE);
              pmeGreensFunction(&cvolw, pssb, pmigw, gpu, &scw);
            }
          }
          break;
        case PrecisionModel::SINGLE:
          {
            ConvolutionWriter<float, float2> cvolw = this->spData();
            if (sc == nullptr) {
              pmeGreensFunction(&cvolw, pssb, pmigw, gpu);
            }
            else {
              ScoreCardWriter scw = sc->data(HybridTargetLevel::DEVICE);
              pmeGreensFunction(&cvolw, pssb, pmigw, gpu, &scw);
            }
          }
          break;
        }
      }
      break;
#endif
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::backwardFFT() {
  for (int i = 0; i < fft_group_count; i++) {
    fft_operations[i].backwardFFT();
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ConvolutionManager::upload() {
  prefactor_offsets.upload();
  double_data.upload();
  float_data.upload();
  frequency_offsets.upload();
  work_units.upload();
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::download() {
  prefactor_offsets.download();
  double_data.download();
  float_data.download();
  frequency_offsets.download();
  work_units.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::allocate() {
  prefactor_offsets.resize(system_count);
  const PMIGridReader pmir = pmig_ptr->data();
  int nreal = 0;
  for (int i = 0; i < system_count; i++) {
    const uint4 i_dims = pmir.dims[i];
    const int maxdim = std::max(std::max(i_dims.x, i_dims.y), i_dims.z);
    prefactor_offsets.putHost(nreal, i);
    nreal += roundUp(maxdim, warp_size_int);
  }
  nreal += roundUp(system_count, warp_size_int);
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
  frequency_offsets.resize(system_count);
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

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::makeFFTGroups(const GpuDetails &gpu) {

  // Devise groups for CPU computations
  std::vector<bool> coverage(system_count, false);
  fft_groups.resize(system_count);
  fft_group_bounds.resize(0);
  PMIGridWriter host_pmigw = pmig_ptr->data();
  PMIGridWriter pmigw = pmig_ptr->data(operating_tier);
  uint ngx, ngy, ngz;
  int n_covered = 0;
#ifdef STORMM_USE_HPC
  uint gpu_gbl_cache_vol;
  HybridTargetLevel fft_op_tier;
  if (gpu != null_gpu) {
    gpu_gbl_cache_vol = static_cast<double>(gpu.getGlobalCacheSize() * 0.1875);
    fft_op_tier = HybridTargetLevel::DEVICE;
  }
  else {
    gpu_gbl_cache_vol = 0;
    fft_op_tier = HybridTargetLevel::HOST;
  }
#else
  const uint gpu_gbl_cache_vol = 0;
  const HybridTargetLevel fft_op_tier = HybridTargetLevel::HOST;
#endif
  const FFTMode fft_op_mode = FFTMode::OUT_OF_PLACE;
  int syspos = 0;
  while (syspos < system_count) {

    // Find the dimensions of this group
    fft_group_count += 1;
    ngx = host_pmigw.dims[syspos].x;
    ngy = host_pmigw.dims[syspos].y;
    ngz = host_pmigw.dims[syspos].z;
    fft_group_bounds.push_back(n_covered);
    uint batch_vol = 0;
    while (syspos < system_count && batch_vol <= gpu_gbl_cache_vol &&
           host_pmigw.dims[syspos].x == ngx && host_pmigw.dims[syspos].y == ngy &&
           host_pmigw.dims[syspos].z == ngz) {
      fft_groups[n_covered] = syspos;
      coverage[syspos] = true;
      n_covered++;
      syspos++;
      batch_vol += (gpu_gbl_cache_vol > 0) ? ngx * ngy * ngz : 0;
    }
  }
  fft_group_bounds.push_back(system_count);

  // Reserve space for the array of FFT operations
  fft_operations.reserve(fft_group_count);

  // Create plans for carrying out FFTs, whether on the CPU or on the GPU
  const uint* frq_ofs_ptr = frequency_offsets.data();
  for (int i = 0; i < fft_group_count; i++) {
    const int example_system_id = fft_groups[fft_group_bounds[i]];
    const size_t ngx_i = host_pmigw.dims[example_system_id].x;
    const size_t ngy_i = host_pmigw.dims[example_system_id].y;
    const size_t ngz_i = host_pmigw.dims[example_system_id].z;
    const size_t q_data_pos = host_pmigw.dims[example_system_id].w;
    const int nset = fft_group_bounds[i + 1] - fft_group_bounds[i];
    const uint frq_ofs = frequency_offsets.readHost(example_system_id);
    size_t gstride, frq_stride;
    if (nset > 1) {
      gstride = host_pmigw.dims[example_system_id + 1].w - q_data_pos;
      frq_stride = frequency_offsets.readHost(example_system_id + 1) - frq_ofs;
    }
    else {
      gstride = ngx_i * ngy_i * ngz_i;
      frq_stride = ((ngx_i / 2) + 1) * ngy_i * ngz_i;
    }
    if (nset > 1) {

      // Get the stride between distinct grids.  Check that both strides are maintained across the
      // entire set of grids.
      for (int j = 0; j < nset - 1; j++) {
        const int jsys_id = example_system_id + j;
        if (host_pmigw.dims[jsys_id + 1].w - host_pmigw.dims[jsys_id].w != gstride ||
            frq_ofs_ptr[jsys_id + 1] - frq_ofs_ptr[jsys_id] != frq_stride) {
          rtErr("Inconsistent grid strides detected in a set of " + std::to_string(nset) +
                " grids with dimensions (" + std::to_string(ngx_i) + " x " +
                std::to_string(ngy_i) + " x " + std::to_string(ngz_i) + ".  The original stride "
                "is " + std::to_string(gstride) + ", but a subsequent stride of " +
                std::to_string(host_pmigw.dims[jsys_id + 1].w - host_pmigw.dims[jsys_id].w) +
                " breaks the pattern and thus the batched FFT.", "ConvolutionManager",
                "makeFFTGroups");
        }
      }
    }

    // Until this point, a host-oriented abstract from the particle mesh interaction grid object
    // has been used, to ensure that its data on the grid boundaries will be accessible to the CPU
    // as it creates FFT groups (these boundaries are identical on the CPU host and the GPU
    // device).  Going forward, the abstract appropriate to the tier at which FFTs will be
    // performed must be used, to ensure that the GPU does not do work involving arrays on the
    // host.  The constants in either abstract are consistent, and still visible to the CPU host.
    switch (pmigw.fftm) {
    case FFTMode::IN_PLACE:
      switch (pmigw.mode) {
      case PrecisionModel::DOUBLE:
        {
          double2* freq_ptr = reinterpret_cast<double2*>(pmigw.ddata);
          fft_operations.emplace_back(&pmigw.ddata[q_data_pos], &freq_ptr[frq_ofs], fft_op_tier,
                                      Normalization::YES, pmigw.fftm, ngx_i, ngy_i, ngz_i, 0, nset,
                                      gstride, frq_stride);
        }
        break;
      case PrecisionModel::SINGLE:
        {
          float2* freq_ptr = reinterpret_cast<float2*>(pmigw.fdata);
          fft_operations.emplace_back(&pmigw.fdata[q_data_pos], &freq_ptr[frq_ofs], fft_op_tier,
                                      Normalization::YES, pmigw.fftm, ngx_i, ngy_i, ngz_i, 0, nset,
                                      gstride, frq_stride);
        }
        break;
      }
      break;
    case FFTMode::OUT_OF_PLACE:
      switch (pmigw.mode) {
      case PrecisionModel::DOUBLE:
        {
          double2* freq_ptr = dp_frequency_data.data(operating_tier);
          fft_operations.emplace_back(&pmigw.ddata[q_data_pos], &freq_ptr[frq_ofs], fft_op_tier,
                                      Normalization::YES, pmigw.fftm, ngx_i, ngy_i, ngz_i, 0, nset,
                                      gstride, frq_stride);
        }
        break;
      case PrecisionModel::SINGLE:
        {
          float2* freq_ptr = sp_frequency_data.data(operating_tier);
          fft_operations.emplace_back(&pmigw.fdata[q_data_pos], &freq_ptr[frq_ofs], fft_op_tier,
                                      Normalization::YES, pmigw.fftm, ngx_i, ngy_i, ngz_i, 0, nset,
                                      gstride, frq_stride);
        }
        break;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConvolutionManager::makeWorkUnits(const GpuDetails &gpu) {

  // Assume that each GPU streaming multiprocessor will support four blocks.  The exact nature of
  // the GPU may not be known at the time the object is created, but the size of the GPU will
  // serve as a reasonable estimate of how large each work unit should be.  It is best to spend as
  // much time as possible on any given system, to maximize the amount of compounding that can
  // happen during energy computations prior to fixed-precision conversion.
  std::vector<uint> fftg_recip_elements(fft_group_count);
  llint total_recip_elements = 0;
  for (int i = 0; i < fft_group_count; i++) {
    const uint idim_a = fft_operations[i].getFrequencyCount(0);
    const uint idim_b = fft_operations[i].getFrequencyCount(1);
    const uint idim_c = fft_operations[i].getFrequencyCount(2);
    const uint nbatch = fft_operations[i].getBatchCount();
    fftg_recip_elements[i] = idim_a * idim_b * idim_c;
    total_recip_elements += fftg_recip_elements[i] * nbatch;
  }

  // The number of work units should reflect the overall size of the synthesis of particle-mesh
  // interaction grids, divide by the number of blocks that will operate on them.
  llint wu_recip_elements = total_recip_elements / static_cast<llint>(gpu.getSMPCount() * 4);
  wu_recip_elements = roundUp<llint>(wu_recip_elements, small_block_size);
  int wu_count = 0;
  for (int i = 0; i < fft_group_count; i++) {
    wu_count += (((fftg_recip_elements[i]) + wu_recip_elements - 1) / wu_recip_elements) *
                fft_operations[i].getBatchCount();
  }
  std::vector<uint4> tmp_wu;
  tmp_wu.reserve(wu_count);
  for (int i = 0; i < fft_group_count; i++) {
    const uint idim_a = fft_operations[i].getFrequencyCount(0);
    const uint idim_b = fft_operations[i].getFrequencyCount(1);
    const uint idim_c = fft_operations[i].getFrequencyCount(2);
    const uint ivol = idim_a * idim_b * idim_c;
    for (int j = fft_group_bounds[i]; j < fft_group_bounds[i + 1]; j++) {
      uint gpos = 0;
      while (gpos < ivol) {
        uint ij_len = (gpos + wu_recip_elements <= ivol) ? wu_recip_elements : ivol - gpos;
        uint4 ij_wu = { static_cast<uint>(fft_groups[j]), gpos, gpos + ij_len, 0 };
        tmp_wu.push_back(ij_wu);
        gpos += ij_len;
      }
    }
  }
  work_units.resize(tmp_wu.size());
  work_units.putHost(tmp_wu, 0, tmp_wu.size());
  hpc_work_unit_count = tmp_wu.size();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis *poly_ps,
                                     const PMIGrid *pmig, ScoreCard *sc) {
  const PsSynthesisBorders pssb = poly_ps->borders();
  const PMIGridReader pmigr = pmig->data();
  std::vector<double> result;
  switch (pmig->getMode()) {
  case PrecisionModel::DOUBLE:
    {
      ConvolutionWriter<double, double2> cvolw = cvol->dpData();
      if (sc == nullptr) {
        result = applyConvolution(&cvolw, pssb, pmigr, nullptr);
      }
      else {
        ScoreCardWriter scw = sc->data();
        result = applyConvolution(&cvolw, pssb, pmigr, &scw);
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    {
      ConvolutionWriter<float, float2> cvolw = cvol->spData();
      if (sc == nullptr) {
        result = applyConvolution(&cvolw, pssb, pmigr, nullptr);
      }
      else {
        ScoreCardWriter scw = sc->data();
        result = applyConvolution(&cvolw, pssb, pmigr, &scw);
      }
      break;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis &poly_ps,
                                     const PMIGrid &pmig, ScoreCard *sc) {
  return applyConvolution(cvol, poly_ps.getSelfPointer(), pmig.getSelfPointer(), sc);
}

} // namespace energy
} // namespace stormm
