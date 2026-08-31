// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
ConvolutionWriter<T, T2>::ConvolutionWriter(const int system_count_in, const int wu_count_in,
                                            const T ew_coeff_in, const int* prf_offsets_in,
                                            const T* self_ecorr_in, const T* bmesh_a_in,
                                            const T* bmesh_b_in, const T* bmesh_c_in,
                                            const T* mval_a_in, const T* mval_b_in,
                                            const T* mval_c_in, const T* msval_a_in,
                                            const T* msval_b_in, const T* msval_c_in,
                                            const T* cmesh_a_in, const T* cmesh_b_in,
                                            const T* cmesh_c_in, T2* freq_data_in,
                                            const uint* frq_offsets_in, const uint4* wu_list_in,
                                            std::vector<FFTStage> *fft_ops_in) :
    system_count{system_count_in}, wu_count{wu_count_in}, ew_coeff{ew_coeff_in},
    prf_offsets{prf_offsets_in}, self_ecorr{self_ecorr_in}, bmesh_a{bmesh_a_in},
    bmesh_b{bmesh_b_in}, bmesh_c{bmesh_c_in}, mval_a{mval_a_in}, mval_b{mval_b_in},
    mval_c{mval_c_in}, msval_a{msval_a_in}, msval_b{msval_b_in}, msval_c{msval_c_in},
    cmesh_a{cmesh_a_in}, cmesh_b{cmesh_b_in}, cmesh_c{cmesh_c_in}, freq_data{freq_data_in},
    frq_offsets{frq_offsets_in}, wu_list{wu_list_in}, fft_ops{fft_ops_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
ConvolutionReader<T, T2>::ConvolutionReader(const int system_count_in, const T ew_coeff_in,
                                            const int* prf_offsets_in, const T* self_ecorr_in,
                                            const T* bmesh_a_in, const T* bmesh_b_in,
                                            const T* bmesh_c_in, const T* mval_a_in,
                                            const T* mval_b_in, const T* mval_c_in,
                                            const T* msval_a_in, const T* msval_b_in,
                                            const T* msval_c_in, const T* cmesh_a_in,
                                            const T* cmesh_b_in, const T* cmesh_c_in,
                                            const T2* freq_data_in, const uint* frq_offsets_in,
                                            const std::vector<FFTStage> *fft_ops_in) :
    system_count{system_count_in}, ew_coeff{ew_coeff_in}, prf_offsets{prf_offsets_in},
    self_ecorr{self_ecorr_in}, bmesh_a{bmesh_a_in}, bmesh_b{bmesh_b_in}, bmesh_c{bmesh_c_in},
    mval_a{mval_a_in}, mval_b{mval_b_in}, mval_c{mval_c_in}, msval_a{msval_a_in},
    msval_b{msval_b_in}, msval_c{msval_c_in}, cmesh_a{cmesh_a_in}, cmesh_b{cmesh_b_in},
    cmesh_c{cmesh_c_in}, freq_data{freq_data_in}, frq_offsets{frq_offsets_in}, fft_ops{fft_ops_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
ConvolutionReader<T, T2>::ConvolutionReader(const ConvolutionWriter<T, T2> *w) :
    system_count{w->system_count}, ew_coeff{w->ew_coeff}, prf_offsets{w->prf_offsets},
    self_ecorr{w->self_ecorr}, bmesh_a{w->bmesh_a}, bmesh_b{w->bmesh_b}, bmesh_c{w->bmesh_c},
    mval_a{w->mval_a}, mval_b{w->mval_b}, mval_c{w->mval_c}, msval_a{w->msval_a},
    msval_b{w->msval_b}, msval_c{w->msval_c}, cmesh_a{w->cmesh_a}, cmesh_b{w->cmesh_b},
    cmesh_c{w->cmesh_c}, freq_data{w->freq_data}, frq_offsets{w->frq_offsets}, fft_ops{w->fft_ops}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
ConvolutionReader<T, T2>::ConvolutionReader(const ConvolutionWriter<T, T2> &w) :
    system_count{w.system_count}, ew_coeff{w.ew_coeff}, prf_offsets{w.prf_offsets},
    self_ecorr{w.self_ecorr}, bmesh_a{w.bmesh_a}, bmesh_b{w.bmesh_b}, bmesh_c{w.bmesh_c},
    mval_a{w.mval_a}, mval_b{w.mval_b}, mval_c{w.mval_c}, msval_a{w.msval_a},
    msval_b{w.msval_b}, msval_c{w.msval_c}, cmesh_a{w.cmesh_a}, cmesh_b{w.cmesh_b},
    cmesh_c{w.cmesh_c}, freq_data{w.freq_data}, frq_offsets{w.frq_offsets}, fft_ops{w.fft_ops}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
const ConvolutionReader<T, T2> ConvolutionManager::data() const {

  // To get automatic templating on the abstract, detect the templated type and make explicit
  // recasts of the proper pointers to the type which they already possess.
  const HybridTargetLevel op_t = operating_tier;
  if (std::type_index(typeid(T)).hash_code() == double_type_index) {
    return ConvolutionReader<T, T2>(system_count, ewald_coefficient, prefactor_offsets.data(op_t),
                                    reinterpret_cast<T*>(self_ecorr.data(op_t)),
                                    reinterpret_cast<T*>(b_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(b_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(b_prefactor_c.data(op_t)),
                                    reinterpret_cast<T*>(m_values_a.data(op_t)),
                                    reinterpret_cast<T*>(m_values_b.data(op_t)),
                                    reinterpret_cast<T*>(m_values_c.data(op_t)),
                                    reinterpret_cast<T*>(mshift_values_a.data(op_t)),
                                    reinterpret_cast<T*>(mshift_values_b.data(op_t)),
                                    reinterpret_cast<T*>(mshift_values_c.data(op_t)),
                                    reinterpret_cast<T*>(c_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(c_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(c_prefactor_c.data(op_t)),
                                    reinterpret_cast<T2*>(dp_frequency_data.data(op_t)),
                                    frequency_offsets.data(op_t), &fft_operations);
  }
  else if (std::type_index(typeid(T)).hash_code() == float_type_index) {
    return ConvolutionReader<T, T2>(system_count, ewald_coefficient, prefactor_offsets.data(op_t),
                                    reinterpret_cast<T*>(sp_self_ecorr.data(op_t)),
                                    reinterpret_cast<T*>(sp_b_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_b_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_b_prefactor_c.data(op_t)),
                                    reinterpret_cast<T*>(sp_m_values_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_m_values_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_m_values_c.data(op_t)),
                                    reinterpret_cast<T*>(sp_mshift_values_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_mshift_values_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_mshift_values_c.data(op_t)),
                                    reinterpret_cast<T*>(sp_c_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_c_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_c_prefactor_c.data(op_t)),
                                    reinterpret_cast<T2*>(sp_frequency_data.data(op_t)),
                                    frequency_offsets.data(op_t), &fft_operations);
  }
  else {
    const std::string t_name = getStormmTypeName<T>();
    rtErr("Type name " + t_name + " is invalid for an abstract of this class.",
          "ConvolutionManager", "data");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
ConvolutionWriter<T, T2> ConvolutionManager::data() {

  // To get automatic templating on the abstract, detect the templated type and make explicit
  // recasts of the proper pointers to the type which they already possess.
  const HybridTargetLevel op_t = operating_tier;
  if (std::type_index(typeid(T)).hash_code() == double_type_index) {
    return ConvolutionWriter<T, T2>(system_count, hpc_work_unit_count, ewald_coefficient,
                                    prefactor_offsets.data(op_t),
                                    reinterpret_cast<T*>(self_ecorr.data(op_t)),
                                    reinterpret_cast<T*>(b_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(b_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(b_prefactor_c.data(op_t)),
                                    reinterpret_cast<T*>(m_values_a.data(op_t)),
                                    reinterpret_cast<T*>(m_values_b.data(op_t)),
                                    reinterpret_cast<T*>(m_values_c.data(op_t)),
                                    reinterpret_cast<T*>(mshift_values_a.data(op_t)),
                                    reinterpret_cast<T*>(mshift_values_b.data(op_t)),
                                    reinterpret_cast<T*>(mshift_values_c.data(op_t)),
                                    reinterpret_cast<T*>(c_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(c_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(c_prefactor_c.data(op_t)),
                                    reinterpret_cast<T2*>(dp_frequency_data.data(op_t)),
                                    frequency_offsets.data(op_t), work_units.data(op_t),
                                    &fft_operations);
  }
  else if (std::type_index(typeid(T)).hash_code() == float_type_index) {
    return ConvolutionWriter<T, T2>(system_count, hpc_work_unit_count, ewald_coefficient,
                                    prefactor_offsets.data(op_t),
                                    reinterpret_cast<T*>(sp_self_ecorr.data(op_t)),
                                    reinterpret_cast<T*>(sp_b_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_b_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_b_prefactor_c.data(op_t)),
                                    reinterpret_cast<T*>(sp_m_values_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_m_values_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_m_values_c.data(op_t)),
                                    reinterpret_cast<T*>(sp_mshift_values_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_mshift_values_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_mshift_values_c.data(op_t)),
                                    reinterpret_cast<T*>(sp_c_prefactor_a.data(op_t)),
                                    reinterpret_cast<T*>(sp_c_prefactor_b.data(op_t)),
                                    reinterpret_cast<T*>(sp_c_prefactor_c.data(op_t)),
                                    reinterpret_cast<T2*>(sp_frequency_data.data(op_t)),
                                    frequency_offsets.data(op_t), work_units.data(op_t),
                                    &fft_operations);
  }
  else {
    const std::string t_name = getStormmTypeName<T>();
    rtErr("Type name " + t_name + " is invalid for an abstract of this class.",
          "ConvolutionManager", "data");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
std::vector<double> pmeGreensFunction(ConvolutionWriter<T, T2> *cvolw,
                                      const PsSynthesisBorders &pssb,
                                      const PMIGridReader &pmigr, ScoreCardWriter *scw) {
  const int xfrm_stride = roundUp(warp_size_int, 32);
  std::vector<double> result(cvolw->system_count, 0.0);
  for (int pos = 0; pos < cvolw->system_count; pos++) {

    // Compute the instantaneous system volume from the inverse transformation matrix.  This
    // quantity is missing even from pre-calculated "C" mesh prefactors.
    const T grid_vol = pmigr.dims[pos].x * pmigr.dims[pos].y * pmigr.dims[pos].z;
    const T pivol = grid_vol / 
                    (pi * pssb.invu[pos * xfrm_stride] * pssb.invu[(pos * xfrm_stride) + 4] *
                     pssb.invu[(pos * xfrm_stride) + 8]);
    const uint npts_a = (pmigr.dims[pos].x / 2) + 1;
    const uint npts_b = pmigr.dims[pos].y;
    const uint npts_c = pmigr.dims[pos].z;
    const uint frq_ofs = cvolw->frq_offsets[pos];
    const T tps = (pi / cvolw->ew_coeff) * (pi / cvolw->ew_coeff);
    const uint prf_ofs = cvolw->prf_offsets[pos];
    for (uint k = 0; k < npts_c; k++) {
      const T dz = cvolw->mval_c[k + prf_ofs];
      const T b_k = cvolw->bmesh_c[k + prf_ofs];
      for (uint j = 0; j < npts_b; j++) {
        const T dy = cvolw->mval_b[j + prf_ofs];
        const T b_jk = b_k * cvolw->bmesh_b[j + prf_ofs];
        for (uint i = 0; i < npts_a; i++) {
          const T dx = cvolw->mval_a[i + prf_ofs];
          const T b_ijk = b_jk * cvolw->bmesh_a[i + prf_ofs];
          const T mmx = (pssb.umat[(pos * xfrm_stride)    ] * dx) +
                        (pssb.umat[(pos * xfrm_stride) + 1] * dy) +
                        (pssb.umat[(pos * xfrm_stride) + 2] * dz);
          const T mmy = (pssb.umat[(pos * xfrm_stride) + 3] * dx) +
                        (pssb.umat[(pos * xfrm_stride) + 4] * dy) +
                        (pssb.umat[(pos * xfrm_stride) + 5] * dz);
          const T mmz = (pssb.umat[(pos * xfrm_stride) + 6] * dx) +
                        (pssb.umat[(pos * xfrm_stride) + 7] * dy) +
                        (pssb.umat[(pos * xfrm_stride) + 8] * dz);
          const T mm_sq = (mmx * mmx) + (mmy * mmy) + (mmz * mmz);
          const T c_ijk = (i == 0 && j == 0 && k == 0) ? 0.0 : pivol * exp(-tps * mm_sq) / mm_sq;
          const uint ijk_idx = frq_ofs + (((k * npts_b) + j) * npts_a) + i;
          T2 fval = cvolw->freq_data[ijk_idx];
          const T qfac = (fval.x * fval.x) + (fval.y * fval.y);
          if (i > 0) {
            result[pos] += qfac * b_ijk * c_ijk;
          }
          else {
            result[pos] += 0.5 * qfac * b_ijk * c_ijk;
          }
          fval.x *= b_ijk * c_ijk;
          fval.y *= b_ijk * c_ijk;
          cvolw->freq_data[ijk_idx] = fval;
        }
      }
    }

    // The grid is normalized here, not during the backward FFT.  Taking this approach makes the
    // convolution simpler to perform on the GPU without an extra kernel call or reliance on a
    // CUDA-specific FFT post-operation feature.
    result[pos] /= grid_vol;
  }
  if (scw != nullptr) {
    switch (pmigr.theme) {
    case NonbondedTheme::ELECTROSTATIC:
      for (int pos = 0; pos < cvolw->system_count; pos++) {
        add(scw, StateVariable::ELECTROSTATIC, static_cast<llint>(result[pos] * scw->nrg_scale_lf),
            pos);
      }
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      for (int pos = 0; pos < cvolw->system_count; pos++) {
        add(scw, StateVariable::VDW, static_cast<llint>(result[pos] * scw->nrg_scale_lf), pos);
      }
      break;
    case NonbondedTheme::ALL:
      
      // This nonsensical case will have been trapped already
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
std::vector<double> applyConvolution(ConvolutionWriter<T, T2> *cvolw,
                                     const PsSynthesisBorders &pssb, const PMIGridReader &pmigr,
                                     ScoreCardWriter *scw) {
  const size_t n_fft = cvolw->fft_ops->size();
  if (n_fft == 0) {
    rtErr("No FFT groups were found in the convolution(s) for " +
          std::to_string(cvolw->system_count) + " systems.\n", "applyConvolution");
  }
  for (size_t i = 0; i < n_fft; i++) {
    cvolw->fft_ops->at(i).forwardFFT();
  }
  std::vector<double> result;
  switch (cvolw->fft_ops->at(0).getTier()) {
  case HybridTargetLevel::HOST:
    result = pmeGreensFunction<T, T2>(cvolw, pssb, pmigr, scw);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }
  for (size_t i = 0; i < n_fft; i++) {
    cvolw->fft_ops->at(i).backwardFFT();
  }
  return result;
}

} // namespace energy
} // namespace stormm
