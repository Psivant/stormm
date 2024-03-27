// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T> BackgroundMeshWriter<T> restoreType(BackgroundMeshWriter<void> *rasa) {
  T* coeff_tmp = reinterpret_cast<T*>(rasa->coeffs);
  return BackgroundMeshWriter<T>(rasa->dims, rasa->kind, rasa->field, rasa->rlrs, coeff_tmp,
                                 rasa->coeff_scale, rasa->probe_radius, rasa->well_depth,
                                 rasa->occ_cost, rasa->mbss);
}

//-------------------------------------------------------------------------------------------------
template <typename T> BackgroundMeshWriter<T> restoreType(const BackgroundMeshWriter<void> &rasa) {
  T* coeff_tmp = reinterpret_cast<T*>(const_cast<void*>(rasa.coeffs));
  return BackgroundMeshWriter<T>(rasa.dims, rasa.kind, rasa.field, rasa.rlrs, coeff_tmp,
                                 rasa.coeff_scale, rasa.probe_radius, rasa.well_depth,
                                 rasa.occ_cost, rasa.mbss);
}

//-------------------------------------------------------------------------------------------------
template <typename T> BackgroundMeshReader<T> restoreType(const BackgroundMeshReader<void> *rasa) {
  const T* coeff_tmp = reinterpret_cast<const T*>(rasa->coeffs);
  return BackgroundMeshReader<T>(rasa->dims, rasa->kind, rasa->field, rasa->rlrs, coeff_tmp,
                                 rasa->coeff_scale, rasa->probe_radius, rasa->well_depth,
                                 rasa->occ_cost, rasa->mbss);
}

//-------------------------------------------------------------------------------------------------
template <typename T> BackgroundMeshReader<T> restoreType(const BackgroundMeshReader<void> &rasa) {
  const T* coeff_tmp = reinterpret_cast<const T*>(rasa.coeffs);
  return BackgroundMeshReader<T>(rasa.dims, rasa.kind, rasa.field, rasa.rlrs, coeff_tmp,
                                 rasa.coeff_scale, rasa.probe_radius, rasa.well_depth,
                                 rasa.occ_cost, rasa.mbss);
}

//-------------------------------------------------------------------------------------------------
template <typename Ttarget, typename Tcontrib>
void checkMeshCompatibility(const BackgroundMeshWriter<Ttarget> &target,
                            const BackgroundMeshReader<Tcontrib> &contrib) {
  if (target.dims.bounds != contrib.dims.bounds) {
    rtErr("Boundary conditions between the accumulator and contribution do not match.",
          "accumulateOcclusionMesh");
  }
  std::vector<double> target_invu(9), contrib_invu(9);
  const double target_na = target.dims.na;
  const double target_nb = target.dims.nb;
  const double target_nc = target.dims.nc;
  const double contrib_na = contrib.dims.na;
  const double contrib_nb = contrib.dims.nb;
  const double contrib_nc = contrib.dims.nc;
  for (int i = 0; i < 3; i++) {
    target_invu[i    ] = target.dims.invu[i    ] * target_na;
    target_invu[i + 3] = target.dims.invu[i + 3] * target_nb;
    target_invu[i + 6] = target.dims.invu[i + 6] * target_nc;
    contrib_invu[i    ] = contrib.dims.invu[i    ] * contrib_na;
    contrib_invu[i + 3] = contrib.dims.invu[i + 3] * contrib_nb;
    contrib_invu[i + 6] = contrib.dims.invu[i + 6] * contrib_nc;
  }
  switch (target.dims.bounds) {
  case BoundaryCondition::PERIODIC:

    // Both meshes must have the same overall size
    if (maxAbsoluteDifference(target_invu, contrib_invu) > 1.0e-5) {
      rtErr("Periodic mesh dimensions must be identical in order to accumulate one into the "
            "other.  Maximum difference in dimensions was " +
            realToString(maxAbsoluteDifference(target_invu, contrib_invu), 9, 5,
                         NumberFormat::STANDARD_REAL) + ".", "checkMeshCompatibility");
    }
    break;
  case BoundaryCondition::ISOLATED:

    // The meshes need not have the same overall size, but the contributing mesh must be at least
    // the size of the target, and no more than two of its elements larger in any direction.
    if (maxAbsoluteDifference(target_invu, contrib_invu) > contrib.dims.max_span + 1.0e-5) {
      rtErr("For isolated meshes, the sizes do not need to align perfectly, but a contributing "
            "mesh that is significantly larger than the accumulating mesh is invalid.",
            "checkMeshCompatibility");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evalElecCSC(const Tcalc4 abcd, const Tcalc r, const Tcalc rswitch, const Tcalc q) {
  if (r < rswitch) {
    Tcalc2 result = evaluateCubicSpline<Tcalc4, Tcalc2, Tcalc>(abcd, r);
    result.x *= q;
    result.y *= q;
    return result;
  }
  else {
    return { q / r, -q / (r * r) };
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc evalElecCSCND(const Tcalc4 abcd, const Tcalc r, const Tcalc rswitch, const Tcalc q) {
  if (r < rswitch) {
    const Tcalc2 result = evaluateCubicSpline<Tcalc4, Tcalc2, Tcalc>(abcd, r);
    return q * result.x;
  }
  else {
    return q / r;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evalLennardJonesCSC(const Tcalc4 abcd, const Tcalc r, const Tcalc rswitch, const Tcalc lja,
                           const Tcalc ljb) {
  if (r < rswitch) {
    return evaluateCubicSpline<Tcalc4, Tcalc2, Tcalc>(abcd, r);
  }
  else {
    Tcalc2 result;
    const Tcalc invr = static_cast<Tcalc>(1.0) / r;
    const Tcalc invr3 = invr * invr * invr;
    const Tcalc invr6 = invr3 * invr3;
    return { ((lja * invr6) - ljb) * invr6,
             ((static_cast<Tcalc>(-12.0) * lja * invr6) +
              (static_cast<Tcalc>(6.0) * ljb)) * invr6 * invr };
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc evalLennardJonesCSCND(const Tcalc4 abcd, const Tcalc r, const Tcalc rswitch, const Tcalc lja,
                            const Tcalc ljb) {
  if (r < rswitch) {
    const Tcalc2 result = evaluateCubicSpline<Tcalc4, Tcalc2, Tcalc>(abcd, r);
    return result.x;
  }
  else {
    Tcalc2 result;
    const Tcalc invr3 = static_cast<Tcalc>(1.0) / (r * r * r);
    const Tcalc invr6 = invr3 * invr3;
    return ((lja * invr6) - ljb) * invr6;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void accumulateNonbondedFieldMesh(BackgroundMeshWriter<llint> *target,
                                  const BackgroundMeshReader<T> &contribution,
                                  const GpuDetails &gpu) {

  // If a GPU is available and the data is available on the device, launch a kernel to compute the
  // derivatives on the GPU.  It is not practical to port the numbers to and from the GPU for what
  // is at most a conversion and accumulation for each value.
  checkMeshCompatibility(target, contribution);

  // Loop over all mesh points and perform the accumulation.
  const size_t ncoeff = static_cast<size_t>(target->dims.na * target->dims.nb * target->dims.nc) *
                        64LLU;
  llint* coeff_ptr = target->coeffs;
  if (isFloatingPointScalarType<T>()) {
    for (size_t i = 0; i < ncoeff; i++) {
      const llint tc = llround(contribution.coeffs[i] * target->coeff_scale);
      coeff_ptr[i] += tc;
    }
  }
  else {
    const int contrib_bits = round(log2(contribution.coeff_scale));
    const int target_bits = round(log2(target->coeff_scale));
    if (contrib_bits > target_bits) {
      const int bit_diff = contrib_bits - target_bits;
      for (size_t i = 0; i < ncoeff; i++) {
        const llint tc = contribution.coeffs[i];
        coeff_ptr[i] += (tc >> bit_diff);
      }
    }
    else if (contrib_bits < target_bits) {
      const int bit_diff = target_bits - contrib_bits;
      for (size_t i = 0; i < ncoeff; i++) {
        const llint tc = contribution.coeffs[i];
        coeff_ptr[i] += (tc << bit_diff);
      }
    }
    else {
      for (size_t i = 0; i < ncoeff; i++) {
        coeff_ptr[i] += contribution.coeffs[i];
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void accumulateNonbondedFieldMesh(BackgroundMesh<llint> *target,
                                  const BackgroundMesh<T> &contribution,
                                  const GpuDetails &gpu, HybridTargetLevel availability) {
  BackgroundMeshWriter targetw = target->data(availability);
  const BackgroundMeshReader<T> contribr = contribution.data(availability);
#ifdef STORMM_USE_HPC
  const size_t ct_contrib = std::type_index(typeid(T)).hash_code();
  if (gpu != null_gpu) {
    switch (availability) {
    case HybridTargetLevel::HOST:
      break;
    case HybridTargetLevel::DEVICE:
      launchAccNonbondedFieldMesh(&target, contribution, ct_contrib, gpu);
      return;
    }
  }
#endif
  accumulateNonbondedFieldMesh(targetw, contribr, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void inferPatchBoundary(std::vector<T> *u_patch, const int dim_a, const int dim_b, const int dim_c,
                        const int start_a, const int start_b, const int start_c, const int end_a,
                        const int end_b, const int end_c, const int inc_a, const int inc_b,
                        const int inc_c) {
  const T value_two = 2.0;
  T* u_data = u_patch->data();
  for (int i = start_a; i < end_a; i += inc_a) {
    int ifulcrum, ilever;
    if (i == 0) {
      ifulcrum = 1;
      ilever   = 2;
    }
    else if (i == dim_a - 1) {
      ifulcrum = dim_a - 2;
      ilever   = dim_a - 3;
    }
    else {
      ifulcrum = i;
      ilever   = i;
    }
    for (int j = start_b; j < end_b; j += inc_b) {
      int jfulcrum, jlever;
      if (j == 0) {
        jfulcrum = 1;
        jlever   = 2;
      }
      else if (j == dim_b - 1) {
        jfulcrum = dim_b - 2;
        jlever   = dim_b - 3;
      }
      else {
        jfulcrum = j;
        jlever   = j;
      }
      for (int k = start_c; k < end_c; k += inc_c) {
        int kfulcrum, klever;
        if (k == 0) {
          kfulcrum = 1;
          klever   = 2;
        }
        else if (k == dim_c - 1) {
          kfulcrum = dim_c - 2;
          klever   = dim_c - 3;
        }
        else {
          kfulcrum = k;
          klever   = k;
        }
        const T fulcrum = u_data[(((kfulcrum * dim_b) + jfulcrum) * dim_a) + ifulcrum];
        const T lever   = u_data[(((klever * dim_b) + jlever) * dim_a) + ilever];
        u_data[(((k * dim_b) + j) * dim_a) + i] = (value_two * fulcrum) - lever;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transferPatchResult(const std::vector<T> &patch, const int patch_dim_a, const int patch_dim_b,
                         const int patch_dim_c, BackgroundMeshWriter<T> *occfieldw,
                         const int a_index, const int b_index, const int c_index,
                         const int readout_a_start, const int readout_b_start,
                         const int readout_c_start, const int element_offset) {
  const size_t sixty_four = 64;
  const int g_na = occfieldw->dims.na;
  const int g_nb = occfieldw->dims.nb;
  const int g_nc = occfieldw->dims.nc;
  for (int i = readout_a_start; i < patch_dim_a - 1; i++) {
    const int output_a_idx = i - readout_a_start + a_index;
    for (int j = readout_b_start; j < patch_dim_b - 1; j++) {
      const int output_b_idx = j - readout_b_start + b_index;
      for (int k = readout_c_start; k < patch_dim_c - 1; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        const int output_c_idx = k - readout_c_start + c_index;

        // There are up to eight places that the output could be applicable.  Test each of them.
        for (int m = 0; m < 8; m++) {
          const size_t m_offset = m + element_offset;
          const int mk = m / 4;
          const int mj = (m - (4 * mk)) / 2;
          const int mi = (m & 1);
          int shift_a = output_a_idx - mi;
          int shift_b = output_b_idx - mj;
          int shift_c = output_c_idx - mk;
          switch (occfieldw->dims.bounds) {
          case BoundaryCondition::PERIODIC:
            {
              shift_a += g_na * ((shift_a < 0) - (shift_a >= g_na));
              shift_b += g_nb * ((shift_b < 0) - (shift_b >= g_nb));
              shift_c += g_nc * ((shift_c < 0) - (shift_c >= g_nc));
              const size_t element_idx = (((shift_c * g_nb) + shift_b) * g_na) + shift_a;
              occfieldw->coeffs[(element_idx * sixty_four) + m_offset] = patch[ijk_idx];
            }
            break;
          case BoundaryCondition::ISOLATED:
            if (shift_a >= 0 && shift_a < g_na && shift_b >= 0 && shift_b < g_nb &&
                shift_c >= 0 && shift_c < g_nc) {
              const size_t element_idx = (((shift_c * g_nb) + shift_b) * g_na) + shift_a;
              occfieldw->coeffs[(element_idx * sixty_four) + m_offset] = patch[ijk_idx];
            }
            break;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void patchNumericDeriv(BackgroundMeshWriter<T> *occfieldw, const int a_index, const int b_index,
                       const int c_index, const T max_occlusion) {

  // Compute the dimensions of the mesh patch.  The patch will be up to eight mesh points in all
  // directions, based on a need to compute derivatives for (na + 1) by (nb + 1) by (nc + 1)
  // points for a mesh of na by nb by nc elements.  Furthermore, in isolated boundary conditions,
  // the patch must include at least two elements in all directions so that the boundary layer can
  // be determined.  In such cases, the second point can be taken from the original mesh element.
  const int g_na = occfieldw->dims.na;
  const int g_nb = occfieldw->dims.nb;
  const int g_nc = occfieldw->dims.nc;
  int actual_a_index = a_index;
  int actual_b_index = b_index;
  int actual_c_index = c_index;
  int bound_offset;
  switch (occfieldw->dims.bounds) {
  case BoundaryCondition::ISOLATED:
    bound_offset = 1;
    actual_a_index -= (a_index == g_na);
    actual_b_index -= (b_index == g_nb);
    actual_c_index -= (c_index == g_nc);
    break;
  case BoundaryCondition::PERIODIC:
    bound_offset = 0;
    break;
  }
  const int patch_dim_a = std::min(g_na + bound_offset - actual_a_index, 8) + 2;
  const int patch_dim_b = std::min(g_nb + bound_offset - actual_b_index, 8) + 2;
  const int patch_dim_c = std::min(g_nc + bound_offset - actual_c_index, 8) + 2;
  const size_t sixty_four = 64;
  const int patch_size = patch_dim_a * patch_dim_b * patch_dim_c;
  std::vector<T> u_patch(patch_size), du_patch(patch_size), ddu_patch(patch_size);
  for (int i = 0; i < patch_dim_a; i++) {
    int gidx_a = actual_a_index - 1 + i;
    bool infer_a = false;
    switch (occfieldw->dims.bounds) {
    case BoundaryCondition::ISOLATED:
      infer_a = (gidx_a < 0 || gidx_a > g_na);
      break;
    case BoundaryCondition::PERIODIC:
      gidx_a += ((gidx_a < 0) - (gidx_a >= g_na)) * g_na;
      break;
    }
    for (int j = 0; j < patch_dim_b; j++) {
      int gidx_b = actual_b_index - 1 + j;
      bool infer_b = false;
      switch (occfieldw->dims.bounds) {
      case BoundaryCondition::ISOLATED:
        infer_b = (gidx_b < 0 || gidx_b > g_nb);
        break;
      case BoundaryCondition::PERIODIC:
        gidx_b += ((gidx_b < 0) - (gidx_b >= g_nb)) * g_nb;
        break;
      }
      for (int k = 0; k < patch_dim_c; k++) {
        int gidx_c = actual_c_index - 1 + k;
        bool infer_c = false;
        switch (occfieldw->dims.bounds) {
        case BoundaryCondition::ISOLATED:
          infer_c = (gidx_c < 0 || gidx_c > g_nc);
          break;
        case BoundaryCondition::PERIODIC:
          gidx_c += ((gidx_c < 0) - (gidx_c >= g_nc)) * g_nc;
          break;
        }

        // Select, or compute, the potentials for the mesh patch.
        if (!(infer_a) && !(infer_b) && !(infer_c)) {
          switch (occfieldw->dims.bounds) {
          case BoundaryCondition::ISOLATED:
            {
              const size_t gelem_idx = (sixty_four * (((((gidx_c - (gidx_c == g_nc)) * g_nb) +
                                                        (gidx_b - (gidx_b == g_nb))) * g_na) +
                                                      (gidx_a - (gidx_a == g_na)))) +
                                       ((((gidx_c == g_nc) * 2) + (gidx_b == g_nb)) * 2) +
                                       (gidx_a == g_na);
              u_patch[(((k * patch_dim_b) + j) * patch_dim_a) + i] = occfieldw->coeffs[gelem_idx];
            }
            break;
          case BoundaryCondition::PERIODIC:
            {
              const size_t gelem_idx = (((gidx_c * g_nb) + gidx_b) * g_na) + gidx_a;
              u_patch[(((k * patch_dim_b) + j) * patch_dim_a) + i] =
                occfieldw->coeffs[sixty_four * gelem_idx];
            }
            break;
          }
        }
      }
    }
  }

  // Make a second and a third pass, if isolated boundary conditions are in effect, to fill out
  // the final layers of the potential patch.
  switch (occfieldw->dims.bounds) {
  case BoundaryCondition::ISOLATED:
    {
      // One layer around the patch filled from the actual grid must be inferred.  Because there
      // are so many conditionals to evaluate to determine which two points in the patch determine
      // the potential at each point on the boundary, it is best done with separate loops.  Start
      // with the AB slabs on the top and bottom, then do other boundary slabs.
      inferPatchBoundary(&u_patch, patch_dim_a, patch_dim_b, patch_dim_c, 0, 0, 0, patch_dim_a,
                         patch_dim_b, patch_dim_c, 1, 1, patch_dim_c - 1);
      inferPatchBoundary(&u_patch, patch_dim_a, patch_dim_b, patch_dim_c, 0, 0, 1, patch_dim_a,
                         patch_dim_b, patch_dim_c - 1, 1, patch_dim_b - 1, 1);
      inferPatchBoundary(&u_patch, patch_dim_a, patch_dim_b, patch_dim_c, 0, 1, 1, patch_dim_a,
                         patch_dim_b - 1, patch_dim_c - 1, patch_dim_a - 1, 1, 1);
    }
    break;
  case BoundaryCondition::PERIODIC:
    break;
  }

  // With the patch assembled, compute du/da, then d2u/dab, then d3u/dabc.  Store the derivatives
  // as appropriate.  Declare first derivatives to be zero wherever the mesh has no occuancy, or
  // near total occupancy.
  const T max_occ_tol = static_cast<T>(max_occlusion) * 1.0e-6;
  const T value_zero = 0.0;
  const T value_half = 0.5;
  const int readout_a_start = (actual_a_index != a_index) + 1;
  const int readout_b_start = (actual_b_index != b_index) + 1;
  const int readout_c_start = (actual_c_index != c_index) + 1;
  for (int i = 1; i < patch_dim_a - 1; i++) {
    for (int j = 0; j < patch_dim_b; j++) {
      for (int k = 0; k < patch_dim_c; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        if (u_patch[ijk_idx] < max_occlusion - max_occ_tol && u_patch[ijk_idx] > max_occ_tol) {
          const int ijk_mi  = ijk_idx - 1;
          const int ijk_pi  = ijk_idx + 1;
          du_patch[ijk_idx] = (u_patch[ijk_pi] - u_patch[ijk_mi]) * value_half;
        }
        else {
          du_patch[ijk_idx] = value_zero;
        }
      }
    }
  }
  transferPatchResult(&du_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 8);
  for (int i = 1; i < patch_dim_a - 1; i++) {
    for (int j = 1; j < patch_dim_b - 1; j++) {
      for (int k = 0; k < patch_dim_c; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        const int ijk_mj  = (((k * patch_dim_b) + j - 1) * patch_dim_a) + i;
        const int ijk_pj  = (((k * patch_dim_b) + j + 1) * patch_dim_a) + i;

        // Store second derivatives in a secondary work array.  Even when computing on a single
        // thread, it can be hard to find an order that guarantees that the first derivatives one
        // wants to use will not have already been overwritten by the second derivative of some
        // previous calculation.
        ddu_patch[ijk_idx] = (du_patch[ijk_pj] - du_patch[ijk_mj]) * value_half;
      }
    }
  }
  transferPatchResult(&ddu_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 32);
  for (int i = 1; i < patch_dim_a - 1; i++) {
    for (int j = 1; j < patch_dim_b - 1; j++) {
      for (int k = 1; k < patch_dim_c - 1; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        const int ijk_mk  = ((((k - 1) * patch_dim_b) + j) * patch_dim_a) + i;
        const int ijk_pk  = ((((k + 1) * patch_dim_b) + j) * patch_dim_a) + i;

        // Place triple derivatives back in the original work array
        du_patch[ijk_idx] = (ddu_patch[ijk_pk] - ddu_patch[ijk_mk]) * value_half;
      }
    }
  }
  transferPatchResult(&du_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 56);

  // Continue to compute du/db, then d2u/dbc.  Store the derivatives as appropriate.
  for (int i = 1; i < patch_dim_a - 1; i++) {
    for (int j = 1; j < patch_dim_b - 1; j++) {
      for (int k = 0; k < patch_dim_c; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        if (u_patch[ijk_idx] < max_occlusion - max_occ_tol && u_patch[ijk_idx] > max_occ_tol) {
          const int ijk_mj  = (((k * patch_dim_b) + j - 1) * patch_dim_a) + i;
          const int ijk_pj  = (((k * patch_dim_b) + j + 1) * patch_dim_a) + i;
          du_patch[ijk_idx] = (u_patch[ijk_pj] - u_patch[ijk_mj]) * value_half;
        }
      }
    }
  }
  transferPatchResult(&du_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 16);
  for (int i = 1; i < patch_dim_a - 1; i++) {
    for (int j = 1; j < patch_dim_b - 1; j++) {
      for (int k = 1; k < patch_dim_c - 1; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        const int ijk_mk  = ((((k - 1) * patch_dim_b) + j) * patch_dim_a) + i;
        const int ijk_pk  = ((((k + 1) * patch_dim_b) + j) * patch_dim_a) + i;
        ddu_patch[ijk_idx] = (du_patch[ijk_pk] - du_patch[ijk_mk]) * value_half;
      }
    }
  }
  transferPatchResult(&ddu_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 48);

  // Finally, compute du/dc, and based on that d2u/dac.  Store the derivatives as appropriate.
  for (int i = 0; i < patch_dim_a; i++) {
    for (int j = 1; j < patch_dim_b - 1; j++) {
      for (int k = 1; k < patch_dim_c - 1; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        if (u_patch[ijk_idx] < max_occlusion - max_occ_tol && u_patch[ijk_idx] > max_occ_tol) {
          const int ijk_mk  = ((((k - 1) * patch_dim_b) + j) * patch_dim_a) + i;
          const int ijk_pk  = ((((k + 1) * patch_dim_b) + j) * patch_dim_a) + i;
          du_patch[ijk_idx] = (u_patch[ijk_pk] - u_patch[ijk_mk]) * value_half;
        }
      }
    }
  }
  transferPatchResult(&du_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 24);
  for (int i = 1; i < patch_dim_a - 1; i++) {
    for (int j = 1; j < patch_dim_b - 1; j++) {
      for (int k = 1; k < patch_dim_c - 1; k++) {
        const int ijk_idx = (((k * patch_dim_b) + j) * patch_dim_a) + i;
        const int ijk_mi  = (((k * patch_dim_b) + j) * patch_dim_a) + i - 1;
        const int ijk_pi  = (((k * patch_dim_b) + j) * patch_dim_a) + i + 1;
        ddu_patch[ijk_idx] = (du_patch[ijk_pi] - du_patch[ijk_mi]) * value_half;
        
      }
    }
  }  
  transferPatchResult(&ddu_patch, patch_dim_a, patch_dim_b, patch_dim_c, occfieldw, a_index,
                      b_index, c_index, readout_a_start, readout_b_start, readout_c_start, 40);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void occlusionFieldDerivatives(BackgroundMesh<T> *occfield, const T max_occlusion,
                               const GpuDetails &gpu, const HybridTargetLevel availability) {
  const size_t ct_mesh = std::type_index(typeid(T)).hash_code();

  // If a GPU is available and the energy data is already staged on the device, launch a kernel
  // to compute the derivatives on the GPU.
#ifdef STORMM_USE_HPC
  if (gpu != null_gpu) {
    switch (availability) {
    case HybridTargetLevel::HOST:
      break;
    case HybridTargetLevel::DEVICE:
      {
        BackgroundMeshWriter<void> bgmw = occfield->templateFreeData(availability);
        launchOccFieldDerivativeCalc(&bgmw, max_occlusion, ct_mesh, gpu);
      }
      return;
    }
  }
#endif
  
  // Compute derivatives on the host
  BackgroundMeshWriter<T> occfieldw = occfield->data();
  const T value_zero = 0.0;
  for (int i = 0; i < occfieldw.dims.na; i += 8) {
    for (int j = 0; j < occfieldw.dims.nb; j += 8) {
      for (int k = 0; k < occfieldw.dims.nc; k += 8) {
        patchNumericDeriv(&occfieldw, i, j, k, max_occlusion);
      }
    }
  }
}
  
} // namespace structure
} // namespace stormm
