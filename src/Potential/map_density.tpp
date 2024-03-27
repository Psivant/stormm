// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tgrid>
void particleAlignment(const Tcalc x, const Tcalc y, const Tcalc z, const Tcalc lpos_inv_scale,
                       const Tcalc* umat, const int cg_mesh_ticks, const int cell_i,
                       const int cell_j, const int cell_k, Tgrid *a_cof, Tgrid *b_cof,
                       Tgrid *c_cof, const int bspline_order, int *grid_a, int *grid_b,
                       int *grid_c) {
  Tcalc rel_a, rel_b, rel_c;
  if (lpos_inv_scale > 0.99) {
    rel_a = (umat[0] * x) + (umat[3] * y) + (umat[6] * z);
    rel_b =                 (umat[4] * y) + (umat[7] * z);
    rel_c =                                 (umat[8] * z);
  }
  else {
    
    // A static_cast is unecessary even if the Cartesian coordinates are expressed in
    // fixed-precision, as this will have happened in the function argument conversion.
    const Tcalc atm_x = x * lpos_inv_scale;
    const Tcalc atm_y = y * lpos_inv_scale;
    const Tcalc atm_z = z * lpos_inv_scale;
    rel_a = (umat[0] * atm_x) + (umat[3] * atm_y) + (umat[6] * atm_z);
    rel_b =                     (umat[4] * atm_y) + (umat[7] * atm_z);
    rel_c =                                         (umat[8] * atm_z);
  }
  
  // There is a very slight chance that the floating point math may work out such that the
  // grid point for the particle lies at or just beyond the cell boundary.  This is an
  // unintended effect but the consequences will be minimal even if mapping is restricted to
  // the cell and its neighbors up / front / right.
  const Tcalc dmt = cg_mesh_ticks;
  Tcalc da = rel_a * dmt;
  Tcalc db = rel_b * dmt;
  Tcalc dc = rel_c * dmt;
  const Tcalc fl_da = floor(da);
  const Tcalc fl_db = floor(db);
  const Tcalc fl_dc = floor(dc);
  da -= fl_da;
  db -= fl_db;
  dc -= fl_dc;
  const int grid_local_a = fl_da;
  const int grid_local_b = fl_db;
  const int grid_local_c = fl_dc;
  *grid_a = grid_local_a + (cg_mesh_ticks * cell_i);
  *grid_b = grid_local_b + (cg_mesh_ticks * cell_j);
  *grid_c = grid_local_c + (cg_mesh_ticks * cell_k);
  
  // The B-spline computation marks a demarcation between computations in the precision of the
  // cell grid and the precision of the particle-mesh interaction grid, which might not be
  // identical (although they are allowed to differ more for purposes of experimentation in the
  // numerics).
  bSpline<Tgrid>(da, bspline_order, a_cof);
  bSpline<Tgrid>(db, bspline_order, b_cof);
  bSpline<Tgrid>(dc, bspline_order, c_cof);
}

//-------------------------------------------------------------------------------------------------
template <typename Tsrc, typename Tcalc, typename Tcalc2>
Tcalc sourceMagnitude(const NonbondedTheme pmig_density, const NonbondedTheme cg_content,
                      const Tsrc q, const bool q_is_real, const int sysid,
                      const SyNonbondedKit<Tcalc, Tcalc2> &synbk) {

  // The switch over the density theme is only necessary to distinguish what to collect if the
  // cell grid serves a non-bonded theme of "ALL."  The correct source of density to be represented
  // on the particle-mesh interaction grid is found by examining the cell grid's contents and
  // then choosing the proper interpretation of its 
  switch (pmig_density) {
  case NonbondedTheme::ELECTROSTATIC:
    switch (cg_content) {
    case NonbondedTheme::ELECTROSTATIC:

      // The charge will be represented in the cell grid data itself.
      if (q_is_real) {
        return q;
      }
      else {
        if (sizeof(Tsrc) == 8) {
          const Ecumenical8 conv = { .lli = static_cast<llint>(q) };
          return conv.d;
        }
        else {
          const Ecumenical4 conv = { .i = static_cast<int>(q) };
          return conv.f;
        }
      }
      break;
    case NonbondedTheme::VAN_DER_WAALS:

      // This erroneous case is trapped earlier in the code path.
      break;
    case NonbondedTheme::ALL:

      // The source data will necessarily have a signed integer format.
      if (sizeof(Tsrc) == 8) {
        const int q_idx = (static_cast<llint>(q) & dp_charge_index_mask);
        return synbk.q_params[q_idx];
      }
      else {
        const int q_idx = (static_cast<int>(q) & sp_charge_index_mask);
        return synbk.q_params[q_idx];
      }
      break;
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    {
      int tlj_idx;
      switch (cg_content) {
      case NonbondedTheme::ELECTROSTATIC:

        // This erroneous case is trapped earlier in the code path.
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        if (q_is_real) {
          if (sizeof(Tsrc) == 8) {
            const Ecumenical8 conv = { .d = static_cast<double>(q) };
            tlj_idx = conv.lli;
          }
          else {
            const Ecumenical4 conv = { .f = static_cast<float>(q) };
            tlj_idx = conv.i;
          }
        }
        else {
          tlj_idx = q;
        }
        break;
      case NonbondedTheme::ALL:

        // The source data will necessarily have a signed integer format.
        if (sizeof(Tsrc) == 8) {
          tlj_idx = (static_cast<llint>(q) >> dp_charge_index_bits);
        }
        else {
          tlj_idx = (static_cast<int>(q) >> sp_charge_index_bits);
        }
        break;
      }
      const Tcalc t_ljb = synbk.ljb_coeff[synbk.ljabc_offsets[sysid] +
                                          (tlj_idx * (synbk.n_lj_types[sysid] + 1))];
      if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
        return sqrt(0.25 * t_ljb);
      }
      else {
        return sqrtf(0.25f * t_ljb);
      }
    }
    break;
  case NonbondedTheme::ALL:

    // No PMI Grid will take both types of non-bonded source density.
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tgrid>
void spreadDensity(const Tcalc* a_cof, const Tcalc* b_cof, const Tcalc* c_cof,
                   const int bspline_order, const int grid_root_a, const int grid_root_b,
                   const int grid_root_c, const uint4 grid_dims, const FFTMode fft_staging,
                   Tgrid* grid_data, int* overflow, const Tcalc mesh_scaling_factor) {
  const int om = bspline_order - 1;
  const bool tgrid_is_llint = (std::type_index(typeid(Tgrid)).hash_code() == llint_type_index);
  uint padded_gdim_x;
  switch (fft_staging) {
  case FFTMode::IN_PLACE:
    padded_gdim_x = 2 * ((grid_dims.x / 2) + 1);
    break;
  case FFTMode::OUT_OF_PLACE:
    padded_gdim_x = grid_dims.x;
    break;
  }
  for (int k = 0; k < bspline_order; k++) {
    int kg_pos = grid_root_c + k;
    kg_pos += ((kg_pos < 0) - (kg_pos >= grid_dims.z)) * grid_dims.z;
    for (int j = 0; j < bspline_order; j++) {
      int jg_pos = grid_root_b + j;
      jg_pos += ((jg_pos < 0) - (jg_pos >= grid_dims.y)) * grid_dims.y;
      const size_t jk_gidx = grid_dims.w + (((kg_pos * grid_dims.y) + jg_pos) * padded_gdim_x);
      const Tcalc jk_contrib = b_cof[j] * c_cof[k];
      for (int i = 0; i < bspline_order; i++) {
        int ig_pos = grid_root_a + i;
        ig_pos += ((ig_pos < 0) - (ig_pos >= grid_dims.x)) * grid_dims.x;
        const size_t gidx = jk_gidx + static_cast<size_t>(ig_pos);
        if (overflow == nullptr) {
          
          // To have no overflow bits array is taken as an indication that accumulation is to
          // occur in real-valued, floating-point numbers.
          grid_data[gidx] += jk_contrib * a_cof[i];
        }
        else {

          // Further branching between the int95_t and int63_t cases
          if (tgrid_is_llint) {
            llint ival_x = grid_data[gidx];
            int ival_y = overflow[gidx];
            hostSplitAccumulation(jk_contrib * a_cof[i] * mesh_scaling_factor, &ival_x, &ival_y);
            grid_data[gidx] = ival_x;
            overflow[gidx]  = ival_y;
          }
          else {
            int ival_x = grid_data[gidx];
            int ival_y = overflow[gidx];
            hostSplitAccumulation(jk_contrib * a_cof[i] * mesh_scaling_factor, &ival_x, &ival_y);
            grid_data[gidx] = ival_x;
            overflow[gidx]  = ival_y;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename Tcalc2, typename T4>
void accumulateCellDensity(PMIGridWriter *pm_wrt, const int sysid, const int cell_i,
                           const int cell_j, const int cell_k,
                           const CellGridReader<T, Tacc, Tcalc, T4> &cgr,
                           const SyNonbondedKit<Tcalc, Tcalc2> &synbk) {

  // Re-derive the cell and cell grid boundaries rather than copying them via input arguments
  const ullint cell_bounds = cgr.system_cell_grids[sysid];
  const int cell_offset = (cell_bounds & 0xfffffffLLU);
  const int cell_na = ((cell_bounds >> 28) & 0xfffLLU);
  const int cell_nb = ((cell_bounds >> 40) & 0xfffLLU);
  const int cell_nc = ((cell_bounds >> 52) & 0xfffLLU);

  // Determine limits and other critical constants
  const bool coord_in_real = (cgr.lpos_scale < 1.01);
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const size_t ijk_cellidx = cell_offset + (((cell_k * cell_nb) + cell_j) * cell_na) + cell_i;
  const uint2 ijk_bounds = cgr.cell_limits[ijk_cellidx];
  const uint mllim = ijk_bounds.x;
  const uint mhlim = mllim + (ijk_bounds.y >> 16);
  const int xfrm_stride = roundUp(9, warp_size_int);
  Tcalc umat[9];
  for (int i = 0; i < 9; i++) {
    umat[i] = cgr.system_cell_umat[(sysid * xfrm_stride) + i];
  }

  // Lay out arrays to collect B-spline coefficients
  const uint4 grid_dims = pm_wrt->dims[sysid];
  switch (pm_wrt->mode) {
  case PrecisionModel::DOUBLE:
    {
      std::vector<double> a_cof(pm_wrt->order), b_cof(pm_wrt->order), c_cof(pm_wrt->order);
      for (uint m = mllim; m < mhlim; m++) {
        const T4 atom_m = cgr.image[m];
        int grid_root_a, grid_root_b, grid_root_c;
        particleAlignment<Tcalc, double>(atom_m.x, atom_m.y, atom_m.z, cgr.lpos_inv_scale, umat,
                                         cgr.mesh_ticks, cell_i, cell_j, cell_k, a_cof.data(),
                                         b_cof.data(), c_cof.data(), pm_wrt->order, &grid_root_a,
                                         &grid_root_b, &grid_root_c);
        const double q = sourceMagnitude<T, Tcalc, Tcalc2>(pm_wrt->theme, cgr.theme, atom_m.w,
                                                           coord_in_real, sysid, synbk);
        for (int i = 0; i < pm_wrt->order; i++) {
          a_cof[i] *= q;
        }
        spreadDensity<double, double>(a_cof.data(), b_cof.data(), c_cof.data(), pm_wrt->order,
                                      grid_root_a, grid_root_b, grid_root_c, pm_wrt->dims[sysid],
                                      pm_wrt->fftm, pm_wrt->ddata);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      std::vector<float> a_cof(pm_wrt->order), b_cof(pm_wrt->order), c_cof(pm_wrt->order);
      for (uint m = mllim; m < mhlim; m++) {
        const T4 atom_m = cgr.image[m];
        int grid_root_a, grid_root_b, grid_root_c;
        particleAlignment<Tcalc, float>(atom_m.x, atom_m.y, atom_m.z, cgr.lpos_inv_scale, umat,
                                        cgr.mesh_ticks, cell_i, cell_j, cell_k, a_cof.data(),
                                        b_cof.data(), c_cof.data(), pm_wrt->order, &grid_root_a,
                                        &grid_root_b, &grid_root_c);
        const float q = sourceMagnitude<T, Tcalc, Tcalc2>(pm_wrt->theme, cgr.theme, atom_m.w,
                                                          coord_in_real, sysid, synbk);
        for (int i = 0; i < pm_wrt->order; i++) {
          a_cof[i] *= q;
        }
        spreadDensity<float, float>(a_cof.data(), b_cof.data(), c_cof.data(), pm_wrt->order,
                                    grid_root_a, grid_root_b, grid_root_c, pm_wrt->dims[sysid],
                                    pm_wrt->fftm, pm_wrt->fdata);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename Tcalc2, typename T4>
void accumulateCellDensity(PMIGridAccumulator *pm_acc, const int sysid, const int cell_i,
                           const int cell_j, const int cell_k,
                           const CellGridReader<T, Tacc, Tcalc, T4> &cgr,
                           const SyNonbondedKit<Tcalc, Tcalc2> &synbk) {

  // Re-derive the cell and cell grid boundaries rather than copying them via input arguments
  const ullint cell_bounds = cgr.system_cell_grids[sysid];
  const int cell_offset = (cell_bounds & 0xfffffffLLU);
  const int cell_na = ((cell_bounds >> 28) & 0xfffLLU);
  const int cell_nb = ((cell_bounds >> 40) & 0xfffLLU);
  const int cell_nc = ((cell_bounds >> 52) & 0xfffLLU);

  // Determine limits and other critical constants
  const bool coord_in_real = (cgr.lpos_scale < 1.01);
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const size_t ijk_cellidx = cell_offset + (((cell_k * cell_nb) + cell_j) * cell_na) + cell_i;
  const uint2 ijk_bounds = cgr.cell_limits[ijk_cellidx];
  const uint mllim = ijk_bounds.x;
  const uint mhlim = mllim + (ijk_bounds.y >> 16);
  const int xfrm_stride = roundUp(9, warp_size_int);
  Tcalc umat[9];
  for (int i = 0; i < 9; i++) {
    umat[i] = cgr.system_cell_umat[(sysid * xfrm_stride) + i];
  }
  const Tcalc ng = cgr.mesh_ticks;

  // Lay out arrays to collect B-spline coefficients
  const uint4 grid_dims = pm_acc->dims[sysid];
  switch (pm_acc->mode) {
  case PrecisionModel::DOUBLE:
    {
      std::vector<double> a_cof(pm_acc->order), b_cof(pm_acc->order), c_cof(pm_acc->order);
      for (uint m = mllim; m < mhlim; m++) {
        const T4 atom_m = cgr.image[m];
        int grid_root_a, grid_root_b, grid_root_c;
        particleAlignment<Tcalc, double>(atom_m.x, atom_m.y, atom_m.z, cgr.lpos_inv_scale, umat,
                                         cgr.mesh_ticks, cell_i, cell_j, cell_k, a_cof.data(),
                                         b_cof.data(), c_cof.data(), pm_acc->order, &grid_root_a,
                                         &grid_root_b, &grid_root_c);
        const double q = sourceMagnitude<T, Tcalc, Tcalc2>(pm_acc->theme, cgr.theme, atom_m.w,
                                                           coord_in_real, sysid, synbk);
        for (int i = 0; i < pm_acc->order; i++) {
          a_cof[i] *= q;
        }
        spreadDensity<double, llint>(a_cof.data(), b_cof.data(), c_cof.data(), pm_acc->order,
                                     grid_root_a, grid_root_b, grid_root_c, pm_acc->dims[sysid],
                                     pm_acc->fftm, pm_acc->lldata, pm_acc->overflow,
                                     pm_acc->fp_scale);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      std::vector<float> a_cof(pm_acc->order), b_cof(pm_acc->order), c_cof(pm_acc->order);
      for (uint m = mllim; m < mhlim; m++) {
        const T4 atom_m = cgr.image[m];
        int grid_root_a, grid_root_b, grid_root_c;
        particleAlignment<Tcalc, float>(atom_m.x, atom_m.y, atom_m.z, cgr.lpos_inv_scale, umat,
                                        cgr.mesh_ticks, cell_i, cell_j, cell_k, a_cof.data(),
                                        b_cof.data(), c_cof.data(), pm_acc->order, &grid_root_a,
                                        &grid_root_b, &grid_root_c);
        const float q = sourceMagnitude<T, Tcalc, Tcalc2>(pm_acc->theme, cgr.theme, atom_m.w,
                                                          coord_in_real, sysid, synbk);
        for (int i = 0; i < pm_acc->order; i++) {
          a_cof[i] *= q;
        }
        spreadDensity<float, int>(a_cof.data(), b_cof.data(), c_cof.data(), pm_acc->order,
                                  grid_root_a, grid_root_b, grid_root_c, pm_acc->dims[sysid],
                                  pm_acc->fftm, pm_acc->idata, pm_acc->overflow, pm_acc->fp_scale);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void mapDensity(PMIGrid *pm, const CellGrid<T, Tacc, Tcalc, T4> *cg,
                const AtomGraphSynthesis *poly_ag) {

  // This form of the function will be restricted to work on the CPU.
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcrd_is_real = isFloatingPointScalarType<T>();

  // Loop over all cells.  Taking the B-spline coefficients of each particle along the mesh a, b,
  // and c axes to be numbered 0, 1, ..., n-1 for nth order interpolation, the mesh point at the
  // cell grid spatial decomposition cell's origin receives the a(0) x b(0) x c(0) contribution
  // from the a-, b-, and c-axis B-spline coefficients while subsequent points i, j, and k along
  // the cell's axes receive the a(i) x b(j) x c(k) contributions for i, j, and k in the range
  // [0, n).  Note that the arrays of B-spline coefficients are output in REVERSE order because of
  // the convenience this affords in calculating them.
  const SyNonbondedKit<double, double2> dsynbk = poly_ag->getDoublePrecisionNonbondedKit();
  const SyNonbondedKit<float, float2>  fsynbk = poly_ag->getSinglePrecisionNonbondedKit();
  const CellGridReader<void, void, void, void> cgr_v = cg->templateFreeData();

  // The fact that the topology synthesis (and individual topologies underneath it) contain
  // explicitly-typed parameter arrays in float and double is the fundamental hurdle for smooth
  // templating in this function.  The problem is not severe: make two versions of the cell grid
  // abstract, one of which will be incorrectly typed but unused.
  const CellGridReader<T, Tacc, double, T4> dcgr = restoreType<T, Tacc, double, T4>(cgr_v);
  const CellGridReader<T, Tacc, float, T4>  fcgr = restoreType<T, Tacc, float, T4>(cgr_v);
  
  const bool acc_density_in_real = (pm->fixedPrecisionEnabled() == false);

  // Initialize the grids.
  pm->initialize();
  
  // Outside of the real-valued parameter arrays, dsynbk and fsynbk will hold equivalent sizing
  // constants and atom indexing.  Use either to manage loops and bounds until the precision
  // model demands a bifurcation of further work.
  for (int sysid = 0; sysid < dsynbk.nsys; sysid++) {
    const ullint cell_bounds = (tcalc_is_double) ? dcgr.system_cell_grids[sysid] :
                                                   fcgr.system_cell_grids[sysid];
    const int cell_offset = (cell_bounds & 0xfffffffLLU);
    const int cell_na = ((cell_bounds >> 28) & 0xfffLLU);
    const int cell_nb = ((cell_bounds >> 40) & 0xfffLLU);
    const int cell_nc = ((cell_bounds >> 52) & 0xfffLLU);
    if (acc_density_in_real) {
      PMIGridWriter pm_wrt = pm->data();
      for (int i = 0; i < cell_na; i++) {
        for (int j = 0; j < cell_nb; j++) {
          for (int k = 0; k < cell_nc; k++) {
          
            // Handle the loops over double- or single-precision calculations, integer or real
            // coordinates, and integer or real accumulation of the resulting density.  Perform
            // this switch outside the innermost loop over all atoms to mitigate some of the cost.
            // The nature of the coordinate representation conveyed by the CellGrid abstract, the
            // calculation precision by the non-bonded parameter abstract, and the method of
            // accumulation by the particle-mesh interaction grid abstract (not the CellGrid
            // abstract) each play a part in selecting the proper overload of the density
            // accumulation function.
            if (tcalc_is_double) {
              accumulateCellDensity<T, Tacc, double, double2, T4>(&pm_wrt, sysid, i, j, k, dcgr,
                                                                  dsynbk);
            }
            else {
              accumulateCellDensity<T, Tacc, float, float2, T4>(&pm_wrt, sysid, i, j, k, fcgr,
                                                                fsynbk);
            }
          }
        }
      }
    }
    else {
      PMIGridAccumulator pm_acc = pm->fpData();
      for (int i = 0; i < cell_na; i++) {
        for (int j = 0; j < cell_nb; j++) {
          for (int k = 0; k < cell_nc; k++) {
            if (tcalc_is_double) {
              accumulateCellDensity<T, Tacc, double, double2, T4>(&pm_acc, sysid, i, j, k, dcgr,
                                                                  dsynbk);
            }
            else {
              accumulateCellDensity<T, Tacc, float, float2, T4>(&pm_acc, sysid, i, j, k, fcgr,
                                                                fsynbk);
            }
          }
        }
      }
    }
  }

  // Convert the grid values back to real, if necessary.
  pm->convertToReal();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void mapDensity(PMIGrid *pm, const CellGrid<T, Tacc, Tcalc, T4> &cg,
                const AtomGraphSynthesis &poly_ag) {
  mapDensity(pm, cg.getSelfPointer(), poly_ag.getSelfPointer());
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void mapDensity(PMIGrid *pm, MolecularMechanicsControls *mm_ctrl,
                const CellGrid<T, Tacc, Tcalc, T4> *cg, const AtomGraphSynthesis *poly_ag,
                const CoreKlManager &launcher, const QMapMethod approach) {
  
  // Handle the automated kernel selection here to avoid work needed for specific cases below.
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    if (mm_ctrl == nullptr) {
      rtErr("A kernel based on asynchronous work units must be launched with a valid "
            "MolecularMechanicsControl object.", "mapDensity");
    }
    break;
  case QMapMethod::AUTOMATIC:
    {
      // Recursively call the function based on availability of a MolecularMechanicsControls object
      // and settings in the two grid representations.
      const QMapMethod revised_approach = (mm_ctrl == nullptr) ? QMapMethod::GENERAL_PURPOSE :
                                                                 pm->getRecommendedMappingMethod();
      mapDensity(pm, mm_ctrl, cg, poly_ag, launcher, revised_approach);
    }
    break;
  case QMapMethod::GENERAL_PURPOSE:
    break;
  }
  
  // This form of the function will launch a kernel to operate on the GPU device.  First, find
  // the appropriate kernel to launch.
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const size_t cg_tmat  = std::type_index(typeid(T)).hash_code();
  const CellGridReader<void, void, void, void> v_cgr = cg->templateFreeData(devc_tier);
  const int block_count = launcher.getGpu().getSMPCount();
  PMIGridWriter pm_wrt = pm->data(devc_tier);
  const PrecisionModel calc_prec = (tcalc_is_double) ? PrecisionModel::DOUBLE :
                                                       PrecisionModel::SINGLE;
  const int2 lp = launcher.getDensityMappingKernelDims(approach, calc_prec, pm_wrt.mode,
                                                       pm->useOverflowAccumulation(), cg_tmat,
                                                       pm_wrt.order);
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    if (tcalc_is_double) {
      const SyNonbondedKit<double,
                           double2> synbk = poly_ag->getDoublePrecisionNonbondedKit(devc_tier);
      MMControlKit<double> ctrl = mm_ctrl->dpData(devc_tier);
      mapDensity(&pm_wrt, nullptr, &ctrl, v_cgr, cg_tmat, synbk, block_count, lp, approach, pm);
    }
    else {
      const SyNonbondedKit<float,
                           float2> synbk = poly_ag->getSinglePrecisionNonbondedKit(devc_tier);
      MMControlKit<float> ctrl = mm_ctrl->spData(devc_tier);
      mapDensity(&pm_wrt, nullptr, &ctrl, v_cgr, cg_tmat, synbk, block_count, lp, approach, pm);
    }
    break;
  case QMapMethod::GENERAL_PURPOSE:
    {
      PMIGridAccumulator pm_acc = pm->fpData(devc_tier);
      if (tcalc_is_double) {
        const SyNonbondedKit<double,
                             double2> synbk = poly_ag->getDoublePrecisionNonbondedKit(devc_tier);
        mapDensity(&pm_wrt, &pm_acc, nullptr, v_cgr, cg_tmat, synbk, block_count, lp, approach,
                   pm);
      }
      else {
        const SyNonbondedKit<float,
                             float2> synbk = poly_ag->getSinglePrecisionNonbondedKit(devc_tier);
        mapDensity(&pm_wrt, &pm_acc, nullptr, v_cgr, cg_tmat, synbk, block_count, lp, approach,
                   pm);
      }
    }
    break;
  case QMapMethod::AUTOMATIC:

    // Automated kernel section was carried out above.
    break;
  }
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
void unrollMapDensityCall(PMIGrid *pm, const size_t cg_tacc, const size_t cg_tcalc,
                          const AtomGraphSynthesis *poly_ag) {
  if (cg_tacc == int_type_index) {
    unrollMapDensityCall<T, int, T4>(pm, cg_tcalc, poly_ag);
  }
  else if (cg_tacc == llint_type_index) {
    unrollMapDensityCall<T, llint, T4>(pm, cg_tcalc, poly_ag);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename T4>
void unrollMapDensityCall(PMIGrid *pm, const size_t cg_tcalc, const AtomGraphSynthesis *poly_ag) {
  if (cg_tcalc == double_type_index) {
    const CellGrid<T, Tacc, double, T4>* cgp = pm->getCellGridPointer<T, Tacc, double, T4>();
    mapDensity<T, Tacc, double, T4>(pm, cgp, poly_ag);
  }
  else if (cg_tcalc == float_type_index) {
    const CellGrid<T, Tacc, float, T4>* cgp = pm->getCellGridPointer<T, Tacc, float, T4>();
    mapDensity<T, Tacc, float, T4>(pm, cgp, poly_ag);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
std::vector<Tcalc> mapDensity(const CoordinateFrameReader &cfr, const NonbondedKit<Tcalc> &nbk,
                              const NonbondedTheme theme, const FFTMode fft_staging,
                              const int grid_dim_a, const int grid_dim_b, const int grid_dim_c,
                              const int order, const BSplineUnity unification) {

  // Trap bad inputs
  if (cfr.natom != nbk.natom) {
    rtErr("The number of atoms in the coordinate system (" + std::to_string(cfr.natom) + ") must "
          "match the number of atoms in the topology (" + std::to_string(nbk.natom) + ").",
          "mapDensity");
  }
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
  case NonbondedTheme::VAN_DER_WAALS:
    break;
  case NonbondedTheme::ALL:
    rtErr("Particles must interact via grid-mediated forces in a sanctioned non-bonded "
          "potential.  " + getEnumerationName(theme) + " is invalid.", "mapDensity");
  }
  
  // Transform the atom into the unit cell space to get its position on the density grid.  The
  // transformation will be carried out in double-precision, as would most molecular dynamics
  // re-imaging steps, but coordinates will then be represented in the calculation type for the
  // purposes of finding the grid indices.
  std::vector<double> frac_x(cfr.natom), frac_y(cfr.natom), frac_z(cfr.natom);
  for (int i = 0; i < cfr.natom; i++) {
    frac_x[i] = cfr.xcrd[i];
    frac_y[i] = cfr.ycrd[i];
    frac_z[i] = cfr.zcrd[i];
  }
  imageCoordinates<double, double>(&frac_x, &frac_y, &frac_z, cfr.umat, cfr.invu, cfr.unit_cell,
                                   ImagingMethod::PRIMARY_UNIT_CELL);
  for (int i = 0; i < cfr.natom; i++) {

    // Computing the fractional coordinates in-place is safe because each step changes one of the
    // Cartesian components which is not needed in subsequent steps.
    frac_x[i] = (cfr.umat[0] * frac_x[i]) + (cfr.umat[3] * frac_y[i]) + (cfr.umat[6] * frac_z[i]);
    frac_y[i] =                             (cfr.umat[4] * frac_y[i]) + (cfr.umat[7] * frac_z[i]);
    frac_z[i] =                                                         (cfr.umat[8] * frac_z[i]);
  }

  // Convert to the calculation type before multiplying through by the number of grid points.
  std::vector<Tcalc> tcx(frac_x.begin(), frac_x.end());
  std::vector<Tcalc> tcy(frac_y.begin(), frac_y.end());
  std::vector<Tcalc> tcz(frac_z.begin(), frac_z.end());
  const Tcalc t_nga = grid_dim_a;
  const Tcalc t_ngb = grid_dim_b;
  const Tcalc t_ngc = grid_dim_c;
  for (int i = 0; i < cfr.natom; i++) {
    tcx[i] *= t_nga;
    tcy[i] *= t_ngb;
    tcz[i] *= t_ngc;
  }

  // Initialize the resulting grid.
  int padded_dim_a;
  switch (fft_staging) {
  case FFTMode::IN_PLACE:
    padded_dim_a = 2 * ((grid_dim_a / 2) + 1);
    break;
  case FFTMode::OUT_OF_PLACE:
    padded_dim_a = grid_dim_a;
    break;
  }
  std::vector<Tcalc> result(padded_dim_a * grid_dim_b * grid_dim_c, 0.0);
  
  // Loop over all particles, compute B-spline coefficients, and map the density to the grid.
  std::vector<Tcalc> bspl_a(order), bspl_b(order), bspl_c(order);
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  for (int pos = 0; pos < cfr.natom; pos++) {
    const int base_ga = tcx[pos];
    const int base_gb = tcy[pos];
    const int base_gc = tcz[pos];
    const Tcalc da = tcx[pos] - static_cast<Tcalc>(base_ga);
    const Tcalc db = tcy[pos] - static_cast<Tcalc>(base_gb);
    const Tcalc dc = tcz[pos] - static_cast<Tcalc>(base_gc);
    switch (unification) {
    case BSplineUnity::CENTER_FILL:
      bSpline(da, order, bspl_a.data());
      bSpline(db, order, bspl_b.data());
      bSpline(dc, order, bspl_c.data());
      break;
    case BSplineUnity::NONE:
      bspl_a = bSplineNoUnity<Tcalc>(da, order);
      bspl_b = bSplineNoUnity<Tcalc>(db, order);
      bspl_c = bSplineNoUnity<Tcalc>(dc, order);
      break;
    }
    Tcalc q;
    switch(theme) {
    case NonbondedTheme::ELECTROSTATIC:
      q = nbk.charge[pos];
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      {
        const Tcalc ljb = nbk.ljb_coeff[(nbk.n_lj_types + 1) * nbk.lj_idx[pos]];
        q = (tcalc_is_double) ? 0.5 * sqrt(ljb) : 0.5f * sqrtf(ljb);
      }
      break;
    case NonbondedTheme::ALL:
      break;
    }
    for (int i = 0; i < order; i++) {
      bspl_a[i] *= q;
    }
    for (int k = 0; k < order; k++) {
      int act_k = base_gc + k;
      act_k += ((act_k < 0) - (act_k >= grid_dim_c)) * grid_dim_c;
      for (int j = 0; j < order; j++) {
        int act_j = base_gb + j;
        act_j += ((act_j < 0) - (act_j >= grid_dim_b)) * grid_dim_b;
        const Tcalc bspl_bc = bspl_b[j] * bspl_c[k];
        const int jk_idx = ((act_k * grid_dim_b) + act_j) * padded_dim_a;
        for (int i = 0; i < order; i++) {
          int act_i = base_ga + i;
          act_i += ((act_i < 0) - (act_i >= grid_dim_a)) * grid_dim_a;
          result[jk_idx + act_i] += bspl_bc * bspl_a[i];
        }
      }
    }
  }
  return result;
}

} // namespace energy
} // namespace stormm
