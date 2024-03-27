// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4>::CellGridWriter(const NonbondedTheme theme_in,
                                                   const int system_count_in,
                                                   const int total_cell_count_in,
                                                   const int total_chain_count_in,
                                                   const int mesh_ticks_in,
                                                   const size_t cell_base_capacity_in,
                                                   const size_t cell_excl_capacity_in,
                                                   const float lpos_scale_in,
                                                   const float lpos_inv_scale_in,
                                                   ullint* system_cell_grids_in,
                                                   Tcalc* system_cell_umat_in,
                                                   T* system_cell_invu_in, T* system_pmig_invu_in,
                                                   uint2* cell_limits_in,
                                                   const uint2* cell_limits_old_in,
                                                   const uint* chain_limits_in,
                                                   const int* system_chain_bounds_in,
                                                   const int* chain_system_owner_in, T4* image_in,
                                                   const T4* image_old_in, uint2* entry_room_in,
                                                   uint2* exit_room_in,
                                                   const int* transit_room_bounds_in,
                                                   int* entry_room_counts_in,
                                                   int* exit_room_counts_in,
                                                   int* nonimg_atom_idx_in, uint* img_atom_idx_in,
                                                   uint* exclusion_maps_in,
                                                   const uint* exclusion_maps_old_in,
                                                   const int* nt_groups_in, Tacc* xfrc_in,
                                                   Tacc* yfrc_in, Tacc* zfrc_in, int* xfrc_ovrf_in,
                                                   int* yfrc_ovrf_in, int* zfrc_ovrf_in,
                                                   Tacc* xfrc_hw_in, Tacc* yfrc_hw_in,
                                                   Tacc* zfrc_hw_in, int* xfrc_hw_ovrf_in,
                                                   int* yfrc_hw_ovrf_in, int* zfrc_hw_ovrf_in) :
    theme{theme_in}, system_count{system_count_in}, total_cell_count{total_cell_count_in},
    total_chain_count{total_chain_count_in}, mesh_ticks{mesh_ticks_in},
    cell_base_capacity{cell_base_capacity_in}, cell_excl_capacity{cell_excl_capacity_in},
    lpos_scale{lpos_scale_in}, lpos_inv_scale{lpos_inv_scale_in},
    system_cell_grids{system_cell_grids_in}, system_cell_umat{system_cell_umat_in},
    system_cell_invu{system_cell_invu_in}, system_pmig_invu{system_pmig_invu_in},
    cell_limits{cell_limits_in}, chain_limits{chain_limits_in},
    system_chain_bounds{system_chain_bounds_in}, chain_system_owner{chain_system_owner_in},
    image{image_in}, image_old{image_old_in}, entry_room{entry_room_in}, exit_room{exit_room_in},
    transit_room_bounds{transit_room_bounds_in}, entry_room_counts{entry_room_counts_in},
    exit_room_counts{exit_room_counts_in}, nonimg_atom_idx{nonimg_atom_idx_in},
    img_atom_idx{img_atom_idx_in}, exclusion_maps{exclusion_maps_in},
    exclusion_maps_old{exclusion_maps_old_in}, nt_groups{nt_groups_in}, xfrc{xfrc_in},
    yfrc{yfrc_in}, zfrc{zfrc_in}, xfrc_ovrf{xfrc_ovrf_in}, yfrc_ovrf{yfrc_ovrf_in},
    zfrc_ovrf{zfrc_ovrf_in}, xfrc_hw{xfrc_hw_in}, yfrc_hw{yfrc_hw_in}, zfrc_hw{zfrc_hw_in},
    xfrc_hw_ovrf{xfrc_hw_ovrf_in}, yfrc_hw_ovrf{yfrc_hw_ovrf_in}, zfrc_hw_ovrf{zfrc_hw_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridReader<T, Tacc, Tcalc, T4>::CellGridReader(const NonbondedTheme theme_in,
                                                   const int system_count_in,
                                                   const int total_cell_count_in,
                                                   const int total_chain_count_in,
                                                   const int mesh_ticks_in,
                                                   const size_t cell_base_capacity_in,
                                                   const size_t cell_excl_capacity_in,
                                                   const float lpos_scale_in,
                                                   const float lpos_inv_scale_in,
                                                   const ullint* system_cell_grids_in,
                                                   const Tcalc* system_cell_umat_in,
                                                   const T* system_cell_invu_in,
                                                   const T* system_pmig_invu_in,
                                                   const uint2* cell_limits_in,
                                                   const uint* chain_limits_in,
                                                   const int* system_chain_bounds_in,
                                                   const int* chain_system_owner_in,
                                                   const T4* image_in,
                                                   const int* nonimg_atom_idx_in,
                                                   const uint* img_atom_idx_in,
                                                   const uint* exclusion_maps_in,
                                                   const int* nt_groups_in,
                                                   const Tacc* xfrc_in, const Tacc* yfrc_in,
                                                   const Tacc* zfrc_in, const int* xfrc_ovrf_in,
                                                   const int* yfrc_ovrf_in,
                                                   const int* zfrc_ovrf_in) :
    theme{theme_in}, system_count{system_count_in}, total_cell_count{total_cell_count_in},
    total_chain_count{total_chain_count_in}, mesh_ticks{mesh_ticks_in},
    cell_base_capacity{cell_base_capacity_in}, cell_excl_capacity{cell_excl_capacity_in},
    lpos_scale{lpos_scale_in}, lpos_inv_scale{lpos_inv_scale_in},
    system_cell_grids{system_cell_grids_in}, system_cell_umat{system_cell_umat_in},
    system_cell_invu{system_cell_invu_in}, system_pmig_invu{system_pmig_invu_in},
    cell_limits{cell_limits_in}, chain_limits{chain_limits_in},
    system_chain_bounds{system_chain_bounds_in}, chain_system_owner{chain_system_owner_in},
    image{image_in}, nonimg_atom_idx{nonimg_atom_idx_in}, img_atom_idx{img_atom_idx_in},
    exclusion_maps{exclusion_maps_in}, nt_groups{nt_groups_in}, xfrc{xfrc_in}, yfrc{yfrc_in},
    zfrc{zfrc_in}, xfrc_ovrf{xfrc_ovrf_in}, yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridReader<T, Tacc, Tcalc, T4>::CellGridReader(const CellGridWriter<T, Tacc, Tcalc, T4> &cgw) :
    theme{cgw.theme}, system_count{cgw.system_count}, total_cell_count{cgw.total_cell_count},
    total_chain_count{cgw.total_chain_count}, mesh_ticks{cgw.mesh_ticks},
    cell_base_capacity{cgw.cell_base_capacity}, cell_excl_capacity{cgw.cell_excl_capacity},
    lpos_scale{cgw.lpos_scale}, lpos_inv_scale{cgw.lpos_inv_scale},
    system_cell_grids{cgw.system_cell_grids}, system_cell_umat{cgw.system_cell_umat},
    system_cell_invu{cgw.system_cell_invu}, system_pmig_invu{cgw.system_pmig_invu},
    cell_limits{cgw.cell_limits}, chain_limits{cgw.chain_limits},
    system_chain_bounds{cgw.system_chain_bounds}, chain_system_owner{cgw.chain_system_owner},
    image{cgw.image}, nonimg_atom_idx{cgw.nonimg_atom_idx}, img_atom_idx{cgw.img_atom_idx},
    exclusion_maps{cgw.exclusion_maps}, nt_groups{cgw.nt_groups}, xfrc{cgw.xfrc}, yfrc{cgw.yfrc},
    zfrc{cgw.zfrc}, xfrc_ovrf{cgw.xfrc_ovrf}, yfrc_ovrf{cgw.yfrc_ovrf}, zfrc_ovrf{cgw.zfrc_ovrf}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridReader<T, Tacc, Tcalc, T4>::CellGridReader(const CellGridWriter<T, Tacc, Tcalc, T4> *cgw) :
    theme{cgw->theme}, system_count{cgw->system_count}, total_cell_count{cgw->total_cell_count},
    total_chain_count{cgw->total_chain_count}, mesh_ticks{cgw->mesh_ticks},
    cell_base_capacity{cgw->cell_base_capacity}, cell_excl_capacity{cgw->cell_excl_capacity},
    lpos_scale{cgw->lpos_scale}, lpos_inv_scale{cgw->lpos_inv_scale},
    system_cell_grids{cgw->system_cell_grids}, system_cell_umat{cgw->system_cell_umat},
    system_cell_invu{cgw->system_cell_invu}, system_pmig_invu{cgw->system_pmig_invu},
    cell_limits{cgw->cell_limits}, chain_limits{cgw->chain_limits},
    system_chain_bounds{cgw->system_chain_bounds}, chain_system_owner{cgw->chain_system_owner},
    image{cgw->image}, nonimg_atom_idx{cgw->nonimg_atom_idx}, img_atom_idx{cgw->img_atom_idx},
    exclusion_maps{cgw->exclusion_maps}, nt_groups{cgw->nt_groups}, xfrc{cgw->xfrc},
    yfrc{cgw->yfrc}, zfrc{cgw->zfrc}, xfrc_ovrf{cgw->xfrc_ovrf}, yfrc_ovrf{cgw->yfrc_ovrf},
    zfrc_ovrf{cgw->zfrc_ovrf}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(const PhaseSpaceSynthesis *poly_ps_ptr_in,
                                       const AtomGraphSynthesis *poly_ag_ptr_in,
                                       const double cutoff, const double padding,
                                       const int mesh_subdivisions_in,
                                       const NonbondedTheme theme_in, const GpuDetails &gpu,
                                       const size_t cell_base_capacity_in,
                                       const ExceptionResponse policy_in) :
    system_count{0}, total_cell_count{0}, total_chain_count{0},
    cell_base_capacity{roundUp<size_t>(cell_base_capacity_in, 32)},
    cell_excl_capacity{(cell_base_capacity * cellgrid_total_excl_allotment *
                        ((cell_base_capacity / 32)) + 1)},
    effective_cutoff{validateEffectiveCutoff(cutoff + padding, policy_in)},
    mesh_subdivisions{mesh_subdivisions_in},
    policy{policy_in},
    theme{theme_in},
    cycle_position{CoordinateCycle::WHITE},
    localpos_scale_bits{0}, localpos_scale{1.0}, localpos_inverse_scale{1.0},
    system_cell_grids{HybridKind::ARRAY, "cg_cell_grids"},
    system_cell_umat{HybridKind::ARRAY, "cg_dp_cell_umat"},
    system_cell_umat_alt{HybridKind::ARRAY, "cg_dp_cell_umat_alt"},
    system_cell_invu{HybridKind::ARRAY, "cg_dp_cell_invu"},
    system_cell_invu_alt{HybridKind::ARRAY, "cg_dp_cell_invu_alt"},
    system_pmig_invu{HybridKind::ARRAY, "cg_dp_pmig_invu"},
    system_pmig_invu_alt{HybridKind::ARRAY, "cg_dp_pmig_invu_alt"},
    image_cell_limits{HybridKind::ARRAY, "cg_cell_bounds"},
    image_cell_limits_alt{HybridKind::ARRAY, "cg_cell_bounds_alt"},
    image_chain_limits{HybridKind::ARRAY, "cg_chain_bounds"},
    system_chain_bounds{HybridKind::ARRAY, "cg_syschn_bounds"},
    chain_system_membership{HybridKind::ARRAY, "cg_chnsys_owner"},
    image{HybridKind::ARRAY, "cg_image"},
    image_alt{HybridKind::ARRAY, "cg_image_alt"},
    image_chain_entering{HybridKind::ARRAY, "cg_entering_atoms"},
    image_chain_leaving{HybridKind::ARRAY, "cg_leaving_atoms"},
    image_vestibule_limits{HybridKind::ARRAY, "cg_vestibule_limits"},
    image_chain_entering_counts{HybridKind::ARRAY, "cg_entering_counts"},
    image_chain_leaving_counts{HybridKind::ARRAY, "cg_leaving_counts"},
    nonimaged_atom_indices{HybridKind::ARRAY, "cg_atom_idx"},
    nonimaged_atom_indices_alt{HybridKind::ARRAY, "cg_atom_idx_alt"},
    image_array_indices{HybridKind::ARRAY, "cg_img_idx"},
    image_array_indices_alt{HybridKind::ARRAY, "cg_img_idx_alt"},
    exclusion_maps{HybridKind::ARRAY, "cg_exclusions"},
    exclusion_maps_alt{HybridKind::ARRAY, "cg_exclusions_alt"},
    nt_work_groups{HybridKind::ARRAY, "cg_nt_work_groups"},
    x_force{HybridKind::ARRAY, "cg_xfrc"},
    y_force{HybridKind::ARRAY, "cg_yfrc"},
    z_force{HybridKind::ARRAY, "cg_zfrc"},
    x_force_overflow{HybridKind::ARRAY, "cg_xfrc_ovrf"},
    y_force_overflow{HybridKind::ARRAY, "cg_yfrc_ovrf"},
    z_force_overflow{HybridKind::ARRAY, "cg_zfrc_ovrf"},
    warp_x_work{HybridKind::ARRAY, "cg_x_warp_work"},
    warp_y_work{HybridKind::ARRAY, "cg_y_warp_work"},
    warp_z_work{HybridKind::ARRAY, "cg_z_warp_work"},
    warp_x_overflow_work{HybridKind::ARRAY, "cg_x_ovrf_work"},
    warp_y_overflow_work{HybridKind::ARRAY, "cg_y_ovrf_work"},
    warp_z_overflow_work{HybridKind::ARRAY, "cg_z_ovrf_work"},
    poly_ps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps_ptr_in)},
    poly_ag_ptr{const_cast<AtomGraphSynthesis*>(poly_ag_ptr_in)}
{
  // Exit with an error if a valid coordinate synthesis was not supplied.
  if (poly_ps_ptr_in == nullptr) {
    rtErr("A valid coordinate synthesis must be provided.", "CellGrid");
  }
  validateCoordinateSynthesis();
  const PsSynthesisReader poly_psr = poly_ps_ptr_in->data();
  system_count = poly_psr.system_count;
  cycle_position = poly_ps_ptr->getCyclePosition();

  // Enforce int or llint data type for accumulation
  const size_t tc_acc = std::type_index(typeid(Tacc)).hash_code();
  if (tc_acc != int_type_index && tc_acc != llint_type_index) {
    rtErr("The only acceptable types for storing neighbor list images are 32- and 64-bit signed "
          "integers.", "CellGrid");
  }

  // Enforce bit scaling depending on the coordinate representation, as well as matches between
  // the cell transformation and the coordinate representation.
  const size_t tc_crd = std::type_index(typeid(T4)).hash_code();
  const size_t tc_mat = std::type_index(typeid(T)).hash_code();
  if (tc_crd == float4_type_index || tc_crd == double4_type_index) {
    if (tc_mat != float_type_index && tc_mat != double_type_index) {
      rtErr("Floating-point coordinate representations (" + getStormmHpcVectorTypeName<T4>() +
            ") must be combined with floating-point cell axis representations (" +
            getStormmScalarTypeName<T>() + ").", "CellGrid");
    }
    if (localpos_scale_bits > 0) {
      localpos_scale_bits = 0;
    }
  }
  else if (tc_crd == int4_type_index) {
    if (tc_mat != int_type_index) {
      rtErr("Coordinates expressed a 32-bit integer format must use a 32-bit integer "
            "representation of the cell axes.", "CellGrid");
    }
  }
  else if (tc_crd == longlong4_type_index) {
    if (tc_mat != llint_type_index) {
      rtErr("Coordinates expressed a 64-bit integer format must use a 64-bit integer "
            "representation of the cell axes.", "CellGrid");
    }
  }
  else {
    rtErr("Local coordinates must be represented in a signed format of minimum size 32 bits.  "
          "Type " + getStormmHpcVectorTypeName<T4>() + " is invalid.", "CellGrid");
  }

  // Enforce float or double for calculations.  The GPU does not perform calculations in long
  // double, which does qualify as a floating point type, but such cases will be trapped by
  // compile-time restrictions.
  const size_t tc_calc = std::type_index(typeid(Tcalc)).hash_code();
  if (isFloatingPointScalarType<Tcalc>() == false) {
    rtErr("Calculations, and the encoding of transformation matrices to take Cartesian "
          "coordinates into the local frames of reference for each spatial decomposition cell, "
          "must take a floating-point scalar data type.", "CellGrid");
  }

  // Enforce bit scaling in terms of force accumulation.
  const int force_scale_bits = poly_ps_ptr->getForceAccumulationBits();
  if (tc_acc == int_type_index && force_scale_bits > force_scale_nonoverflow_bits) {
    rtErr("Force accumulation cannot take place with " + std::to_string(force_scale_bits) +
          " bits of unit precision if the primary accumulators are 32-bit integers.", "CellGrid");
  }
  
  // Determine the total number of cells in each system, and overall
  const int xfrm_stride = roundUp(9, warp_size_int);
  std::vector<int> cell_na(system_count), cell_nb(system_count), cell_nc(system_count);
  llint tmp_total_cell_count = 0LL;
  int tmp_total_chain_count = 0;
  const std::vector<uint> prime_factors = { 2, 3, 5, 7 };
  const ullint big_product = ipowl(2, 14) * ipowl(3, 10) * ipowl(5, 6) * ipowl(7, 4);
  for (int i = 0; i < poly_psr.system_count; i++) {
    double cell_a, cell_b, cell_c;
    hessianNormalWidths(&poly_psr.invu[i * xfrm_stride], &cell_a, &cell_b, &cell_c);
    cell_na[i] = floor(cell_a / effective_cutoff);
    cell_nb[i] = floor(cell_b / effective_cutoff);
    cell_nc[i] = floor(cell_c / effective_cutoff);
    const int3 best_abc = optimizeCellConfiguration(cell_na[i], cell_nb[i], cell_nc[i],
                                                    mesh_subdivisions);
    cell_na[i] = best_abc.x;
    cell_nb[i] = best_abc.y;
    cell_nc[i] = best_abc.z;
    if (cell_na[i] > maximum_spatial_decomposition_cells ||
        cell_nb[i] > maximum_spatial_decomposition_cells ||
        cell_nc[i] > maximum_spatial_decomposition_cells) {
      rtErr("The system spatial decomposition would require a grid of [ " +
            std::to_string(cell_na[i]) + " x " + std::to_string(cell_nb[i]) + " x " +
            std::to_string(cell_nc[i]) + " ] cells.  A maximum of " +
            std::to_string(maximum_spatial_decomposition_cells) + " is enforced along all axes.",
            "CellGrid");
    }
    tmp_total_cell_count += cell_na[i] * cell_nb[i] * cell_nc[i];
    tmp_total_chain_count += cell_nb[i] * cell_nc[i];
  }
  if (tmp_total_cell_count > static_cast<llint>(maximum_cellgrid_cell_count)) {
    rtErr("Too many cells would be required by the systems at hand (" +
          std::to_string(tmp_total_cell_count) + " versus a maximum of " +
          std::to_string(maximum_cellgrid_cell_count) + ").", "CellGrid");
  }
  total_cell_count = tmp_total_cell_count;
  total_chain_count = tmp_total_chain_count;

  // Estimate the proportion of particles that might wander into or out of each cell in any given
  // step.  
  const double migratory_factor = computeMigrationRate(effective_cutoff, 0.1);

  // Allocate space for cells and supporting arrays.  Initialize the partition cell transformation
  // matrices.
  system_cell_grids.resize(system_count);
  system_cell_umat.resize(system_count * xfrm_stride);
  system_cell_umat_alt.resize(system_count * xfrm_stride);
  system_cell_invu.resize(system_count * xfrm_stride);
  system_cell_invu_alt.resize(system_count * xfrm_stride);
  system_pmig_invu.resize(system_count * xfrm_stride);
  system_pmig_invu_alt.resize(system_count * xfrm_stride);
  image_cell_limits.resize(total_cell_count);
  image_cell_limits_alt.resize(total_cell_count);
  image_chain_limits.resize(total_chain_count + 1);
  system_chain_bounds.resize(system_count + 1);
  chain_system_membership.resize(total_chain_count);
  std::vector<double> umat_cell_stage(system_count * xfrm_stride, 0.0);
  std::vector<double> umat_cell_stage_alt(system_count * xfrm_stride, 0.0);
  std::vector<double> invu_cell_stage(system_count * xfrm_stride, 0.0);
  std::vector<double> invu_cell_stage_alt(system_count * xfrm_stride, 0.0);
  std::vector<double> invu_pmig_stage(system_count * xfrm_stride, 0.0);
  std::vector<double> invu_pmig_stage_alt(system_count * xfrm_stride, 0.0);
  double* umat_stage_ptr = umat_cell_stage.data();
  double* umat_stage_ptr_alt = umat_cell_stage_alt.data();
  double* cell_stage_ptr = invu_cell_stage.data();
  double* cell_stage_ptr_alt = invu_cell_stage_alt.data();
  double* pmig_stage_ptr = invu_pmig_stage.data();
  double* pmig_stage_ptr_alt = invu_pmig_stage_alt.data();
  int cell_counter = 0;
  int chain_counter = 0;
  for (int i = 0; i < poly_psr.system_count; i++) {
    const double cdim_a = 1.0 / static_cast<double>(cell_na[i]);
    const double cdim_b = 1.0 / static_cast<double>(cell_nb[i]);
    const double cdim_c = 1.0 / static_cast<double>(cell_nc[i]);
    const std::vector<double> cell_isclr = { static_cast<double>(cell_na[i]), 0.0, 0.0,
                                             0.0, static_cast<double>(cell_nb[i]), 0.0,
                                             0.0, 0.0, static_cast<double>(cell_nc[i]) };
    const std::vector<double> cell_sclr = { cdim_a, 0.0, 0.0, 0.0, cdim_b, 0.0, 0.0, 0.0, cdim_c };
    matrixMultiply<double>(&poly_psr.umat[i * xfrm_stride], 3, 3, cell_isclr.data(), 3, 3,
                           &umat_stage_ptr[i * xfrm_stride]);
    matrixMultiply<double>(&poly_psr.umat_alt[i * xfrm_stride], 3, 3, cell_isclr.data(), 3, 3,
                           &umat_stage_ptr_alt[i * xfrm_stride]);
    matrixMultiply<double>(&poly_psr.invu[i * xfrm_stride], 3, 3, cell_sclr.data(), 3, 3,
                           &cell_stage_ptr[i * xfrm_stride]);
    matrixMultiply<double>(&poly_psr.invu_alt[i * xfrm_stride], 3, 3, cell_sclr.data(), 3, 3,
                           &cell_stage_ptr_alt[i * xfrm_stride]);
    const double gdim_a = 1.0 / static_cast<double>(cell_na[i] * mesh_subdivisions);
    const double gdim_b = 1.0 / static_cast<double>(cell_nb[i] * mesh_subdivisions);
    const double gdim_c = 1.0 / static_cast<double>(cell_nc[i] * mesh_subdivisions);
    const std::vector<double> pmig_sclr = { gdim_a, 0.0, 0.0, 0.0, gdim_b, 0.0, 0.0, 0.0, gdim_c };
    matrixMultiply<double>(&poly_psr.invu[i * xfrm_stride], 3, 3, pmig_sclr.data(), 3, 3,
                           &pmig_stage_ptr[i * xfrm_stride]);
    matrixMultiply<double>(&poly_psr.invu_alt[i * xfrm_stride], 3, 3, pmig_sclr.data(), 3, 3,
                           &pmig_stage_ptr_alt[i * xfrm_stride]);

    // Pack the cell sizes and simulation starting points.
    ullint icells = cell_counter;
    icells |= (static_cast<ullint>(cell_na[i]) << 28);
    icells |= (static_cast<ullint>(cell_nb[i]) << 40);
    icells |= (static_cast<ullint>(cell_nc[i]) << 52);
    system_cell_grids.putHost(icells, i);
    cell_counter += cell_na[i] * cell_nb[i] * cell_nc[i];
    system_chain_bounds.putHost(chain_counter, i);
    for (int j = 0; j < cell_nb[i] * cell_nc[i]; j++) {
      chain_system_membership.putHost(i, chain_counter + j);
    }
    chain_counter += cell_nb[i] * cell_nc[i];
  }
  system_chain_bounds.putHost(chain_counter, poly_psr.system_count);
  computeFixedPrecisionModel(invu_cell_stage);
  std::vector<T> cell_imatrix(9), cell_imatrix_alt(9), pmig_imatrix(9), pmig_imatrix_alt(9);
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int ipos = i * xfrm_stride;
    if (localpos_scale_bits > 0) {
      for (int j = 0; j < 9; j++) {
        cell_imatrix[j]     = llround(invu_cell_stage[ipos + j] * localpos_scale);
        cell_imatrix_alt[j] = llround(invu_cell_stage_alt[ipos + j] * localpos_scale);
        pmig_imatrix[j]     = llround(invu_pmig_stage[ipos + j] * localpos_scale);
        pmig_imatrix_alt[j] = llround(invu_pmig_stage_alt[ipos + j] * localpos_scale);
      }
    }
    else {
      for (int j = 0; j < 9; j++) {
        cell_imatrix[j]     = invu_cell_stage[ipos + j];
        cell_imatrix_alt[j] = invu_cell_stage_alt[ipos + j];
        pmig_imatrix[j]     = invu_pmig_stage[ipos + j];
        pmig_imatrix_alt[j] = invu_pmig_stage_alt[ipos + j];
      }
    }
    system_cell_invu.putHost(cell_imatrix, ipos, 9);
    system_cell_invu_alt.putHost(cell_imatrix_alt, ipos, 9);
    system_pmig_invu.putHost(pmig_imatrix, ipos, 9);
    system_pmig_invu_alt.putHost(pmig_imatrix_alt, ipos, 9);
    for (int j = 0; j < 9; j++) {
      system_cell_umat.putHost(umat_cell_stage[ipos + j], ipos + j);
      system_cell_umat_alt.putHost(umat_cell_stage_alt[ipos + j], ipos + j);
    }
  }

  // Check that the cell base capacity is large enough.  This requires its own pass over all atoms
  // in all systems.  If the cell base capacity is already large enough, then its value will not be
  // changed.  This allows user input to set the cell base capacity to an arbitrarily high value,
  // in the event that a simulation is known to need more than the capacity which might be expected
  // by this algorithm.
  std::vector<int> cell_dest, cell_populations;
  int trial_base_capacity = cell_base_capacity;
  for (int pos = 0; pos < poly_psr.system_count; pos++) {
    tallyCellPopulations(&cell_populations, pos, cell_na[pos], cell_nb[pos], cell_nc[pos]);
    const int max_pop = maxValue(cell_populations);
    cell_base_capacity = std::max(roundUp<size_t>(max_pop * 3 / 2, 32), cell_base_capacity);
  }
  
  // Allocate the bulk of the space that the object will use in partitioning atoms.
  const size_t image_size = static_cast<size_t>(total_cell_count) * cell_base_capacity;
  image.resize(image_size);
  image_alt.resize(image_size);
  const int total_atoms = poly_psr.atom_starts[poly_psr.system_count - 1] +
                          roundUp(poly_psr.atom_counts[poly_psr.system_count - 1], warp_size_int);
  image_array_indices.resize(total_atoms);
  image_array_indices_alt.resize(total_atoms);
  x_force.resize(image_size);
  y_force.resize(image_size);
  z_force.resize(image_size);
  x_force_overflow.resize(image_size);
  y_force_overflow.resize(image_size);
  z_force_overflow.resize(image_size);
  const size_t workspace_size = (gpu.getSMPCount() * gpu.getMaxThreadsPerSMP() / warp_size_int) *
                                cell_base_capacity;
  warp_x_work.resize(workspace_size);
  warp_y_work.resize(workspace_size);
  warp_z_work.resize(workspace_size);
  warp_x_overflow_work.resize(workspace_size);
  warp_y_overflow_work.resize(workspace_size);
  warp_z_overflow_work.resize(workspace_size);

  // Initialize atoms in the cells by looping over all systems on the CPU.  No GPU equivalent is
  // economical as the imaging calculations are expensive relative to the information uploads and
  // downloads that would be required.
  std::vector<double> atom_x, atom_y, atom_z;
  uint icl_atom_counter = 0U;
  for (int pos = 0; pos < poly_psr.system_count; pos++) {

    // Each chain of the system will need its boundaries defined, and they are intended not to
    // change throughout the life of the object.  If it needs to change, the object might be
    // destroyed so that a new cell grid can be built in its place (the atoms and their official
    // positions are stored in the higher-precision PhaseSpaceSynthesis), or the program will
    // raise an exception.
    chain_counter = system_chain_bounds.readHost(pos);
    const uint max_atoms_per_chain = static_cast<size_t>(cell_na[pos]) * cell_base_capacity;
    for (int k = 0; k < cell_nc[pos]; k++) {
      for (int j = 0; j < cell_nb[pos]; j++) {
        image_chain_limits.putHost(icl_atom_counter, chain_counter + (k * cell_nb[pos]) + j);
        icl_atom_counter += max_atoms_per_chain;
      }
    }

    // Atoms cannot be placed into the images until the images have been allocated, and the
    // images cannot be allocated until the total padded size is known.  Storing the atoms in
    // their imaged states, however, would require a great deal of extra memory.  The
    // conservative option is to run through the systems once more and recompute the results
    // once the arrays have been allocated.
  }

  // Set the upper bounds for various arrays
  image_chain_limits.putHost(icl_atom_counter,
                             system_chain_bounds.readHost(poly_psr.system_count));

  // Allocate the images and related bounds arrays
  image.resize(icl_atom_counter);
  image_alt.resize(icl_atom_counter);
  nonimaged_atom_indices.resize(icl_atom_counter);
  nonimaged_atom_indices_alt.resize(icl_atom_counter);
    
  // Loop back over all structures, reimage atoms, and pack the image arrays.
  populateImage(CoordinateCycle::WHITE);
  populateImage(CoordinateCycle::BLACK);
  
  // Create work units to support the "tower and plate" neutral territory decomposition
  prepareWorkGroups();
  
  // Handle the exclusions maps.  This is integral to the manner in which non-bonded interactions
  // are processed.  The tower-plate neutral territory method organizes the interactions, with the
  // "tower" consisting of five cells and the "plate" another twelve.  The "plate" is aligned with
  // cell chains, in the AB plane of the system's particular cell grid, while the "tower" runs
  // along the C axis.  Atoms of the smaller tower will be cached in their properly imaged states,
  // positioned relative to the home cell, which is also the intersection of the tower and plate.
  // Exclusions are stored for atoms of the central three cells of the tower interacting with each
  // other and with the central four cells unique to the plate: a total of fourteen cell:cell pairs
  // with demarcations in the associated exclusion_demaractions array.  Each cell gets an array of
  // 32 byte-aligned uint values indicating the starting point and number of relevant exclusion bit
  // masks for one of the fourteen pairs, as well as the offset for the cell's exclusions overall.
  // In the first tranch of exclusion masks, one bit mask is present for each atom of the four
  // central cells of the plate, containing bits for each atom of the central cells of the tower.
  // This covers twelve of the cell:cell interactions.  In the second tranch of exclusion masks,
  // each atom of the lower two of the central cells of the tower (home cell and the -c cell) gets
  // a bit mask and each mask contains one bit for every atom of the home cell.
  initializeExclusionMasks(CoordinateCycle::WHITE);
  initializeExclusionMasks(CoordinateCycle::BLACK);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(const PhaseSpaceSynthesis &poly_ps_in,
                                       const AtomGraphSynthesis &poly_ag_in, const double cutoff,
                                       const double padding, const int mesh_subdivisions_in,
                                       const NonbondedTheme theme_in, const GpuDetails &gpu,
                                       const size_t base_capacity_in,
                                       const ExceptionResponse policy_in) :
    CellGrid(poly_ps_in.getSelfPointer(), poly_ag_in.getSelfPointer(), cutoff, padding,
             mesh_subdivisions_in, theme_in, gpu, base_capacity_in, policy_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getTotalCellCount() const {
  return total_cell_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getCellCount(const int index) const {
  return (system_cell_grids.readHost(index + 1) & 0xfffffffLLU) -
         (system_cell_grids.readHost(index    ) & 0xfffffffLLU);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getCellCount(const int index, const UnitCellAxis axis) const {
  const ullint st_val = system_cell_grids.readHost(index);
  switch (axis) {
  case UnitCellAxis::A:
    return ((st_val >> 28) & 0xfff);
  case UnitCellAxis::B:
    return ((st_val >> 40) & 0xfff);
  case UnitCellAxis::C:
    return ((st_val >> 52) & 0xfff);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getCellCount(const int index,
                                               const CartesianDimension axis) const {
  switch (axis) {
  case CartesianDimension::X:
    return getCellCount(index, UnitCellAxis::A);
  case CartesianDimension::Y:
    return getCellCount(index, UnitCellAxis::B);
  case CartesianDimension::Z:
    return getCellCount(index, UnitCellAxis::C);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
size_t CellGrid<T, Tacc, Tcalc, T4>::getCellBaseCapacity() const {
  return cell_base_capacity;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
double CellGrid<T, Tacc, Tcalc, T4>::getEffectiveCutoff() const {
  return effective_cutoff;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
NonbondedTheme CellGrid<T, Tacc, Tcalc, T4>::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CoordinateCycle CellGrid<T, Tacc, Tcalc, T4>::getCyclePosition() const {
  return cycle_position;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getPositionScalingBits() const {
  return localpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int CellGrid<T, Tacc, Tcalc, T4>::getMeshSubdivisions() const {
  return mesh_subdivisions;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
uint CellGrid<T, Tacc, Tcalc, T4>::getImageIndex(const int system_index, const int atom_index,
                                                 const CoordinateCycle orientation) const {
  validateAtomIndex(system_index, atom_index, "getImageIndex");
  const int synth_idx = poly_ps_ptr->getAtomOffset(system_index) + atom_index;
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return image_array_indices_alt.readHost(synth_idx);
  case CoordinateCycle::WHITE:
    return image_array_indices.readHost(synth_idx);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
uint CellGrid<T, Tacc, Tcalc, T4>::getImageIndex(const int system_index,
                                                 const int atom_index) const {
  return getImageIndex(system_index, atom_index, cycle_position);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int4 CellGrid<T, Tacc, Tcalc, T4>::getCellLocation(const int system_index, const int atom_index,
                                                   const CoordinateCycle orientation) const {
  validateAtomIndex(system_index, atom_index, "getCellLocation");
  const int synth_idx = poly_ps_ptr->getAtomOffset(system_index) + atom_index;
  const ullint sys_cdims = system_cell_grids.readHost(system_index);
  const uint cell_init = (sys_cdims & 0xfffffff);
  const uint cell_na = ((sys_cdims >> 28) & 0xfff);
  const uint cell_nb = ((sys_cdims >> 40) & 0xfff);
  const uint cell_nc = (sys_cdims >> 52);
  uint image_idx;
  switch (orientation) {
  case CoordinateCycle::BLACK:
    image_idx = image_array_indices_alt.readHost(synth_idx);
    break;
  case CoordinateCycle::WHITE:
    image_idx = image_array_indices.readHost(synth_idx);
    break;
  }
  const int chn_low_bnd = system_chain_bounds.readHost(system_index);
  const uint rel_image_idx = image_idx - image_chain_limits.readHost(chn_low_bnd);
  const uint chain_length = cell_na * cell_base_capacity;
  const int chn_rel_idx = (rel_image_idx / chain_length);
  const int chain_cpos = chn_rel_idx / cell_nb;
  const int chain_bpos = chn_rel_idx - (chain_cpos * cell_nb);
  const int chain_idx = chn_rel_idx + chn_low_bnd;
  const int chain_loc = rel_image_idx - (chain_length * chn_rel_idx);
  int low_cell_est = cell_init + (chn_rel_idx * cell_na);
  int hgh_cell_est = cell_init + ((chn_rel_idx + 1) * cell_na);
  int mid_cell_est = low_cell_est + ((hgh_cell_est - low_cell_est) >> 1);
  bool found;
  do {
    uint2 mid_lims;
    switch (orientation) {
    case CoordinateCycle::BLACK:
      mid_lims = image_cell_limits_alt.readHost(mid_cell_est);
      break;
    case CoordinateCycle::WHITE:
      mid_lims = image_cell_limits.readHost(mid_cell_est);
      break;
    }
    const uint mid_min = mid_lims.x;
    const uint mid_max = mid_lims.x + (mid_lims.y >> 16);
    if (image_idx < mid_min) {
      hgh_cell_est = mid_cell_est;
      found = false;
    }
    else if (image_idx < mid_max) {
      found = true;
    }
    else {
      low_cell_est = mid_cell_est;
      found = false;
    }
    mid_cell_est = ((low_cell_est + hgh_cell_est) >> 1);
  } while (! found);
  const int chain_apos = mid_cell_est - (cell_init + (chn_rel_idx * cell_na));
  return { chain_apos, chain_bpos, chain_cpos, mid_cell_est };
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int4 CellGrid<T, Tacc, Tcalc, T4>::getCellLocation(const int system_index,
                                                   const int atom_index) const {
  return getCellLocation(system_index, atom_index, cycle_position);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int2 CellGrid<T, Tacc, Tcalc, T4>::getChainLocation(const int system_index, const int atom_index,
                                                    const CoordinateCycle orientation) const {
  validateAtomIndex(system_index, atom_index, "getCellLocation");
  const int synth_idx = poly_ps_ptr->getAtomOffset(system_index) + atom_index;
  const ullint sys_cdims = system_cell_grids.readHost(system_index);
  const uint cell_init = (sys_cdims & 0xfffffff);
  const uint cell_na = ((sys_cdims >> 28) & 0xfff);
  const uint cell_nb = ((sys_cdims >> 40) & 0xfff);
  const uint cell_nc = (sys_cdims >> 52);
  uint image_idx;
  switch (orientation) {
  case CoordinateCycle::BLACK:
    image_idx = image_array_indices_alt.readHost(synth_idx);
    break;
  case CoordinateCycle::WHITE:
    image_idx = image_array_indices.readHost(synth_idx);
    break;
  }
  const int chn_low_bnd = system_chain_bounds.readHost(system_index);
  const uint rel_image_idx = image_idx - image_chain_limits.readHost(chn_low_bnd);
  const uint chain_length = cell_na * cell_base_capacity;
  const int chn_rel_idx = (rel_image_idx / chain_length);
  const int chain_idx = chn_rel_idx + chn_low_bnd;
  const int chain_loc = rel_image_idx - (chain_length * chn_rel_idx);
  return { chain_idx, chain_loc };
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
int2 CellGrid<T, Tacc, Tcalc, T4>::getChainLocation(const int system_index,
                                                    const int atom_index) const {
  return getChainLocation(system_index, atom_index, cycle_position);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const PhaseSpaceSynthesis* CellGrid<T, Tacc, Tcalc, T4>::getCoordinateSynthesisPointer() const {
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const AtomGraphSynthesis* CellGrid<T, Tacc, Tcalc, T4>::getTopologySynthesisPointer() const {
  return poly_ag_ptr;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4>
CellGrid<T, Tacc, Tcalc, T4>::data(const CoordinateCycle orientation,
                                   const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return CellGridWriter<T, Tacc,
                          Tcalc, T4>(theme, system_count, total_cell_count, total_chain_count,
                                     mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                     localpos_scale, localpos_inverse_scale,
                                     system_cell_grids.data(tier), system_cell_umat_alt.data(tier),
                                     system_cell_invu_alt.data(tier),
                                     system_pmig_invu_alt.data(tier),
                                     image_cell_limits_alt.data(tier),
                                     image_cell_limits.data(tier), image_chain_limits.data(tier),
                                     system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image_alt.data(tier),
                                     image.data(tier), image_chain_entering.data(tier),
                                     image_chain_leaving.data(tier),
                                     image_vestibule_limits.data(tier),
                                     image_chain_entering_counts.data(tier),
                                     image_chain_leaving_counts.data(tier),
                                     nonimaged_atom_indices_alt.data(tier),
                                     image_array_indices_alt.data(tier),
                                     exclusion_maps_alt.data(tier), exclusion_maps.data(tier),
                                     nt_work_groups.data(tier), x_force.data(tier),
                                     y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier), warp_x_work.data(tier),
                                     warp_y_work.data(tier), warp_z_work.data(tier),
                                     warp_x_overflow_work.data(tier),
                                     warp_y_overflow_work.data(tier),
                                     warp_z_overflow_work.data(tier));
  case CoordinateCycle::WHITE:
    return CellGridWriter<T, Tacc,
                          Tcalc, T4>(theme, system_count, total_cell_count, total_chain_count,
                                     mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                     localpos_scale, localpos_inverse_scale,
                                     system_cell_grids.data(tier), system_cell_umat.data(tier),
                                     system_cell_invu.data(tier), system_pmig_invu.data(tier),
                                     image_cell_limits.data(tier),
                                     image_cell_limits_alt.data(tier),
                                     image_chain_limits.data(tier), system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image.data(tier),
                                     image_alt.data(tier), image_chain_entering.data(tier),
                                     image_chain_leaving.data(tier),
                                     image_vestibule_limits.data(tier),
                                     image_chain_entering_counts.data(tier),
                                     image_chain_leaving_counts.data(tier),
                                     nonimaged_atom_indices.data(tier),
                                     image_array_indices.data(tier), exclusion_maps.data(tier),
                                     exclusion_maps_alt.data(tier), nt_work_groups.data(tier),
                                     x_force.data(tier), y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier), warp_x_work.data(tier),
                                     warp_y_work.data(tier), warp_z_work.data(tier),
                                     warp_x_overflow_work.data(tier),
                                     warp_y_overflow_work.data(tier),
                                     warp_z_overflow_work.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4>
CellGrid<T, Tacc, Tcalc, T4>::data(const HybridTargetLevel tier) {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
CellGrid<T, Tacc, Tcalc, T4>::data(const CoordinateCycle orientation,
                                   const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return CellGridReader<T, Tacc,
                          Tcalc, T4>(theme, system_count, total_cell_count, total_chain_count,
                                     mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                     localpos_scale, localpos_inverse_scale,
                                     system_cell_grids.data(tier), system_cell_umat_alt.data(tier),
                                     system_cell_invu_alt.data(tier),
                                     system_pmig_invu_alt.data(tier),
                                     image_cell_limits_alt.data(tier),
                                     image_chain_limits.data(tier), system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image_alt.data(tier),
                                     nonimaged_atom_indices_alt.data(tier),
                                     image_array_indices_alt.data(tier),
                                     exclusion_maps_alt.data(tier), nt_work_groups.data(tier),
                                     x_force.data(tier), y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier));
  case CoordinateCycle::WHITE:
    return CellGridReader<T, Tacc,
                          Tcalc, T4>(theme, system_count, total_cell_count, total_chain_count,
                                     mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                     localpos_scale, localpos_inverse_scale,
                                     system_cell_grids.data(tier), system_cell_umat.data(tier),
                                     system_cell_invu.data(tier), system_pmig_invu.data(tier),
                                     image_cell_limits.data(tier), image_chain_limits.data(tier),
                                     system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image.data(tier),
                                     nonimaged_atom_indices.data(tier),
                                     image_array_indices.data(tier), exclusion_maps.data(tier),
                                     nt_work_groups.data(tier), x_force.data(tier),
                                     y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
CellGrid<T, Tacc, Tcalc, T4>::data(const HybridTargetLevel tier) const {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<void, void, void, void>
CellGrid<T, Tacc, Tcalc, T4>::templateFreeData(const CoordinateCycle orientation,
                                               const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return CellGridWriter<void, void,
                          void, void>(theme, system_count, total_cell_count, total_chain_count,
                                      mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                      localpos_scale, localpos_inverse_scale,
                                      system_cell_grids.data(tier),
                                      reinterpret_cast<void*>(system_cell_umat_alt.data(tier)),
                                      reinterpret_cast<void*>(system_cell_invu_alt.data(tier)),
                                      reinterpret_cast<void*>(system_pmig_invu_alt.data(tier)),
                                      image_cell_limits_alt.data(tier),
                                      image_cell_limits.data(tier), image_chain_limits.data(tier),
                                      system_chain_bounds.data(tier),
                                      chain_system_membership.data(tier),
                                      reinterpret_cast<void*>(image_alt.data(tier)),
                                      reinterpret_cast<void*>(image.data(tier)),
                                      image_chain_entering.data(tier),
                                      image_chain_leaving.data(tier),
                                      image_vestibule_limits.data(tier),
                                      image_chain_entering_counts.data(tier),
                                      image_chain_leaving_counts.data(tier),
                                      nonimaged_atom_indices_alt.data(tier),
                                      image_array_indices_alt.data(tier),
                                      exclusion_maps.data(tier), exclusion_maps_alt.data(tier),
                                      nt_work_groups.data(tier),
                                      reinterpret_cast<void*>(x_force.data(tier)),
                                      reinterpret_cast<void*>(y_force.data(tier)),
                                      reinterpret_cast<void*>(z_force.data(tier)),
                                      x_force_overflow.data(tier), y_force_overflow.data(tier),
                                      z_force_overflow.data(tier),
                                      reinterpret_cast<void*>(warp_x_work.data(tier)),
                                      reinterpret_cast<void*>(warp_y_work.data(tier)),
                                      reinterpret_cast<void*>(warp_z_work.data(tier)),
                                      warp_x_overflow_work.data(tier),
                                      warp_y_overflow_work.data(tier),
                                      warp_z_overflow_work.data(tier));
  case CoordinateCycle::WHITE:
    return CellGridWriter<void, void,
                          void, void>(theme, system_count, total_cell_count, total_chain_count,
                                      mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                      localpos_scale, localpos_inverse_scale,
                                      system_cell_grids.data(tier),
                                      reinterpret_cast<void*>(system_cell_umat.data(tier)),
                                      reinterpret_cast<void*>(system_cell_invu.data(tier)),
                                      reinterpret_cast<void*>(system_pmig_invu.data(tier)),
                                      image_cell_limits.data(tier),
                                      image_cell_limits_alt.data(tier),
                                      image_chain_limits.data(tier),
                                      system_chain_bounds.data(tier),
                                      chain_system_membership.data(tier),
                                      reinterpret_cast<void*>(image.data(tier)),
                                      reinterpret_cast<void*>(image_alt.data(tier)),
                                      image_chain_entering.data(tier),
                                      image_chain_leaving.data(tier),
                                      image_vestibule_limits.data(tier),
                                      image_chain_entering_counts.data(tier),
                                      image_chain_leaving_counts.data(tier),
                                      nonimaged_atom_indices.data(tier),
                                      image_array_indices.data(tier), exclusion_maps.data(tier),
                                      exclusion_maps_alt.data(tier), nt_work_groups.data(tier),
                                      reinterpret_cast<void*>(x_force.data(tier)),
                                      reinterpret_cast<void*>(y_force.data(tier)),
                                      reinterpret_cast<void*>(z_force.data(tier)),
                                      x_force_overflow.data(tier), y_force_overflow.data(tier),
                                      z_force_overflow.data(tier),
                                      reinterpret_cast<void*>(warp_x_work.data(tier)),
                                      reinterpret_cast<void*>(warp_y_work.data(tier)),
                                      reinterpret_cast<void*>(warp_z_work.data(tier)),
                                      warp_x_overflow_work.data(tier),
                                      warp_y_overflow_work.data(tier),
                                      warp_z_overflow_work.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<void, void, void, void>
CellGrid<T, Tacc, Tcalc, T4>::templateFreeData(const HybridTargetLevel tier) {
  return templateFreeData(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<void, void, void, void>
CellGrid<T, Tacc, Tcalc, T4>::templateFreeData(const CoordinateCycle orientation,
                                               const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    {
      const void* cell_umat_ptr = reinterpret_cast<const void*>(system_cell_umat_alt.data(tier));
      const void* cell_invu_ptr = reinterpret_cast<const void*>(system_cell_invu_alt.data(tier));
      const void* pmig_invu_ptr = reinterpret_cast<const void*>(system_pmig_invu_alt.data(tier));
      return CellGridReader<void, void,
                            void, void>(theme, system_count, total_cell_count, total_chain_count,
                                        mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                        localpos_scale, localpos_inverse_scale,
                                        system_cell_grids.data(tier), cell_umat_ptr, cell_invu_ptr,
                                        pmig_invu_ptr, image_cell_limits_alt.data(tier),
                                        image_chain_limits.data(tier),
                                        system_chain_bounds.data(tier),
                                        chain_system_membership.data(tier),
                                        reinterpret_cast<const void*>(image_alt.data(tier)),
                                        nonimaged_atom_indices_alt.data(tier),
                                        image_array_indices_alt.data(tier),
                                        exclusion_maps_alt.data(tier), nt_work_groups.data(tier),
                                        reinterpret_cast<const void*>(x_force.data(tier)),
                                        reinterpret_cast<const void*>(y_force.data(tier)),
                                        reinterpret_cast<const void*>(z_force.data(tier)),
                                        x_force_overflow.data(tier), y_force_overflow.data(tier),
                                        z_force_overflow.data(tier));
    }
    break;
  case CoordinateCycle::WHITE:
    return CellGridReader<void, void,
                          void, void>(theme, system_count, total_cell_count, total_chain_count,
                                      mesh_subdivisions, cell_base_capacity, cell_excl_capacity,
                                      localpos_scale, localpos_inverse_scale,
                                      system_cell_grids.data(tier),
                                      reinterpret_cast<const void*>(system_cell_umat.data(tier)),
                                      reinterpret_cast<const void*>(system_cell_invu.data(tier)),
                                      reinterpret_cast<const void*>(system_pmig_invu.data(tier)),
                                      image_cell_limits.data(tier),
                                      image_chain_limits.data(tier),
                                      system_chain_bounds.data(tier),
                                      chain_system_membership.data(tier),
                                      reinterpret_cast<const void*>(image.data(tier)),
                                      nonimaged_atom_indices.data(tier),
                                      image_array_indices.data(tier), exclusion_maps.data(tier),
                                      nt_work_groups.data(tier),
                                      reinterpret_cast<const void*>(x_force.data(tier)),
                                      reinterpret_cast<const void*>(y_force.data(tier)),
                                      reinterpret_cast<const void*>(z_force.data(tier)),
                                      x_force_overflow.data(tier),
                                      y_force_overflow.data(tier), z_force_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<void, void, void, void>
CellGrid<T, Tacc, Tcalc, T4>::templateFreeData(const HybridTargetLevel tier) const {
  return templateFreeData(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGrid<T, Tacc, Tcalc, T4>* CellGrid<T, Tacc, Tcalc, T4>::getSelfPointer() const {
  return this;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::upload() {
  system_cell_grids.upload();
  system_cell_umat.upload();
  system_cell_umat_alt.upload();
  system_cell_invu.upload();
  system_cell_invu_alt.upload();
  system_pmig_invu.upload();
  system_pmig_invu_alt.upload();
  image_cell_limits.upload();
  image_cell_limits_alt.upload();
  image_chain_limits.upload();
  system_chain_bounds.upload();
  chain_system_membership.upload();
  image.upload();
  image_alt.upload();
  image_chain_entering.upload();
  image_chain_leaving.upload();
  image_vestibule_limits.upload();
  image_chain_entering_counts.upload();
  image_chain_leaving_counts.upload();
  nonimaged_atom_indices.upload();
  nonimaged_atom_indices_alt.upload();
  image_array_indices.upload();
  image_array_indices_alt.upload();
  exclusion_maps.upload();
  exclusion_maps_alt.upload();
  nt_work_groups.upload();
  x_force.upload();
  y_force.upload();
  z_force.upload();
  warp_x_work.upload();
  warp_y_work.upload();
  warp_z_work.upload();
  x_force_overflow.upload();
  y_force_overflow.upload();
  z_force_overflow.upload();
  warp_x_overflow_work.upload();
  warp_y_overflow_work.upload();
  warp_z_overflow_work.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::download() {
  system_cell_grids.download();
  system_cell_umat.download();
  system_cell_umat_alt.download();
  system_cell_invu.download();
  system_cell_invu_alt.download();
  system_pmig_invu.download();
  system_pmig_invu_alt.download();
  image_cell_limits.download();
  image_cell_limits_alt.download();
  image_chain_limits.download();
  system_chain_bounds.download();
  chain_system_membership.download();
  image.download();
  image_alt.download();
  image_chain_entering.download();
  image_chain_leaving.download();
  image_vestibule_limits.download();
  image_chain_entering_counts.download();
  image_chain_leaving_counts.download();
  nonimaged_atom_indices.download();
  nonimaged_atom_indices_alt.download();
  image_array_indices.download();
  image_array_indices_alt.download();
  exclusion_maps.download();
  exclusion_maps_alt.download();
  nt_work_groups.download();
  x_force.download();
  y_force.download();
  z_force.download();
  warp_x_work.download();
  warp_y_work.download();
  warp_z_work.download();
  x_force_overflow.download();
  y_force_overflow.download();
  z_force_overflow.download();
  warp_x_overflow_work.download();
  warp_y_overflow_work.download();
  warp_z_overflow_work.download();
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::initializeForces(const HybridTargetLevel tier,
                                                    const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const Tacc value_zero = 0;
      CellGridWriter<T, Tacc, Tcalc, T4> cgw = data(cycle_position);
      for (int i = 0; i < cgw.total_cell_count; i++) {
        const uint2 cell_lims = cgw.cell_limits[i];
        const uint hlim = cell_lims.x + (cell_lims.y >> 16);
        for (uint j = cell_lims.x; j < hlim; j++) {
          cgw.xfrc[j] = value_zero;
          cgw.yfrc[j] = value_zero;
          cgw.zfrc[j] = value_zero;
          cgw.xfrc_ovrf[j] = 0;
          cgw.yfrc_ovrf[j] = 0;
          cgw.zfrc_ovrf[j] = 0;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    if (gpu != null_gpu) {
      CellGridWriter<void, void, void, void> cgw = templateFreeData(cycle_position, tier);
      launchCellGridAction(&cgw, std::type_index(typeid(T)).hash_code(),
                           std::type_index(typeid(Tacc)).hash_code(), gpu,
                           CellGridAction::INIT_FORCES);
    }
    else {
      rtErr("A valid GPU must be provided.", "CellGrid", "initializeForces");
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::contributeForces(PhaseSpaceSynthesis *dest,
                                                    const HybridTargetLevel tier,
                                                    const GpuDetails &gpu) const {

  // Check that the destination synthesis and the internally referenced object have identical
  // fixed precision force accumulation.
  if (dest->getForceAccumulationBits() != poly_ps_ptr->getForceAccumulationBits()) {
    rtErr("The force scaling model in the provided coordinate synthesis does not match that of "
          "the internally referenced synthesis (" +
          std::to_string(dest->getForceAccumulationBits()) + " vs. " +
          std::to_string(poly_ps_ptr->getForceAccumulationBits()) + " bits of unit precision).",
          "CellGrid", "contributeForces");
  }
  
  // The coordinate synthesis will receive forces at its own current cycle position.  Taking the
  // abstract with no cycle position argument defaults to the PhaseSpaceSynthesis's own setting.
  PsSynthesisWriter destw = dest->data(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      CellGridWriter<T, Tacc, Tcalc, T4> cgw = data(cycle_position);
      for (int pos = 0; pos < destw.system_count; pos++) {
        if (destw.atom_counts[pos] != poly_ps_ptr->getAtomCount(pos)) {
          rtErr("The atom counts for system index " + std::to_string(pos) + " do not match "
                "between the destination and internally referenced coordinate syntheses (" +
                std::to_string(destw.atom_counts[pos]) + ", " +
                std::to_string(poly_ps_ptr->getAtomCount(pos)) + ").", "CellGrid",
                "contributeForces");
        }
        if (destw.atom_starts[pos] != poly_ps_ptr->getAtomOffset(pos)) {
          rtErr("The atom offsets for system index " + std::to_string(pos) + " do not match "
                "between the destination and internally referenced coordinate syntheses (" +
                std::to_string(destw.atom_starts[pos]) + ", " +
                std::to_string(poly_ps_ptr->getAtomOffset(pos)) + ").", "CellGrid",
                "contributeForces");
        }
        const int hlim = destw.atom_starts[pos] + destw.atom_counts[pos];
        if (destw.frc_bits <= force_scale_nonoverflow_bits) {
          for (int i = destw.atom_starts[pos]; i < hlim; i++) {
            const uint image_idx = cgw.img_atom_idx[i];
            destw.xfrc[i] += cgw.xfrc[image_idx];
            destw.yfrc[i] += cgw.yfrc[image_idx];
            destw.zfrc[i] += cgw.zfrc[image_idx];
          }
        }
        else {
          for (int i = destw.atom_starts[pos]; i < hlim; i++) {
            const uint image_idx = cgw.img_atom_idx[i];
            const int95_t nfx = hostInt95Sum(destw.xfrc[i], destw.xfrc_ovrf[i],
                                             cgw.xfrc[image_idx], cgw.xfrc_ovrf[image_idx]);
            const int95_t nfy = hostInt95Sum(destw.yfrc[i], destw.yfrc_ovrf[i],
                                             cgw.yfrc[image_idx], cgw.yfrc_ovrf[image_idx]);
            const int95_t nfz = hostInt95Sum(destw.zfrc[i], destw.zfrc_ovrf[i],
                                             cgw.zfrc[image_idx], cgw.zfrc_ovrf[image_idx]);
            destw.xfrc[i] = nfx.x;
            destw.yfrc[i] = nfy.x;
            destw.zfrc[i] = nfz.x;
            destw.xfrc_ovrf[i] = nfx.y;
            destw.yfrc_ovrf[i] = nfy.y;
            destw.zfrc_ovrf[i] = nfz.y;
          }
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    if (gpu != null_gpu) {
      const CellGridReader<void, void, void, void> cgr = templateFreeData(cycle_position, tier);
      launchCellGridAction(cgr, std::type_index(typeid(T)).hash_code(),
                           std::type_index(typeid(Tacc)).hash_code(), &destw, gpu,
                           CellGridAction::XFER_FORCES);
    }
    else {
      rtErr("A valid GPU must be provided.", "CellGrid", "initializeForces");
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::contributeForces(const HybridTargetLevel tier,
                                                    const GpuDetails &gpu) const {
  contributeForces(poly_ps_ptr, tier, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::validateCoordinateSynthesis() const {
  if (poly_ps_ptr->getSystemCount() >= maximum_cellgrid_systems) {
    rtErr("A system count of " + std::to_string(poly_ps_ptr->getSystemCount()) + " is too large.  "
          "A maximum of " + std::to_string(maximum_cellgrid_systems) + " periodic simulations are "
          "supported at one time.", "CellGrid", "validateCoordinateSynthesis");
  }
  switch (poly_ps_ptr->getUnitCellType()) {
  case UnitCellType::NONE:
    rtErr("Periodic boundary conditions must apply to all systems.", "CellGrid",
          "validateCoordinateSynthesis");
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
double
CellGrid<T, Tacc, Tcalc, T4>::validateEffectiveCutoff(const double eff_cut_in,
                                                      const ExceptionResponse policy_in) const {
  if (eff_cut_in < minimum_cell_width) {
    switch (policy_in) {
    case ExceptionResponse::DIE:
      rtErr("An effective cutoff of " +
            realToString(eff_cut_in, 8, 4, NumberFormat::STANDARD_REAL) + " is too short and may "
            "result in missed bonded exclusions.", "CellGrid", "validateEffectiveCutoff");
    case ExceptionResponse::WARN:
      rtWarn("An effective cutoff of " +
             realToString(eff_cut_in, 8, 4, NumberFormat::STANDARD_REAL) + " is too short and may "
             "result in missed bonded exclusions.  The minimum value of " +
             realToString(minimum_cell_width, 8, 4, NumberFormat::STANDARD_REAL) +
             " will be applied instead.", "CellGrid", "validateEffectiveCutoff");
      return minimum_cell_width;
    case ExceptionResponse::SILENT:
      return minimum_cell_width;
    }
  }
  else if (eff_cut_in >= maximum_cell_width) {
    switch (policy_in) {
    case ExceptionResponse::DIE:
      rtErr("An effective cutoff of " +
            realToString(eff_cut_in, 8, 4, NumberFormat::STANDARD_REAL) + " is too long and may "
            "result in cells with too many particles.", "CellGrid", "validateEffectiveCutoff");
    case ExceptionResponse::WARN:
      rtWarn("An effective cutoff of " +
             realToString(eff_cut_in, 8, 4, NumberFormat::STANDARD_REAL) + " is too long and may "
             "result in cells with too many particles.  The maximum value of " +
             realToString(maximum_cell_width, 8, 4, NumberFormat::STANDARD_REAL) + "will be "
             "applied instead.", "CellGrid", "validateEffectiveCutoff");
      return maximum_cell_width;
    case ExceptionResponse::SILENT:
      return maximum_cell_width;
    }
  }
  return eff_cut_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::validateAtomIndex(const int system_index, const int atom_index,
                                                     const char* caller) const {
  if (system_index < 0 || system_index >= poly_ps_ptr->getSystemCount()) {
    rtErr("System index " + std::to_string(system_index) + " is ivalid for a synthesis of " +
          std::to_string(poly_ps_ptr->getSystemCount()) + " systems.", "CellGrid", caller);
  }
  if (atom_index < 0 || atom_index >= poly_ps_ptr->getAtomCount(system_index)) {
    rtErr("Atom index " + std::to_string(atom_index) + " is invalid for system index " +
          std::to_string(system_index) + " (total atoms " +
          std::to_string(poly_ps_ptr->getAtomCount(system_index)) + ".", "CellGrid",
          caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void
CellGrid<T, Tacc, Tcalc, T4>::computeFixedPrecisionModel(const std::vector<double> &invu_samples) {

  // Set the position bits to zero for floating-point coordinate representations
  if (isFloatingPointScalarType<T>()) {
    localpos_scale_bits = 0;
    localpos_scale = 1.0;
    localpos_inverse_scale = 1.0;
    return;
  }
  else {

    // Compute the largest decomposition cell across all systems.
    double max_x = 0.0;
    double max_y = 0.0;
    double max_z = 0.0;
    double min_x = 0.0;
    double min_y = 0.0;
    double min_z = 0.0;
    const double* invu_ptr = invu_samples.data();
    const int xfrm_stride = roundUp(9, warp_size_int);
    for (int pos = 0; pos < system_count; pos++) {
      const double* p_invu = &invu_ptr[pos * xfrm_stride];
      for (int i = -2; i <= 3; i++) {
        const double di = i;
        for (int j = -2; j <= 3; j++) {
          const double dj = j;
          for (int k = -2; k <= 3; k++) {
            const double dk = k;
            const double x = (p_invu[0] * di) + (p_invu[3] * dj) + (p_invu[6] * dk);
            const double y =                    (p_invu[4] * dj) + (p_invu[7] * dk);
            const double z =                                       (p_invu[8] * dk);
            max_x = std::max(x, max_x);
            max_y = std::max(y, max_y);
            max_z = std::max(z, max_z);
            min_x = std::min(x, min_x);
            min_y = std::min(y, min_y);
            min_z = std::min(z, min_z);
          }
        }
      }
    }

    // The fixed precision model must accommodate the maximum dimensions
    double max_xyz = std::max(fabs(max_x), fabs(min_x));
    max_xyz = std::max(max_xyz, std::max(fabs(max_y), fabs(min_y)));
    max_xyz = std::max(max_xyz, std::max(fabs(max_z), fabs(min_z)));
    if (max_xyz < constants::tiny) {
      rtErr("No unit cell appears to have any dimensions.", "CellGrid",
            "computeFixedPrecisionModel");
    }
    const double ltwo_xyz = log2(max_xyz);
    const int bits_ahead = ceil(ltwo_xyz) + static_cast<double>(ceil(ltwo_xyz) - ltwo_xyz < 0.2);
    localpos_scale_bits = (sizeof(T) * 8) - bits_ahead - 1;
    localpos_scale = pow(2.0, localpos_scale_bits);
    localpos_inverse_scale = 1.0 / localpos_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template<typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::tallyCellPopulations(std::vector<int> *cell_populations,
                                                        const int system_index, const int na,
                                                        const int nb, const int nc) {
  CoordinateFrame cf = poly_ps_ptr->exportCoordinates(system_index);
  CoordinateFrameWriter cfw = cf.data();
  const AtomGraph *ag_ptr = poly_ps_ptr->getSystemTopologyPointer(system_index);
  const NonbondedKit nbk = ag_ptr->getDoublePrecisionNonbondedKit();
  const VdwCombiningRule lj_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff,
                                                              nbk.n_lj_types);
  std::vector<int> cell_dest(cfw.natom, 0);

  // Intercept the fractional coordinates and image them into the primary unit cell, then
  // calculate the spatial decomposition cell before returning to real space.
  locateDecompositionCell(cfw.xcrd, cfw.ycrd, cfw.zcrd, cell_dest.data(), cfw.natom, cfw.umat,
                          cfw.invu, na, nb, nc, cfw.unit_cell);
  const int icell_count = na * nb * nc;
  cell_populations->resize(icell_count);
  int* pop_ptr = cell_populations->data();
  for (int i = 0; i < icell_count; i++) {
    pop_ptr[i] = 0;
  }
  for (int i = 0; i < cfw.natom; i++) {
    switch (theme) {
    case NonbondedTheme::ELECTROSTATIC:
      if (fabs(nbk.charge[i]) > constants::small) {
        pop_ptr[cell_dest[i]] += 1;
      }
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      if (hasVdwProperties<double>(nbk, i, lj_rule)) {
        pop_ptr[cell_dest[i]] += 1;
      }
      break;
    case NonbondedTheme::ALL:
      pop_ptr[cell_dest[i]] += 1;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::populateImage(const CoordinateCycle cyc) {

  // Recall the type setting for coordinates--this will determine the manner in which data is
  // translated from the coordinate frame to the local axes, as well as atom properties from the
  // topology to the final members of the image tuples.
  const size_t tc_crd = std::type_index(typeid(T4)).hash_code();
  
  // Set pointers for the image of interest
  const PsSynthesisReader poly_psr = poly_ps_ptr->data(cyc);
  const SyNonbondedKit<double, double2> poly_nbk = poly_ag_ptr->getDoublePrecisionNonbondedKit();
  uint2* cell_limits_ptr;
  T4* image_ptr;
  int* nonimg_atom_idx_ptr;
  uint* image_idx_ptr;
  switch (cyc) {
  case CoordinateCycle::WHITE:
    cell_limits_ptr = image_cell_limits.data();
    image_ptr = image.data();
    nonimg_atom_idx_ptr = nonimaged_atom_indices.data();
    image_idx_ptr = image_array_indices.data();
    break;
  case CoordinateCycle::BLACK:
    cell_limits_ptr = image_cell_limits.data();
    image_ptr = image_alt.data();
    nonimg_atom_idx_ptr = nonimaged_atom_indices_alt.data();
    image_idx_ptr = image_array_indices_alt.data();
    break;
  }

  // Loop over all systems, reallocating temporary arrays and filling the pointers set above
  std::vector<int> cell_dest, cell_contents, pcell_limits;
  
  for (int pos = 0; pos < poly_psr.system_count; pos++) {
    const ullint pcell_dims = system_cell_grids.readHost(pos);
    const int pcell_offset = (pcell_dims & 0xfffffffLLU);
    const int cell_na = ((pcell_dims >> 28) & 0xfffLLU);
    const int cell_nb = ((pcell_dims >> 40) & 0xfffLLU);
    const int cell_nc = ((pcell_dims >> 52) & 0xfffLLU);
    CoordinateFrame cf = poly_ps_ptr->exportCoordinates(pos, cyc);
    CoordinateFrameWriter cfw = cf.data();

    // Check that the coordinates obtained from the synthesis at this point in the time cycle are
    // valid.  If all coordinates are zero, then the synthesis was probably not initialized with
    // real data at this point in its own time cycle.  Swap to the alternate time point in such a
    // case to have valid coordinates and avoid overflowing the cell containing the origin of each
    // system's coordinates.
    if (cfw.natom >= 2 && mean(cfw.xcrd, cfw.natom) < constants::small &&
                          mean(cfw.ycrd, cfw.natom) < constants::small &&
                          mean(cfw.zcrd, cfw.natom) < constants::small) {
      switch (cyc) {
      case CoordinateCycle::WHITE:
        coordCopy(&cf, *poly_ps_ptr, pos, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK);
        break;
      case CoordinateCycle::BLACK:
        coordCopy(&cf, *poly_ps_ptr, pos, TrajectoryKind::POSITIONS, CoordinateCycle::WHITE);
        break;
      }
    }
    const AtomGraph *ag_ptr = poly_ps_ptr->getSystemTopologyPointer(pos);
    const NonbondedKit nbk = ag_ptr->getDoublePrecisionNonbondedKit();
    const VdwCombiningRule lj_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff,
                                                                nbk.n_lj_types);

    // Intercept the fractional coordinates and image them into the primary unit cell, then
    // calculate the spatial decomposition cell before returning to real space.  On this second
    // pass, the decomposition occurs with both time cycle points from the original coordinate
    // synthesis.  Furthermore, the points in the time cycle are explicily spelled out, rather than
    // just taking the current time cycle position from the coordinate synthesis as was done in the
    // first pass.
    cell_dest.resize(cfw.natom);
    locateDecompositionCell(cfw.xcrd, cfw.ycrd, cfw.zcrd, cell_dest.data(), cfw.natom, cfw.umat,
                            cfw.invu, cell_na, cell_nb, cell_nc, cfw.unit_cell);
    
    // Place atoms in their destinations.  The information can now be added to the image array.
    // Fill the appropriate image, non-imaged atom indices, and cell limits in the synthesis-wide
    // cell grid bounds array.
    const int patom_offset = poly_psr.atom_starts[pos];
    const int pcell_count  = cell_na * cell_nb * cell_nc;
    pcell_limits.resize(pcell_count + 1);
    cell_contents.resize(cfw.natom);
    indexingArray(cell_dest, &cell_contents, &pcell_limits, pcell_count);
    const double dcell_na = cell_na;
    const double dcell_nb = cell_nb;
    const double dcell_nc = cell_nc;
    for (int j = 0; j < cell_nb; j++) {
      const double cp_b = static_cast<double>(j) / dcell_nb;
      for (int k = 0; k < cell_nc; k++) {
        const int chain_idx = system_chain_bounds.readHost(pos) + (k * cell_nb) + j;
        const uint chain_start = image_chain_limits.readHost(chain_idx);
        uint chain_pos = 0;
        const double cp_c = static_cast<double>(k) / dcell_nc;
        for (int i = 0; i < cell_na; i++) {

          // Compute the origin of cell {i, j, k}
          const double cp_a = static_cast<double>(i) / dcell_na;
          const double orig_x = (cfw.invu[0] * cp_a) + (cfw.invu[3] * cp_b) + (cfw.invu[6] * cp_c);
          const double orig_y =                        (cfw.invu[4] * cp_b) + (cfw.invu[7] * cp_c);
          const double orig_z =                                               (cfw.invu[8] * cp_c);

          // The cell limits of the ith cell in this, the {j, k}th chain of system pos, begin at
          // the current local chain position plus the chain offset in the image.
          uint2 tmp_cell_limits;
          tmp_cell_limits.x = chain_start + chain_pos;
          
          // Compute the indexing limits of the cell within the locally prepared cell_contents
          // array and loop over these limits with the counter m.  Compute the location of each
          // atom relative to the cell origin, scale as appropriate, and store the result.
          const size_t cell_ijk_idx = (((k * cell_nb) + j) * cell_na) + i;
          const uint mllim = pcell_limits[cell_ijk_idx];
          const uint mhlim = pcell_limits[cell_ijk_idx + 1];
          for (int m = mllim; m < mhlim; m++) {

            // As obtained here, the "topological atom index" refers to the atom's place in the one
            // system.  The cell_contents array is local to this routine, created and resized as
            // needed for each system of the synthesis.  As such the topological atom index is
            // useful for lookups within the specific system's AtomGraph.  However, at the end of
            // this process, the "non-imaged atom index" must be recorded with respect to the
            // system's atom offset within the synthesis.  This entire process implies that the
            // indexing of atoms in each system of the synthesis must not change order with respect
            // to the original AtomGraph objects, but there are other aspects of the code that have
            // implicit reliance on this and it is not restrictive to keep the ordering intact.
            const int topl_atom_idx = cell_contents[m];

            // Form the coordinate / property tuple and place this in the proper image.  Due to the
            // consistency of the ordering in individual topologies and the topology synthesis, the
            // single topology's abstract can be used to make these distinctions.
            bool valid_atom;
            switch (theme) {
            case NonbondedTheme::ELECTROSTATIC:
              valid_atom = (fabs(nbk.charge[topl_atom_idx]) > constants::small);
              break;
            case NonbondedTheme::VAN_DER_WAALS:
              valid_atom = hasVdwProperties<double>(nbk, topl_atom_idx, lj_rule);
              break;
            case NonbondedTheme::ALL:
              valid_atom = true;
              break;
            }
            if (valid_atom) {

              // Form the coordinate / property tuple.
              T4 crdv;
              if (localpos_scale_bits > 0) {
                crdv.x = (cfw.xcrd[topl_atom_idx] - orig_x) * localpos_scale;
                crdv.y = (cfw.ycrd[topl_atom_idx] - orig_y) * localpos_scale;
                crdv.z = (cfw.zcrd[topl_atom_idx] - orig_z) * localpos_scale;
              }
              else {
                crdv.x = cfw.xcrd[topl_atom_idx] - orig_x;
                crdv.y = cfw.ycrd[topl_atom_idx] - orig_y;
                crdv.z = cfw.zcrd[topl_atom_idx] - orig_z;
              }
              switch (theme) {
              case NonbondedTheme::ELECTROSTATIC:
                if (tc_crd == int4_type_index) {
                  Ecumenical4 qconv = { .f = static_cast<float>(poly_nbk.charge[patom_offset +
                                                                                topl_atom_idx]) };
                  crdv.w = qconv.i;
                }
                else if (tc_crd == longlong4_type_index) {
                  Ecumenical8 qconv = { .d = poly_nbk.charge[patom_offset + topl_atom_idx] };
                  crdv.w = qconv.lli;
                }
                else {
                  crdv.w = poly_nbk.charge[patom_offset + topl_atom_idx];
                }
                break;
              case NonbondedTheme::VAN_DER_WAALS:

                // To record the Lennard-Jones type index within the local topology would be
                // acceptable, as each system comes with a special type offset (ljabc_offsets in
                // the SyNonbondedKit abstract) which jogs forward to a table in the synthesis
                // which is equivalent to the Lennard-Jones parameter table in the individual
                // topology.  However, for consistency, use the index as it appears in the
                // synthesis array. (While it may seem wasteful to have the synthesis store
                // multiple tables which might contain mostly the same parameters, it is preferable
                // to combining all of the types into one giant table as the size of the matrix
                // would grow as the square of the number of unique Lennard-Jones atom types that
                // possibly interact with one another.  Furthermore, the synthesis does reduce the
                // number of tables it stores to the unique list of tables.)
                if (tc_crd == float4_type_index) {
                  Ecumenical4 tconv = { .i = poly_nbk.lj_idx[patom_offset + topl_atom_idx] };
                  crdv.w = tconv.f;
                }
                else if (tc_crd == double4_type_index) {
                  Ecumenical8 tconv = { .lli = poly_nbk.lj_idx[patom_offset + topl_atom_idx] };
                  crdv.w = tconv.d;
                }
                else {
                  crdv.w = poly_nbk.lj_idx[patom_offset + topl_atom_idx];
                }
                break;
              case NonbondedTheme::ALL:

                // The Lennard-Jones index from the individual topology is valid, and the charge
                // index can likewise be valid by jogging forward using the proper system offset in
                // the topology synthesis.  Record the system-specific charge index alongside the
                // Lennard-Jones index in the bit string.  Charge parameters are not ordered the
                // same way in the synthesis as they are in the underlying individual topologies.
                if (tc_crd == float4_type_index) {
                  const uint tlj_idx = poly_nbk.lj_idx[patom_offset + topl_atom_idx];
                  const uint tc_arg = ((tlj_idx << sp_charge_index_bits) |
                                       poly_nbk.q_idx[patom_offset + topl_atom_idx]);
                  Ecumenical4 tconv = { .ui = tc_arg };
                  crdv.w = tconv.f;
                }
                else if (tc_crd == double4_type_index) {
                  const ullint tlj_idx = poly_nbk.lj_idx[patom_offset + topl_atom_idx];
                  const ullint tc_arg = ((tlj_idx << dp_charge_index_bits) |
                                         poly_nbk.q_idx[patom_offset + topl_atom_idx]);
                  Ecumenical8 tconv = { .ulli = tc_arg };
                  crdv.w = tconv.d;
                }
                else if (tc_crd == int4_type_index) {
                  const uint tlj_idx = poly_nbk.lj_idx[patom_offset + topl_atom_idx];
                  const uint tc_arg = ((tlj_idx << sp_charge_index_bits) |
                                       poly_nbk.q_idx[patom_offset + topl_atom_idx]);
                  const Ecumenical4 tconv = { .ui = tc_arg };
                  crdv.w = tconv.i;
                }
                else {
                  
                  // The format is large enough that the Lennard-Jones type index will not push it
                  // into negative territory.  There is no need to pass through an unsigned integer
                  // intermediate to protect the highest bits across platforms.
                  const llint tlj_idx = poly_nbk.lj_idx[patom_offset + topl_atom_idx];
                  crdv.w = ((tlj_idx << dp_charge_index_bits) |
                            poly_nbk.q_idx[patom_offset + topl_atom_idx]);
                }
                break;
              }

              // Place the atom's coordinate / property tuple in the proper index of image.  Set
              // supplemental arrays at the same index to record the atom's relationship to the
              // original coordinate synthesis.
              const size_t image_idx = chain_start + chain_pos;
              nonimg_atom_idx_ptr[image_idx] = topl_atom_idx + patom_offset;
              image_idx_ptr[topl_atom_idx + patom_offset] = image_idx;
              image_ptr[image_idx] = crdv;
              chain_pos++;
            }
            else {
              image_idx_ptr[topl_atom_idx + patom_offset] = 0xffffffffU;
            }
          }

          // The cell limits of the ith cell in chain {j, k} of this system end at the current
          // chain position.  The possibility that some atoms in the system may not be included in
          // the image (due to lack of electrostatic or van-der Waals properties) forces the tuple
          // to be completed here rather than when it is initialized.  Store the completed tuple.
          const uint natom_covered = chain_start + chain_pos - tmp_cell_limits.x;
          tmp_cell_limits.y = (pos | (natom_covered << 16));
          cell_limits_ptr[pcell_offset + (((k * cell_nb) + j) * cell_na) + i] = tmp_cell_limits;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::prepareWorkGroups() {
  std::vector<int> twu(16);
  std::vector<int> rel_a = {  0,  0,  0,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1,  0,  1,  2 };
  std::vector<int> rel_b = {  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2 };
  std::vector<int> rel_c = { -2, -1,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
  nt_work_groups.resize(static_cast<size_t>(total_cell_count) * 16LLU);
  for (int i = 0; i < total_cell_count; i++) {
    const uint2 icell_lims = image_cell_limits.readHost(i);
    const int system_idx = (icell_lims.y & 0xffff);
    const ullint sys_cg_dims = system_cell_grids.readHost(system_idx);
    const int icell_offset = (sys_cg_dims & 0xfffffff);
    const int icell_na = ((sys_cg_dims >> 28) & 0xfff);
    const int icell_nb = ((sys_cg_dims >> 40) & 0xfff);
    const int icell_nc = ((sys_cg_dims >> 52) & 0xfff);

    // Determine the location of this cell in the system grid overall
    const int abc_idx = i - icell_offset;
    const int c_idx = abc_idx / (icell_na * icell_nb);
    const int b_idx = (abc_idx - (c_idx * icell_na * icell_nb)) / icell_na;
    const int a_idx = abc_idx - (((c_idx * icell_nb) + b_idx) * icell_na);

    // List the other cells that must be brought into play, by their own indices within the grid
    // as a whole.  Place them in a standard order according to the description of the work_groups
    // array, encoded in rel_a, rel_b, and rel_c above.
    for (int j = 0; j < 16; j++) {
      int ra_idx = a_idx + rel_a[j];
      int rb_idx = b_idx + rel_b[j];
      int rc_idx = c_idx + rel_c[j];
      ra_idx += ((ra_idx < 0) + (ra_idx >= icell_na)) * icell_na;
      rb_idx += ((rb_idx < 0) + (rb_idx >= icell_nb)) * icell_nb;
      rc_idx += ((rc_idx < 0) + (rc_idx >= icell_nc)) * icell_nc;
      twu[j] = icell_offset + (((rc_idx * icell_nb) + rb_idx) * icell_na) + ra_idx;
    }
    nt_work_groups.putHost(twu, 16LLU * static_cast<ullint>(i), 16);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::initializeExclusionMasks(const CoordinateCycle cyc) {

  // Determine the necessary size of the array.  Each cell-cell interaction that might hold
  // exclusions to be processed by the non-bonded loop must have a number of bits sufficient to
  // hold the maximum number of atoms it might hold times the maximum number of atoms that the
  // other cell might hold.  Some cells may exceed the maximum atom capacity, but in the same way
  // that pooling the cells in a chain along the system's A axis buffers against catastrophe,
  // pooling the space for exclusion lists covered by fourteen cell:cell interactions protects
  // against catastrophe in tracking exclusions.
  cell_excl_capacity = cell_base_capacity * cellgrid_total_excl_allotment *
                       ((cell_base_capacity / 32) + 1);
  const size_t total_exclusion_capacity = static_cast<size_t>(total_cell_count) *
                                          cell_excl_capacity;
#if 0
  exclusion_maps.resize(total_exclusion_capacity);
  exclusion_maps_alt.resize(total_exclusion_capacity);
#endif
  
  // Loop over all cells, proceeding from a given cell to obtain the system index, the cell grid
  // dimensions, and finally the bounds of all cells in the neutral-territory work unit.  A similar
  // process will be performed by each warp when computing the non-bonded interactions.  The goals
  // are, in order:
  //
  // - Load the cell limits for the home cell
  // - Determine, from the "y" member of the home cell's limits entry, the system index
  // - Load the system's system_cell_grids entry, holding the offset for all cells in the system
  //   and the number of cells it holds along all three of its unit cell axes
  // - Load information pertaining to additional cells of the tower and plate, recording their
  //   offsets as well as the numbers of atoms in each cell
  // - Compute the number of atoms in the whole core of the tower, the core of the plate, and the
  //   lower core of the tower.  These will function as keys into the exclusions array.
  //
  // Once the landmarks are established, lay out the exclusions for that home cell's interactions
  // as an array of blank unsigned integers.  Record the landmarks in temporary arrays so that the
  // exclusions in each system can be sprinkled into the masks in the next stage.
  std::vector<int> tower_whole_core_atoms(total_cell_count);
  std::vector<int> tower_lower_core_atoms(total_cell_count);
  std::vector<int> plate_core_atoms(total_cell_count);
  CellGridWriter<T, Tacc, Tcalc, T4> cgw = data(cyc);
  for (int i = 0; i < total_cell_count; i++) {

    // Load the work group for this home cell
    std::vector<int> twu = nt_work_groups.readHost(static_cast<ullint>(i) * 16LLU, 16);

    // Loop through each topology using the synthesis of forward exclusion masks
    
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4> restoreType(CellGridWriter<void, void, void, void> *rasa) {
  return CellGridWriter<T, Tacc,
                        Tcalc, T4>(rasa->theme, rasa->system_count, rasa->total_cell_count,
                                   rasa->total_chain_count, rasa->mesh_ticks,
                                   rasa->cell_base_capacity, rasa->cell_excl_capacity,
                                   rasa->lpos_scale, rasa->lpos_inv_scale, rasa->system_cell_grids,
                                   reinterpret_cast<Tcalc*>(rasa->system_cell_umat),
                                   reinterpret_cast<T*>(rasa->system_cell_invu),
                                   reinterpret_cast<T*>(rasa->system_pmig_invu),
                                   rasa->cell_limits, rasa->cell_limits_old, rasa->chain_limits,
                                   rasa->system_chain_bounds, rasa->chain_system_owner,
                                   reinterpret_cast<T4*>(rasa->image),
                                   reinterpret_cast<const T4*>(rasa->image_old),
                                   rasa->entry_room, rasa->exit_room, rasa->transit_room_bounds,
                                   rasa->entry_room_counts, rasa->exit_room_counts,
                                   rasa->nonimg_atom_idx, rasa->img_atom_idx,
                                   rasa->exclusion_maps, rasa->exclusion_maps_old,
                                   rasa->nt_groups, reinterpret_cast<Tacc*>(rasa->xfrc),
                                   reinterpret_cast<Tacc*>(rasa->yfrc),
                                   reinterpret_cast<Tacc*>(rasa->zfrc), rasa->xfrc_ovrf,
                                   rasa->yfrc_ovrf, rasa->zfrc_ovrf,
                                   reinterpret_cast<Tacc*>(rasa->xfrc_hw),
                                   reinterpret_cast<Tacc*>(rasa->yfrc_hw),
                                   reinterpret_cast<Tacc*>(rasa->zfrc_hw), rasa->xfrc_hw_ovrf,
                                   rasa->yfrc_hw_ovrf, rasa->zfrc_hw_ovrf);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4> restoreType(CellGridWriter<void, void, void, void> &rasa) {
  return CellGridWriter<T, Tacc,
                        Tcalc, T4>(rasa.theme, rasa.system_count, rasa.total_cell_count,
                                   rasa.total_chain_count, rasa.mesh_ticks,
                                   rasa.cell_base_capacity, rasa.cell_excl_capacity,
                                   rasa.lpos_scale, rasa.lpos_inv_scale, rasa.system_cell_grids,
                                   reinterpret_cast<Tcalc*>(rasa.system_cell_umat),
                                   reinterpret_cast<T*>(rasa.system_cell_invu),
                                   reinterpret_cast<T*>(rasa.system_pmig_invu),
                                   rasa.cell_limits, rasa.cell_limits_old, rasa.chain_limits,
                                   rasa.system_chain_bounds, rasa.chain_system_owner,
                                   reinterpret_cast<T4*>(rasa.image),
                                   reinterpret_cast<const T4*>(rasa.image_old),
                                   rasa.entry_room, rasa.exit_room, rasa.transit_room_bounds,
                                   rasa.entry_room_counts, rasa.exit_room_counts,
                                   rasa.nonimg_atom_idx, rasa.img_atom_idx, rasa.exclusion_maps,
                                   rasa.exclusion_maps_old, rasa.nt_groups,
                                   reinterpret_cast<Tacc*>(rasa.xfrc),
                                   reinterpret_cast<Tacc*>(rasa.yfrc),
                                   reinterpret_cast<Tacc*>(rasa.zfrc), rasa.xfrc_ovrf,
                                   rasa.yfrc_ovrf, rasa.zfrc_ovrf,
                                   reinterpret_cast<Tacc*>(rasa.xfrc_hw),
                                   reinterpret_cast<Tacc*>(rasa.yfrc_hw),
                                   reinterpret_cast<Tacc*>(rasa.zfrc_hw), rasa.xfrc_hw_ovrf,
                                   rasa.yfrc_hw_ovrf, rasa.zfrc_hw_ovrf);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
restoreType(const CellGridReader<void, void, void, void> *rasa) {
  return CellGridReader<T, Tacc,
                        Tcalc, T4>(rasa->theme, rasa->system_count, rasa->total_cell_count,
                                   rasa->total_chain_count, rasa->mesh_ticks,
                                   rasa->cell_base_capacity, rasa->cell_excl_capacity,
                                   rasa->lpos_scale, rasa->lpos_inv_scale, rasa->system_cell_grids,
                                   reinterpret_cast<const Tcalc*>(rasa->system_cell_umat),
                                   reinterpret_cast<const T*>(rasa->system_cell_invu),
                                   reinterpret_cast<const T*>(rasa->system_pmig_invu),
                                   rasa->cell_limits, rasa->chain_limits,
                                   rasa->system_chain_bounds, rasa->chain_system_owner,
                                   reinterpret_cast<const T4*>(rasa->image),
                                   rasa->nonimg_atom_idx, rasa->img_atom_idx,
                                   rasa->exclusion_maps, rasa->nt_groups,
                                   reinterpret_cast<const Tacc*>(rasa->xfrc),
                                   reinterpret_cast<const Tacc*>(rasa->yfrc),
                                   reinterpret_cast<const Tacc*>(rasa->zfrc), rasa->xfrc_ovrf,
                                   rasa->yfrc_ovrf, rasa->zfrc_ovrf);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
restoreType(const CellGridReader<void, void, void, void> &rasa) {
  return CellGridReader<T, Tacc,
                        Tcalc, T4>(rasa.theme, rasa.system_count, rasa.total_cell_count,
                                   rasa.total_chain_count, rasa.mesh_ticks,
                                   rasa.cell_base_capacity, rasa.cell_excl_capacity,
                                   rasa.lpos_scale, rasa.lpos_inv_scale, rasa.system_cell_grids,
                                   reinterpret_cast<const Tcalc*>(rasa.system_cell_umat),
                                   reinterpret_cast<const T*>(rasa.system_cell_invu),
                                   reinterpret_cast<const T*>(rasa.system_pmig_invu),
                                   rasa.cell_limits, rasa.chain_limits, rasa.system_chain_bounds,
                                   rasa.chain_system_owner,
                                   reinterpret_cast<const T4*>(rasa.image), rasa.nonimg_atom_idx,
                                   rasa.img_atom_idx, rasa.exclusion_maps, rasa.nt_groups,
                                   reinterpret_cast<const Tacc*>(rasa.xfrc),
                                   reinterpret_cast<const Tacc*>(rasa.yfrc),
                                   reinterpret_cast<const Tacc*>(rasa.zfrc), rasa.xfrc_ovrf,
                                   rasa.yfrc_ovrf, rasa.zfrc_ovrf);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int3 locateDecompositionCell(T *x, T *y, T *z, const double* umat, const double* invu,
                             const double cell_na, const double cell_nb, const double cell_nc,
                             const UnitCellType unit_cell, const double scale) {
  int3 result;
  double frac_x, frac_y, frac_z;
  switch (unit_cell) {
  case UnitCellType::NONE:
    rtErr("Periodic boundary conditions are essential for cell grid placement.",
          "locateDecompositionCell");
  case UnitCellType::ORTHORHOMBIC:
    if (scale > 1.01) {
      frac_x = umat[0] * (static_cast<double>(*x) / scale);
      frac_y = umat[4] * (static_cast<double>(*y) / scale);
      frac_z = umat[8] * (static_cast<double>(*z) / scale);
    }
    else {
      frac_x = umat[0] * (*x);
      frac_y = umat[4] * (*y);
      frac_z = umat[8] * (*z);
    }
    break;
  case UnitCellType::TRICLINIC:
    if (scale > 1.01) {
      const double dx = static_cast<double>(*x) / scale;
      const double dy = static_cast<double>(*y) / scale;
      const double dz = static_cast<double>(*z) / scale;
      frac_x = (umat[0] * dx) + (umat[3] * dy) + (umat[6] * dz);
      frac_y =                  (umat[4] * dy) + (umat[7] * dz);
      frac_z =                                   (umat[8] * dz);
    }
    else {
      frac_x = (umat[0] * (*x)) + (umat[3] * (*y)) + (umat[6] * (*z));
      frac_y =                    (umat[4] * (*y)) + (umat[7] * (*z));
      frac_z =                                       (umat[8] * (*z));
    }
    break;
  }
  frac_x -= floor(frac_x);
  frac_y -= floor(frac_y);
  frac_z -= floor(frac_z);
  result.x = floor(frac_x * cell_na);
  result.y = floor(frac_y * cell_nb);
  result.z = floor(frac_z * cell_nc);
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    if (scale > 1.01) {
      *x = llround(invu[0] * frac_x * scale);
      *y = llround(invu[4] * frac_y * scale);
      *z = llround(invu[8] * frac_z * scale);
    }
    else {
      *x = invu[0] * frac_x;
      *y = invu[4] * frac_y;
      *z = invu[8] * frac_z;
    }
    break;
  case UnitCellType::TRICLINIC:
    if (scale > 1.01) {
      *x = llround(((invu[0] * frac_x) + (invu[3] * frac_y) + (invu[6] * frac_z)) * scale);
      *y =                      llround(((invu[4] * frac_y) + (invu[7] * frac_z)) * scale);
      *z =                                             llround(invu[8] * frac_z * scale);
    }
    else {
      *x = (invu[0] * frac_x) + (invu[3] * frac_y) + (invu[6] * frac_z);
      *y =                      (invu[4] * frac_y) + (invu[7] * frac_z);
      *z =                                           (invu[8] * frac_z);
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void locateDecompositionCell(T* x, T* y, T* z, int* cell_idx, const int natom, const double* umat,
                             const double* invu, const int cell_na, const int cell_nb,
                             const int cell_nc, const UnitCellType unit_cell, const double scale) {
  const double dcell_na = cell_na;
  const double dcell_nb = cell_nb;
  const double dcell_nc = cell_nc;
  const bool fixed_precision = (scale > 1.01);
  switch (unit_cell) {
  case UnitCellType::NONE:
    rtErr("Periodic boundary conditions are essential for cell grid placement.",
          "locateDecompositionCell");
  case UnitCellType::ORTHORHOMBIC:
    for (int i = 0; i < natom; i++) {
      double frac_x, frac_y, frac_z;
      if (fixed_precision) {
        frac_x = umat[0] * (static_cast<double>(x[i]) / scale);
        frac_y = umat[4] * (static_cast<double>(y[i]) / scale);
        frac_z = umat[8] * (static_cast<double>(z[i]) / scale);
      }
      else {
        frac_x = umat[0] * (x[i]);
        frac_y = umat[4] * (y[i]);
        frac_z = umat[8] * (z[i]);
      }
      frac_x -= floor(frac_x);
      frac_y -= floor(frac_y);
      frac_z -= floor(frac_z);
      const int tc_a = floor(frac_x * dcell_na);
      const int tc_b = floor(frac_y * dcell_nb);
      const int tc_c = floor(frac_z * dcell_nc);
      cell_idx[i] = (((tc_c * cell_nb) + tc_b) * cell_na) + tc_a;
      if (fixed_precision) {
        x[i] = llround(invu[0] * frac_x * scale);
        y[i] = llround(invu[4] * frac_y * scale);
        z[i] = llround(invu[8] * frac_z * scale);
      }
      else {
        x[i] = invu[0] * frac_x;
        y[i] = invu[4] * frac_y;
        z[i] = invu[8] * frac_z;
      }
    }
    break;
  case UnitCellType::TRICLINIC:
    for (int i = 0; i < natom; i++) {
      double frac_x, frac_y, frac_z;
      if (fixed_precision) {
        const double dx = static_cast<double>(x[i]) / scale;
        const double dy = static_cast<double>(y[i]) / scale;
        const double dz = static_cast<double>(z[i]) / scale;
        frac_x = (umat[0] * dx) + (umat[3] * dy) + (umat[6] * dz);
        frac_y =                  (umat[4] * dy) + (umat[7] * dz);
        frac_z =                                   (umat[8] * dz);
      }
      else {
        frac_x = (umat[0] * x[i]) + (umat[3] * y[i]) + (umat[6] * z[i]);
        frac_y =                    (umat[4] * y[i]) + (umat[7] * z[i]);
        frac_z =                                       (umat[8] * z[i]);
      }
      frac_x -= floor(frac_x);
      frac_y -= floor(frac_y);
      frac_z -= floor(frac_z);
      const int tc_a = floor(frac_x * dcell_na);
      const int tc_b = floor(frac_y * dcell_nb);
      const int tc_c = floor(frac_z * dcell_nc);
      cell_idx[i] = (((tc_c * cell_nb) + tc_b) * cell_na) + tc_a;
      if (fixed_precision) {
        x[i] = llround(((invu[0] * frac_x) + (invu[3] * frac_y) + (invu[6] * frac_z)) * scale);
        y[i] =                      llround(((invu[4] * frac_y) + (invu[7] * frac_z)) * scale);
        z[i] =                                             llround(invu[8] * frac_z * scale);
      }
      else {
        x[i] = (invu[0] * frac_x) + (invu[3] * frac_y) + (invu[6] * frac_z);
        y[i] =                      (invu[4] * frac_y) + (invu[7] * frac_z);
        z[i] =                                           (invu[8] * frac_z);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void locateDecompositionCell(std::vector<T> *x, std::vector<T> *y, std::vector<T> *z,
                             std::vector<int> *cell_idx, const double* umat, const double* invu,
                             const double cell_na, const double cell_nb, const double cell_nc,
                             const UnitCellType unit_cell) {
  const size_t natom = x->size();
  if (y->size() != natom || z->size() != natom || cell_idx->size() != natom) {
    rtErr("Consistent numbers of Cartesian X, Y, and Z coordinates, as well as space to store the "
          "result cell indices, are required.", "locateDecompositionCell");
  }
  locateDecompositionCell(x->data(), y->data(), z->data(), cell_idx->data(), natom, umat, invu,
                          cell_na, cell_nb, cell_nc, unit_cell);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void locateDecompositionCell(Hybrid<T> *x, Hybrid<T> *y, Hybrid<T> *z, Hybrid<int> *cell_idx,
                             const double* umat, const double* invu, const double cell_na,
                             const double cell_nb, const double cell_nc,
                             const UnitCellType unit_cell) {
  const size_t natom = x->size();
  if (y->size() != natom || z->size() != natom || cell_idx->size() != natom) {
    rtErr("Consistent numbers of Cartesian X, Y, and Z coordinates, as well as space to store the "
          "result cell indices, are required.", "locateDecompositionCell");
  }
  locateDecompositionCell(x->data(), y->data(), z->data(), cell_idx->data(), natom, umat, invu,
                          cell_na, cell_nb, cell_nc, unit_cell);
}

} // Namespace energy
} // namespace stormm
