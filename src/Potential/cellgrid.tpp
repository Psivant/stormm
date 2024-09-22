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
                                                   const uint cell_base_capacity_in,
                                                   const float lpos_scale_in,
                                                   const float inv_lpos_scale_in,
                                                   const float frc_scale_in,
                                                   const ullint* system_cell_grids_in,
                                                   Tcalc* system_cell_umat_in,
                                                   T* system_cell_invu_in, T* system_pmig_invu_in,
                                                   uint2* cell_limits_in,
                                                   uint2* cell_limits_alt_in,
                                                   const uint* chain_limits_in,
                                                   const int* system_chain_bounds_in,
                                                   const int* chain_system_owner_in, T4* image_in,
                                                   T4* image_alt_in, uchar* migration_keys_in,
                                                   int* wander_count_in, int* wander_count_alt_in,
                                                   uint2* wanderers_in, int* nonimg_atom_idx_in,
                                                   int* nonimg_atom_idx_alt_in,
                                                   uint* img_atom_idx_in,
                                                   uint* img_atom_idx_alt_in,
                                                   ushort* img_atom_chn_cell_in,
                                                   ushort* img_atom_chn_cell_alt_in,
                                                   const int* nt_groups_in, Tacc* xfrc_in,
                                                   Tacc* yfrc_in, Tacc* zfrc_in, int* xfrc_ovrf_in,
                                                   int* yfrc_ovrf_in, int* zfrc_ovrf_in) :
    theme{theme_in}, system_count{system_count_in}, total_cell_count{total_cell_count_in},
    twice_cell_count{total_cell_count_in * 2}, total_chain_count{total_chain_count_in},
    mesh_ticks{mesh_ticks_in}, cell_base_capacity{cell_base_capacity_in},
    lpos_scale{lpos_scale_in}, inv_lpos_scale{inv_lpos_scale_in}, frc_scale{frc_scale_in},
    system_cell_grids{system_cell_grids_in}, system_cell_umat{system_cell_umat_in},
    system_cell_invu{system_cell_invu_in}, system_pmig_invu{system_pmig_invu_in},
    cell_limits{cell_limits_in}, cell_limits_alt{cell_limits_alt_in},
    chain_limits{chain_limits_in}, system_chain_bounds{system_chain_bounds_in},
    chain_system_owner{chain_system_owner_in}, image{image_in}, image_alt{image_alt_in},
    migration_keys{migration_keys_in}, wander_count{wander_count_in},
    wander_count_alt{wander_count_alt_in}, wanderers{wanderers_in},
    nonimg_atom_idx{nonimg_atom_idx_in}, nonimg_atom_idx_alt{nonimg_atom_idx_alt_in},
    img_atom_idx{img_atom_idx_in}, img_atom_idx_alt{img_atom_idx_alt_in},
    img_atom_chn_cell{img_atom_chn_cell_in}, img_atom_chn_cell_alt{img_atom_chn_cell_alt_in},
    nt_groups{nt_groups_in}, xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in}, xfrc_ovrf{xfrc_ovrf_in},
    yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridReader<T, Tacc, Tcalc, T4>::CellGridReader(const NonbondedTheme theme_in,
                                                   const int system_count_in,
                                                   const int total_cell_count_in,
                                                   const int total_chain_count_in,
                                                   const int mesh_ticks_in,
                                                   const uint cell_base_capacity_in,
                                                   const float lpos_scale_in,
                                                   const float inv_lpos_scale_in,
                                                   const float inv_frc_scale_in,
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
                                                   const ushort* img_atom_chn_cell_in,
                                                   const int* nt_groups_in,
                                                   const Tacc* xfrc_in, const Tacc* yfrc_in,
                                                   const Tacc* zfrc_in, const int* xfrc_ovrf_in,
                                                   const int* yfrc_ovrf_in,
                                                   const int* zfrc_ovrf_in) :
    theme{theme_in}, system_count{system_count_in}, total_cell_count{total_cell_count_in},
    twice_cell_count{total_cell_count_in * 2}, total_chain_count{total_chain_count_in},
    mesh_ticks{mesh_ticks_in}, cell_base_capacity{cell_base_capacity_in},
    lpos_scale{lpos_scale_in}, inv_lpos_scale{inv_lpos_scale_in}, inv_frc_scale{inv_frc_scale_in},
    system_cell_grids{system_cell_grids_in}, system_cell_umat{system_cell_umat_in},
    system_cell_invu{system_cell_invu_in}, system_pmig_invu{system_pmig_invu_in},
    cell_limits{cell_limits_in}, chain_limits{chain_limits_in},
    system_chain_bounds{system_chain_bounds_in}, chain_system_owner{chain_system_owner_in},
    image{image_in}, nonimg_atom_idx{nonimg_atom_idx_in}, img_atom_idx{img_atom_idx_in},
    img_atom_chn_cell{img_atom_chn_cell_in}, nt_groups{nt_groups_in}, xfrc{xfrc_in},
    yfrc{yfrc_in}, zfrc{zfrc_in}, xfrc_ovrf{xfrc_ovrf_in}, yfrc_ovrf{yfrc_ovrf_in},
    zfrc_ovrf{zfrc_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridReader<T, Tacc, Tcalc, T4>::CellGridReader(const CellGridWriter<T, Tacc, Tcalc, T4> &cgw) :
    theme{cgw.theme}, system_count{cgw.system_count}, total_cell_count{cgw.total_cell_count},
    twice_cell_count{cgw.twice_cell_count}, total_chain_count{cgw.total_chain_count},
    mesh_ticks{cgw.mesh_ticks}, cell_base_capacity{cgw.cell_base_capacity},
    lpos_scale{cgw.lpos_scale}, inv_lpos_scale{cgw.inv_lpos_scale},
    inv_frc_scale{static_cast<float>(1.0 / cgw.frc_scale)},
    system_cell_grids{cgw.system_cell_grids}, system_cell_umat{cgw.system_cell_umat},
    system_cell_invu{cgw.system_cell_invu}, system_pmig_invu{cgw.system_pmig_invu},
    cell_limits{cgw.cell_limits}, chain_limits{cgw.chain_limits},
    system_chain_bounds{cgw.system_chain_bounds}, chain_system_owner{cgw.chain_system_owner},
    image{cgw.image}, nonimg_atom_idx{cgw.nonimg_atom_idx}, img_atom_idx{cgw.img_atom_idx},
    img_atom_chn_cell{cgw.img_atom_chn_cell}, nt_groups{cgw.nt_groups}, xfrc{cgw.xfrc},
    yfrc{cgw.yfrc}, zfrc{cgw.zfrc}, xfrc_ovrf{cgw.xfrc_ovrf}, yfrc_ovrf{cgw.yfrc_ovrf},
    zfrc_ovrf{cgw.zfrc_ovrf}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridReader<T, Tacc, Tcalc, T4>::CellGridReader(const CellGridWriter<T, Tacc, Tcalc, T4> *cgw) :
    theme{cgw->theme}, system_count{cgw->system_count}, total_cell_count{cgw->total_cell_count},
    twice_cell_count{cgw->twice_cell_count}, total_chain_count{cgw->total_chain_count},
    mesh_ticks{cgw->mesh_ticks}, cell_base_capacity{cgw->cell_base_capacity},
    lpos_scale{cgw->lpos_scale}, inv_lpos_scale{cgw->inv_lpos_scale},
    inv_frc_scale{static_cast<float>(1.0 / cgw->frc_scale)},
    system_cell_grids{cgw->system_cell_grids}, system_cell_umat{cgw->system_cell_umat},
    system_cell_invu{cgw->system_cell_invu}, system_pmig_invu{cgw->system_pmig_invu},
    cell_limits{cgw->cell_limits}, chain_limits{cgw->chain_limits},
    system_chain_bounds{cgw->system_chain_bounds}, chain_system_owner{cgw->chain_system_owner},
    image{cgw->image}, nonimg_atom_idx{cgw->nonimg_atom_idx}, img_atom_idx{cgw->img_atom_idx},
    img_atom_chn_cell{cgw->img_atom_chn_cell}, nt_groups{cgw->nt_groups},
    xfrc{cgw->xfrc}, yfrc{cgw->yfrc}, zfrc{cgw->zfrc}, xfrc_ovrf{cgw->xfrc_ovrf},
    yfrc_ovrf{cgw->yfrc_ovrf}, zfrc_ovrf{cgw->zfrc_ovrf}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(const PhaseSpaceSynthesis *poly_ps_ptr_in,
                                       const AtomGraphSynthesis *poly_ag_ptr_in,
                                       const double cutoff, const double padding,
                                       const int mesh_subdivisions_in,
                                       const NonbondedTheme theme_in,
                                       const uint cell_base_capacity_in,
                                       const ExceptionResponse policy_in) :
    system_count{0}, total_cell_count{0}, total_chain_count{0},
    cell_base_capacity{cell_base_capacity_in},
    effective_cutoff{validateEffectiveCutoff(cutoff + padding, policy_in)},
    mesh_subdivisions{mesh_subdivisions_in},
    policy{policy_in},
    theme{theme_in},
    has_tiny_box{TinyBoxPresence::NO},
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
    nonimaged_atom_indices{HybridKind::ARRAY, "cg_atom_idx"},
    nonimaged_atom_indices_alt{HybridKind::ARRAY, "cg_atom_idx_alt"},
    image_array_indices{HybridKind::ARRAY, "cg_img_idx"},
    image_array_indices_alt{HybridKind::ARRAY, "cg_img_idx_alt"},
    image_chain_cell_indices{HybridKind::ARRAY, "cg_img_chn_cell"},
    image_chain_cell_indices_alt{HybridKind::ARRAY, "cg_img_chn_cell_alt"},
    cell_migrations{HybridKind::ARRAY, "cg_migrations"},
    wandering_atom_count{HybridKind::ARRAY, "cg_wander_count"},
    wandering_atom_count_alt{HybridKind::ARRAY, "cg_wander_count_alt"},
    wanderers{HybridKind::ARRAY, "cg_wanderers"},
    nt_work_groups{HybridKind::ARRAY, "cg_nt_work_groups"},
    x_force{HybridKind::ARRAY, "cg_xfrc"},
    y_force{HybridKind::ARRAY, "cg_yfrc"},
    z_force{HybridKind::ARRAY, "cg_zfrc"},
    x_force_overflow{HybridKind::ARRAY, "cg_xfrc_ovrf"},
    y_force_overflow{HybridKind::ARRAY, "cg_yfrc_ovrf"},
    z_force_overflow{HybridKind::ARRAY, "cg_zfrc_ovrf"},
    origin_offset_stride{warp_size_int},
    cell_origins_ax{HybridKind::POINTER, "cg_ax_orig"},
    cell_origins_bx{HybridKind::POINTER, "cg_bx_orig"},
    cell_origins_by{HybridKind::POINTER, "cg_by_orig"},
    cell_origins_cx{HybridKind::POINTER, "cg_cx_orig"},
    cell_origins_cy{HybridKind::POINTER, "cg_cy_orig"},
    cell_origins_cz{HybridKind::POINTER, "cg_cy_orig"},
    cell_origins_ax_overflow{HybridKind::POINTER, "cg_ax_orig_ovrf"},
    cell_origins_bx_overflow{HybridKind::POINTER, "cg_bx_orig_ovrf"},
    cell_origins_by_overflow{HybridKind::POINTER, "cg_by_orig_ovrf"},
    cell_origins_cx_overflow{HybridKind::POINTER, "cg_cx_orig_ovrf"},
    cell_origins_cy_overflow{HybridKind::POINTER, "cg_cy_orig_ovrf"},
    cell_origins_cz_overflow{HybridKind::POINTER, "cg_cy_orig_ovrf"},
    alt_cell_origins_ax{HybridKind::POINTER, "cg_ax_alt_orig"},
    alt_cell_origins_bx{HybridKind::POINTER, "cg_bx_alt_orig"},
    alt_cell_origins_by{HybridKind::POINTER, "cg_by_alt_orig"},
    alt_cell_origins_cx{HybridKind::POINTER, "cg_cx_alt_orig"},
    alt_cell_origins_cy{HybridKind::POINTER, "cg_cy_alt_orig"},
    alt_cell_origins_cz{HybridKind::POINTER, "cg_cy_alt_orig"},
    alt_cell_origins_ax_overflow{HybridKind::POINTER, "cg_ax_alt_orig_ovrf"},
    alt_cell_origins_bx_overflow{HybridKind::POINTER, "cg_bx_alt_orig_ovrf"},
    alt_cell_origins_by_overflow{HybridKind::POINTER, "cg_by_alt_orig_ovrf"},
    alt_cell_origins_cx_overflow{HybridKind::POINTER, "cg_cx_alt_orig_ovrf"},
    alt_cell_origins_cy_overflow{HybridKind::POINTER, "cg_cy_alt_orig_ovrf"},
    alt_cell_origins_cz_overflow{HybridKind::POINTER, "cg_cy_alt_orig_ovrf"},
    origin_llint_data{HybridKind::ARRAY, "cg_orig_llint"},
    origin_int_data{HybridKind::ARRAY, "cg_orig_int"},
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
    if (cell_na[i] > maximum_cellgrid_span || cell_nb[i] > maximum_cellgrid_span ||
        cell_nc[i] > maximum_cellgrid_span) {
      rtErr("The system spatial decomposition would require a grid of [ " +
            std::to_string(cell_na[i]) + " x " + std::to_string(cell_nb[i]) + " x " +
            std::to_string(cell_nc[i]) + " ] cells.  A maximum of " +
            std::to_string(maximum_cellgrid_span) + " is enforced along all axes.", "CellGrid");
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

  // Check for the existence of a tiny box, or one with a short side which would require special
  // treatment of interparticle displacement calculations.
  for (int i = 0; i < poly_psr.system_count; i++) {
    if (cell_na[i] < 5 || cell_nb[i] < 5 || cell_nc[i] < 5) {
      has_tiny_box = TinyBoxPresence::YES;
    }
  }

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
  const int minimum_chain_length = minValue(cell_na);
  for (int pos = 0; pos < poly_psr.system_count; pos++) {
    tallyCellPopulations(&cell_populations, pos, cell_na[pos], cell_nb[pos], cell_nc[pos]);
    const int max_pop = maxValue(cell_populations);
    uint trial_base_capacity;
    if (minimum_chain_length < 8) {
      trial_base_capacity = roundUp<uint>(max_pop * 3 / 2, 64);
    }
    else if (minimum_chain_length < 12) {
      trial_base_capacity = roundUp<uint>(max_pop * 3 / 2, 32);
    }
    else if (minimum_chain_length < 16) {
      trial_base_capacity = roundUp<uint>(max_pop * 3 / 2, 16);
    }
    else if (minimum_chain_length < 24) {
      trial_base_capacity = roundUp<uint>(max_pop * 5 / 4, 8);
    }
    else {
      trial_base_capacity = roundUp<uint>(max_pop, 8);
    }
    cell_base_capacity = std::max(trial_base_capacity, cell_base_capacity);
  }
  
  // Allocate the bulk of the space that the object will use in partitioning atoms.
  const uint image_size = static_cast<uint>(total_cell_count) * cell_base_capacity;
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
    const uint max_atoms_per_chain = static_cast<uint>(cell_na[pos]) * cell_base_capacity;
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
  image_chain_cell_indices.resize(icl_atom_counter);
  image_chain_cell_indices_alt.resize(icl_atom_counter);  

  // Allocate space for moving atoms between chains and cells within each chain
  cell_migrations.resize(image_size);
  wandering_atom_count.resize(warp_size_int);
  wandering_atom_count_alt.resize(warp_size_int);
  wanderers.resize(maximum_wandering_atoms);
  
  // Loop back over all structures, reimage atoms, and pack the image arrays.
  populateImage(CoordinateCycle::WHITE);
  populateImage(CoordinateCycle::BLACK);

  // Arrange the fixed-precision cell origins
  drawNeighborListRulers();  
  
  // Create work units to support the "tower and plate" neutral territory decomposition
  prepareWorkGroups();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(const PhaseSpaceSynthesis *poly_ps_in,
                                       const AtomGraphSynthesis &poly_ag_in, const double cutoff,
                                       const double padding, const int mesh_subdivisions_in,
                                       const NonbondedTheme theme_in, const uint base_capacity_in,
                                       const ExceptionResponse policy_in) :
    CellGrid(poly_ps_in, poly_ag_in.getSelfPointer(), cutoff, padding, mesh_subdivisions_in,
             theme_in, base_capacity_in, policy_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(const PhaseSpaceSynthesis &poly_ps_in,
                                       const AtomGraphSynthesis &poly_ag_in, const double cutoff,
                                       const double padding, const int mesh_subdivisions_in,
                                       const NonbondedTheme theme_in, const uint base_capacity_in,
                                       const ExceptionResponse policy_in) :
    CellGrid(poly_ps_in.getSelfPointer(), poly_ag_in.getSelfPointer(), cutoff, padding,
             mesh_subdivisions_in, theme_in, base_capacity_in, policy_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(const CellGrid &original) :
    system_count{original.system_count},
    total_cell_count{original.total_cell_count},
    total_chain_count{original.total_chain_count},
    cell_base_capacity{original.cell_base_capacity},
    effective_cutoff{original.effective_cutoff},
    mesh_subdivisions{original.mesh_subdivisions},
    policy{original.policy},
    theme{original.theme},
    has_tiny_box{original.has_tiny_box},
    cycle_position{original.cycle_position},
    localpos_scale_bits{original.localpos_scale_bits},
    localpos_scale{original.localpos_scale},
    localpos_inverse_scale{original.localpos_inverse_scale},
    system_cell_grids{original.system_cell_grids},
    system_cell_umat{original.system_cell_umat},
    system_cell_umat_alt{original.system_cell_umat_alt},
    system_cell_invu{original.system_cell_invu},
    system_cell_invu_alt{original.system_cell_invu_alt},
    system_pmig_invu{original.system_pmig_invu},
    system_pmig_invu_alt{original.system_pmig_invu_alt},
    particle_mesh_grid_invu{original.particle_mesh_grid_invu},
    particle_mesh_grid_invu_alt{original.particle_mesh_grid_invu_alt},
    image_cell_limits{original.image_cell_limits},
    image_cell_limits_alt{original.image_cell_limits_alt},
    image_chain_limits{original.image_chain_limits},
    system_chain_bounds{original.system_chain_bounds},
    chain_system_membership{original.chain_system_membership},
    image{original.image},
    image_alt{original.image_alt},
    nonimaged_atom_indices{original.nonimaged_atom_indices},
    nonimaged_atom_indices_alt{original.nonimaged_atom_indices_alt},
    image_array_indices{original.image_array_indices},
    image_array_indices_alt{original.image_array_indices_alt},
    image_chain_cell_indices{original.image_chain_cell_indices},
    image_chain_cell_indices_alt{original.image_chain_cell_indices_alt},
    cell_migrations{original.cell_migrations},
    wandering_atom_count{original.wandering_atom_count},
    wandering_atom_count_alt{original.wandering_atom_count_alt},
    wanderers{original.wanderers},
    nt_work_groups{original.nt_work_groups},
    x_force{original.x_force},
    y_force{original.y_force},
    z_force{original.z_force},
    x_force_overflow{original.x_force_overflow},
    y_force_overflow{original.y_force_overflow},
    z_force_overflow{original.z_force_overflow},
    origin_offset_stride{original.origin_offset_stride},
    cell_origins_ax{original.cell_origins_ax},
    cell_origins_bx{original.cell_origins_bx},
    cell_origins_by{original.cell_origins_by},
    cell_origins_cx{original.cell_origins_cx},
    cell_origins_cy{original.cell_origins_cy},
    cell_origins_cz{original.cell_origins_cz},
    cell_origins_ax_overflow{original.cell_origins_ax_overflow},
    cell_origins_bx_overflow{original.cell_origins_bx_overflow},
    cell_origins_by_overflow{original.cell_origins_by_overflow},
    cell_origins_cx_overflow{original.cell_origins_cx_overflow},
    cell_origins_cy_overflow{original.cell_origins_cy_overflow},
    cell_origins_cz_overflow{original.cell_origins_cz_overflow},
    alt_cell_origins_ax{original.alt_cell_origins_ax},
    alt_cell_origins_bx{original.alt_cell_origins_bx},
    alt_cell_origins_by{original.alt_cell_origins_by},
    alt_cell_origins_cx{original.alt_cell_origins_cx},
    alt_cell_origins_cy{original.alt_cell_origins_cy},
    alt_cell_origins_cz{original.alt_cell_origins_cz},
    alt_cell_origins_ax_overflow{original.alt_cell_origins_ax_overflow},
    alt_cell_origins_bx_overflow{original.alt_cell_origins_bx_overflow},
    alt_cell_origins_by_overflow{original.alt_cell_origins_by_overflow},
    alt_cell_origins_cx_overflow{original.alt_cell_origins_cx_overflow},
    alt_cell_origins_cy_overflow{original.alt_cell_origins_cy_overflow},
    alt_cell_origins_cz_overflow{original.alt_cell_origins_cz_overflow},
    origin_llint_data{original.origin_llint_data},
    origin_int_data{original.origin_int_data},
    poly_ps_ptr{original.poly_ps_ptr},
    poly_ag_ptr{original.poly_ag_ptr}
{
  // Repair the pointers in the cell origins rulers
  rebaseRulers();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>::CellGrid(CellGrid &&original) :
    system_count{original.system_count},
    total_cell_count{original.total_cell_count},
    total_chain_count{original.total_chain_count},
    cell_base_capacity{original.cell_base_capacity},
    effective_cutoff{original.effective_cutoff},
    mesh_subdivisions{original.mesh_subdivisions},
    policy{original.policy},
    theme{original.theme},
    has_tiny_box{original.has_tiny_box},
    cycle_position{original.cycle_position},
    localpos_scale_bits{original.localpos_scale_bits},
    localpos_scale{original.localpos_scale},
    localpos_inverse_scale{original.localpos_inverse_scale},
    system_cell_grids{std::move(original.system_cell_grids)},
    system_cell_umat{std::move(original.system_cell_umat)},
    system_cell_umat_alt{std::move(original.system_cell_umat_alt)},
    system_cell_invu{std::move(original.system_cell_invu)},
    system_cell_invu_alt{std::move(original.system_cell_invu_alt)},
    system_pmig_invu{std::move(original.system_pmig_invu)},
    system_pmig_invu_alt{std::move(original.system_pmig_invu_alt)},
    particle_mesh_grid_invu{std::move(original.particle_mesh_grid_invu)},
    particle_mesh_grid_invu_alt{std::move(original.particle_mesh_grid_invu_alt)},
    image_cell_limits{std::move(original.image_cell_limits)},
    image_cell_limits_alt{std::move(original.image_cell_limits_alt)},
    image_chain_limits{std::move(original.image_chain_limits)},
    system_chain_bounds{std::move(original.system_chain_bounds)},
    chain_system_membership{std::move(original.chain_system_membership)},
    image{std::move(original.image)},
    image_alt{std::move(original.image_alt)},
    nonimaged_atom_indices{std::move(original.nonimaged_atom_indices)},
    nonimaged_atom_indices_alt{std::move(original.nonimaged_atom_indices_alt)},
    image_array_indices{std::move(original.image_array_indices)},
    image_array_indices_alt{std::move(original.image_array_indices_alt)},
    image_chain_cell_indices{std::move(original.image_chain_cell_indices)},
    image_chain_cell_indices_alt{std::move(original.image_chain_cell_indices_alt)},
    cell_migrations{std::move(original.cell_migrations)},
    wandering_atom_count{std::move(original.wandering_atom_count)},
    wandering_atom_count_alt{std::move(original.wandering_atom_count_alt)},
    wanderers{std::move(original.wanderers)},
    nt_work_groups{std::move(original.nt_work_groups)},
    x_force{std::move(original.x_force)},
    y_force{std::move(original.y_force)},
    z_force{std::move(original.z_force)},
    x_force_overflow{std::move(original.x_force_overflow)},
    y_force_overflow{std::move(original.y_force_overflow)},
    z_force_overflow{std::move(original.z_force_overflow)},
    origin_offset_stride{original.origin_offset_stride},
    cell_origins_ax{std::move(original.cell_origins_ax)},
    cell_origins_bx{std::move(original.cell_origins_bx)},
    cell_origins_by{std::move(original.cell_origins_by)},
    cell_origins_cx{std::move(original.cell_origins_cx)},
    cell_origins_cy{std::move(original.cell_origins_cy)},
    cell_origins_cz{std::move(original.cell_origins_cz)},
    cell_origins_ax_overflow{std::move(original.cell_origins_ax_overflow)},
    cell_origins_bx_overflow{std::move(original.cell_origins_bx_overflow)},
    cell_origins_by_overflow{std::move(original.cell_origins_by_overflow)},
    cell_origins_cx_overflow{std::move(original.cell_origins_cx_overflow)},
    cell_origins_cy_overflow{std::move(original.cell_origins_cy_overflow)},
    cell_origins_cz_overflow{std::move(original.cell_origins_cz_overflow)},
    alt_cell_origins_ax{std::move(original.alt_cell_origins_ax)},
    alt_cell_origins_bx{std::move(original.alt_cell_origins_bx)},
    alt_cell_origins_by{std::move(original.alt_cell_origins_by)},
    alt_cell_origins_cx{std::move(original.alt_cell_origins_cx)},
    alt_cell_origins_cy{std::move(original.alt_cell_origins_cy)},
    alt_cell_origins_cz{std::move(original.alt_cell_origins_cz)},
    alt_cell_origins_ax_overflow{std::move(original.alt_cell_origins_ax_overflow)},
    alt_cell_origins_bx_overflow{std::move(original.alt_cell_origins_bx_overflow)},
    alt_cell_origins_by_overflow{std::move(original.alt_cell_origins_by_overflow)},
    alt_cell_origins_cx_overflow{std::move(original.alt_cell_origins_cx_overflow)},
    alt_cell_origins_cy_overflow{std::move(original.alt_cell_origins_cy_overflow)},
    alt_cell_origins_cz_overflow{std::move(original.alt_cell_origins_cz_overflow)},
    origin_llint_data{std::move(original.origin_llint_data)},
    origin_int_data{std::move(original.origin_int_data)},
    poly_ps_ptr{original.poly_ps_ptr},
    poly_ag_ptr{original.poly_ag_ptr}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>& CellGrid<T, Tacc, Tcalc, T4>::operator=(const CellGrid &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  system_count = other.system_count;
  total_cell_count = other.total_cell_count;
  total_chain_count = other.total_chain_count;
  cell_base_capacity = other.cell_base_capacity;
  effective_cutoff = other.effective_cutoff;
  mesh_subdivisions = other.mesh_subdivisions;
  policy = other.policy;
  theme = other.theme;
  has_tiny_box = other.has_tiny_box;
  cycle_position = other.cycle_position;
  localpos_scale_bits = other.localpos_scale_bits;
  localpos_scale = other.localpos_scale;
  localpos_inverse_scale = other.localpos_inverse_scale;
  system_cell_grids = other.system_cell_grids;
  system_cell_umat = other.system_cell_umat;
  system_cell_umat_alt = other.system_cell_umat_alt;
  system_cell_invu = other.system_cell_invu;
  system_cell_invu_alt = other.system_cell_invu_alt;
  system_pmig_invu = other.system_pmig_invu;
  system_pmig_invu_alt = other.system_pmig_invu_alt;
  particle_mesh_grid_invu = other.particle_mesh_grid_invu;
  particle_mesh_grid_invu_alt = other.particle_mesh_grid_invu_alt;
  image_cell_limits = other.image_cell_limits;
  image_cell_limits_alt = other.image_cell_limits_alt;
  image_chain_limits = other.image_chain_limits;
  system_chain_bounds = other.system_chain_bounds;
  chain_system_membership = other.chain_system_membership;
  image = other.image;
  image_alt = other.image_alt;
  nonimaged_atom_indices = other.nonimaged_atom_indices;
  nonimaged_atom_indices_alt = other.nonimaged_atom_indices_alt;
  image_array_indices = other.image_array_indices;
  image_array_indices_alt = other.image_array_indices_alt;
  image_chain_cell_indices = other.image_chain_cell_indices;
  image_chain_cell_indices_alt = other.image_chain_cell_indices_alt;
  cell_migrations = other.cell_migrations;
  wandering_atom_count = other.wandering_atom_count;
  wandering_atom_count_alt = other.wandering_atom_count_alt;
  wanderers = other.wanderers;
  nt_work_groups = other.nt_work_groups;
  x_force = other.x_force;
  y_force = other.y_force;
  z_force = other.z_force;
  x_force_overflow = other.x_force_overflow;
  y_force_overflow = other.y_force_overflow;
  z_force_overflow = other.z_force_overflow;
  origin_offset_stride = other.origin_offset_stride;
  cell_origins_ax = other.cell_origins_ax;
  cell_origins_bx = other.cell_origins_bx;
  cell_origins_by = other.cell_origins_by;
  cell_origins_cx = other.cell_origins_cx;
  cell_origins_cy = other.cell_origins_cy;
  cell_origins_cz = other.cell_origins_cz;
  cell_origins_ax_overflow = other.cell_origins_ax_overflow;
  cell_origins_bx_overflow = other.cell_origins_bx_overflow;
  cell_origins_by_overflow = other.cell_origins_by_overflow;
  cell_origins_cx_overflow = other.cell_origins_cx_overflow;
  cell_origins_cy_overflow = other.cell_origins_cy_overflow;
  cell_origins_cz_overflow = other.cell_origins_cz_overflow;
  alt_cell_origins_ax = other.alt_cell_origins_ax;
  alt_cell_origins_bx = other.alt_cell_origins_bx;
  alt_cell_origins_by = other.alt_cell_origins_by;
  alt_cell_origins_cx = other.alt_cell_origins_cx;
  alt_cell_origins_cy = other.alt_cell_origins_cy;
  alt_cell_origins_cz = other.alt_cell_origins_cz;
  alt_cell_origins_ax_overflow = other.alt_cell_origins_ax_overflow;
  alt_cell_origins_bx_overflow = other.alt_cell_origins_bx_overflow;
  alt_cell_origins_by_overflow = other.alt_cell_origins_by_overflow;
  alt_cell_origins_cx_overflow = other.alt_cell_origins_cx_overflow;
  alt_cell_origins_cy_overflow = other.alt_cell_origins_cy_overflow;
  alt_cell_origins_cz_overflow = other.alt_cell_origins_cz_overflow;
  origin_llint_data = other.origin_llint_data;
  origin_int_data = other.origin_int_data;
  poly_ps_ptr = other.poly_ps_ptr;
  poly_ag_ptr = other.poly_ag_ptr;

  // Repair the pointers in the cell origins rulers before returning the result
  rebaseRulers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGrid<T, Tacc, Tcalc, T4>& CellGrid<T, Tacc, Tcalc, T4>::operator=(CellGrid &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  system_count = other.system_count;
  total_cell_count = other.total_cell_count;
  total_chain_count = other.total_chain_count;
  cell_base_capacity = other.cell_base_capacity;
  effective_cutoff = other.effective_cutoff;
  mesh_subdivisions = other.mesh_subdivisions;
  policy = other.policy;
  theme = other.theme;
  has_tiny_box = other.has_tiny_box;
  cycle_position = other.cycle_position;
  localpos_scale_bits = other.localpos_scale_bits;
  localpos_scale = other.localpos_scale;
  localpos_inverse_scale = other.localpos_inverse_scale;
  system_cell_grids = std::move(other.system_cell_grids);
  system_cell_umat = std::move(other.system_cell_umat);
  system_cell_umat_alt = std::move(other.system_cell_umat_alt);
  system_cell_invu = std::move(other.system_cell_invu);
  system_cell_invu_alt = std::move(other.system_cell_invu_alt);
  system_pmig_invu = std::move(other.system_pmig_invu);
  system_pmig_invu_alt = std::move(other.system_pmig_invu_alt);
  particle_mesh_grid_invu = std::move(other.particle_mesh_grid_invu);
  particle_mesh_grid_invu_alt = std::move(other.particle_mesh_grid_invu_alt);
  image_cell_limits = std::move(other.image_cell_limits);
  image_cell_limits_alt = std::move(other.image_cell_limits_alt);
  image_chain_limits = std::move(other.image_chain_limits);
  system_chain_bounds = std::move(other.system_chain_bounds);
  chain_system_membership = std::move(other.chain_system_membership);
  image = std::move(other.image);
  image_alt = std::move(other.image_alt);
  nonimaged_atom_indices = std::move(other.nonimaged_atom_indices);
  nonimaged_atom_indices_alt = std::move(other.nonimaged_atom_indices_alt);
  image_array_indices = std::move(other.image_array_indices);
  image_array_indices_alt = std::move(other.image_array_indices_alt);
  image_chain_cell_indices = std::move(other.image_chain_cell_indices);
  image_chain_cell_indices_alt = std::move(other.image_chain_cell_indices_alt);
  cell_migrations = std::move(other.cell_migrations);
  wandering_atom_count = std::move(other.wandering_atom_count);
  wandering_atom_count_alt = std::move(other.wandering_atom_count_alt);
  wanderers = std::move(other.wanderers);
  nt_work_groups = std::move(other.nt_work_groups);
  x_force = std::move(other.x_force);
  y_force = std::move(other.y_force);
  z_force = std::move(other.z_force);
  x_force_overflow = std::move(other.x_force_overflow);
  y_force_overflow = std::move(other.y_force_overflow);
  z_force_overflow = std::move(other.z_force_overflow);
  origin_offset_stride = other.origin_offset_stride;
  cell_origins_ax = std::move(other.cell_origins_ax);
  cell_origins_bx = std::move(other.cell_origins_bx);
  cell_origins_by = std::move(other.cell_origins_by);
  cell_origins_cx = std::move(other.cell_origins_cx);
  cell_origins_cy = std::move(other.cell_origins_cy);
  cell_origins_cz = std::move(other.cell_origins_cz);
  cell_origins_ax_overflow = std::move(other.cell_origins_ax_overflow);
  cell_origins_bx_overflow = std::move(other.cell_origins_bx_overflow);
  cell_origins_by_overflow = std::move(other.cell_origins_by_overflow);
  cell_origins_cx_overflow = std::move(other.cell_origins_cx_overflow);
  cell_origins_cy_overflow = std::move(other.cell_origins_cy_overflow);
  cell_origins_cz_overflow = std::move(other.cell_origins_cz_overflow);
  alt_cell_origins_ax = std::move(other.alt_cell_origins_ax);
  alt_cell_origins_bx = std::move(other.alt_cell_origins_bx);
  alt_cell_origins_by = std::move(other.alt_cell_origins_by);
  alt_cell_origins_cx = std::move(other.alt_cell_origins_cx);
  alt_cell_origins_cy = std::move(other.alt_cell_origins_cy);
  alt_cell_origins_cz = std::move(other.alt_cell_origins_cz);
  alt_cell_origins_ax_overflow = std::move(other.alt_cell_origins_ax_overflow);
  alt_cell_origins_bx_overflow = std::move(other.alt_cell_origins_bx_overflow);
  alt_cell_origins_by_overflow = std::move(other.alt_cell_origins_by_overflow);
  alt_cell_origins_cx_overflow = std::move(other.alt_cell_origins_cx_overflow);
  alt_cell_origins_cy_overflow = std::move(other.alt_cell_origins_cy_overflow);
  alt_cell_origins_cz_overflow = std::move(other.alt_cell_origins_cz_overflow);
  origin_llint_data = std::move(other.origin_llint_data);
  origin_int_data = std::move(other.origin_int_data);
  poly_ps_ptr = other.poly_ps_ptr;
  poly_ag_ptr = other.poly_ag_ptr;

  // Return the result.  The std::move function obviates the need for pointer repair.
  return *this;
}

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
uint CellGrid<T, Tacc, Tcalc, T4>::getCellBaseCapacity() const {
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
TinyBoxPresence CellGrid<T, Tacc, Tcalc, T4>::getTinyBoxPresence() const {
  return has_tiny_box;
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
const PsSynthesisBorders
CellGrid<T, Tacc, Tcalc, T4>::getUnitCellTransforms(const CoordinateCycle orientation,
                                                    const HybridTargetLevel tier) const {
  return poly_ps_ptr->borders(orientation, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const PsSynthesisBorders
CellGrid<T, Tacc, Tcalc, T4>::getUnitCellTransforms(const HybridTargetLevel tier) const {
  return poly_ps_ptr->borders(cycle_position, tier);
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
                                     mesh_subdivisions, cell_base_capacity, localpos_scale,
                                     localpos_inverse_scale, poly_ps_ptr->getForceScalingFactor(),
                                     system_cell_grids.data(tier), system_cell_umat_alt.data(tier),
                                     system_cell_invu_alt.data(tier),
                                     system_pmig_invu_alt.data(tier),
                                     image_cell_limits_alt.data(tier),
                                     image_cell_limits.data(tier), image_chain_limits.data(tier),
                                     system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image_alt.data(tier),
                                     image.data(tier), cell_migrations.data(tier),
                                     wandering_atom_count_alt.data(tier),
                                     wandering_atom_count.data(tier), wanderers.data(tier),
                                     nonimaged_atom_indices_alt.data(tier),
                                     nonimaged_atom_indices.data(tier),
                                     image_array_indices_alt.data(tier),
                                     image_array_indices.data(tier),
                                     image_chain_cell_indices_alt.data(tier),
                                     image_chain_cell_indices.data(tier),
                                     nt_work_groups.data(tier), x_force.data(tier),
                                     y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier));
  case CoordinateCycle::WHITE:
    return CellGridWriter<T, Tacc,
                          Tcalc, T4>(theme, system_count, total_cell_count, total_chain_count,
                                     mesh_subdivisions, cell_base_capacity, localpos_scale,
                                     localpos_inverse_scale, poly_ps_ptr->getForceScalingFactor(),
                                     system_cell_grids.data(tier), system_cell_umat.data(tier),
                                     system_cell_invu.data(tier), system_pmig_invu.data(tier),
                                     image_cell_limits.data(tier),
                                     image_cell_limits_alt.data(tier),
                                     image_chain_limits.data(tier), system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image.data(tier),
                                     image_alt.data(tier), cell_migrations.data(tier),
                                     wandering_atom_count.data(tier),
                                     wandering_atom_count_alt.data(tier), wanderers.data(tier),
                                     nonimaged_atom_indices.data(tier),
                                     nonimaged_atom_indices_alt.data(tier),
                                     image_array_indices.data(tier),
                                     image_array_indices_alt.data(tier),
                                     image_chain_cell_indices.data(tier),
                                     image_chain_cell_indices_alt.data(tier),
                                     nt_work_groups.data(tier), x_force.data(tier),
                                     y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier));
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
                                     mesh_subdivisions, cell_base_capacity, localpos_scale,
                                     localpos_inverse_scale,
                                     poly_ps_ptr->getInverseForceScalingFactor(),
                                     system_cell_grids.data(tier), system_cell_umat_alt.data(tier),
                                     system_cell_invu_alt.data(tier),
                                     system_pmig_invu_alt.data(tier),
                                     image_cell_limits_alt.data(tier),
                                     image_chain_limits.data(tier), system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image_alt.data(tier),
                                     nonimaged_atom_indices_alt.data(tier),
                                     image_array_indices_alt.data(tier),
                                     image_chain_cell_indices_alt.data(tier),
                                     nt_work_groups.data(tier), x_force.data(tier),
                                     y_force.data(tier), z_force.data(tier),
                                     x_force_overflow.data(tier), y_force_overflow.data(tier),
                                     z_force_overflow.data(tier));
  case CoordinateCycle::WHITE:
    return CellGridReader<T, Tacc,
                          Tcalc, T4>(theme, system_count, total_cell_count, total_chain_count,
                                     mesh_subdivisions, cell_base_capacity, localpos_scale,
                                     localpos_inverse_scale,
                                     poly_ps_ptr->getInverseForceScalingFactor(),
                                     system_cell_grids.data(tier), system_cell_umat.data(tier),
                                     system_cell_invu.data(tier), system_pmig_invu.data(tier),
                                     image_cell_limits.data(tier), image_chain_limits.data(tier),
                                     system_chain_bounds.data(tier),
                                     chain_system_membership.data(tier), image.data(tier),
                                     nonimaged_atom_indices.data(tier),
                                     image_array_indices.data(tier),
                                     image_chain_cell_indices.data(tier),
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
                                      mesh_subdivisions, cell_base_capacity, localpos_scale,
                                      localpos_inverse_scale, poly_ps_ptr->getForceScalingFactor(),
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
                                      cell_migrations.data(tier),
                                      wandering_atom_count_alt.data(tier),
                                      wandering_atom_count.data(tier), wanderers.data(tier),
                                      nonimaged_atom_indices_alt.data(tier),
                                      nonimaged_atom_indices.data(tier),
                                      image_array_indices_alt.data(tier),
                                      image_array_indices.data(tier),
                                      image_chain_cell_indices_alt.data(tier),
                                      image_chain_cell_indices.data(tier),
                                      nt_work_groups.data(tier),
                                      reinterpret_cast<void*>(x_force.data(tier)),
                                      reinterpret_cast<void*>(y_force.data(tier)),
                                      reinterpret_cast<void*>(z_force.data(tier)),
                                      x_force_overflow.data(tier), y_force_overflow.data(tier),
                                      z_force_overflow.data(tier));
  case CoordinateCycle::WHITE:
    return CellGridWriter<void, void,
                          void, void>(theme, system_count, total_cell_count, total_chain_count,
                                      mesh_subdivisions, cell_base_capacity, localpos_scale,
                                      localpos_inverse_scale, poly_ps_ptr->getForceScalingFactor(),
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
                                      cell_migrations.data(tier),
                                      wandering_atom_count.data(tier),
                                      wandering_atom_count_alt.data(tier), wanderers.data(tier),
                                      nonimaged_atom_indices.data(tier),
                                      nonimaged_atom_indices_alt.data(tier),
                                      image_array_indices.data(tier),
                                      image_array_indices_alt.data(tier),
                                      image_chain_cell_indices.data(tier),
                                      image_chain_cell_indices_alt.data(tier),
                                      nt_work_groups.data(tier),
                                      reinterpret_cast<void*>(x_force.data(tier)),
                                      reinterpret_cast<void*>(y_force.data(tier)),
                                      reinterpret_cast<void*>(z_force.data(tier)),
                                      x_force_overflow.data(tier), y_force_overflow.data(tier),
                                      z_force_overflow.data(tier));
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
                                        mesh_subdivisions, cell_base_capacity, localpos_scale,
                                        localpos_inverse_scale,
                                        poly_ps_ptr->getInverseForceScalingFactor(),
                                        system_cell_grids.data(tier),
                                        cell_umat_ptr, cell_invu_ptr, pmig_invu_ptr,
                                        image_cell_limits_alt.data(tier),
                                        image_chain_limits.data(tier),
                                        system_chain_bounds.data(tier),
                                        chain_system_membership.data(tier),
                                        reinterpret_cast<const void*>(image_alt.data(tier)),
                                        nonimaged_atom_indices_alt.data(tier),
                                        image_array_indices_alt.data(tier),
                                        image_chain_cell_indices_alt.data(tier),
                                        nt_work_groups.data(tier),
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
                                      mesh_subdivisions, cell_base_capacity, localpos_scale,
                                      localpos_inverse_scale,
                                      poly_ps_ptr->getInverseForceScalingFactor(),
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
                                      image_array_indices.data(tier),
                                      image_chain_cell_indices.data(tier),
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
const CellOriginsReader
CellGrid<T, Tacc, Tcalc, T4>::getRulers(const CoordinateCycle orientation,
                                        const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::WHITE:
    return CellOriginsReader(origin_offset_stride, cell_origins_ax.data(tier),
                             cell_origins_ax_overflow.data(tier), cell_origins_bx.data(tier),
                             cell_origins_bx_overflow.data(tier), cell_origins_by.data(tier),
                             cell_origins_by_overflow.data(tier), cell_origins_cx.data(tier),
                             cell_origins_cx_overflow.data(tier), cell_origins_cy.data(tier),
                             cell_origins_cy_overflow.data(tier), cell_origins_cz.data(tier),
                             cell_origins_cz_overflow.data(tier));
  case CoordinateCycle::BLACK:
    return CellOriginsReader(origin_offset_stride, cell_origins_ax.data(tier),
                             alt_cell_origins_ax_overflow.data(tier),
                             alt_cell_origins_bx.data(tier),
                             alt_cell_origins_bx_overflow.data(tier),
                             alt_cell_origins_by.data(tier),
                             alt_cell_origins_by_overflow.data(tier),
                             alt_cell_origins_cx.data(tier),
                             alt_cell_origins_cx_overflow.data(tier),
                             alt_cell_origins_cy.data(tier),
                             alt_cell_origins_cy_overflow.data(tier),
                             alt_cell_origins_cz.data(tier),
                             alt_cell_origins_cz_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellOriginsReader
CellGrid<T, Tacc, Tcalc, T4>::getRulers(const HybridTargetLevel tier) const {
  return getRulers(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellOriginsWriter CellGrid<T, Tacc, Tcalc, T4>::getRulers(const CoordinateCycle orientation,
                                                          const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::WHITE:
    return CellOriginsWriter(origin_offset_stride, cell_origins_ax.data(tier),
                             cell_origins_ax_overflow.data(tier), cell_origins_bx.data(tier),
                             cell_origins_bx_overflow.data(tier), cell_origins_by.data(tier),
                             cell_origins_by_overflow.data(tier), cell_origins_cx.data(tier),
                             cell_origins_cx_overflow.data(tier), cell_origins_cy.data(tier),
                             cell_origins_cy_overflow.data(tier), cell_origins_cz.data(tier),
                             cell_origins_cz_overflow.data(tier));
  case CoordinateCycle::BLACK:
    return CellOriginsWriter(origin_offset_stride, cell_origins_ax.data(tier),
                             alt_cell_origins_ax_overflow.data(tier),
                             alt_cell_origins_bx.data(tier),
                             alt_cell_origins_bx_overflow.data(tier),
                             alt_cell_origins_by.data(tier),
                             alt_cell_origins_by_overflow.data(tier),
                             alt_cell_origins_cx.data(tier),
                             alt_cell_origins_cx_overflow.data(tier),
                             alt_cell_origins_cy.data(tier),
                             alt_cell_origins_cy_overflow.data(tier),
                             alt_cell_origins_cz.data(tier),
                             alt_cell_origins_cz_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellOriginsWriter CellGrid<T, Tacc, Tcalc, T4>::getRulers(const HybridTargetLevel tier) {
  return getRulers(cycle_position, tier);
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
  image_array_indices.upload();
  image_array_indices_alt.upload();
  image_chain_cell_indices.upload();
  image_chain_cell_indices_alt.upload();
  cell_migrations.upload();
  wandering_atom_count.upload();
  wandering_atom_count_alt.upload();
  wanderers.upload();
  nonimaged_atom_indices.upload();
  nonimaged_atom_indices_alt.upload();
  nt_work_groups.upload();
  x_force.upload();
  y_force.upload();
  z_force.upload();
  x_force_overflow.upload();
  y_force_overflow.upload();
  z_force_overflow.upload();
  origin_llint_data.upload();
  origin_int_data.upload();
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
  image_array_indices.download();
  image_array_indices_alt.download();
  image_chain_cell_indices.download();
  image_chain_cell_indices_alt.download();
  cell_migrations.download();
  wandering_atom_count.download();
  wandering_atom_count_alt.download();
  wanderers.download();
  nonimaged_atom_indices.download();
  nonimaged_atom_indices_alt.download();
  nt_work_groups.download();
  x_force.download();
  y_force.download();
  z_force.download();
  x_force_overflow.download();
  y_force_overflow.download();
  z_force_overflow.download();
  origin_llint_data.download();
  origin_int_data.download();
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
  const uint image_size = total_cell_count * cell_base_capacity;
  
  // The coordinate synthesis will receive forces at its own current cycle position.  Taking the
  // abstract with no cycle position argument defaults to the PhaseSpaceSynthesis's own setting.
  PsSynthesisWriter destw = dest->data(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const CellGridReader<T, Tacc, Tcalc, T4> cgr = data(cycle_position);
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
        if (std::type_index(typeid(Tacc)).hash_code() == int_type_index) {
          for (int i = destw.atom_starts[pos]; i < hlim; i++) {
            const uint image_idx = cgr.img_atom_idx[i];

            // Atoms from the synthesis that are not represented on the neighbor list grid will
            // have indices of -1, which becomes UINT_MAX in an unsigned int.
            if (image_idx < image_size) {
              const llint nfx = hostInt63ToLongLong(cgr.xfrc[image_idx], cgr.xfrc_ovrf[image_idx]);
              const llint nfy = hostInt63ToLongLong(cgr.yfrc[image_idx], cgr.yfrc_ovrf[image_idx]);
              const llint nfz = hostInt63ToLongLong(cgr.zfrc[image_idx], cgr.zfrc_ovrf[image_idx]);
              if (destw.frc_bits <= force_scale_nonoverflow_bits) {
                destw.xfrc[i] += nfx;
                destw.yfrc[i] += nfy;
                destw.zfrc[i] += nfz;
              }
              else {
                const int95_t snfx = hostInt95Sum(destw.xfrc[i], destw.xfrc_ovrf[i], nfx, 0);
                const int95_t snfy = hostInt95Sum(destw.yfrc[i], destw.yfrc_ovrf[i], nfy, 0);
                const int95_t snfz = hostInt95Sum(destw.zfrc[i], destw.zfrc_ovrf[i], nfz, 0);
                destw.xfrc[i] = snfx.x;
                destw.yfrc[i] = snfy.x;
                destw.zfrc[i] = snfz.x;
                destw.xfrc_ovrf[i] = snfx.y;
                destw.yfrc_ovrf[i] = snfy.y;
                destw.zfrc_ovrf[i] = snfz.y;
              }
            }
          }
        }
        else {
          for (int i = destw.atom_starts[pos]; i < hlim; i++) {
            const uint image_idx = cgr.img_atom_idx[i];
            if (image_idx < image_size) {
              const int95_t nfx = hostInt95Sum(destw.xfrc[i], destw.xfrc_ovrf[i],
                                               cgr.xfrc[image_idx], cgr.xfrc_ovrf[image_idx]);
              const int95_t nfy = hostInt95Sum(destw.yfrc[i], destw.yfrc_ovrf[i],
                                               cgr.yfrc[image_idx], cgr.yfrc_ovrf[image_idx]);
              const int95_t nfz = hostInt95Sum(destw.zfrc[i], destw.zfrc_ovrf[i],
                                               cgr.zfrc[image_idx], cgr.zfrc_ovrf[image_idx]);
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
void CellGrid<T, Tacc, Tcalc, T4>::updatePositions(PhaseSpaceSynthesis *dest,
                                                   const HybridTargetLevel tier,
                                                   const GpuDetails &gpu) {
  CellGridWriter<T, Tacc, Tcalc, T4> cgw = this->data(tier);
  PsSynthesisWriter destw = dest->data(tier);

  // Allocate space for tracking cell migrations.  The GPU will do this in chip cache on a
  // chain-by-chain basis.
  std::vector<int> flux(cgw.total_cell_count, 0);
  std::vector<uint> fill_counters(cgw.total_cell_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const int xfrm_stride = roundUp(9, warp_size_int);
      for (int sys_idx = 0; sys_idx < destw.system_count; sys_idx++) {
        const ullint c_lims = cgw.system_cell_grids[sys_idx];
        const int ncell_a = ((c_lims >> 28) & 0xfff);
        const int ncell_b = ((c_lims >> 40) & 0xfff);
        const int ncell_c = (c_lims >> 52);
        const int cell_start = (c_lims & 0xfffffff);
        const uint chain_capacity = ncell_a * static_cast<int>(cgw.cell_base_capacity);
        const uint base_img_idx = static_cast<uint>(cell_start) * cgw.cell_base_capacity;
        const double dncell_a = ncell_a;
        const double dncell_b = ncell_b;
        const double dncell_c = ncell_c;

        // Phase 1: Compute the new particle positions and cell residencies in the CellGrid.  In
        //          the coordinate synthesis, the new particle positions will be expected to
        //          already be in the "alternate" arrays of the coordinate synthesis.  Compare to
        //          the current particle positions and cell residencies.  Record the new local
        //          positions in the current image (it is an out-of-place sort aligned to the
        //          "tick-tock" mechanics of the PhaseSpaceSynthesis).  Record the cell migration
        //          as a relative move in the parallel array of one-byte keys.  This relies on an
        //          assumption that particles will not move more than one cell in any given time
        //          step.  On the GPU, this operation can be fused with the valence interactions
        //          kernel, which may read forces from the CellGrid but does not rely on the local
        //          particle positions.  Here, each phase is performed on an entire system before
        //          moving on because it is convenient to do so in serial CPU code, without
        //          repeating the derivation of the length constants above.
        const int llim = destw.atom_starts[sys_idx];
        const int hlim = llim + destw.atom_counts[sys_idx];
        const double* umat = &destw.umat_alt[xfrm_stride * sys_idx];
        const T* invu_ptr = &cgw.system_cell_invu[xfrm_stride * sys_idx];
        for (int i = llim; i < hlim; i++) {

          // Re-image the particle into the primary unit cell.  Find its current neighbor list
          // cell.
          double x, y, z;
          if (destw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
            x = static_cast<double>(destw.xalt[i]) * destw.inv_gpos_scale;
            y = static_cast<double>(destw.yalt[i]) * destw.inv_gpos_scale;
            z = static_cast<double>(destw.zalt[i]) * destw.inv_gpos_scale;
          }
          else {
            x = hostInt95ToDouble(destw.xalt[i], destw.xalt_ovrf[i]) * destw.inv_gpos_scale;
            y = hostInt95ToDouble(destw.yalt[i], destw.yalt_ovrf[i]) * destw.inv_gpos_scale;
            z = hostInt95ToDouble(destw.zalt[i], destw.zalt_ovrf[i]) * destw.inv_gpos_scale;
          }
          double ximg = (umat[0] * x) + (umat[3] * y) + (umat[6] * z);
          double yimg =                 (umat[4] * y) + (umat[7] * z);
          double zimg =                                 (umat[8] * z);
          ximg -= floor(ximg);
          yimg -= floor(yimg);
          zimg -= floor(zimg);
          ximg *= dncell_a;
          yimg *= dncell_b;
          zimg *= dncell_c;
          const int cell_aidx = ximg;
          const int cell_bidx = yimg;
          const int cell_cidx = zimg;
          ximg -= cell_aidx;
          yimg -= cell_bidx;
          zimg -= cell_cidx;
          const T xn_cell = (invu_ptr[0] * ximg) + (invu_ptr[3] * yimg) + (invu_ptr[6] * zimg);
          const T yn_cell =                        (invu_ptr[4] * yimg) + (invu_ptr[7] * zimg);
          const T zn_cell =                                               (invu_ptr[8] * zimg);
          const uint curr_img_idx = cgw.img_atom_idx[i];
          T4 xyz_prop = cgw.image[curr_img_idx];
          xyz_prop.x = xn_cell;
          xyz_prop.y = yn_cell;
          xyz_prop.z = zn_cell;
          cgw.image[curr_img_idx] = xyz_prop;
          
          // Look up the particle's index in the current image.  Determine its current cell
          // position based on this index.  The cell B- and C-axis indices can be obtained by
          // calculating the current chain residence, while the A-axis index can be obtained by
          // assuming that the particle moves no more than one cell width from the current
          // position and testing the limits accordingly.
          const uint del_img_idx = curr_img_idx - base_img_idx;
          const int curr_chain_idx = (del_img_idx / chain_capacity);
          const int curr_cell_c = (curr_chain_idx / ncell_b);
          const int curr_cell_b = curr_chain_idx - (curr_cell_c * ncell_b);
          int curr_cell_a = cell_aidx;
          int cell_guess = cell_start + (curr_chain_idx * ncell_a) + curr_cell_a;
          uint2 test_limits = cgw.cell_limits[cell_guess];
          if (curr_img_idx > test_limits.x) {
            int guess_del_idx = static_cast<int>(curr_img_idx - test_limits.x) -
                                static_cast<int>(test_limits.y >> 16);
            while (guess_del_idx >= 0) {
              curr_cell_a++;
              cell_guess++;
              test_limits = cgw.cell_limits[cell_guess];
              guess_del_idx -= static_cast<int>(test_limits.y >> 16);
            }
          }
          else {

            // If the current cell grid index is equal to the lower bound of a particular cell's
            // holdings, then it would appear the correct cell has been found.  However, this is
            // only true if the cell has holdings at all.  Otherwise, advance to the next cell that
            // does have a nonzero particle population.
            while ((test_limits.y >> 16) == 0U) {
              curr_cell_a++;
              cell_guess++;
              test_limits = cgw.cell_limits[cell_guess];
            }

            // Otherwise, the current cell is too far along in the chain.  The first cell that
            // has a lower bound not less than the current cell grid image index will have a
            // nonzero particle population and be the current home cell.
            while (curr_img_idx < test_limits.x) {
              curr_cell_a--;
              cell_guess--;
              test_limits = cgw.cell_limits[cell_guess];
            }
          }
          curr_cell_a += ((curr_cell_a < 0) - (curr_cell_a >= ncell_a)) * ncell_a;
          
          // Compute the atom migration and construct the atom's migration key.  If the key is
          // nonzero, mark the flux arrays to show that the atom is transiting from its original
          // cell into a new one.
          int migration_a = cell_aidx - curr_cell_a;
          int migration_b = cell_bidx - curr_cell_b;
          int migration_c = cell_cidx - curr_cell_c;
          if (abs(migration_a) > 1) {
            migration_a = (migration_a < 0) - (migration_a > 0);
          }
          if (abs(migration_b) > 1) {
            migration_b = (migration_b < 0) - (migration_b > 0);
          }
          if (abs(migration_c) > 1) {
            migration_c = (migration_c < 0) - (migration_c > 0);
          }
          const int mig_key =  (migration_a < 0)       + ((migration_a > 0) << 1) +
                              ((migration_b < 0) << 2) + ((migration_b > 0) << 3) +
                              ((migration_c < 0) << 4) + ((migration_c > 0) << 5);
          cgw.migration_keys[curr_img_idx] = static_cast<uchar>(mig_key);
          if (mig_key != 0) {
            flux[cell_start +
                 (((curr_cell_c * ncell_b) + curr_cell_b) * ncell_a) + curr_cell_a] -= 1;
            flux[cell_start + (((cell_cidx * ncell_b) + cell_bidx) * ncell_a) + cell_aidx] += 1;
          }
        }

        // Phase 2: Compute the new cell indexing boundaries.  This will populate the next cell
        //          limits array, first in terms of the total number of particles in each cell and
        //          then by prefix sums over each chain.  This will require a new kernel launch on
        //          the GPU, as the total numbers of atoms in each cell cannot be known until all
        //          particles' moves and new positions are computed.
        for (int k = 0; k < ncell_c; k++) {
          for (int j = 0; j < ncell_b; j++) {
            size_t cell_idx = cell_start + (((k * ncell_b) + j) * ncell_a);
            uint img_running_idx = cgw.cell_limits[cell_idx].x;
            for (int i = 0; i < ncell_a; i++) {
              const uint2 ccli = cgw.cell_limits[cell_idx];
              const uint next_particle_count = (ccli.y >> 16) + flux[cell_idx];
              const uint sysid = (ccli.y & 0xffff);
              cgw.cell_limits_alt[cell_idx] = { img_running_idx,
                                                (sysid | (next_particle_count << 16)) };
              fill_counters[cell_idx] = img_running_idx;
              img_running_idx += next_particle_count;
              flux[cell_idx] = 0;
              cell_idx++;
            }
          }
        }

        // Phase 3: Populate the new cells for the next image.  On the GPU, this will be
        //          accomplished by devoting one warp to each cell, looping over particles in the
        //          cognate cell of the original image.  Particles that do not move between cells
        //          (the majority of all cases) will get their positions in the new image based on
        //          a warp-wide reduction of all non-migrating particles in a warp vote.  The
        //          complete set of all non-migrating particles for the warp's batch will get
        //          their new indices based on a single atomicAdd() operation carried out by the
        //          first lane.  Next, migrating particles will be placed in their new cells based
        //          on indices found with atomicAdd() operations carried out by individual threads.
        //          This also requires a new kernel launch on the GPU.
        for (int k = 0; k < ncell_c; k++) {
          for (int j = 0; j < ncell_b; j++) {
            size_t cell_idx = cell_start + (((k * ncell_b) + j) * ncell_a);
            for (int i = 0; i < ncell_a; i++) {
              const uint2 ccli = cgw.cell_limits[cell_idx];
              const uint natom = (ccli.y >> 16);
              for (uint m = 0; m < natom; m++) {
                const uint curr_img_idx = ccli.x + m;
                const T4 xyz_prop = cgw.image[curr_img_idx];
                const int nonimg_idx = cgw.nonimg_atom_idx[curr_img_idx];
                const int mig_key = static_cast<int>(cgw.migration_keys[curr_img_idx]);
                int cell_maidx = i -  (mig_key &  0x1)       + ((mig_key &  0x2) >> 1);
                int cell_mbidx = j - ((mig_key &  0x4) >> 2) + ((mig_key &  0x8) >> 3);
                int cell_mcidx = k - ((mig_key & 0x10) >> 4) + ((mig_key & 0x20) >> 5);
                cell_maidx += ((cell_maidx < 0) - (cell_maidx >= ncell_a)) * ncell_a;
                cell_mbidx += ((cell_mbidx < 0) - (cell_mbidx >= ncell_b)) * ncell_b;
                cell_mcidx += ((cell_mcidx < 0) - (cell_mcidx >= ncell_c)) * ncell_c;
                const uint next_cell_idx = cell_start +
                                           (((cell_mcidx * ncell_b) + cell_mbidx) * ncell_a) +
                                           cell_maidx;
                const uint next_img_idx = fill_counters[next_cell_idx];
                cgw.image_alt[next_img_idx] = xyz_prop;
                cgw.img_atom_idx_alt[nonimg_idx] = next_img_idx;
                cgw.nonimg_atom_idx_alt[next_img_idx] = nonimg_idx;
                fill_counters[next_cell_idx] = next_img_idx + 1;
              }
              cell_idx++;
            }
          }
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::updatePositions(const HybridTargetLevel tier,
                                                   const GpuDetails &gpu) {
  updatePositions(poly_ps_ptr, tier, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::updateCyclePosition() {
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
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::updateCyclePosition(CoordinateCycle time_point) {
  cycle_position = time_point;
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
             realToString(maximum_cell_width, 8, 4, NumberFormat::STANDARD_REAL) + " will be "
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
  ushort* image_chain_cell_ptr;
  switch (cyc) {
  case CoordinateCycle::WHITE:
    cell_limits_ptr = image_cell_limits.data();
    image_ptr = image.data();
    nonimg_atom_idx_ptr = nonimaged_atom_indices.data();
    image_idx_ptr = image_array_indices.data();
    image_chain_cell_ptr = image_chain_cell_indices.data();
    break;
  case CoordinateCycle::BLACK:
    cell_limits_ptr = image_cell_limits_alt.data();
    image_ptr = image_alt.data();
    nonimg_atom_idx_ptr = nonimaged_atom_indices_alt.data();
    image_idx_ptr = image_array_indices_alt.data();
    image_chain_cell_ptr = image_chain_cell_indices_alt.data();
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
              image_chain_cell_ptr[image_idx] = i;
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
void CellGrid<T, Tacc, Tcalc, T4>::rebaseRulers() {
  const int ruler_stride = origin_offset_stride * system_count;
  cell_origins_ax.setPointer(&origin_llint_data,                 0, ruler_stride);
  cell_origins_bx.setPointer(&origin_llint_data,      ruler_stride, ruler_stride);
  cell_origins_by.setPointer(&origin_llint_data,  2 * ruler_stride, ruler_stride);
  cell_origins_cx.setPointer(&origin_llint_data,  3 * ruler_stride, ruler_stride);
  cell_origins_cy.setPointer(&origin_llint_data,  4 * ruler_stride, ruler_stride);
  cell_origins_cz.setPointer(&origin_llint_data,  5 * ruler_stride, ruler_stride);
  alt_cell_origins_ax.setPointer(&origin_llint_data,  6 * ruler_stride, ruler_stride);
  alt_cell_origins_bx.setPointer(&origin_llint_data,  7 * ruler_stride, ruler_stride);
  alt_cell_origins_by.setPointer(&origin_llint_data,  8 * ruler_stride, ruler_stride);
  alt_cell_origins_cx.setPointer(&origin_llint_data,  9 * ruler_stride, ruler_stride);
  alt_cell_origins_cy.setPointer(&origin_llint_data, 10 * ruler_stride, ruler_stride);
  alt_cell_origins_cz.setPointer(&origin_llint_data, 11 * ruler_stride, ruler_stride);
  cell_origins_ax_overflow.setPointer(&origin_int_data,                 0, ruler_stride);
  cell_origins_bx_overflow.setPointer(&origin_int_data,      ruler_stride, ruler_stride);
  cell_origins_by_overflow.setPointer(&origin_int_data,  2 * ruler_stride, ruler_stride);
  cell_origins_cx_overflow.setPointer(&origin_int_data,  3 * ruler_stride, ruler_stride);
  cell_origins_cy_overflow.setPointer(&origin_int_data,  4 * ruler_stride, ruler_stride);
  cell_origins_cz_overflow.setPointer(&origin_int_data,  5 * ruler_stride, ruler_stride);
  alt_cell_origins_ax_overflow.setPointer(&origin_int_data,  6 * ruler_stride, ruler_stride);
  alt_cell_origins_bx_overflow.setPointer(&origin_int_data,  7 * ruler_stride, ruler_stride);
  alt_cell_origins_by_overflow.setPointer(&origin_int_data,  8 * ruler_stride, ruler_stride);
  alt_cell_origins_cx_overflow.setPointer(&origin_int_data,  9 * ruler_stride, ruler_stride);
  alt_cell_origins_cy_overflow.setPointer(&origin_int_data, 10 * ruler_stride, ruler_stride);
  alt_cell_origins_cz_overflow.setPointer(&origin_int_data, 11 * ruler_stride, ruler_stride);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::allocateRulers() {

  // Get the total number of cells along each direction across all systems
  int tmp_data_stride = 0;
  for (int i = 0; i < system_count; i++) {
    const ullint sys_grids = system_cell_grids.readHost(i);
    const int cell_na = ((sys_grids >> 28) & 0xfff);
    const int cell_nb = ((sys_grids >> 40) & 0xfff);
    const int cell_nc = (sys_grids >> 52);
    tmp_data_stride = std::max(tmp_data_stride, roundUp(cell_na + 1, warp_size_int));
    tmp_data_stride = std::max(tmp_data_stride, roundUp(cell_nb + 1, warp_size_int));
    tmp_data_stride = std::max(tmp_data_stride, roundUp(cell_nc + 1, warp_size_int));
  }
  origin_offset_stride = tmp_data_stride;
  origin_llint_data.resize(12 * tmp_data_stride * system_count);
  origin_int_data.resize(12 * tmp_data_stride * system_count);
  rebaseRulers();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::drawNeighborListRulers() {
  allocateRulers();
  const std::vector<CoordinateCycle> stages = { CoordinateCycle::WHITE, CoordinateCycle::BLACK };
  for (size_t icyc = 0; icyc < stages.size(); icyc++) {
    PsSynthesisWriter poly_psw = poly_ps_ptr->data(stages[icyc]);
    const int xfrm_stride = roundUp(9, warp_size_int);
    for (int i = 0; i < system_count; i++) {
      const ullint sys_grids = system_cell_grids.readHost(i);
      const int cell_na = ((sys_grids >> 28) & 0xfff);
      const int cell_nb = ((sys_grids >> 40) & 0xfff);
      const int cell_nc = (sys_grids >> 52);
      const double cell_dna = cell_na;
      const double cell_dnb = cell_nb;
      const double cell_dnc = cell_nc;

      // Build the markers along the unit cell A axis.  Cartesian X coordinates are all that is
      // required.  The integer representation of the box vectors is authoritative, but the real
      // representation is needed to determine an initial guess for the length.
      double* invu_ptr;
      llint* boxv_ptr;
      int* boxv_ovrf_ptr;
      switch (stages[icyc]) {
      case CoordinateCycle::WHITE:
        invu_ptr = &poly_psw.invu[(i * xfrm_stride)];
        boxv_ptr = &poly_psw.boxvecs[(i * xfrm_stride)];
        boxv_ovrf_ptr = &poly_psw.boxvec_ovrf[(i * xfrm_stride)];
        break;
      case CoordinateCycle::BLACK:
        invu_ptr = &poly_psw.invu_alt[(i * xfrm_stride)];
        boxv_ptr = &poly_psw.alt_boxvecs[(i * xfrm_stride)];
        boxv_ovrf_ptr = &poly_psw.alt_boxvec_ovrf[(i * xfrm_stride)];
        break;
      }
      const double cell_ax = invu_ptr[0] / cell_dna;
      const double cell_bx = invu_ptr[3] / cell_dna;
      const double cell_by = invu_ptr[4] / cell_dnb;
      const double cell_cx = invu_ptr[6] / cell_dna;
      const double cell_cy = invu_ptr[7] / cell_dnb;
      const double cell_cz = invu_ptr[8] / cell_dnc;
      const int95_t cell_iax = hostDoubleToInt95(cell_ax * poly_psw.gpos_scale);
      const int95_t cell_ibx = hostDoubleToInt95(cell_bx * poly_psw.gpos_scale);
      const int95_t cell_iby = hostDoubleToInt95(cell_by * poly_psw.gpos_scale);
      const int95_t cell_icx = hostDoubleToInt95(cell_cx * poly_psw.gpos_scale);
      const int95_t cell_icy = hostDoubleToInt95(cell_cy * poly_psw.gpos_scale);
      const int95_t cell_icz = hostDoubleToInt95(cell_cz * poly_psw.gpos_scale);
      const int95_t sys_iax = { boxv_ptr[0], boxv_ovrf_ptr[0] };
      const int95_t sys_ibx = { boxv_ptr[3], boxv_ovrf_ptr[3] };
      const int95_t sys_iby = { boxv_ptr[4], boxv_ovrf_ptr[4] };
      const int95_t sys_icx = { boxv_ptr[6], boxv_ovrf_ptr[6] };
      const int95_t sys_icy = { boxv_ptr[7], boxv_ovrf_ptr[7] };
      const int95_t sys_icz = { boxv_ptr[8], boxv_ovrf_ptr[8] };
      const int95_t sys_iax_est = hostSplitFPMult(cell_iax, cell_na);
      const int95_t sys_ibx_est = hostSplitFPMult(cell_ibx, cell_na);
      const int95_t sys_iby_est = hostSplitFPMult(cell_iby, cell_nb);
      const int95_t sys_icx_est = hostSplitFPMult(cell_icx, cell_na);
      const int95_t sys_icy_est = hostSplitFPMult(cell_icy, cell_nb);
      const int95_t sys_icz_est = hostSplitFPMult(cell_icz, cell_nc);
      const int95_t iax_err = hostSplitFPSubtract(sys_iax, sys_iax_est);
      const int95_t ibx_err = hostSplitFPSubtract(sys_ibx, sys_ibx_est);
      const int95_t iby_err = hostSplitFPSubtract(sys_iby, sys_iby_est);
      const int95_t icx_err = hostSplitFPSubtract(sys_icx, sys_icx_est);
      const int95_t icy_err = hostSplitFPSubtract(sys_icy, sys_icy_est);
      const int95_t icz_err = hostSplitFPSubtract(sys_icz, sys_icz_est);

      // Set pointers based on the stage of the coordinate cycle;
      llint *ax_ptr, *bx_ptr, *by_ptr, *cx_ptr, *cy_ptr, *cz_ptr;
      int *ax_ovrf_ptr, *bx_ovrf_ptr, *by_ovrf_ptr, *cx_ovrf_ptr, *cy_ovrf_ptr, *cz_ovrf_ptr;
      switch (stages[icyc]) {
      case CoordinateCycle::WHITE:
        ax_ptr = cell_origins_ax.data();
        bx_ptr = cell_origins_bx.data();
        by_ptr = cell_origins_by.data();
        cx_ptr = cell_origins_cx.data();
        cy_ptr = cell_origins_cy.data();
        cz_ptr = cell_origins_cz.data();
        ax_ovrf_ptr = cell_origins_ax_overflow.data();
        bx_ovrf_ptr = cell_origins_bx_overflow.data();
        by_ovrf_ptr = cell_origins_by_overflow.data();
        cx_ovrf_ptr = cell_origins_cx_overflow.data();
        cy_ovrf_ptr = cell_origins_cy_overflow.data();
        cz_ovrf_ptr = cell_origins_cz_overflow.data();
        break;
      case CoordinateCycle::BLACK:
        ax_ptr = alt_cell_origins_ax.data();
        bx_ptr = alt_cell_origins_bx.data();
        by_ptr = alt_cell_origins_by.data();
        cx_ptr = alt_cell_origins_cx.data();
        cy_ptr = alt_cell_origins_cy.data();
        cz_ptr = alt_cell_origins_cz.data();
        ax_ovrf_ptr = alt_cell_origins_ax_overflow.data();
        bx_ovrf_ptr = alt_cell_origins_bx_overflow.data();
        by_ovrf_ptr = alt_cell_origins_by_overflow.data();
        cx_ovrf_ptr = alt_cell_origins_cx_overflow.data();
        cy_ovrf_ptr = alt_cell_origins_cy_overflow.data();
        cz_ovrf_ptr = alt_cell_origins_cz_overflow.data();
        break;
      }
      
      // The error will be completely contained in the primary component of the fixed-precision
      // representation.  Double-precision transformations, and the fact that the format cannot
      // extend the primary unit cell markers very far into the overflow accumulator even at the
      // limits of the allowed precision, ensure that the remainder is small.
      loadRuler(cell_na, cell_iax, { boxv_ptr[0], boxv_ovrf_ptr[0] }, iax_err.x,
                &ax_ptr[(i * origin_offset_stride)], &ax_ovrf_ptr[(i * origin_offset_stride)]);
      loadRuler(cell_nb, cell_ibx, { boxv_ptr[3], boxv_ovrf_ptr[3] }, ibx_err.x,
                &bx_ptr[(i * origin_offset_stride)], &bx_ovrf_ptr[(i * origin_offset_stride)]);
      loadRuler(cell_nb, cell_iby, { boxv_ptr[4], boxv_ovrf_ptr[4] }, iby_err.x,
                &by_ptr[(i * origin_offset_stride)], &by_ovrf_ptr[(i * origin_offset_stride)]);
      loadRuler(cell_nc, cell_icx, { boxv_ptr[6], boxv_ovrf_ptr[6] }, icx_err.x,
                &cx_ptr[(i * origin_offset_stride)], &cx_ovrf_ptr[(i * origin_offset_stride)]);
      loadRuler(cell_nc, cell_icy, { boxv_ptr[7], boxv_ovrf_ptr[7] }, icy_err.x,
                &cy_ptr[(i * origin_offset_stride)], &cy_ovrf_ptr[(i * origin_offset_stride)]);
      loadRuler(cell_nc, cell_icz, { boxv_ptr[8], boxv_ovrf_ptr[8] }, icz_err.x,
                &cz_ptr[(i * origin_offset_stride)], &cz_ovrf_ptr[(i * origin_offset_stride)]);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::prepareWorkGroups() {
  std::vector<int> twu(32, 0);
  std::vector<int> rel_a = {  0,  0,  0,  0,  0, -2, -1,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1 };
  std::vector<int> rel_b = {  0,  0,  0,  0,  0, -2, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0 };
  std::vector<int> rel_c = { -2, -1,  0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
  nt_work_groups.resize(static_cast<size_t>(total_cell_count) * 32LLU);
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
    for (int j = 0; j < 17; j++) {
      int ra_idx = a_idx + rel_a[j];
      int rb_idx = b_idx + rel_b[j];
      int rc_idx = c_idx + rel_c[j];
      ra_idx += ((ra_idx < 0) - (ra_idx >= icell_na)) * icell_na;
      rb_idx += ((rb_idx < 0) - (rb_idx >= icell_nb)) * icell_nb;
      rc_idx += ((rc_idx < 0) - (rc_idx >= icell_nc)) * icell_nc;
      twu[j + ((j >= 5) * 11)] = icell_offset + (((rc_idx * icell_nb) + rb_idx) * icell_na) +
                                 ra_idx;
    }

    // Indices 6-15 (and 29-31) are available to store additional information on the work unit.
    twu[29] = system_idx;
    nt_work_groups.putHost(twu, 32LLU * static_cast<ullint>(i), 32);
  }
}

//-------------------------------------------------------------------------------------------------
#if 0
template <typename T, typename Tacc, typename Tcalc, typename T4>
void CellGrid<T, Tacc, Tcalc, T4>::prepareCullingField() {

}
#endif

//-------------------------------------------------------------------------------------------------
template <typename Tsrc, typename Tcalc, typename Tcalc2>
Tcalc sourceMagnitude(const NonbondedTheme requested_property, const NonbondedTheme cg_content,
                      const Tsrc q, const bool q_is_real, const int sysid,
                      const SyNonbondedKit<Tcalc, Tcalc2> &synbk) {

  // The switch over the density theme is only necessary to distinguish what to collect if the
  // cell grid serves a non-bonded theme of "ALL."  The correct source of density to be represented
  // on the particle-mesh interaction grid is found by examining the cell grid's contents and
  // then choosing the proper interpretation of its contents.
  switch (requested_property) {
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
      rtErr(getEnumerationName(requested_property) + " properties cannot be extracted from a "
            "source value with " + getEnumerationName(cg_content) + " information.",
            "sourceMagnitude");
    case NonbondedTheme::ALL:

      // The source data will necessarily have a signed integer format, even if it is recast as a
      // floating point number to fit the CellGrid's T4 (Tcoord x 4) tuple type.
      if (sizeof(Tsrc) == 8) {
        int q_idx;
        if (q_is_real) {
          const Ecumenical8 conv = { .d = static_cast<double>(q) };
          q_idx = (conv.lli & dp_charge_index_mask);
        }
        else {
          q_idx = (static_cast<llint>(q) & dp_charge_index_mask);
        }
        return synbk.q_params[q_idx];
      }
      else {
        int q_idx;
        if (q_is_real) {
          const Ecumenical4 conv = { .f = static_cast<float>(q) };
          q_idx = (conv.i & sp_charge_index_mask);
        }
        else {
          q_idx = (static_cast<int>(q) & sp_charge_index_mask);
        }
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
        rtErr(getEnumerationName(requested_property) + " properties cannot be extracted from a "
              "source value with " + getEnumerationName(cg_content) + " information.",
              "sourceMagnitude");
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
          if (q_is_real) {
            const Ecumenical8 conv = { .d = static_cast<double>(q) };
            tlj_idx = (conv.lli >> dp_charge_index_bits);
          }
          else {
            tlj_idx = (static_cast<llint>(q) >> dp_charge_index_bits);
          }
        }
        else {
          if (q_is_real) {
            const Ecumenical4 conv = { .f = static_cast<float>(q) };
            tlj_idx = (conv.i >> sp_charge_index_bits);
          }
          else {
            tlj_idx = (static_cast<int>(q) >> sp_charge_index_bits);
          }
        }
        break;
      }
      const Tcalc2 t_ljab = synbk.ljab_coeff[synbk.ljabc_offsets[sysid] +
                                             (tlj_idx * (synbk.n_lj_types[sysid] + 1))];
      if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
        return sqrt(0.25 * t_ljab.y);
      }
      else {
        return sqrtf(0.25f * t_ljab.y);
      }
    }
    break;
  case NonbondedTheme::ALL:
    rtErr("Only one non-bonded output property is allowed.", "sourceMagnitude");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tsrc>
int sourceIndex(const NonbondedTheme requested_property, const NonbondedTheme cg_content,
                const Tsrc q, const bool q_is_real) {

  // The requested property must be accessed through the CellGrid object by an index, which will
  // be the case for any van-der Waals property but may or may not be the case for an electrostatic
  // property (charge magnitude).
  int prop_idx;
  switch (requested_property) {
  case NonbondedTheme::ELECTROSTATIC:
    switch (cg_content) {
    case NonbondedTheme::ELECTROSTATIC:
      if (q_is_real) {
        rtErr("Real-valued charge data in an " + getEnumerationName(cg_content) + " CellGrid "
              "object cannot produce an index.", "sourceIndex");
      }
      else {
        prop_idx = q;
      }
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      rtErr(getEnumerationName(requested_property) + " parameter indices cannot be extracted "
            "from a source value with " + getEnumerationName(cg_content) + " information.",
            "sourceIndex");
    case NonbondedTheme::ALL:
      if (q_is_real) {
        if (sizeof(Tsrc) == 8) {
          if (q_is_real) {
            const Ecumenical8 conv = { .d = static_cast<double>(q) };
            prop_idx = (conv.lli & dp_charge_index_mask);
          }
          else {
            prop_idx = (static_cast<llint>(q) & dp_charge_index_mask);
          }
        }
        else {
          if (q_is_real) {
            const Ecumenical4 conv = { .f = static_cast<float>(q) };
            prop_idx = (conv.i & sp_charge_index_mask);
          }
          else {
            prop_idx = (static_cast<int>(q) & sp_charge_index_mask);
          }
        }
      }
      else {
        if (sizeof(Tsrc) == 8) {
          prop_idx = (static_cast<llint>(q) & dp_charge_index_mask);
        }
        else {
          prop_idx = (static_cast<int>(q) & sp_charge_index_mask);
        }
      }
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    {
      switch (cg_content) {
      case NonbondedTheme::ELECTROSTATIC:
        rtErr(getEnumerationName(requested_property) + " parameter indices cannot be extracted "
              "from a source value with " + getEnumerationName(cg_content) + " information.",
              "sourceIndex");
      case NonbondedTheme::VAN_DER_WAALS:
        if (q_is_real) {
          if (sizeof(Tsrc) == 8) {
            const Ecumenical8 conv = { .d = static_cast<double>(q) };
            prop_idx = conv.lli;
          }
          else {
            const Ecumenical4 conv = { .f = static_cast<float>(q) };
            prop_idx = conv.i;
          }
        }
        else {
          prop_idx = q;
        }
        break;
      case NonbondedTheme::ALL:

        // The source data will necessarily have a signed integer format.
        if (sizeof(Tsrc) == 8) {
          if (q_is_real) {
            const Ecumenical8 conv = { .d = static_cast<double>(q) };
            prop_idx = (conv.lli >> dp_charge_index_bits);
          }
          else {
            prop_idx = (static_cast<llint>(q) >> dp_charge_index_bits);
          }
        }
        else {
          if (q_is_real) {
            const Ecumenical4 conv = { .f = static_cast<float>(q) };
            prop_idx = (conv.i >> sp_charge_index_bits);
          }
          else {
            prop_idx = (static_cast<int>(q) >> sp_charge_index_bits);
          }
        }
        break;
      }
    }
    break;
  case NonbondedTheme::ALL:
    rtErr("Only one non-bonded output property is allowed.", "sourceMagnitude");
  }
  return prop_idx;
}
                  
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4> restoreType(CellGridWriter<void, void, void, void> *rasa) {
  return CellGridWriter<T, Tacc,
                        Tcalc, T4>(rasa->theme, rasa->system_count, rasa->total_cell_count,
                                   rasa->total_chain_count, rasa->mesh_ticks,
                                   rasa->cell_base_capacity, rasa->lpos_scale,
                                   rasa->inv_lpos_scale, rasa->frc_scale, rasa->system_cell_grids,
                                   reinterpret_cast<Tcalc*>(rasa->system_cell_umat),
                                   reinterpret_cast<T*>(rasa->system_cell_invu),
                                   reinterpret_cast<T*>(rasa->system_pmig_invu),
                                   rasa->cell_limits, rasa->cell_limits_alt, rasa->chain_limits,
                                   rasa->system_chain_bounds, rasa->chain_system_owner,
                                   reinterpret_cast<T4*>(rasa->image),
                                   reinterpret_cast<T4*>(rasa->image_alt),
                                   rasa->migration_keys, rasa->wander_count,
                                   rasa->wander_count_alt, rasa->wanderers, rasa->nonimg_atom_idx,
                                   rasa->nonimg_atom_idx_alt, rasa->img_atom_idx,
                                   rasa->img_atom_idx_alt, rasa->img_atom_chn_cell,
                                   rasa->img_atom_chn_cell_alt, rasa->nt_groups,
                                   reinterpret_cast<Tacc*>(rasa->xfrc),
                                   reinterpret_cast<Tacc*>(rasa->yfrc),
                                   reinterpret_cast<Tacc*>(rasa->zfrc), rasa->xfrc_ovrf,
                                   rasa->yfrc_ovrf, rasa->zfrc_ovrf);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4> restoreType(CellGridWriter<void, void, void, void> &rasa) {
  return CellGridWriter<T, Tacc,
                        Tcalc, T4>(rasa.theme, rasa.system_count, rasa.total_cell_count,
                                   rasa.total_chain_count, rasa.mesh_ticks,
                                   rasa.cell_base_capacity, rasa.lpos_scale, rasa.inv_lpos_scale,
                                   rasa.frc_scale, rasa.system_cell_grids,
                                   reinterpret_cast<Tcalc*>(rasa.system_cell_umat),
                                   reinterpret_cast<T*>(rasa.system_cell_invu),
                                   reinterpret_cast<T*>(rasa.system_pmig_invu),
                                   rasa.cell_limits, rasa.cell_limits_alt, rasa.chain_limits,
                                   rasa.system_chain_bounds, rasa.chain_system_owner,
                                   reinterpret_cast<T4*>(rasa.image),
                                   reinterpret_cast<T4*>(rasa.image_alt),
                                   rasa.migration_keys, rasa.wander_count, rasa.wander_count_alt,
                                   rasa.wanderers, rasa.nonimg_atom_idx, rasa.nonimg_atom_idx_alt,
                                   rasa.img_atom_idx, rasa.img_atom_idx_alt,
                                   rasa.img_atom_chn_cell, rasa.img_atom_chn_cell_alt,
                                   rasa.nt_groups, reinterpret_cast<Tacc*>(rasa.xfrc),
                                   reinterpret_cast<Tacc*>(rasa.yfrc),
                                   reinterpret_cast<Tacc*>(rasa.zfrc), rasa.xfrc_ovrf,
                                   rasa.yfrc_ovrf, rasa.zfrc_ovrf);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
restoreType(const CellGridReader<void, void, void, void> *rasa) {
  return CellGridReader<T, Tacc,
                        Tcalc, T4>(rasa->theme, rasa->system_count, rasa->total_cell_count,
                                   rasa->total_chain_count, rasa->mesh_ticks,
                                   rasa->cell_base_capacity, rasa->lpos_scale,
                                   rasa->inv_lpos_scale, rasa->inv_frc_scale,
                                   rasa->system_cell_grids,
                                   reinterpret_cast<const Tcalc*>(rasa->system_cell_umat),
                                   reinterpret_cast<const T*>(rasa->system_cell_invu),
                                   reinterpret_cast<const T*>(rasa->system_pmig_invu),
                                   rasa->cell_limits, rasa->chain_limits,
                                   rasa->system_chain_bounds, rasa->chain_system_owner,
                                   reinterpret_cast<const T4*>(rasa->image),
                                   rasa->nonimg_atom_idx, rasa->img_atom_idx,
                                   rasa->img_atom_chn_cell, rasa->nt_groups,
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
                                   rasa.cell_base_capacity, rasa.lpos_scale, rasa.inv_lpos_scale,
                                   rasa.inv_frc_scale, rasa.system_cell_grids,
                                   reinterpret_cast<const Tcalc*>(rasa.system_cell_umat),
                                   reinterpret_cast<const T*>(rasa.system_cell_invu),
                                   reinterpret_cast<const T*>(rasa.system_pmig_invu),
                                   rasa.cell_limits, rasa.chain_limits, rasa.system_chain_bounds,
                                   rasa.chain_system_owner,
                                   reinterpret_cast<const T4*>(rasa.image), rasa.nonimg_atom_idx,
                                   rasa.img_atom_idx, rasa.img_atom_chn_cell, rasa.nt_groups,
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
