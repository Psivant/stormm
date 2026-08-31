// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tgrid, typename Tcalc, typename Tacc>
void pullPMEForces(const Tgrid* a_cof, const Tgrid* b_cof, const Tgrid* c_cof, const Tgrid* da_cof,
                   const Tgrid* db_cof, const Tgrid* dc_cof, int bspline_order,
                   const Tgrid* pme_potential, const int grid_root_a, const int grid_root_b,
                   const int grid_root_c, const uint4 grid_dims, const FFTMode fft_staging,
                   const Tcalc* umat, const uint cg_img_index, Tacc* cg_xfrc, Tacc* cg_yfrc,
                   Tacc* cg_zfrc, int* cg_xovrf, int* cg_yovrf, int* cg_zovrf, Tacc* cg_net_frc,
                   int* cg_net_ovrf, const Tcalc cg_frc_scl) {
  const bool tacc_is_llint = (std::type_index(typeid(Tacc)).hash_code() == llint_type_index);
  uint padded_gdim_x;
  switch (fft_staging) {
  case FFTMode::IN_PLACE:
    padded_gdim_x = 2 * ((grid_dims.x / 2) + 1);
    break;
  case FFTMode::OUT_OF_PLACE:
    padded_gdim_x = grid_dims.x;
    break;
  }
  const Tcalc value_zero = 0.0;
  Tcalc fa = value_zero;
  Tcalc fb = value_zero;
  Tcalc fc = value_zero;
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
        fa -= da_cof[i] * b_cof[j] * c_cof[k] * pme_potential[gidx];
        fb -= a_cof[i] * db_cof[j] * c_cof[k] * pme_potential[gidx];
        fc -= a_cof[i] * b_cof[j] * dc_cof[k] * pme_potential[gidx];
      }
    }
  }
  fa *= static_cast<Tcalc>(grid_dims.x);
  fb *= static_cast<Tcalc>(grid_dims.y);
  fc *= static_cast<Tcalc>(grid_dims.z);
  const Tcalc fx = (umat[0] * fa);
  const Tcalc fy = (umat[3] * fa) + (umat[4] * fb);
  const Tcalc fz = (umat[6] * fa) + (umat[7] * fb) + (umat[8] * fc);
      
  // Accumulate the results
  if (tacc_is_llint) {
    const int95_t ifx = hostDoubleToInt95(fx * cg_frc_scl);
    const int95_t ify = hostDoubleToInt95(fy * cg_frc_scl);
    const int95_t ifz = hostDoubleToInt95(fz * cg_frc_scl);
    const int95_t itfx = hostSplitFPSum(ifx, cg_xfrc[cg_img_index], cg_xovrf[cg_img_index]);
    const int95_t itfy = hostSplitFPSum(ify, cg_yfrc[cg_img_index], cg_yovrf[cg_img_index]);
    const int95_t itfz = hostSplitFPSum(ifz, cg_zfrc[cg_img_index], cg_zovrf[cg_img_index]);
    cg_xfrc[cg_img_index] = itfx.x;
    cg_yfrc[cg_img_index] = itfy.x;
    cg_zfrc[cg_img_index] = itfz.x;
    cg_xovrf[cg_img_index] = itfx.y;
    cg_yovrf[cg_img_index] = itfy.y;
    cg_zovrf[cg_img_index] = itfz.y;
    const int95_t inetfx = hostSplitFPSum(ifx, cg_net_frc[0], cg_net_ovrf[0]);
    const int95_t inetfy = hostSplitFPSum(ify, cg_net_frc[1], cg_net_ovrf[1]);
    const int95_t inetfz = hostSplitFPSum(ifz, cg_net_frc[2], cg_net_ovrf[2]);
    cg_net_frc[0] = inetfx.x;
    cg_net_frc[1] = inetfy.x;
    cg_net_frc[2] = inetfz.x;
    cg_net_ovrf[0] = inetfx.y;
    cg_net_ovrf[1] = inetfy.y;
    cg_net_ovrf[2] = inetfz.y;
  }
  else {
    const int2 ifx = hostDoubleToInt63(fx * cg_frc_scl);
    const int2 ify = hostDoubleToInt63(fy * cg_frc_scl);
    const int2 ifz = hostDoubleToInt63(fz * cg_frc_scl);
    const int2 itfx = hostSplitFPSum(ifx, cg_xfrc[cg_img_index], cg_xovrf[cg_img_index]);
    const int2 itfy = hostSplitFPSum(ify, cg_yfrc[cg_img_index], cg_yovrf[cg_img_index]);
    const int2 itfz = hostSplitFPSum(ifz, cg_zfrc[cg_img_index], cg_zovrf[cg_img_index]);
    cg_xfrc[cg_img_index] = itfx.x;
    cg_yfrc[cg_img_index] = itfy.x;
    cg_zfrc[cg_img_index] = itfz.x;
    cg_xovrf[cg_img_index] = itfx.y;
    cg_yovrf[cg_img_index] = itfy.y;
    cg_zovrf[cg_img_index] = itfz.y;
    const int2 inetfx = hostSplitFPSum(ifx, cg_net_frc[0], cg_net_ovrf[0]);
    const int2 inetfy = hostSplitFPSum(ify, cg_net_frc[1], cg_net_ovrf[1]);
    const int2 inetfz = hostSplitFPSum(ifz, cg_net_frc[2], cg_net_ovrf[2]);
    cg_net_frc[0] = inetfx.x;
    cg_net_frc[1] = inetfy.x;
    cg_net_frc[2] = inetfz.x;
    cg_net_ovrf[0] = inetfx.y;
    cg_net_ovrf[1] = inetfy.y;
    cg_net_ovrf[2] = inetfz.y;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename Tcalc2, typename T4>
void gatherCellForces(CellGridWriter<T, Tacc, Tcalc, T4> *cgw, const PMIGridReader &pmigr,
                      const PsSynthesisBorders &pssb, const int sysid, const int cell_i,
                      const int cell_j, const int cell_k,
                      const SyNonbondedKit<Tcalc, Tcalc2> &synbk) {

  // Re-derive the cell and cell grid boundaries rather than copying them via input arguments
  const ullint cell_bounds = cgw->system_cell_grids[sysid];
  const int cell_offset = (cell_bounds & 0xfffffffLLU);
  const int cell_na = ((cell_bounds >> 28) & 0xfffLLU);
  const int cell_nb = ((cell_bounds >> 40) & 0xfffLLU);

  // Determine limits and other critical constants
  const bool coord_in_real = (cgw->lpos_scale < 1.01);
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const size_t ijk_cellidx = cell_offset + (((cell_k * cell_nb) + cell_j) * cell_na) + cell_i;
  const uint2 ijk_bounds = cgw->cell_limits[ijk_cellidx];
  const uint mllim = ijk_bounds.x;
  const uint mhlim = mllim + (ijk_bounds.y >> 16);
  const int xfrm_stride = roundUp(9, warp_size_int);

  // Lay out arrays to collect B-spline coefficients
  const uint4 grid_dims = pmigr.dims[sysid];
  Tcalc system_box_umat[9];
  for (int i = 0; i < 9; i++) {
    system_box_umat[i] = pssb.umat[(sysid * xfrm_stride) + i];
  }
  switch (pmigr.mode) {
  case PrecisionModel::DOUBLE:
    {
      std::vector<double> a_cof(pmigr.order), b_cof(pmigr.order), c_cof(pmigr.order);
      std::vector<double> da_cof(pmigr.order), db_cof(pmigr.order), dc_cof(pmigr.order);
      for (uint m = mllim; m < mhlim; m++) {
        const T4 atom_m = cgw->image[m];
        int grid_root_a, grid_root_b, grid_root_c;
        particleAlignment<Tcalc, double>(atom_m.x, atom_m.y, atom_m.z, cgw->inv_lpos_scale,
                                         &cgw->system_cell_umat[sysid * xfrm_stride],
                                         cgw->mesh_ticks, cell_i, cell_j, cell_k, a_cof.data(),
                                         b_cof.data(), c_cof.data(), pmigr.order, &grid_root_a,
                                         &grid_root_b, &grid_root_c, da_cof.data(), db_cof.data(),
                                         dc_cof.data());
        const double q = sourceMagnitude<T, Tcalc, Tcalc2>(pmigr.theme, cgw->theme, atom_m.w,
                                                           coord_in_real, sysid, synbk);
        for (int i = 0; i < pmigr.order; i++) {
          a_cof[i] *= q;
          da_cof[i] *= q;
        }
        pullPMEForces<double, Tcalc, Tacc>(a_cof.data(), b_cof.data(), c_cof.data(), da_cof.data(),
                                           db_cof.data(), dc_cof.data(), pmigr.order, pmigr.ddata,
                                           grid_root_a, grid_root_b, grid_root_c, grid_dims,
                                           pmigr.fftm, system_box_umat, m, cgw->xfrc, cgw->yfrc,
                                           cgw->zfrc, cgw->xfrc_ovrf, cgw->yfrc_ovrf,
                                           cgw->zfrc_ovrf, cgw->net_frc, cgw->net_frc_ovrf,
                                           cgw->frc_scale);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      std::vector<float> a_cof(pmigr.order), b_cof(pmigr.order), c_cof(pmigr.order);
      std::vector<float> da_cof(pmigr.order), db_cof(pmigr.order), dc_cof(pmigr.order);
      for (uint m = mllim; m < mhlim; m++) {
        const T4 atom_m = cgw->image[m];
        int grid_root_a, grid_root_b, grid_root_c;
        particleAlignment<Tcalc, float>(atom_m.x, atom_m.y, atom_m.z, cgw->inv_lpos_scale,
                                        &cgw->system_cell_umat[(sysid * xfrm_stride)],
                                        cgw->mesh_ticks, cell_i, cell_j, cell_k, a_cof.data(),
                                        b_cof.data(), c_cof.data(), pmigr.order, &grid_root_a,
                                        &grid_root_b, &grid_root_c, da_cof.data(), db_cof.data(),
                                        dc_cof.data());
        const float q = sourceMagnitude<T, Tcalc, Tcalc2>(pmigr.theme, cgw->theme, atom_m.w,
                                                          coord_in_real, sysid, synbk);
        for (int i = 0; i < pmigr.order; i++) {
          a_cof[i] *= q;
          da_cof[i] *= q;
        }
        pullPMEForces<float, Tcalc, Tacc>(a_cof.data(), b_cof.data(), c_cof.data(), da_cof.data(),
                                          db_cof.data(), dc_cof.data(), pmigr.order, pmigr.fdata,
                                          grid_root_a, grid_root_b, grid_root_c, grid_dims,
                                          pmigr.fftm, system_box_umat, m, cgw->xfrc, cgw->yfrc,
                                          cgw->zfrc, cgw->xfrc_ovrf, cgw->yfrc_ovrf,
                                          cgw->zfrc_ovrf, cgw->net_frc, cgw->net_frc_ovrf,
                                          cgw->frc_scale);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tnb_calc, typename Tnb_calc2, typename T4>
void gatherForces(CellGridWriter<T, Tacc, Tnb_calc, T4> *cgw, const PMIGridReader &pmigr,
                  const PsSynthesisBorders &pssb,
                  const SyNonbondedKit<Tnb_calc, Tnb_calc2> &synbk) {
  for (int sysid = 0; sysid < synbk.nsys; sysid++) {
    const ullint cell_bounds = cgw->system_cell_grids[sysid];
    const int cell_na = ((cell_bounds >> 28) & 0xfffLLU);
    const int cell_nb = ((cell_bounds >> 40) & 0xfffLLU);
    const int cell_nc = ((cell_bounds >> 52) & 0xfffLLU);
    for (int i = 0; i < cell_na; i++) {
      for (int j = 0; j < cell_nb; j++) {
        for (int k = 0; k < cell_nc; k++) {
          gatherCellForces<T, Tacc, Tnb_calc, Tnb_calc2, T4>(cgw, pmigr, pssb, sysid, i, j, k,
                                                             synbk);
        }
      }
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void gatherForces(CellGrid<T, Tacc, Tcalc, T4> *cg, const PMIGrid *pm,
                  const AtomGraphSynthesis *poly_ag) {

  // As in the density mapping API, the need to "unroll" the templating of the CellGrid object
  // arises from a choice to build the API around the PMIGrid and topology synthesis.  The PMIGrid
  // object contains a pointer to a CellGrid object, but its template parameters must be recovered.
  // This form of the function will be restricted to work on the CPU.
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcrd_is_real = isFloatingPointScalarType<T>();
  const SyNonbondedKit<double, double2> dsynbk = poly_ag->getDoublePrecisionNonbondedKit();
  const SyNonbondedKit<float, float2>  fsynbk = poly_ag->getSinglePrecisionNonbondedKit();
  CellGridWriter<void, void, void, void> cgw_v = cg->templateFreeData();
  CellGridWriter<T, Tacc, double, T4> dcgw = restoreType<T, Tacc, double, T4>(cgw_v);
  CellGridWriter<T, Tacc, float, T4>  fcgw = restoreType<T, Tacc, float, T4>(cgw_v);
  const PsSynthesisBorders pssb = cg->getUnitCellTransforms();
  const PMIGridReader pmigr = pm->data();

  // Outside of the real-valued parameter arrays, dsynbk and fsynbk will hold equivalent sizing
  // constants and atom indexing.  Use either to manage loops and bounds until the precision
  // model demands a bifurcation of further work.
  if (tcalc_is_double) {
    gatherForces<T, Tacc, double, double2, T4>(&dcgw, pmigr, pssb, dsynbk);
  }
  else {
    gatherForces<T, Tacc, float, float2, T4>(&fcgw, pmigr, pssb, fsynbk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void gatherForces(CellGrid<T, Tacc, Tcalc, T4> *cg, const PMIGrid &pm,
                  const AtomGraphSynthesis &poly_ag) {
  gatherForces(cg, pm.getSelfPointer(), poly_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
std::vector<double3> gatherForces(const CoordinateFrameReader &cfr, const Tdata* potential_field,
                                  const NonbondedKit<Tcalc> &nbk, const NonbondedTheme theme,
                                  const FFTMode fft_staging, const int grid_dim_a,
                                  const int grid_dim_b, const int grid_dim_c, const int order,
                                  const BSplineUnity unification) {

  // As with the prior density mapping operation, trap bad inputs.
  if (cfr.natom != nbk.natom) {
    rtErr("The number of atoms in the coordinate system (" + std::to_string(cfr.natom) + ") must "
          "match the number of atoms in the topology (" + std::to_string(nbk.natom) + ").",
          "gatherForces");
  }
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
  case NonbondedTheme::VAN_DER_WAALS:
    break;
  case NonbondedTheme::ALL:
    rtErr("Particles must interact via grid-mediated forces in a sanctioned non-bonded "
          "potential. \"" + getEnumerationName(theme) + "\" is invalid.", "gatherForces");
  }

  // Transform the particle positions into the unit cell space to position each on the grid.
  std::vector<double> frac_x(cfr.natom), frac_y(cfr.natom), frac_z(cfr.natom);
  for (int i = 0; i < cfr.natom; i++) {
    frac_x[i] = cfr.xcrd[i];
    frac_y[i] = cfr.ycrd[i];
    frac_z[i] = cfr.zcrd[i];
  }
  for (int i = 0; i < cfr.natom; i++) {
    frac_x[i] = (cfr.umat[0] * frac_x[i]) + (cfr.umat[3] * frac_y[i]) + (cfr.umat[6] * frac_z[i]);
    frac_y[i] =                             (cfr.umat[4] * frac_y[i]) + (cfr.umat[7] * frac_z[i]);
    frac_z[i] =                                                         (cfr.umat[8] * frac_z[i]);
    frac_x[i] -= floor(frac_x[i]);
    frac_y[i] -= floor(frac_y[i]);
    frac_z[i] -= floor(frac_z[i]);
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

  // Calculate the number of padded grid indices.
  int padded_dim_a;
  switch (fft_staging) {
  case FFTMode::IN_PLACE:
    padded_dim_a = 2 * ((grid_dim_a / 2) + 1);
    break;
  case FFTMode::OUT_OF_PLACE:
    padded_dim_a = grid_dim_a;
    break;
  }

  // Compute the particles' B-spline coefficients and derivatives thereof, then pull forces from
  // the grid.
  const Tcalc value_zero = 0.0;
  std::vector<double3> result(nbk.natom);
  std::vector<Tcalc> bspl_a(order), bspl_b(order), bspl_c(order);
  std::vector<Tcalc> dbspl_a(order), dbspl_b(order), dbspl_c(order);
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
      bSpline(da, order, bspl_a.data(), dbspl_a.data());
      bSpline(db, order, bspl_b.data(), dbspl_b.data());
      bSpline(dc, order, bspl_c.data(), dbspl_c.data());
      break;
    case BSplineUnity::NONE:
      bspl_a = bSplineNoUnity<Tcalc>(da, order);
      bspl_b = bSplineNoUnity<Tcalc>(db, order);
      bspl_c = bSplineNoUnity<Tcalc>(dc, order);
      dbspl_a = dBSpline(da, order, false);
      dbspl_b = dBSpline(db, order, false);
      dbspl_c = dBSpline(dc, order, false);
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
      dbspl_a[i] *= q;
    }
    Tcalc tmp_rslt_x = value_zero;
    Tcalc tmp_rslt_y = value_zero;
    Tcalc tmp_rslt_z = value_zero;
    for (int k = 0; k < order; k++) {
      int act_k = base_gc + k;
      act_k += ((act_k < 0) - (act_k >= grid_dim_c)) * grid_dim_c;
      for (int j = 0; j < order; j++) {
	int act_j = base_gb + j;
        act_j += ((act_j < 0) - (act_j >= grid_dim_b)) * grid_dim_b;
        const int jk_idx = ((act_k * grid_dim_b) + act_j) * padded_dim_a;
        for (int i = 0; i < order; i++) {
          int act_i = base_ga + i;
          act_i += ((act_i < 0) - (act_i >= grid_dim_a)) * grid_dim_a;
          tmp_rslt_x += dbspl_a[i] * bspl_b[j] * bspl_c[k];
          tmp_rslt_y += bspl_a[i] * dbspl_b[j] * bspl_c[k];
          tmp_rslt_z += bspl_a[i] * bspl_b[j] * dbspl_c[k];
        }
      }
    }
    result[pos] = { static_cast<double>(tmp_rslt_x), static_cast<double>(tmp_rslt_y),
                    static_cast<double>(tmp_rslt_z) };
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tdata>
std::vector<double3> gatherForces(const CoordinateFrame *cf, const Tdata* potential_field,
                                  const AtomGraph *ag, const NonbondedTheme theme,
                                  const FFTMode fft_staging, const int grid_dim_a,
                                  const int grid_dim_b, const int grid_dim_c, const int order,
                                  const PrecisionModel prec, const BSplineUnity unification) {
  const CoordinateFrameReader cfr = cf->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
      return gatherForces<Tdata, double>(cfr, potential_field, nbk, theme, fft_staging,
                                         grid_dim_a, grid_dim_b, grid_dim_c, order, unification);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const NonbondedKit<float> nbk = ag->getSinglePrecisionNonbondedKit();
      return gatherForces<Tdata, float>(cfr, potential_field, nbk, theme, fft_staging, grid_dim_a,
                                        grid_dim_b, grid_dim_c, order, unification);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata>
std::vector<double3> gatherForces(const CoordinateFrame &cf, const Tdata* potential_field,
                                  const AtomGraph &ag, const NonbondedTheme theme,
                                  const FFTMode fft_staging, const int grid_dim_a,
                                  const int grid_dim_b, const int grid_dim_c, const int order,
                                  const PrecisionModel prec, const BSplineUnity unification) {
  return gatherForces(cf.getSelfPointer(), potential_field, ag.getSelfPointer(), theme,
                      fft_staging, grid_dim_a, grid_dim_b, grid_dim_c, order, prec);
}


//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
void unrollGatherForcesCall(PMIGrid *pm, const size_t cg_tacc, const size_t cg_tcalc,
                            const AtomGraphSynthesis *poly_ag) {
  if (cg_tacc == int_type_index) {
    unrollGatherForcesCall<T, int, T4>(pm, cg_tcalc, poly_ag);
  }
  else if (cg_tacc == llint_type_index) {
    unrollGatherForcesCall<T, llint, T4>(pm, cg_tcalc, poly_ag);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename T4>
void unrollGatherForcesCall(PMIGrid *pm, const size_t cg_tcalc,
                            const AtomGraphSynthesis *poly_ag) {
  if (cg_tcalc == double_type_index) {
    const CellGrid<T, Tacc, double, T4> *c_cgp = pm->getCellGridPointer<T, Tacc, double, T4>();
    CellGrid<T, Tacc, double, T4> *cgp = const_cast<CellGrid<T, Tacc, double, T4>*>(c_cgp);
    gatherForces<T, Tacc, double, T4>(cgp, pm, poly_ag);
  }
  else if (cg_tcalc == float_type_index) {
    const CellGrid<T, Tacc, float, T4> *c_cgp = pm->getCellGridPointer<T, Tacc, float, T4>();
    CellGrid<T, Tacc, float, T4> *cgp = const_cast<CellGrid<T, Tacc, float, T4>*>(c_cgp);
    gatherForces<T, Tacc, float, T4>(cgp, pm, poly_ag);
  }
}

} // namespace energy
} // namespace stormm
