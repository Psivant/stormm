// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc maxClashingDistance(const NonbondedKit<Tcalc> &nbk, const Tcalc elec_limit,
                          const Tcalc vdw_ratio) {
  Tcalc result = elec_limit;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    for (int j = 0; j <= i; j++) {
      result = std::max(result, nbk.lj_sigma[(i * nbk.n_lj_types) + j]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool trivialClashCheck(const std::vector<Tcalc> &cachi_xcrd, const std::vector<Tcalc> &cachj_xcrd,
                       const std::vector<Tcalc> &cachi_ycrd, const std::vector<Tcalc> &cachj_ycrd,
                       const std::vector<Tcalc> &cachi_zcrd, const std::vector<Tcalc> &cachj_zcrd,
                       const int ni_atoms, const int nj_atoms, const Tcalc max_clash) {
  Tcalc min_ix = cachi_xcrd[0];
  Tcalc min_iy = cachi_ycrd[0];
  Tcalc min_iz = cachi_zcrd[0];
  Tcalc max_ix = cachi_xcrd[0];
  Tcalc max_iy = cachi_ycrd[0];
  Tcalc max_iz = cachi_zcrd[0];
  for (int i = 1; i < ni_atoms; i++) {
    min_ix = std::min(cachi_xcrd[i], min_ix);
    min_iy = std::min(cachi_ycrd[i], min_iy);
    min_iz = std::min(cachi_zcrd[i], min_iz);
    max_ix = std::max(cachi_xcrd[i], max_ix);
    max_iy = std::max(cachi_ycrd[i], max_iy);
    max_iz = std::max(cachi_zcrd[i], max_iz);
  }
  Tcalc min_jx = cachj_xcrd[0];
  Tcalc min_jy = cachj_ycrd[0];
  Tcalc min_jz = cachj_zcrd[0];
  Tcalc max_jx = cachj_xcrd[0];
  Tcalc max_jy = cachj_ycrd[0];
  Tcalc max_jz = cachj_zcrd[0];
  for (int j = 1; j < nj_atoms; j++) {
    min_jx = std::min(cachj_xcrd[j], min_jx);
    min_jy = std::min(cachj_ycrd[j], min_jy);
    min_jz = std::min(cachj_zcrd[j], min_jz);
    max_jx = std::max(cachj_xcrd[j], max_jx);
    max_jy = std::max(cachj_ycrd[j], max_jy);
    max_jz = std::max(cachj_zcrd[j], max_jz);
  }
  min_ix -= max_clash;
  min_iy -= max_clash;
  min_iz -= max_clash;
  max_ix += max_clash;
  max_iy += max_clash;
  max_iz += max_clash;
  return (max_ix > min_jx && min_ix < max_jx && max_iy > min_jy && min_iy < max_jy &&
          max_iz > min_jz && min_iz < max_jz);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool directClashTesting(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &maskr,
                        const Tcalc elec_limit, const Tcalc vdw_ratio, const Tcalc inv_scale,
                        ClashReport *summary) {

  // Determine the calculation type and a maximum distance parameter
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc max_clash = maxClashingDistance<Tcalc>(nbk, elec_limit, vdw_ratio);
  
  // Lay out a workspace and begin searching
  bool clash_found = false;
  std::vector<Tcalc> cachi_xcrd(tile_length), cachi_ycrd(tile_length), cachi_zcrd(tile_length);
  std::vector<Tcalc> cachj_xcrd(tile_length), cachj_ycrd(tile_length), cachj_zcrd(tile_length);
  std::vector<int> cachi_ljidx(tile_length), cachj_ljidx(tile_length);
  for (int sti = 0; sti < maskr.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {
      const int stni_atoms = std::min(maskr.natom - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(maskr.natom - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;
      const int diag_supertile = (sti == stj);
      const int stij_map_index = maskr.supertile_map_idx[(stj * maskr.supertile_stride_count) +
                                                         sti];
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);
          const int tij_map_index = maskr.tile_map_idx[stij_map_index +
                                                       (tj * tile_lengths_per_supertile) + ti];
          if (isFloatingPointScalarType<Tcoord>()) {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              cachi_xcrd[i]  = xcrd[atom_i];
              cachi_ycrd[i]  = ycrd[atom_i];
              cachi_zcrd[i]  = zcrd[atom_i];
              cachi_ljidx[i] = nbk.lj_idx[atom_i];
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              cachj_xcrd[j]  = xcrd[atom_j];
              cachj_ycrd[j]  = ycrd[atom_j];
              cachj_zcrd[j]  = zcrd[atom_j];
              cachj_ljidx[j] = nbk.lj_idx[atom_j];
            }
          }
          else {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              cachi_xcrd[i]  = xcrd[atom_i] * inv_scale;
              cachi_ycrd[i]  = ycrd[atom_i] * inv_scale;
              cachi_zcrd[i]  = zcrd[atom_i] * inv_scale;
              cachi_ljidx[i] = nbk.lj_idx[atom_i];
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              cachj_xcrd[j]  = xcrd[atom_j] * inv_scale;
              cachj_ycrd[j]  = ycrd[atom_j] * inv_scale;
              cachj_zcrd[j]  = zcrd[atom_j] * inv_scale;
              cachj_ljidx[j] = nbk.lj_idx[atom_j];
            }
          }

          // Perform a trivial clash check
          if (diag_tile || (diag_supertile && ti - tj < 3) ||
              trivialClashCheck<Tcalc>(cachi_xcrd, cachj_xcrd, cachi_ycrd, cachj_ycrd, cachi_zcrd,
                                       cachj_zcrd, ni_atoms, nj_atoms, max_clash)) {
            for (int i = 0; i < ni_atoms; i++) {
              const int ljt_i = cachi_ljidx[i];
              const uint mask_i = maskr.mask_data[tij_map_index + i];
              const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
              if (tcalc_is_double) {
                for (int j = 0; j < jlim; j++) {
                  if ((mask_i >> j) & 0x1) {
                    continue;
                  }
                  const int ljt_ij = (nbk.n_lj_types * ljt_i) + cachj_ljidx[j];
                  const Tcalc dx = cachj_xcrd[j] - cachi_xcrd[i];
                  const Tcalc dy = cachj_ycrd[j] - cachi_ycrd[i];
                  const Tcalc dz = cachj_zcrd[j] - cachi_zcrd[i];
                  const Tcalc r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                  if (r < nbk.lj_sigma[ljt_ij] * vdw_ratio || r < elec_limit) {
                    clash_found = true;
                    if (summary != nullptr) {
                      const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
                      const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
                      summary->addClash(atom_i, atom_j, r);
                    }
                  }
                }
              }
              else {
                for (int j = 0; j < jlim; j++) {
                  if ((mask_i >> j) & 0x1) {
                    continue;
                  }
                  const int ljt_ij = (nbk.n_lj_types * ljt_i) + cachj_ljidx[j];
                  const Tcalc dx = cachj_xcrd[j] - cachi_xcrd[i];
                  const Tcalc dy = cachj_ycrd[j] - cachi_ycrd[i];
                  const Tcalc dz = cachj_zcrd[j] - cachi_zcrd[i];
                  const Tcalc r = sqrtf((dx * dx) + (dy * dy) + (dz * dz));
                  if (r < nbk.lj_sigma[ljt_ij] * vdw_ratio || r < elec_limit) {
                    clash_found = true;
                    if (summary != nullptr) {
                      const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
                      const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
                      summary->addClash(atom_i, atom_j, r);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return clash_found;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
int3 clashGridDecomposition(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                            const int natom, double3 *grid_origin, double3 *grid_lengths) {
  
  // Use the atom count to set an upper limit as to the number of grid cells.  The minimum is 512
  // (a grid of 8 x 8 x 8, or some rearrangement of that many cells), the maximum will be
  // determined by placing up to sixteen atoms in each cell.
  const llint max_cells = std::max(natom / 16, 216);
  Tcoord min_x = xcrd[0];
  Tcoord max_x = xcrd[0];
  Tcoord min_y = ycrd[0];
  Tcoord max_y = ycrd[0];
  Tcoord min_z = zcrd[0];
  Tcoord max_z = zcrd[0];
  for (int i = 1; i < natom; i++) {
    min_x = std::min(xcrd[i], min_x);
    max_x = std::max(xcrd[i], max_x);
    min_y = std::min(ycrd[i], min_y);
    max_y = std::max(ycrd[i], max_y);
    min_z = std::min(zcrd[i], min_z);
    max_z = std::max(zcrd[i], max_z);
  }
  grid_origin->x = min_x;
  grid_origin->y = min_y;
  grid_origin->z = min_z;
  grid_lengths->x = max_x - min_x;
  grid_lengths->y = max_y - min_y;
  grid_lengths->z = max_z - min_z;

  // Find the manner in which the grid cells can be partitioned which makes them as close as
  // possible to cubes.
  const Tcoord len_x = max_x - min_x;
  const Tcoord len_y = max_y - min_y;
  const Tcoord len_z = max_z - min_z;
  const Tcoord ratio_xy = len_y / len_x;
  const Tcoord ratio_xz = len_z / len_x;
  const Tcoord prod_ratio_xyz = ratio_xy * ratio_xz;
  const Tcoord dmax_cells = max_cells;
  const int guess_xdim = round(cbrt(dmax_cells / prod_ratio_xyz));
  llint min_diff = max_cells;
  int3 result;
  for (int i = guess_xdim - 2; i < guess_xdim + 2; i++) {
    if (i < 1) {
      continue;
    }
    llint nx = i;
    llint ny = round(i * ratio_xy);
    llint nz = round(i * ratio_xz);
    const llint nxyz  = nx * ny * nz;
    const llint ndiff = llabs(max_cells - nxyz);
    if (ndiff < min_diff || min_diff == max_cells) {
      min_diff = ndiff;
      result.x = nx;
      result.y = ny;
      result.z = nz;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool clashCellToCell(const int icell_othr, const int jcell_othr, const int kcell_othr,
                     const int icell_home, const int jcell_home, const int kcell_home,
                     const int3 grid_cells, const Tcoord* xcrd, const Tcoord* ycrd,
                     const Tcoord* zcrd, const std::vector<int> &cell_contents,
                     const std::vector<int> &cell_content_bounds, const NonbondedKit<Tcalc> &nbk,
                     const StaticExclusionMask *mask, const Tcalc elec_limit,
                     const Tcalc vdw_ratio, const Tcalc inv_scale, ClashReport *summary) {

  // Bail out immediately if the "other" cell is off the grid (non-periodic boundary conditions are
  // assumed).  The home cell is assumed to be on the grid.
  if (icell_othr < 0 || icell_othr >= grid_cells.x || jcell_othr < 0 ||
      jcell_othr >= grid_cells.y || kcell_othr < 0 || kcell_othr >= grid_cells.z) {
    return false;
  }

  // Determine the two cell indices
  const int home_idx = (((kcell_home * grid_cells.y) + jcell_home) * grid_cells.x) + icell_home;
  const int othr_idx = (((kcell_othr * grid_cells.y) + jcell_othr) * grid_cells.x) + icell_othr;

  // Loop over all atoms in each cell
  bool result = false;
  bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const int ilim = cell_content_bounds[home_idx + 1];
  const Tcalc sq_elec = elec_limit * elec_limit;
  for (int i = cell_content_bounds[home_idx]; i < ilim; i++) {
    const int iatom  = cell_contents[i];
    const int iljidx = nbk.lj_idx[iatom];
    const Tcalc ixloc = xcrd[iatom];
    const Tcalc iyloc = ycrd[iatom];
    const Tcalc izloc = zcrd[iatom];
    const int jlim = cell_content_bounds[othr_idx + 1];
    for (int j = cell_content_bounds[othr_idx]; j < jlim; j++) {
      const int jatom = cell_contents[j];
      const Tcalc dij_xloc = xcrd[jatom] - ixloc;
      const Tcalc dij_yloc = ycrd[jatom] - iyloc;
      const Tcalc dij_zloc = zcrd[jatom] - izloc;
      const Tcalc sq_dist = (dij_xloc * dij_xloc) + (dij_yloc * dij_yloc) + (dij_zloc * dij_zloc);
      const int jljidx = nbk.lj_idx[jatom];
      const Tcalc ij_sigma = vdw_ratio * nbk.lj_sigma[(iljidx * nbk.n_lj_types) + jljidx];
      if (sq_dist < ij_sigma * ij_sigma || sq_dist < sq_elec) {

        // Check that this non-bonded interaction, which appears to be clashing, is not in fact an
        // exclusion.
        if (mask->testExclusion(iatom, jatom) == false) {
          result = true;
          if (summary != nullptr) {
            summary->addClash(iatom, jatom, (tcalc_is_double) ? sqrt(sq_dist) : sqrtf(sq_dist));
          }
        }
      }
    }
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, const Tcalc elec_limit,
                 const Tcalc vdw_ratio, const Tcalc inv_scale, ClashReport *summary) {

  // Set the summary to make use of the supplied clash limits.  In many cases this will overwrite
  // the summary's internal parameters with values that it supplied in the first place.
  if (summary != nullptr) {
    summary->setMinimumDistance(elec_limit);
    summary->setMinimumSigmaRatio(vdw_ratio);
    summary->setTopologyPointer(mask->getTopologyPointer());
    summary->clear();
  }
  
  // Determine the calculation type
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  bool clash_found = false;
  
  // Do the direct calculation if the size is very small
  if (nbk.natom < clash_direct_calculation_size_limit) {
    clash_found = directClashTesting(xcrd, ycrd, zcrd, nbk, mask->data(), vdw_ratio, elec_limit,
                                     inv_scale, summary);
  }
  else {
    
    // Find the largest clashing distance.
    const Tcalc max_clash = maxClashingDistance<Tcalc>(nbk, elec_limit, vdw_ratio);
    if (max_clash >= constants::tiny) {

      // Design a grid based on the maximum clash distance.  In order to be worthwhile, the grid
      // cells must be at least max_clash thick between all parallel faces, and six or more
      // grid cells in all directions.
      double3 grid_origin, grid_lengths;
      const int3 grid_cells = clashGridDecomposition<Tcoord>(xcrd, ycrd, zcrd, nbk.natom,
                                                             &grid_origin, &grid_lengths);
      const double3 dgrid_cells = { static_cast<double>(grid_cells.x),
                                    static_cast<double>(grid_cells.y),
                                    static_cast<double>(grid_cells.z) };
      const int ncell = grid_cells.x * grid_cells.y * grid_cells.z;
      std::vector<int> cell_content_bounds(ncell + 1, 0);
      std::vector<int> cell_homes(nbk.natom);
      std::vector<int> cell_contents(nbk.natom);
      for (int i = 0; i < nbk.natom; i++) {
        const int icellx = ((xcrd[i] - grid_origin.x) / grid_lengths.x) * dgrid_cells.x;
        const int icelly = ((ycrd[i] - grid_origin.y) / grid_lengths.y) * dgrid_cells.y;
        const int icellz = ((zcrd[i] - grid_origin.z) / grid_lengths.z) * dgrid_cells.z;
        cell_homes[i] = (((icellz * grid_cells.y) + icelly) * grid_cells.x) + icellx;
      }
      indexingArray(cell_homes, &cell_contents, &cell_content_bounds);
      for (int i = 0; i < grid_cells.x; i++) {
        for (int j = 0; j < grid_cells.y; j++) {
          for (int k = 0; k < grid_cells.z; k++) {
            clash_found = (clash_found ||
                           clashCellToCell(i - 1, j - 1, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i - 1,     j, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i - 1, j + 1, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(    i, j - 1, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(    i,     j, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(    i, j + 1, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i + 1, j - 1, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i + 1,     j, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i + 1, j + 1, k - 1, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i - 1, j - 1,     k, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i - 1,     j,     k, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(i - 1, j + 1,     k, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(    i, j - 1,     k, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
            clash_found = (clash_found ||
                           clashCellToCell(    i,     j,     k, i, j, k, grid_cells, xcrd, ycrd,
                                           zcrd, cell_contents, cell_content_bounds, nbk, mask,
                                           vdw_ratio, elec_limit, inv_scale, summary));
          }
        }
      }
    }
  }

  // Loop over all 1:4 interactions to clean up possible clashes among these pairs.
  for (int pos = 0; pos < vk.ndihe; pos++) {
    if (vk.dihe14_param_idx[pos] == 0) {
      continue;
    }
    const int atom_i = vk.dihe_i_atoms[pos];
    const int atom_l = vk.dihe_l_atoms[pos];
    const int ljt_il = (nbk.n_lj_types * nbk.lj_idx[atom_i]) + nbk.lj_idx[atom_l];
    const Tcalc dx = xcrd[atom_l] - xcrd[atom_i];
    const Tcalc dy = ycrd[atom_l] - ycrd[atom_i];
    const Tcalc dz = zcrd[atom_l] - zcrd[atom_i];
    const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                        sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    if (r < nbk.lj_14_sigma[ljt_il] * vdw_ratio) {
      clash_found = true;
      if (summary != nullptr) {
        summary->addClash(atom_i, atom_l, r);
      }
    }
  }
  for (int pos = 0; pos < vk.ninfr14; pos++) {
    if (vk.infr14_param_idx[pos] == 0) {
      continue;
    }
    const int atom_i = vk.infr14_i_atoms[pos];
    const int atom_l = vk.infr14_l_atoms[pos];
    const int ljt_il = (nbk.n_lj_types * nbk.lj_idx[atom_i]) + nbk.lj_idx[atom_l];
    const Tcalc dx = xcrd[atom_l] - xcrd[atom_i];
    const Tcalc dy = ycrd[atom_l] - ycrd[atom_i];
    const Tcalc dz = zcrd[atom_l] - zcrd[atom_i];
    const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                        sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    if (r < nbk.lj_14_sigma[ljt_il] * vdw_ratio) {
      clash_found = true;
      if (summary != nullptr) {
        summary->addClash(atom_i, atom_l, r);
      }
    }
  }
  
  // If this point has been reached, no clash was detected.
  return clash_found;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeriesReader<Tcoord> &csr, const size_t frame,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, const Tcalc elec_limit,
                 const Tcalc vdw_ratio, ClashReport *summary) {
  const size_t padded_atoms = roundUp(csr.natom, warp_size_int);
  const size_t atom_start = frame * padded_atoms;
  const size_t xfrm_start = frame * roundUp<size_t>(9, warp_size_zu);
  return detectClash<Tcoord, Tcalc>(&csr.xcrd[atom_start], &csr.ycrd[atom_start],
                                    &csr.zcrd[atom_start], vk, nbk, mask, elec_limit, vdw_ratio,
                                    csr.inv_gpos_scale, summary);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> *cs, const int frame, const AtomGraph *ag,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary) {
  const size_t ct = std::type_index(typeid(Tcalc)).hash_code();
  if (ct == float_type_index) {
    return detectClash<Tcoord, float>(cs->data(), frame, ag->getSinglePrecisionValenceKit(),
                                      ag->getSinglePrecisionNonbondedKit(), mask, elec_limit,
                                      vdw_ratio, summary);
  }
  else {
    return detectClash<Tcoord, double>(cs->data(), frame, ag->getDoublePrecisionValenceKit(),
                                       ag->getDoublePrecisionNonbondedKit(), mask, elec_limit,
                                       vdw_ratio, summary);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> &cs, const int frame, const AtomGraph &ag,
                 const StaticExclusionMask &mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary) {
  return detectClash<Tcoord, Tcalc>(cs.getSelfPointer(), frame, ag.getSelfPointer(),
                                    mask.getSelfPointer(), elec_limit, vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> *cs, const int frame, const AtomGraph *ag,
                 const StaticExclusionMask *mask, ClashReport *summary) {
  if (summary == nullptr) {
    return detectClash<Tcoord, Tcalc>(cs, frame, ag, mask, default_minimize_clash_r0,
                                      default_minimize_clash_ratio);
  }
  else {
    return detectClash<Tcoord, Tcalc>(cs, frame, ag, mask, summary->getMinimumDistance(),
                                      summary->getMinimumSigmaRatio(), summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> &cs, const int frame, const AtomGraph &ag,
                 const StaticExclusionMask &mask, ClashReport *summary) {
  if (summary == nullptr) {
    return detectClash<Tcoord, Tcalc>(cs.getSelfPointer(), frame, ag.getSelfPointer(),
                                      mask.getSelfPointer(), default_minimize_clash_r0,
                                      default_minimize_clash_ratio);
  }
  else {
    return detectClash<Tcoord, Tcalc>(cs.getSelfPointer(), frame, ag.getSelfPointer(),
                                      mask.getSelfPointer(), summary->getMinimumDistance(),
                                      summary->getMinimumSigmaRatio(), summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const PsSynthesisReader &poly_psr, const int system_index,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary) {
  const int atom_start = poly_psr.atom_starts[system_index];
  const int natom = poly_psr.atom_counts[system_index];
  const int atom_limit = atom_start + natom;
  if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    std::vector<float> xcrd(natom);
    std::vector<float> ycrd(natom);
    std::vector<float> zcrd(natom);
    for (int i = atom_start; i < atom_limit; i++) {
      xcrd[i - atom_start] = static_cast<float>(poly_psr.xcrd[i]) * poly_psr.inv_gpos_scale_f;
      ycrd[i - atom_start] = static_cast<float>(poly_psr.ycrd[i]) * poly_psr.inv_gpos_scale_f;
      zcrd[i - atom_start] = static_cast<float>(poly_psr.zcrd[i]) * poly_psr.inv_gpos_scale_f;
    }
    return detectClash<float, Tcalc>(xcrd.data(), ycrd.data(), zcrd.data(), vk, nbk, mask,
                                     elec_limit, vdw_ratio, 1.0, summary);
  }
  else {
    std::vector<double> xcrd(natom);
    std::vector<double> ycrd(natom);
    std::vector<double> zcrd(natom);
    for (int i = atom_start; i < atom_limit; i++) {
      xcrd[i - atom_start] = hostInt95ToDouble(poly_psr.xcrd[i], poly_psr.xcrd_ovrf[i]) *
                             poly_psr.inv_gpos_scale;
      ycrd[i - atom_start] = hostInt95ToDouble(poly_psr.ycrd[i], poly_psr.ycrd_ovrf[i]) *
                             poly_psr.inv_gpos_scale;
      zcrd[i - atom_start] = hostInt95ToDouble(poly_psr.zcrd[i], poly_psr.zcrd_ovrf[i]) *
                             poly_psr.inv_gpos_scale;
    }
    return detectClash<double, Tcalc>(xcrd.data(), ycrd.data(), zcrd.data(), vk, nbk, mask,
                                      elec_limit, vdw_ratio, 1.0, summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis *poly_ps, const int system_index,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary) {
  const AtomGraph *ag_ptr = poly_ps->getSystemTopologyPointer(system_index);
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    return detectClash<double>(poly_ps->data(), system_index,
                               ag_ptr->getDoublePrecisionValenceKit(),
                               ag_ptr->getDoublePrecisionNonbondedKit(), mask, elec_limit,
                               vdw_ratio, summary);
  }
  else {
    return detectClash<float>(poly_ps->data(), system_index,
                              ag_ptr->getSinglePrecisionValenceKit(),
                              ag_ptr->getSinglePrecisionNonbondedKit(), mask, elec_limit,
                              vdw_ratio, summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis &poly_ps, const int system_index,
                 const StaticExclusionMask &mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary) {
  return detectClash<Tcalc>(poly_ps.getSelfPointer(), system_index, mask.getSelfPointer(),
                            elec_limit, vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis *poly_ps, const int system_index,
                 const StaticExclusionMask *mask, ClashReport *summary) {
  if (summary == nullptr) {
    return detectClash<Tcalc>(poly_ps, system_index, mask, default_minimize_clash_r0,
                              default_minimize_clash_ratio);
  }
  else {
    return detectClash<Tcalc>(poly_ps, system_index, mask, summary->getMinimumDistance(),
                              summary->getMinimumSigmaRatio(), summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis &poly_ps, const int system_index,
                 const StaticExclusionMask &mask, ClashReport *summary) {
  if (summary == nullptr) {
    return detectClash<Tcalc>(poly_ps.getSelfPointer(), system_index, mask.getSelfPointer(),
                              default_minimize_clash_r0, default_minimize_clash_ratio);
  }
  else {
    return detectClash<Tcalc>(poly_ps.getSelfPointer(), system_index, mask.getSelfPointer(),
                              summary->getMinimumDistance(), summary->getMinimumSigmaRatio(),
                              summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const CondensateReader &cdnsr, const int system_index,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary) {
  const size_t atom_start = cdnsr.atom_starts[system_index];
  switch (cdnsr.mode) {
  case PrecisionModel::DOUBLE:
    return detectClash<double, Tcalc>(&cdnsr.xcrd[atom_start], &cdnsr.ycrd[atom_start],
                                      &cdnsr.zcrd[atom_start], vk, nbk, mask, elec_limit,
                                      vdw_ratio, 1.0, summary);
  case PrecisionModel::SINGLE:
    return detectClash<float, Tcalc>(&cdnsr.xcrd_sp[atom_start], &cdnsr.ycrd_sp[atom_start],
                                     &cdnsr.zcrd_sp[atom_start], vk, nbk, mask, elec_limit,
                                     vdw_ratio, 1.0, summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const Condensate *cdns, const int system_index, const AtomGraph *ag,
                 const StaticExclusionMask *mask, const double elec_limit,
                 const double vdw_ratio, ClashReport *summary) {
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    return detectClash<double>(cdns->data(), system_index, ag->getDoublePrecisionValenceKit(),
                               ag->getDoublePrecisionNonbondedKit(), mask, elec_limit, vdw_ratio,
                               summary);
  }
  else {
    return detectClash<float>(cdns->data(), system_index, ag->getSinglePrecisionValenceKit(),
                              ag->getSinglePrecisionNonbondedKit(), mask, elec_limit, vdw_ratio,
                              summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const Condensate &cdns, const int system_index, const AtomGraph &ag,
                 const StaticExclusionMask &mask, const double elec_limit,
                 const double vdw_ratio, ClashReport *summary) {
  return detectClash<Tcalc>(cdns.getSelfPointer(), system_index, ag.getSelfPointer(),
                            mask.getSelfPointer(), elec_limit, vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const Condensate *cdns, const int system_index, const AtomGraph *ag,
                 const StaticExclusionMask *mask, ClashReport *summary) {
  if (summary == nullptr) {
    return detectClash<Tcalc>(cdns, system_index, ag, mask, default_minimize_clash_r0,
                              default_minimize_clash_ratio);
  }
  else {
    return detectClash<Tcalc>(cdns, system_index, ag, mask, summary->getMinimumDistance(),
                              summary->getMinimumSigmaRatio(), summary);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
bool detectClash(const Condensate &cdns, const int system_index, const AtomGraph &ag,
                 const StaticExclusionMask &mask, ClashReport *summary) {
  return detectClash<Tcalc>(cdns.getSelfPointer(), system_index, ag.getSelfPointer(),
                            mask.getSelfPointer(), summary);
}

} // namespace structure
} // namespace stormm
