#include <algorithm>
#include "copyright.h"
#include "Math/rounding.h"
#include "Reporting/error_format.h"
#include "brickwork.h"

namespace stormm {
namespace synthesis {

using stmath::roundUp;
using stmath::partition;
  
//-------------------------------------------------------------------------------------------------
VolumePartition::VolumePartition(const int a_dim_in, const int b_dim_in, const int c_dim_in,
                                 const int a_orig_in, const int b_orig_in, const int c_orig_in,
                                 const int halo_under_in, const int halo_over_in,
                                 const int system_index_in) :
    a_dim{a_dim_in}, b_dim{b_dim_in}, c_dim{c_dim_in}, a_orig{a_orig_in}, b_orig{b_orig_in},
    c_orig{c_orig_in}, halo_under{halo_under_in}, halo_over{halo_over_in},
    system_index{system_index_in}
{}

//-------------------------------------------------------------------------------------------------
int3 VolumePartition::getLength() const {
  return { a_dim, b_dim, c_dim };
}

//-------------------------------------------------------------------------------------------------
int VolumePartition::getLength(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return a_dim;
  case UnitCellAxis::B:
    return b_dim;
  case UnitCellAxis::C:
    return c_dim;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int VolumePartition::getLength(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return a_dim;
  case CartesianDimension::Y:
    return b_dim;
  case CartesianDimension::Z:
    return c_dim;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int3 VolumePartition::getOrigin() const {
  return { a_orig, b_orig, c_orig };
}

//-------------------------------------------------------------------------------------------------
int3 VolumePartition::getLengthPlusHalo() const {
  const int all_halo = halo_over + halo_under;
  return { a_dim + all_halo, b_dim + all_halo, c_dim + all_halo };
}

//-------------------------------------------------------------------------------------------------
int VolumePartition::getLengthPlusHalo(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return a_dim + halo_under + halo_over;
  case UnitCellAxis::B:
    return b_dim + halo_under + halo_over;
  case UnitCellAxis::C:
    return c_dim + halo_under + halo_over;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int VolumePartition::getLengthPlusHalo(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return a_dim + halo_under + halo_over;
  case CartesianDimension::Y:
    return b_dim + halo_under + halo_over;
  case CartesianDimension::Z:
    return c_dim + halo_under + halo_over;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
  int VolumePartition::getVolumeWithHalo() const {
  const int all_halo = halo_under + halo_over;
  return (a_dim + all_halo) * (b_dim + all_halo) * (c_dim + all_halo);
}

//-------------------------------------------------------------------------------------------------
int VolumePartition::getSystemIndex() const {
  return system_index;
}

//-------------------------------------------------------------------------------------------------
VolumePartition VolumePartition::split(const int a_dim_update) {
  const int cleaved_adim = a_dim - a_dim_update;
  a_dim = a_dim_update;
  return VolumePartition(cleaved_adim, b_dim, c_dim, a_orig + a_dim_update, b_orig, c_orig,
                         halo_under, halo_over, system_index);
}
  
//-------------------------------------------------------------------------------------------------
Brickwork::Brickwork(const std::vector<int3> &system_dimensions_in, const int a_span_max_in,
                     const int bc_cross_section_max_in, const int halo_under_in,
                     const int halo_over_in, const int max_nonhalo_volume_in,
                     const int target_multiple, const std::vector<int> &preferred_a_lengths,
                     const std::vector<int> &discouraged_a_lengths) :
    a_span_max{a_span_max_in}, bc_cross_section_max{bc_cross_section_max_in},
    halo_under{halo_under_in}, halo_over{halo_over_in}, max_nonhalo_volume{max_nonhalo_volume_in},
    system_dimensions{system_dimensions_in}, bricks{}
{
  // Check that the maximum non-halo volume can accommodate the halo-inclusive cross section.
  const int all_halo = halo_under + halo_over;
  if ((max_nonhalo_volume + all_halo) * (1 + all_halo) < bc_cross_section_max) {
    rtErr("A maximum non-halo volume of " + std::to_string(max_nonhalo_volume) + " is "
          "incompatible with a maximum (halo-inclusive) cross section of " +
          std::to_string(bc_cross_section_max) + ".", "Brickwork");
  }
  
  // Validate the halos in light of the maximum cross-section
  if ((all_halo + 1) * (all_halo + 1) > bc_cross_section_max) {
    rtErr("Halo spacings of " + std::to_string(halo_under) + " (under) and " +
          std::to_string(halo_over) + " (over) are too large for a maximum cross sectional area "
          "of " + std::to_string(bc_cross_section_max) + ".", "Brickwork");
  }
  subdivide(target_multiple, preferred_a_lengths, discouraged_a_lengths);
}

//-------------------------------------------------------------------------------------------------
int Brickwork::getBrickCount() const {
  return bricks.size();
}

//-------------------------------------------------------------------------------------------------
int Brickwork::getBrickCount(const int system_index) const {
  validateSystemIndex(system_index, "getBrickCount");
  const int nbricks = bricks.size();
  int result;
  for (int i = 0; i < nbricks; i++) {
    result += (bricks[i].getSystemIndex() == system_index);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int Brickwork::getSystemCount() const {
  return system_dimensions.size();
}

//-------------------------------------------------------------------------------------------------
double2 Brickwork::getAverageBrickLength(const UnitCellAxis dim) const {
  double2 result = { 0.0, 0.0 };
  const size_t nbricks = bricks.size();
  double running_sq = 0.0;
  for (size_t i = 0; i < nbricks; i++) {
    const double len_i = bricks[i].getLength(dim);
    result.x += len_i;
    running_sq += len_i * len_i;
  }
  const double dn = nbricks;
  result.y = sqrt((dn * running_sq) - (result.x * result.x)) / sqrt(dn * (dn - 1.0));
  result.x /= dn;
  return result;
}

//-------------------------------------------------------------------------------------------------
double2 Brickwork::getAverageBrickLength(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return getAverageBrickLength(UnitCellAxis::A);
  case CartesianDimension::Y:
    return getAverageBrickLength(UnitCellAxis::B);
  case CartesianDimension::Z:
    return getAverageBrickLength(UnitCellAxis::C);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int3 Brickwork::getBrickOrigin(const int brick_index) const {
  validateBrickIndex(brick_index, "getBrickOrigin");
  return bricks[brick_index].getOrigin();
}

//-------------------------------------------------------------------------------------------------
int3 Brickwork::getBrickOriginWithHalo(const int brick_index) const {
  validateBrickIndex(brick_index, "getBrickLengthsWithHalo");
  int3 result = bricks[brick_index].getOrigin();
  result.x -= halo_under;
  result.y -= halo_under;
  result.z -= halo_under;
  const size_t system_index = bricks[brick_index].getSystemIndex();
  result.x += (result.x < 0) * system_dimensions[system_index].x;
  result.y += (result.y < 0) * system_dimensions[system_index].y;
  result.z += (result.z < 0) * system_dimensions[system_index].z;
  return result;
}

//-------------------------------------------------------------------------------------------------
int3 Brickwork::getBrickLengths(const int brick_index) const {
  validateBrickIndex(brick_index, "getBrickLengths");
  return bricks[brick_index].getLength();
}

//-------------------------------------------------------------------------------------------------
int3 Brickwork::getBrickLengthsWithHalo(const int brick_index) const {
  validateBrickIndex(brick_index, "getBrickLengthsWithHalo");
  return bricks[brick_index].getLengthPlusHalo();
}

//-------------------------------------------------------------------------------------------------
int Brickwork::getSystemMembership(const int brick_index) const {
  validateBrickIndex(brick_index, "getSystemMembership");
  return bricks[brick_index].getSystemIndex();
}

//-------------------------------------------------------------------------------------------------
int3 Brickwork::getSystemDimensions(const int system_index) const {
  validateSystemIndex(system_index, "getSystemDimensions");
  return system_dimensions[system_index];
}

//-------------------------------------------------------------------------------------------------
void Brickwork::subdivide(const int target_work_unit_multiple,
                          const std::vector<int> &preferred_a_lengths,
                          const std::vector<int> &discouraged_a_lengths) {
  const int nsys = system_dimensions.size();

  // Lay out a search range to find the optimal spread of work units.  Determine the largest cross
  // section, given the constraints of the halo.
  const int all_halo = halo_under + halo_over;
  int best_bdim_wh, best_cdim_wh;
  const int low_search_range = 1 + all_halo;
  const int hgh_search_range = bc_cross_section_max / (1 + all_halo);

  // Initialize the number of work units to be that obtained by putting each cell of the synthesis
  // of systems into its own work unit.
  int best_nwu = 0;
  for (int i = 0; i < nsys; i++) {
    const int3 idims = system_dimensions[i];
    best_nwu += idims.x * idims.y * idims.z;
  }

  // Scan a range of different possibilities for the basic work unit.  The B/C cross section limit
  // includes the halo region.
  int best_bdim = 1 + all_halo;
  int best_cdim = 1 + all_halo;
  for (int i = low_search_range; i <= hgh_search_range; i++) {
    const int max_cdim = bc_cross_section_max / i;
    for (int j = 1 + all_halo; j < max_cdim; j++) {
      const int nwu = countMinimalWorkUnits(i - all_halo, j - all_halo);
      if (nwu < best_nwu) {
        best_nwu = nwu;
        best_bdim = i;
        best_cdim = j;
      }
    }
  }
  best_bdim -= all_halo;
  best_cdim -= all_halo;
  
  // Create work units based on the most efficient layout
  bricks.reserve(best_nwu);
  bricks.resize(0);
  std::vector<int> a_lengths(16);
  for (int i = 0; i < nsys; i++) {
    const int3 idims = system_dimensions[i];

    const int nlrg_b = idims.y / best_bdim;
    const int nlrg_c = idims.z / best_cdim;
    const int remainder_b = idims.y - (nlrg_b * best_bdim);
    const int remainder_c = idims.z - (nlrg_c * best_cdim);

    // Lay out the largest possible work units for this layer with respect to the A axis
    int ta_span = std::min(max_nonhalo_volume / (best_bdim * best_cdim),
                           a_span_max - all_halo);
    if (ta_span == 0) {
      rtErr("A work unit non-halo region with non-halo cross section " +
            std::to_string(best_bdim) + " x " + std::to_string(best_cdim) + " and maximum "
            "non-halo volume " + std::to_string(max_nonhalo_volume) + " results in a unit cell A "
            "axis dimension of zero.", "Brickwork", "subdivide");
    }

    // Adjust the best length downwards to meet the highest possible preferred length.  Determine
    // a suitable combination of lengths that are all in the "preferred" category, or otherwise
    // avoid the "discouraged" category.  Obtain a list of lengths for work units of this
    // cross section.
    partition(idims.x, ta_span, preferred_a_lengths, discouraged_a_lengths, &a_lengths);
    const size_t na_lim = a_lengths.size();
    for (int nb = 0; nb < nlrg_b; nb++) {
      for (int nc = 0; nc < nlrg_c; nc++) {
        int aprog = 0;
        for (int na = 0; na < na_lim; na++) {
          bricks.emplace_back(a_lengths[na], best_bdim, best_cdim, aprog, nb * best_bdim,
                              nc * best_cdim, halo_under, halo_over, i);
          aprog += a_lengths[na];
        }
      }
    }

    // Fill out work units to complete the far right side, moving along the C axis
    if (remainder_b > 0) {
      const int height = (bc_cross_section_max / (remainder_b + all_halo)) - all_halo;
      const int va_span = std::min(max_nonhalo_volume / (height * remainder_b),
                                   a_span_max - all_halo);
      if (va_span == 0) {
        rtErr("A work unit non-halo region with non-halo cross section " +
              std::to_string(remainder_b) + " x " + std::to_string(height) + " and maximum "
              "non-halo volume " + std::to_string(max_nonhalo_volume) + " results in a unit cell "
              "A axis dimension of zero.  This occurred when back-filling the unit cell B axis.",
              "Brickwork", "subdivide");
      }
      partition(idims.x, va_span, preferred_a_lengths, discouraged_a_lengths, &a_lengths);
      const size_t na_lim = a_lengths.size();
      for (int nc = 0; nc < idims.z; nc += height) {
        const int choose_cdim = (nc + height <= idims.z) ? height : idims.z - nc;
        int aprog = 0;
        for (int na = 0; na < na_lim; na++) {
          bricks.emplace_back(a_lengths[na], remainder_b, choose_cdim, aprog,
                              idims.y - remainder_b, nc, halo_under, halo_over, i);
          aprog += a_lengths[na];
        }
      }
    }

    // Create out work units to fill out the top, moving along the B axis
    if (remainder_c > 0) {
      const int total_b_span = nlrg_b * best_bdim;
      const int width = (bc_cross_section_max / (remainder_c + all_halo)) - all_halo;
      const int va_span = std::min(max_nonhalo_volume / (width * remainder_c),
                                   a_span_max - all_halo);
      if (va_span == 0) {
        rtErr("A work unit non-halo region with non-halo cross section " +
              std::to_string(width) + " x " + std::to_string(remainder_c) + " and maximum "
              "non-halo volume " + std::to_string(max_nonhalo_volume) + " results in a unit cell "
              "A axis dimension of zero.  This occurred when back-filling the unit cell C axis.",
              "Brickwork", "subdivide");
      }
      partition(idims.x, va_span, preferred_a_lengths, discouraged_a_lengths, &a_lengths);
      const size_t na_lim = a_lengths.size();
      for (int nb = 0; nb < total_b_span; nb += width) {
        const int choose_bdim = (nb + width <= total_b_span) ? width : total_b_span - nb;
        int aprog = 0;
        for (int na = 0; na < na_lim; na++) {
          bricks.emplace_back(a_lengths[na], choose_bdim, remainder_c, aprog, nb,
                              idims.z - remainder_c, halo_under, halo_over, i);
          aprog += a_lengths[na];
        }
      }
    }
  }

  // Order the work units
  std::sort(bricks.begin(), bricks.end(),
            [](const VolumePartition &a,
               const VolumePartition &b) { return (a.getVolumeWithHalo() >
                                                   b.getVolumeWithHalo()); });
  
  // If the number of work units if less than four times the target and far from a multiple of it,
  // split the largest work units until the number can be brought in line with the goal.
  const size_t next_bp = roundUp(bricks.size(), static_cast<size_t>(target_work_unit_multiple));
  if (next_bp >= 4 * target_work_unit_multiple) {
    return;
  }
  else {
    bricks.reserve(next_bp);
  }
  bool more_bricks = true;
  while (bricks.size() < next_bp && more_bricks) {
    const int nbricks = bricks.size();
    int i = 0;
    while (i < nbricks && bricks.size() < next_bp) {
      const int ilen = bricks[i].getLength(UnitCellAxis::A);
      if (ilen >= 2) {
        VolumePartition nvp = bricks[i].split(ilen / 2);
        bricks.push_back(nvp);
      }
      i++;
    }
    more_bricks = (bricks.size() >= static_cast<size_t>(nbricks));
    
    // Order the work units, before exiting and before possibly trying again
    std::sort(bricks.begin(), bricks.end(),
              [](const VolumePartition &a,
                 const VolumePartition &b) { return (a.getVolumeWithHalo() >
                                                     b.getVolumeWithHalo()); });
  }
}

//-------------------------------------------------------------------------------------------------
void Brickwork::validateSystemIndex(const int system_index, const char* caller) const {
  if (system_index < 0 || system_index >= system_dimensions.size()) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for a collection of " +
          std::to_string(system_dimensions.size()) + ".", "Brickwork", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void Brickwork::validateBrickIndex(const int index, const char* caller) const {
  if (index < 0 || index >= bricks.size()) {
    rtErr("Brick index " + std::to_string(index) + " is invalid.", "Brickwork", caller);
  }
}

//-------------------------------------------------------------------------------------------------
int Brickwork::countMinimalWorkUnits(const int trimmed_bdim, const int trimmed_cdim) const {
  int result = 0;
  const int all_halo = halo_under + halo_over;
  const int nsys = system_dimensions.size();
  for (int i = 0; i < nsys; i++) {
    const int3 idims = system_dimensions[i];

    // Lay out work units with the most efficient cross sections first.
    const int nlrg_b = idims.y / trimmed_bdim;
    const int nlrg_c = idims.z / trimmed_cdim;
    const int remainder_b = idims.y - (nlrg_b * trimmed_bdim);
    const int remainder_c = idims.z - (nlrg_c * trimmed_cdim);
    int i_result = nlrg_b * nlrg_c;
    
    // Add work units to fill in the final part of the B-axis range
    if (remainder_b > 0) {
      const int height = (bc_cross_section_max / (remainder_b + all_halo)) - all_halo;
      i_result += (idims.z + height - 1) / height;
    }

    // Add work units to fill in the final part of the C-axis range
    if (remainder_c > 0) {
      const int total_b_span = nlrg_b * trimmed_bdim;
      const int width = (bc_cross_section_max / (remainder_c + all_halo)) - all_halo;
      i_result += (total_b_span + width - 1) / width;
    }
    
    // Assume that the bricks shear evenly across the BC plane
    result += ((idims.x + a_span_max - all_halo - 1) / (a_span_max - all_halo)) * i_result;
  }
  return result;
}

} // namespace synthesis
} // namespace stormm
