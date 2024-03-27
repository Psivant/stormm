#include <algorithm>
#include "copyright.h"
#include "static_exclusionmask.h"

namespace stormm {
namespace energy {

using topology::NonbondedKit;

//-------------------------------------------------------------------------------------------------
StaticExclusionMaskReader::StaticExclusionMaskReader(const int natom_in,
                                                     const int supertile_stride_count_in,
                                                     const int* supertile_map_idx_in,
                                                     const int* tile_map_idx_in,
                                                     const uint* mask_data_in) :
    natom{natom_in}, supertile_stride_count{supertile_stride_count_in},
    supertile_map_idx{supertile_map_idx_in}, tile_map_idx{tile_map_idx_in}, mask_data{mask_data_in}
{}

//-------------------------------------------------------------------------------------------------
StaticExclusionMask::StaticExclusionMask(const AtomGraph *ag_in) :
    atom_count{(ag_in == nullptr) ? 0 : ag_in->getAtomCount()},
    supertile_stride_count{(atom_count + supertile_length - 1) / supertile_length},
    tile_stride_count{(atom_count + tile_length - 1) / tile_length},
    unique_supertile_count{0},
    unique_tile_count{0},
    all_masks{0, "exclmask_data"},
    supertile_map_indices{0, "exclmask_supertiles"},
    tile_map_indices{0, "exclmask_tiles"}
{
  // Return immediately if the nullptr was fed in
  if (ag_in == nullptr) {
    return;
  }
  
  // Check the system's idea of an int
  if (sizeof(int) < 4LLU) {
    rtErr("This architecture stores an int in " + std::to_string(sizeof(int) * 8) + "-bit format, "
          "but a static, all-to-all exclusion mask requires that an int have at least 32 bits.",
          "StaticExclusionMask");
  }

  // Check the total number of atoms
  if (supertile_stride_count > max_supertile_count) {
    rtErr("This system is too large to store a static exclusion list.  Due to the need to store "
          "at least one int index for the interaction of every group of 256 atoms with any other, "
          "with int assumed to be 32 bit, the maximum number of atoms that can be indexed in "
          "this manner is about sqrt(2^31) * 256 ~ 46,300 * 256 = 11.8 million.",
          "StaticExclusionMask");
  }

  // Initialize the first (blank) tile: 32 unsigned integers set to zero, to which all tiles in
  // the blank supertile will index, indicating that there are no exclusions.
  std::vector<uint> tmp_masks(2 * tile_length, 0);
  std::vector<int> tmp_tile_map_indices(tiles_per_supertile, 0);
  int n_unique_supertiles = 1;
  int n_unique_tiles = 1;

  // Step over all supertiles with a list of the pair exclusions in hand
  const NonbondedKit<double> nbk = ag_in->getDoublePrecisionNonbondedKit();
  std::vector<int> tmp_supertile_map_indices(supertile_stride_count * supertile_stride_count);
  for (int sti = 0; sti < supertile_stride_count; sti++) {
    const int isptl_start = sti * supertile_length;
    const int isptl_end = std::min((sti + 1) * supertile_length, atom_count);
    for (int stj = 0; stj < supertile_stride_count; stj++) {
      const int jsptl_start = stj * supertile_length;
      const int jsptl_end = std::min((stj + 1) * supertile_length, atom_count);

      // Identify any exclusions within this tile.
      bool excl_found = ((isptl_end - isptl_start) % tile_length > 0 ||
                         (jsptl_end - jsptl_start) % tile_length > 0 ||
                         isptl_start == jsptl_start);
      for (int i = isptl_start; i < isptl_end; i++) {
        
        // Check 1:1 exclusions
        for (int j = nbk.nb11_bounds[i]; j < nbk.nb11_bounds[i + 1]; j++) {
          excl_found = (excl_found || (nbk.nb11x[j] >= jsptl_start && nbk.nb11x[j] < jsptl_end));
        }

        // Check 1:2 exclusions
        for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
          excl_found = (excl_found || (nbk.nb12x[j] >= jsptl_start && nbk.nb12x[j] < jsptl_end));
        }

        // Check 1:3 exclusions
        for (int j = nbk.nb13_bounds[i]; j < nbk.nb13_bounds[i + 1]; j++) {
          excl_found = (excl_found || (nbk.nb13x[j] >= jsptl_start && nbk.nb13x[j] < jsptl_end));
        }

        // Check 1:4 exclusions
        for (int j = nbk.nb14_bounds[i]; j < nbk.nb14_bounds[i + 1]; j++) {
          excl_found = (excl_found || (nbk.nb14x[j] >= jsptl_start && nbk.nb14x[j] < jsptl_end));
        }

        // Break if an exclusion has been identified
        if (excl_found) {
          break;
        }
      }

      // If an exclusion was not found, this supertile is blank and indexes to the blank map.
      if (excl_found == false) {
        tmp_supertile_map_indices[(stj * supertile_stride_count) + sti] = 0;
        continue;
      }

      // Otherwise, this supertile has some exclusion somewhere within it, which means that
      // at least 256 indices to individual tile exclusion maps must be supplied.
      tmp_supertile_map_indices[(stj * supertile_stride_count) + sti] = n_unique_supertiles * 256;
      std::vector<int> this_supertile_map(tiles_per_supertile, 0);

      // Within this supertile, count the number of unique tiles and log them in the master
      // store of exclusion masks.  By default, every tile's map points to the blank exclusion
      // mask, tile 0.
      const int ni_tiles = (isptl_end - isptl_start + tile_length - 1) / tile_length;
      const int nj_tiles = (jsptl_end - jsptl_start + tile_length - 1) / tile_length;
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int itl_start = isptl_start + (ti * tile_length);
        const int itl_end = std::min(itl_start + tile_length, isptl_end);
        for (int tj = 0; tj < nj_tiles; tj++) {
          const int jtl_start = jsptl_start + (tj * tile_length);
          const int jtl_end = std::min(jtl_start + tile_length, jsptl_end);
          std::vector<uint> mask_buffer(2 * tile_length, 0);
          excl_found = (itl_end - itl_start < tile_length ||
                        jtl_end - jtl_start < tile_length || itl_start == jtl_start);
          
          // If the tile is incomplete due to running off the end of the number of system atoms,
          // the extra interactions must be listed as exclusions for the GPU code to understand
          // not to count them.
          uint jmask = 0U;
          for (int i = itl_end - itl_start; i < tile_length; i++) {
            jmask |= (0x1 << i);
          }
          for (int i = itl_end - itl_start; i < tile_length; i++) {
            mask_buffer[i] = 0xffffffff;
          }
          for (int i = tile_length; i < 2 * tile_length; i++) {
            mask_buffer[i] |= jmask;
          }
          int imask = 0U;
          for (int j = jtl_end - jtl_start; j < tile_length; j++) {
            imask |= (0x1 << j);
          }
          for (int j = jtl_end - jtl_start; j < tile_length; j++) {
            mask_buffer[tile_length + j] = 0xffffffff;
          }
          for (int j = 0; j < tile_length; j++) {
            mask_buffer[j] |= imask;
          }

          // Diagonal tiles must exclude all atoms for which j >= i
          if (itl_start == jtl_start) {
            uint diag_mask = 0xffff;
            for (int i = 0; i < tile_length; i++) {
              mask_buffer[i] |= diag_mask;
              diag_mask ^= (0x1 << i);
            }
            diag_mask = 0U;
            for (int i = tile_length; i < 2 * tile_length; i++) {
              diag_mask |= (0x1 << (i - tile_length));
              mask_buffer[i] |= diag_mask;
            }
          }

          // Scan all atoms along the abscissa to see if their known 1:1, 1:2, 1:3, or 1:4
          // exclusions cover any atoms along the ordinate for this tile.
          for (int i = itl_start; i < itl_end; i++) {

            // The self-exclusion is not handled in the context of this list (it would introduce
            // more tiles in some situations, which are not necessary).  Code that special-cases
            // tiles with self interactions will be necessary just to prevent double-counting,
            // and this code can incorporate prohibitions against self interactions.

            // Check 1:1 exclusions
            for (int j = nbk.nb11_bounds[i]; j < nbk.nb11_bounds[i + 1]; j++) {
              const int jatom = nbk.nb11x[j];
              if (jatom >= jtl_start && jatom < jtl_end) {
                mask_buffer[i - itl_start] |= (0x1 << (jatom - jtl_start));
                mask_buffer[tile_length + jatom - jtl_start] |= (0x1 << (i - itl_start));
                excl_found = true;
              }
            }

            // Check 1:2 exclusions
            for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
              const int jatom = nbk.nb12x[j];
              if (jatom >= jtl_start && jatom < jtl_end) {
                mask_buffer[i - itl_start] |= (0x1 << (jatom - jtl_start));
                mask_buffer[tile_length + jatom - jtl_start] |= (0x1 << (i - itl_start));
                excl_found = true;
              }
            }

            // Check 1:3 exclusions
            for (int j = nbk.nb13_bounds[i]; j < nbk.nb13_bounds[i + 1]; j++) {
              const int jatom = nbk.nb13x[j];
              if (jatom >= jtl_start && jatom < jtl_end) {
                mask_buffer[i - itl_start] |= (0x1 << (jatom - jtl_start));
                mask_buffer[tile_length + jatom - jtl_start] |= (0x1 << (i - itl_start));
                excl_found = true;
              }
            }

            // Check 1:4 exclusions
            for (int j = nbk.nb14_bounds[i]; j < nbk.nb14_bounds[i + 1]; j++) {
              const int jatom = nbk.nb14x[j];
              if (jatom >= jtl_start && jatom < jtl_end) {
                mask_buffer[i - itl_start] |= (0x1 << (jatom - jtl_start));
                mask_buffer[tile_length + jatom - jtl_start] |= (0x1 << (i - itl_start));
                excl_found = true;
              }
            }
          }

          // Log this tile mask if an exclusion was found
          if (excl_found) {
            tmp_masks.insert(tmp_masks.end(), mask_buffer.begin(), mask_buffer.end());
            this_supertile_map[(tj * tile_lengths_per_supertile) + ti] = n_unique_tiles *
                                                                         2 * tile_length;
            n_unique_tiles++;
            if (n_unique_tiles > max_unique_tiles) {
              rtErr("The maximum number of unique 16 x 16 atom tile masks has been exceeded.  A "
                    "system may create no more than (2^26) - 1 unique tile masks containing "
                    "exclusions, about 67 million or enough to accommodate a system of 11.8 "
                    "million atoms (the maximum system size) with every group of 16 atoms "
                    "sharing exclusions with up to 90 others.", "StaticExclusionMask");
            }
          }
        }
      }

      // Store the present supertile's map in the stack.  This will add a new set of 256
      // indices into the growing list of 16 x 16 atom tile masks.  Most of the supertile's
      // 16 x 16 atom tiles will still index to the zero position, as they will still have
      // no exclusions.  The process for an HPC warp crunching through the non-bonded loop
      // is, then, four-fold to obtain its relevant exclusion mask:
      //   1.) Get the tile number to work on next (atomicAdd 1 to some counter and take the
      //       result returned by the atomic operation as the current tile assignment)
      //   2.) Compute, based on the atom count, what supertile this tile assignment comes from
      //   3.) Read the supertile's mapping index from the StaticExclusionMask's member variable
      //       supertile_map_indices.  If that index is zero, then the supertile has no mask with
      //       any exclusions in it (skip step 4).
      //   4.) If the supertile's map index is nonzero, compute the specific tile index from within
      //       the supertile and read that element within the tile_map_indices member variable.
      //       If that result is zero, then there are no exclusions.  If that result (let it be
      //       called tile_mask_index) is nonzero, read
      //       all_masks[32 * tile_mask_index : 32 * (tile_mask_index + 1) as the exclusion mask
      //       for the current tile assignment.
      tmp_tile_map_indices.insert(tmp_tile_map_indices.end(), this_supertile_map.begin(),
                                  this_supertile_map.end());
      n_unique_supertiles++;
    }
  }

  // Allocate memory for the object (although the supertile_map_indices array could be allocated
  // during inline initialization, that is done here to first go through a check on the overall
  // system size).
  unique_supertile_count = n_unique_supertiles;
  unique_tile_count = n_unique_tiles;
  all_masks.resize(unique_tile_count * 2 * tile_length);
  supertile_map_indices.resize(supertile_stride_count * supertile_stride_count);
  tile_map_indices.resize(unique_supertile_count * tiles_per_supertile);

  // Load the object
  all_masks.putHost(tmp_masks);
  supertile_map_indices.putHost(tmp_supertile_map_indices);
  tile_map_indices.putHost(tmp_tile_map_indices);
  ag_pointer = ag_in;
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMask::StaticExclusionMask(const AtomGraph &ag_in) :
    StaticExclusionMask(ag_in.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMask::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMask::getSuperTileStrideCount() const {
  return supertile_stride_count;
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMask::getTileStrideCount() const {
  return tile_stride_count;
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMask::getUniqueSuperTileCount() const {
  return unique_supertile_count;
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMask::getUniqueTileCount() const {
  return unique_tile_count;
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMask::getUniqueTileCount(const int supertile_i_index,
                                            const int supertile_j_index) const {

  // Get the super-tile index from the master table of all super-tiles
  const int sptl_index = (supertile_stride_count * supertile_j_index) + supertile_i_index;
  const int sptl_map_index = supertile_map_indices.readHost(sptl_index);

  // The first index is the special case of zero exclusions in any tile
  if (sptl_map_index == 0) {
    return 0;
  }
  else {
    const int offset = tiles_per_supertile * sptl_map_index;
    const int* tilemap_ptr = tile_map_indices.data();
    int result = 0;
    for (int i = 0; i < tiles_per_supertile; i++) {
      result += (tilemap_ptr[i] > 0);
    }
    return result;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<uint>
StaticExclusionMask::getTileExclusions(const int supertile_i_index, const int supertile_j_index,
                                       const int tile_i_index, const int tile_j_index) const {
  const int sptl_index = (supertile_stride_count * supertile_j_index) + supertile_i_index;
  const int sptl_map_index = supertile_map_indices.readHost(sptl_index);
  const int tl_index = (tile_j_index * tile_lengths_per_supertile) + tile_i_index;
  const int tl_map_index = tile_map_indices.readHost(sptl_map_index + tl_index);
  return all_masks.readHost(tl_map_index, 2 * tile_length);
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* StaticExclusionMask::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const StaticExclusionMaskReader StaticExclusionMask::data(const HybridTargetLevel tier) const {
  return StaticExclusionMaskReader(atom_count, supertile_stride_count,
                                   supertile_map_indices.data(tier), tile_map_indices.data(tier),
                                   all_masks.data(tier));
}

//-------------------------------------------------------------------------------------------------
const StaticExclusionMask* StaticExclusionMask::getSelfPointer() const {
  return this;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void StaticExclusionMask::upload() {
  all_masks.upload();
  supertile_map_indices.upload();
  tile_map_indices.upload();
}

//-------------------------------------------------------------------------------------------------
void StaticExclusionMask::download() {
  all_masks.download();
  supertile_map_indices.download();
  tile_map_indices.download();
}
#endif

//-------------------------------------------------------------------------------------------------
bool StaticExclusionMask::testExclusion(const int atom_i, const int atom_j) const {

  // Get the supertile map index
  const int sptl_i = atom_i / supertile_length;
  const int sptl_j = atom_j / supertile_length;
  const int sptl_index = (supertile_stride_count * sptl_j) + sptl_i;
  const int sptl_map_index = supertile_map_indices.readHost(sptl_index);
  if (sptl_map_index == 0) {
    return false;
  }
  const int tl_i = (atom_i - (supertile_length * sptl_i)) / tile_length;
  const int tl_j = (atom_j - (supertile_length * sptl_j)) / tile_length;
  const int tl_index = (tl_j * tile_lengths_per_supertile) + tl_i;
  const int tl_map_index = tile_map_indices.readHost(sptl_map_index + tl_index);
  if (tl_map_index == 0) {
    return false;
  }
  const int atom_i_offset = atom_i - (supertile_length * sptl_i) - (tile_length * tl_i);
  const int atom_j_offset = atom_j - (supertile_length * sptl_j) - (tile_length * tl_j);
  const uint imask = all_masks.readHost(tl_map_index + atom_i_offset);
  return ((imask >> atom_j_offset) & 0x1);
}

} // namespace energy
} // namespace stormm
