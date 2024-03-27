//-*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename Tmask> std::vector<NonbondedWorkUnit>
enumerateNonbondedWorkUnits(const Tmask &se, const int target_tile_count,
                            const int nbwu_estimate, const std::vector<int> &atom_counts,
                            const std::vector<int> &atom_offsets) {
  std::vector<NonbondedWorkUnit> result;
  result.reserve(nbwu_estimate);
  const int system_count = atom_counts.size();
  std::vector<int3> tile_list(target_tile_count);
  const int max_imports = small_block_max_atoms / tile_length;
  std::vector<int> import_list(max_imports);

  // The tile length must be less than the warp size (the intention is that it be half the
  // warp size), but the tiles for each system must not overlap in the list of atoms from the
  // topology or coordinate syntheses.
  int total_tile_lengths = 0;
  std::vector<int> system_tile_starts(system_count);
  for (int i = 0; i < system_count; i++) {
    system_tile_starts[i] = total_tile_lengths;
    total_tile_lengths += (atom_counts[i] + tile_length - 1) / tile_length;
  }

  // Assemble work units
  int current_tile_count = 0;
  int import_count = 0;
  std::vector<int> import_coverage(total_tile_lengths, 0);
  for (int sysid = 0; sysid < system_count; sysid++) {
    const int system_tile_lengths = (atom_counts[sysid] + tile_length - 1) / tile_length;

    // The strategy is to list tiles in the typical order, keeping track of import requirements.
    // Get as close as possible to the target work unit size, keeping to multiples of eight if at
    // all possible, then commit work units to the list as they reach the target size or it
    // becomes impossible to grow the work unit without exceeding the import limits.
    const int ordi_stride = (target_tile_count ==   tiny_nbwu_tiles) ? 2 :
                            (target_tile_count ==  small_nbwu_tiles) ? 4 :
                            (target_tile_count == medium_nbwu_tiles) ? 6 :
                            (target_tile_count ==  large_nbwu_tiles) ? 8 :
                            tile_lengths_per_supertile;
    int ordi_llim = 0;
    int ordi_hlim = std::min(ordi_llim + ordi_stride, system_tile_lengths);
    int absc_increment = 1;
    while (ordi_llim < system_tile_lengths) {
      int istart, ilimit;
      if (absc_increment == 1) {
        istart = 0;
        ilimit = system_tile_lengths;
      }
      else {
        istart = system_tile_lengths - 1;
        ilimit = -1;
      }
      for (int i = istart; i != ilimit; i += absc_increment) {
        for (int j = ordi_llim; j < ordi_hlim; j++) {
          if (j > i) {
            continue;
          }

          // Test whether this new tile can be added within the limits of the work unit.
          if (addTileToWorkUnitList(tile_list.data(), import_coverage.data(), system_tile_starts,
                                    &import_count, &current_tile_count, i, j, sysid)) {
            if (current_tile_count == target_tile_count) {
              
              // This list of tiles has reached its target size and should be converted into a
              // work unit.
              result.emplace_back(se, tile_list);
              for (int k = 0 ; k < current_tile_count; k++) {
                import_coverage[tile_list[k].x + system_tile_starts[tile_list[k].z]] = 0;
                import_coverage[tile_list[k].y + system_tile_starts[tile_list[k].z]] = 0;
              }
              import_count = 0;
              current_tile_count = 0;
            }
          }
          else {
            
            // Backtrack to the most recent multiple of the target batch size.
            const int nbatch = current_tile_count / small_block_tile_width;
            const int fallback_size = nbatch * small_block_tile_width;
            std::vector<int3> tmp_tile_list(fallback_size);
            for (int k = 0; k < fallback_size; k++) {
              tmp_tile_list[k] = tile_list[k];
            }
            result.emplace_back(se, tmp_tile_list);
            for (int k = 0 ; k < current_tile_count; k++) {
              import_coverage[tile_list[k].x + system_tile_starts[tile_list[k].z]] = 0;
              import_coverage[tile_list[k].y + system_tile_starts[tile_list[k].z]] = 0;
            }
            for (int k = fallback_size; k < current_tile_count; k++) {
              tmp_tile_list[k - fallback_size] = tile_list[k];
            }
            const int tmp_tile_count = current_tile_count - fallback_size;
            current_tile_count = 0;
            import_count = 0;
            
            // Add the tiles that could not be included in the work unit just commited to the
            // start of a new work unit.  There will be room, although the new work unit may
            // need to be committed sooner than usual.
            for (int k = 0; k < tmp_tile_count; k++) {
              addTileToWorkUnitList(tile_list.data(), import_coverage.data(), system_tile_starts,
                                    &import_count, &current_tile_count, tmp_tile_list[k].x,
                                    tmp_tile_list[k].y, tmp_tile_list[k].z);
            }
            
            // Add the tile at hand, the one which was too much to add to the previous work unit
            // and thus triggered the backtracking.
            addTileToWorkUnitList(tile_list.data(), import_coverage.data(), system_tile_starts,
                                  &import_count, &current_tile_count, i, j, sysid);
          }
        }
      }
      absc_increment *= -1;
      ordi_llim = ordi_hlim;
      ordi_hlim = std::min(ordi_llim + ordi_stride, system_tile_lengths);
    }
  }
  
  // Commit the final tile group, then return the result
  tile_list.resize(current_tile_count);
  result.emplace_back(se, tile_list);

  // Check that no pollution was left the in the tile coverage array--this would indicate possible
  // problems in the construction and could lead to work units importing too many tiles.
  for (int i = 0; i < current_tile_count; i++) {
    import_coverage[tile_list[i].x + system_tile_starts[tile_list[i].z]] = 0;
    import_coverage[tile_list[i].y + system_tile_starts[tile_list[i].z]] = 0;    
  }
  const int stray_imports = sum<int>(import_coverage);
  if (stray_imports != 0) {
    rtErr("A total of " + std::to_string(stray_imports) + " tiles were not cleaned in the import "
          "coverage array.", "enumerateNonbondedWorkUnits");
  }

  // Set initialization masks for all work units, then return the result
  if (target_tile_count == huge_nbwu_tiles) {
    setPropertyInitializationMasks(&result, NbwuKind::SUPERTILES);
  }
  else {
    setPropertyInitializationMasks(&result, NbwuKind::TILE_GROUPS);
  }
  return result;
}

} // namespace synthesis
} // namespace stormm
