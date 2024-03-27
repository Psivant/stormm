#include "copyright.h"
#include "Constants/scaling.h"
#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "nonbonded_workunit.h"

namespace stormm {
namespace synthesis {

using constants::kilo;
using constants::mega;
using energy::StaticExclusionMaskReader;
using energy::supertile_length;
using energy::tile_length;
using energy::tile_lengths_per_supertile;
using stmath::accumulateBitmask;
using stmath::findBin;
using stmath::reduceUniqueValues;
using stmath::roundUp;
  
//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMask &se,
                                     const std::vector<int3> &tile_list) :
    tile_count{static_cast<int>(tile_list.size())},
    import_count{0},
    kind{NbwuKind::TILE_GROUPS},
    score{0}, init_accumulation{0}, refresh_atom_start{0}, refresh_code{0},
    imports{std::vector<int>(small_block_max_imports)},
    import_system_indices{std::vector<int>(small_block_max_imports, 0)},
    import_size_keys{std::vector<int>(small_block_max_imports / 4, 0)},
    tile_instructions{std::vector<uint2>(tile_count)}
{
  // In this form of the constructor, one system is present with a single StaticExclusionMask.
  // Get the total number of atoms and inflate the tile starting positions to obtain the atom
  // import limits.
  const StaticExclusionMaskReader ser = se.data();
  const int n_supert = se.getSuperTileStrideCount();
  std::vector<int> tmp_imports(2 * tile_count);
  for (int i = 0; i < tile_count; i++) {
    if (tile_list[i].x < 0 || tile_list[i].y < 0 || tile_list[i].x * tile_length >= ser.natom ||
        tile_list[i].y * tile_length >= ser.natom) {
      rtErr("A tile with lower left corner (" + std::to_string(tile_list[i].x * tile_length) +
            ", " + std::to_string(tile_list[i].y * tile_length) + ") cannot be part of a system "
            "with " + std::to_string(ser.natom) + " atoms.", "NonbondedWorkUnit");
    }
    const int sti = tile_list[i].x / tile_lengths_per_supertile;
    const int stj = tile_list[i].y / tile_lengths_per_supertile;
    const int stij_map_index = ser.supertile_map_idx[(stj * n_supert) + sti];
    const int ti = tile_list[i].x - (sti * tile_lengths_per_supertile);
    const int tj = tile_list[i].y - (stj * tile_lengths_per_supertile);
    const int tij_map_index = ser.tile_map_idx[stij_map_index +
                                               (tj * tile_lengths_per_supertile) + ti];
    tile_instructions[i].y = tij_map_index;
    tmp_imports[(2 * i)    ] = tile_list[i].x;
    tmp_imports[(2 * i) + 1] = tile_list[i].y;
  }
  reduceUniqueValues(&tmp_imports);
  import_count = tmp_imports.size();
  if (import_count >= small_block_max_imports) {
    rtErr("Too many imports (" + std::to_string(import_count) + ", with " +
          std::to_string(small_block_max_imports) + " being the maximum) would be required in a "
          "non-bonded work unit for topology " + se.getTopologyPointer()->getFileName() + ".",
          "NonbondedWorkUnit");
  }
  for (int i = 0; i < import_count; i++) {
    tmp_imports[i] *= tile_length;
    imports[i] = tmp_imports[i];
  }
  for (int i = import_count; i < small_block_max_imports; i++) {
    imports[i] = -1;
  }
  
  // Compute the workunit effort score, an estimate of how long one work unit will take relative
  // to others.  Work units with the highest effort will be placed first, to backfill idle
  // processes with shorter work units.
  score = import_count + (tile_count * 8);
  
  // The limited number of imports, and the manipulations that need to happen to assemble the
  // work unit instructions, likely makes it preferable to just search them all rather than
  // calling a binary search function.  Each instruction's y member has already been filled out.
  for (int i = 0; i < import_count; i++) {
    const int import_idx = imports[i] / tile_length;
    const uint absc_mask = tile_length * i;
    const uint ordi_mask = (absc_mask << 16);
    for (int j = 0; j < tile_count; j++) {
      if (tile_list[j].x == import_idx) {
        tile_instructions[j].x |= absc_mask;
      }
      if (tile_list[j].y == import_idx) {
        tile_instructions[j].x |= ordi_mask;
      }
    }
  }

  // Fill out the size keys for each import, a series of bit-packed integers with one import
  // batch size (up to tile_length atoms, i.e. 16) per eight bits.
  if (constants::int_bit_count_int < 32) {
    rtErr("Descriptors for non-bonded work units require that signed integers be of size >= 32 "
          "bits.", "NonbondedWorkUnit");
  }
  for (int i = 0; i < import_count; i++) {
    const int key_idx = i / 4;
    const int key_pos = i - (key_idx * 4);
    int import_batch_size;
    if (imports[i] >= 0) {
      import_batch_size = std::min(ser.natom - imports[i], tile_length);
    }
    else {
      import_batch_size = 0;
    }
    import_size_keys[key_idx] |= (import_batch_size << (8 * key_pos));
  }
}

//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMask &se,
                                     const int abscissa_start, const int ordinate_start) :
    tile_count{0},
    import_count{0},
    kind{NbwuKind::SUPERTILES},
    score{0}, init_accumulation{0}, refresh_atom_start{0}, refresh_code{0},
    imports{},
    import_system_indices{std::vector<int>(1, 0)},
    import_size_keys{},
    tile_instructions{}
{
  imports.resize(3);
  imports[0] = abscissa_start;
  imports[1] = ordinate_start;
  
  // Compute the tile count and the work unit score
  const int atom_limit = se.getAtomCount();
  imports[2] = (abscissa_start + supertile_length > atom_limit) ? atom_limit - abscissa_start :
                                                                  supertile_length;
  imports[3] = (ordinate_start + supertile_length > atom_limit) ? atom_limit - ordinate_start :
                                                                  supertile_length;
  if (abscissa_start == ordinate_start) {
    if (abscissa_start < atom_limit - supertile_length &&
        ordinate_start < atom_limit - supertile_length) {
      tile_count = tile_lengths_per_supertile * (tile_lengths_per_supertile + 1) / 2;
      import_count = tile_lengths_per_supertile;
    }
    else {
      const int ntile_rem = (atom_limit - abscissa_start + tile_length - 1) / tile_length;  
      tile_count = ntile_rem * (ntile_rem + 1) / 2;
      import_count = ntile_rem;
    }
  }
  else {
    if (abscissa_start < atom_limit - supertile_length &&
        ordinate_start < atom_limit - supertile_length) {
      tile_count = tile_lengths_per_supertile * tile_lengths_per_supertile;
      import_count = 2 * tile_lengths_per_supertile;
    }
    else if (abscissa_start < atom_limit - supertile_length) {
      const int ntile_rem = (atom_limit - ordinate_start + tile_length - 1) / tile_length;
      tile_count = tile_lengths_per_supertile * ntile_rem;
      import_count = tile_lengths_per_supertile + ntile_rem;
    }
    else if (ordinate_start < atom_limit - supertile_length) {
      const int ntile_rem = (atom_limit - abscissa_start + tile_length - 1) / tile_length;
      tile_count = tile_lengths_per_supertile * ntile_rem;
      import_count = tile_lengths_per_supertile + ntile_rem;
    }
    else {
      rtErr("A trapezoidal tile should not exist in a supertile-based work unit.",
            "NonbondedWorkUnit");
    }
  }
  score = import_count + (tile_count * 8);
}

//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMaskSynthesis &se,
                                     const std::vector<int3> &tile_list) :
    tile_count{static_cast<int>(tile_list.size())},
    import_count{0},
    kind{NbwuKind::TILE_GROUPS},
    score{0}, init_accumulation{0}, refresh_atom_start{0}, refresh_code{0},
    imports{std::vector<int>(small_block_max_imports)},
    import_system_indices{std::vector<int>(small_block_max_imports)},
    import_size_keys{std::vector<int>(small_block_max_imports / 4, 0)},
    tile_instructions{std::vector<uint2>(tile_count)}
{
  // As before, loop over tiles to determine the array of imported atoms.  Work from the static
  // exclusion mask compilation rather than a solitary system's exclusions.  The imports will be
  // determined not just for tile numbers but also for system indices--the fourth tile edge's
  // worth of atoms in the first system is not the fourth tile edge's worth of atoms in the third
  // system.  Long long ints and standard sorting to the rescue.
  SeMaskSynthesisReader ser = se.data();
  std::vector<llint> tmp_imports(2 * tile_count);
  for (int i = 0; i < tile_count; i++) {

    // Get the atom count, and thus the supertile count, for the particular tile's system
    const int natom = ser.atom_counts[tile_list[i].z];
    const int n_supert = (natom + supertile_length - 1) / supertile_length;
    if (tile_list[i].x < 0 || tile_list[i].y < 0 || tile_list[i].x * tile_length >= natom ||
        tile_list[i].y * tile_length >= natom) {
      rtErr("A tile with lower left corner (" + std::to_string(tile_list[i].x * tile_length) +
            ", " + std::to_string(tile_list[i].y * tile_length) + ") cannot be part of a system "
            "(index " + std::to_string(tile_list[i].z) + ") with " + std::to_string(natom) +
            " atoms.", "NonbondedWorkUnit");
    }
    const int sti = tile_list[i].x / tile_lengths_per_supertile;
    const int stj = tile_list[i].y / tile_lengths_per_supertile;
    const int stm_bound = ser.supertile_map_bounds[tile_list[i].z];
    const int stij_map_index = ser.supertile_map_idx[stm_bound + (stj * n_supert) + sti];
    const int ti = tile_list[i].x - (sti * tile_lengths_per_supertile);
    const int tj = tile_list[i].y - (stj * tile_lengths_per_supertile);
    const int tij_map_index = ser.tile_map_idx[stij_map_index +
                                               (tj * tile_lengths_per_supertile) + ti];
    tile_instructions[i].y = tij_map_index;
    tmp_imports[(2 * i)    ] = (tile_list[i].x |
                                (static_cast<llint>(tile_list[i].z) << int_bit_count_int));
    tmp_imports[(2 * i) + 1] = (tile_list[i].y |
                                (static_cast<llint>(tile_list[i].z) << int_bit_count_int));
  }

  // As in the simple constructor involving just one system, we reduce the import list to its
  // unique values, then inflate the imports to read as atom (rather than tile) indices.  This
  // time, however, we add the system lower bound as well.  The imports array will then have
  // absolute atom indices to obtain from the coordinate arrays pertaining to the compilation of
  // systems, but it is still of interest what system they come from, as this will reveal the
  // maximum numbers of atoms in each system and thus the extent of the import (whether it goes
  // all the way to tile_length, or stops at the upper limit of system atoms).
  reduceUniqueValues(&tmp_imports);
  import_count = tmp_imports.size();
  for (int i = 0; i < import_count; i++) {
    import_system_indices[i] = static_cast<int>((tmp_imports[i] >> int_bit_count_int) &
                                                0x00000000ffffffffLL);
    imports[i] = (static_cast<int>(tmp_imports[i] & 0x00000000ffffffffLL) * tile_length) +
                 ser.atom_offsets[import_system_indices[i]];
  }
  for (int i = import_count; i < small_block_max_imports; i++) {
    imports[i] = -1;
  }
  score = import_count + (tile_count * 8);
  
  // Assemble work unit instructions
  for (int i = 0; i < import_count; i++) {
    const int system_idx = import_system_indices[i];
    const int tile_idx = (imports[i] - ser.atom_offsets[system_idx]) / tile_length;
    const uint absc_mask = 16 * i;
    const uint ordi_mask = (absc_mask << 16);
    for (int j = 0; j < tile_count; j++) {
      if (tile_list[j].x == tile_idx && tile_list[j].z == system_idx) {
        tile_instructions[j].x |= absc_mask;
      }
      if (tile_list[j].y == tile_idx && tile_list[j].z == system_idx) {
        tile_instructions[j].x |= ordi_mask;
      }
    }
  }

  // Fill out the size keys for each import, a series of bit-packed integers with one import batch
  // size (up to tile_length atoms, i.e. 16) per eight bits.
  if (constants::int_bit_count_int < 32) {
    rtErr("Descriptors for non-bonded work units require that signed integers be of size >= 32 "
          "bits.", "NonbondedWorkUnit");
  }
  for (int i = 0; i < import_count; i++) {
    const int key_idx = i / 4;
    const int key_pos = i - (key_idx * 4);
    const int system_idx = import_system_indices[i];
    const int atom_limit = ser.atom_offsets[system_idx] + ser.atom_counts[system_idx];
    int import_batch_size;
    if (imports[i] >= 0) {
      import_batch_size = std::min(atom_limit - imports[i], tile_length);
    }
    else {
      import_batch_size = 0;
    }
    import_size_keys[key_idx] |= (import_batch_size << (8 * key_pos));
  }
}

//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMaskSynthesis &se,
                                     const int abscissa_start, const int ordinate_start,
                                     const int system_index) :
    tile_count{0},
    import_count{0},
    kind{NbwuKind::SUPERTILES},
    score{0}, init_accumulation{0}, refresh_atom_start{0}, refresh_code{0},
    imports{},
    import_system_indices{std::vector<int>(1)},
    import_size_keys{},
    tile_instructions{}
{
  imports.resize(3);
  imports[0] = abscissa_start;
  imports[1] = ordinate_start;

  /// Obtain the system index
  const SeMaskSynthesisReader ser = se.data();
  import_system_indices[0] = getImportSystemIndex(ser, abscissa_start);
  if (import_system_indices[0] != getImportSystemIndex(ser, ordinate_start)) {
    rtErr("The abscissa and ordinate starting indices do not correspond to the same atoms.",
          "NonbondedWorkUnit");
  }
  
  // Compute the tile count and the work unit score
  const int atom_limit = se.getAtomCount(system_index);
  imports[2] = (abscissa_start + supertile_length > atom_limit) ? atom_limit - abscissa_start :
                                                                  supertile_length;
  imports[3] = (ordinate_start + supertile_length > atom_limit) ? atom_limit - ordinate_start :
                                                                  supertile_length;
  if (abscissa_start == ordinate_start) {
    if (abscissa_start < atom_limit - supertile_length &&
        ordinate_start < atom_limit - supertile_length) {
      tile_count = tile_lengths_per_supertile * (tile_lengths_per_supertile + 1) / 2;
      import_count = tile_lengths_per_supertile;
    }
    else {
      const int ntile_rem = (atom_limit - abscissa_start + tile_length - 1) / tile_length;  
      tile_count = ntile_rem * (ntile_rem + 1) / 2;
      import_count = ntile_rem;
    }
  }
  else {
    if (abscissa_start < atom_limit - supertile_length &&
        ordinate_start < atom_limit - supertile_length) {
      tile_count = tile_lengths_per_supertile * tile_lengths_per_supertile;
      import_count = 2 * tile_lengths_per_supertile;
    }
    else if (abscissa_start < atom_limit - supertile_length) {
      const int ntile_rem = (atom_limit - ordinate_start + tile_length - 1) / tile_length;
      tile_count = tile_lengths_per_supertile * ntile_rem;
      import_count = tile_lengths_per_supertile + ntile_rem;
    }
    else if (ordinate_start < atom_limit - supertile_length) {
      const int ntile_rem = (atom_limit - abscissa_start + tile_length - 1) / tile_length;
      tile_count = tile_lengths_per_supertile * ntile_rem;
      import_count = tile_lengths_per_supertile + ntile_rem;
    }
    else {
      rtErr("A trapezoidal tile should not exist in a supertile-based work unit.",
            "NonbondedWorkUnit");
    }
  }
  score = import_count + (tile_count * 8);
}
  
//-------------------------------------------------------------------------------------------------
int NonbondedWorkUnit::getTileCount() const {
  return tile_count;
}

//-------------------------------------------------------------------------------------------------
int NonbondedWorkUnit::getImportCount() const {
  return import_count;
}

//-------------------------------------------------------------------------------------------------
int NonbondedWorkUnit::getInitializationMask() const {
  return init_accumulation;
}

//-------------------------------------------------------------------------------------------------
int4 NonbondedWorkUnit::getTileLimits(const int index) const {
  const int absc_idx = (tile_instructions[index].x & 0xffff) / tile_length;
  const int ordi_idx = ((tile_instructions[index].x >> 16) & 0xffff) / tile_length;
  const int absc_key = absc_idx / 4;
  const int absc_pos = absc_idx - (absc_key * 4);
  const int ordi_key = ordi_idx / 4;
  const int ordi_pos = ordi_idx - (ordi_key * 4);
  const int absc_length = ((import_size_keys[absc_key] >> (8 * absc_pos)) & 0xff);
  const int ordi_length = ((import_size_keys[ordi_key] >> (8 * ordi_pos)) & 0xff);
  return { imports[absc_idx], imports[absc_idx] + absc_length,
           imports[ordi_idx], imports[ordi_idx] + ordi_length };
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& NonbondedWorkUnit::getTileInstructions() const {
  return tile_instructions;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> NonbondedWorkUnit::getAbstract(const int instruction_start) const {
  std::vector<int> result;
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    {
      result.resize(tile_groups_wu_abstract_length, -1);
      result[0] = import_count;
      for (int i = 0; i < import_count; i++) {
        result[i + 1] = imports[i];
      }
      for (int i = 0; i < 5; i++) {
        result[1 + small_block_max_imports + i] = import_size_keys[i];
      }
      result[small_block_max_imports + 6] = instruction_start;
      result[small_block_max_imports + 7] = instruction_start + tile_count;
      for (int i = 0; i < import_count; i++) {
        result[small_block_max_imports + 8 + i] = import_system_indices[i];
      }
      result[(2 * small_block_max_imports) +  8] = init_accumulation;
      result[(2 * small_block_max_imports) +  9] = refresh_atom_start;
      result[(2 * small_block_max_imports) + 10] = refresh_code;
    }
    return result;
  case NbwuKind::SUPERTILES:
    result.resize(supertile_wu_abstract_length);
    for (int i = 0; i < 4; i++) {
      result[i] = imports[i];
    }
    result[4] = (import_system_indices[0]);
    result[5] = init_accumulation;
    result[6] = refresh_atom_start;
    result[7] = refresh_code;
    return result;
  case NbwuKind::HONEYCOMB:
  case NbwuKind::UNKNOWN:
    return result;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int NonbondedWorkUnit::getImportSystemIndex(const SeMaskSynthesisReader &ser,
                                            const int atom_index) {
  if (atom_index >= ser.atom_offsets[ser.nsys - 1]) {
    if (atom_index >= ser.atom_offsets[ser.nsys - 1] + ser.atom_counts[ser.nsys - 1]) {
      rtErr("Atom index " + std::to_string(atom_index) + " is out of range in a synthesis of " +
            std::to_string(ser.nsys) + " systems with a maximum limit of " +
            std::to_string(ser.atom_offsets[ser.nsys - 1]) + " + " +
            std::to_string(ser.atom_counts[ser.nsys - 1]) + " atoms.", "NonbondedWorkUnit",
            "getImportSystemIndex");
    }
    return ser.nsys - 1;
  }
  else {
    return findBin(ser.atom_counts, atom_index, ser.nsys - 1);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void NonbondedWorkUnit::setInitializationMask(const int mask_in) {
  init_accumulation = mask_in;
}

//-------------------------------------------------------------------------------------------------
void NonbondedWorkUnit::setRefreshAtomIndex(int index_in) {
  refresh_atom_start = index_in;
}

//-------------------------------------------------------------------------------------------------
void NonbondedWorkUnit::setRefreshWorkCode(int code_in) {
  refresh_code = code_in;
}

//-------------------------------------------------------------------------------------------------
size_t estimateNonbondedWorkUnitCount(const std::vector<int> &atom_counts, const int nbwu_tile) {
  const size_t nsys = atom_counts.size();
  size_t result = 0LLU;
  for (size_t i = 0; i < nsys; i++) {
    const size_t ntile_side = (nbwu_tile == huge_nbwu_tiles) ?
                              (atom_counts[i] + supertile_length - 1) / supertile_length :
                              (atom_counts[i] + tile_length - 1) / tile_length;
    result += ntile_side * (ntile_side + 1LLU) / 2LLU;
  }
  if (nbwu_tile == tiny_nbwu_tiles || nbwu_tile == small_nbwu_tiles ||
      nbwu_tile == medium_nbwu_tiles || nbwu_tile == large_nbwu_tiles) {
    const size_t nbwu_tile_zu = static_cast<size_t>(nbwu_tile);
    return (result + nbwu_tile_zu - 1LLU) / nbwu_tile_zu;
  }
  else if (nbwu_tile == huge_nbwu_tiles) {
    return result;
  }
  else {
    rtErr("Only 8, 16, 32, 64, or " + std::to_string(huge_nbwu_tiles) + " tiles may be the target "
          "size of any one NonbondedWorkUnit.", "estimateNonbondedWorkUnitCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool addTileToWorkUnitList(int3* tile_list, int* import_coverage,
                           const std::vector<int> &system_tile_starts, int *import_count,
                           int *current_tile_count, const int ti, const int tj, const int sysid) {
  const int c_ic  = *import_count;
  if (ti == tj) {
    if (c_ic + 1 - import_coverage[ti + system_tile_starts[sysid]] > small_block_max_imports) {
      return false;
    }
    else {
      const size_t c_tile = *current_tile_count;
      tile_list[c_tile].x = ti;
      tile_list[c_tile].y = tj;
      tile_list[c_tile].z = sysid;
      *import_count = c_ic + 1 - import_coverage[ti + system_tile_starts[sysid]];
      import_coverage[ti + system_tile_starts[sysid]] = 1;
      *current_tile_count += 1;
      return true;
    }
  }
  else {
    if (c_ic + 2 - import_coverage[ti + system_tile_starts[sysid]] -
                   import_coverage[tj + system_tile_starts[sysid]] > small_block_max_imports) {
      return false;
    }
    else {
      const size_t c_tile = *current_tile_count;
      tile_list[c_tile].x = ti;
      tile_list[c_tile].y = tj;
      tile_list[c_tile].z = sysid;
      *import_count = c_ic + 2 - import_coverage[ti + system_tile_starts[sysid]] -
                                 import_coverage[tj + system_tile_starts[sysid]];
      import_coverage[ti + system_tile_starts[sysid]] = 1;
      import_coverage[tj + system_tile_starts[sysid]] = 1;
      *current_tile_count += 1;
      return true;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void setPropertyInitializationMasks(std::vector<NonbondedWorkUnit> *all_wu,
                                    const NbwuKind kind) {
  const size_t small_block_max_imports_zu = static_cast<size_t>(small_block_max_imports);
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    {
      const size_t nwu = all_wu->size();
      NonbondedWorkUnit* all_wu_data = all_wu->data();
      std::vector<int> import_counts(nwu);
      std::vector<int> import_list(nwu * small_block_max_imports_zu);
      int max_atom_idx = 0;
      for (size_t i = 0; i < nwu; i++) {
        import_counts[i] = all_wu_data[i].getImportCount();
        const std::vector<int> tile_abstract = all_wu_data[i].getAbstract();
        for (int j = 0; j < tile_abstract[0]; j++) {
          max_atom_idx = std::max(max_atom_idx, tile_abstract[j + 1]);
        }
        for (size_t j = 0LLU; j < small_block_max_imports_zu; j++) {
          import_list[j + (i * small_block_max_imports_zu)] = tile_abstract[j + 1LLU];
        }
      }
      std::vector<bool> coverage(max_atom_idx + 1, false);
      for (size_t i = 0; i < nwu; i++) {
        uint init_mask = 0;
        for (int j = 0; j < import_counts[i]; j++) {
          const size_t j_zu = j;
          const int jimp_start = import_list[j_zu + (i * small_block_max_imports_zu)];

          // Check that there are no duplicate imports in any work unit.
          for (int k = j + 1; k < import_counts[i]; k++) {
            const size_t k_zu = k;
            if (jimp_start == import_list[k_zu + (i * small_block_max_imports_zu)]) {
              rtErr("Duplicate imports were detected in non-bonded work unit " +
                    std::to_string(i) + ", slots " + std::to_string(j) + " and " +
                    std::to_string(k) + ".", "setPropertyInitializationMasks");
            }
          }
          if (coverage[jimp_start] == false) { 
            coverage[jimp_start] = true;
            init_mask |= (0x1 << j);
          }
        }
        all_wu_data[i].setInitializationMask(init_mask);
      }
    }
    break;
  case NbwuKind::SUPERTILES:
  case NbwuKind::HONEYCOMB:
  case NbwuKind::UNKNOWN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
int nonbondedWorkUnitInitCode(const InitializationTask init_request, const int cache_depth) {

  // The following codes use the format X / Y / Z force accumulators (bits 0, 1, or 2,
  // respectively), Generalized Born radii and derivative accumulators (bits 3 and 4), and the
  // depth of the random numbers cache (bits 8 through 11).  This code will be shifted 16 bits to
  // the left to accommodate the number of atoms that each work unit's initialization applies to.
  switch (init_request) {
  case InitializationTask::GENERAL_DYNAMICS:
    return 7;
  case InitializationTask::GB_DYNAMICS:
    return 31;
  case InitializationTask::GB_MINIMIZATION:
    return 24;
  case InitializationTask::LANGEVIN_DYNAMICS:
    return (7 | (cache_depth << 8));
  case InitializationTask::GB_LANGEVIN_DYNAMICS:
    return (31 | (cache_depth << 8));
  case InitializationTask::GENERAL_MINIMIZATION:
  case InitializationTask::NONE:
    return 0;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void distributeInitializationRanges(std::vector<NonbondedWorkUnit> *result, const int natom,
                                    const int init_code, const GpuDetails &gpu) {

  // Handle blank refresh instructions
  NonbondedWorkUnit* res_ptr = result->data();
  int nnbwu = result->size();
  if (init_code == 0) {
    for (int i = 0; i < nnbwu; i++) {
      res_ptr[i].setRefreshAtomIndex(0);
      res_ptr[i].setRefreshWorkCode(0);
    }
    return;
  }

  // Distribute the work of refreshing accumulators over the non-bonded work units.
  int n_task = 0;
  for (int i = 0; i < 8; i++) {
    if ((init_code >> i) & 0x1) {
      n_task++;
    }
  }
  const int n_atom_blocks = (natom + small_block_size - 1) / small_block_size;
  const bool compute_rng = ((init_code >> 8) > 0);
  if (compute_rng && gpu.getGpuSupported()) {
    
    // A non-null GPU was detected and will affect the profile of the workload.  If there are
    // enough blocks in the launch grid relative to the number of work units, expand the list of
    // non-bonded work units with a complement devoted to computing random numbers.
    if (nnbwu + n_atom_blocks < 2 * gpu.getSMPCount()) {
      const StaticExclusionMaskSynthesis blank_se;
      NonbondedWorkUnit blank_nbwu(blank_se, {});
      result->resize(nnbwu + n_atom_blocks, blank_nbwu);
      res_ptr = result->data();
      for (int i = nnbwu + n_atom_blocks - 1; i >= nnbwu; i--) {
        res_ptr[i] = res_ptr[i - nnbwu];
      }
      for (int i = 0; i < n_atom_blocks; i++) {
        res_ptr[i] = blank_nbwu;
      }
      nnbwu += n_atom_blocks;
    }
  }
  n_task += compute_rng;
  if (n_task * n_atom_blocks < nnbwu) {
    
    // If there are enough non-bonded work units to distribute just one set of initializations to
    // each, do each task individually, starting with the random number computations (placing these
    // at the front of the kernel gives more time for whatever queues the write instructions go
    // into to drain by the end of the kernel).
    int atom_idx = 0;
    int nbwu_idx = 0;
    int work_code = (init_code & 0xff00);
    if (work_code > 0) {
      while (atom_idx < natom) {
        res_ptr[nbwu_idx].setRefreshAtomIndex(atom_idx);
        const int nbatch = std::min(small_block_size, natom - atom_idx);
        res_ptr[nbwu_idx].setRefreshWorkCode(work_code | (nbatch << 16));
        atom_idx += small_block_size;
        nbwu_idx++;
      }
    }
    
    // Now distribute the initialization / refresh operations.
    work_code = 1;
    for (int i = 0; i < 8; i++) {
      if ((init_code >> i) & 0x1) {
        atom_idx = 0;
        while (atom_idx < natom) {
          res_ptr[nbwu_idx].setRefreshAtomIndex(atom_idx);
          const int nbatch = std::min(small_block_size, natom - atom_idx);
          res_ptr[nbwu_idx].setRefreshWorkCode(work_code | (nbatch << 16));
          atom_idx += small_block_size;
          nbwu_idx++;
        }
      }
      work_code *= 2;
    }
  }
  else {
    
    // Otherwise, do not consider the middle-ground cases where the work might be more evenly
    // distributed by assigning some combinations of tasks to various work units.  Every work
    // unit gets an assignment to initialize a stretch of atoms according to the same work code.
    // The overwhelming odds are that there are more than enough work units to spread the
    // initialization work thin, not taking this path.
    int work_code = (init_code & 0xffff);
    if (work_code > 0) {
      const int atom_blocks_per_wu = std::max((((natom + small_block_size - 1) /
                                                small_block_size) + nnbwu - 1) / nnbwu, 1);
      const int atoms_per_wu = atom_blocks_per_wu * small_block_size;
      int atom_idx = 0;
      int nbwu_idx = 0;
      while (atom_idx < natom) {
        res_ptr[nbwu_idx].setRefreshAtomIndex(atom_idx);
        const int nbatch = std::min(atoms_per_wu, natom - atom_idx);
        res_ptr[nbwu_idx].setRefreshWorkCode(work_code | (nbatch << 16));
        atom_idx += atoms_per_wu;
        nbwu_idx++;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const StaticExclusionMaskSynthesis &poly_se,
                        const InitializationTask init_request, const int random_cache_depth,
                        const GpuDetails &gpu) {
  
  // Determine the optimal overall size for work units.  Given that this process is guided by
  // static (as opposed to forward) exclusion masks, this is a matter of how many atoms are
  // present in all topologies and the number of tiles it would take to cover them all.
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  std::vector<int> atom_counts(poly_ser.nsys);
  std::vector<int> atom_offsets(poly_ser.nsys);
  for (int i = 0; i < poly_ser.nsys; i++) {
    atom_counts[i]  = poly_ser.atom_counts[i];
    atom_offsets[i] = poly_ser.atom_offsets[i];
  }
  const int padded_natom = atom_offsets[poly_ser.nsys - 1] +
                           roundUp(atom_counts[poly_ser.nsys - 1], warp_size_int);

  // Try work units of eight, sixteen, 32, and 64 tiles before going to the gargantuan,
  // 256-tile "huge" work units.
  const size_t tiny_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, tiny_nbwu_tiles);
  const size_t small_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, small_nbwu_tiles);
  const size_t medium_wu_count = estimateNonbondedWorkUnitCount(atom_counts, medium_nbwu_tiles);
  const size_t large_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, large_nbwu_tiles);
  const size_t huge_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, huge_nbwu_tiles);
  std::vector<NonbondedWorkUnit> result;
  if (tiny_wu_count < 2 * kilo) {
    result = enumerateNonbondedWorkUnits(poly_se, tiny_nbwu_tiles, tiny_wu_count, atom_counts,
                                         atom_offsets);
  }
  else if (small_wu_count < 4 * kilo) {
    result = enumerateNonbondedWorkUnits(poly_se, small_nbwu_tiles, small_wu_count, atom_counts,
                                         atom_offsets);
  }
  else if (medium_wu_count < 8 * kilo) {
    result = enumerateNonbondedWorkUnits(poly_se, medium_nbwu_tiles, medium_wu_count, atom_counts,
                                         atom_offsets);
  }
  else if (large_wu_count < 128 * kilo) {
    result = enumerateNonbondedWorkUnits(poly_se, large_nbwu_tiles, large_wu_count, atom_counts,
                                         atom_offsets);
  }
  else {    
    result = enumerateNonbondedWorkUnits(poly_se, huge_nbwu_tiles, huge_wu_count, atom_counts,
                                         atom_offsets);
  }
  
  // Add initialization instructions and return the result
  const int init_code = nonbondedWorkUnitInitCode(init_request, random_cache_depth);
  distributeInitializationRanges(&result, padded_natom, init_code, gpu);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const StaticExclusionMask &se, const InitializationTask init_request,
                        const int random_cache_depth) {
  const std::vector<int> atom_counts = { se.getAtomCount() };
  const std::vector<int> atom_offsets(1, 0);
  const size_t tiny_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, tiny_nbwu_tiles);
  const size_t small_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, small_nbwu_tiles);
  const size_t medium_wu_count = estimateNonbondedWorkUnitCount(atom_counts, medium_nbwu_tiles);
  const size_t large_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, large_nbwu_tiles);
  const size_t huge_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, huge_nbwu_tiles);
  std::vector<NonbondedWorkUnit> result;
  if (tiny_wu_count < 2 * kilo) {
    result = enumerateNonbondedWorkUnits(se, tiny_nbwu_tiles, tiny_wu_count, atom_counts,
                                         atom_offsets);
  }
  else if (small_wu_count < 4 * kilo) {
    result = enumerateNonbondedWorkUnits(se, small_nbwu_tiles, small_wu_count, atom_counts,
                                         atom_offsets);
  }
  else if (medium_wu_count < 16 * kilo) {
    result = enumerateNonbondedWorkUnits(se, medium_nbwu_tiles, medium_wu_count, atom_counts,
                                         atom_offsets);
  }
  else if (large_wu_count < 128 * kilo) {
    result = enumerateNonbondedWorkUnits(se, large_nbwu_tiles, large_wu_count, atom_counts,
                                         atom_offsets);
  }
  else {
    result = enumerateNonbondedWorkUnits(se, huge_nbwu_tiles, huge_wu_count, atom_counts,
                                         atom_offsets);
  }

  // Add initialization instructions and return the result
  const int init_code = nonbondedWorkUnitInitCode(init_request, random_cache_depth);
  distributeInitializationRanges(&result, se.getAtomCount(), init_code);
  return result;
}

} // namespace synthesis
} // namespace stormm
