#include <algorithm>
#include "copyright.h"
#include "Parsing/parse.h"
#include "forward_exclusionmask.h"

namespace stormm {
namespace energy {

using parse::char4ToString;
using topology::NonbondedKit;

//-------------------------------------------------------------------------------------------------
ForwardExclusionMaskReader::ForwardExclusionMaskReader(const int natom_in, const int mask2_size_in,
                                                       const int* primary_idx_in,
                                                       const uint2* mask1_in,
                                                       const uint2* mask2_in) :
  natom{natom_in},
  mask2_size{mask2_size_in},
  primary_idx{primary_idx_in},
  mask1{mask1_in},
  mask2{mask2_in}
{}

//-------------------------------------------------------------------------------------------------
ForwardExclusionMask::ForwardExclusionMask(const AtomGraph *ag_in) :
  atom_count{(ag_in == nullptr) ? 0 : ag_in->getAtomCount()},
  primary_mask_indices{static_cast<size_t>(atom_count), "fwdmask_idx"},
  primary_masks{static_cast<size_t>(atom_count), "fwdmask_primary"},
  secondary_masks{0, "fwdmask_secondary"},
  ag_pointer{const_cast<AtomGraph*>(ag_in)}
{
  // Return immediately if the pointer was null.
  if (ag_in == nullptr) {
    return;
  }
  
  // Step through each atom, assembling a list of exclusions whos atomic indices are ahead of it
  // in the topology list.  Determine whether the list will fit within a single 32-bit unsigned int
  // mask (the primary mask), and if not make additional masks (secondary masks).
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  std::vector<int> excl_atoms;
  std::vector<uint2> tmp_secondary_masks;
  int n_secondary_masks = 0;
  std::vector<uint2> tmp_primary_masks(atom_count);
  for (int i = 0; i < atom_count; i++) {

    // Determine an upper bound on the number of atoms this may involve (the actual number
    // should be exactly half of this).  Allocate enough memory to hold everything.
    const int dbl_excl_est = (nbk.nb11_bounds[i + 1] - nbk.nb11_bounds[i] +
                              nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i] +
                              nbk.nb13_bounds[i + 1] - nbk.nb13_bounds[i] +
                              nbk.nb14_bounds[i + 1] - nbk.nb14_bounds[i]);
    excl_atoms.resize(dbl_excl_est);

    // Step through each exclusion list and accumulate atoms.
    int nx = 0;
    for (int j = nbk.nb11_bounds[i]; j < nbk.nb11_bounds[i + 1]; j++) {
      if (nbk.nb11x[j] > i) {
        excl_atoms[nx] = nbk.nb11x[j];
        nx++;
      }
    }
    for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
      if (nbk.nb12x[j] > i) {
        excl_atoms[nx] = nbk.nb12x[j];
        nx++;
      }
    }
    for (int j = nbk.nb13_bounds[i]; j < nbk.nb13_bounds[i + 1]; j++) {
      if (nbk.nb13x[j] > i) {
        excl_atoms[nx] = nbk.nb13x[j];
        nx++;
      }
    }
    for (int j = nbk.nb14_bounds[i]; j < nbk.nb14_bounds[i + 1]; j++) {
      if (nbk.nb14x[j] > i) {
        excl_atoms[nx] = nbk.nb14x[j];
        nx++;
      }
    }
    excl_atoms.resize(nx);

    // Sort the exclusions in ascending order.
    std::sort(excl_atoms.begin(), excl_atoms.end(), [](int a, int b) { return a < b; });

    // Fill out the primary mask and commit it to the object's array
    uint pmask = 0U;
    int k = 0;
    while (k < nx && excl_atoms[k] <= i + 32) {
      pmask |= (0x1 << (excl_atoms[k] - i - 1));
      k++;
    }
    tmp_primary_masks[i].x = pmask;

    // Determine whether the primary mask covered all of the exclusions
    uint pextra = 0U;
    if (k < nx) {

      // Attempt to express the remaining atoms as an extension of the primary mask: the reference
      // atom index must be < i + 32768 and any additional exclusions must lie within 16 atoms of
      // the reference index.
      if (excl_atoms[k] - i < 32768 && nx - k <= 17 && excl_atoms[nx - 1] - excl_atoms[k] <= 17) {
        const int ref_idx = excl_atoms[k];
        pextra = ref_idx - i;
        for (int m = k + 1; m < nx; m++) {

          // The action is to take the index of the mth exclusion relative to the reference atom
          // (the kth exclusion), again subtracting 1 as the mask elements refere to the first,
          // second, third, ..., up to the 16th atom after the reference.  However, the mask is
          // moved forward by 15 bits to make room for the number of the reference atom, so
          // -1 + 15 = 14.  Add 14 to the bit shift.
          pextra |= (0x1 << (excl_atoms[m] - ref_idx + 14));
        }
      }
      else {

        // If there are already too many secondary masks, trap the case
        if (n_secondary_masks >= 4194304) {
          rtErr("The topology described by " + ag_pointer->getFileName() + " already has " +
                std::to_string(n_secondary_masks) + " secondary masks.  This is a very high "
                "number and the format cannot handle any more.", "ForwardExclusionMask");
        }

        // Accumulate a list of secondary masks for this atom.
        pextra = n_secondary_masks;
        bool cluster_engaged = false;
        uint2 contrib;
        while (k < nx) {
          if (cluster_engaged == false) {
            contrib.x = excl_atoms[k];
            contrib.y = 0U;
            k++;
            cluster_engaged = true;
          }
          else {
            if (excl_atoms[k] > contrib.x + 32) {
              tmp_secondary_masks.push_back(contrib);
              n_secondary_masks++;
              cluster_engaged = false;
            }
            else {
              contrib.y |= (0x1 << (excl_atoms[k] - contrib.x - 1));
              k++;
            }
          }
        }
        if (cluster_engaged) {
          tmp_secondary_masks.push_back(contrib);
          n_secondary_masks++;
        }

        // Trap cases where some atom requires too many secondary masks (this should be impossible
        // with any conventional force field describing a real chemical system, and extraordinarily
        // difficult even in a theoretical framework designed to break the format).
        if (n_secondary_masks - pextra >= 512) {
          rtErr("A total of " + std::to_string(n_secondary_masks) + " secondary exclusion masks "
                "were required to cover the exclusions of atom index " + std::to_string(i + 1) + 
                "(" + char4ToString(ag_pointer->getAtomName(i)) + " of residue " +
                char4ToString(ag_pointer->getResidueName(ag_pointer->getResidueIndex(i))) +
                ") in topology " + ag_pointer->getFileName() + ".", "ForwardExclusionMask");
        }

        // Record the number of secondary masks needed by this atom
        pextra |= ((n_secondary_masks - pextra) << 22);

        // Note that this atom has secondary masks, not just a short mask stuffed into the
        // primary mask array.
        pextra |= (0x1 << 31);
      }
    }

    // Commit the result
    tmp_primary_masks[i].y = pextra;
  }

  // The primary masks array now holds two unsigned integers arranged in a tuple for each atom.
  // However, there may be a great deal of repetition in this array, and the goal should be to
  // compress this to the smallest amount of information possible, as it is a reasonable goal to
  // cache all of the unique forward exclusion masks.  The secondary masks array might be
  // compressible as well, but it is so unlikely to be populated that there is likely no need.
  std::vector<bool> mask_catalogged(atom_count, false);
  std::vector<uint2> compacted_primary_masks;
  std::vector<int> tmp_primary_mask_indices(atom_count);
  int n_unique_masks = 0;
  for (int i = 0; i < atom_count; i++) {
    if (mask_catalogged[i]) {
      continue;
    }
    compacted_primary_masks.push_back(tmp_primary_masks[i]);
    for (int j = i; j < atom_count; j++) {
      if (tmp_primary_masks[j].x == tmp_primary_masks[i].x &&
          tmp_primary_masks[j].y == tmp_primary_masks[i].y) {
        mask_catalogged[j] = true;
        tmp_primary_mask_indices[j] = n_unique_masks;
      }
    }
    n_unique_masks++;
  }

  // Load the mask arrays in the object
  primary_mask_indices.putHost(tmp_primary_mask_indices);
  primary_masks.resize(compacted_primary_masks.size());
  primary_masks.putHost(compacted_primary_masks);
  secondary_masks.resize(tmp_secondary_masks.size());
  secondary_masks.putHost(tmp_secondary_masks);
}

//-------------------------------------------------------------------------------------------------
ForwardExclusionMask::ForwardExclusionMask(const AtomGraph &ag_in) :
    ForwardExclusionMask(ag_in.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
int ForwardExclusionMask::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int ForwardExclusionMask::getPrimaryMaskCount() const {
  return static_cast<int>(primary_masks.size());
}

//-------------------------------------------------------------------------------------------------
int ForwardExclusionMask::getExtendedMaskCount() const {
  int result = 0;
  const uint2* primary_ptr = primary_masks.data();
  const int nprim = primary_masks.size();
  for (int i = 0; i < nprim; i++) {
    const uint xmask = primary_ptr[i].y;
    result += (xmask > 0 && (xmask & 0x80000000) == 0x0);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int ForwardExclusionMask::getSecondaryMaskCount() const {
  return static_cast<int>(secondary_masks.size());
}

//-------------------------------------------------------------------------------------------------
int ForwardExclusionMask::getTotalExclusionsCount() const {
  const uint2* primary_ptr = primary_masks.data();
  const uint2* secondary_ptr = secondary_masks.data();
  const int nprim = primary_masks.size();
  const int nscnd = secondary_masks.size();
  std::vector<int> prim_excl(nprim);
  std::vector<int> scnd_excl(nscnd);

  // Compute secondary mask bit counts first.
  for (int i = 0; i < nscnd; i++) {
    const uint pmask = primary_ptr[i].y;
    int nexcl = 0;
    for (int j = 0; j < 32; j++) {
      nexcl += ((pmask >> j) & 0x1);
    }
    scnd_excl[i] = nexcl;
  }

  // Compute the number of exclusions in each primary mask.  Any secondary mask, including an
  // extended "half" mask, appended to the primary mask implies one additional exclusion plus
  // however many bits are checked in the actual mask component.
  for (int i = 0; i < nprim; i++) {
    const uint pmask = primary_ptr[i].x;
    int nexcl = 0;
    for (int j = 0; j < 32; j++) {
      nexcl += ((pmask >> j) & 0x1);
    }
    const uint xmask = primary_ptr[i].y;
    if (xmask & 0x80000000) {
      const int mask2_start = (xmask & 0x003fffff);
      const int mask2_length = ((xmask & 0x7fc00000) >> 22);
      nexcl++;
      for (int j = mask2_start; j < mask2_start + mask2_length; j++) {
        nexcl += scnd_excl[j];
      }
    }
    else if (xmask > 0) {
      nexcl++;
      for (int j = 15; j < 31; j++) {
        nexcl += ((xmask >> j) & 0x1);
      }
    }
    prim_excl[i] = nexcl;
  }

  // Sum the results for all primary masks
  const int* index_ptr = primary_mask_indices.data();
  int result = 0;
  for (int i = 0; i < atom_count; i++) {
    result += prim_excl[index_ptr[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
uint2 ForwardExclusionMask::getPrimaryMask(const int index) const {
  return primary_masks.readHost(primary_mask_indices.readHost(index));
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2>
ForwardExclusionMask::getSecondaryMasks(const int start, const int length) const {
  return secondary_masks.readHost(start, length);
}

//-------------------------------------------------------------------------------------------------
const ForwardExclusionMaskReader ForwardExclusionMask::data(HybridTargetLevel tier) const {
  return ForwardExclusionMaskReader(atom_count, static_cast<int>(secondary_masks.size()),
                                    primary_mask_indices.data(tier), primary_masks.data(tier),
                                    secondary_masks.data(tier));
}

//-------------------------------------------------------------------------------------------------
bool ForwardExclusionMask::testExclusion(const int atom_i, const int atom_j) const {

  // Ensure that atom_j > atom_i
  if (atom_j == atom_i) {
    return false;
  }
  if (atom_j < atom_i) {
    return testExclusion(atom_j, atom_i);
  }

  // Get the primary mask
  uint2 mask1 = primary_masks.readHost(primary_mask_indices.readHost(atom_i));

  // If atom J is within 32 indices of atom I, the primary mask will identify an exclusion
  const int dji_idx = atom_j - atom_i;
  if (dji_idx <= 32) {
    return ((mask1.x >> (dji_idx - 1)) & 0x1);
  }

  // If atom J is greater than 32 indices from atom I, there are more things to check.  First, are
  // there any additional exclusions listed in the second part of the primary mask?  If so, are
  // they a compact extra mask stuffed in the second part of the primary mask, or is there a whole
  // list of additional masks?
  if (mask1.y == 0) {
    return false;
  }
  else if (mask1.y & 0x80000000) {

    // Take in a list of additional secondary masks
    const int mask2_start = (mask1.y & 0x003fffff);

    // Use an integer divide to take some pressure off the bitwise operation pipelines
    const int mask2_length = ((mask1.y & 0x7fc00000) >> 22);
    const std::vector<uint2> mask2 = secondary_masks.readHost(mask2_start, mask2_length);

    // For each secondary mask, determine whether atom J falls within range of the reference
    // atom.  Check for an exclusion if it does.
    for (int i = 0; i < mask2_length; i++) {
      const int reference_atom = mask2[i].x;
      const int drj_idx = atom_j - reference_atom;
      if (drj_idx == 0) {

        // The reference atom itself must be an exclusion
        return true;
      }
      else if (drj_idx > 0 && drj_idx <= 32) {

        // Use this secondary mask to identify an exclusion
        return ((mask2[i].y >> (drj_idx - 1)) & 0x1);
      }
    }

    // No secondary mask identified an exclusion for atom J
    return false;
  }
  else {

    // Interpret the second part of the primary mask as half exclusion mask, half relative index
    // for atom I
    const int reference_atom = atom_i + (mask1.y & 0x00007fff);
    const int drj_idx = atom_j - reference_atom;
    if (drj_idx == 0) {
      return true;
    }
    else if (drj_idx > 0 && drj_idx <= 16) {

      // Use the half mask to identify an exclusion
      const uint mini_mask = ((mask1.y >> 15) & 0x0000ffff);
      return ((mini_mask >> (drj_idx - 1)) & 0x1);
    }
    else {

      // Atom J did not fall within the range of the primary mask or the secondary half mask
      return false;
    }
  }
  __builtin_unreachable();
}

} // namespace energy
} // namespace stormm
