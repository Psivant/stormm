#include <algorithm>
#include "copyright.h"
#include "Math/formulas.h"
#include "Topology/atomgraph_analysis.h"
#include "Topology/topology_util.h"
#include "local_exclusionmask.h"

namespace stormm {
namespace energy {

using card::HybridKind;
using card::setPointer;
using stmath::ipow;
using topology::hasVdwProperties;
using topology::inferCombiningRule;

//-------------------------------------------------------------------------------------------------
LocalExclusionMaskReader::LocalExclusionMaskReader(const int* prof_idx_in,
                                                   const ullint* profiles_in,
                                                   const uint2* aux_masks_in) :
    prof_idx{prof_idx_in}, profiles{profiles_in}, aux_masks{aux_masks_in}
{}
  
//-------------------------------------------------------------------------------------------------
LocalExclusionMask::LocalExclusionMask(const AtomGraph *ag_in, const NonbondedTheme theme_in) :
    total_atom_count{(ag_in == nullptr) ? 0 : ag_in->getAtomCount()},
    theme{theme_in},
    mode_a_profile_count{0}, mode_b_profile_count{0}, mode_c_profile_count{0},
    mode_d_profile_count{0}, mode_e_profile_count{0}, mode_f_profile_count{0},
    atom_profile_codes{static_cast<size_t>(total_atom_count), "locmask_idx"},
    atom_profiles{HybridKind::ARRAY, "locmask_atom_prof"},
    secondary_masks{HybridKind::ARRAY, "locmask_atom_bkp"},
    ag_pointer{const_cast<AtomGraph*>(ag_in)},
    poly_ag_pointer{nullptr}
{
  setTopology(ag_in);
  countProfileKinds();
}

//-------------------------------------------------------------------------------------------------
LocalExclusionMask::LocalExclusionMask(const AtomGraph &ag_in, const NonbondedTheme theme_in) :
  LocalExclusionMask(ag_in.getSelfPointer(), theme_in)
{}

//-------------------------------------------------------------------------------------------------
LocalExclusionMask::LocalExclusionMask(const AtomGraphSynthesis *poly_ag_in,
                                       const NonbondedTheme theme_in) :
    total_atom_count{(poly_ag_in == nullptr) ? 0 : poly_ag_in->getPaddedAtomCount()},
    theme{theme_in},
    mode_a_profile_count{0}, mode_b_profile_count{0}, mode_c_profile_count{0},
    mode_d_profile_count{0}, mode_e_profile_count{0}, mode_f_profile_count{0},
    atom_profile_codes{static_cast<size_t>(total_atom_count), "locmask_idx"},
    atom_profiles{HybridKind::ARRAY, "locmask_atom_prof"},
    secondary_masks{HybridKind::ARRAY, "locmask_atom_bkp"},
    ag_pointer{nullptr},
    poly_ag_pointer{const_cast<AtomGraphSynthesis*>(poly_ag_in)}
{
  setTopology(poly_ag_in);
  countProfileKinds();
}

//-------------------------------------------------------------------------------------------------
LocalExclusionMask::LocalExclusionMask(const AtomGraphSynthesis &poly_ag_in,
                                       const NonbondedTheme theme_in) :
    LocalExclusionMask(poly_ag_in.getSelfPointer(), theme_in)
{}

//-------------------------------------------------------------------------------------------------
int LocalExclusionMask::getAtomCount() const {
  return total_atom_count;
}

//-------------------------------------------------------------------------------------------------
int LocalExclusionMask::getProfileCount() const {
  return atom_profiles.size();
}

//-------------------------------------------------------------------------------------------------
int LocalExclusionMask::getProfileCount(const int mode_index) const {
  switch (mode_index) {
  case 0:
    return mode_a_profile_count;
  case 1:
    return mode_b_profile_count;
  case 2:
    return mode_c_profile_count;
  case 3:
    return mode_d_profile_count;
  case 4:
    return mode_e_profile_count;
  case 5:
    return mode_f_profile_count;
  default:
    rtErr("Invalid profile mode index " + std::to_string(mode_index) + ".  Indices 0-6 are "
          "valid.", "LocalExclusionMask", "getProfileCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int LocalExclusionMask::getMode(const int atom_index, const int system_index) const {
  validateAtomIndex(atom_index, system_index, "getMode");
  const int actual_index = (system_index > 0 && poly_ag_pointer != nullptr) ?
                           atom_index + poly_ag_pointer->getAtomOffset(system_index) : atom_index;
  const ullint prof = (atom_profiles.readHost(atom_profile_codes.readHost(actual_index)) &
                       lmask_mode_bitmask);
  return static_cast<int>(prof >> (64 - lmask_mode_bits));
}

//-------------------------------------------------------------------------------------------------
size_t LocalExclusionMask::getSecondaryMaskCount() const {
  return secondary_masks.size();
}

//-------------------------------------------------------------------------------------------------
ullint LocalExclusionMask::getProfile(const int atom_index, const int system_index) const {
  validateAtomIndex(atom_index, system_index, "getProfile");
  const int atom_offset = poly_ag_pointer->getAtomOffset(system_index);
  return atom_profiles.readHost(atom_profile_codes.readHost(atom_offset));
}

//-------------------------------------------------------------------------------------------------
const Hybrid<uint2> LocalExclusionMask::getSecondaryMaskView(const int atom_index,
                                                             const int system_index) {
  validateAtomIndex(atom_index, system_index, "getSecondaryMaskView");
  ullint prof;
  if (poly_ag_pointer != nullptr) {
    prof = getProfile(atom_index, system_index);
  }
  else if (ag_pointer != nullptr) {
    prof = getProfile(atom_index);
  }
  switch (prof & lmask_mode_bitmask) {
  case lmask_mode_a:
  case lmask_mode_b:
  case lmask_mode_c:
  case lmask_mode_f:
    rtErr("Atom " + std::to_string(atom_index) + " of system " + std::to_string(system_index) +
          " does not use a profile mode that entails secondary masks.");
  case lmask_mode_d:
    {
      // The resulting Hybrid object will point to a place in the secondary masks array given by
      // some middle bits the profile with a length determined by yet higher bits.  The exclusion
      // mask includes some direct, near indices in the lower  bits of the profile.
      const size_t offset = ((prof & lmask_d_array_idx) >> lmask_d_array_idx_pos);
      const size_t length = ((prof & lmask_d_array_cnt) >> lmask_d_array_cnt_pos);
      return setPointer(&secondary_masks, offset, length, "lmask_viewer");
    }
    break;
  case lmask_mode_e:
    {
      // The resulting Hybrid object will point to a place in the secondary masks array given by
      // the low bits of the profile with a length determined by higher bits.  There is no
      // exclusion mask incorporated into a profile of this type.
      const size_t offset = (prof & lmask_e_array_idx);
      const size_t length = ((prof & lmask_e_array_cnt) >> lmask_e_array_cnt_pos);
      return setPointer<uint2>(&secondary_masks, offset, length, "lmask_viewer");
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool LocalExclusionMask::testExclusion(const int atom_i, const int atom_j) const {
  const ullint prof = atom_profiles.readHost(atom_profile_codes.readHost(atom_i));
  const int del_ij = atom_j - atom_i;
  switch (prof & lmask_mode_bitmask) {
  case lmask_mode_a:
    {
      return (abs(del_ij) <= lmask_long_local_span &&
              ((prof >> (lmask_long_local_span + del_ij)) & 0x1));
    }
    break;
  case lmask_mode_b:
    {
      const int abs_dij = abs(del_ij);
      if (abs_dij > lmask_b_max_reach) {
        return false;
      }
      else if (abs_dij <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }
      else if (del_ij > 0) {
        const int upper_shft = ((prof & lmask_b_upper_shft) >> lmask_b_upper_shft_pos);
        const int upper_mask_start = lmask_short_local_span + upper_shft;
        const int rel_ij = del_ij - upper_mask_start;
        return (rel_ij >= 0 && rel_ij < lmask_short_extra_span &&
                ((prof >> (rel_ij + lmask_b_upper_mask_pos)) & 0x1));
      }
      else {

        // The only remaining case is that del_ij < 0
        const int lower_shft = ((prof & lmask_b_lower_shft) >> lmask_b_lower_shft_pos);
        const int lower_mask_start = -lmask_short_local_span - lower_shft - lmask_short_extra_span;
        const int rel_ij = del_ij - lower_mask_start;
        return (rel_ij >= 0 && rel_ij < lmask_short_extra_span &&
                ((prof >> (rel_ij + lmask_b_lower_mask_pos)) & 0x1));
      }
    }
    break;
  case lmask_mode_c:
    {
      if (abs(del_ij) <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }

      // Forming the unsigned long long int on the r.h.s. and then converting it to a signed
      // short int will translate the bit string appropriately.
      const int alt_mask_shft = static_cast<short int>((prof & lmask_c_shft) >> lmask_c_shft_pos);
      
      // Run the shift in terms of the index atom
      const int rel_ij = del_ij - alt_mask_shft;
      return (rel_ij >= 0 && rel_ij < lmask_long_extra_span &&
              ((prof >> (rel_ij + lmask_c_alt_mask_pos)) & 0x1));
    }
    break;
  case lmask_mode_d:
    {
      if (abs(del_ij) <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }

      // This is the best possible path.  Obtain the number of masks and loop over all of them.
      const size_t nmasks = ((prof & lmask_d_array_cnt) >> lmask_d_array_cnt_pos);
      const size_t start_idx = ((prof & lmask_d_array_idx) >> lmask_d_array_idx_pos);
      const uint2* secondary_ptr = secondary_masks.data();
      for (size_t i = 0; i < nmasks; i++) {
        const uint2 tmask = secondary_ptr[start_idx + i];
        const int tmask_x = tmask.x;
        if (del_ij >= tmask_x && del_ij < tmask_x + 32 &&
            ((tmask.y >> (del_ij - tmask_x)) & 0x1)) {
          return true;
        }
      }
      return false;
    }
    break;
  case lmask_mode_e:
    {
      // This is the bet possible path and there is no local exclusion arrangement to test.  Loop
      // over all the masks.
      const size_t nmasks = ((prof & lmask_e_array_cnt) >> lmask_e_array_cnt_pos);
      const size_t start_idx = (prof & lmask_e_array_idx);
      const uint2* secondary_ptr = secondary_masks.data();
      for (size_t i = 0; i < nmasks; i++) {
        const uint2 tmask = secondary_ptr[start_idx + i];
        const int tmask_x = tmask.x;
        if (del_ij >= tmask_x && del_ij < tmask_x + 32 &&
            ((tmask.y >> (del_ij - tmask_x)) & 0x1)) {
          return true;
        }
      }
      return false;
    }
    break;
  case lmask_mode_f:
    break;
  default:
    rtErr("No known profile code was matched to " + lMaskModeToString(prof) + ".",
          "LocalExclusionMask", "testExclusion");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool LocalExclusionMask::testExclusion(const int atom_i, const int atom_j,
                                       const int system_index) const {
  if (poly_ag_pointer == nullptr) {
    rtErr("A synthesis of topologies must be associated in order to compute exclusions within a "
          "specific system.", "LocalExclusionMask", "testExclusion");
  }
  const int sys_offset = poly_ag_pointer->getAtomOffset(system_index);
  const int atom_count = poly_ag_pointer->getAtomCount(system_index);
  if (atom_i < 0 || atom_i >= atom_count) {
    rtErr("Atom index " + std::to_string(atom_i) + " is invalid for a system of " +
          std::to_string(atom_count) + " atoms.", "LocalExclusionMask", "testExclusion");
  }
  if (atom_j < 0 || atom_j >= atom_count) {
    rtErr("Atom index " + std::to_string(atom_j) + " is invalid for a system of " +
          std::to_string(atom_count) + " atoms.", "LocalExclusionMask", "testExclusion");
  }
  return testExclusion(atom_i + sys_offset, atom_j + sys_offset);
}

//-------------------------------------------------------------------------------------------------
LocalExclusionMaskReader LocalExclusionMask::data(const HybridTargetLevel tier) const {
  return LocalExclusionMaskReader(atom_profile_codes.data(tier), atom_profiles.data(tier),
                                  secondary_masks.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::upload() {
  atom_profile_codes.upload();
  atom_profiles.upload();
  secondary_masks.upload();
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::download() {
  atom_profile_codes.download();
  atom_profiles.download();
  secondary_masks.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::setTopology(const AtomGraph* ag_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_in);
  if (ag_pointer != nullptr) {
    total_atom_count = ag_pointer->getAtomCount();
    atom_profile_codes.resize(total_atom_count);
    const std::vector<int> tmp_prof_codes = extendMasks(ag_pointer);
    atom_profile_codes.putHost(tmp_prof_codes, 0, total_atom_count);
  }
  else {

    // If the object is not already set to serve a topology synthesis, associating it with a null
    // topology will wipe it clean.
    if (poly_ag_pointer == nullptr) {
      nullify();
    }
  }
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::setTopology(const AtomGraph& ag_in) {
  setTopology(ag_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::setTopology(const AtomGraphSynthesis* poly_ag_in) {
  poly_ag_pointer = const_cast<AtomGraphSynthesis*>(poly_ag_in);
  if (poly_ag_pointer != nullptr) {

    // Loop over all topologies, taking the offset and getting exclusions from the system-specific
    // non-bonded kit.
    const int nsys = poly_ag_pointer->getSystemCount();
    std::vector<bool> coverage(nsys, false);
    const std::vector<int> pag_indices = poly_ag_pointer->getTopologyIndices();
    total_atom_count = poly_ag_pointer->getPaddedAtomCount();
    atom_profile_codes.resize(total_atom_count);
    for (int pos = 0; pos < nsys; pos++) {
      if (coverage[pos]) {
        continue;
      }
      const int pos_topology_idx = pag_indices[pos];
      const AtomGraph *ag_ptr = poly_ag_pointer->getSystemTopologyPointer(pos);
      std::vector<int> tmp_prof_codes = extendMasks(ag_ptr);
      for (int i = pos; i < nsys; i++) {
        if (pag_indices[i] == pos_topology_idx) {
          const int system_offset = poly_ag_pointer->getAtomOffset(i);
          atom_profile_codes.putHost(tmp_prof_codes, system_offset, ag_ptr->getAtomCount());
        }
      }
    }
  }
  else {

    // If the object is not already associated with a valid topology, setting the synthesis pointer
    // to null will wipe the object clean.
    if (ag_pointer == nullptr) {
      nullify();
    }
  }
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::setTopology(const AtomGraphSynthesis& poly_ag_in) {
  setTopology(poly_ag_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
std::vector<int> LocalExclusionMask::extendMasks(const NonbondedKit<double> &nbk,
                                                 const ChemicalDetailsKit &cdk) {
  VdwCombiningRule lj_rule;
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
  case NonbondedTheme::ALL:
    lj_rule = VdwCombiningRule::GEOMETRIC;
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    lj_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types);
    break;
  }
  std::vector<int> result(nbk.natom);

  // Create an array of residue names corresponding to each atom, to expedite finding new examples
  // of any given atom and residue name once a profile is set.
  std::vector<char4> atomized_residue_names(cdk.natom);
  for (int i = 0; i < cdk.nres; i++) {
    for (int j = cdk.res_limits[i]; j < cdk.res_limits[i + 1]; j++) {
      atomized_residue_names[j] = cdk.res_names[i];
    }
  }
  
  // Seed an A-mode profile (no exclusions) if none yet exists.
  std::vector<ullint> tmp_profiles = atom_profiles.readHost();
  std::vector<uint2> tmp_secondary_masks = secondary_masks.readHost();
  std::vector<bool> atoms_mapped(nbk.natom, false);
  if (tmp_profiles.size() == 0) {
    tmp_profiles.push_back(lmask_mode_a);
  }

  // Find profiles for other atoms
  for (int pos = 0; pos < nbk.natom; pos++) {

    // Skip atoms that have already been mapped.
    if (atoms_mapped[pos]) {
      continue;
    }

    // The first A mode profile applies if the atom itself is not relevant to the non-bonded list.
    if (hasRelevantProperties(nbk, pos, theme, lj_rule) == false) {
      result[pos] = 0;
      atoms_mapped[pos] = true;
      continue;
    }
    std::vector<int> excluded_atoms = compileLocalExclusionList(nbk, pos, theme, lj_rule);

    // The first A mode profile applies if there are no exclusions.
    if (excluded_atoms.size() == 0) {
      result[pos] = 0;
      atoms_mapped[pos] = true;
      continue;
    }
    
    // Determine the type of profile that this will require.  Begin by checking the range of
    // relevant exclusions.  The switch below will add at most one entry to the the temporary
    // atom profiles array, through one and only one branch, so it is safe to record its present
    // size here.
    const size_t nprof = tmp_profiles.size();
    ullint tprof_mode;
    bool mode_found = false;
    if (lMaskFitsModeA(pos, excluded_atoms)) {
      mode_found = true;
      tprof_mode = lmask_mode_a;
    }
    else if (lMaskFitsModeB(pos, excluded_atoms)) {
      mode_found = true;
      tprof_mode = lmask_mode_b;
    }
    else if (lMaskFitsModeC(pos, excluded_atoms)) {
      mode_found = true;
      tprof_mode = lmask_mode_c;
    }
    if (mode_found) {
      
      // One of the mode A, B, or C profiles will work and is preferred.  Create the profile.
      const ullint tprof = lMaskCreateProfile(pos, tprof_mode, excluded_atoms);
      setProfile(pos, tprof, &result, &tmp_profiles);
      atoms_mapped[pos] = true;
      continue;
    }

    // The D profile can mop up other cases but involves references to the array of secondary
    // masks.  It is very unlikely that any simulation will require even one of these profiles.
    // However, because the profile indexes into the array of secondary masks with a lower bound
    // than the F profile, all possible Ds must be evaluated, until either no more atoms fit the
    // description or there is no more room for D profiles in the accessible space of secondary
    // masks.
    if (lMaskFitsModeD(pos, excluded_atoms, tmp_secondary_masks)) {
      ullint tprof;
      const std::vector<uint2> tprof_secondary_masks = lMaskCreateProfile(pos, lmask_mode_d,
                                                                          excluded_atoms, &tprof);
      
      // Check whether the profile has already been included, along with its complement of
      // secondary masks.
      setProfile(pos, tprof, lmask_mode_d, tprof_secondary_masks, &result, &tmp_profiles,
                 &tmp_secondary_masks);
      atoms_mapped[pos] = true;
    }
  }
  
  // If nothing else, any set of exclusions will fit mode E.  There are too many possible
  // enumerations for it not to work.  Implement this solution for any atoms that do not yet have
  // profiles.
  for (int pos = 0; pos < nbk.natom; pos++) {

    // An evaluation of relevant properties was done in the previous loop.  Atoms without relevant
    // properties were marked as "mapped" and will be skipped in this loop as well.
    if (atoms_mapped[pos]) {
      continue;
    }
    std::vector<int> excluded_atoms = compileLocalExclusionList(nbk, pos, theme, lj_rule);
    ullint tprof;
    const std::vector<uint2> tprof_secondary_masks = lMaskCreateProfile(pos, lmask_mode_e,
                                                                        excluded_atoms, &tprof);
    setProfile(pos, tprof, lmask_mode_e, tprof_secondary_masks, &result, &tmp_profiles,
               &tmp_secondary_masks);
    atoms_mapped[pos] = true;
  }
  
  // Contribute all profiles and secondary masks to the object's main arrays.
  atom_profiles.resize(tmp_profiles.size());
  atom_profiles.putHost(tmp_profiles, 0, tmp_profiles.size());
  secondary_masks.resize(tmp_secondary_masks.size());
  secondary_masks.putHost(tmp_secondary_masks, 0, tmp_secondary_masks.size());
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> LocalExclusionMask::extendMasks(const AtomGraph *ag) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  return extendMasks(nbk, cdk);
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::nullify() {
  total_atom_count = 0;
  mode_a_profile_count = 0;
  mode_b_profile_count = 0;
  mode_c_profile_count = 0;
  mode_d_profile_count = 0;
  mode_e_profile_count = 0;
  mode_f_profile_count = 0;
  atom_profile_codes.resize(0);
  atom_profiles.resize(0);
  secondary_masks.resize(0);
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::validateAtomIndex(const int atom_index, const int system_index,
                                           const char* caller) const {
  if (atom_index < 0 || atom_index >= total_atom_count) {
    rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a collection of " +
          std::to_string(total_atom_count) + " atoms.", "LocalExclusionMask", caller);
  }
  else if (ag_pointer != nullptr && system_index > 0) {
    rtErr("A nonzero system index (" + std::to_string(system_index) + ") is incompatible with a "
          "mask associated with only a single topology.", "LocalExclusionMask", caller);
  }
  if (poly_ag_pointer != nullptr) {
    if (system_index > 0) {
      const int sys_atoms = poly_ag_pointer->getAtomCount(system_index);
      if (atom_index >= sys_atoms) {
        rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a system (index " +
              std::to_string(system_index) + ") with " + std::to_string(sys_atoms) + " atoms.",
              "LocalExclusionMask", caller);
      }
    }
    else if (atom_index >= total_atom_count) {

      // Allow an atom index anywhere in the synthesis to be submitted with system zero, even if
      // it would violate the atom bounds of system zero.
      rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a collection of systems "
            "with " + std::to_string(total_atom_count) + " atoms in total.", "LocalExclusionMask",
            caller);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::validateAtomIndex(const int atom_index, const char* caller) const {
  validateAtomIndex(atom_index, 0, caller);
}

//-------------------------------------------------------------------------------------------------
void LocalExclusionMask::countProfileKinds() {

  // Accumulate the numbers of each profile type.
  const ullint* atom_profile_ptr = atom_profiles.data();
  const size_t n_profiles = atom_profiles.size();
  mode_a_profile_count = 0;
  mode_b_profile_count = 0;
  mode_c_profile_count = 0;
  mode_d_profile_count = 0;
  mode_e_profile_count = 0;
  mode_f_profile_count = 0;
  for (size_t i = 0; i < n_profiles; i++) {
    switch (atom_profile_ptr[i] & lmask_mode_bitmask) {
    case lmask_mode_a:
      mode_a_profile_count += 1;
      break;
    case lmask_mode_b:
      mode_b_profile_count += 1;
      break;
    case lmask_mode_c:
      mode_c_profile_count += 1;
      break;
    case lmask_mode_d:
      mode_d_profile_count += 1;
      break;
    case lmask_mode_e:
      mode_e_profile_count += 1;
      break;
    case lmask_mode_f:
      mode_f_profile_count += 1;
      break;
      break;
    default:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> compileLocalExclusionList(const NonbondedKit<double> &nbk, const int atom_index,
                                           const NonbondedTheme theme,
                                           const VdwCombiningRule lj_rule) {

  // Loop over all 1:1, 1:2, 1:3, and 1:4 non-bonded exclusions, accumulating a list that counts
  // all forward and backward interactions.  
  const int llim11 = nbk.nb11_bounds[atom_index];
  const int hlim11 = nbk.nb11_bounds[atom_index + 1];
  const int llim12 = nbk.nb12_bounds[atom_index];
  const int hlim12 = nbk.nb12_bounds[atom_index + 1];
  const int llim13 = nbk.nb13_bounds[atom_index];
  const int hlim13 = nbk.nb13_bounds[atom_index + 1];
  const int llim14 = nbk.nb14_bounds[atom_index];
  const int hlim14 = nbk.nb14_bounds[atom_index + 1];
  const int total_i_excl = (hlim11 - llim11) + (hlim12 - llim12) + (hlim13 - llim13) +
                           (hlim14 - llim14);
  std::vector<int> result(total_i_excl);
  int nx = 0;
  for (int i = llim11; i < hlim11; i++) {

    // Depending on the non-bonded theme, check the properties of each atom to ensure that it
    // should be counted as an exclusion.  The atom index will remain the same, and so exclusion
    // masks for a non-bonded neighbor list specific to van-der Waals interactions will remain
    // as far-flung as they would when serving a neighbor list of all atoms.  However, in the
    // overwhelming majority of cases the amount of data transmission will remain the same.
    if (hasRelevantProperties(nbk, nbk.nb11x[i], theme, lj_rule)) {
      result[nx] = nbk.nb11x[i];
      nx++;
    }
  }
  for (int i = llim12; i < hlim12; i++) {
    if (hasRelevantProperties(nbk, nbk.nb12x[i], theme, lj_rule)) {
      result[nx] = nbk.nb12x[i];
      nx++;
    }
  }
  for (int i = llim13; i < hlim13; i++) {
    if (hasRelevantProperties(nbk, nbk.nb13x[i], theme, lj_rule)) {
      result[nx] = nbk.nb13x[i];
      nx++;
    }
  }
  for (int i = llim14; i < hlim14; i++) {
    if (hasRelevantProperties(nbk, nbk.nb14x[i], theme, lj_rule)) {
      result[nx] = nbk.nb14x[i];
      nx++;
    }
  }
  result.resize(nx);

  // As with the forward mask, sort the atoms in ascending order.  Duplicates will have been pruned
  // at the stage when non-bonded 1:1, 1:2, 1:3, and 1:4 exclusions were computed for the topology.
  std::sort(result.begin(), result.end(), [](int a, int b) { return a < b; });
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string lMaskModeToString(const ullint mode) {
  const uint mode_rep = (mode >> lmask_mode_pos);
  return std::to_string((mode_rep & 0x4) >> 2) + std::to_string((mode_rep & 0x2) >> 1) +
         std::to_string(mode_rep & 0x1);
}

//-------------------------------------------------------------------------------------------------
bool lMaskFitsModeA(const int atom_index, const std::vector<int> &excluded_atoms) {
  const int nexcl = excluded_atoms.size();
  if (nexcl == 0) {
    return true;
  }
  else {
    return (excluded_atoms[0]         >= atom_index - lmask_long_local_span &&
            excluded_atoms[nexcl - 1] <= atom_index + lmask_long_local_span);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool lMaskFitsModeB(const int atom_index, const std::vector<int> &excluded_atoms) {
  const int nexcl = excluded_atoms.size();
  const int max_reach = lmask_short_local_span + lmask_short_extra_span +
                        ipow(2, lmask_b_shft_bits);
  if (nexcl == 0 || atom_index - excluded_atoms[0] > max_reach ||
      excluded_atoms[nexcl - 1] - atom_index > max_reach) {
    return false;
  }
  else {
    
    // Mark off all exclusions near enough to the atom as to fall under the proximate mask.
    std::vector<bool> coverage(nexcl);
    for (int i = 0; i < nexcl; i++) {
      coverage[i] = (abs(excluded_atoms[i] - atom_index) <= lmask_short_local_span);
    }

    // The two additional masks must fall below and above the base atom index, respectively--there
    // is no way for the profile configuration to resolve the indexing otherwise.
    bool lower_occ = false;
    if (coverage[0] == false && excluded_atoms[0] < atom_index) {
      const int lower_mask_anchor = excluded_atoms[0];
      int i = 0;
      while (i < nexcl && coverage[i] == false &&
             excluded_atoms[i] < lower_mask_anchor + lmask_short_extra_span) {
        coverage[i] = true;
        lower_occ = true;
        i++;
      }
    }
    bool upper_occ = false;
    if (coverage[nexcl - 1] == false && excluded_atoms[nexcl - 1] > atom_index) {
      const int upper_mask_ceiling = excluded_atoms[nexcl - 1];
      int i = nexcl - 1;
      while (i >= 0 && coverage[i] == false &&
             excluded_atoms[i] > upper_mask_ceiling - lmask_short_extra_span) {
        coverage[i] = true;
        upper_occ = true;
        i--;
      }
    }
    bool all_covered = true;
    for (int i = 0; i < nexcl; i++) {
      all_covered = (all_covered && coverage[i]);
    }
    return (all_covered && lower_occ && upper_occ);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool lMaskFitsModeC(const int atom_index, const std::vector<int> &excluded_atoms) {
  const int nexcl = excluded_atoms.size();
  if (nexcl == 0) {
    return false;
  }

  // Mark off all exclusions which can fall into the proximate mask
  std::vector<bool> coverage(nexcl);
  for (int i = 0; i < nexcl; i++) {
    coverage[i] = (abs(excluded_atoms[i] - atom_index) <= lmask_short_local_span);
  }

  // The additional mask must be able to cover all atoms within its allotted 14 bits.
  bool seek_bounds = true;
  int min_uncovered_loc, max_uncovered_loc;
  for (int i = 0; i < nexcl; i++) {
    if (coverage[i] == false) {
      if (seek_bounds) {
        min_uncovered_loc = i;
        seek_bounds = false;
      }
      max_uncovered_loc = i;
    }
  }

  // The mask of other atoms must begin with an atom that lies within 2^15 of the central atom
  // (atom_index) on the topological list.  Furthermore, if all atoms were covered, then the mask
  // would fit mode A, not C.
  const int max_shift = ipow(2, lmask_c_shft_bits - 1);
  if ((atom_index > excluded_atoms[min_uncovered_loc] &&
       atom_index - excluded_atoms[min_uncovered_loc] > max_shift) ||
      (atom_index < excluded_atoms[min_uncovered_loc] &&
       excluded_atoms[min_uncovered_loc] - atom_index > max_shift - 1) || seek_bounds) {
    return false;
  }
  return (excluded_atoms[max_uncovered_loc] - excluded_atoms[min_uncovered_loc] <
          lmask_long_extra_span);
}

//-------------------------------------------------------------------------------------------------
bool lMaskFitsModeD(const int atom_index, const std::vector<int> &excluded_atoms,
                    const std::vector<uint2> &tmp_secondary_masks) {
  const int nexcl = excluded_atoms.size();
  if (nexcl == 0) {
    return false;
  }

  // Mark off all exclusions which can fall into the proximate mask
  std::vector<bool> coverage(nexcl);
  for (int i = 0; i < nexcl; i++) {
    coverage[i] = (abs(excluded_atoms[i] - atom_index) <= lmask_short_local_span);
  }

  // Count the number of groups that would be needed to cover all remaining atoms.
  int n_extra_group = 0;
  bool on_group = false;
  int group_start;
  for (int i = 0; i < nexcl; i++) {
    if (coverage[i] == false) {

      // The atom is not covered.  If tracking a new group, try adding the atom to that group.
      // If the atom fits, continue, otherwise close the group.
      if (on_group) {
        if (excluded_atoms[i] >= group_start + 32) {
          n_extra_group++;
          group_start = excluded_atoms[i];
        }
      }
      else {
        on_group = true;
        group_start = excluded_atoms[i];
      }
    }
  }

  // Add one extra group for the final excluded atom(s).  If no groups were ever found, then all
  // exclusions in fact fit within the local mask and this would be a mode A profile.
  if (on_group) {
    n_extra_group++;
  }
  else {
    return false;
  }
  return (n_extra_group < ipow(2, lmask_d_array_cnt_bits) &&
          tmp_secondary_masks.size() + n_extra_group < ipow(2, lmask_d_array_idx_bits));
}

//-------------------------------------------------------------------------------------------------
ullint lMaskCreateProfile(const int atom_index, const ullint mode,
                          const std::vector<int> &excluded_atoms) {
  const int nx = excluded_atoms.size();
  ullint tprof = mode;
  switch (mode) {
  case lmask_mode_a:
    tprof |= (0x1LLU << lmask_long_local_span);
    for (int i = 0; i < nx; i++) {
      tprof |= (0x1LLU << (lmask_long_local_span + excluded_atoms[i] - atom_index));
    }
    break;
  case lmask_mode_b:
    {
      std::vector<bool> coverage(nx, false);
      tprof |= (0x1LLU << lmask_short_local_span);
      for (int i = 0; i < nx; i++) {
        if (abs(excluded_atoms[i] - atom_index) <= lmask_short_local_span) {
          tprof |= (0x1LLU << (lmask_short_local_span + excluded_atoms[i] - atom_index));
          coverage[i] = true;
        }
      }

      // If there is a lower mask to include, it will be evident in the first element of the
      // excluded atoms array.
      if (coverage[0] == false && excluded_atoms[0] < atom_index) {

        // Check whether the shift might be zero.  This can happen in cases where the extended
        // mask on one side of the central atom requires a nonzero shift but the other does not.
        const ullint lower_shift = std::max(atom_index - lmask_short_local_span -
                                            lmask_short_extra_span - excluded_atoms[0], 0);
        ullint lower_mask = 0LLU;
        int i = 0;
        while (coverage[i] == false &&
               excluded_atoms[i] < excluded_atoms[0] + lmask_short_extra_span) {
          lower_mask |= (0x1LLU << (excluded_atoms[i] - excluded_atoms[0]));
          i++;
        }
        tprof |= (lower_shift << ((2 * lmask_short_local_span) + 1 +
                                  (2 * lmask_short_extra_span)));
        tprof |= (lower_mask << ((2 * lmask_short_local_span) + 1));
      }

      // If there is an upper mask to include, it will be evident in the last element of the
      // excluded atoms array.
      if (coverage[nx - 1] == false && excluded_atoms[nx - 1] > atom_index) {
        const ullint upper_shift = std::max(excluded_atoms[nx - 1] - atom_index -
                                            lmask_short_extra_span - lmask_short_local_span, 0);
        tprof |= (upper_shift << ((2 * lmask_short_local_span) + 1 +
                                  (2 * lmask_short_extra_span) + lmask_b_shft_bits));

        // All of the remaining atoms will fit in the upper mask.
        ullint upper_mask = 0LLU;
        for (int i = 0; i < nx; i++) {
          upper_mask |= (0x1LLU << (excluded_atoms[i] - atom_index - lmask_short_local_span -
                                    upper_shift));
          i++;
        }
        tprof |= (upper_mask << ((2 * lmask_short_local_span) + 1 + lmask_short_extra_span));
      }
    }
    break;
  case lmask_mode_c:
    {
      // The C profile is also more complex that the A profile, not as complicated to evaluate but
      // in some ways more limited than the B profile.  The C mode is, overall, more common than
      // B.  Begin building the mask with local exclusions.
      std::vector<bool> coverage(nx, false);
      tprof |= (0x1LLU << lmask_short_local_span);
      for (int i = 0; i < nx; i++) {
        if (abs(excluded_atoms[i] - atom_index) <= lmask_short_local_span) {
          tprof |= (0x1LLU << (lmask_short_local_span + excluded_atoms[i] - atom_index));
          coverage[i] = true;
        }
      }

      // Loop over all atoms until encountering the unaccounted atom with the lowest topological
      // index.  Record the necessary offset and continue building the additional mask.
      int iseek = 0;
      while (iseek < nx && coverage[iseek]) {
        iseek++;
      }
      const ushort offset = excluded_atoms[iseek] - atom_index;
      const ullint u_ofs = offset;
      tprof |= (u_ofs << lmask_c_shft_pos);
      ullint tmask = 0LLU;
      for (int i = iseek; i < nx; i++) {
        if (coverage[i] == false) {
          tmask |= (0x1LLU << (excluded_atoms[i] - excluded_atoms[iseek]));
        }
      }
      tprof |= (tmask << ((2 * lmask_short_local_span) + 1));
    }
    break;
  case lmask_mode_d:
  case lmask_mode_e:
    rtErr("Use an alternative overload to return the list of secondary masks associated with "
          "mode " + lMaskModeToString(mode) + ".", "lMaskCreateProfile");
  case lmask_mode_f:
    rtErr("This profile cannot be created in the typical manner.  It is created as a "
          "memory-intensive but rapid resolution to the most difficult exclusion cases, if needed "
          "at all.", "lMaskCreateProfile");
  }
  return tprof;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2> lMaskCreateProfile(const int atom_index, const ullint mode,
                                      const std::vector<int> &excluded_atoms, ullint *t_profile) {
  std::vector<uint2> result;
  const int nx = excluded_atoms.size();
  switch (mode) {
  case lmask_mode_a:
  case lmask_mode_b:
  case lmask_mode_c:
    *t_profile = lMaskCreateProfile(atom_index, mode, excluded_atoms);
    break;
  case lmask_mode_d:
  case lmask_mode_e:
    {
      *t_profile = mode;

      // The D profile can mop up other cases but involves references to the array of secondary
      // masks.  It is very unlikely that any simulation will require even one of these profiles.
      // However, because the profile indexes into the array of secondary masks with a lower bound
      // than the E profile, all possible Ds must be evaluated, until either no more atoms fit the
      // description or there is no more room for E profiles in the accessible space of secondary
      // masks.
      std::vector<bool> coverage(nx, false);
      if (mode == lmask_mode_d) {
        ullint tmp_prof = 0LLU;
        tmp_prof |= (0x1LLU << lmask_short_local_span);
        for (int i = 0; i < nx; i++) {
          if (abs(excluded_atoms[i] - atom_index) <= lmask_short_local_span) {
            tmp_prof |= (0x1LLU << (lmask_short_local_span + excluded_atoms[i] - atom_index));
            coverage[i] = true;
          }
        }
        *t_profile |= tmp_prof;
      }
      
      // Loop over all atoms.  As different groups of atoms are encountered, they will fill the
      // array of secondary masks.  The check leading into this branch ensures that the secondary
      // masks array (the working copy being tmp_secondary_masks, a Standard Template Library
      // object rather than a Hybrid object to expedite reallocation and growth) does not exceed
      // the format limits.
      int iseek = 0;
      while (coverage[iseek]) {
        iseek++;
      }
      int group_start_topl = excluded_atoms[iseek];
      int group_start_list = iseek;
      for (int i = iseek; i < nx; i++) {
        if (coverage[i] == false) {
          if (excluded_atoms[i] >= group_start_topl + 32) {
            uint tmask = 0U;
            for (int j = group_start_list; j < i; j++) {
              if (coverage[i] == false) {
                tmask |= (0x1LLU << (excluded_atoms[j] - group_start_topl));
              }
            }
            result.push_back({ static_cast<uint>(group_start_topl - atom_index), tmask });
            group_start_list = i;
            group_start_topl = excluded_atoms[i];
          }
        }
      }
      uint final_tmask = 0U;
      for (int i = group_start_list; i < nx; i++) {
        if (coverage[i] == false) {
          final_tmask |= (0x1LLU << (excluded_atoms[i] - group_start_topl));
        }
      }
      result.push_back({ static_cast<uint>(group_start_topl - atom_index), final_tmask });
      const ullint n_masks = result.size();
      if (mode == lmask_mode_d) {
        *t_profile |= (n_masks << lmask_d_array_cnt_pos);
      }
      else {
        *t_profile |= (n_masks << lmask_e_array_cnt_pos);        
      }
    }
    break;
  case lmask_mode_f:
    rtErr("This profile cannot be created in the typical manner.  It is created as a "
          "memory-intensive but rapid resolution to the most difficult exclusion cases, if needed "
          "at all.", "lMaskCreateProfile");
  }
  return result;    
}

//-------------------------------------------------------------------------------------------------
void setProfile(const int pos, const ullint tprof, std::vector<int> *result,
                std::vector<ullint> *tmp_profiles) {

  // The zero index will always be the one (and only) A mode profile with no non-self exclusions.
  // Check whether the atom's profile is already known or add it to the repository if not, and
  // assign the proper index to the result.
  const ullint* tmp_profile_ptr = tmp_profiles->data();
  const size_t nprof = tmp_profiles->size();
  size_t i = 0;
  while (i < nprof && tmp_profile_ptr[i] != tprof) {
    i++;
  }

  // Set the result for the specific atom.
  if (i == nprof) {
    tmp_profiles->push_back(tprof);
    result->at(pos) = nprof;
  }
  else {
    result->at(pos) = i;
  }
}

//-------------------------------------------------------------------------------------------------
void setProfile(const int pos, const ullint tprof, const ullint mode,
                const std::vector<uint2> &tprof_secondary_masks, std::vector<int> *result,
                std::vector<ullint> *tmp_profiles, std::vector<uint2> *tmp_secondary_masks) {

  // Lay out parsing parameters depending on the mode.  The D mode involves a local exclusion list,
  // but the E mode does not.
  ullint count_mask, index_mask, count_pos, index_pos;
  if (mode == lmask_mode_d) {
    count_mask = lmask_d_array_cnt;
    index_mask = lmask_d_array_idx;
    count_pos  = lmask_d_array_cnt_pos;
    index_pos  = lmask_d_array_idx_pos;
  }
  else if (mode == lmask_mode_e) {
    count_mask = lmask_e_array_cnt;
    index_mask = lmask_e_array_idx;
    count_pos  = lmask_e_array_cnt_pos;
    index_pos  = 0;
  }
  else {
    rtErr("The only modes using this profile setting function are E and F.", "setProfile");
  }
  const ullint mode_plus_count_mask = (mode == lmask_mode_d) ?
                                      (count_mask | lmask_mode_bitmask | lmask_d_excl) :
                                      (count_mask | lmask_mode_bitmask);
  
  // Check to see whether there are any profiles that might match this one--by scanning for
  // matches of the profile itself first, the search can be greatly abbreviated.
  size_t i = 0;
  bool array_unique = true;
  const size_t nprof = tmp_profiles->size();
  const ullint* tmp_prof_ptr = tmp_profiles->data();
  const uint2* tmp_scnd_mask_ptr = tmp_secondary_masks->data();
  while (i < nprof && array_unique) {
    array_unique = ((tprof & mode_plus_count_mask) != (tmp_prof_ptr[i] & mode_plus_count_mask));
    if (array_unique == false) {
      
      // This profile has the same number of additional masks associated with it as another
      // profile which has already been recorded.  Attempt to recover the uniqueness of the
      // current mask array by comparing the individual masks to what has already been found.
      const size_t seq_length = ((tprof & count_mask) >> count_pos);
      const size_t cmp_start  = ((tmp_prof_ptr[i] & index_mask) >> index_pos);
      for (size_t i = 0; i < seq_length; i++) {
        const uint2 seq_mask = tprof_secondary_masks[i];
        const uint2 cmp_mask = tmp_scnd_mask_ptr[cmp_start + i];
        array_unique = (array_unique || seq_mask.x != cmp_mask.x || seq_mask.y != cmp_mask.y);
      }
    }

    // If the present mask array is indeed unique, continue the search.
    if (array_unique) {
      i++;
    }
  }
  if (i == nprof) {

    // Complete the new profile with the index of the first additional mask.  Add the profile and
    // all of its associated masks to the proper arrays.
    const ullint current_mask_count = tmp_secondary_masks->size();
    const int index_shift = (mode == lmask_mode_d) ? lmask_d_array_idx_pos : 0;
    tmp_profiles->push_back(tprof | (current_mask_count << index_shift));
    tmp_secondary_masks->insert(tmp_secondary_masks->end(), tprof_secondary_masks.begin(),
                                tprof_secondary_masks.end());
    result->at(pos) = nprof;
  }
  else {
    result->at(pos) = i;
  }
}

} // namespace energy
} // namespace stormm
