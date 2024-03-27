#include <cmath>
#include <cstdio>
#include <climits>
#include "copyright.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "atomgraph.h"
#include "topology_util.h"

namespace stormm {
namespace topology {

using card::HybridTargetLevel;
using stmath::findBin;
using stmath::locateValue;
using stmath::reduceUniqueValues;
using parse::char4ToString;
  
//-------------------------------------------------------------------------------------------------
std::string AtomGraph::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
std::string AtomGraph::getFileName() const {
  return source;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getResidueCount() const {
  return residue_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getMoleculeCount() const {
  return molecule_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLargestResidueSize() const {
  return largest_residue_size;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLastSoluteResidue() const {
  return last_solute_residue;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLastSoluteAtom() const {
  return last_solute_atom;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getFirstSolventMolecule() const {
  return first_solvent_molecule;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getFirstSolventAtom() const {

  // The first solvent molecule may be outside of the system bounds, indicating that there is, in
  // fact, no solvent.  Trap that case and return -1.
  if (first_solvent_molecule < 0 || first_solvent_molecule >= molecule_count) {
    return -1;
  }
  return molecule_contents.readHost(molecule_limits.readHost(first_solvent_molecule));
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLargestMoleculeSize() const {
  return largest_molecule_size;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getTotalMass() const {
  double tmass = 0.0;
  const double* mass_ptr = atomic_masses.data();
  for (int i = 0; i < atom_count; i++) {
    tmass += mass_ptr[i];
  }
  return tmass;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDegreesOfFreedom() const {
  return unconstrained_dof;
}
  
//-------------------------------------------------------------------------------------------------
int AtomGraph::getConstrainedDegreesOfFreedom() const {
  return constrained_dof;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDescriptor(const TopologyDescriptor choice) const {
  return descriptors.readHost(static_cast<ulint>(choice));
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDescriptor(const SanderDescriptor choice) const {
  return descriptors.readHost(static_cast<ulint>(choice));
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int>& AtomGraph::getResidueLimits() const {
  return residue_limits;
}

//-------------------------------------------------------------------------------------------------
int2 AtomGraph::getResidueLimits(const int index) const {
  int2 tmp = {residue_limits.readHost(index), residue_limits.readHost(index + 1)};
  return tmp;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueIndex() const {

  // Fill in an entire vector based on the residue limits array
  std::vector<int> result(atom_count);
  for (int i = 0; i < residue_count; i++) {
    const int llim = residue_limits.readHost(i);
    const int hlim = residue_limits.readHost(i + 1);
    for (int j = llim; j < hlim; j++) {
      result[j] = i;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getResidueIndex(const int atom_index) const {

  // This will still happen on the fly, rather than storing a long list of numbers.  It's simply
  // not common, or critical in performant code, to access the residue number.
  return findBin(residue_limits.data(), atom_index, residue_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getStructuralAtomNumber() const {
  return atom_struc_numbers.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getStructuralAtomNumber(const int low_index,
                                                    const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getStructuralAtomNumber");
  return atom_struc_numbers.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getStructuralAtomNumber(const int index) const {
  return atom_struc_numbers.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueNumber() const {
  return residue_numbers.readHost(0, residue_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueNumber(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getResidueNumber");
  return residue_numbers.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getResidueNumber(const int index) const {
  return residue_numbers.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeLimits() const {
  return molecule_limits.readHost(0, molecule_count + 1);
}

//-------------------------------------------------------------------------------------------------
int2 AtomGraph::getMoleculeLimits(const int index) const {
  int2 tmp = {molecule_limits.readHost(index), molecule_limits.readHost(index + 1)};
  return tmp;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getParticlesPerMolecule() const {
  std::vector<int> atoms_per_molecule = molecule_limits.readHost(0, molecule_count + 1);
  for (int i = 0; i < molecule_count; i++) {
    atoms_per_molecule[i] = atoms_per_molecule[i + 1] - atoms_per_molecule[i];
  }
  atoms_per_molecule.pop_back();
  return atoms_per_molecule;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getParticlesPerMolecule(const int index) const {
  return molecule_limits.readHost(index + 1) - molecule_limits.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getAtomicNumber() const {
  return atomic_numbers.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getAtomicNumber(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getAtomicNumber");
  return atomic_numbers.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAtomicNumber(const int index) const {
  return atomic_numbers.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> AtomGraph::getAtomMobility() const {
  return getAtomMobility(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> AtomGraph::getAtomMobility(const int low_index, const int high_index) const {

  // Range check as this will use the pointer
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getAtomMobility");
  std::vector<bool> mobiles(high_index - low_index, true);
  const int* m_ptr = mobile_atoms.data();
  const int int_bits = sizeof(int) * 8;
  for (int i = low_index; i < high_index; i++) {
    const int access_index = i / int_bits;
    mobiles[i - low_index] = ((static_cast<uint>(m_ptr[access_index]) >>
                               (i - (access_index * int_bits))) & 0x1);
  }
  return mobiles;
}

//-------------------------------------------------------------------------------------------------
bool AtomGraph::getAtomMobility(const int index) const {
  const int int_bits = sizeof(int) * 8;
  const int access_index = index / int_bits;
  const uint m_val = static_cast<uint>(mobile_atoms.readHost(access_index));
  return ((m_val >> (index - (access_index * int_bits))) & 0x1);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> AtomGraph::getAtomMobilityMask() const {
  return getAtomMobilityMask(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> AtomGraph::getAtomMobilityMask(const int low_index, const int high_index) const {
  const int uint_bits = sizeof(uint) * 8;
  std::vector<uint> result((high_index - low_index + uint_bits - 1) / uint_bits, 0U);
  int result_pos = 0;
  int result_bit = 0;
  const int* mobility_ptr = mobile_atoms.data();
  int mask_pos = low_index / uint_bits;
  int mask_bit = low_index - mask_pos * uint_bits;
  for (int i = low_index; i < high_index; i++) {
    result[result_pos] |= (((static_cast<uint>(mobility_ptr[mask_pos]) >> mask_bit) & 0x1) <<
                           result_bit);
    result_bit++;
    if (result_bit == uint_bits) {
      result_bit = 0;
      result_pos++;
    }
    mask_bit++;
    if (mask_bit == uint_bits) {
      mask_bit = 0;
      mask_pos++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeMembership() const {
  return molecule_membership.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeMembership(const int low_index,
                                                  const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getMoleculeMembership");
  return molecule_membership.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getMoleculeMembership(const int index) const {
  return molecule_membership.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeContents() const {
  return molecule_contents.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeContents(const int index) const {
  return molecule_contents.readHost(molecule_limits.readHost(index),
                                    molecule_limits.readHost(index + 1));
}

//-------------------------------------------------------------------------------------------------
const Hybrid<char4>& AtomGraph::getAtomName() const {
  return atom_names;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getAtomName(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getAtomName");
  return atom_names.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
char4 AtomGraph::getAtomName(const int index) const {
  return atom_names.readHost(index);
}

//-------------------------------------------------------------------------------------------------
const Hybrid<char4>& AtomGraph::getAtomType() const {
  return atom_types;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getAtomType(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getAtomType");
  return atom_types.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
char4 AtomGraph::getAtomType(const int index) const {
  return atom_types.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<char4>> AtomGraph::getAtomTypeNameTable() const {
  std::vector<std::vector<char4>> result(lj_type_count);
  std::vector<bool> coverage(lj_type_count, false);
  int types_located = 0;
  int i = 0;
  const int* ljidx_ptr = lennard_jones_indices.data();
  const char4* atyp_ptr = atom_types.data();
  while (i < atom_count && types_located < lj_type_count) {
    if (ljidx_ptr[i] < 0 || ljidx_ptr[i] >= lj_type_count) {
      rtErr("Lennard-Jones type index " + std::to_string(ljidx_ptr[i]) + " is invalid.",
            "AtomGraph", "getAtomTypeNameTable");
    }
    if (coverage[ljidx_ptr[i]] == false) {
      coverage[ljidx_ptr[i]] = true;
      std::vector<uint> tmp_tname_list;
      int n_sharing = 1;
      for (int j = i + 1; j < atom_count; j++) {
        n_sharing += (ljidx_ptr[j] == ljidx_ptr[i]);
      }
      tmp_tname_list.resize(n_sharing);
      tmp_tname_list[0] = char4ToUint(atyp_ptr[i]);
      n_sharing = 1;
      for (int j = i + 1; j < atom_count; j++) {
        if (ljidx_ptr[j] == ljidx_ptr[i]) {
          tmp_tname_list[n_sharing] = char4ToUint(atyp_ptr[j]);
          n_sharing++;
        }
      }
      reduceUniqueValues(&tmp_tname_list);
      result[ljidx_ptr[i]] = std::vector<char4>(tmp_tname_list.size());
      const size_t n_tname = tmp_tname_list.size();
      for (size_t j = 0; j < n_tname; j++) {
        result[ljidx_ptr[i]][j] = uintToChar4(tmp_tname_list[j]);
      }
      types_located++;
    }
    i++;
  }
  if (types_located < lj_type_count) {
    rtErr("The topology contains only " + std::to_string(types_located) + " unique Lennard-Jones "
          "atom types, whereas it is stated to contain " + std::to_string(lj_type_count) + ".",
          "AtomGraph", "getAtomTypeNameTable");
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
const Hybrid<char4>& AtomGraph::getResidueName() const {
  return residue_names;
}

//-------------------------------------------------------------------------------------------------
char4 AtomGraph::getResidueName(const int index) const {
  return residue_names.readHost(index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getUreyBradleyTermCount() const {
  return urey_bradley_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCharmmImprTermCount() const {
  return charmm_impr_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCmapTermCount() const {
  return cmap_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getUreyBradleyParameterCount() const {
  return urey_bradley_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCharmmImprParameterCount() const {
  return charmm_impr_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCmapSurfaceCount() const {
  return cmap_surface_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCmapDimension(const int index) const {
  return cmap_surface_dimensions.readHost(index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getBondTermCount() const {
  return bond_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAngleTermCount() const {
  return angl_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDihedralTermCount() const {
  return dihe_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getBondParameterCount() const {
  return bond_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAngleParameterCount() const {
  return angl_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDihedralParameterCount() const {
  return dihe_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteCount() const {
  return virtual_site_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteParameterSetCount() const {
  return virtual_site_parameter_set_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::findVirtualSites() const {
  return virtual_site_atoms.readHost(0, virtual_site_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> AtomGraph::findVirtualSites(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "findVirtualSites");
  std::vector<int2> result;
  const int* vstmp = virtual_site_atoms.data();
  const int nvs = virtual_site_count;
  for (int i = 0; i < nvs; i++) {
    if (vstmp[i] >= low_index && vstmp[i] < high_index) {
      result.push_back({vstmp[i], i});
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::findVirtualSites(const int index) const {
  return virtual_site_atoms.readHost(index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteIndex(const int atom_index, const ExceptionResponse policy) const {
  const int result = locateValue(virtual_site_atoms, atom_index);

  // Check that the atom found is, indeed, a virtual site
  if (virtual_site_atoms.readHost(result) != atom_index) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Search for the virtual site index of atom index " + std::to_string(atom_index) +
            "yielded a virtual site with atom index " +
            std::to_string(virtual_site_atoms.readHost(result)) +
            ".  The atomic number of atom index " + std::to_string(atom_index) + " is " +
            std::to_string(atomic_numbers.readHost(atom_index)) + ".", "AtomGraph",
            "getVirtualSiteIndex");
    case ExceptionResponse::WARN:
      return -1;
    case ExceptionResponse::SILENT:
      return -1;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKind AtomGraph::getVirtualSiteFrameType(const int index) const {
  const int vs_pidx = virtual_site_parameter_indices.readHost(index);
  return static_cast<VirtualSiteKind>(virtual_site_frame_types.readHost(vs_pidx));
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteFrameAtom(const int index, const int nfrm) const {
  switch (nfrm) {
  case 1:
    return virtual_site_frame1_atoms.readHost(index);
  case 2:
    return virtual_site_frame2_atoms.readHost(index);
  case 3:
    return virtual_site_frame3_atoms.readHost(index);
  case 4:
    return virtual_site_frame4_atoms.readHost(index);
  default:
    rtErr("Virtual sites cannot have a frame atom number " + std::to_string(nfrm) + ".",
          "AtomGraph", "getVirtualSiteFrameAtom");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getChargeTypeCount() const {
  return charge_type_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLJTypeCount() const {
  return lj_type_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getTotalExclusions() const {
  return total_exclusions;
}

//-------------------------------------------------------------------------------------------------
UnitCellType AtomGraph::getUnitCellType() const {
  return periodic_box_class;
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventModel AtomGraph::getImplicitSolventModel() const {
  return gb_style;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getDielectricConstant() const {
  return dielectric_constant;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getSaltConcentration() const {
  return salt_concentration;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getCoulombConstant() const {
  return coulomb_constant;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getElectrostatic14Screening() const {
  return elec14_screening_factor;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getVanDerWaals14Screening() const {
  return vdw14_screening_factor;
}
  
//-------------------------------------------------------------------------------------------------
std::string AtomGraph::getPBRadiiSet() const {
  return pb_radii_set;
}

//-------------------------------------------------------------------------------------------------
std::string AtomGraph::getWaterResidueName() const {
  return char4ToString(water_residue_name);
}
  
//-------------------------------------------------------------------------------------------------
int AtomGraph::getRigidWaterCount() const {
  int result = 0;
  const int* ox_ptr = settle_oxygen_atoms.data();
  const int* h1_ptr = settle_hydro1_atoms.data();
  const int* h2_ptr = settle_hydro2_atoms.data();
  const int* zn_ptr = atomic_numbers.data();
  for (int i = 0; i < settle_group_count; i++) {
    result += (zn_ptr[ox_ptr[i]] == 8 && zn_ptr[h1_ptr[i]] == 1 && zn_ptr[h2_ptr[i]] == 1);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getSettleGroupCount() const {
  return settle_group_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getBondConstraintCount() const {
  return bond_constraint_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getConstraintGroupCount() const {
  return constraint_group_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getConstraintGroupAtoms(const int index) const {
  if (index < 0 || index >= constraint_group_count) {
    rtErr("Constraint group index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(constraint_group_count) + " constraint groups.", "AtomGraph",
          "getConstraintGroupAtoms");
  }
  std::vector<int> result;
  const int hlim = constraint_group_bounds.readHost(index + 1);
  for (int i = constraint_group_bounds.readHost(index); i < hlim; i++) {
    result.push_back(constraint_group_atoms.readHost(i));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getConstraintGroupTotalSize() const {
  return constraint_group_bounds.readHost(constraint_group_count);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getNonrigidParticleCount() const {
  return nonrigid_particle_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getChargeIndex() const {
  return charge_indices.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getChargeIndex(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getChargeIndex");
  return charge_indices.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getChargeIndex(const int index) const {
  return charge_indices.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getLennardJonesIndex() const {
  return lennard_jones_indices.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getLennardJonesIndex(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getLennardJonesIndex");
  return lennard_jones_indices.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLennardJonesIndex(const int index) const {
  return lennard_jones_indices.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getAtomExclusions(const int index) const {
  
  // Assemble 1:1, then 1:2, 1:3, and finally 1:4 exclusions.  They will all get returned.
  std::vector<int> result = getNonbonded11Exclusions(index);
  std::vector<int> result2 = getNonbonded12Exclusions(index);
  std::vector<int> result3 = getNonbonded13Exclusions(index);
  std::vector<int> result4 = getNonbonded14Exclusions(index);
  result.insert(result.end(), result2.begin(), result2.end());
  result.insert(result.end(), result3.begin(), result3.end());
  result.insert(result.end(), result4.begin(), result4.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded11Exclusions(const int index) const {
  std::vector<int> result;
  const int nb11_low = nb11_exclusion_bounds.readHost(index);
  const int nb11_high = nb11_exclusion_bounds.readHost(index + 1);
  const int* nb11_list = nb11_exclusion_list.data();
  for (int i = nb11_low; i < nb11_high; i++) {
    result.push_back(nb11_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded12Exclusions(const int index) const {
  std::vector<int> result;
  const int nb12_low = nb12_exclusion_bounds.readHost(index);
  const int nb12_high = nb12_exclusion_bounds.readHost(index + 1);
  const int* nb12_list = nb12_exclusion_list.data();
  for (int i = nb12_low; i < nb12_high; i++) {
    result.push_back(nb12_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded13Exclusions(const int index) const {
  std::vector<int> result;
  const int nb13_low = nb13_exclusion_bounds.readHost(index);
  const int nb13_high = nb13_exclusion_bounds.readHost(index + 1);
  const int* nb13_list = nb13_exclusion_list.data();
  for (int i = nb13_low; i < nb13_high; i++) {
    result.push_back(nb13_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded14Exclusions(const int index) const {
  std::vector<int> result;
  const int nb14_low = nb14_exclusion_bounds.readHost(index);
  const int nb14_high = nb14_exclusion_bounds.readHost(index + 1);
  const int* nb14_list = nb14_exclusion_list.data();
  for (int i = nb14_low; i < nb14_high; i++) {
    result.push_back(nb14_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ValenceKit<double> AtomGraph::getDoublePrecisionValenceKit(const HybridTargetLevel tier) const {
  return ValenceKit<double>(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                            bond_parameter_count, angl_parameter_count, dihe_parameter_count,
                            inferred_14_attenuations, attenuated_14_type_count,
                            urey_bradley_term_count, charmm_impr_term_count, cmap_term_count,
                            urey_bradley_parameter_count, charmm_impr_parameter_count,
                            cmap_surface_count, bond_stiffnesses.data(tier),
                            bond_equilibria.data(tier), angl_stiffnesses.data(tier),
                            angl_equilibria.data(tier), dihe_amplitudes.data(tier),
                            dihe_periodicities.data(tier), dihe_phase_angles.data(tier),
                            attn14_elec_factors.data(tier), attn14_vdw_factors.data(tier),
                            bond_i_atoms.data(tier), bond_j_atoms.data(tier),
                            bond_parameter_indices.data(tier), bond_modifiers.data(tier),
                            angl_i_atoms.data(tier), angl_j_atoms.data(tier),
                            angl_k_atoms.data(tier), angl_parameter_indices.data(tier),
                            angl_modifiers.data(tier), dihe_i_atoms.data(tier),
                            dihe_j_atoms.data(tier), dihe_k_atoms.data(tier),
                            dihe_l_atoms.data(tier), dihe_parameter_indices.data(tier),
                            dihe14_parameter_indices.data(tier), dihe_modifiers.data(tier),
                            infr14_i_atoms.data(tier), infr14_l_atoms.data(tier),
                            infr14_parameter_indices.data(tier), urey_bradley_i_atoms.data(tier),
                            urey_bradley_k_atoms.data(tier),
                            urey_bradley_parameter_indices.data(tier),
                            charmm_impr_i_atoms.data(tier), charmm_impr_j_atoms.data(tier),
                            charmm_impr_k_atoms.data(tier), charmm_impr_l_atoms.data(tier),
                            charmm_impr_parameter_indices.data(tier), cmap_i_atoms.data(tier),
                            cmap_j_atoms.data(tier), cmap_k_atoms.data(tier),
                            cmap_l_atoms.data(tier), cmap_m_atoms.data(tier),
                            cmap_surface_dimensions.data(tier), cmap_surface_bounds.data(tier),
                            cmap_patch_bounds.data(tier), cmap_surface_indices.data(tier),
                            urey_bradley_stiffnesses.data(tier),
                            urey_bradley_equilibria.data(tier), charmm_impr_stiffnesses.data(tier),
                            charmm_impr_phase_angles.data(tier), cmap_surfaces.data(tier),
                            cmap_phi_derivatives.data(tier), cmap_psi_derivatives.data(tier),
                            cmap_phi_psi_derivatives.data(tier), cmap_patches.data(tier),
                            bond_assigned_atoms.data(tier), bond_assigned_index.data(tier),
                            bond_assigned_terms.data(tier), bond_assigned_bounds.data(tier),
                            angl_assigned_atoms.data(tier), angl_assigned_index.data(tier),
                            angl_assigned_terms.data(tier), angl_assigned_bounds.data(tier),
                            dihe_assigned_atoms.data(tier), dihe_assigned_index.data(tier),
                            dihe_assigned_terms.data(tier), dihe_assigned_bounds.data(tier),
                            urey_bradley_assigned_atoms.data(tier),
                            urey_bradley_assigned_index.data(tier),
                            urey_bradley_assigned_terms.data(tier),
                            urey_bradley_assigned_bounds.data(tier),
                            charmm_impr_assigned_atoms.data(tier),
                            charmm_impr_assigned_index.data(tier),
                            charmm_impr_assigned_terms.data(tier),
                            charmm_impr_assigned_bounds.data(tier), cmap_assigned_atoms.data(tier),
                            cmap_assigned_index.data(tier), cmap_assigned_terms.data(tier),
                            cmap_assigned_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
ValenceKit<float> AtomGraph::getSinglePrecisionValenceKit(const HybridTargetLevel tier) const {
  return ValenceKit<float>(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                           bond_parameter_count, angl_parameter_count, dihe_parameter_count,
                           inferred_14_attenuations, attenuated_14_type_count,
                           urey_bradley_term_count, charmm_impr_term_count, cmap_term_count,
                           urey_bradley_parameter_count, charmm_impr_parameter_count,
                           cmap_surface_count, sp_bond_stiffnesses.data(tier),
                           sp_bond_equilibria.data(tier), sp_angl_stiffnesses.data(tier),
                           sp_angl_equilibria.data(tier), sp_dihe_amplitudes.data(tier),
                           sp_dihe_periodicities.data(tier), sp_dihe_phase_angles.data(tier),
                           sp_attn14_elec_factors.data(tier), sp_attn14_vdw_factors.data(tier),
                           bond_i_atoms.data(tier), bond_j_atoms.data(tier),
                           bond_parameter_indices.data(tier), bond_modifiers.data(tier),
                           angl_i_atoms.data(tier), angl_j_atoms.data(tier),
                           angl_k_atoms.data(tier), angl_parameter_indices.data(tier),
                           angl_modifiers.data(tier), dihe_i_atoms.data(tier),
                           dihe_j_atoms.data(tier), dihe_k_atoms.data(tier),
                           dihe_l_atoms.data(tier), dihe_parameter_indices.data(tier),
                           dihe14_parameter_indices.data(tier), dihe_modifiers.data(tier),
                           infr14_i_atoms.data(tier), infr14_l_atoms.data(tier),
                           infr14_parameter_indices.data(tier), urey_bradley_i_atoms.data(tier),
                           urey_bradley_k_atoms.data(tier),
                           urey_bradley_parameter_indices.data(tier),
                           charmm_impr_i_atoms.data(tier), charmm_impr_j_atoms.data(tier),
                           charmm_impr_k_atoms.data(tier), charmm_impr_l_atoms.data(tier),
                           charmm_impr_parameter_indices.data(tier), cmap_i_atoms.data(tier),
                           cmap_j_atoms.data(tier), cmap_k_atoms.data(tier),
                           cmap_l_atoms.data(tier), cmap_m_atoms.data(tier),
                           cmap_surface_dimensions.data(tier), cmap_surface_bounds.data(tier),
                           cmap_patch_bounds.data(tier), cmap_surface_indices.data(tier),
                           sp_urey_bradley_stiffnesses.data(tier),
                           sp_urey_bradley_equilibria.data(tier),
                           sp_charmm_impr_stiffnesses.data(tier),
                           sp_charmm_impr_phase_angles.data(tier), sp_cmap_surfaces.data(tier),
                           sp_cmap_phi_derivatives.data(tier), sp_cmap_psi_derivatives.data(tier),
                           sp_cmap_phi_psi_derivatives.data(tier), sp_cmap_patches.data(tier),
                           bond_assigned_atoms.data(tier), bond_assigned_index.data(tier),
                           bond_assigned_terms.data(tier), bond_assigned_bounds.data(tier),
                           angl_assigned_atoms.data(tier), angl_assigned_index.data(tier),
                           angl_assigned_terms.data(tier), angl_assigned_bounds.data(tier),
                           dihe_assigned_atoms.data(tier), dihe_assigned_index.data(tier),
                           dihe_assigned_terms.data(tier), dihe_assigned_bounds.data(tier),
                           urey_bradley_assigned_atoms.data(tier),
                           urey_bradley_assigned_index.data(tier),
                           urey_bradley_assigned_terms.data(tier),
                           urey_bradley_assigned_bounds.data(tier),
                           charmm_impr_assigned_atoms.data(tier),
                           charmm_impr_assigned_index.data(tier),
                           charmm_impr_assigned_terms.data(tier),
                           charmm_impr_assigned_bounds.data(tier), cmap_assigned_atoms.data(tier),
                           cmap_assigned_index.data(tier), cmap_assigned_terms.data(tier),
                           cmap_assigned_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
NonbondedKit<double>
AtomGraph::getDoublePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  return NonbondedKit<double>(atom_count, charge_type_count, lj_type_count, coulomb_constant,
                              atomic_charges.data(tier), charge_indices.data(tier),
                              lennard_jones_indices.data(tier), charge_parameters.data(tier),
                              lj_a_values.data(tier), lj_b_values.data(tier),
                              lj_c_values.data(tier), lj_14_a_values.data(tier),
                              lj_14_b_values.data(tier), lj_14_c_values.data(tier),
                              lj_sigma_values.data(tier), lj_14_sigma_values.data(tier),
                              nb11_exclusion_list.data(tier), nb11_exclusion_bounds.data(tier),
                              nb12_exclusion_list.data(tier), nb12_exclusion_bounds.data(tier),
                              nb13_exclusion_list.data(tier), nb13_exclusion_bounds.data(tier),
                              nb14_exclusion_list.data(tier), nb14_exclusion_bounds.data(tier),
                              lj_type_corrections.data(tier));
}

//-------------------------------------------------------------------------------------------------
NonbondedKit<float> AtomGraph::getSinglePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  return NonbondedKit<float>(atom_count, charge_type_count, lj_type_count, coulomb_constant,
                             sp_atomic_charges.data(tier), charge_indices.data(tier),
                             lennard_jones_indices.data(tier), sp_charge_parameters.data(tier),
                             sp_lj_a_values.data(tier), sp_lj_b_values.data(tier),
                             sp_lj_c_values.data(tier), sp_lj_14_a_values.data(tier),
                             sp_lj_14_b_values.data(tier), sp_lj_14_c_values.data(tier),
                             sp_lj_sigma_values.data(tier), sp_lj_14_sigma_values.data(tier),
                             nb11_exclusion_list.data(tier), nb11_exclusion_bounds.data(tier),
                             nb12_exclusion_list.data(tier), nb12_exclusion_bounds.data(tier),
                             nb13_exclusion_list.data(tier), nb13_exclusion_bounds.data(tier),
                             nb14_exclusion_list.data(tier), nb14_exclusion_bounds.data(tier),
                             sp_lj_type_corrections.data(tier));
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventKit<double>
AtomGraph::getDoublePrecisionImplicitSolventKit(const HybridTargetLevel tier) const {
  return ImplicitSolventKit<double>(atom_count, gb_style, dielectric_constant, salt_concentration,
                                    neck_gb_indices.data(tier), atomic_pb_radii.data(tier),
                                    gb_screening_factors.data(tier),
                                    gb_alpha_parameters.data(tier), gb_beta_parameters.data(tier),
                                    gb_gamma_parameters.data(tier));
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventKit<float>
AtomGraph::getSinglePrecisionImplicitSolventKit(const HybridTargetLevel tier) const {
  return ImplicitSolventKit<float>(atom_count, gb_style, dielectric_constant, salt_concentration,
                                   neck_gb_indices.data(tier), sp_atomic_pb_radii.data(tier),
                                   sp_gb_screening_factors.data(tier),
                                   sp_gb_alpha_parameters.data(tier),
                                   sp_gb_beta_parameters.data(tier),
                                   sp_gb_gamma_parameters.data(tier));
}

//-------------------------------------------------------------------------------------------------
ChemicalDetailsKit AtomGraph::getChemicalDetailsKit(const HybridTargetLevel tier) const {
  return ChemicalDetailsKit(atom_count, residue_count, molecule_count, unconstrained_dof,
                            constrained_dof, atom_names.data(tier), residue_names.data(tier),
                            atom_types.data(tier), atomic_numbers.data(tier),
                            residue_limits.data(tier), atom_struc_numbers.data(tier),
                            residue_numbers.data(tier), molecule_membership.data(tier),
                            molecule_contents.data(tier), molecule_limits.data(tier),
                            atomic_masses.data(tier), sp_atomic_masses.data(tier),
                            inverse_atomic_masses.data(tier), sp_inverse_atomic_masses.data(tier));
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKit<double>
AtomGraph::getDoublePrecisionVirtualSiteKit(const HybridTargetLevel tier) const {
  return VirtualSiteKit<double>(virtual_site_count, virtual_site_parameter_set_count,
                                virtual_site_atoms.data(tier),
                                virtual_site_frame1_atoms.data(tier),
                                virtual_site_frame2_atoms.data(tier),
                                virtual_site_frame3_atoms.data(tier),
                                virtual_site_frame4_atoms.data(tier),
                                virtual_site_parameter_indices.data(tier),
                                virtual_site_frame_types.data(tier),
                                virtual_site_frame_dim1.data(tier),
                                virtual_site_frame_dim2.data(tier),
                                virtual_site_frame_dim3.data(tier));
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKit<float>
AtomGraph::getSinglePrecisionVirtualSiteKit(const HybridTargetLevel tier) const {
  return VirtualSiteKit<float>(virtual_site_count, virtual_site_parameter_set_count,
                               virtual_site_atoms.data(tier),
                               virtual_site_frame1_atoms.data(tier),
                               virtual_site_frame2_atoms.data(tier),
                               virtual_site_frame3_atoms.data(tier),
                               virtual_site_frame4_atoms.data(tier),
                               virtual_site_parameter_indices.data(tier),
                               virtual_site_frame_types.data(tier),
                               sp_virtual_site_frame_dim1.data(tier),
                               sp_virtual_site_frame_dim2.data(tier),
                               sp_virtual_site_frame_dim3.data(tier));
}

//-------------------------------------------------------------------------------------------------
ConstraintKit<double>
AtomGraph::getDoublePrecisionConstraintKit(const HybridTargetLevel tier) const {
  return ConstraintKit<double>(settle_group_count, settle_parameter_count, constraint_group_count,
                               constraint_parameter_count, settle_oxygen_atoms.data(tier),
                               settle_hydro1_atoms.data(tier), settle_hydro2_atoms.data(tier),
                               settle_parameter_indices.data(tier),
                               constraint_group_atoms.data(tier),
                               constraint_group_bounds.data(tier),
                               constraint_parameter_indices.data(tier),
                               constraint_parameter_bounds.data(tier),
                               settle_mormt.data(tier), settle_mhrmt.data(tier),
                               settle_ra.data(tier), settle_rb.data(tier), settle_rc.data(tier),
                               settle_invra.data(tier), constraint_squared_lengths.data(tier),
                               constraint_inverse_masses.data(tier));
}

//-------------------------------------------------------------------------------------------------
ConstraintKit<float>
AtomGraph::getSinglePrecisionConstraintKit(const HybridTargetLevel tier) const {
  return ConstraintKit<float>(settle_group_count, settle_parameter_count, constraint_group_count,
                              constraint_parameter_count, settle_oxygen_atoms.data(tier),
                              settle_hydro1_atoms.data(tier), settle_hydro2_atoms.data(tier),
                              settle_parameter_indices.data(tier),
                              constraint_group_atoms.data(tier),
                              constraint_group_bounds.data(tier),
                              constraint_parameter_indices.data(tier),
                              constraint_parameter_bounds.data(tier),
                              sp_settle_mormt.data(tier), sp_settle_mhrmt.data(tier),
                              sp_settle_ra.data(tier), sp_settle_rb.data(tier),
                              sp_settle_rc.data(tier), sp_settle_invra.data(tier),
                              sp_constraint_squared_lengths.data(tier),
                              sp_constraint_inverse_masses.data(tier));
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
ValenceKit<double> AtomGraph::getDeviceViewToHostDPValenceKit() const {
  return ValenceKit<double>(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                            bond_parameter_count, angl_parameter_count, dihe_parameter_count,
                            inferred_14_attenuations, attenuated_14_type_count,
                            urey_bradley_term_count, charmm_impr_term_count, cmap_term_count,
                            urey_bradley_parameter_count, charmm_impr_parameter_count,
                            cmap_surface_count, bond_stiffnesses.getDeviceValidHostPointer(),
                            bond_equilibria.getDeviceValidHostPointer(),
                            angl_stiffnesses.getDeviceValidHostPointer(),
                            angl_equilibria.getDeviceValidHostPointer(),
                            dihe_amplitudes.getDeviceValidHostPointer(),
                            dihe_periodicities.getDeviceValidHostPointer(),
                            dihe_phase_angles.getDeviceValidHostPointer(),
                            attn14_elec_factors.getDeviceValidHostPointer(),
                            attn14_vdw_factors.getDeviceValidHostPointer(),
                            bond_i_atoms.getDeviceValidHostPointer(),
                            bond_j_atoms.getDeviceValidHostPointer(),
                            bond_parameter_indices.getDeviceValidHostPointer(),
                            bond_modifiers.getDeviceValidHostPointer(),
                            angl_i_atoms.getDeviceValidHostPointer(),
                            angl_j_atoms.getDeviceValidHostPointer(),
                            angl_k_atoms.getDeviceValidHostPointer(),
                            angl_parameter_indices.getDeviceValidHostPointer(),
                            angl_modifiers.getDeviceValidHostPointer(),
                            dihe_i_atoms.getDeviceValidHostPointer(),
                            dihe_j_atoms.getDeviceValidHostPointer(),
                            dihe_k_atoms.getDeviceValidHostPointer(),
                            dihe_l_atoms.getDeviceValidHostPointer(),
                            dihe_parameter_indices.getDeviceValidHostPointer(),
                            dihe14_parameter_indices.getDeviceValidHostPointer(),
                            dihe_modifiers.getDeviceValidHostPointer(),
                            infr14_i_atoms.getDeviceValidHostPointer(),
                            infr14_l_atoms.getDeviceValidHostPointer(),
                            infr14_parameter_indices.getDeviceValidHostPointer(),
                            urey_bradley_i_atoms.getDeviceValidHostPointer(),
                            urey_bradley_k_atoms.getDeviceValidHostPointer(),
                            urey_bradley_parameter_indices.getDeviceValidHostPointer(),
                            charmm_impr_i_atoms.getDeviceValidHostPointer(),
                            charmm_impr_j_atoms.getDeviceValidHostPointer(),
                            charmm_impr_k_atoms.getDeviceValidHostPointer(),
                            charmm_impr_l_atoms.getDeviceValidHostPointer(),
                            charmm_impr_parameter_indices.getDeviceValidHostPointer(),
                            cmap_i_atoms.getDeviceValidHostPointer(),
                            cmap_j_atoms.getDeviceValidHostPointer(),
                            cmap_k_atoms.getDeviceValidHostPointer(),
                            cmap_l_atoms.getDeviceValidHostPointer(),
                            cmap_m_atoms.getDeviceValidHostPointer(),
                            cmap_surface_dimensions.getDeviceValidHostPointer(),
                            cmap_surface_bounds.getDeviceValidHostPointer(),
                            cmap_patch_bounds.getDeviceValidHostPointer(),
                            cmap_surface_indices.getDeviceValidHostPointer(),
                            urey_bradley_stiffnesses.getDeviceValidHostPointer(),
                            urey_bradley_equilibria.getDeviceValidHostPointer(),
                            charmm_impr_stiffnesses.getDeviceValidHostPointer(),
                            charmm_impr_phase_angles.getDeviceValidHostPointer(),
                            cmap_surfaces.getDeviceValidHostPointer(),
                            cmap_phi_derivatives.getDeviceValidHostPointer(),
                            cmap_psi_derivatives.getDeviceValidHostPointer(),
                            cmap_phi_psi_derivatives.getDeviceValidHostPointer(),
                            cmap_patches.getDeviceValidHostPointer(),
                            bond_assigned_atoms.getDeviceValidHostPointer(),
                            bond_assigned_index.getDeviceValidHostPointer(),
                            bond_assigned_terms.getDeviceValidHostPointer(),
                            bond_assigned_bounds.getDeviceValidHostPointer(),
                            angl_assigned_atoms.getDeviceValidHostPointer(),
                            angl_assigned_index.getDeviceValidHostPointer(),
                            angl_assigned_terms.getDeviceValidHostPointer(),
                            angl_assigned_bounds.getDeviceValidHostPointer(),
                            dihe_assigned_atoms.getDeviceValidHostPointer(),
                            dihe_assigned_index.getDeviceValidHostPointer(),
                            dihe_assigned_terms.getDeviceValidHostPointer(),
                            dihe_assigned_bounds.getDeviceValidHostPointer(),
                            urey_bradley_assigned_atoms.getDeviceValidHostPointer(),
                            urey_bradley_assigned_index.getDeviceValidHostPointer(),
                            urey_bradley_assigned_terms.getDeviceValidHostPointer(),
                            urey_bradley_assigned_bounds.getDeviceValidHostPointer(),
                            charmm_impr_assigned_atoms.getDeviceValidHostPointer(),
                            charmm_impr_assigned_index.getDeviceValidHostPointer(),
                            charmm_impr_assigned_terms.getDeviceValidHostPointer(),
                            charmm_impr_assigned_bounds.getDeviceValidHostPointer(),
                            cmap_assigned_atoms.getDeviceValidHostPointer(),
                            cmap_assigned_index.getDeviceValidHostPointer(),
                            cmap_assigned_terms.getDeviceValidHostPointer(),
                            cmap_assigned_bounds.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
ValenceKit<float> AtomGraph::getDeviceViewToHostSPValenceKit() const {
  return ValenceKit<float>(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                           bond_parameter_count, angl_parameter_count, dihe_parameter_count,
                           inferred_14_attenuations, attenuated_14_type_count,
                           urey_bradley_term_count, charmm_impr_term_count, cmap_term_count,
                           urey_bradley_parameter_count, charmm_impr_parameter_count,
                           cmap_surface_count, sp_bond_stiffnesses.getDeviceValidHostPointer(),
                           sp_bond_equilibria.getDeviceValidHostPointer(),
                           sp_angl_stiffnesses.getDeviceValidHostPointer(),
                           sp_angl_equilibria.getDeviceValidHostPointer(),
                           sp_dihe_amplitudes.getDeviceValidHostPointer(),
                           sp_dihe_periodicities.getDeviceValidHostPointer(),
                           sp_dihe_phase_angles.getDeviceValidHostPointer(),
                           sp_attn14_elec_factors.getDeviceValidHostPointer(),
                           sp_attn14_vdw_factors.getDeviceValidHostPointer(),
                           bond_i_atoms.getDeviceValidHostPointer(),
                           bond_j_atoms.getDeviceValidHostPointer(),
                           bond_parameter_indices.getDeviceValidHostPointer(),
                           bond_modifiers.getDeviceValidHostPointer(),
                           angl_i_atoms.getDeviceValidHostPointer(),
                           angl_j_atoms.getDeviceValidHostPointer(),
                           angl_k_atoms.getDeviceValidHostPointer(),
                           angl_parameter_indices.getDeviceValidHostPointer(),
                           angl_modifiers.getDeviceValidHostPointer(),
                           dihe_i_atoms.getDeviceValidHostPointer(),
                           dihe_j_atoms.getDeviceValidHostPointer(),
                           dihe_k_atoms.getDeviceValidHostPointer(),
                           dihe_l_atoms.getDeviceValidHostPointer(),
                           dihe_parameter_indices.getDeviceValidHostPointer(),
                           dihe14_parameter_indices.getDeviceValidHostPointer(),
                           dihe_modifiers.getDeviceValidHostPointer(),
                           infr14_i_atoms.getDeviceValidHostPointer(),
                           infr14_l_atoms.getDeviceValidHostPointer(),
                           infr14_parameter_indices.getDeviceValidHostPointer(),
                           urey_bradley_i_atoms.getDeviceValidHostPointer(),
                           urey_bradley_k_atoms.getDeviceValidHostPointer(),
                           urey_bradley_parameter_indices.getDeviceValidHostPointer(),
                           charmm_impr_i_atoms.getDeviceValidHostPointer(),
                           charmm_impr_j_atoms.getDeviceValidHostPointer(),
                           charmm_impr_k_atoms.getDeviceValidHostPointer(),
                           charmm_impr_l_atoms.getDeviceValidHostPointer(),
                           charmm_impr_parameter_indices.getDeviceValidHostPointer(),
                           cmap_i_atoms.getDeviceValidHostPointer(),
                           cmap_j_atoms.getDeviceValidHostPointer(),
                           cmap_k_atoms.getDeviceValidHostPointer(),
                           cmap_l_atoms.getDeviceValidHostPointer(),
                           cmap_m_atoms.getDeviceValidHostPointer(),
                           cmap_surface_dimensions.getDeviceValidHostPointer(),
                           cmap_surface_bounds.getDeviceValidHostPointer(),
                           cmap_patch_bounds.getDeviceValidHostPointer(),
                           cmap_surface_indices.getDeviceValidHostPointer(),
                           sp_urey_bradley_stiffnesses.getDeviceValidHostPointer(),
                           sp_urey_bradley_equilibria.getDeviceValidHostPointer(),
                           sp_charmm_impr_stiffnesses.getDeviceValidHostPointer(),
                           sp_charmm_impr_phase_angles.getDeviceValidHostPointer(),
                           sp_cmap_surfaces.getDeviceValidHostPointer(),
                           sp_cmap_phi_derivatives.getDeviceValidHostPointer(),
                           sp_cmap_psi_derivatives.getDeviceValidHostPointer(),
                           sp_cmap_phi_psi_derivatives.getDeviceValidHostPointer(),
                           sp_cmap_patches.getDeviceValidHostPointer(),
                           bond_assigned_atoms.getDeviceValidHostPointer(),
                           bond_assigned_index.getDeviceValidHostPointer(),
                           bond_assigned_terms.getDeviceValidHostPointer(),
                           bond_assigned_bounds.getDeviceValidHostPointer(),
                           angl_assigned_atoms.getDeviceValidHostPointer(),
                           angl_assigned_index.getDeviceValidHostPointer(),
                           angl_assigned_terms.getDeviceValidHostPointer(),
                           angl_assigned_bounds.getDeviceValidHostPointer(),
                           dihe_assigned_atoms.getDeviceValidHostPointer(),
                           dihe_assigned_index.getDeviceValidHostPointer(),
                           dihe_assigned_terms.getDeviceValidHostPointer(),
                           dihe_assigned_bounds.getDeviceValidHostPointer(),
                           urey_bradley_assigned_atoms.getDeviceValidHostPointer(),
                           urey_bradley_assigned_index.getDeviceValidHostPointer(),
                           urey_bradley_assigned_terms.getDeviceValidHostPointer(),
                           urey_bradley_assigned_bounds.getDeviceValidHostPointer(),
                           charmm_impr_assigned_atoms.getDeviceValidHostPointer(),
                           charmm_impr_assigned_index.getDeviceValidHostPointer(),
                           charmm_impr_assigned_terms.getDeviceValidHostPointer(),
                           charmm_impr_assigned_bounds.getDeviceValidHostPointer(),
                           cmap_assigned_atoms.getDeviceValidHostPointer(),
                           cmap_assigned_index.getDeviceValidHostPointer(),
                           cmap_assigned_terms.getDeviceValidHostPointer(),
                           cmap_assigned_bounds.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
NonbondedKit<double>
AtomGraph::getDeviceViewToHostDPNonbondedKit() const {
  return NonbondedKit<double>(atom_count, charge_type_count, lj_type_count, coulomb_constant,
                              atomic_charges.getDeviceValidHostPointer(),
                              charge_indices.getDeviceValidHostPointer(),
                              lennard_jones_indices.getDeviceValidHostPointer(),
                              charge_parameters.getDeviceValidHostPointer(),
                              lj_a_values.getDeviceValidHostPointer(),
                              lj_b_values.getDeviceValidHostPointer(),
                              lj_c_values.getDeviceValidHostPointer(),
                              lj_14_a_values.getDeviceValidHostPointer(),
                              lj_14_b_values.getDeviceValidHostPointer(),
                              lj_14_c_values.getDeviceValidHostPointer(),
                              lj_sigma_values.getDeviceValidHostPointer(),
                              lj_14_sigma_values.getDeviceValidHostPointer(),
                              nb11_exclusion_list.getDeviceValidHostPointer(),
                              nb11_exclusion_bounds.getDeviceValidHostPointer(),
                              nb12_exclusion_list.getDeviceValidHostPointer(),
                              nb12_exclusion_bounds.getDeviceValidHostPointer(),
                              nb13_exclusion_list.getDeviceValidHostPointer(),
                              nb13_exclusion_bounds.getDeviceValidHostPointer(),
                              nb14_exclusion_list.getDeviceValidHostPointer(),
                              nb14_exclusion_bounds.getDeviceValidHostPointer(),
                              lj_type_corrections.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
NonbondedKit<float> AtomGraph::getDeviceViewToHostSPNonbondedKit() const {
  return NonbondedKit<float>(atom_count, charge_type_count, lj_type_count, coulomb_constant,
                             sp_atomic_charges.getDeviceValidHostPointer(),
                             charge_indices.getDeviceValidHostPointer(),
                             lennard_jones_indices.getDeviceValidHostPointer(),
                             sp_charge_parameters.getDeviceValidHostPointer(),
                             sp_lj_a_values.getDeviceValidHostPointer(),
                             sp_lj_b_values.getDeviceValidHostPointer(),
                             sp_lj_c_values.getDeviceValidHostPointer(),
                             sp_lj_14_a_values.getDeviceValidHostPointer(),
                             sp_lj_14_b_values.getDeviceValidHostPointer(),
                             sp_lj_14_c_values.getDeviceValidHostPointer(),
                             sp_lj_sigma_values.getDeviceValidHostPointer(),
                             sp_lj_14_sigma_values.getDeviceValidHostPointer(),
                             nb11_exclusion_list.getDeviceValidHostPointer(),
                             nb11_exclusion_bounds.getDeviceValidHostPointer(),
                             nb12_exclusion_list.getDeviceValidHostPointer(),
                             nb12_exclusion_bounds.getDeviceValidHostPointer(),
                             nb13_exclusion_list.getDeviceValidHostPointer(),
                             nb13_exclusion_bounds.getDeviceValidHostPointer(),
                             nb14_exclusion_list.getDeviceValidHostPointer(),
                             nb14_exclusion_bounds.getDeviceValidHostPointer(),
                             sp_lj_type_corrections.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventKit<double>
AtomGraph::getDeviceViewToHostDPImplicitSolventKit() const {
  return ImplicitSolventKit<double>(atom_count, gb_style, dielectric_constant, salt_concentration,
                                    neck_gb_indices.getDeviceValidHostPointer(),
                                    atomic_pb_radii.getDeviceValidHostPointer(),
                                    gb_screening_factors.getDeviceValidHostPointer(),
                                    gb_alpha_parameters.getDeviceValidHostPointer(),
                                    gb_beta_parameters.getDeviceValidHostPointer(),
                                    gb_gamma_parameters.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventKit<float>
AtomGraph::getDeviceViewToHostSPImplicitSolventKit() const {
  return ImplicitSolventKit<float>(atom_count, gb_style, dielectric_constant, salt_concentration,
                                   neck_gb_indices.getDeviceValidHostPointer(),
                                   sp_atomic_pb_radii.getDeviceValidHostPointer(),
                                   sp_gb_screening_factors.getDeviceValidHostPointer(),
                                   sp_gb_alpha_parameters.getDeviceValidHostPointer(),
                                   sp_gb_beta_parameters.getDeviceValidHostPointer(),
                                   sp_gb_gamma_parameters.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
ChemicalDetailsKit AtomGraph::getDeviceViewToHostChemicalDetailsKit() const {
  return ChemicalDetailsKit(atom_count, residue_count, molecule_count,
                            unconstrained_dof, constrained_dof,
                            atom_names.getDeviceValidHostPointer(),
                            residue_names.getDeviceValidHostPointer(),
                            atom_types.getDeviceValidHostPointer(),
                            atomic_numbers.getDeviceValidHostPointer(),
                            residue_limits.getDeviceValidHostPointer(),
                            atom_struc_numbers.getDeviceValidHostPointer(),
                            residue_numbers.getDeviceValidHostPointer(),
                            molecule_membership.getDeviceValidHostPointer(),
                            molecule_contents.getDeviceValidHostPointer(),
                            molecule_limits.getDeviceValidHostPointer(),
                            atomic_masses.getDeviceValidHostPointer(),
                            sp_atomic_masses.getDeviceValidHostPointer(),
                            inverse_atomic_masses.getDeviceValidHostPointer(),
                            sp_inverse_atomic_masses.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKit<double>
AtomGraph::getDeviceViewToHostDPVirtualSiteKit() const {
  return VirtualSiteKit<double>(virtual_site_count, virtual_site_parameter_set_count,
                                virtual_site_atoms.getDeviceValidHostPointer(),
                                virtual_site_frame1_atoms.getDeviceValidHostPointer(),
                                virtual_site_frame2_atoms.getDeviceValidHostPointer(),
                                virtual_site_frame3_atoms.getDeviceValidHostPointer(),
                                virtual_site_frame4_atoms.getDeviceValidHostPointer(),
                                virtual_site_parameter_indices.getDeviceValidHostPointer(),
                                virtual_site_frame_types.getDeviceValidHostPointer(),
                                virtual_site_frame_dim1.getDeviceValidHostPointer(),
                                virtual_site_frame_dim2.getDeviceValidHostPointer(),
                                virtual_site_frame_dim3.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKit<float>
AtomGraph::getDeviceViewToHostSPVirtualSiteKit() const {
  return VirtualSiteKit<float>(virtual_site_count, virtual_site_parameter_set_count,
                               virtual_site_atoms.getDeviceValidHostPointer(),
                               virtual_site_frame1_atoms.getDeviceValidHostPointer(),
                               virtual_site_frame2_atoms.getDeviceValidHostPointer(),
                               virtual_site_frame3_atoms.getDeviceValidHostPointer(),
                               virtual_site_frame4_atoms.getDeviceValidHostPointer(),
                               virtual_site_parameter_indices.getDeviceValidHostPointer(),
                               virtual_site_frame_types.getDeviceValidHostPointer(),
                               sp_virtual_site_frame_dim1.getDeviceValidHostPointer(),
                               sp_virtual_site_frame_dim2.getDeviceValidHostPointer(),
                               sp_virtual_site_frame_dim3.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
ConstraintKit<double>
AtomGraph::getDeviceViewToHostDPConstraintKit() const {
  return ConstraintKit<double>(settle_group_count, settle_parameter_count, constraint_group_count,
                               constraint_parameter_count,
                               settle_oxygen_atoms.getDeviceValidHostPointer(),
                               settle_hydro1_atoms.getDeviceValidHostPointer(),
                               settle_hydro2_atoms.getDeviceValidHostPointer(),
                               settle_parameter_indices.getDeviceValidHostPointer(),
                               constraint_group_atoms.getDeviceValidHostPointer(),
                               constraint_group_bounds.getDeviceValidHostPointer(),
                               constraint_parameter_indices.getDeviceValidHostPointer(),
                               constraint_parameter_bounds.getDeviceValidHostPointer(),
                               settle_mormt.getDeviceValidHostPointer(),
                               settle_mhrmt.getDeviceValidHostPointer(),
                               settle_ra.getDeviceValidHostPointer(),
                               settle_rb.getDeviceValidHostPointer(),
                               settle_rc.getDeviceValidHostPointer(),
                               settle_invra.getDeviceValidHostPointer(),
                               constraint_squared_lengths.getDeviceValidHostPointer(),
                               constraint_inverse_masses.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
ConstraintKit<float>
AtomGraph::getDeviceViewToHostSPConstraintKit() const {
  return ConstraintKit<float>(settle_group_count, settle_parameter_count, constraint_group_count,
                              constraint_parameter_count,
                              settle_oxygen_atoms.getDeviceValidHostPointer(),
                              settle_hydro1_atoms.getDeviceValidHostPointer(),
                              settle_hydro2_atoms.getDeviceValidHostPointer(),
                              settle_parameter_indices.getDeviceValidHostPointer(),
                              constraint_group_atoms.getDeviceValidHostPointer(),
                              constraint_group_bounds.getDeviceValidHostPointer(),
                              constraint_parameter_indices.getDeviceValidHostPointer(),
                              constraint_parameter_bounds.getDeviceValidHostPointer(),
                              sp_settle_mormt.getDeviceValidHostPointer(),
                              sp_settle_mhrmt.getDeviceValidHostPointer(),
                              sp_settle_ra.getDeviceValidHostPointer(),
                              sp_settle_rb.getDeviceValidHostPointer(),
                              sp_settle_rc.getDeviceValidHostPointer(),
                              sp_settle_invra.getDeviceValidHostPointer(),
                              sp_constraint_squared_lengths.getDeviceValidHostPointer(),
                              sp_constraint_inverse_masses.getDeviceValidHostPointer());
}
#  endif
#endif

//-------------------------------------------------------------------------------------------------
const AtomGraph* AtomGraph::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
AtomGraph* AtomGraph::getSelfPointer() {
  return this;
}

} // namespace topology
} // namespace stormm
