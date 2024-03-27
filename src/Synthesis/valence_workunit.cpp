#include <algorithm>
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Math/statistics.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Topology/atomgraph_enumerators.h"
#include "Topology/topology_util.h"
#include "valence_workunit.h"

namespace stormm {
namespace synthesis {

using constants::small_block_size;
using stmath::accumulateBitmask;
using stmath::DataOrder;
using stmath::locateValue;
using stmath::maxValue;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::reduceUniqueValues;
using stmath::roundUp;
using stmath::sum;
using parse::char4ToString;
using topology::extractBoundedListEntries;
using topology::markAffectorAtoms;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using topology::TorsionKind;
using topology::VirtualSiteKind;
using topology::writeAtomList;
  
//-------------------------------------------------------------------------------------------------
ValenceDelegator::ValenceDelegator(const AtomGraph *ag_in, const RestraintApparatus *ra_in) :
    atom_count{ag_in->getAtomCount()},
    first_unassigned_atom{0}, max_presence_allocation{2},
    bond_affector_list{}, bond_affector_bounds{}, angl_affector_list{}, angl_affector_bounds{},
    dihe_affector_list{}, dihe_affector_bounds{}, ubrd_affector_list{}, ubrd_affector_bounds{},
    cimp_affector_list{}, cimp_affector_bounds{}, cmap_affector_list{}, cmap_affector_bounds{},
    infr_affector_list{}, infr_affector_bounds{}, vste_affector_list{}, vste_affector_bounds{},
    cnst_affector_list{}, cnst_affector_bounds{}, sett_affector_list{}, sett_affector_bounds{},
    work_unit_assignment_count{}, work_unit_presence{}, assigned_update_work_units{},
    assigned_bond_acc_work_units{}, assigned_angl_acc_work_units{}, assigned_dihe_acc_work_units{},
    assigned_ubrd_acc_work_units{}, assigned_cimp_acc_work_units{}, assigned_cmap_acc_work_units{},
    assigned_infr14_acc_work_units{}, assigned_rposn_acc_work_units{},
    assigned_rbond_acc_work_units{}, assigned_rangl_acc_work_units{},
    assigned_rdihe_acc_work_units{}, ag_pointer{ag_in}, ra_pointer{ra_in}
{
  // Get relevant abstracts
  const ValenceKit<double> vk = ag_in->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_in->getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag_in->getDoublePrecisionConstraintKit();
  const RestraintKit<double, double2, double4> rar = ra_in->dpData();
  
  // Allocate and fill the arrays
  allocate();
  fillAffectorArrays(vk, vsk, cnk, rar);
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getAtomAssignmentCount(const int atom_index) const {
  return work_unit_assignment_count[atom_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getFirstUnassignedAtom() const {
  return first_unassigned_atom;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getBondAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(bond_affector_list, bond_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getAngleAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(angl_affector_list, angl_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getDihedralAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(dihe_affector_list, dihe_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getUreyBradleyAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(ubrd_affector_list, ubrd_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getCharmmImproperAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(cimp_affector_list, cimp_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getCmapAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(cmap_affector_list, cmap_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getInferred14Affectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(infr_affector_list, infr_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getPositionalRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rposn_affector_list, rposn_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getDistanceRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rbond_affector_list, rbond_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getAngleRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rangl_affector_list, rangl_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getDihedralRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rdihe_affector_list, rdihe_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getVirtualSiteAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(vste_affector_list, vste_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getSettleGroupAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(sett_affector_list, sett_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getConstraintGroupAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(cnst_affector_list, cnst_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::findMovementPartners(const int atom_idx,
                                       const std::vector<int> &caller_stack) const {
  std::vector<int> result;
  result.reserve(64);
  const ConstraintKit<double> cnk  = ag_pointer->getDoublePrecisionConstraintKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  
  // Form a new call stack from the input.  Check that the call stack has not grown too large,
  // as this might indicate an infinite loop and some nonsensical feature of the topology.
  std::vector<int> updated_caller_stack(caller_stack);
  updated_caller_stack.push_back(atom_idx);
  if (updated_caller_stack.size() > max_atom_search_stacks) {
    const std::string atom_list = writeAtomList(updated_caller_stack,
                                                ag_pointer->getChemicalDetailsKit());
    rtErr("Call stack for finding movement partners is too deep.  This likely indicates that, "
          "somehow, a constraint group and virtual site in the topology are interacting in a "
          "way that creates a paradox.  Calls include atoms " + atom_list + ".",
          "findMovementPartners");
  }

  // The movement partners must necessarily include the atom itself, even if there are no
  // constraint or virtual sites affecting it.
  result.push_back(atom_idx);

  // If the atom is a virtual site, then all of its frame atoms are likewise required to move
  // before its new position can be determined.  Furthermore, any constraint groups which they
  // are a part of must be included.  Use a recursive call to sweep up the frame atoms'
  // constraint groups.
  if (ag_pointer->getAtomicNumber(atom_idx) == 0) {

    // Find the precise virtual site frame index to include only those atoms
    const std::vector<int> relevant_vstes = extractBoundedListEntries(vste_affector_list,
                                                                      vste_affector_bounds,
                                                                      atom_idx);
    const int nvsite = relevant_vstes.size();
    for (int i = 0; i < nvsite; i++) {
      const int vste_idx = relevant_vstes[i];
      if (vsk.vs_atoms[vste_idx] == atom_idx) {
        const std::vector<int> tmp1 = findMovementPartners(vsk.frame1_idx[vste_idx],
                                                           updated_caller_stack);
        result.insert(result.end(), tmp1.begin(), tmp1.end());
        const std::vector<int> tmp2 = findMovementPartners(vsk.frame2_idx[vste_idx],
                                                           updated_caller_stack);
        result.insert(result.end(), tmp2.begin(), tmp2.end());
        switch (static_cast<VirtualSiteKind>(vsk.vs_types[vste_idx])) {
        case VirtualSiteKind::FLEX_2:
        case VirtualSiteKind::FIXED_2:
        case VirtualSiteKind::NONE:
          break;
        case VirtualSiteKind::FLEX_3:
        case VirtualSiteKind::FIXED_3:
        case VirtualSiteKind::FAD_3:
        case VirtualSiteKind::OUT_3:
          {
            const std::vector<int> tmp3 = findMovementPartners(vsk.frame3_idx[vste_idx],
                                                               updated_caller_stack);
            result.insert(result.end(), tmp3.begin(), tmp3.end());
          }
          break;
        case VirtualSiteKind::FIXED_4:
          {
            const std::vector<int> tmp3 = findMovementPartners(vsk.frame3_idx[vste_idx],
                                                               updated_caller_stack);
            result.insert(result.end(), tmp3.begin(), tmp3.end());
            const std::vector<int> tmp4 = findMovementPartners(vsk.frame4_idx[vste_idx],
                                                               updated_caller_stack);
            result.insert(result.end(), tmp4.begin(), tmp4.end());
          }
          break;
        }
      }
    }
  }
  
  // Any constraint groups that affect the atom must have all of their atoms moved along
  // with the atom.
  const std::vector<int> relevant_setts = extractBoundedListEntries(sett_affector_list,
                                                                    sett_affector_bounds,
                                                                    atom_idx);
  const int nsett = relevant_setts.size();
  for (int i = 0; i < nsett; i++) {
    const int sett_index = relevant_setts[i];
    result.push_back(cnk.settle_ox_atoms[sett_index]);
    result.push_back(cnk.settle_h1_atoms[sett_index]);
    result.push_back(cnk.settle_h2_atoms[sett_index]);
  }
  const std::vector<int> relevant_cnsts = extractBoundedListEntries(cnst_affector_list,
                                                                    cnst_affector_bounds,
                                                                    atom_idx);
  const int ncnst = relevant_cnsts.size();
  for (int i = 0; i < ncnst; i++) {
    const int cnst_index = relevant_cnsts[i];
    for (int j = cnk.group_bounds[cnst_index]; j < cnk.group_bounds[cnst_index + 1]; j++) {
      result.push_back(cnk.group_list[j]);
    }
  }

  // Sort the list of required atoms and prune duplicate entries
  reduceUniqueValues(&result);  
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::findForcePartners(const int atom_idx,
                                                     const std::vector<int> &caller_stack) const {
  std::vector<int> result;
  result.reserve(32);
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();
  
  // Form a new call stack from the input.  Check that the call stack has not grown too large,
  // as this might indicate an infinite loop and some nonsensical feature of the topology.
  std::vector<int> updated_caller_stack(caller_stack);
  updated_caller_stack.push_back(atom_idx);
  if (updated_caller_stack.size() > max_atom_search_stacks) {
    const std::string atom_list = writeAtomList(updated_caller_stack,
                                                ag_pointer->getChemicalDetailsKit());
    rtErr("Call stack for finding force partners is too deep.  This likely indicates that, "
          "somehow, a virtual site in the topology is including itself as one of its frame atoms "
          "or some similar paradox.  Calls include atoms " + atom_list + ".", "findForcePartners");
  }

  // Push the atom itself onto the list first.  Even if there are no interactions, having the
  // atom is, by definition, critical to evaluating its movement.
  result.push_back(atom_idx);

  // Push atoms relevant to all valence force field terms onto the list
  std::vector<int> afl_extractions(32);
  extractBoundedListEntries(&afl_extractions, bond_affector_list, bond_affector_bounds, atom_idx);
  const int nbond = afl_extractions.size();
  for (int i = 0; i < nbond; i++) {
    const int bond_term_index = afl_extractions[i];
    result.push_back(vk.bond_i_atoms[bond_term_index]);
    result.push_back(vk.bond_j_atoms[bond_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, angl_affector_list, angl_affector_bounds, atom_idx);
  const int nangl = afl_extractions.size();
  for (int i = 0; i < nangl; i++) {
    const int angl_term_index = afl_extractions[i];
    result.push_back(vk.angl_i_atoms[angl_term_index]);
    result.push_back(vk.angl_j_atoms[angl_term_index]);
    result.push_back(vk.angl_k_atoms[angl_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, dihe_affector_list, dihe_affector_bounds, atom_idx);
  const int ndihe = afl_extractions.size();
  for (int i = 0; i < ndihe; i++) {
    const int dihe_term_index = afl_extractions[i];
    result.push_back(vk.dihe_i_atoms[dihe_term_index]);
    result.push_back(vk.dihe_j_atoms[dihe_term_index]);
    result.push_back(vk.dihe_k_atoms[dihe_term_index]);
    result.push_back(vk.dihe_l_atoms[dihe_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, ubrd_affector_list, ubrd_affector_bounds, atom_idx);
  const int nubrd = afl_extractions.size();
  for (int i = 0; i < nubrd; i++) {
    const int ubrd_term_index = afl_extractions[i];
    result.push_back(vk.ubrd_i_atoms[ubrd_term_index]);
    result.push_back(vk.ubrd_k_atoms[ubrd_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, cimp_affector_list, cimp_affector_bounds, atom_idx);
  const int ncimp = afl_extractions.size();
  for (int i = 0; i < ncimp; i++) {
    const int cimp_term_index = afl_extractions[i];
    result.push_back(vk.cimp_i_atoms[cimp_term_index]);
    result.push_back(vk.cimp_j_atoms[cimp_term_index]);
    result.push_back(vk.cimp_k_atoms[cimp_term_index]);
    result.push_back(vk.cimp_l_atoms[cimp_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, cmap_affector_list, cmap_affector_bounds, atom_idx);
  const int ncmap = afl_extractions.size();
  for (int i = 0; i < ncmap; i++) {
    const int cmap_term_index = afl_extractions[i];
    result.push_back(vk.cmap_i_atoms[cmap_term_index]);
    result.push_back(vk.cmap_j_atoms[cmap_term_index]);
    result.push_back(vk.cmap_k_atoms[cmap_term_index]);
    result.push_back(vk.cmap_l_atoms[cmap_term_index]);
    result.push_back(vk.cmap_m_atoms[cmap_term_index]);
  }

  // Most 1:4 attenuated interactions are implicitly covered by dihedral interactions to which
  // they can be linked.  Most virtual sites will have additional 1:4 attenuated interactions
  // that are not covered by dihedral interactions, as the virtual sites do not typically
  // participate in dihedral terms.  Push these onto the list.
  extractBoundedListEntries(&afl_extractions, infr_affector_list, infr_affector_bounds, atom_idx);
  const int ninfr = afl_extractions.size();
  for (int i = 0; i < ninfr; i++) {
    const int infr_term_index = afl_extractions[i];
    result.push_back(vk.infr14_i_atoms[infr_term_index]);
    result.push_back(vk.infr14_l_atoms[infr_term_index]);
  }

  // Push atoms relevant to NMR restraints onto the list.
  extractBoundedListEntries(&afl_extractions, rposn_affector_list, rposn_affector_bounds,
                            atom_idx);
  const int nrposn = afl_extractions.size();
  for (int i = 0; i < nrposn; i++) {
    const int rposn_term_index = afl_extractions[i];
    result.push_back(rar.rposn_atoms[rposn_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, rbond_affector_list, rbond_affector_bounds,
                            atom_idx);
  const int nrbond = afl_extractions.size();
  for (int i = 0; i < nrbond; i++) {
    const int rbond_term_index = afl_extractions[i];
    result.push_back(rar.rbond_i_atoms[rbond_term_index]);
    result.push_back(rar.rbond_j_atoms[rbond_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, rangl_affector_list, rangl_affector_bounds,
                            atom_idx);
  const int nrangl = afl_extractions.size();
  for (int i = 0; i < nrangl; i++) {
    const int rangl_term_index = afl_extractions[i];
    result.push_back(rar.rangl_i_atoms[rangl_term_index]);
    result.push_back(rar.rangl_j_atoms[rangl_term_index]);
    result.push_back(rar.rangl_k_atoms[rangl_term_index]);
  }
  extractBoundedListEntries(&afl_extractions, rdihe_affector_list, rdihe_affector_bounds,
                            atom_idx);
  const int nrdihe = afl_extractions.size();
  for (int i = 0; i < nrdihe; i++) {
    const int rdihe_term_index = afl_extractions[i];
    result.push_back(rar.rdihe_i_atoms[rdihe_term_index]);
    result.push_back(rar.rdihe_j_atoms[rdihe_term_index]);
    result.push_back(rar.rdihe_k_atoms[rdihe_term_index]);
    result.push_back(rar.rdihe_l_atoms[rdihe_term_index]);
  }

  // Virtual sites transfer the non-bonded forces they have accumuated (including contributions
  // from attenuated 1:4 non-bonded interactions, possibly standard valence terms if the force
  // field is very abstract) to their frame atoms.  The virtual site, all of its frame atoms, and
  // any atoms making interactions with the virtual site must therefore be included.
  extractBoundedListEntries(&afl_extractions, vste_affector_list, vste_affector_bounds, atom_idx);
  const int nvsite = afl_extractions.size();
  for (int i = 0; i < nvsite; i++) {
    const int vste_index = afl_extractions[i];

    // If the virtual site is atom_idx, then the primary and inferred 1:4 interactions have already
    // been counted in the work above.  This part of the step can be skipped, and since the atom
    // itself has already been added to the list of relevant atoms, it is not necessary to add it
    // again.  If the virtual site is not atom_idx, this step will cover adding the virtual site
    // itself, and all of the particles that participate in valence or restraint forces with it, to
    // the list of relevant atoms.
    const int vs_idx = vsk.vs_atoms[vste_index];
    if (vs_idx != atom_idx) {
      const std::vector<int> tmpv = findForcePartners(vs_idx, updated_caller_stack);
      result.insert(result.end(), tmpv.begin(), tmpv.end());
    }

    // Regardless of whether atoms relevant to additional interactions to the virtual site were
    // included, the frame atoms must be added to the list of relevant atoms so that it is known
    // how forces from the virtual site will be split up, some of which will land upon atom_idx
    // if atom_idx was not the virtual site itself.
    result.push_back(vsk.frame1_idx[vste_index]);
    result.push_back(vsk.frame2_idx[vste_index]);
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[vste_index])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      result.push_back(vsk.frame3_idx[vste_index]);
      break;
    case VirtualSiteKind::FIXED_4:
      result.push_back(vsk.frame3_idx[vste_index]);
      result.push_back(vsk.frame4_idx[vste_index]);
      break;
    }
  }
    
  // Sort the list of required atoms and prune duplicate entries
  reduceUniqueValues(&result);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getUpdateDependencies(const int atom_index) const {
  std::vector<int> result;
  std::vector<int> move_req = findMovementPartners(atom_index);
  const size_t nmove = move_req.size();
  for (size_t i = 0; i < nmove; i++) {

    // Each of the atoms in move_req will be included in the tmpj vector, so no need to
    // explicitly add the contents of move_req to the list.
    const std::vector<int> tmpi = findForcePartners(move_req[i]);
    result.insert(result.end(), tmpi.begin(), tmpi.end());
  }
  reduceUniqueValues(&result);
  return result;
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getUpdateWorkUnit(const int atom_index) const {
  return assigned_update_work_units[atom_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getBondAccumulatorWorkUnit(const int bond_index) const {
  return assigned_bond_acc_work_units[bond_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getAngleAccumulatorWorkUnit(const int angl_index) const {
  return assigned_angl_acc_work_units[angl_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getDihedralAccumulatorWorkUnit(const int dihe_index) const {
  return assigned_dihe_acc_work_units[dihe_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getUreyBradleyAccumulatorWorkUnit(const int ubrd_index) const {
  return assigned_ubrd_acc_work_units[ubrd_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getCharmmImproperAccumulatorWorkUnit(const int cimp_index) const {
  return assigned_cimp_acc_work_units[cimp_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getCmapAccumulatorWorkUnit(const int cmap_index) const {
  return assigned_cmap_acc_work_units[cmap_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getInferred14AccumulatorWorkUnit(const int infr14_index) const {
  return assigned_infr14_acc_work_units[infr14_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getPositionalRestraintAccumulatorWorkUnit(const int rposn_index) const {
  return assigned_rposn_acc_work_units[rposn_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getDistanceRestraintAccumulatorWorkUnit(const int rbond_index) const {
  return assigned_rbond_acc_work_units[rbond_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getAngleRestraintAccumulatorWorkUnit(const int rangl_index) const {
  return assigned_rangl_acc_work_units[rangl_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getDihedralRestraintAccumulatorWorkUnit(const int rdihe_index) const {
  return assigned_rdihe_acc_work_units[rdihe_index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ValenceDelegator::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus* ValenceDelegator::getRestraintApparatusPointer() const {
  return ra_pointer;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::checkPresence(const int atom_index, const int vwu_index) const {
  bool result = false;
  const int offset = atom_index * max_presence_allocation;
  for (int i = 0; i < work_unit_assignment_count[atom_index]; i++) {
    result = (result || work_unit_presence[offset + i] == vwu_index);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markAtomAddition(const int vwu_index, const int atom_index) {
  const int atom_wua = work_unit_assignment_count[atom_index];
  
  // If necessary, resize and repaginate the record of which work units each atom takes part in
  if (atom_wua == max_presence_allocation) {
    work_unit_presence.resize(atom_count * (max_presence_allocation + 1));
    for (int i = atom_count - 1; i >= 0; i--) {
      const int new_offset = (max_presence_allocation + 1) * i;
      const int old_offset = max_presence_allocation * i;
      for (int j = 0; j < max_presence_allocation; j++) {
        work_unit_presence[new_offset + j] = work_unit_presence[old_offset + j];
      }
    }
    max_presence_allocation += 1;
  }

  // Mark the new work unit to which the indexed atom has been assigned
  work_unit_presence[(max_presence_allocation * atom_index) + atom_wua] = vwu_index;
  work_unit_assignment_count[atom_index] = atom_wua + 1;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setUpdateWorkUnit(const int atom_index, const int vwu_index) {
  if (assigned_update_work_units[atom_index] >= 0) {
    return false;
  }
  assigned_update_work_units[atom_index] = vwu_index;

  // Update the first unassigned atom, if appropriate
  while (first_unassigned_atom < atom_count &&
         assigned_update_work_units[first_unassigned_atom] >= 0) {
    first_unassigned_atom += 1;
  }
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setBondAccumulatorWorkUnit(const int bond_index, const int vwu_index) {
  if (assigned_bond_acc_work_units[bond_index] >= 0) {
    return false;
  }
  assigned_bond_acc_work_units[bond_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setAngleAccumulatorWorkUnit(const int angl_index, const int vwu_index) {
  if (assigned_angl_acc_work_units[angl_index] >= 0) {
    return false;
  }
  assigned_angl_acc_work_units[angl_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setDihedralAccumulatorWorkUnit(const int dihe_index, const int vwu_index) {
  if (assigned_dihe_acc_work_units[dihe_index] >= 0) {
    return false;
  }
  assigned_dihe_acc_work_units[dihe_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setUreyBradleyAccumulatorWorkUnit(const int ubrd_index,
                                                         const int vwu_index) {
  if (assigned_ubrd_acc_work_units[ubrd_index] >= 0) {
    return false;
  }
  assigned_ubrd_acc_work_units[ubrd_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setCharmmImproperAccumulatorWorkUnit(const int cimp_index,
                                                            const int vwu_index) {
  if (assigned_cimp_acc_work_units[cimp_index] >= 0) {
    return false;
  }
  assigned_cimp_acc_work_units[cimp_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setCmapAccumulatorWorkUnit(const int cmap_index, const int vwu_index) {
  if (assigned_cmap_acc_work_units[cmap_index] >= 0) {
    return false;
  }
  assigned_cmap_acc_work_units[cmap_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setInferred14AccumulatorWorkUnit(const int infr14_index,
                                                        const int vwu_index) {
  if (assigned_infr14_acc_work_units[infr14_index] >= 0) {
    return false;
  }
  assigned_infr14_acc_work_units[infr14_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setPositionalRestraintAccumulatorWorkUnit(const int rposn_index,
                                                                 const int vwu_index) {
  if (assigned_rposn_acc_work_units[rposn_index] >= 0) {
    return false;
  }
  assigned_rposn_acc_work_units[rposn_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setDistanceRestraintAccumulatorWorkUnit(const int rbond_index,
                                                               const int vwu_index) {
  if (assigned_rbond_acc_work_units[rbond_index] >= 0) {
    return false;
  }
  assigned_rbond_acc_work_units[rbond_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setAngleRestraintAccumulatorWorkUnit(const int rangl_index,
                                                            const int vwu_index) {
  if (assigned_rangl_acc_work_units[rangl_index] >= 0) {
    return false;
  }
  assigned_rangl_acc_work_units[rangl_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setDihedralRestraintAccumulatorWorkUnit(const int rdihe_index,
                                                               const int vwu_index) {
  if (assigned_rdihe_acc_work_units[rdihe_index] >= 0) {
    return false;
  }
  assigned_rdihe_acc_work_units[rdihe_index] = vwu_index;
  return true;
}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::allocate() {
 
  // Allocate memory as needed
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag_pointer->getDoublePrecisionConstraintKit();
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();
  bond_affector_list.resize(2 * vk.nbond);
  bond_affector_bounds.resize(atom_count + 1, 0);
  angl_affector_list.resize(3 * vk.nangl);
  angl_affector_bounds.resize(atom_count + 1, 0);
  dihe_affector_list.resize(4 * vk.ndihe);
  dihe_affector_bounds.resize(atom_count + 1, 0);
  ubrd_affector_list.resize(2 * vk.nubrd);
  ubrd_affector_bounds.resize(atom_count + 1, 0);
  cimp_affector_list.resize(4 * vk.ncimp);
  cimp_affector_bounds.resize(atom_count + 1, 0);
  cmap_affector_list.resize(5 * vk.ncmap);
  cmap_affector_bounds.resize(atom_count + 1, 0);
  infr_affector_list.resize(2 * vk.ninfr14);
  infr_affector_bounds.resize(atom_count + 1, 0);
  rposn_affector_list.resize(rar.nposn);
  rposn_affector_bounds.resize(atom_count + 1, 0);
  rbond_affector_list.resize(2 * rar.nbond);
  rbond_affector_bounds.resize(atom_count + 1, 0);
  rangl_affector_list.resize(3 * rar.nangl);
  rangl_affector_bounds.resize(atom_count + 1, 0);
  rdihe_affector_list.resize(4 * rar.ndihe);
  rdihe_affector_bounds.resize(atom_count + 1, 0);
  vste_affector_list.resize(5 * vsk.nsite);
  vste_affector_bounds.resize(atom_count + 1, 0);
  sett_affector_list.resize(3 * cnk.nsettle);
  sett_affector_bounds.resize(atom_count + 1, 0);
  cnst_affector_list.resize(cnk.group_bounds[cnk.ngroup], -1);
  cnst_affector_bounds.resize(atom_count + 1, 0);  
  work_unit_assignment_count.resize(atom_count, 0);
  work_unit_presence.resize(max_presence_allocation * atom_count);
  assigned_update_work_units.resize(atom_count, -1);
  assigned_bond_acc_work_units.resize(vk.nbond, -1);
  assigned_angl_acc_work_units.resize(vk.nangl, -1);
  assigned_dihe_acc_work_units.resize(vk.ndihe, -1);
  assigned_ubrd_acc_work_units.resize(vk.nubrd, -1);
  assigned_cimp_acc_work_units.resize(vk.ncimp, -1);
  assigned_cmap_acc_work_units.resize(vk.ncmap, -1);
  assigned_infr14_acc_work_units.resize(vk.ninfr14, -1);
  assigned_rposn_acc_work_units.resize(rar.nposn, -1);
  assigned_rbond_acc_work_units.resize(rar.nbond, -1);
  assigned_rangl_acc_work_units.resize(rar.nangl, -1);
  assigned_rdihe_acc_work_units.resize(rar.ndihe, -1);
}
  
//-------------------------------------------------------------------------------------------------
void ValenceDelegator::fillAffectorArrays(const ValenceKit<double> &vk,
                                          const VirtualSiteKit<double> &vsk,
                                          const ConstraintKit<double> &cnk,
                                          const RestraintKit<double, double2, double4> &rar) {

  // Pass through the topology, filling out the valence term affector arrays and the virtual site
  // frame atom arrays.
  markAffectorAtoms(&bond_affector_bounds, &bond_affector_list, vk.nbond, vk.bond_i_atoms,
                    vk.bond_j_atoms);
  markAffectorAtoms(&angl_affector_bounds, &angl_affector_list, vk.nangl, vk.angl_i_atoms,
                    vk.angl_j_atoms, vk.angl_k_atoms);
  markAffectorAtoms(&dihe_affector_bounds, &dihe_affector_list, vk.ndihe, vk.dihe_i_atoms,
                    vk.dihe_j_atoms, vk.dihe_k_atoms, vk.dihe_l_atoms);
  markAffectorAtoms(&ubrd_affector_bounds, &ubrd_affector_list, vk.nubrd, vk.ubrd_i_atoms,
                    vk.ubrd_k_atoms);
  markAffectorAtoms(&cimp_affector_bounds, &cimp_affector_list, vk.ncimp, vk.cimp_i_atoms,
                    vk.cimp_j_atoms, vk.cimp_k_atoms, vk.cimp_l_atoms);
  markAffectorAtoms(&cmap_affector_bounds, &cmap_affector_list, vk.ncmap, vk.cmap_i_atoms,
                    vk.cmap_j_atoms, vk.cmap_k_atoms, vk.cmap_l_atoms, vk.cmap_m_atoms);
  markAffectorAtoms(&infr_affector_bounds, &infr_affector_list, vk.ninfr14, vk.infr14_i_atoms,
                    vk.infr14_l_atoms);
  markAffectorAtoms(&sett_affector_bounds, &sett_affector_list, cnk.nsettle, cnk.settle_ox_atoms,
                    cnk.settle_h1_atoms, cnk.settle_h2_atoms);
  markAffectorAtoms(&rposn_affector_bounds, &rposn_affector_list, rar.nposn, rar.rposn_atoms);
  markAffectorAtoms(&rbond_affector_bounds, &rbond_affector_list, rar.nbond, rar.rbond_i_atoms,
                    rar.rbond_j_atoms);
  markAffectorAtoms(&rangl_affector_bounds, &rangl_affector_list, rar.nangl, rar.rangl_i_atoms,
                    rar.rangl_j_atoms, rar.rangl_k_atoms);
  markAffectorAtoms(&rdihe_affector_bounds, &rdihe_affector_list, rar.ndihe, rar.rdihe_i_atoms,
                    rar.rdihe_j_atoms, rar.rdihe_k_atoms, rar.rdihe_l_atoms);

  // Virtual sites can contain -1 in the atom arrays and therefore cannot go through the
  // general procedure.  Best to handle them with a filter over the frame type, anyway.
  for (int pos = 0; pos < vsk.nsite; pos++) {

    // There will be a minimum of two frame atoms, but the virtual site itself is part of the
    // "term" that puts forces on particles.  The non-bonded forces thus far determined to act
    // on the virtual site trasmit onto as many as four frame atoms.
    vste_affector_bounds[vsk.vs_atoms[pos]] += 1;
    vste_affector_bounds[vsk.frame1_idx[pos]] += 1;
    vste_affector_bounds[vsk.frame2_idx[pos]] += 1;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[vsk.vs_param_idx[pos]])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
      break;
    case VirtualSiteKind::FIXED_4:
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
      vste_affector_bounds[vsk.frame4_idx[pos]] += 1;
      break;
    }
  }

  // Constraint groups are not of any pre-defined size and therefore cannot go through the
  // general procedure.
  for (int pos = 0; pos < cnk.ngroup; pos++) {
    for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
      cnst_affector_bounds[cnk.group_list[j]] += 1;
    }
  }
  prefixSumInPlace(&vste_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace(&cnst_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  
  // Finish up the virtual site mapping
  for (int pos = 0; pos < vsk.nsite; pos++) {
    const int virtual_atom = vsk.vs_atoms[pos];
    const int parent_atom = vsk.frame1_idx[pos];
    const int frame2_atom = vsk.frame2_idx[pos];
    const int list_vs_idx = vste_affector_bounds[virtual_atom];
    const int list_p_idx  = vste_affector_bounds[parent_atom];
    const int list_f2_idx = vste_affector_bounds[frame2_atom];
    vste_affector_list[list_vs_idx] = pos;
    vste_affector_list[list_p_idx ] = pos;
    vste_affector_list[list_f2_idx] = pos;
    vste_affector_bounds[virtual_atom] += 1;
    vste_affector_bounds[parent_atom]  += 1;
    vste_affector_bounds[frame2_atom]  += 1;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[vsk.vs_param_idx[pos]])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom = vsk.frame3_idx[pos];
        const int list_f3_idx = vste_affector_bounds[frame3_atom];
        vste_affector_list[list_f3_idx] = pos;
        vste_affector_bounds[frame3_atom] += 1;
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom = vsk.frame3_idx[pos];
        const int frame4_atom = vsk.frame4_idx[pos];
        const int list_f3_idx = vste_affector_bounds[frame3_atom];
        const int list_f4_idx = vste_affector_bounds[frame4_atom];
        vste_affector_list[list_f3_idx] = pos;
        vste_affector_list[list_f4_idx] = pos;
        vste_affector_bounds[frame3_atom] += 1;
        vste_affector_bounds[frame4_atom] += 1;
      }
      break;
    }
  }

  // Finish up the constraint group mapping
  for (int pos = 0; pos < cnk.ngroup; pos++) {
    for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
      const int c_atom = cnk.group_list[j];
      const int list_c_idx = cnst_affector_bounds[c_atom];
      cnst_affector_list[list_c_idx] = pos;
      cnst_affector_bounds[c_atom] += 1;
    }
  }

  // Rewind the final prefix sums after using them to populate the virtual site and constraint
  // affector lists.
  const int natom = atom_count;
  for (int i = natom; i > 0; i--) {
    vste_affector_bounds[i]  = vste_affector_bounds[i - 1];
    cnst_affector_bounds[i]  = cnst_affector_bounds[i - 1];
  }
  vste_affector_bounds[0] = 0;
  cnst_affector_bounds[0] = 0;
}
  
//-------------------------------------------------------------------------------------------------
ValenceWorkUnit::ValenceWorkUnit(ValenceDelegator *vdel_in, std::vector<int> *tvwu_coverage,
                                 const int list_index_in, const int seed_atom_in,
                                 const int max_atoms_in) :
    imported_atom_count{0}, moved_atom_count{0}, updated_atom_count{0}, bond_term_count{0},
    angl_term_count{0}, dihe_term_count{0}, ubrd_term_count{0}, cbnd_term_count{0},
    cimp_term_count{0}, cdhe_term_count{0}, cmap_term_count{0}, infr14_term_count{0},
    rposn_term_count{0}, rbond_term_count{0}, rangl_term_count{0}, rdihe_term_count{0},
    cnst_group_count{0}, sett_group_count{0}, vste_count{0}, list_index{list_index_in},
    min_atom_index{-1}, max_atom_index{-1}, atom_limit{max_atoms_in}, atom_import_list{},
    atom_update_mask{}, bond_term_list{}, angl_term_list{}, dihe_term_list{}, ubrd_term_list{},
    cimp_term_list{}, cmap_term_list{}, infr14_term_list{}, bond_i_atoms{}, bond_j_atoms{},
    angl_i_atoms{}, angl_j_atoms{}, angl_k_atoms{}, dihe_i_atoms{}, dihe_j_atoms{}, dihe_k_atoms{},
    dihe_l_atoms{}, ubrd_i_atoms{}, ubrd_k_atoms{}, cimp_i_atoms{}, cimp_j_atoms{}, cimp_k_atoms{},
    cimp_l_atoms{}, cmap_i_atoms{}, cmap_j_atoms{}, cmap_k_atoms{}, cmap_l_atoms{}, cmap_m_atoms{},
    infr14_i_atoms{}, infr14_l_atoms{}, cbnd_term_list{}, cbnd_is_ubrd{}, cbnd_i_atoms{},
    cbnd_jk_atoms{}, cdhe_term_list{}, cdhe_is_cimp{}, cdhe_i_atoms{}, cdhe_j_atoms{},
    cdhe_k_atoms{}, cdhe_l_atoms{}, rposn_term_list{}, rbond_term_list{}, rangl_term_list{},
    rdihe_term_list{}, rposn_atoms{}, rbond_i_atoms{}, rbond_j_atoms{}, rangl_i_atoms{},
    rangl_j_atoms{}, rangl_k_atoms{}, rdihe_i_atoms{}, rdihe_j_atoms{}, rdihe_k_atoms{},
    rdihe_l_atoms{}, acc_bond_energy{}, acc_angl_energy{}, acc_dihe_energy{}, acc_ubrd_energy{},
    acc_cimp_energy{}, acc_cmap_energy{}, acc_infr14_energy{}, acc_rposn_energy{},
    acc_rbond_energy{}, acc_rangl_energy{}, acc_rdihe_energy{}, acc_cbnd_energy{},
    acc_cdhe_energy{}, cnst_group_list{}, sett_group_list{}, cnst_group_atoms{},
    cnst_group_bounds{}, sett_ox_atoms{}, sett_h1_atoms{}, sett_h2_atoms{}, virtual_site_list{},
    vsite_atoms{}, vsite_frame1_atoms{}, vsite_frame2_atoms{}, vsite_frame3_atoms{},
    vsite_frame4_atoms{}, cbnd_instructions{}, angl_instructions{}, cdhe_instructions{},
    cmap_instructions{}, infr14_instructions{}, rposn_instructions{}, rbond_instructions{},
    rangl_instructions{}, rdihe_instructions{}, vste_instructions{}, sett_instructions{},
    cnst_instructions{}, vdel_pointer{vdel_in}, ag_pointer{vdel_in->getTopologyPointer()},
    ra_pointer{vdel_in->getRestraintApparatusPointer()}
{
  // Check the atom bounds
  if (atom_limit < maximum_valence_work_unit_atoms / 16 ||
      atom_limit > maximum_valence_work_unit_atoms) {
    const std::string err_msg = (atom_limit > maximum_valence_work_unit_atoms) ?
                                "could lead to work units that do not fit in HPC resources." :
                                "is too small to make effective work units.";
    rtErr("The maximum allowed number of atoms should be between " +
          std::to_string(minimum_valence_work_unit_atoms) + " and " +
          std::to_string(maximum_valence_work_unit_atoms) + ".  A value of " +
          std::to_string(atom_limit) + err_msg, "ValenceWorkUnit");
  }

  // Shortcut to the data array for this valence work unit's coverage in the topology as a whole
  int* tvwu_cov_ptr = tvwu_coverage->data();
  
  // Reserve space for atoms
  atom_import_list.reserve(atom_limit);

  // Unpack the original topology
  ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  
  // Starting with the seed atom, branch out and add new atoms until approaching the maximum
  // allowed number of atoms.  If the atom is part of some molecule that can all be accommodated
  // by the one work unit, include all atoms.
  std::vector<int> candidate_additions(1, seed_atom_in);
  std::vector<int> growth_points;
  growth_points.reserve(32);
  candidate_additions.reserve(32);
  int fua_atom = vdel_pointer->getFirstUnassignedAtom();

  // If no atoms have yet been assigned and all atoms of this system will fit within the work
  // unit's bounds, simply log them all and return.
  if (fua_atom == 0&& cdk.natom <= atom_limit) {
    for (int i = 0; i < cdk.natom; i++) {
      addNewAtomImport(i);
      addNewAtomUpdate(i);
    }
    return;
  }
  
  if (fua_atom != seed_atom_in && fua_atom < cdk.natom) {
    candidate_additions.push_back(fua_atom);
  }
  
  // Proceed to add atoms layer by layer in the current molecule, or if that has been totally
  // covered, jump to the next molecule and continue the process.
  int ncandidate = candidate_additions.size();
  bool capacity_reached = false;
  while (ncandidate > 0 && capacity_reached == false) {

    // A single successful addition is needed to prove that the capacity has not yet been reached.
    capacity_reached = true;
    
    // Transfer candidate atoms to become growth points.  In construction, valence work units
    // do not overlap.  Atoms included in some previous work unit do not become candidates.
    growth_points.resize(0);
    for (int i = 0; i < ncandidate; i++) {

      // This is an atom that the work unit will be responsible for updating.  Check the number
      // of its dependencies and add them to the import list, if there is room.  Mark the atom
      // as one for the work unit to update.
      const std::vector<int> up_deps = vdel_pointer->getUpdateDependencies(candidate_additions[i]);
      
      const int ndeps = up_deps.size();
      std::vector<bool> incl_up_deps(ndeps, false);
      int n_new_atoms = 0;
      for (int j = 0; j < ndeps; j++) {
        if (vdel_pointer->checkPresence(up_deps[j], list_index) == false &&
            tvwu_cov_ptr[up_deps[j]] == 0) {
          n_new_atoms++;
          incl_up_deps[j] = true;
        }
      }
      if (imported_atom_count + n_new_atoms <= atom_limit) {

        // The candidate atom will be part of its own dependencies list.
        for (int j = 0; j < ndeps; j++) {
          if (incl_up_deps[j]) {
            addNewAtomImport(up_deps[j]);
            tvwu_cov_ptr[up_deps[j]] = 1;
          }
        }
        growth_points.push_back(candidate_additions[i]);
        addNewAtomUpdate(candidate_additions[i]);
        capacity_reached = false;
      }
    }
    candidate_additions.resize(0);

    // Loop over the growth points and determine new candidate atoms.  During construction,
    // valence work units do not overlap.  Atoms included in some previous work unit do not
    // become new candidates.
    const int ngrow = growth_points.size();
    for (int i = 0; i < ngrow; i++) {
      const int grow_atom = growth_points[i];
      for (int j = nbk.nb12_bounds[grow_atom]; j < nbk.nb12_bounds[grow_atom + 1]; j++) {
        if (vdel_pointer->getUpdateWorkUnit(nbk.nb12x[j]) == -1) {
          candidate_additions.push_back(nbk.nb12x[j]);
        }
      }
    }
    reduceUniqueValues(&candidate_additions);
    ncandidate = candidate_additions.size();

    // If no candidate atoms have yet been found, try jumping to the first unassigned atom.
    // In all likelihood, this will be on another molecule within the same topology.  That will
    // be the seed for the next round of additions.
    if (ncandidate == 0) {
      fua_atom = vdel_pointer->getFirstUnassignedAtom();
      if (fua_atom < cdk.natom) {
        candidate_additions.push_back(fua_atom);
        ncandidate = 1;
      }
    }
  }

  // Erase the coverage for this valence work unit so that the coverage array is clean for the
  // next valence work unit to use without reallocating.
  for (int i = 0; i < imported_atom_count; i++) {
    tvwu_cov_ptr[atom_import_list[i]] = 0;
  }
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getImportedAtomCount() const {
  return imported_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMovedAtomCount() const {
  return moved_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getUpdatedAtomCount() const {
  return updated_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getListIndex() const {
  return list_index;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMinAtomIndex() const {
  return min_atom_index;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMaxAtomIndex() const {
  return max_atom_index;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMaxAtoms() const {
  return atom_limit;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceWorkUnit::getAtomImportList(const int atom_offset) const {
  if (atom_offset == 0) {
    return atom_import_list;
  }
  else {
    std::vector<int> result(imported_atom_count);
    for (int i = 0; i < imported_atom_count; i++) {
      result[i] = atom_import_list[i] + atom_offset;
    }
    return result;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getImportedAtomIndex(const int slot, const int atom_offset) const {
  return atom_import_list[slot] + atom_offset;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2> ValenceWorkUnit::getAtomManipulationMasks() const {
  const int nbits = static_cast<int>(sizeof(uint)) * 8;
  const int n_segment = roundUp((imported_atom_count + nbits - 1) / nbits, warp_size_int);
  std::vector<uint2> result(n_segment, {0U, 0U});
  const int nfull = imported_atom_count / nbits;
  for (int i = 0; i < nfull; i++) {
    result[i].x = 0xffffffff;
  }
  for (int i = nfull * nbits; i < imported_atom_count; i++) {
    const int seg_idx = (i / nbits);
    const int bit_idx = i - (seg_idx * nbits);
    result[seg_idx].x |= (0x1 << bit_idx);
  }
  const int mask_length = (imported_atom_count + nbits - 1) / nbits;
  if (atom_update_mask.size() != static_cast<size_t>(mask_length)) {
    rtErr("The atom update mask is not of the anticipated length (" +
          std::to_string(atom_update_mask.size()) + ", " + std::to_string(mask_length) + ").",
          "ValenceWorkUnit", "getAtomManipulationMasks");
  }
  for (int i = 0; i < mask_length; i++) {
    result[i].y = atom_update_mask[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getPaddedConstraintInstructionCount() const {
  int result = 0;
  for (int i = 0; i < cnst_group_count; i++) {

    // Subtract 1 to get the number of constrained bonds--the central atom adds 1 to the total
    // atom count in the group.
    const int total_bonds = cnst_group_bounds[i + 1] - cnst_group_bounds[i] - 1;
    if (total_bonds >= 16) {
      rtErr("The largest constraint group size was determined to involve " +
            std::to_string(total_bonds) + " bonds, which exceeds the maximum allowed size of " +
            std::to_string(max_constraint_group_size) + ".", "ValenceWorkUnit",
            "getPaddedConstraintInstructionCount");
    }
    const int rsiz_limit = roundUp(result + 1, warp_size_int);
    const int rbnd_limit = roundUp(result + total_bonds, warp_size_int);
    if (rsiz_limit != rbnd_limit) {
      result = rsiz_limit;
    }
    result += total_bonds;
  }
  return roundUp(result, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceWorkUnit::getTaskCounts() const {
  std::vector<int> result(static_cast<int>(VwuTask::ALL_TASKS));
  result[static_cast<int>(VwuTask::BOND)]   = bond_term_count;
  result[static_cast<int>(VwuTask::ANGL)]   = angl_term_count;
  result[static_cast<int>(VwuTask::DIHE)]   = dihe_term_count;
  result[static_cast<int>(VwuTask::UBRD)]   = ubrd_term_count;
  result[static_cast<int>(VwuTask::CBND)]   = cbnd_term_count;
  result[static_cast<int>(VwuTask::CIMP)]   = cimp_term_count;
  result[static_cast<int>(VwuTask::CDHE)]   = cdhe_term_count;
  result[static_cast<int>(VwuTask::CMAP)]   = cmap_term_count;
  result[static_cast<int>(VwuTask::INFR14)] = infr14_term_count;
  result[static_cast<int>(VwuTask::RPOSN)]  = rposn_term_count;
  result[static_cast<int>(VwuTask::RBOND)]  = rbond_term_count;
  result[static_cast<int>(VwuTask::RANGL)]  = rangl_term_count;
  result[static_cast<int>(VwuTask::RDIHE)]  = rdihe_term_count;

  // The constraint group count is a standout--the number of tasks is depends on padding and warp
  // sizes.  Each constraint task is one of the members of a group.
  result[static_cast<int>(VwuTask::CGROUP)] = getPaddedConstraintInstructionCount();
  result[static_cast<int>(VwuTask::SETTLE)] = sett_group_count;
  result[static_cast<int>(VwuTask::VSITE)]  = vste_count;
  return result;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeCompositeBondInstructions(const std::vector<int> &bond_param_map,
                                                     const std::vector<int> &ubrd_param_map) {
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  cbnd_instructions.resize(bond_term_count + ubrd_term_count);
  const bool bond_mapped = (bond_param_map.size() > 0LLU);
  const bool ubrd_mapped = (ubrd_param_map.size() > 0LLU);
  if ((vk.nbond > 0 && bond_mapped && vk.nubrd > 0 && ubrd_mapped == false) ||
      (vk.nbond > 0 && bond_mapped == false && vk.nubrd > 0 && ubrd_mapped)) {
    rtErr("Parameter correspondence tables must be supplied for both harmonic bond and "
          "Urey-bradley terms, or neither.", "ValenceWorkUnit", "storeCompositeBondInstructions");
  }
  const bool use_map = (bond_mapped || ubrd_mapped);
  if (use_map && static_cast<int>(bond_param_map.size()) != vk.nbond_param) {
    rtErr("A parameter correspondence table with " + std::to_string(bond_param_map.size()) +
          " entries cannot serve a topology with " + std::to_string(vk.nbond_param) +
          " bond parameter sets.", "ValenceWorkUnit", "storeCompositeBondInstructions");
  }
  if (use_map && static_cast<int>(ubrd_param_map.size()) != vk.nubrd_param) {
    rtErr("A parameter correspondence table with " + std::to_string(ubrd_param_map.size()) +
          " entries cannot serve a topology with " + std::to_string(vk.nbond_param) +
          " Urey-Bradley parameter sets.", "ValenceWorkUnit", "storeCompositeBondInstructions");
  }
  for (int pos = 0; pos < cbnd_term_count; pos++) {
    cbnd_instructions[pos].x = ((cbnd_jk_atoms[pos] << 10) | cbnd_i_atoms[pos]);
    if (cbnd_is_ubrd[pos]) {
      const int ut_idx = vk.ubrd_param_idx[cbnd_term_list[pos]];
      cbnd_instructions[pos].x |= (0x1 << 20);
      cbnd_instructions[pos].y = (use_map) ? ubrd_param_map[ut_idx] : ut_idx;
    }
    else {
      const int bt_idx = vk.bond_param_idx[cbnd_term_list[pos]];
      cbnd_instructions[pos].y = (use_map) ? bond_param_map[bt_idx] : bt_idx;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeAngleInstructions(const std::vector<int> &parameter_map) {
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  angl_instructions.resize(angl_term_count);
  const bool use_map = (parameter_map.size() > 0LLU);
  if (use_map && static_cast<int>(parameter_map.size()) != vk.nangl_param) {
    rtErr("A parameter correspondence table with " + std::to_string(parameter_map.size()) +
          "entries cannot serve a topology with " + std::to_string(vk.nangl_param) +
          " bond angle parameter sets.", "ValenceWorkUnit", "storeAngleInstructions");
  }
  for (int pos = 0; pos < angl_term_count; pos++) {
    const int at_idx = vk.angl_param_idx[angl_term_list[pos]];
    angl_instructions[pos].x = ((angl_k_atoms[pos] << 20) | (angl_j_atoms[pos] << 10) |
                                angl_i_atoms[pos]);
    angl_instructions[pos].y = (use_map) ? parameter_map[at_idx] : at_idx;
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeCompositeDihedralInstructions(const std::vector<int> &dihe_param_map,
                                                         const std::vector<int> &dihe14_param_map,
                                                         const std::vector<int> &cimp_param_map) {
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  cdhe_instructions.resize(cdhe_term_count);
  const bool dihe_mapped = (dihe_param_map.size() > 0LLU);
  const bool cimp_mapped = (cimp_param_map.size() > 0LLU);
  if ((vk.ndihe > 0 && dihe_mapped && vk.ncimp > 0 && cimp_mapped == false) ||
      (vk.ndihe > 0 && dihe_mapped == false && vk.ncimp > 0 && cimp_mapped)) {
    rtErr("Parameter correspondence tables must be supplied for both dihedral and CHARMM improper "
          "terms, or neither.", "ValenceWorkUnit", "storeCompositeDihedralInstructions");
  }
  const bool use_map = (dihe_mapped || cimp_mapped); 
  if (use_map && static_cast<int>(dihe_param_map.size()) != vk.ndihe_param) {
    rtErr("A parameter correspondence table with " + std::to_string(dihe_param_map.size()) +
          "entries cannot serve a topology with " + std::to_string(vk.ndihe_param) +
          " torsion parameter sets.", "ValenceWorkUnit", "storeCompositeDihedralInstructions");
  }
  if (use_map && static_cast<int>(cimp_param_map.size()) != vk.ncimp_param) {
    rtErr("A parameter correspondence table with " + std::to_string(cimp_param_map.size()) +
          "entries cannot serve a topology with " + std::to_string(vk.ncimp_param) +
          " torsion parameter sets.", "ValenceWorkUnit", "storeCompositeDihedralInstructions");
  }
  for (int pos = 0; pos < cdhe_term_count; pos++) {

    // Write out the atom indices
    cdhe_instructions[pos].x = ((cdhe_k_atoms[pos] << 20) | (cdhe_j_atoms[pos] << 10) |
                                cdhe_i_atoms[pos]);
    cdhe_instructions[pos].y = cdhe_l_atoms[pos];

    // Branch based on whether this is a CHARMM improper
    if (cdhe_is_cimp[pos]) {
      cdhe_instructions[pos].x |= (0x1 << 30);

      // Record the parameter index of the CHARMM improper.  Do not permit more than 65536 unique
      // CHARMM improper torsion parameter sets.
      const int ht_idx = vk.cimp_param_idx[cdhe_term_list[pos].x];
      const int ht_mapped_idx = (use_map) ? cimp_param_map[ht_idx] : ht_idx;
      if (ht_mapped_idx >= 65536) {
        const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
        rtErr("In valence work unit " + std::to_string(list_index) + ", the CHARMM improper "
              "dihedral involving atoms " +
              char4ToString(cdk.atom_names[atom_import_list[cdhe_i_atoms[pos]]]) + ", " +
              char4ToString(cdk.atom_names[atom_import_list[cdhe_j_atoms[pos]]]) + ", " +
              char4ToString(cdk.atom_names[atom_import_list[cdhe_k_atoms[pos]]]) + ", and " +
              char4ToString(cdk.atom_names[atom_import_list[cdhe_l_atoms[pos]]]) + " requires a "
              "parameter set with index " + std::to_string(ht_mapped_idx) + ", which exceeds the "
              "allowed number of unique parameter pairs indexed in the first composite dihedral "
              "position.", "ValenceWorkUnit", "storeCompositeDihedralInstructions");        
      }
      cdhe_instructions[pos].y |= (ht_mapped_idx << 16);
      cdhe_instructions[pos].z = 0U;
    }
    else {

      // Construct the first (of possibly two) dihedral instructions

      // Indicate whether this interaction contributes to the cosine proper or improper
      // diehdral energy (this setting will determine the contributions of the second dihedral
      // in the composite pair, if one exists)
      switch (static_cast<TorsionKind>(vk.dihe_modifiers[cdhe_term_list[pos].x].w)) {
      case TorsionKind::PROPER:
      case TorsionKind::PROPER_NO_14:
        break;
      case TorsionKind::IMPROPER:
      case TorsionKind::IMPROPER_NO_14:
        cdhe_instructions[pos].x |= (0x1 << 31);
        break;
      }

      // Record the 1:4 scaling factors for the I:L interaction (up to 32 unique pairs of
      // electrostatic and Lennard-Jones scaling factors are available, with the 0th being no
      // interaction)
      const int dt14_idx = vk.dihe14_param_idx[cdhe_term_list[pos].x];
      const int dt14_mapped_idx = (use_map) ? dihe14_param_map[dt14_idx] : dt14_idx;
      if (dt14_mapped_idx >= 32) {
        const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
        rtErr("In valence work unit " + std::to_string(list_index) + ", the dihedral involving "
              "atoms " + char4ToString(cdk.atom_names[atom_import_list[cdhe_i_atoms[pos]]]) +
              ", " + char4ToString(cdk.atom_names[atom_import_list[cdhe_j_atoms[pos]]]) + ", " +
              char4ToString(cdk.atom_names[atom_import_list[cdhe_k_atoms[pos]]]) + ", and " +
              char4ToString(cdk.atom_names[atom_import_list[cdhe_l_atoms[pos]]]) + " requires a "
              "1:4 scaling factor pair with index " + std::to_string(dt14_mapped_idx) +
              ", which exceeds the allowed number of unique parameter pairs in the first "
              "composite dihedral 1:4 position.", "ValenceWorkUnit",
              "storeCompositeDihedralInstructions");
      }
      cdhe_instructions[pos].y |= (dt14_mapped_idx << 10);
      
      // Indicate whether there is a second dihedral in this composite term
      if (cdhe_term_list[pos].y >= 0) {
        cdhe_instructions[pos].y |= (0x1 << 15);
      }
      
      // Indicate the parameter index of the first dihedral (up to 65535 unique parameter sets).
      // If the index turns out to be 65535 and there are more than 65535 total parameter sets,
      // then 65535 is a valid parameter set but it will not do anything when the code later
      // evaluates this first dihedral.  65535 is the new zero: in this context, it indicates no
      // dihedral term to evaluate.  
      const int dt_idx = vk.dihe_param_idx[cdhe_term_list[pos].x];
      cdhe_instructions[pos].y |= (use_map) ? (std::min(dihe_param_map[dt_idx], 65536) << 16) :
                                   (std::min(dt_idx, 65536) << 16);

      // Record the parameter index for the second dihedral.
      if (cdhe_term_list[pos].y >= 0) {
        const int dt2_idx = vk.dihe_param_idx[cdhe_term_list[pos].y];
        const int dt2_14_idx = vk.dihe14_param_idx[cdhe_term_list[pos].y];
        const int dt2_mapped_idx = (use_map) ? dihe_param_map[dt2_idx] : dt2_idx;
        const int dt2_14_mapped_idx = (use_map) ? dihe14_param_map[dt2_14_idx] : dt2_14_idx;

        // There is no more room to go up--if the index of the 1:4 scaling parameters exceeds 4096
        // or the index of the dihedral term itself exceeds 1048576, that should have been
        // intercepted during the initial construction.  Trap it here for added protection.
        if (dt2_mapped_idx >= 1048576 || dt2_14_mapped_idx >= 4096) {
          const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
          rtErr("In valence work unit " + std::to_string(list_index) + ", the dihedral involving "
                "atoms " + char4ToString(cdk.atom_names[atom_import_list[cdhe_i_atoms[pos]]]) +
                ", " + char4ToString(cdk.atom_names[atom_import_list[cdhe_j_atoms[pos]]]) + ", " +
                char4ToString(cdk.atom_names[atom_import_list[cdhe_k_atoms[pos]]]) + ", and " +
                char4ToString(cdk.atom_names[atom_import_list[cdhe_l_atoms[pos]]]) + " requires a "
                "1:4 scaling factor pair with index " + std::to_string(dt2_14_idx) + " and a "
                "dihedral parameter set with index " + std::to_string(dt2_idx) + ", which exceeds "
                "the allowed number of unique parameters.", "ValenceWorkUnit",
                "storeCompositeDihedralInstructions");
        }
        cdhe_instructions[pos].z = dt2_mapped_idx;
        cdhe_instructions[pos].z |= (dt2_14_mapped_idx << 20);
      }
      else {
        cdhe_instructions[pos].z = 0U;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeCmapInstructions(const std::vector<int> &parameter_map) {
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  cmap_instructions.resize(cmap_term_count);
  const bool use_map = (parameter_map.size() > 0LLU);
  if (use_map && static_cast<int>(parameter_map.size()) != vk.ncmap_surf) {
    rtErr("A parameter correspondence table with " + std::to_string(parameter_map.size()) +
          "entries cannot serve a topology with " + std::to_string(vk.ncmap_surf) +
          " CMAP surfaces.", "ValenceWorkUnit", "storeCmapInstructions");
  }
  for (int pos = 0; pos < cmap_term_count; pos++) {
    const int mt_idx = vk.cmap_surf_idx[cmap_term_list[pos]];
    cmap_instructions[pos].x = ((cmap_k_atoms[pos] << 20) | (cmap_j_atoms[pos] << 10) |
                                cmap_i_atoms[pos]);
    cmap_instructions[pos].y = ((cmap_m_atoms[pos] << 10) | cmap_l_atoms[pos]);
    cmap_instructions[pos].y |= (use_map) ? (parameter_map[mt_idx] << 20) : (mt_idx << 20);
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeInferred14Instructions(const std::vector<int> &parameter_map) {
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  infr14_instructions.resize(infr14_term_count);
  const bool use_map = (parameter_map.size() > 0LLU);
  if (use_map && static_cast<int>(parameter_map.size()) != vk.nattn14_param) {
    rtErr("A parameter correspondence table with " + std::to_string(parameter_map.size()) +
          "entries cannot serve a topology with " + std::to_string(vk.nattn14_param) +
          " pairs of 1:4 scaling factors.", "ValenceWorkUnit", "storeInferred14Instructions");
  }
  for (int pos = 0; pos < infr14_term_count; pos++) {
    const int ft_idx = vk.infr14_param_idx[infr14_term_list[pos]];
    infr14_instructions[pos] = ((infr14_l_atoms[pos] << 10) | infr14_i_atoms[pos]);
    infr14_instructions[pos] |= (use_map) ? (parameter_map[ft_idx] << 20) : (ft_idx << 20);
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storePositionalRestraintInstructions(const std::vector<int> &kr_param_map,
                                                           const std::vector<int> &xyz_param_map) {
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();
  rposn_instructions.resize(rposn_term_count);
  const bool use_map = (kr_param_map.size() > 0LLU);
  if (kr_param_map.size() != xyz_param_map.size()) {
    rtErr("Parameter correspondence tables must be supplied for both the k(1,2) / r(1,2,3,4) "
          "parameter set and the x / y / z target positions, or neither.", "ValenceWorkUnit",
          "storePositionalRestraintInstructions");
  }
  if (use_map && static_cast<int>(kr_param_map.size()) != rar.nposn) {
    rtErr("A parameter correspondence table with " + std::to_string(kr_param_map.size()) +
          "entries cannot serve a RestraintApparatus with " + std::to_string(rar.nposn) +
          "terms.", "ValenceWorkUnit", "storePositionalRestraintInstructions");
  }
  for (int pos = 0; pos < rposn_term_count; pos++) {

    // The RestraintApparatus object does not have a concept of restraints as parameters with
    // parameter indices into a reduced table of unique values.  There, the k(2,3) and r(1,2,3,4)
    // parameters of the ith restraint are the ith entries of the k and r data arrays.  In the
    // AtomGraphSynthesis, however, there is a concept of restraint parameter uniqueness.  The
    // instructions are tailored either way, depending on whether the parameter correspondence
    // from an AtomGraphSynthesis has been supplied.  It's analogous to what happens with each of
    // the valence terms above, but for those terms the AtomGraph itself stores parameter tables
    // for each term, so there is an extra layer of parameter referencing to unroll in order to
    // get to the AtomGraphSynthesis parameter indexing.
    const int kr_idx  = (use_map) ? kr_param_map[rposn_term_list[pos]] : rposn_term_list[pos];
    const int xyz_idx = (use_map) ? xyz_param_map[rposn_term_list[pos]] : rposn_term_list[pos];
    rposn_instructions[pos].x = rposn_atoms[pos];
    rposn_instructions[pos].x |= (kr_idx << 10);
    if (rar.time_dependence) {
      rposn_instructions[pos].x |= (0x1 << 31);
    }
    rposn_instructions[pos].y = xyz_idx;
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeDistanceRestraintInstructions(const std::vector<int> &kr_param_map) {
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();
  rbond_instructions.resize(rbond_term_count);
  const bool use_map = (kr_param_map.size() > 0LLU);
  if (use_map && static_cast<int>(kr_param_map.size()) != rar.nbond) {
    rtErr("A parameter correspondence table with " + std::to_string(kr_param_map.size()) +
          "entries cannot serve a RestraintApparatus with " + std::to_string(rar.nbond) +
          "terms.", "ValenceWorkUnit", "storeDistanceRestraintInstructions");
  }
  for (int pos = 0; pos < rbond_term_count; pos++) {
    const int kr_idx  = (use_map) ? kr_param_map[rbond_term_list[pos]] : rbond_term_list[pos];
    rbond_instructions[pos].x = ((rbond_j_atoms[pos] << 10) | rbond_i_atoms[pos]);
    rbond_instructions[pos].y = kr_idx;
    if (rar.time_dependence) {
      rbond_instructions[pos].x |= (0x1 << 31);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeAngleRestraintInstructions(const std::vector<int> &kr_param_map) {
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();
  rangl_instructions.resize(rangl_term_count);
  const bool use_map = (kr_param_map.size() > 0LLU);
  if (use_map && static_cast<int>(kr_param_map.size()) != rar.nangl) {
    rtErr("A parameter correspondence table with " + std::to_string(kr_param_map.size()) +
          "entries cannot serve a RestraintApparatus with " + std::to_string(rar.nangl) +
          "terms.", "ValenceWorkUnit", "storeAngleRestraintInstructions");
  }
  for (int pos = 0; pos < rangl_term_count; pos++) {
    const int kr_idx  = (use_map) ? kr_param_map[rangl_term_list[pos]] : rangl_term_list[pos];
    rangl_instructions[pos].x = ((rangl_k_atoms[pos] << 20) | (rangl_j_atoms[pos] << 10) |
                                 rangl_i_atoms[pos]);
    rangl_instructions[pos].y = kr_idx;
    if (rar.time_dependence) {
      rangl_instructions[pos].x |= (0x1 << 31);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeDihedralRestraintInstructions(const std::vector<int> &kr_param_map) {
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();
  rdihe_instructions.resize(rdihe_term_count);
  const bool use_map = (kr_param_map.size() > 0LLU);
  if (use_map && static_cast<int>(kr_param_map.size()) != rar.ndihe) {
    rtErr("A parameter correspondence table with " + std::to_string(kr_param_map.size()) +
          "entries cannot serve a RestraintApparatus with " + std::to_string(rar.ndihe) +
          "terms.", "ValenceWorkUnit", "storeDihedralRestraintInstructions");
  }
  for (int pos = 0; pos < rdihe_term_count; pos++) {
    const int kr_idx  = (use_map) ? kr_param_map[rdihe_term_list[pos]] : rdihe_term_list[pos];
    rdihe_instructions[pos].x = ((rdihe_k_atoms[pos] << 20) | (rdihe_j_atoms[pos] << 10) |
                                 rdihe_i_atoms[pos]);
    rdihe_instructions[pos].y = ((kr_idx << 10) | rdihe_l_atoms[pos]);
    if (rar.time_dependence) {
      rdihe_instructions[pos].x |= (0x1 << 31);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeVirtualSiteInstructions(const std::vector<int> &parameter_map) {
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  vste_instructions.resize(vste_count);
  const bool use_map = (parameter_map.size() > 0LLU);
  if (use_map && static_cast<int>(parameter_map.size()) != vsk.nframe_set) {
    rtErr("A parameter correspondence table with " + std::to_string(parameter_map.size()) +
          "entries cannot serve a topology with " + std::to_string(vsk.nframe_set) +
          " virtual sites.", "ValenceWorkUnit", "storeVirtualSiteInstructions");
  }
  for (int pos = 0; pos < vste_count; pos++) {
    const int vp_idx = vsk.vs_param_idx[virtual_site_list[pos]];
    vste_instructions[pos].x = ((vsite_frame2_atoms[pos] << 20) | (vsite_frame1_atoms[pos] << 10) |
                                vsite_atoms[pos]);
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[vp_idx])) {
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::NONE:

      // The y member must be initialized here, to ensure that information written in previous
      // passes cannot linger and pollute subsequent expressions given to the y member.
      vste_instructions[pos].y = 0;
      break;
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      vste_instructions[pos].y = vsite_frame3_atoms[pos];
      break;
    case VirtualSiteKind::FIXED_4:
      vste_instructions[pos].y = ((vsite_frame4_atoms[pos] << 10) | vsite_frame3_atoms[pos]);
      break;
    }
    vste_instructions[pos].y |= (use_map) ? (parameter_map[vp_idx] << 20) : (vp_idx << 20);
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::storeSettleGroupInstructions(const std::vector<int> &parameter_map) {
  const ConstraintKit<double> cnk = ag_pointer->getDoublePrecisionConstraintKit();
  sett_instructions.resize(sett_group_count);
  const bool use_map = (parameter_map.size() > 0LLU);
  if (use_map && static_cast<int>(parameter_map.size()) != cnk.nsett_param) {
    rtErr("A parameter correspondence table with " + std::to_string(parameter_map.size()) +
          "entries cannot serve a topology with " + std::to_string(cnk.nsett_param) +
          " SETTLE configurations.", "ValenceWorkUnit", "storeVirtualSiteInstructions");
  }
  for (int pos = 0; pos < sett_group_count; pos++) {
    const int st_idx = cnk.settle_param_idx[sett_group_list[pos]];
    sett_instructions[pos].x = ((sett_h2_atoms[pos] << 20) | (sett_h1_atoms[pos] << 10) |
                                sett_ox_atoms[pos]);
    sett_instructions[pos].y = (use_map) ? parameter_map[st_idx] : st_idx;
  }
}

//-------------------------------------------------------------------------------------------------
void
ValenceWorkUnit::storeConstraintGroupInstructions(const std::vector<int> &parameter_map,
                                                  const std::vector<int> &group_param_bounds) {
  const ConstraintKit<double> cnk = ag_pointer->getDoublePrecisionConstraintKit();
  const int result_size = getPaddedConstraintInstructionCount();
  cnst_instructions.resize(result_size);
  const bool use_map = (parameter_map.size() > 0LLU);
  if (use_map && static_cast<int>(parameter_map.size()) != cnk.ncnst_param) {
    rtErr("A parameter correspondence table with " + std::to_string(parameter_map.size()) +
          "entries cannot serve a topology with " + std::to_string(cnk.ncnst_param) +
          " constraint configurations.", "ValenceWorkUnit", "storeConstraintGroupInstructions");
  }
  if (use_map && maxValue(parameter_map) >= static_cast<int>(group_param_bounds.size())) {
    rtErr("If a parameter correspondence table is provided, a bounds array that accommodates each "
          "of those constraint group parameter sets must be provided.  The correspondence table "
          "lists parameters up to " + std::to_string(maxValue(parameter_map)) + ", whereas the "
          "constraint group parameter bounds array only accommodates " +
          std::to_string(group_param_bounds.size() - 1LLU) + " unique parameter sets.",
          "ValenceWorkUnit", "storeConstraintGroupInstructions");
  }
  int ridx = 0;
  for (int pos = 0; pos < cnst_group_count; pos++) {
    const int ct_idx = cnk.group_param_idx[cnst_group_list[pos]];
    const int actual_param_idx = (use_map) ? parameter_map[ct_idx] : ct_idx;
    const int param_llim = (use_map) ? group_param_bounds[actual_param_idx] :
                                       cnk.group_param_bounds[actual_param_idx];
    const int total_bonds = cnst_group_bounds[pos + 1] - cnst_group_bounds[pos] - 1;

    // Check that the constraint group will fit within the warp.  If not, pad the rest of the warp
    // with blank interactions.
    const int ridx_limit = roundUp(ridx + 1, warp_size_int);
    const int rbnd_limit = roundUp(ridx + total_bonds, warp_size_int);
    if (ridx_limit != rbnd_limit) {
      while (ridx < ridx_limit) {
        
        // Record atom indices as zero, and show the total number of bonds and participating
        // threads in this group.  Threads that do not have a valid bond to constrain will detect
        // that by comparing their position in the warp with these two numbers.
        cnst_instructions[ridx].x = ((ridx & warp_bits_mask_int) << 20);
        cnst_instructions[ridx].y = 0U;
        ridx++;
      }
    }
    const int ridx_base = ridx;
    for (int i = cnst_group_bounds[pos] + 1; i < cnst_group_bounds[pos + 1]; i++) {
      cnst_instructions[ridx].x = ((cnst_group_atoms[i] << 10) |
                                   cnst_group_atoms[cnst_group_bounds[pos]]);
      cnst_instructions[ridx].x |= ((total_bonds << 28) |
                                    ((ridx_base & warp_bits_mask_int) << 20));

      // Use the loop control variable as the offset here, but do not subtract 1 as the
      // length parameter array will have a blank for the central atom.
      cnst_instructions[ridx].y = param_llim + i - cnst_group_bounds[pos];
      ridx++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getCompositeBondInstructions() const {
  return cbnd_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getAngleInstructions() const {
  return angl_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint3>& ValenceWorkUnit::getCompositeDihedralInstructions() const {
  return cdhe_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getCmapInstructions() const {
  return cmap_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint>& ValenceWorkUnit::getInferred14Instructions() const {
  return infr14_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getPositionalRestraintInstructions() const {
  return rposn_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getDistanceRestraintInstructions() const {
  return rbond_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getAngleRestraintInstructions() const {
  return rangl_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getDihedralRestraintInstructions() const {
  return rdihe_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getVirtualSiteInstructions() const {
  return vste_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getSettleGroupInstructions() const {
  return sett_instructions;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint2>& ValenceWorkUnit::getConstraintGroupInstructions() const {
  return cnst_instructions;
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getCompositeBondInstruction(const int index) const {
  return cbnd_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getAngleInstruction(const int index) const {
  return angl_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint3 ValenceWorkUnit::getCompositeDihedralInstruction(const int index) const {
  return cdhe_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getCmapInstruction(const int index) const {
  return cmap_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint ValenceWorkUnit::getInferred14Instruction(const int index) const {
  return infr14_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getPositionalRestraintInstruction(const int index) const {
  return rposn_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getDistanceRestraintInstruction(const int index) const {
  return rbond_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getAngleRestraintInstruction(const int index) const {
  return rangl_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getDihedralRestraintInstruction(const int index) const {
  return rdihe_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getVirtualSiteInstruction(const int index) const {
  return vste_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getSettleGroupInstruction(const int index) const {
  return sett_instructions[index];
}

//-------------------------------------------------------------------------------------------------
uint2 ValenceWorkUnit::getConstraintGroupInstruction(const int index) const {
  return cnst_instructions[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint>& ValenceWorkUnit::getAccumulationFlags(const VwuTask vtask) const {
  switch (vtask) {
  case VwuTask::BOND:
    return acc_bond_energy;
  case VwuTask::ANGL:
    return acc_angl_energy;
  case VwuTask::DIHE:
    return acc_dihe_energy;
  case VwuTask::UBRD:
    return acc_ubrd_energy;
  case VwuTask::CBND:
    return acc_cbnd_energy;
  case VwuTask::CIMP:
    return acc_cimp_energy;
  case VwuTask::CDHE:
    return acc_cdhe_energy;
  case VwuTask::CMAP:
    return acc_cmap_energy;
  case VwuTask::INFR14:
    return acc_infr14_energy;
  case VwuTask::RPOSN:
    return acc_rposn_energy;
  case VwuTask::RBOND:
    return acc_rbond_energy;
  case VwuTask::RANGL:
    return acc_rangl_energy;
  case VwuTask::RDIHE:
    return acc_rdihe_energy;
  case VwuTask::SETTLE:
  case VwuTask::CGROUP:
  case VwuTask::VSITE:
    rtErr("Geometry constraints and virtual sites are not terms pertaining to the energy "
          "accumulation.", "ValenceWorkUnit", "getAccumulationFlags");
  case VwuTask::ALL_TASKS:
    rtErr("Accumulation flags are only defined for a single, specific task.");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint>& ValenceWorkUnit::getAtomUpdateFlags() const {
  return atom_update_mask;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& ValenceWorkUnit::getSimpleTaskList(const VwuTask vtask) const {
  switch (vtask) {
  case VwuTask::BOND:
    return bond_term_list;
  case VwuTask::ANGL:
    return angl_term_list;
  case VwuTask::DIHE:
    return dihe_term_list;
  case VwuTask::UBRD:
    return ubrd_term_list;
  case VwuTask::CIMP:
    return cimp_term_list;
  case VwuTask::CBND:
  case VwuTask::CDHE:
    rtErr("Composite terms cannot be returned by this function.  Access them via the "
          "getComposite<name>TaskList functions.", "ValenceWorkUnit", "getSimpleTaskList");
  case VwuTask::CMAP:
    return cmap_term_list;
  case VwuTask::INFR14:
    return infr14_term_list;
  case VwuTask::RPOSN:
    return rposn_term_list;
  case VwuTask::RBOND:
    return rbond_term_list;
  case VwuTask::RANGL:
    return rangl_term_list;
  case VwuTask::RDIHE:
    return rdihe_term_list;
  case VwuTask::SETTLE:
    return sett_group_list;
  case VwuTask::CGROUP:
    return cnst_group_list;
  case VwuTask::VSITE:
    return virtual_site_list;
  case VwuTask::ALL_TASKS:
    rtErr("A task list is only defined for a single, specific task.");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& ValenceWorkUnit::getCompositeBondTaskList() const {
  return cbnd_term_list;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int2>& ValenceWorkUnit::getCompositeDihedralTaskList() const {
  return cdhe_term_list;
}

//-------------------------------------------------------------------------------------------------
ValenceDelegator* ValenceWorkUnit::getDelegatorPointer() {
  return vdel_pointer;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ValenceWorkUnit::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus* ValenceWorkUnit::getRestraintApparatusPointer() const {
  return ra_pointer;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::setListIndex(const int list_index_in) {
  list_index = list_index_in;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::setAtomLimit(const int new_limit) {
  if (new_limit < imported_atom_count) {
    rtErr("The atom limit cannot be set below the number of atoms currently in a work unit.  "
          "ValenceWorkUnit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + " contains " + std::to_string(imported_atom_count) +
          " atoms and cannot reduce its limit to " + std::to_string(new_limit) + ".",
          "ValenceWorkUnit", "setAtomLimit");
  }
  atom_limit = new_limit;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAtomImport(const int atom_index) {

  // Check that the atom is not already in the work unit
  if (vdel_pointer->checkPresence(atom_index, list_index)) {
    return;
  }
  
  // Add the new atom to the list of atom imports.
  atom_import_list.push_back(atom_index);
  imported_atom_count += 1;

  // Check the appropriate boxes in the delegator.
  vdel_pointer->markAtomAddition(list_index, atom_index);  
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAtomUpdate(const int atom_index) {

  // Check that the atom is on the import list
  if (vdel_pointer->checkPresence(atom_index, list_index) == false) {
    rtErr("Atom with topological index " + std::to_string(atom_index) + " is not imported by "
          "ValenceWorkUnit " + std::to_string(list_index) + ".", "ValenceWorkUnit",
          "addNewAtomUpdate");
  }
  
  // Check the appropriate boxes in the delegator.
  if (vdel_pointer->setUpdateWorkUnit(atom_index, list_index)) {
    atom_update_list.push_back(atom_index);
    updated_atom_count += 1;
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::makeAtomMoveList() {
  std::vector<int> required_atoms;
  for (int i = 0; i < updated_atom_count; i++) {

    // The list of movement partners comprises the particle itself.
    const std::vector<int> tmpv = vdel_pointer->findMovementPartners(atom_update_list[i]);
    required_atoms.insert(required_atoms.end(), tmpv.begin(), tmpv.end());
  }
  atom_move_list = reduceUniqueValues(required_atoms);
  moved_atom_count = atom_move_list.size();
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::sortAtomSets() {
  std::sort(atom_import_list.begin(), atom_import_list.end(), [](int a, int b) { return a < b; });
  std::sort(atom_move_list.begin(), atom_move_list.end(), [](int a, int b) { return a < b; });
  std::sort(atom_update_list.begin(), atom_update_list.end(), [](int a, int b) { return a < b; });

  // Check the sizing of each list
  if (static_cast<int>(atom_import_list.size()) != imported_atom_count) {
    rtErr("Atom import count (" + std::to_string(imported_atom_count) + ") does not reflect the "
          "true length of the imported atom list (" + std::to_string(atom_import_list.size()) +
          ") in work unit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + ".", "ValenceWorkUnit", "sortAtomSets");
  }
  if (static_cast<int>(atom_move_list.size()) != moved_atom_count) {
    rtErr("Moving atom count (" + std::to_string(moved_atom_count) + ") does not reflect the "
          "true length of the atom move list (" + std::to_string(atom_move_list.size()) +
          ") in work unit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + ".", "ValenceWorkUnit", "sortAtomSets");    
  }
  if (static_cast<int>(atom_update_list.size()) != updated_atom_count) {
    rtErr("Atom update count (" + std::to_string(updated_atom_count) + ") does not reflect the "
          "true length of the atom move list (" + std::to_string(atom_update_list.size()) +
          ") in work unit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + ".", "ValenceWorkUnit", "sortAtomSets");    
  }
  
  // Check the atom import list to ensure that all entries are unique.
  for (int i = 1; i < imported_atom_count; i++) {
    if (atom_import_list[i] == atom_import_list[i - 1]) {
      const int ires = ag_pointer->getResidueIndex(atom_import_list[i]);
      rtErr("A duplicate entry is present in the atom import list of work unit " +
            std::to_string(list_index) + " serving topology " + ag_pointer->getFileName() +
            ".  The topological atom index is " + std::to_string(atom_import_list[i]) + ", name " +
            char4ToString(ag_pointer->getAtomName(atom_import_list[i])) + ", residue " +
            char4ToString(ag_pointer->getResidueName(ires)) + " " +
            std::to_string(ag_pointer->getResidueNumber(atom_import_list[i])) + ".",
            "ValenceWorkUnit", "sortAtomSets");
    }
  }

  // Check the movement list to ensure that each atom is present in the import list.
  for (int i = 0; i < moved_atom_count; i++) {
    if (locateValue(atom_import_list, atom_move_list[i],
                    DataOrder::ASCENDING) == imported_atom_count) {
      rtErr("Atom index " + std::to_string(atom_move_list[i]) + "(" +
            char4ToString(ag_pointer->getAtomName(atom_move_list[i])) + ") in topology " +
            ag_pointer->getFileName() + " is scheduled to be moved by work unit " +
            std::to_string(list_index) + " but not present in the import list of " +
            std::to_string(imported_atom_count) + " atoms.", "ValenceWorkUnit", "sortAtomSets");
    }
  }
  for (int i = 0; i < updated_atom_count; i++) {
    if (locateValue(atom_move_list, atom_update_list[i],
                    DataOrder::ASCENDING) == moved_atom_count) {
      rtErr("Atom index " + std::to_string(atom_update_list[i]) + "(" +
            char4ToString(ag_pointer->getAtomName(atom_update_list[i])) + ") in topology " +
            ag_pointer->getFileName() + " is scheduled to be updated and logged by work unit " +
            std::to_string(list_index) + " but not present in the movement list of " +
            std::to_string(moved_atom_count) + " atoms.", "ValenceWorkUnit", "sortAtomSets");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::makeAtomUpdateMask() {
  atom_update_mask.resize((imported_atom_count + uint_bit_count_int - 1) / uint_bit_count_int, 0U);
  size_t import_pos = 0;
  const size_t uac_zu = updated_atom_count;
  for (size_t i = 0; i < uac_zu; i++) {

    // With imported and updated atom lists sorted, it is possible to seek the position of each
    // successive atom update within the list of atom imports by simply advancing a counter.
    while (atom_import_list[import_pos] != atom_update_list[i]) {
      import_pos++;
    }
    accumulateBitmask(&atom_update_mask, import_pos);
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::logActivities() {

  // Extract information from the topology
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag_pointer->getDoublePrecisionConstraintKit();
  const RestraintKit<double, double2, double4> rar = ra_pointer->dpData();

  // Make a "straight table" based on the minimum and maximum imported atom indices.  Most work
  // units will involve a confined sequence of atoms from within the topology.  While it might
  // be inefficient to have each work unit allocate a list mapping the entire topology to one of
  // its internal atoms (a bounded memory requirement in that ValenceWorkUnits are constructed
  // serially, but memory allocation effort still growing as the number of work units times the
  // number of atoms), it is feasible to take an offset and then allocate an array to map the
  // relevant part of the topology to the list of atoms that each work unit does hold.
  const int mapping_offset = atom_import_list[0];
  const int map_length = atom_import_list[imported_atom_count - 1] + 1 - mapping_offset;
  std::vector<int> import_map(map_length, -1);
  for (int i = 0; i < imported_atom_count; i++) {
    import_map[atom_import_list[i] - mapping_offset] = i;
  }
  
  // Detail valence interactions for this work unit.  Add instructions to accumulate various
  // energy terms if this work unit is the first in the list to evaluate a given term.
  const DataOrder order = DataOrder::ASCENDING;
  const int uint_bits = sizeof(uint) * 8;
  const int uint_bits_m1 = uint_bits - 1;
  const std::vector<int> relevant_bonds = vdel_pointer->getBondAffectors(atom_move_list);
  const size_t nbond = relevant_bonds.size();
  acc_bond_energy.resize((nbond + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nbond; pos++) {
    const int bond_idx = relevant_bonds[pos];
    bond_term_list.push_back(bond_idx);
    bond_i_atoms.push_back(import_map[vk.bond_i_atoms[bond_idx] - mapping_offset]);
    bond_j_atoms.push_back(import_map[vk.bond_j_atoms[bond_idx] - mapping_offset]);
    if (vdel_pointer->setBondAccumulatorWorkUnit(bond_idx, list_index)) {
      accumulateBitmask(&acc_bond_energy, pos);
    }
  }
  const std::vector<int> relevant_angls = vdel_pointer->getAngleAffectors(atom_move_list);
  const size_t nangl = relevant_angls.size();
  acc_angl_energy.resize((nangl + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nangl; pos++) {
    const int angl_idx = relevant_angls[pos];
    angl_term_list.push_back(angl_idx);
    angl_i_atoms.push_back(import_map[vk.angl_i_atoms[angl_idx] - mapping_offset]);
    angl_j_atoms.push_back(import_map[vk.angl_j_atoms[angl_idx] - mapping_offset]);
    angl_k_atoms.push_back(import_map[vk.angl_k_atoms[angl_idx] - mapping_offset]);
    if (vdel_pointer->setAngleAccumulatorWorkUnit(angl_idx, list_index)) {
      accumulateBitmask(&acc_angl_energy, pos);
    }
  }
  const std::vector<int> relevant_dihes = vdel_pointer->getDihedralAffectors(atom_move_list);
  const size_t ndihe = relevant_dihes.size();
  acc_dihe_energy.resize((ndihe + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < ndihe; pos++) {
    const int dihe_idx = relevant_dihes[pos];
    dihe_term_list.push_back(dihe_idx);
    dihe_i_atoms.push_back(import_map[vk.dihe_i_atoms[dihe_idx] - mapping_offset]);
    dihe_j_atoms.push_back(import_map[vk.dihe_j_atoms[dihe_idx] - mapping_offset]);
    dihe_k_atoms.push_back(import_map[vk.dihe_k_atoms[dihe_idx] - mapping_offset]);
    dihe_l_atoms.push_back(import_map[vk.dihe_l_atoms[dihe_idx] - mapping_offset]);
    if (vdel_pointer->setDihedralAccumulatorWorkUnit(dihe_idx, list_index)) {
      accumulateBitmask(&acc_dihe_energy, pos);
    }
  }
  const std::vector<int> relevant_ubrds = vdel_pointer->getUreyBradleyAffectors(atom_move_list);
  const size_t nubrd = relevant_ubrds.size();
  acc_ubrd_energy.resize((nubrd + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nubrd; pos++) {
    const int ubrd_idx = relevant_ubrds[pos];
    ubrd_term_list.push_back(ubrd_idx);
    ubrd_i_atoms.push_back(import_map[vk.ubrd_i_atoms[ubrd_idx] - mapping_offset]);
    ubrd_k_atoms.push_back(import_map[vk.ubrd_k_atoms[ubrd_idx] - mapping_offset]);
    if (vdel_pointer->setUreyBradleyAccumulatorWorkUnit(ubrd_idx, list_index)) {
      accumulateBitmask(&acc_ubrd_energy, pos);
    }
  }
  const std::vector<int> relevant_cimps = vdel_pointer->getCharmmImproperAffectors(atom_move_list);
  const size_t ncimp = relevant_cimps.size();
  acc_cimp_energy.resize((ncimp + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < ncimp; pos++) {
    const int cimp_idx = relevant_cimps[pos];
    cimp_term_list.push_back(cimp_idx);
    cimp_i_atoms.push_back(import_map[vk.cimp_i_atoms[cimp_idx] - mapping_offset]);
    cimp_j_atoms.push_back(import_map[vk.cimp_j_atoms[cimp_idx] - mapping_offset]);
    cimp_k_atoms.push_back(import_map[vk.cimp_k_atoms[cimp_idx] - mapping_offset]);
    cimp_l_atoms.push_back(import_map[vk.cimp_l_atoms[cimp_idx] - mapping_offset]);
    if (vdel_pointer->setCharmmImproperAccumulatorWorkUnit(cimp_idx, list_index)) {
      accumulateBitmask(&acc_cimp_energy, pos);
    }
  }
  const std::vector<int> relevant_cmaps = vdel_pointer->getCmapAffectors(atom_move_list);
  const size_t ncmap = relevant_cmaps.size();
  acc_cmap_energy.resize((ncmap + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < ncmap; pos++) {
    const int cmap_idx = relevant_cmaps[pos];
    cmap_term_list.push_back(cmap_idx);
    cmap_i_atoms.push_back(import_map[vk.cmap_i_atoms[cmap_idx] - mapping_offset]);
    cmap_j_atoms.push_back(import_map[vk.cmap_j_atoms[cmap_idx] - mapping_offset]);
    cmap_k_atoms.push_back(import_map[vk.cmap_k_atoms[cmap_idx] - mapping_offset]);
    cmap_l_atoms.push_back(import_map[vk.cmap_l_atoms[cmap_idx] - mapping_offset]);
    cmap_m_atoms.push_back(import_map[vk.cmap_m_atoms[cmap_idx] - mapping_offset]);
    if (vdel_pointer->setCmapAccumulatorWorkUnit(cmap_idx, list_index)) {
      accumulateBitmask(&acc_cmap_energy, pos);
    }
  }
  const std::vector<int> relevant_infrs = vdel_pointer->getInferred14Affectors(atom_move_list);
  const size_t ninfr = relevant_infrs.size();
  acc_infr14_energy.resize((ninfr + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < ninfr; pos++) {
    const int infr_idx = relevant_infrs[pos];
    infr14_term_list.push_back(infr_idx);
    infr14_i_atoms.push_back(import_map[vk.infr14_i_atoms[infr_idx] - mapping_offset]);
    infr14_l_atoms.push_back(import_map[vk.infr14_l_atoms[infr_idx] - mapping_offset]);
    if (vdel_pointer->setInferred14AccumulatorWorkUnit(infr_idx, list_index)) {
      accumulateBitmask(&acc_infr14_energy, pos);
    }
  }

  // Detail the restraint terms for this work unit to manage
  const std::vector<int> relevant_rposns =
    vdel_pointer->getPositionalRestraintAffectors(atom_move_list);
  const size_t nrposn = relevant_rposns.size();
  acc_rposn_energy.resize((nrposn + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nrposn; pos++) {
    const int rposn_idx = relevant_rposns[pos];
    rposn_term_list.push_back(rposn_idx);
    rposn_atoms.push_back(import_map[rar.rposn_atoms[rposn_idx] - mapping_offset]);
    if (vdel_pointer->setPositionalRestraintAccumulatorWorkUnit(rposn_idx, list_index)) {
      accumulateBitmask(&acc_rposn_energy, pos);
    }
  }
  const std::vector<int> relevant_rbonds =
    vdel_pointer->getDistanceRestraintAffectors(atom_move_list);
  const size_t nrbond = relevant_rbonds.size();
  acc_rbond_energy.resize((nrbond + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nrbond; pos++) {
    const int rbond_idx = relevant_rbonds[pos];
    rbond_term_list.push_back(rbond_idx);
    rbond_i_atoms.push_back(import_map[rar.rbond_i_atoms[rbond_idx] - mapping_offset]);
    rbond_j_atoms.push_back(import_map[rar.rbond_j_atoms[rbond_idx] - mapping_offset]);
    if (vdel_pointer->setDistanceRestraintAccumulatorWorkUnit(rbond_idx, list_index)) {
      accumulateBitmask(&acc_rbond_energy, pos);
    }
  }
  const std::vector<int> relevant_rangls =
    vdel_pointer->getAngleRestraintAffectors(atom_move_list);
  const size_t nrangl = relevant_rangls.size();
  acc_rangl_energy.resize((nrangl + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nrangl; pos++) {
    const int rangl_idx = relevant_rangls[pos];
    rangl_term_list.push_back(rangl_idx);
    rangl_i_atoms.push_back(import_map[rar.rangl_i_atoms[rangl_idx] - mapping_offset]);
    rangl_j_atoms.push_back(import_map[rar.rangl_j_atoms[rangl_idx] - mapping_offset]);
    rangl_k_atoms.push_back(import_map[rar.rangl_k_atoms[rangl_idx] - mapping_offset]);
    if (vdel_pointer->setAngleRestraintAccumulatorWorkUnit(rangl_idx, list_index)) {
      accumulateBitmask(&acc_rangl_energy, pos);
    }
  }
  const std::vector<int> relevant_rdihes =
    vdel_pointer->getDihedralRestraintAffectors(atom_move_list);
  const size_t nrdihe = relevant_rdihes.size();
  acc_rdihe_energy.resize((nrdihe + uint_bits_m1) / uint_bits);
  for (size_t pos = 0; pos < nrdihe; pos++) {
    const int rdihe_idx = relevant_rdihes[pos];
    rdihe_term_list.push_back(rdihe_idx);
    rdihe_i_atoms.push_back(import_map[rar.rdihe_i_atoms[rdihe_idx] - mapping_offset]);
    rdihe_j_atoms.push_back(import_map[rar.rdihe_j_atoms[rdihe_idx] - mapping_offset]);
    rdihe_k_atoms.push_back(import_map[rar.rdihe_k_atoms[rdihe_idx] - mapping_offset]);
    rdihe_l_atoms.push_back(import_map[rar.rdihe_l_atoms[rdihe_idx] - mapping_offset]);
    if (vdel_pointer->setDihedralRestraintAccumulatorWorkUnit(rdihe_idx, list_index)) {
      accumulateBitmask(&acc_rdihe_energy, pos);
    }
  }

  // Detail virtual sites for this work unit to manage
  const std::vector<int> relevant_vstes = vdel_pointer->getVirtualSiteAffectors(atom_move_list);
  const size_t nvste = relevant_vstes.size();
  for (size_t pos = 0; pos < nvste; pos++) {
    const int vste_idx = relevant_vstes[pos];
    virtual_site_list.push_back(vste_idx);
    vsite_atoms.push_back(import_map[vsk.vs_atoms[vste_idx] - mapping_offset]);
    vsite_frame1_atoms.push_back(import_map[vsk.frame1_idx[vste_idx] - mapping_offset]);
    vsite_frame2_atoms.push_back(import_map[vsk.frame2_idx[vste_idx] - mapping_offset]);
    const int frame_param_idx = vsk.vs_param_idx[vste_idx];
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[frame_param_idx])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      vsite_frame3_atoms.push_back(-1);
      vsite_frame4_atoms.push_back(-1);
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      vsite_frame3_atoms.push_back(import_map[vsk.frame3_idx[vste_idx] - mapping_offset]);
      vsite_frame4_atoms.push_back(-1);
      break;
    case VirtualSiteKind::FIXED_4:
      vsite_frame3_atoms.push_back(import_map[vsk.frame3_idx[vste_idx] - mapping_offset]);
      vsite_frame4_atoms.push_back(import_map[vsk.frame4_idx[vste_idx] - mapping_offset]);
      break;
    }
  }

  // Detail the constraint groups for this work unit
  const std::vector<int> relevant_setts = vdel_pointer->getSettleGroupAffectors(atom_move_list);
  const size_t nsett = relevant_setts.size();
  for (size_t pos = 0; pos < nsett; pos++) {
    const int sett_idx = relevant_setts[pos];
    sett_group_list.push_back(sett_idx);
    sett_ox_atoms.push_back(import_map[cnk.settle_ox_atoms[sett_idx] - mapping_offset]);
    sett_h1_atoms.push_back(import_map[cnk.settle_h1_atoms[sett_idx] - mapping_offset]);
    sett_h2_atoms.push_back(import_map[cnk.settle_h2_atoms[sett_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_cnsts =
    vdel_pointer->getConstraintGroupAffectors(atom_move_list);
  const size_t ncnst = relevant_cnsts.size();
  cnst_group_bounds.push_back(0);
  for (size_t pos = 0; pos < ncnst; pos++) {
    const int cnst_idx = relevant_cnsts[pos];
    cnst_group_list.push_back(cnst_idx);
    for (int i = cnk.group_bounds[cnst_idx]; i < cnk.group_bounds[cnst_idx + 1]; i++) {
      cnst_group_atoms.push_back(import_map[cnk.group_list[i] - mapping_offset]);
    }
    cnst_group_bounds.push_back(cnst_group_atoms.size());
  }

  // Log the counts of each basic task
  bond_term_count   = bond_term_list.size();
  angl_term_count   = angl_term_list.size();
  dihe_term_count   = dihe_term_list.size();
  ubrd_term_count   = ubrd_term_list.size();
  cimp_term_count   = cimp_term_list.size();
  cmap_term_count   = cmap_term_list.size();
  infr14_term_count = infr14_term_list.size();
  rposn_term_count  = rposn_term_list.size();
  rbond_term_count  = rbond_term_list.size();
  rangl_term_count  = rangl_term_list.size();
  rdihe_term_count  = rdihe_term_list.size();
  cnst_group_count  = cnst_group_list.size();
  sett_group_count  = sett_group_list.size();
  vste_count        = virtual_site_list.size();

  // Concatenate the bond and Urey-Bradley term lists into the composite bond list.
  cbnd_term_list = bond_term_list;
  cbnd_term_list.insert(cbnd_term_list.end(), ubrd_term_list.begin(), ubrd_term_list.end());
  cbnd_i_atoms = bond_i_atoms;
  cbnd_jk_atoms = bond_j_atoms;
  cbnd_i_atoms.insert(cbnd_i_atoms.end(), ubrd_i_atoms.begin(), ubrd_i_atoms.end());
  cbnd_jk_atoms.insert(cbnd_jk_atoms.end(), ubrd_k_atoms.begin(), ubrd_k_atoms.end());
  cbnd_is_ubrd.resize(cbnd_term_list.size(), true);
  for (int i = 0; i < bond_term_count; i++) {
    cbnd_is_ubrd[i] = false;
  }
  cbnd_term_count = cbnd_term_list.size();
  acc_cbnd_energy.resize((cbnd_term_count + uint_bits_m1) / uint_bits);
  for (int i = 0; i < cbnd_term_count; i++) {
    if (cbnd_is_ubrd[i]) {
      if (vdel_pointer->getUreyBradleyAccumulatorWorkUnit(cbnd_term_list[i]) == list_index) {
        accumulateBitmask(&acc_cbnd_energy, i);        
      }
    }
    else {
      if (vdel_pointer->getBondAccumulatorWorkUnit(cbnd_term_list[i]) == list_index) {
        accumulateBitmask(&acc_cbnd_energy, i);        
      }
    }
  }
  
  // Compute the number of composite dihedrals: the dihedrals already subsume a great deal of the
  // 1:4 attenuated interactions (probably the entirety of these interactions if the system has no
  // virtual sites).  However, many dihedrals have additional cosine terms overlaid on the same
  // four atoms.  Bundle these terms in pairs, if possible, and eat the blank calculation if there
  // is no secondary dihedral to pair in any one case.  For a typical simulation with an Amber
  // force field, this will cut down on the total number of dihedral computations by about 40%,
  // and bundling the 1:4 interactions will likewise reduce the memory traffic.  The overall
  // reduction in memory bandwidth for reading instructions (relative to an approach that reads
  // each dihedral and 1:4 term separately, with a 64-bit instruction for each dihedral and a
  // 32-bit instruction for each 1:4 term) is roughly 30%, on top of other optimizations that
  // streamline the information access.  The composite dihedrals will replace the lists of standard
  // dihedral instructions and CHARMM impropers, fusing them into one.
  std::vector<int2> dihe_presence_bounds(imported_atom_count);
  std::vector<bool> presence_found(imported_atom_count, false);
  std::vector<bool> is_improper(dihe_term_count, false);
  for (int pos = 0; pos < dihe_term_count; pos++) {
    const int term_idx = dihe_term_list[pos];
    switch (static_cast<TorsionKind>(vk.dihe_modifiers[term_idx].w)) {
    case TorsionKind::PROPER_NO_14:
    case TorsionKind::PROPER:
      break;
    case TorsionKind::IMPROPER_NO_14:
    case TorsionKind::IMPROPER:
      is_improper[pos] = true;
    }
    const int i_atom = dihe_i_atoms[pos];
    const int j_atom = dihe_j_atoms[pos];
    const int k_atom = dihe_k_atoms[pos];
    const int l_atom = dihe_l_atoms[pos];
    if (presence_found[i_atom]) {
      dihe_presence_bounds[i_atom].y = pos;
    }
    else {
      dihe_presence_bounds[i_atom].x = pos;
      presence_found[i_atom] = true;
    }
    if (presence_found[j_atom]) {
      dihe_presence_bounds[j_atom].y = pos;
    }
    else {
      dihe_presence_bounds[j_atom].x = pos;
      presence_found[j_atom] = true;
    }
    if (presence_found[k_atom]) {
      dihe_presence_bounds[k_atom].y = pos;
    }
    else {
      dihe_presence_bounds[k_atom].x = pos;
      presence_found[k_atom] = true;
    }
    if (presence_found[l_atom]) {
      dihe_presence_bounds[l_atom].y = pos;
    }
    else {
      dihe_presence_bounds[l_atom].x = pos;
      presence_found[l_atom] = true;
    }
  }

  // Find all dihedrals that can be paired when acting on the same atoms
  std::vector<bool> coverage(dihe_term_count, false);
  for (int pos = 0; pos < dihe_term_count; pos++) {

    // Filter out terms with extermely large parameter indices--these must appear in their own
    // composite terms, as the only dihedral.  This filtering applies only to proper dihedrals,
    // which are more diverse in their parameter combinations than CHARMM improper dihedrals.
    if (vk.dihe_param_idx[dihe_term_list[pos]] >= 65535) {
      continue;
    }

    // Filter out impropers and skip terms that have already been included in a combined item.
    if (coverage[pos] || is_improper[pos]) {
      continue;
    }
    const int i_atom = dihe_i_atoms[pos];
    const int j_atom = dihe_j_atoms[pos];
    const int k_atom = dihe_k_atoms[pos];
    const int l_atom = dihe_l_atoms[pos];
    int min_np = std::max(dihe_presence_bounds[i_atom].x, dihe_presence_bounds[j_atom].x);
    min_np = std::max(min_np, dihe_presence_bounds[k_atom].x);
    min_np = std::max(min_np, dihe_presence_bounds[l_atom].x);
    min_np = std::max(min_np, pos + 1);
    int max_np = std::min(dihe_presence_bounds[i_atom].y, dihe_presence_bounds[j_atom].y);
    max_np = std::min(max_np, dihe_presence_bounds[k_atom].y);
    max_np = std::min(max_np, dihe_presence_bounds[l_atom].y);
    for (int npos = min_np; npos < max_np; npos++) {
      if (coverage[npos] || is_improper[npos]) {
        continue;
      }
      if ((dihe_i_atoms[npos] == dihe_i_atoms[pos] && dihe_j_atoms[npos] == dihe_j_atoms[pos] &&
           dihe_k_atoms[npos] == dihe_k_atoms[pos] && dihe_l_atoms[npos] == dihe_l_atoms[pos]) ||
          (dihe_i_atoms[npos] == dihe_l_atoms[pos] && dihe_j_atoms[npos] == dihe_k_atoms[pos] &&
           dihe_k_atoms[npos] == dihe_j_atoms[pos] && dihe_l_atoms[npos] == dihe_i_atoms[pos])) {

        // This dihedral applies to the same atoms.  Make a composite dihedral, then bail out.
        cdhe_term_list.push_back({dihe_term_list[pos], dihe_term_list[npos]});
        cdhe_is_cimp.push_back(false);
        cdhe_i_atoms.push_back(dihe_i_atoms[pos]);
        cdhe_j_atoms.push_back(dihe_j_atoms[pos]);
        cdhe_k_atoms.push_back(dihe_k_atoms[pos]);
        cdhe_l_atoms.push_back(dihe_l_atoms[pos]);
        coverage[pos] = true;
        coverage[npos] = true;
        break;
      }
    }
  }

  // Add the remaining dihedrals to the composite dihedral list as singlets
  for (int pos = 0; pos < dihe_term_count; pos++) {
    if (coverage[pos]) {
      continue;
    }
    cdhe_term_list.push_back({dihe_term_list[pos], -1});
    cdhe_is_cimp.push_back(false);
    cdhe_i_atoms.push_back(dihe_i_atoms[pos]);
    cdhe_j_atoms.push_back(dihe_j_atoms[pos]);
    cdhe_k_atoms.push_back(dihe_k_atoms[pos]);
    cdhe_l_atoms.push_back(dihe_l_atoms[pos]);
    coverage[pos] = true;
  }
  
  // Add CHARMM improper dihedrals to the composite list as more singlet terms
  for (int pos = 0; pos < cimp_term_count; pos++) {
    cdhe_term_list.push_back({cimp_term_list[pos], -1});
    cdhe_is_cimp.push_back(true);
    cdhe_i_atoms.push_back(cimp_i_atoms[pos]);
    cdhe_j_atoms.push_back(cimp_j_atoms[pos]);
    cdhe_k_atoms.push_back(cimp_k_atoms[pos]);
    cdhe_l_atoms.push_back(cimp_l_atoms[pos]);
    coverage[pos] = true;
  }
  cdhe_term_count = cdhe_term_list.size();

  // Mark which of the composite dihedrals this work unit will accumulate into the energy totals.
  // If the composite term involves a pair of dihedrals, the work unit must accumulate both of them
  // in order to accumulate the composite term.  If the dihedrals somehow are accumulated by
  // different work units (which should be impossible, as work unit are built sequentially and all
  // dihedrals involving any four atoms will be present in any work unit that moves one of those
  // atoms), set this work unit to accumulate the potential if it controls the dihedral with the
  // lower topological term index.  
  acc_cdhe_energy.resize((cdhe_term_count + uint_bits_m1) / uint_bits);
  for (int pos = 0; pos < cdhe_term_count; pos++) {

    // Check that no composite term lists the same dihedral twice
    const int xterm = cdhe_term_list[pos].x;
    const int yterm = cdhe_term_list[pos].y;
    if (xterm == yterm) {
      rtErr("Composite term " + std::to_string(pos) + " in ValenceWorkUnit " +
            std::to_string(list_index) + " has both dihedrals referencing term index " +
            std::to_string(xterm) + ".", "ValenceWorkUnit", "logActivities");
    }
    if (cdhe_is_cimp[pos]) {
      if (vdel_pointer->getCharmmImproperAccumulatorWorkUnit(xterm) == list_index) {
        accumulateBitmask(&acc_cdhe_energy, pos);
      }
    }
    else if (xterm >= 0 && yterm >= 0) {
      if (vdel_pointer->getDihedralAccumulatorWorkUnit(xterm) == list_index &&
          vdel_pointer->getDihedralAccumulatorWorkUnit(yterm) == list_index) {
        accumulateBitmask(&acc_cdhe_energy, pos);
      }
      else if (xterm < yterm &&
               vdel_pointer->getDihedralAccumulatorWorkUnit(xterm) == list_index) {
        accumulateBitmask(&acc_cdhe_energy, pos);
      }
      else if (yterm < xterm &&
               vdel_pointer->getDihedralAccumulatorWorkUnit(xterm) == list_index) {
        accumulateBitmask(&acc_cdhe_energy, pos);
      }
    }
    else if (xterm >= 0  &&
             vdel_pointer->getDihedralAccumulatorWorkUnit(xterm) == list_index) {
      accumulateBitmask(&acc_cdhe_energy, pos);
    }
  }

  // Store instruction sets based on the original topology only.  If the instructions need to be
  // realigned to the parameter tables in an AtomGraphSynthesis, that can be done with external
  // calls to the same functions, supplying the necessary parameter translation tables.
  storeCompositeBondInstructions();
  storeAngleInstructions();
  storeCompositeDihedralInstructions();
  storeCmapInstructions();
  storeInferred14Instructions();
  storePositionalRestraintInstructions();
  storeDistanceRestraintInstructions();
  storeAngleRestraintInstructions();
  storeDihedralRestraintInstructions();
  storeVirtualSiteInstructions();
  storeSettleGroupInstructions();
  storeConstraintGroupInstructions();
}

//-------------------------------------------------------------------------------------------------
int calculateValenceWorkUnitSize(const int* atom_counts, const int system_count) {

  // The maximum topology size is relevant--valence work units pertain to one and only one
  // topology due to considerations for importing atoms and accumulating energy.
  int min_inefficiency, best_size;
  for (int vs = minimum_valence_work_unit_atoms; vs <= maximum_valence_work_unit_atoms; vs *= 2) {
    const int unique_atom_coverage = vs - 12;
    int inefficiency = 0;
    for (int i = 0; i < system_count; i++) {
      if (atom_counts[i] < vs) {
        inefficiency = vs - atom_counts[i];
      }
      else {
        const int tnvwu = (atom_counts[i] + unique_atom_coverage - 1) / unique_atom_coverage;
        inefficiency += (tnvwu * vs) - atom_counts[i];
      }
    }
    if (vs == minimum_valence_work_unit_atoms || inefficiency < min_inefficiency) {
      min_inefficiency = inefficiency;
      best_size = vs;
    }
  }
  return best_size;
}

//-------------------------------------------------------------------------------------------------
int calculateValenceWorkUnitSize(const std::vector<int> &atom_counts) {
  return calculateValenceWorkUnitSize(atom_counts.data(), atom_counts.size());
}

//-------------------------------------------------------------------------------------------------
int calculateValenceWorkUnitSize(const Hybrid<int> &atom_counts) {
  return calculateValenceWorkUnitSize(atom_counts.data(), atom_counts.size());
}

//-------------------------------------------------------------------------------------------------
int calculateValenceWorkUnitSize(const int* atom_counts, const int system_count,
                                 const int sm_count, ValenceKernelSize *kwidth) {
  int best_size;
  ValenceKernelSize best_layout;
  double best_efficiency = 0.0;
  const int vwu_size_incr = maximum_valence_work_unit_atoms / 32;
  const double total_atoms = sum<int>(atom_counts, system_count);
  for (int vs = 2 * vwu_size_incr; vs <= maximum_valence_work_unit_atoms; vs += vwu_size_incr) {
    
    // Estimate the number of valence work units that would be needed to cover the workload with
    // this size of work unit.  The overlap is conservative, implying slightly more work units than
    // may be necessary, to ensure that the algorithm is unlikely to overrun break points at which
    // a given card might end up with one extra block to run, e.g. slots for 320 blocks and 321 in
    // the workload.
    int nvwu_est = 0;
    const double dvs = vs;
    double average_vwu_atoms = 0.0;
    const int unique_atom_coverage = vs - 12;
    for (int i = 0; i < system_count; i++) {
      if (atom_counts[i] < vs) {
        nvwu_est++;
        average_vwu_atoms += static_cast<double>(atom_counts[i]);
      }
      else {
        const int tnvwu = (atom_counts[i] + unique_atom_coverage - 1) / unique_atom_coverage;
        nvwu_est += tnvwu;

        // The intrinsic efficiency incorporates the new number of work units, less one, at the
        // estimated unique atom coverage per work unit and one work unit with the remaining number
        // of atoms.
        const int last_vwu_atoms_est = atom_counts[i] - ((tnvwu - 1) * unique_atom_coverage);
        average_vwu_atoms += static_cast<double>(((tnvwu - 1) * vs) + last_vwu_atoms_est + 12);
      }
    }
    const double intrinsic_eff = total_atoms / average_vwu_atoms;
    average_vwu_atoms /= static_cast<double>(nvwu_est);
    
    // Loop over the possible kernel block sizes that could accommodate work units of this size.
    // Find the optimal efficiency.
    const int vks_max = static_cast<int>(ValenceKernelSize::SM);
    ValenceKernelSize best_kwidth = ValenceKernelSize::XL;
    for (int i = 0; i <= vks_max; i++) {

      // Pair each enumeration to a specific kernel size, to future-proof against changes in the
      // order of the list of ValenceKernelSize values.
      int test_katoms, block_multiplier;
      const ValenceKernelSize iwidth  = static_cast<ValenceKernelSize>(i);
      switch (iwidth) {
      case ValenceKernelSize::XL:
        test_katoms = maximum_valence_work_unit_atoms;
        block_multiplier = 2;
        break;
      case ValenceKernelSize::LG:
        test_katoms = half_valence_work_unit_atoms;
        block_multiplier = 4;
        break;
      case ValenceKernelSize::MD:
        test_katoms = quarter_valence_work_unit_atoms;
        block_multiplier = 8;
        break;
      case ValenceKernelSize::SM:
        test_katoms = eighth_valence_work_unit_atoms;
        block_multiplier = 16;
        break;
      }
      if (vs <= test_katoms) {
        const int nblocks = sm_count * block_multiplier;
        const int ncyc = (nvwu_est + nblocks - 1) / nblocks;

        // Use the sqrt function to represent the fact that most atoms will have more than one
        // associated dihedral angle, and other terms.  Underfull work units will still obtain
        // thread utilization closer to the optimum than the overall percentage of the atom
        // capacity might suggest.
        const double block_eff = sqrt(average_vwu_atoms / test_katoms) *
                                 (static_cast<double>(nvwu_est) /
                                  static_cast<double>(ncyc * nblocks));
        const double overall_eff = block_eff * intrinsic_eff;
        if (overall_eff > best_efficiency) {
          best_efficiency = overall_eff;
          best_size = vs;
          best_layout = iwidth;
        }
      }
    }
  }
  *kwidth = best_layout;
  return best_size;
}

//-------------------------------------------------------------------------------------------------
int calculateValenceWorkUnitSize(const std::vector<int> &atom_counts, const int sm_count,
                                 ValenceKernelSize *kwidth) {
  return calculateValenceWorkUnitSize(atom_counts.data(), atom_counts.size(), sm_count, kwidth);
}

//-------------------------------------------------------------------------------------------------
int calculateValenceWorkUnitSize(const Hybrid<int> &atom_counts, const int sm_count,
                                 ValenceKernelSize *kwidth) {
  return calculateValenceWorkUnitSize(atom_counts.data(), atom_counts.size(), sm_count, kwidth);
}

//-------------------------------------------------------------------------------------------------
std::vector<ValenceWorkUnit> buildValenceWorkUnits(const AtomGraph *ag,
                                                   const RestraintApparatus *ra,
                                                   const int max_atoms_per_vwu) {
  ValenceDelegator vdel(ag, ra);
  return buildValenceWorkUnits(&vdel, max_atoms_per_vwu);
}

//-------------------------------------------------------------------------------------------------
std::vector<ValenceWorkUnit> buildValenceWorkUnits(ValenceDelegator *vdel,
                                                   const int max_atoms_per_vwu) {
  std::vector<ValenceWorkUnit> result;

  // Assign one and only one unit to log the moves to each atom, including updating the positions
  // of constrained atoms and virtual sites.  By construction, the ValenceWorkUnit objects will
  // have a modicum of extra space to store additional atoms.  In order to determine the position
  // of an atom that is part of a constraint group, the work unit must possess and move all atoms
  // in the constraint group.  Likewise, in order to place a virtual site, the work unit must have
  // current locations of all frame atoms for the virtual site in their updated positions.  In
  // order to move any atom, all forces on that atom must be known.  The hierarchy is therefore:
  //
  // Atom to move and log
  //  - Atoms that contribute directly to forces on the atom to move
  //  - Other atoms in constraint group
  //     - Atoms that contribute to forces on the constraint group atoms
  //  - Frame atoms, if the atom is a virtual site
  //     - Atoms that contribute to forces on the frame atoms
  //
  // Atom A can contribute to forces on some other atom B by being part of a valence term, NMR
  // restraint, or a virtual site frame that also involves B.  The hierarchy cannot chain
  // indefinitely, as only one cycle of force computation (including transmission from a virtual
  // site) can occur before moving atoms.
  const AtomGraph *ag_ptr = vdel->getTopologyPointer();
  const RestraintApparatus *ra_ptr = vdel->getRestraintApparatusPointer();
  std::vector<int> tvwu_coverage(ag_ptr->getAtomCount(), 0);
  while (vdel->getFirstUnassignedAtom() < ag_ptr->getAtomCount()) {
    const int n_units = result.size();
    result.emplace_back(vdel, &tvwu_coverage, n_units, vdel->getFirstUnassignedAtom(),
                        max_atoms_per_vwu);
  }
  
  // Loop once more over the update list and construct the atom movement list.  The movement list
  // is a superset of all the atoms that the work unit will update, and the import list is a
  // superset of the movement list.  This completes the atom sets that define any work unit.
  // Sort each of the sets.
  const int nvwu = result.size();
  for (int i = 0; i < nvwu; i++) {
    result[i].makeAtomMoveList();
    result[i].sortAtomSets();
    result[i].makeAtomUpdateMask();
  }
  
  // With the atom update assignments of each work unit known and the import list furnishing any
  // additional required dependencies (halo atoms), the task is now to loop over all updates and
  // trace back the force terms that must be computed.
  for (int i = 0; i < nvwu; i++) {
    result[i].logActivities();
  }  
  return result;
}

} // namespace synthesis
} // namespace stormm
