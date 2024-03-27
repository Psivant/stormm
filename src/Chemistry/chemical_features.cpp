#include <algorithm>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/scaling.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/mixed_types.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Math/summation.h"
#include "Parsing/parse.h"
#include "Topology/atomgraph_analysis.h"
#include "Topology/topology_util.h"
#include "UnitTesting/approx.h"
#include "chemical_features.h"
#include "indigo.h"

namespace stormm {
namespace chemistry {
  
using card::HybridKind;
using constants::PrecisionModel;
using diskutil::getBaseName;
using stmath::accumulateBitmask;
using stmath::addScalarToVector;
using stmath::numberSeriesToBitMask;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::readBitFromMask;
using stmath::reduceUniqueValues;
using stmath::sum;
using stmath::unsetBitInMask;
using numerics::globalpos_scale_nonoverflow_bits;
using parse::char4ToString;
using testing::Approx;
using topology::colorConnectivity;
using topology::isBonded;
using topology::selectRotatingAtoms;
using topology::TorsionKind;
using trajectory::CoordinateFrameReader;

//-------------------------------------------------------------------------------------------------
BondedNode::BondedNode() :
    previous_atom_index{-1}, previous_node_index{-1}, atom_index{-1}, layer_index{-1},
    root_bond_order{0.0}, branch_count{0}, branch_atoms{nullptr}, rings_completed{0}
{}

//-------------------------------------------------------------------------------------------------
void BondedNode::setBranchPointer(std::vector<int> *vi, const size_t pos,
                                  const size_t max_branches) {

  // Check that the vector has enough data to accept the pointer
  if (vi->size() < (pos + 1) * max_branches) {
    rtErr("Storage vector does not have sufficient space.", "BondedNode");
  }
  branch_atoms = &(vi->data()[pos * max_branches]);
}

//-------------------------------------------------------------------------------------------------
void BondedNode::addToTree(const int previous_in, const int current_atom, const int current_layer,
                           const NonbondedKit<double> &nbk,
                           const std::vector<bool> &valid_atom_mask) {
  previous_atom_index = previous_in;
  atom_index = current_atom;
  layer_index = current_layer;
  const int excl_start = nbk.nb12_bounds[current_atom];
  const int link_candidate_count = nbk.nb12_bounds[current_atom + 1] - excl_start;
  int j = 0;
  for (int i = 0; i < link_candidate_count; i++) {
    const int candidate_atom = nbk.nb12x[i + excl_start];
    if (candidate_atom != previous_in && valid_atom_mask[candidate_atom]) {
      branch_atoms[j] = candidate_atom;
      j++;
    }
  }
  branch_count = j;
}

//-------------------------------------------------------------------------------------------------
void BondedNode::addToTree(const int previous_in, const int current_atom, const int current_layer,
                           const int previous_node_in, const NonbondedKit<double> &nbk,
                           const std::vector<bool> &valid_atom_mask) {
  previous_node_index = previous_node_in;
  addToTree(previous_in, current_atom, current_layer, nbk, valid_atom_mask);
}

//-------------------------------------------------------------------------------------------------
void BondedNode::addBondOrder(const ValenceKit<double> &vk, const Hybrid<double> &bond_orders) {
  for (int i = vk.bond_asgn_bounds[atom_index]; i < vk.bond_asgn_bounds[atom_index + 1]; i++) {
    if (vk.bond_asgn_atoms[i] == previous_atom_index) {
      root_bond_order = bond_orders.readHost(vk.bond_asgn_terms[i]);
      return;
    }
  }
  for (int i = vk.bond_asgn_bounds[previous_atom_index];
       i < vk.bond_asgn_bounds[previous_atom_index + 1]; i++) {
    if (vk.bond_asgn_atoms[i] == atom_index) {
      root_bond_order = bond_orders.readHost(vk.bond_asgn_terms[i]);
      return;
    }
  }
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getPreviousAtom() const {
  return previous_atom_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getPreviousNode() const {
  return previous_node_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getAtom() const {
  return atom_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getLayer() const {
  return layer_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getBranchCount() const {
  return branch_count;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getBranchAtom(const int index) const {
  return branch_atoms[index];
}

//-------------------------------------------------------------------------------------------------
int BondedNode::findBranchAtom(const int search_index) const {
  for (int i = 0; i < branch_count; i++) {
    if (branch_atoms[i] == search_index) {
      return i;
    }
  }
  rtErr("The atom with topology index " + std::to_string(search_index + 1) + " is no branch of " +
        "atom " + std::to_string(atom_index + 1) + ".", "BondedNode", "findBranchAtom");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
uint BondedNode::getRingCompletion(const int branch_index) const {
  return ((rings_completed >> branch_index) & 0x1);
}

//-------------------------------------------------------------------------------------------------
double BondedNode::getRootBondOrder() const {
  return root_bond_order;
}
  
//-------------------------------------------------------------------------------------------------
void BondedNode::setRingCompletion(const int branch_index) {
  rings_completed |= (0x1 << branch_index);
}

//-------------------------------------------------------------------------------------------------
void BondedNode::wipeRingCompletion() {
  rings_completed = 0U;
}

//-------------------------------------------------------------------------------------------------
IsomerPlan::IsomerPlan(const ConformationEdit motion_in,
                       const ChiralInversionProtocol chiral_plan_in, const int root_atom_in,
                       const int pivot_atom_in, const std::vector<int> &moving_atoms_in,
                       const AtomGraph* ag_pointer_in) :
    motion{motion_in}, chiral_plan{chiral_plan_in}, root_atom{root_atom_in},
    pivot_atom{pivot_atom_in}, root_handle{-1}, pivot_handle{-1}, moving_atoms{moving_atoms_in},
    ag_pointer{const_cast<AtomGraph*>(ag_pointer_in)}
{
  switch (motion) {
  case ConformationEdit::BOND_ROTATION:
  case ConformationEdit::CIS_TRANS_FLIP:
    {
      // Find the root atom's handle
      const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
      const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
      const int root_min_bound = nbk.nb12_bounds[root_atom];
      const int root_max_bound = nbk.nb12_bounds[root_atom + 1];
      int max_root_z = -1;
      for (int i = root_min_bound; i < root_max_bound; i++) {
        if (nbk.nb12x[i] != pivot_atom && cdk.z_numbers[nbk.nb12x[i]] > max_root_z) {
          max_root_z = cdk.z_numbers[nbk.nb12x[i]];
          root_handle = nbk.nb12x[i];
        }
      }

      // Find the pivot atom's handle
      const int pivot_min_bound = nbk.nb12_bounds[pivot_atom];
      const int pivot_max_bound = nbk.nb12_bounds[pivot_atom + 1];
      int max_pivot_z = -1;
      for (int i = pivot_min_bound; i < pivot_max_bound; i++) {
        if (nbk.nb12x[i] != root_atom && cdk.z_numbers[nbk.nb12x[i]] > max_pivot_z) {
          max_pivot_z = cdk.z_numbers[nbk.nb12x[i]];
          pivot_handle = nbk.nb12x[i];
        }
      }
    }
    break;
  case ConformationEdit::CHIRAL_INVERSION:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
IsomerPlan::IsomerPlan(const ConformationEdit motion_in, const int root_atom_in,
                       const int pivot_atom_in, const std::vector<int> &moving_atoms_in,
                       const AtomGraph* ag_pointer_in) :
    IsomerPlan(motion_in, ChiralInversionProtocol::DO_NOT_INVERT, root_atom_in, pivot_atom_in,
               moving_atoms_in, ag_pointer_in)
{}

//-------------------------------------------------------------------------------------------------
IsomerPlan::IsomerPlan(const ConformationEdit motion_in,
                       const ChiralInversionProtocol chiral_plan_in, const int root_atom_in,
                       const int pivot_atom_in, const AtomGraph* ag_pointer_in) :
    IsomerPlan(motion_in, chiral_plan_in, root_atom_in, pivot_atom_in, {}, ag_pointer_in)
{}

//-------------------------------------------------------------------------------------------------
IsomerPlan::IsomerPlan(const ConformationEdit motion_in, const int root_atom_in,
                       const int pivot_atom_in, const AtomGraph* ag_pointer_in) :
    IsomerPlan(motion_in, ChiralInversionProtocol::DO_NOT_INVERT, root_atom_in, pivot_atom_in, {},
               ag_pointer_in)
{}

//-------------------------------------------------------------------------------------------------
ConformationEdit IsomerPlan::getMotion() const {
  return motion;
}

//-------------------------------------------------------------------------------------------------
ChiralInversionProtocol IsomerPlan::getChiralPlan() const {
  return chiral_plan;
}

//-------------------------------------------------------------------------------------------------
int IsomerPlan::getRootAtom() const {
  return root_atom;
}

//-------------------------------------------------------------------------------------------------
int IsomerPlan::getPivotAtom() const {
  return pivot_atom;
}

//-------------------------------------------------------------------------------------------------
int IsomerPlan::getRootHandle() const {
  return root_handle;
}

//-------------------------------------------------------------------------------------------------
int IsomerPlan::getPivotHandle() const {
  return pivot_handle;
}

//-------------------------------------------------------------------------------------------------
int IsomerPlan::getMovingAtomCount() const {
  return moving_atoms.size();
}

//-------------------------------------------------------------------------------------------------
int IsomerPlan::getMovingAtom(const size_t index) const {
  if (index > moving_atoms.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a plan containing " +
          std::to_string(moving_atoms.size()) + " moving atoms.", "IsomerPlan", "getMovingAtom");
  }
  return moving_atoms[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& IsomerPlan::getMovingAtoms() const {
  return moving_atoms;
}

//-------------------------------------------------------------------------------------------------
void IsomerPlan::addMovingAtoms(const std::vector<int> &new_atom_idx) {
  moving_atoms.insert(moving_atoms.end(), new_atom_idx.begin(), new_atom_idx.end());
}

//-------------------------------------------------------------------------------------------------
void IsomerPlan::eraseMovingAtoms() {
  moving_atoms.resize(0);
}

//-------------------------------------------------------------------------------------------------
CoupledEdit::CoupledEdit() :
    edit{ConformationEdit::BOND_ROTATION}, index{0}
{}

//-------------------------------------------------------------------------------------------------
CoupledEdit::CoupledEdit(const ConformationEdit edit_in, const int index_in) :
    edit{edit_in}, index{index_in}
{}

//-------------------------------------------------------------------------------------------------
CoupledEdit::CoupledEdit(const int2 data_in) :
    edit{static_cast<ConformationEdit>(data_in.x)}, index{data_in.y}
{}

//-------------------------------------------------------------------------------------------------
CoupledEdit CoupledEdit::operator=(const int2 &other) {
  CoupledEdit result(other);
  return result;
}

//-------------------------------------------------------------------------------------------------
ChemicalFeaturesReader::ChemicalFeaturesReader(const int atom_count_in,
                                               const int planar_atom_count_in,
                                               const int ring_count_in,
                                               const int fused_ring_count_in,
                                               const int twistable_ring_count_in,
                                               const int conjugated_group_count_in,
                                               const int aromatic_group_count_in,
                                               const int polar_hydrogen_count_in,
                                               const int hbond_donor_count_in,
                                               const int hbond_acceptor_count_in,
                                               const int chiral_center_count_in,
                                               const int rotatable_bond_count_in,
                                               const int cis_trans_bond_count_in,
                                               const int double_bond_count_in,
                                               const int triple_bond_count_in,
                                               const int max_ring_size_in,
                                               const double temperature_in,
                                               const bool rotating_groups_mapped_in,
                                               const bool chiralities_computed_in,
                                               const int* planar_centers_in,
                                               const int* ring_atoms_in,
                                               const int* ring_atom_bounds_in,
                                               const int* aromatic_pi_electrons_in,
                                               const int* aromatic_groups_in,
                                               const int* aromatic_group_bounds_in,
                                               const int* polar_hydrogens_in,
                                               const int* hydrogen_bond_donors_in,
                                               const int* hydrogen_bond_acceptors_in,
                                               const int* chiral_centers_in,
                                               const int* chiral_inversion_methods_in,
                                               const int* rotatable_groups_in,
                                               const int* rotatable_group_bounds_in,
                                               const int* cis_trans_groups_in,
                                               const int* cis_trans_group_bounds_in,
                                               const int* invertible_groups_in,
                                               const int* invertible_group_bounds_in,
                                               const int4* chiral_arm_atoms_in,
                                               const double* formal_charges_in,
                                               const double* bond_orders_in,
                                               const double* free_electrons_in,
                                               const ullint* ring_inclusion_in,
                                               const AtomGraph *ag_pointer_in) :
    atom_count{atom_count_in}, planar_atom_count{planar_atom_count_in}, ring_count{ring_count_in},
    fused_ring_count{fused_ring_count_in}, twistable_ring_count{twistable_ring_count_in},
    conjugated_group_count{conjugated_group_count_in},
    aromatic_group_count{aromatic_group_count_in}, polar_hydrogen_count{polar_hydrogen_count_in},
    hbond_donor_count{hbond_donor_count_in}, hbond_acceptor_count{hbond_acceptor_count_in},
    chiral_center_count{chiral_center_count_in}, rotatable_bond_count{rotatable_bond_count_in},
    cis_trans_bond_count{cis_trans_bond_count_in}, double_bond_count{double_bond_count_in},
    triple_bond_count{triple_bond_count_in}, max_ring_size{max_ring_size_in},
    temperature{temperature_in}, rotating_groups_mapped{rotating_groups_mapped_in},
    chiralities_computed{chiralities_computed_in}, planar_centers{planar_centers_in},
    ring_atoms{ring_atoms_in}, ring_atom_bounds{ring_atom_bounds_in},
    aromatic_pi_electrons{aromatic_pi_electrons_in}, aromatic_groups{aromatic_groups_in},
    aromatic_group_bounds{aromatic_group_bounds_in}, polar_hydrogens{polar_hydrogens_in},
    hydrogen_bond_donors{hydrogen_bond_donors_in},
    hydrogen_bond_acceptors{hydrogen_bond_acceptors_in}, chiral_centers{chiral_centers_in},
    chiral_inversion_methods{chiral_inversion_methods_in}, rotatable_groups{rotatable_groups_in},
    rotatable_group_bounds{rotatable_group_bounds_in}, cis_trans_groups{cis_trans_groups_in},
    cis_trans_group_bounds{cis_trans_group_bounds_in}, invertible_groups{invertible_groups_in},
    invertible_group_bounds{invertible_group_bounds_in}, chiral_arm_atoms{chiral_arm_atoms_in},
    formal_charges{formal_charges_in}, bond_orders{bond_orders_in},
    free_electrons{free_electrons_in}, ring_inclusion{ring_inclusion_in}, ag_pointer{ag_pointer_in}
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    atom_count{(ag_in == nullptr) ? 0 : ag_in->getAtomCount()}, planar_atom_count{0},
    ring_count{0}, fused_ring_count{0}, twistable_ring_count{0}, conjugated_group_count{0},
    aromatic_group_count{0}, polar_hydrogen_count{0}, hbond_donor_count{0},
    hbond_acceptor_count{0}, chiral_center_count{0}, rotatable_bond_count{0},
    cis_trans_bond_count{0}, double_bond_count{0}, triple_bond_count{0},
    max_ring_size{8 * sizeof(ullint)}, temperature{temperature_in},
    rotating_groups_mapped{false},
    chiralities_computed{false},
    planar_centers{HybridKind::POINTER, "chemfe_planarity"},
    ring_inclusion{HybridKind::ARRAY, "chemfe_rings"},
    ring_atom_bounds{HybridKind::POINTER, "chemfe_ring_bounds"},
    ring_atoms{HybridKind::POINTER, "chemfe_ring_atoms"},
    aromatic_group_bounds{HybridKind::POINTER, "chemfe_arom_bounds"},
    aromatic_pi_electrons{HybridKind::POINTER, "chemfe_pi_elec"},
    aromatic_groups{HybridKind::POINTER, "chemfe_arom_groups"},
    polar_hydrogens{HybridKind::POINTER, "chemfe_polar_h"},
    hydrogen_bond_donors{HybridKind::POINTER, "chemfe_hbond_donor"},
    hydrogen_bond_acceptors{HybridKind::POINTER, "chemfe_hbond_accpt"},
    chiral_arm_atoms{HybridKind::ARRAY, "chemfe_chiral_arms"},
    chiral_centers{HybridKind::POINTER, "chemfe_chirals"},
    chiral_inversion_methods{HybridKind::POINTER, "chemfe_chiral_meth"},
    rotatable_groups{HybridKind::POINTER, "chemfe_rotators"},
    rotatable_group_bounds{HybridKind::POINTER, "chemfe_rotator_bounds"},
    cis_trans_groups{HybridKind::POINTER, "chemfe_cistrans"},
    cis_trans_group_bounds{HybridKind::POINTER, "chemfe_cistrans_bounds"},
    invertible_groups{HybridKind::POINTER, "chemfe_invertors"},
    invertible_group_bounds{HybridKind::POINTER, "chemfe_invertor_bounds"},
    anchor_a_branches{HybridKind::POINTER, "chemfe_anchor_i"},
    anchor_b_branches{HybridKind::POINTER, "chemfe_anchor_ii"},
    formal_charges{HybridKind::POINTER, "chemfe_formal_charges"},
    bond_orders{HybridKind::POINTER, "chemfe_bond_orders"},
    free_electrons{HybridKind::POINTER, "chemfe_free_e"},
    int_data{HybridKind::ARRAY, "chemfe_int"},
    mutable_group_data{HybridKind::ARRAY, "chemfe_muta_int"},
    double_data{HybridKind::ARRAY, "chemfe_double"},
    zerok_formal_charges{},
    zerok_bond_orders{},
    zerok_free_electrons{},
    bond_in_ring{},
    ag_pointer{const_cast<AtomGraph*>(ag_in)},
    timer{timer_in}
{
  if (ag_in == nullptr) {
    return;
  }
  analyzeTopology(map_groups_in);
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph &ag_in, const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    ChemicalFeatures(ag_in.getSelfPointer(), map_groups_in, temperature_in, timer_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrameReader &cfr,
                                   const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    ChemicalFeatures(ag_in, map_groups_in, temperature_in, timer_in)
{
  // Correct the chiral centers' orientations
  findChiralOrientations(cfr);
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrame &cf,
                                   const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    ChemicalFeatures(ag_in, cf.data(), map_groups_in, temperature_in, timer_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const PhaseSpace &ps,
                                   const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    ChemicalFeatures(ag_in, CoordinateFrameReader(ps), map_groups_in, temperature_in, timer_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph &ag_in, const CoordinateFrame &cf,
                                   const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    ChemicalFeatures(ag_in.getSelfPointer(), cf, map_groups_in, temperature_in, timer_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph &ag_in, const PhaseSpace &ps,
                                   const MapRotatableGroups map_groups_in,
                                   const double temperature_in, StopWatch *timer_in) :
    ChemicalFeatures(ag_in.getSelfPointer(), CoordinateFrameReader(ps), map_groups_in,
                     temperature_in, timer_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const ChemicalFeatures &original) :
    atom_count{original.atom_count},
    planar_atom_count{original.planar_atom_count},
    ring_count{original.ring_count},
    fused_ring_count{original.fused_ring_count},
    twistable_ring_count{original.twistable_ring_count},
    conjugated_group_count{original.conjugated_group_count},
    aromatic_group_count{original.aromatic_group_count},
    polar_hydrogen_count{original.polar_hydrogen_count},
    hbond_donor_count{original.hbond_donor_count},
    hbond_acceptor_count{original.hbond_acceptor_count},
    chiral_center_count{original.chiral_center_count},
    rotatable_bond_count{original.rotatable_bond_count},
    cis_trans_bond_count{original.cis_trans_bond_count},
    double_bond_count{original.double_bond_count},
    triple_bond_count{original.triple_bond_count},
    max_ring_size{original.max_ring_size},
    temperature{original.temperature},
    rotating_groups_mapped{original.rotating_groups_mapped},
    chiralities_computed{original.chiralities_computed},
    planar_centers{original.planar_centers},
    ring_inclusion{original.ring_inclusion},
    ring_atom_bounds{original.ring_atom_bounds},
    ring_atoms{original.ring_atoms},
    aromatic_group_bounds{original.aromatic_group_bounds},
    aromatic_pi_electrons{original.aromatic_pi_electrons},
    aromatic_groups{original.aromatic_groups},
    polar_hydrogens{original.polar_hydrogens},
    hydrogen_bond_donors{original.hydrogen_bond_donors},
    hydrogen_bond_acceptors{original.hydrogen_bond_acceptors},
    chiral_arm_atoms{original.chiral_arm_atoms},
    chiral_centers{original.chiral_centers},
    chiral_inversion_methods{original.chiral_inversion_methods},
    rotatable_groups{original.rotatable_groups},
    rotatable_group_bounds{original.rotatable_group_bounds},
    cis_trans_groups{original.cis_trans_groups},
    cis_trans_group_bounds{original.cis_trans_group_bounds},
    invertible_groups{original.invertible_groups},
    invertible_group_bounds{original.invertible_group_bounds},
    anchor_a_branches{original.anchor_a_branches},
    anchor_b_branches{original.anchor_b_branches},
    formal_charges{original.formal_charges},
    bond_orders{original.bond_orders},
    free_electrons{original.free_electrons},
    int_data{original.int_data},
    mutable_group_data{original.mutable_group_data},
    double_data{original.double_data},
    zerok_formal_charges{original.zerok_formal_charges},
    zerok_bond_orders{original.zerok_bond_orders},
    zerok_free_electrons{original.zerok_free_electrons},
    bond_in_ring{original.bond_in_ring},
    ag_pointer{original.ag_pointer},
    timer{original.timer}
{
  repairPointers();
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(ChemicalFeatures &&original) :
    atom_count{original.atom_count},
    planar_atom_count{original.planar_atom_count},
    ring_count{original.ring_count},
    fused_ring_count{original.fused_ring_count},
    twistable_ring_count{original.twistable_ring_count},
    conjugated_group_count{original.conjugated_group_count},
    aromatic_group_count{original.aromatic_group_count},
    polar_hydrogen_count{original.polar_hydrogen_count},
    hbond_donor_count{original.hbond_donor_count},
    hbond_acceptor_count{original.hbond_acceptor_count},
    chiral_center_count{original.chiral_center_count},
    rotatable_bond_count{original.rotatable_bond_count},
    cis_trans_bond_count{original.cis_trans_bond_count},
    double_bond_count{original.double_bond_count},
    triple_bond_count{original.triple_bond_count},
    max_ring_size{original.max_ring_size},
    temperature{original.temperature},
    rotating_groups_mapped{original.rotating_groups_mapped},
    chiralities_computed{original.chiralities_computed},
    planar_centers{std::move(original.planar_centers)},
    ring_inclusion{std::move(original.ring_inclusion)},
    ring_atom_bounds{std::move(original.ring_atom_bounds)},
    ring_atoms{std::move(original.ring_atoms)},
    aromatic_group_bounds{std::move(original.aromatic_group_bounds)},
    aromatic_pi_electrons{std::move(original.aromatic_pi_electrons)},
    aromatic_groups{std::move(original.aromatic_groups)},
    polar_hydrogens{std::move(original.polar_hydrogens)},
    hydrogen_bond_donors{std::move(original.hydrogen_bond_donors)},
    hydrogen_bond_acceptors{std::move(original.hydrogen_bond_acceptors)},
    chiral_arm_atoms{std::move(original.chiral_arm_atoms)},
    chiral_centers{std::move(original.chiral_centers)},
    chiral_inversion_methods{std::move(original.chiral_inversion_methods)},
    rotatable_groups{std::move(original.rotatable_groups)},
    rotatable_group_bounds{std::move(original.rotatable_group_bounds)},
    cis_trans_groups{std::move(original.cis_trans_groups)},
    cis_trans_group_bounds{std::move(original.cis_trans_group_bounds)},
    invertible_groups{std::move(original.invertible_groups)},
    invertible_group_bounds{std::move(original.invertible_group_bounds)},
    anchor_a_branches{std::move(original.anchor_a_branches)},
    anchor_b_branches{std::move(original.anchor_b_branches)},
    formal_charges{std::move(original.formal_charges)},
    bond_orders{std::move(original.bond_orders)},
    free_electrons{std::move(original.free_electrons)},
    int_data{std::move(original.int_data)},
    mutable_group_data{std::move(original.mutable_group_data)},
    double_data{std::move(original.double_data)},
    zerok_formal_charges{std::move(original.zerok_formal_charges)},
    zerok_bond_orders{std::move(original.zerok_bond_orders)},
    zerok_free_electrons{std::move(original.zerok_free_electrons)},
    bond_in_ring{std::move(original.bond_in_ring)},
    ag_pointer{original.ag_pointer},
    timer{original.timer}
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures& ChemicalFeatures::operator=(const ChemicalFeatures &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy elements of the other object
  atom_count = other.atom_count;
  planar_atom_count = other.planar_atom_count;
  ring_count = other.ring_count;
  fused_ring_count = other.fused_ring_count;
  twistable_ring_count = other.twistable_ring_count;
  conjugated_group_count = other.conjugated_group_count;
  aromatic_group_count = other.aromatic_group_count;
  polar_hydrogen_count = other.polar_hydrogen_count;
  hbond_donor_count = other.hbond_donor_count;
  hbond_acceptor_count = other.hbond_acceptor_count;
  chiral_center_count = other.chiral_center_count;
  rotatable_bond_count = other.rotatable_bond_count;
  cis_trans_bond_count = other.cis_trans_bond_count;
  double_bond_count = other.double_bond_count;
  triple_bond_count = other.triple_bond_count;
  max_ring_size = other.max_ring_size;
  temperature = other.temperature;
  rotating_groups_mapped = other.rotating_groups_mapped;
  chiralities_computed = other.chiralities_computed;
  planar_centers = other.planar_centers;
  ring_inclusion = other.ring_inclusion;
  ring_atom_bounds = other.ring_atom_bounds;
  ring_atoms = other.ring_atoms;
  aromatic_group_bounds = other.aromatic_group_bounds;
  aromatic_pi_electrons = other.aromatic_pi_electrons;
  aromatic_groups = other.aromatic_groups;
  polar_hydrogens = other.polar_hydrogens;
  hydrogen_bond_donors = other.hydrogen_bond_donors;
  hydrogen_bond_acceptors = other.hydrogen_bond_acceptors;
  chiral_arm_atoms = other.chiral_arm_atoms;
  chiral_centers = other.chiral_centers;
  chiral_inversion_methods = other.chiral_inversion_methods;
  rotatable_groups = other.rotatable_groups;
  rotatable_group_bounds = other.rotatable_group_bounds;
  cis_trans_groups = other.cis_trans_groups;
  cis_trans_group_bounds = other.cis_trans_group_bounds;
  invertible_groups = other.invertible_groups;
  invertible_group_bounds = other.invertible_group_bounds;
  anchor_a_branches = other.anchor_a_branches;
  anchor_b_branches = other.anchor_b_branches;
  formal_charges = other.formal_charges;
  bond_orders = other.bond_orders;
  free_electrons = other.free_electrons;
  int_data = other.int_data;
  mutable_group_data = other.mutable_group_data;
  double_data = other.double_data;
  zerok_formal_charges = other.zerok_formal_charges;
  zerok_bond_orders = other.zerok_bond_orders;
  zerok_free_electrons = other.zerok_free_electrons;
  bond_in_ring = other.bond_in_ring;
  ag_pointer = other.ag_pointer;
  timer = other.timer;
  repairPointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures& ChemicalFeatures::operator=(ChemicalFeatures &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy elements of the other object
  atom_count = other.atom_count;
  planar_atom_count = other.planar_atom_count;
  ring_count = other.ring_count;
  fused_ring_count = other.fused_ring_count;
  twistable_ring_count = other.twistable_ring_count;
  conjugated_group_count = other.conjugated_group_count;
  polar_hydrogen_count = other.polar_hydrogen_count;
  hbond_donor_count = other.hbond_donor_count;
  hbond_acceptor_count = other.hbond_acceptor_count;
  aromatic_group_count = other.aromatic_group_count;
  chiral_center_count = other.chiral_center_count;
  rotatable_bond_count = other.rotatable_bond_count;
  cis_trans_bond_count = other.cis_trans_bond_count;
  double_bond_count = other.double_bond_count;
  triple_bond_count = other.triple_bond_count;
  max_ring_size = other.max_ring_size;
  temperature = other.temperature;
  rotating_groups_mapped = other.rotating_groups_mapped;
  chiralities_computed = other.chiralities_computed;
  planar_centers = std::move(other.planar_centers);
  ring_inclusion = std::move(other.ring_inclusion);
  ring_atom_bounds = std::move(other.ring_atom_bounds);
  ring_atoms = std::move(other.ring_atoms);
  aromatic_group_bounds = std::move(other.aromatic_group_bounds);
  aromatic_pi_electrons = std::move(other.aromatic_pi_electrons);
  aromatic_groups = std::move(other.aromatic_groups);
  polar_hydrogens = std::move(other.polar_hydrogens);
  hydrogen_bond_donors = std::move(other.hydrogen_bond_donors);
  hydrogen_bond_acceptors = std::move(other.hydrogen_bond_acceptors);
  chiral_arm_atoms = std::move(other.chiral_arm_atoms);
  chiral_centers = std::move(other.chiral_centers);
  chiral_inversion_methods = std::move(other.chiral_inversion_methods);
  rotatable_groups = std::move(other.rotatable_groups);
  rotatable_group_bounds = std::move(other.rotatable_group_bounds);
  cis_trans_groups = std::move(other.cis_trans_groups);
  cis_trans_group_bounds = std::move(other.cis_trans_group_bounds);
  invertible_groups = std::move(other.invertible_groups);
  invertible_group_bounds = std::move(other.invertible_group_bounds);
  anchor_a_branches = std::move(other.anchor_a_branches);
  anchor_b_branches = std::move(other.anchor_b_branches);
  formal_charges = std::move(other.formal_charges);
  bond_orders = std::move(other.bond_orders);
  free_electrons = std::move(other.free_electrons);
  int_data = std::move(other.int_data);
  mutable_group_data = std::move(other.mutable_group_data);
  double_data = std::move(other.double_data);
  zerok_formal_charges = std::move(other.zerok_formal_charges);
  zerok_bond_orders = std::move(other.zerok_bond_orders);
  zerok_free_electrons = std::move(other.zerok_free_electrons);
  bond_in_ring = std::move(other.bond_in_ring);
  ag_pointer = other.ag_pointer;
  timer = other.timer;
  return *this;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::analyzeTopology(const MapRotatableGroups map_groups_in) {

  // Set up the StopWatch and obtain indices of the relevant timing bins
  int chemfe_misc_timing, ring_trace_timing, aromatic_group_timing, lewis_timing;
  int chiral_center_timing;
  if (timer != nullptr) {
    timer->assignTime(0);
    chemfe_misc_timing = timer->addCategory("[ChemicalFeatures] Miscellaneous");
    ring_trace_timing = timer->addCategory("[ChemicalFeatures] Ring tracing");
    lewis_timing = timer->addCategory("[ChemicalFeatures] Lewis structure drawing");
    aromatic_group_timing = timer->addCategory("[ChemicalFeatures] Aromaticity detection");
    chiral_center_timing = timer->addCategory("[ChemicalFeatures] Chiral center detection");
  }  
  
  // Obtain abstracts from the topology
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();

  // Detect improper torsion terms and the atoms that make up their centers.  This is a first pass
  // at what is inhibited from rotating.
  std::vector<int> tmp_planar_centers = findPlanarAtoms(vk);
  
  // Detect rings by the bonding pattern
  std::vector<ullint> tmp_ring_inclusion(atom_count, 0LLU);
  std::vector<int> tmp_ring_atom_bounds(1, 0);
  std::vector<int> tmp_ring_atoms;
  if (timer != nullptr) timer->assignTime(chemfe_misc_timing);
  traceTopologicalRings(nbk, cdk, &tmp_ring_inclusion, &tmp_ring_atoms, &tmp_ring_atom_bounds);
  if (timer != nullptr) timer->assignTime(ring_trace_timing);
  ring_inclusion.resize(atom_count);
  ring_inclusion.putHost(tmp_ring_inclusion);
  ring_count = static_cast<int>(tmp_ring_atom_bounds.size()) - 1;

  // Prepare a table of atoms that are part of rings
  bond_in_ring.resize(vk.nbond, false);
  for (int i = 0; i < ring_count; i++) {
    for (int j = tmp_ring_atom_bounds[i]; j < tmp_ring_atom_bounds[i + 1]; j++) {
      const int jatom = tmp_ring_atoms[j];
      for (int k = vk.bond_asgn_bounds[jatom]; k < vk.bond_asgn_bounds[jatom + 1]; k++) {
        const int katom = vk.bond_asgn_atoms[k];
        bool tb_in_ring = false;
        for (int m = tmp_ring_atom_bounds[i]; m < tmp_ring_atom_bounds[i + 1]; m++) {
          tb_in_ring = (tb_in_ring || katom == tmp_ring_atoms[m]);
        }
        bond_in_ring[vk.bond_asgn_terms[k]] = (bond_in_ring[vk.bond_asgn_terms[k]] || tb_in_ring);
      }
    }
  }

  // Allocate the double-precision real data that will result from Lewis structure determination
  const int padded_atom_count = roundUp(atom_count, warp_size_int);
  const int ndbl = 2 * padded_atom_count + roundUp(vk.nbond, warp_size_int);
  double_data.resize(ndbl);
  double_data.shrinkToFit();
  formal_charges.setPointer(&double_data, 0, padded_atom_count);
  free_electrons.setPointer(&double_data, padded_atom_count, padded_atom_count);
  bond_orders.setPointer(&double_data, 2 * padded_atom_count, roundUp(vk.nbond, warp_size_int));
  zerok_formal_charges.resize(padded_atom_count);
  zerok_free_electrons.resize(padded_atom_count);
  zerok_bond_orders.resize(bond_orders.size());

  // Draw a Lewis structure and record the results.  The formal_charges, bond orders, and
  // free_electrons arrays can now be allocated.  Lewis structures will be drawn for each unique
  // molecule and copied otherwise, but that requires a list of all unique molecules.
  if (timer != nullptr) timer->assignTime(chemfe_misc_timing);
  drawLewisStructures(vk, nbk, cdk);
  if (timer != nullptr) timer->assignTime(lewis_timing);
  
  // Mark all aromatic atoms after dissecting numerous details about atoms in rings.
  std::vector<int> tmp_aromatic_group_bounds(1, 0);
  std::vector<int> tmp_aromatic_pi_electrons;
  std::vector<int> tmp_aromatic_groups;
  findAromaticGroups(cdk, vk, tmp_ring_atoms, tmp_ring_atom_bounds, &tmp_aromatic_group_bounds,
                     &tmp_aromatic_pi_electrons, &tmp_aromatic_groups);
  if (timer != nullptr) timer->assignTime(aromatic_group_timing);

  // Mark polar hydrogens and hydrogen bond donors and acceptors
  std::vector<int> tmp_polar_hydrogens;
  std::vector<int> tmp_hydrogen_bond_donors;
  std::vector<int> tmp_hydrogen_bond_acceptors;
  findHydrogenBondElements(nbk, cdk, &tmp_polar_hydrogens, &tmp_hydrogen_bond_donors,
                           &tmp_hydrogen_bond_acceptors);
  polar_hydrogen_count = tmp_polar_hydrogens.size();
  hbond_donor_count = tmp_hydrogen_bond_donors.size();
  hbond_acceptor_count = tmp_hydrogen_bond_acceptors.size();
  if (timer != nullptr) timer->assignTime(chemfe_misc_timing);
  
  // Find chiral centers
  const std::vector<int> tmp_chiral_centers = detailChiralCenters(nbk, vk, cdk);
  chiral_center_count = tmp_chiral_centers.size();
  const std::vector<int>
    tmp_chiral_inversion_methods = findChiralInversionMethods(tmp_chiral_centers, tmp_ring_atoms,
                                                              tmp_ring_atom_bounds);
  if (timer != nullptr) timer->assignTime(chiral_center_timing);

  // Store the integer results
  const size_t nint = roundUp(tmp_planar_centers.size(), warp_size_zu) +
                      roundUp(tmp_ring_atom_bounds.size(), warp_size_zu) +
                      roundUp(tmp_ring_atoms.size(), warp_size_zu) +
                      roundUp(tmp_aromatic_group_bounds.size(), warp_size_zu) +
                      roundUp(tmp_aromatic_pi_electrons.size(), warp_size_zu) +
                      roundUp(tmp_aromatic_groups.size(), warp_size_zu) +
                      (2 * roundUp(static_cast<size_t>(chiral_center_count), warp_size_zu)) +
                      roundUp(tmp_polar_hydrogens.size(), warp_size_zu) +
                      roundUp(tmp_hydrogen_bond_donors.size(), warp_size_zu) +
                      roundUp(tmp_hydrogen_bond_acceptors.size(), warp_size_zu);
  int_data.resize(nint);
  int_data.shrinkToFit();
  size_t ic = planar_centers.putHost(&int_data, tmp_planar_centers, 0, warp_size_zu);
  ic = ring_atom_bounds.putHost(&int_data, tmp_ring_atom_bounds, ic, warp_size_zu);
  ic = ring_atoms.putHost(&int_data, tmp_ring_atoms, ic, warp_size_zu);
  ic = aromatic_group_bounds.putHost(&int_data, tmp_aromatic_group_bounds, ic, warp_size_zu);
  ic = aromatic_pi_electrons.putHost(&int_data, tmp_aromatic_pi_electrons, ic, warp_size_zu);
  ic = aromatic_groups.putHost(&int_data, tmp_aromatic_groups, ic, warp_size_zu);
  ic = polar_hydrogens.putHost(&int_data, tmp_polar_hydrogens, ic, warp_size_zu);
  ic = hydrogen_bond_donors.putHost(&int_data, tmp_hydrogen_bond_donors, ic, warp_size_zu);
  ic = hydrogen_bond_acceptors.putHost(&int_data, tmp_hydrogen_bond_acceptors, ic, warp_size_zu);
  ic = chiral_centers.putHost(&int_data, tmp_chiral_centers, ic, warp_size_zu);
  ic = chiral_inversion_methods.putHost(&int_data, tmp_chiral_inversion_methods, ic, warp_size_zu);
  if (timer != nullptr) timer->assignTime(chemfe_misc_timing);
  
  // Find rotatable bonds, if group mapping is active, and map the invertible groups that will
  // flip the chirality of already detected chiral centers.
  switch (map_groups_in) {
  case MapRotatableGroups::YES:
    findRotatableBondGroups(timer);
    break;
  case MapRotatableGroups::NO:
    allocateMutableData({}, {}, {}, {}, {}, {}, {}, {});
    break;
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::findPlanarAtoms(const ValenceKit<double> &vk) const {
  std::vector<int> result;
  for (int i = 0; i < vk.ndihe; i++) {
    const TorsionKind ikind = static_cast<TorsionKind>(vk.dihe_modifiers[i].w);
    if (ikind == TorsionKind::IMPROPER_NO_14 || ikind == TorsionKind::IMPROPER) {

      // The third atom of an improper is the central atom, about which a plane is enforced
      result.push_back(vk.dihe_k_atoms[i]);
    }
  }
  for (int i = 0; i < vk.ncimp; i++) {
    result.push_back(vk.cimp_k_atoms[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::traceTopologicalRings(const NonbondedKit<double> &nbk,
                                             const ChemicalDetailsKit &cdk,
                                             std::vector<ullint> *tmp_ring_inclusion,
                                             std::vector<int> *tmp_ring_atoms,
                                             std::vector<int> *tmp_ring_atom_bounds) {
  int max_branches = 0;
  for (int i = 0; i < atom_count; i++) {
    max_branches = std::max(nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i], max_branches);
  }
  const int max_molecule = ag_pointer->getLargestMoleculeSize();
  std::vector<int> all_branch_atoms(max_branches * max_molecule);
  std::vector<BondedNode> links(atom_count);
  std::vector<int> tree_positions(atom_count, -1);
  for (int i = 0; i < max_molecule; i++) {
    links[i].setBranchPointer(&all_branch_atoms, i, max_branches);
  }

  // Loop over all atoms, proceeding to explore all available bonds, until all atoms have been
  // either used as the start of chain / ring exploration or have been included in the exploration
  // initiated by some other atom.
  std::vector<bool> valid_atom_mask(atom_count);
  for (int i = 0; i < atom_count; i++) {
    valid_atom_mask[i] = (cdk.z_numbers[i] > 0);
  }
  std::vector<bool> atom_touched(atom_count, false);
  for (int i = 0; i < atom_count; i++) {
    if (cdk.z_numbers[i] == 0) {
      atom_touched[i] = true;
    }
    if (atom_touched[i]) {
      continue;
    }
    
    // Initiate the chain
    atom_touched[i] = true;
    links[0].addToTree(-1, i, 0, nbk, valid_atom_mask);
    tree_positions[i] = 0;
    int node_count = 1;
    int current_layer = 1;
    int layer_llim = 0;
    int layer_hlim = 1;
    while (layer_hlim > layer_llim) {
      const int next_layer_llim = node_count;
      for (int j = layer_llim; j < layer_hlim; j++) {
        const int j_atom = links[j].getAtom();
        const int j_branch_count = links[j].getBranchCount();
        for (int k = 0; k < j_branch_count; k++) {
          const int k_atom = links[j].getBranchAtom(k);
          
          // Check the status of the next atom: has it already been incorporated into this chain
          // or an earlier one?  If so, this signifies the completion of some loop.  Determine that
          // loop based on the histories of the current atom and the atom it touches.  Otherwise,
          // add the next atom to the chain.
          if (atom_touched[k_atom] && cdk.z_numbers[k_atom] != 0) {
            markRingAtoms(j_atom, k_atom, tree_positions, node_count, &links, tmp_ring_inclusion,
                          tmp_ring_atoms, tmp_ring_atom_bounds, cdk);
          }
          else {
            atom_touched[k_atom] = true;
            links[node_count].addToTree(j_atom, k_atom, current_layer, nbk, valid_atom_mask);
            tree_positions[k_atom] = node_count;
            node_count++;
          }
        }
      }
      layer_llim = next_layer_llim;
      layer_hlim = node_count;
      current_layer++;
    }
    
    // Wipe any rings that were found, to allow future trees to evaluate new rings
    for (int j = 0; j < layer_hlim; j++) {
      links[j].wipeRingCompletion();
    }
  }
  
  // Fused rings will imply additional cyclic structures than have already been detected.
  // Further analyze any rings with overlapping sets of atoms.
  std::vector<int> ring_participation_bounds(atom_count + 1, 0);
  int* ring_atoms_ptr  = tmp_ring_atoms->data();
  int* ring_bounds_ptr = tmp_ring_atom_bounds->data();
  const int nring = tmp_ring_atom_bounds->size() - 1;
  const int nring_atoms = ring_bounds_ptr[nring];
  for (int i = 0; i < nring; i++) {
    for (int j = ring_bounds_ptr[i]; j < ring_bounds_ptr[i + 1]; j++) {
      ring_participation_bounds[ring_atoms_ptr[j]] += 1;
    }
  }
  prefixSumInPlace(&ring_participation_bounds, PrefixSumType::EXCLUSIVE, "ChemicalFeatures");
  std::vector<int> ring_participation(ring_participation_bounds[atom_count]);
  for (int i = 0; i < nring; i++) {
    for (int j = ring_bounds_ptr[i]; j < ring_bounds_ptr[i + 1]; j++) {
      const size_t atom_idx = ring_atoms_ptr[j];
      const int partic_idx = ring_participation_bounds[atom_idx];
      ring_participation[partic_idx] = i;
      ring_participation_bounds[atom_idx] = partic_idx + 1;
    }
  }
  for (int i = atom_count; i > 0; i--) {
    ring_participation_bounds[i] = ring_participation_bounds[i - 1];    
  }
  ring_participation_bounds[0] = 0;

  // Search all bonds coming out of atoms in fused rings: do they imply other rings that have
  // not yet been detected?
  std::vector<int> fused_rings(16);
  std::vector<int> fused_ring_atoms(32);
  std::vector<int> all_fused_rings;
  std::vector<int> fused_ring_atom_connections;
  std::vector<int> ranked_connections;
  std::vector<bool> fusion_coverage(atom_count, false);
  std::vector<int> layer_bounds(32);
  std::vector<int> jk_paths(32);
  std::vector<int> jk_path_bounds(32);
  std::vector<int> subsystem_ring_atoms;
  std::vector<int> subsystem_ring_atom_bounds(1, 0);
  std::vector<int2> subsystem_ring_ranks;
  std::vector<int> treaded_atoms(atom_count, 0);
  for (int i = 0; i < atom_count; i++) {
    valid_atom_mask[i] = false;
  }
  for (int i = 0; i < atom_count; i++) {
    if (fusion_coverage[i]) {
      continue;
    }
    if (ring_participation_bounds[i + 1] - ring_participation_bounds[i] <= 1) {
      fusion_coverage[i] = true;
      continue;
    }
    
    // This atom participates in two or more rings and has not yet been analyzed, implying an
    // unattended multi-ring system.  Map the extent of that system.
    fused_rings.resize(0);
    for (int j = ring_participation_bounds[i]; j < ring_participation_bounds[i + 1]; j++) {
      fused_rings.push_back(ring_participation[j]);
    }
    int prev_ring_count = 0;
    int curr_ring_count = fused_rings.size();
    while (curr_ring_count > prev_ring_count) {

      // Loop over all newly added rings and see if the atoms they contain point to other rings.
      for (int pos = prev_ring_count; pos < curr_ring_count; pos++) {
        const int pring = fused_rings[pos];        
        for (int j = ring_bounds_ptr[pring]; j < ring_bounds_ptr[pring + 1]; j++) {
          const int atom_idx = ring_atoms_ptr[j];
          for (int k = ring_participation_bounds[atom_idx];
               k < ring_participation_bounds[atom_idx + 1]; k++) {
            if (ring_participation[k] != pring) {
              fused_rings.push_back(ring_participation[k]);
            }
          }
        }
      }

      // Cull duplicate entries in the list of fused rings.
      reduceUniqueValues(&fused_rings);
      prev_ring_count = curr_ring_count;
      curr_ring_count = fused_rings.size();
    }
    
    // Mark all atoms in the complete fused ring system to avoid investigating them again.  Take
    // this opportunity to make a list of rings already present in the fused system, to accumulate
    // secondary ring findings.  Store the indices of unique fused rings in this subsystem for
    // culling later.
    const int nfused_rings = fused_rings.size();
    int s2_est_size = 0;
    for (int pos = 0; pos < nfused_rings; pos++) {
      const int pring = fused_rings[pos];
      all_fused_rings.push_back(fused_rings[pos]);
      s2_est_size += ring_bounds_ptr[pring + 1] - ring_bounds_ptr[pring];
    }
    std::vector<int> s2_ring_atoms(s2_est_size * 2);
    std::vector<int> s2_ring_bounds(nfused_rings * 2);
    s2_ring_atoms.resize(0);
    s2_ring_bounds.resize(0);
    fused_ring_atoms.resize(0);
    for (int pos = 0; pos < nfused_rings; pos++) {
      const int pring = fused_rings[pos];
      for (int j = ring_bounds_ptr[pring]; j < ring_bounds_ptr[pring + 1]; j++) {
        fused_ring_atoms.push_back(ring_atoms_ptr[j]);
        s2_ring_atoms.push_back(ring_atoms_ptr[j]);
      }
      s2_ring_bounds[pos + 1] = s2_ring_atoms.size();
    }
    reduceUniqueValues(&fused_ring_atoms);
    const int nfused_atoms = fused_ring_atoms.size();
    for (int j = 0; j < nfused_atoms; j++) {
      fusion_coverage[fused_ring_atoms[j]] = true;
      valid_atom_mask[fused_ring_atoms[j]] = true;
    }
    
    // Search the 1:2 exclusions of each atom in the fused ring system to see if it connects to
    // other atoms in the fused ring system.  Find all such atoms that connect to three or more
    // other atoms.  These are the atoms in the fused ring system that will form smaller rings,
    // and must be searched in O(N^2) fashion for the smallest rings that encompass both atoms.
    fused_ring_atom_connections.resize(nfused_atoms);
    ranked_connections.resize(nfused_atoms);
    for (int pos = 0; pos < nfused_atoms; pos++) {
      fused_ring_atom_connections[pos] = 0;
      const int atom_idx = fused_ring_atoms[pos];
      for (int j = nbk.nb12_bounds[atom_idx]; j < nbk.nb12_bounds[atom_idx + 1]; j++) {
        if (valid_atom_mask[nbk.nb12x[j]] && cdk.z_numbers[nbk.nb12x[j]] > 0) {
          fused_ring_atom_connections[pos] += 1;
        }
      }
      ranked_connections[pos] = fused_ring_atom_connections[pos];
    }
    std::sort(ranked_connections.begin(), ranked_connections.end(),
              [](int a, int b) { return a > b; });
    addScalarToVector(&ranked_connections, -1);
    
    // Calculate the maximum size of a tree for this fused ring system.  Keep the search bounded
    // by searching only for secondary rings of up to a certain size starting at each atom j.
    const int max_splits = ranked_connections[0];
    int max_tree_depth = 1;
    int prospective_tree_size = max_splits + 1;
    int layer_size = max_splits;
    int ranked_pos = 1;
    while (prospective_tree_size < max_fused_ring_tree_size && max_tree_depth < nfused_atoms - 1) {
      layer_size *= ranked_connections[ranked_pos];
      prospective_tree_size += layer_size;
      ranked_pos++;
      max_tree_depth++;
    }
    max_tree_depth -= (prospective_tree_size > max_fused_ring_tree_size);
    
    // Reallocate the tree search arrays if more memory is needed
    if (static_cast<size_t>(prospective_tree_size) > tree_positions.size()) {
      tree_positions.resize(prospective_tree_size);
    }
    if (static_cast<size_t>(prospective_tree_size) > links.size()) {
      links.resize(prospective_tree_size);
    }
    if (static_cast<size_t>(prospective_tree_size * max_branches) > all_branch_atoms.size()) {
      all_branch_atoms.resize(prospective_tree_size * max_branches);
      for (int i = 0; i < prospective_tree_size; i++) {
        links[i].setBranchPointer(&all_branch_atoms, i, max_branches);
      }
    }
    int n_multi_split = 0;
    for (int j = 1; j < nfused_atoms; j++) {
      if (fused_ring_atom_connections[j] == 2) {
        continue;
      }

      // Increment the count of nodes that might branch into new rings.  If only one has thus
      // far been encountered, keep going, as it requires two such nodes to complete a ring.
      n_multi_split += fused_ring_atom_connections[j] - 2;
      if (n_multi_split == 1) {
        continue;
      }
      const int jatom = fused_ring_atoms[j];
      
      // Make a tree, starting at jatom and extending through all available combinations to
      // arrive back at jatom.
      links[0].addToTree(-1, jatom, 0, -1, nbk, valid_atom_mask);
      int node_count = 1;
      int layer_count = 1;
      int layer_llim = 0;
      int layer_hlim = 1;
      layer_bounds.resize(1);
      layer_bounds[0] = 0;
      layer_bounds.push_back(node_count);
      while (layer_hlim > layer_llim && layer_count < max_tree_depth) {
        const int next_layer_llim = node_count;
        for (int pos = layer_llim; pos < layer_hlim; pos++) {
          const int pos_atom = links[pos].getAtom();
          const int pos_branch_count = links[pos].getBranchCount();
          for (int pos2 = 0; pos2 < pos_branch_count; pos2++) {
            const int pos2_atom = links[pos].getBranchAtom(pos2);

            // Check this atom against the history in the tree.  If it has already been visited
            // by this branch of the tree, do not incorporate it into the next layer.
            bool found = false;
            int prev_node = links[pos].getPreviousNode();
            while ((! found) && prev_node != -1) {
              const int prev_atom = links[prev_node].getAtom();
              prev_node = links[prev_node].getPreviousNode();
              found = (found || (prev_atom == pos2_atom &&
                                 (fused_ring_atom_connections[j] < 4 ||
                                  prev_node != -1)));
            }
            if (! found) {
              links[node_count].addToTree(pos_atom, pos2_atom, layer_count, pos, nbk,
                                          valid_atom_mask);
              node_count++;
            }
          }
        }
        layer_llim = next_layer_llim;
        layer_hlim = node_count;
        layer_bounds.push_back(node_count);
        layer_count++;
      }
      const int klim = j + (fused_ring_atom_connections[j] > 3);
      for (int k = 0; k < klim; k++) {
        if (fused_ring_atom_connections[k] == 2) {
          continue;
        }
        const int katom = fused_ring_atoms[k];
        const int path_quota = std::min(fused_ring_atom_connections[j],
                                        fused_ring_atom_connections[k]);

        // Find the shortest and second-shortest paths between atoms jatom and katom, staying
        // within the fused ring atom system.  This occurs when a tree starting at jatom includes
        // katom twice.  If katom is included more than twice in any layer, this implies additional
        // rings of equal size.  Log them all.
        int njk_paths = 0;
        jk_path_bounds.resize(1);
        jk_path_bounds[0] = 0;
        jk_paths.resize(0);
        std::vector<int> current_path;
        current_path.reserve(32);
        for (int lcon = 1; lcon < layer_count; lcon++) {
          for (int node_idx = layer_bounds[lcon]; node_idx < layer_bounds[lcon + 1]; node_idx++) {
            if (links[node_idx].getAtom() == katom) {
              current_path.resize(0);
              current_path.push_back(links[node_idx].getAtom());
              int prev_node = links[node_idx].getPreviousNode();

              // Record the new path, all the way back to jatom
              while (prev_node >= 0) {
                const int atom_on_path = links[prev_node].getAtom();
                current_path.push_back(atom_on_path);
                prev_node = links[prev_node].getPreviousNode();
              }

              // Check that the path is unique, by testing previously found paths in both
              // directions.  If so, add it to the growing list of new paths.
              bool path_independent = true;
              const int current_path_length = current_path.size();
              for (int m = 0; m < njk_paths; m++) {
                if (jk_path_bounds[m + 1] - jk_path_bounds[m] == current_path_length) {
                  int forward_tracking = 0;
                  size_t cpos = 0;
                  const size_t nlim = jk_path_bounds[m + 1];
                  for (size_t n = jk_path_bounds[m]; n < nlim; n++) {
                    forward_tracking += (jk_paths[n] == current_path[cpos]);
                    cpos++;
                  }
                  int backward_tracking = 0;
                  cpos = current_path_length - 1;
                  for (size_t n = jk_path_bounds[m]; n < nlim; n++) {
                    backward_tracking += (jk_paths[n] == current_path[cpos]);
                    cpos--;
                  }
                  if (forward_tracking == current_path_length ||
                      backward_tracking == current_path_length) {
                    path_independent = false;
                  }
                }
              }
              if (path_independent) {
                jk_paths.insert(jk_paths.end(), current_path.begin(), current_path.end());
                jk_path_bounds.push_back(jk_paths.size());
                njk_paths++;
              }
            }
          }
        }

        // With the possible paths computed, list out the rings for this subsystem.
        if (jatom == katom) {
          for (int m = 0; m < njk_paths; m++) {

            // Incorporate the first path, minus the final atom to avoid duplication
            const int mpos_lim = jk_path_bounds[m + 1] - 1;
            for (int mpos = jk_path_bounds[m]; mpos < mpos_lim; mpos++) {
              subsystem_ring_atoms.push_back(jk_paths[mpos]);
            }

            // Update the bounds
            const int2 tmp_rank = { mpos_lim - jk_path_bounds[m],
                                    static_cast<int>(subsystem_ring_atom_bounds.size()) - 1 };
            subsystem_ring_atom_bounds.push_back(subsystem_ring_atoms.size());
            subsystem_ring_ranks.push_back(tmp_rank);
          }
        }
        else {
          for (int m = 1; m < njk_paths; m++) {
            for (int n = 0; n < m; n++) {

              // Incorporate the first path, in its entirety to include both jatom and katom
              for (int npos = jk_path_bounds[n]; npos < jk_path_bounds[n + 1]; npos++) {
                subsystem_ring_atoms.push_back(jk_paths[npos]);
              }
            
              // Incorporate the second path, in reverse, omitting the endpoints, to complete the
              // ring in a sensible order (all paths proceed from katom to jatom, and include both
              // katom and jatom).
              for (int mpos = jk_path_bounds[m + 1] - 2; mpos > jk_path_bounds[m]; mpos--) {
                subsystem_ring_atoms.push_back(jk_paths[mpos]);
              }

              // Update the bounds
              const int2 tmp_rank = { jk_path_bounds[n + 1] - jk_path_bounds[n] +
                                      jk_path_bounds[m + 1] - 2 - jk_path_bounds[m],
                                      static_cast<int>(subsystem_ring_atom_bounds.size()) - 1 };
              subsystem_ring_atom_bounds.push_back(subsystem_ring_atoms.size());
              subsystem_ring_ranks.push_back(tmp_rank);
            }
          }
        }
      }
    }

    // Clear the valid atom mask
    for (int j = 0; j < nfused_atoms; j++) {
      valid_atom_mask[fused_ring_atoms[j]] = false;
    }
  }

  // Return immediately if there were no fused rings
  if (all_fused_rings.size() == 0LLU) {
    return;
  }

  // Cull larger rings if they can be found to be a superset of two or more smaller rings.  This
  // is done by starting with the smallest rings in the set and marking the atoms of the fused
  // ring subsystem that they cover.  As larger and larger rings are added, those that do not
  // cover any new atoms can be deemed irrelevant.  Use the valid_atom_mask array a second time,
  // marking it TRUE for each atom touched by one of the discovered rings.
  std::sort(subsystem_ring_ranks.begin(), subsystem_ring_ranks.end(),
            [](int2 a, int2 b) { return a.x < b.x; });
  const int nsubsys_rings = subsystem_ring_atom_bounds.size() - 1LLU;
  ullint* ring_inclusion_ptr = tmp_ring_inclusion->data();
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  std::vector<bool> valid_bond_mask(vk.nbond, false);
  for (int i = 0; i < nsubsys_rings; i++) {
    const int ring_idx = subsystem_ring_ranks[i].y;
    bool ring_is_essential = false;
    for (int j = subsystem_ring_atom_bounds[ring_idx];
         j < subsystem_ring_atom_bounds[ring_idx + 1]; j++) {

      // Mark the atoms covered by the ring--any ring that covers a new atom is essential.
      const int atom_idx = subsystem_ring_atoms[j];
      ring_is_essential = (ring_is_essential || valid_atom_mask[atom_idx] == false);
      valid_atom_mask[atom_idx] = true;

      // Mark the bonds covered by the ring--any ring that covers a new bonds is essential.
      const int other_idx = (j == subsystem_ring_atom_bounds[ring_idx + 1] - 1) ?
                            subsystem_ring_atoms[subsystem_ring_atom_bounds[ring_idx]] :
                            subsystem_ring_atoms[j + 1];
      for (int k = vk.bond_asgn_bounds[atom_idx]; k < vk.bond_asgn_bounds[atom_idx + 1]; k++) {
        if (vk.bond_asgn_atoms[k] == other_idx) {
          ring_is_essential = (ring_is_essential ||
                               valid_bond_mask[vk.bond_asgn_terms[k]] == false);
          valid_bond_mask[vk.bond_asgn_terms[k]] = true;
        }
      }
      for (int k = vk.bond_asgn_bounds[other_idx]; k < vk.bond_asgn_bounds[other_idx + 1]; k++) {
        if (vk.bond_asgn_atoms[k] == atom_idx) {
          ring_is_essential = (ring_is_essential ||
                               valid_bond_mask[vk.bond_asgn_terms[k]] == false);
          valid_bond_mask[vk.bond_asgn_terms[k]] = true;
        }
      }
    }
    if (ring_is_essential) {
      const int tring_size = subsystem_ring_atom_bounds[ring_idx + 1] -
                             subsystem_ring_atom_bounds[ring_idx];
      if (tring_size < max_ring_size) {
        const ullint incl_mask = (0x1LLU << tring_size);
        for (int j = subsystem_ring_atom_bounds[ring_idx];
             j < subsystem_ring_atom_bounds[ring_idx + 1]; j++) {
          ring_inclusion_ptr[subsystem_ring_atoms[j]] |= incl_mask;
        }
      }
      for (int j = subsystem_ring_atom_bounds[ring_idx];
         j < subsystem_ring_atom_bounds[ring_idx + 1]; j++) {
        tmp_ring_atoms->push_back(subsystem_ring_atoms[j]);
      }
      tmp_ring_atom_bounds->push_back(tmp_ring_atoms->size());
    }
  }

  // Cull the original fused ring systems--they are no longer necessary.  Reset the ring_bounds_ptr
  // and ring_atoms_ptr pointers, as the target arrays have been resized.
  const int tmp_overall_ring_count = static_cast<int>(tmp_ring_atom_bounds->size()) - 1;
  const int original_fused_ring_count = all_fused_rings.size();
  std::vector<bool> cull_rings(tmp_overall_ring_count, false);
  for (int i = 0; i < original_fused_ring_count; i++) {
    cull_rings[all_fused_rings[i]] = true;
  }
  int ratoms_counter = 0;
  int rbounds_counter = 1;
  ring_bounds_ptr = tmp_ring_atom_bounds->data();
  ring_atoms_ptr  = tmp_ring_atoms->data();
  for (int i = 0; i < tmp_overall_ring_count; i++) {
    if (cull_rings[i] == false) {      
      for (int j = ring_bounds_ptr[i]; j < ring_bounds_ptr[i + 1]; j++) {
        ring_atoms_ptr[ratoms_counter] = ring_atoms_ptr[j];
        ratoms_counter++;
      }
      ring_bounds_ptr[rbounds_counter] = ratoms_counter;
      rbounds_counter++;
    }
  }
  tmp_ring_atoms->resize(ratoms_counter);
  tmp_ring_atoms->shrink_to_fit();
  tmp_ring_atom_bounds->resize(rbounds_counter);
  tmp_ring_atom_bounds->shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::markRingAtoms(const int j_atom, const int k_atom,
                                     const std::vector<int> &tree_positions, const int node_count,
                                     std::vector<BondedNode> *links,
                                     std::vector<ullint> *tmp_ring_inclusion,
                                     std::vector<int> *tmp_ring_atoms,
                                     std::vector<int> *tmp_ring_atom_bounds,
                                     const ChemicalDetailsKit &cdk) {

  // Check that the ring has not already been marked
  BondedNode* linkdata = links->data();
  const int jpos = tree_positions[j_atom];
  const int kpos = tree_positions[k_atom];
  const int kbranch_in_j_atom = linkdata[jpos].findBranchAtom(k_atom);
  const int jbranch_in_k_atom = linkdata[kpos].findBranchAtom(j_atom);
  const bool j_marked = (linkdata[jpos].getRingCompletion(kbranch_in_j_atom));
  const bool k_marked = (linkdata[kpos].getRingCompletion(jbranch_in_k_atom));
  if (j_marked != k_marked) {
    rtErr("Inconsistency detected in ring formation for atoms " + std::to_string(j_atom + 1) +
          " and " + std::to_string(k_atom + 1) + ".", "ChemicalFeatures", "markRingAtoms");
  }
  if (j_marked) {
    return;
  }

  // Search backwards from each atom to find the point at which they reach a common atom
  std::vector<int> j_history, k_history;
  j_history.push_back(j_atom);
  k_history.push_back(k_atom);
  int j_length = 1;
  int k_length = 1;
  bool common_point_found = false;
  bool jgrow = true;
  bool kgrow = true;
  int jtrack = jpos;
  int ktrack = kpos;
  int k_pivot, j_pivot;
  while (common_point_found == false) {
    const int jprev = linkdata[jtrack].getPreviousAtom();
    if (jprev >= 0) {
      jtrack = tree_positions[jprev];
    }
    jgrow = (jgrow && (jprev >= 0));
    if (jgrow) {
      j_history.push_back(jprev);
      j_length++;
      for (int k = 0; k < k_length; k++) {
        if (jprev == k_history[k]) {
          common_point_found = true;
          j_pivot = j_length;
          k_pivot = k + 1;
          break;
        }
      }
    }
    if (common_point_found == false) {
      const int kprev = linkdata[ktrack].getPreviousAtom();
      if (kprev >= 0) {
        ktrack = tree_positions[kprev];
      }
      kgrow = (kgrow && (kprev >= 0));
      if (kgrow) {
        k_history.push_back(kprev);
        k_length++;
        for (int j = 0; j < j_length; j++) {
          if (kprev == j_history[j]) {
            common_point_found = true;
            j_pivot = j + 1;
            k_pivot = k_length;
            break;
          }
        }
      }
    }
    if (jgrow == false && kgrow == false) {
      common_point_found = true;
    }
  }
  
  // Set the appropriate bit on the various ring atoms.  Store the atoms of this ring in a
  // growing list.  This needs to be done with push_back as the ring computation is somewhat
  // involved and would be inefficient to repeat.
  const int ring_size = j_pivot + k_pivot - 1;
  if (ring_size < max_ring_size && ring_size >= 3) {
    
    // Step through the J history until the pivot point, then backwards through the K history
    // starting at the pivot point.  This will list all atoms in the ring.    
    int n_planar = 0;
    const ullint ring_mask = (0x1LLU << ring_size);
    ullint* ring_ptr = tmp_ring_inclusion->data();
    for (int j = 0; j < j_pivot; j++) {
      ring_ptr[j_history[j]] |= ring_mask;
      tmp_ring_atoms->push_back(j_history[j]);
    }
    for (int k = k_pivot - 2; k >= 0; k--) {
      ring_ptr[k_history[k]] |= ring_mask;
      tmp_ring_atoms->push_back(k_history[k]);
    }
    tmp_ring_atom_bounds->push_back(tmp_ring_atoms->size());
  }
  
  // Mark the ring completion in each node
  linkdata[jpos].setRingCompletion(kbranch_in_j_atom);
  linkdata[kpos].setRingCompletion(jbranch_in_k_atom);
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::drawLewisStructures(const ValenceKit<double> &vk,
                                           const NonbondedKit<double> &nbk,
                                           const ChemicalDetailsKit &cdk) {

  // Loop over all molecules and find those which are similar
  const int3 init_midx = {-1, 0, 0};
  std::vector<int3> mol_index(cdk.nmol, init_midx);
  int nmi = 0;
  for (int i = 0; i < cdk.nmol; i++) {    
    if (mol_index[i].x >= 0) {
      continue;
    }
    mol_index[i].x = nmi;
    mol_index[i].y = 0;
    mol_index[i].z = 0;
    const int i_llim = cdk.mol_limits[i];
    const int i_hlim = cdk.mol_limits[i + 1];
    const int ni_atom = i_hlim - i_llim;
    for (int j = i + 1; j < cdk.nmol; j++) {
      const int j_llim = cdk.mol_limits[j];
      const int j_hlim = cdk.mol_limits[j + 1];
      if (j_hlim - j_llim != ni_atom) {
        continue;
      }

      // Check the Z numbers
      bool znums_differ = false;
      for (int k = 0; k < ni_atom; k++) {
        const int moli_atom = cdk.mol_contents[i_llim + k];
        const int molj_atom = cdk.mol_contents[j_llim + k];
        znums_differ = (znums_differ || cdk.z_numbers[moli_atom] != cdk.z_numbers[molj_atom]);
      }
      if (znums_differ) {
        continue;
      }
      
      // Check the charges
      bool charges_differ = false;
      for (int k = 0; k < ni_atom; k++) {
        const int moli_atom = cdk.mol_contents[i_llim + k];
        const int molj_atom = cdk.mol_contents[j_llim + k];
        charges_differ = (charges_differ ||
                          fabs(nbk.charge[moli_atom] - nbk.charge[molj_atom]) > 1.0e-4);
      }
      if (charges_differ) {
        continue;
      }

      // Check the bonds
      bool bonds_similar = true;
      for (int k = 0; k < ni_atom; k++) {
        const int moli_atom = cdk.mol_contents[i_llim + k];
        const int molj_atom = cdk.mol_contents[j_llim + k];
        bonds_similar = (bonds_similar &&
                        (vk.bond_asgn_bounds[moli_atom + 1] - vk.bond_asgn_bounds[moli_atom] == 
                         vk.bond_asgn_bounds[molj_atom + 1] - vk.bond_asgn_bounds[molj_atom]));
        if (bonds_similar) {
          const int offset = vk.bond_asgn_atoms[vk.bond_asgn_bounds[molj_atom]] -
                             vk.bond_asgn_atoms[vk.bond_asgn_bounds[moli_atom]];
          int jm = vk.bond_asgn_bounds[molj_atom];
          for (int im = vk.bond_asgn_bounds[moli_atom]; im < vk.bond_asgn_bounds[moli_atom + 1];
               im++) {
              bonds_similar = (bonds_similar &&
                               vk.bond_asgn_atoms[jm] - vk.bond_asgn_atoms[im] == offset);
            jm++;
          }
        }
      }
      if (bonds_similar == false) {
        continue;
      }
      
      // The molecules have been found to be the same.
      mol_index[j].x = nmi;
      mol_index[j].y = j_llim - i_llim;
      mol_index[j].z = vk.bond_asgn_terms[vk.bond_asgn_bounds[cdk.mol_contents[j_llim]]] -
                       vk.bond_asgn_terms[vk.bond_asgn_bounds[cdk.mol_contents[i_llim]]];
    }

    // Increment the number of molecules
    nmi++;
  }
  
  // Proceed through the list of molecules, using the pre-determined molecular matches.
  std::vector<bool> ls_covered(cdk.nmol, false);
  for (int i = 0; i < cdk.nmol; i++) {
    if (ls_covered[i]) {
      continue;
    }
    IndigoTable idg_tab(ag_pointer, i, temperature);
    const std::vector<CombineIDp> fc = idg_tab.getGroundStateFormalCharges();
    const std::vector<CombineIDp> bo = idg_tab.getGroundStateBondOrders();
    const std::vector<CombineIDp> fe = idg_tab.getGroundStateFreeElectrons();
    const std::vector<int2> zk_fc = idg_tab.getZeroKelvinFormalCharges();
    const std::vector<int2> zk_bo = idg_tab.getZeroKelvinBondOrders();
    const std::vector<int2> zk_fe = idg_tab.getZeroKelvinFreeElectrons();
    const int match_idx = mol_index[i].x;
    for (int j = i; j < cdk.nmol; j++) {
      if (mol_index[j].x == match_idx) {
        for (size_t k = 0; k < fc.size(); k++) {
          formal_charges.putHost(fc[k].y, fc[k].x + mol_index[j].y);
          zerok_formal_charges[zk_fc[k].x + mol_index[j].y] = zk_fc[k].y;
        }
        for (size_t k = 0; k < bo.size(); k++) {

          // Find the atoms to which the original instance of the bond applied, then the atoms to
          // which the bond's cognate in the current residue applies.
          const int orig_i_atom = vk.bond_i_atoms[bo[k].x];
          const int orig_j_atom = vk.bond_j_atoms[bo[k].x];
          const int next_i_atom = orig_i_atom + mol_index[j].y;
          const int next_j_atom = orig_j_atom + mol_index[j].y;

          // Find the index of the new bond
          int next_bond_index = -1;
          for (int m = vk.bond_asgn_bounds[next_i_atom]; m < vk.bond_asgn_bounds[next_i_atom + 1];
               m++) {
            if (vk.bond_asgn_atoms[m] == next_j_atom) {
              next_bond_index = vk.bond_asgn_terms[m];
              break;
            }
          }
          if (next_bond_index == -1) {
            for (int m = vk.bond_asgn_bounds[next_j_atom];
                 m < vk.bond_asgn_bounds[next_j_atom + 1]; m++) {
              if (vk.bond_asgn_atoms[m] == next_i_atom) {
                next_bond_index = vk.bond_asgn_terms[m];
                break;
              }
            }
          }
          if (next_bond_index == -1) {
            rtErr("No bond in molecule " + std::to_string(j + 1) + " could be identified as a "
                  "cognate to the bond between atoms " +
                  char4ToString(ag_pointer->getAtomName(orig_i_atom)) + " and " +
                  char4ToString(ag_pointer->getAtomName(orig_j_atom)) + " of molecule " +
                  std::to_string(i + 1) + ".", "ChemicalFeatures", "DrawLewisStructures");
          }
          bond_orders.putHost(bo[k].y, next_bond_index);
          zerok_bond_orders[next_bond_index] = zk_bo[k].y;
        }
        for (size_t k = 0; k < fe.size(); k++) {
          free_electrons.putHost(fe[k].y, fe[k].x + mol_index[j].y);
          zerok_free_electrons[zk_fe[k].x + mol_index[j].y] = zk_fe[k].y;
        }
        ls_covered[j] = true;
      }      
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findAromaticGroups(const ChemicalDetailsKit &cdk,
                                          const ValenceKit<double> &vk,
                                          const std::vector<int> &tmp_ring_atoms,
                                          const std::vector<int> &tmp_ring_atom_bounds,
                                          std::vector<int> *tmp_aromatic_group_bounds,
                                          std::vector<int> *tmp_aromatic_pi_electrons,
                                          std::vector<int> *tmp_aromatic_groups) {

  // Identify sp2- and possible sp2-centers within the ring atoms
  const int nring_atoms = tmp_ring_atoms.size();
  std::vector<bool> ring_sp2_character(nring_atoms, false);
  const double* bo_ptr = bond_orders.data();
  for (int pos = 0; pos < vk.nbond; pos++) {

    // Use a number just shy of 1 1/4 to accommodate double bonds distributed over up to four
    // atoms--it's hard to imagine such a case, but perhaps a sulfur atom making a tetrahedral
    // center between two aromatic rings at 90-degree angles to one another is physcially possible?
    // In the guanidino head group of arginine the double bond is distributed over three atoms.
    if (bo_ptr[pos] > 1.24999) {
      const int bi_atom = vk.bond_i_atoms[pos];
      const int bj_atom = vk.bond_j_atoms[pos];
      for (int i = 0; i < nring_atoms; i++) {
        ring_sp2_character[i] = (ring_sp2_character[i] || tmp_ring_atoms[i] == bi_atom ||
                                 tmp_ring_atoms[i] == bj_atom);
      }
    }
  }
  for (int i = 0; i < nring_atoms; i++) {
    const int atom_idx = tmp_ring_atoms[i];
    switch (cdk.z_numbers[atom_idx]) {
    case 1:
    case 9:
    case 17:
    case 35:

      // Hydrogen and each of the halogens are assumed to not participate in rings
      break;
    case 6:
    case 7:
    case 8:
    case 15:
    case 16:

      // Allow for a tiny amount of imprecision in whatever floating point arithmetic led to
      // the calculated free electron content.  Two free electrons can switch into the pi system
      // if this atom makes two or more bonds to other atoms, which will be the case for any atom
      // in a ring.  Even a nitrogen making two bonds to other ring atoms, with a hydrogen hanging
      // off of it (three bonds to other atoms in all) can donate a pair of electrons to the pi
      // system in the ring.  Allow one free electron to contribute, as a protonated histidine
      // ring will distribute a single positive formal charge and bonds of order 1.5 across both
      // nitrogens, leaving a lone free electron, on average, for each of them.  The individual
      // resonance states have either a full positive charge or a lone pair on either nitrogen,
      // and in either case the lone pair can contribute to an aromatic ring with a net positive
      // charge and six pi electrons.  Making inferences from the averages in this manner may
      // create some way for errors to creep in, but roll with it for now.
      ring_sp2_character[i] = (ring_sp2_character[i] ||
                               free_electrons.readHost(atom_idx) >= 0.999);
      break;
    default:
      break;
    }
  }
  
  // Determine fused rings, in preparation for aromaticity checks
  std::vector<int> ring_participation_bounds(atom_count + 1, 0);
  for (int i = 0; i < ring_count; i++) {
    for (int j = tmp_ring_atom_bounds[i]; j < tmp_ring_atom_bounds[i+1]; j++) {
      ring_participation_bounds[tmp_ring_atoms[j]] += 1;
    }
  }
  prefixSumInPlace(&ring_participation_bounds, PrefixSumType::EXCLUSIVE, "ChemicalFeatures");
  std::vector<int> ring_participation(ring_participation_bounds[atom_count]);
  for (int i = 0; i < ring_count; i++) {
    for (int j = tmp_ring_atom_bounds[i]; j < tmp_ring_atom_bounds[i+1]; j++) {
      ring_participation[ring_participation_bounds[tmp_ring_atoms[j]]] = i;
      ring_participation_bounds[tmp_ring_atoms[j]] += 1;
    }
  }
  for (int i = atom_count; i > 0; i--) {
    ring_participation_bounds[i] = ring_participation_bounds[i - 1];
  }
  ring_participation_bounds[0] = 0;

  // Assess fused ring systems
  std::vector<int> fused_rings(ring_count);
  std::vector<bool> ring_covered(ring_count, false);
  std::vector<int> fused_ring_bounds(1, 0);
  int nr = 0;
  for (int i = 0; i < ring_count; i++) {
    if (ring_covered[i]) {
      continue;
    }
    ring_covered[i] = true;
    std::vector<int> new_rings(1, i);
    int current_ring_count;
    int previous_ring_count = 0;
    do {

      // Snapshot the current number of rings
      current_ring_count = new_rings.size();

      // Loop over all current rings and add more
      for (int rpos = previous_ring_count; rpos < current_ring_count; rpos++) {
        for (int j = tmp_ring_atom_bounds[new_rings[rpos]];
             j < tmp_ring_atom_bounds[new_rings[rpos] + 1]; j++) {
          const int tring_atom = tmp_ring_atoms[j];

          // Check for atoms that are also part of an additional ring
          for (int k = ring_participation_bounds[tring_atom] + 1;
               k < ring_participation_bounds[tring_atom + 1]; k++) {
            bool extra_ring_accounted = ring_covered[ring_participation[k]];
            const int total_ring_count = new_rings.size();
            for (int m = 0; m < total_ring_count; m++) {
              extra_ring_accounted = (extra_ring_accounted ||
                                      new_rings[m] == ring_participation[k]);
            }
            if (extra_ring_accounted == false) {
              new_rings.push_back(ring_participation[k]);
              ring_covered[ring_participation[k]] = true;
            }
          }
        }
      }

      // Update the number of rings previously known, in preparation for another pass
      previous_ring_count = current_ring_count;
    } while (new_rings.size() > current_ring_count);

    // Systems with only one ring are not fused rings
    if (new_rings.size() == 1) {
      continue;
    }

    // Contribute all contiguous rings to the fused ring list
    for (size_t j = 0; j < new_rings.size(); j++) {
      ring_covered[new_rings[j]] = true;
      fused_rings[nr] = new_rings[j];
      nr++;
    }
    fused_ring_bounds.push_back(nr);
  }

  // Record the number of fused ring systems, for convenience
  fused_ring_count = fused_ring_bounds.size() - 1;
  
  // Determine fused ring atom and bond content
  std::vector<int> fused_ring_atom_counts(fused_ring_count, 0);
  std::vector<int> fused_ring_bond_counts(fused_ring_count, 0);
  std::vector<int> fused_ring_atoms;
  std::vector<int> fused_ring_atom_bounds(1, 0);
  for (int i = 0; i < fused_ring_count; i++) {

    // Make a list of all atoms
    std::vector<int> atom_list;
    for (int j = fused_ring_bounds[i]; j < fused_ring_bounds[i + 1]; j++) {
      const int jring = fused_rings[j];
      for (int k = tmp_ring_atom_bounds[jring]; k < tmp_ring_atom_bounds[jring + 1]; k++) {
        bool found = false;
        const int atom_list_length = atom_list.size();
        const int tkr_atom = tmp_ring_atoms[k];
        for (int m = 0; m < atom_list_length; m++) {
          found = (found || atom_list[m] == tkr_atom);
        }
        if (found == false) {
          atom_list.push_back(tkr_atom);
        }
      }
    }
    fused_ring_atom_counts[i] = atom_list.size();
    fused_ring_atoms.insert(fused_ring_atoms.end(), atom_list.begin(), atom_list.end());
    fused_ring_atom_bounds.push_back(fused_ring_atoms.size());

    // Referencing each atom, make a list of all bonds
    for (int j = 0; j < fused_ring_atom_counts[i]; j++) {
      const int atom_idx = atom_list[j];
      for (int k = vk.bond_asgn_bounds[atom_idx]; k < vk.bond_asgn_bounds[atom_idx + 1]; k++) {
        bool found = false;
        const int partner_atom = vk.bond_asgn_atoms[k];
        for (int m = 0; m < fused_ring_atom_counts[i]; m++) {
          found = (found || atom_list[m] == partner_atom);
        }
        if (found) {
          fused_ring_bond_counts[i] += 1;
        }
      }
    }
  }
  
  // Mark rings as aromatic or not, based on the Lewis structure, rings, and fused_rings
  std::vector<int> pi_electron_count(ring_count);
  std::vector<int> lp_electron_count(ring_count);
  for (int i = 0; i < ring_count; i++) {
    const int rb_llim = tmp_ring_atom_bounds[i];
    const int rb_hlim = tmp_ring_atom_bounds[i + 1];
    double pi_electrons = 0.0;
    double lp_electrons = 0.0;
    for (int j = rb_llim; j < rb_hlim; j++) {

      // Compute the degree of pi bonding for bonds relating to this atom lying along the ring.
      const int atom_idx = tmp_ring_atoms[j];
      for (int k = vk.bond_asgn_bounds[atom_idx]; k < vk.bond_asgn_bounds[atom_idx + 1]; k++) {
        const int partner_idx = vk.bond_asgn_atoms[k];
        bool partner_in_ring = false;
        for (int m = rb_llim; m < rb_hlim; m++) {
          partner_in_ring = (partner_in_ring || tmp_ring_atoms[m] == partner_idx);
        }
        if (partner_in_ring) {
          const int bond_idx = vk.bond_asgn_terms[k];
          pi_electrons += 2.0 * (bond_orders.readHost(bond_idx) - 1.0);
        }
      }

      // Add the free electron content
      const double nfree_e = free_electrons.readHost(atom_idx);
      if (nfree_e >= 2.0) {
        lp_electrons += 2.0;
      }
      else if (nfree_e >= 0.0) {
        lp_electrons += nfree_e;
      }
    }
    pi_electron_count[i] = round(pi_electrons);
    lp_electron_count[i] = round(lp_electrons);
  }
  for (int i = 0; i < fused_ring_count; i++) {

    // Immediately eliminate fused rings that have less than five atoms, or that have one or more
    // atoms which cannot take on sp2 character.
    if (fused_ring_atom_counts[i] < 5) {
      continue;
    }
    const int frb_llim = fused_ring_bounds[i];
    const int frb_hlim = fused_ring_bounds[i + 1];
    bool ring_qualifies = true;
    for (int j = frb_llim; j < frb_hlim; j++) {
      const int jring = fused_rings[j];
      for (int k = tmp_ring_atom_bounds[jring]; k < tmp_ring_atom_bounds[jring + 1]; k++) {
        ring_qualifies = (ring_qualifies && ring_sp2_character[k]);
      }
    }
    if (ring_qualifies == false) {
      continue;
    }

    // The numbers of pi and lone pair electrons are the sum of each of the constituent rings,
    // without double-counting bonds or atoms shared by both rings.
    int pi_electrons = 0;
    int lp_electrons = 0;
    std::vector<int> ring_occupancy(fused_ring_atom_counts[i], 0);
    for (int j = frb_llim; j < frb_hlim; j++) {
      const int jring = fused_rings[j];      
      pi_electrons += pi_electron_count[jring];
      lp_electrons += lp_electron_count[jring];
      for (int k = tmp_ring_atom_bounds[jring]; k < tmp_ring_atom_bounds[jring + 1]; k++) {
        const int m_start = fused_ring_atom_bounds[i];
        for (int m = m_start; m < fused_ring_atom_bounds[i + 1]; m++) {
          ring_occupancy[m - m_start] += (fused_ring_atoms[m] == tmp_ring_atoms[k]);
        }
      }
    }
    for (int j = 0; j < fused_ring_atom_counts[i]; j++) {
      if (ring_occupancy[j] == 1) {
        continue;
      }
      const int n_extra = ring_occupancy[j] - 1;
      const int atomj_idx = fused_ring_atoms[fused_ring_atom_bounds[i] + j];
      const int nfree_e = static_cast<int>(round(free_electrons.readHost(atomj_idx)));
      if (nfree_e >= 2) {
        lp_electrons -= 2 * n_extra;
      }
      else if (nfree_e >= 0) {
        lp_electrons -= nfree_e * n_extra;
      }

      // Seek out bonds between atoms that are included in similar numbers of rings within
      // this fused system.
      for (int k = j + 1; k < fused_ring_atom_counts[i]; k++) {
        if (ring_occupancy[k] == ring_occupancy[j]) {
          const int atomk_idx = fused_ring_atoms[fused_ring_atom_bounds[i] + k];
          for (int m = vk.bond_asgn_bounds[atomk_idx]; m < vk.bond_asgn_bounds[atomk_idx + 1];
               m++) {
            if (vk.bond_asgn_atoms[m] == atomj_idx) {
              const double tbo = bond_orders.readHost(vk.bond_asgn_terms[m]);
              pi_electrons -= static_cast<int>(round(2.0 * (tbo - 1.0))) * n_extra;
            }
          }
          for (int m = vk.bond_asgn_bounds[atomj_idx]; m < vk.bond_asgn_bounds[atomj_idx + 1];
               m++) {
            const double tbo = bond_orders.readHost(vk.bond_asgn_terms[m]);
            pi_electrons -= static_cast<int>(round(2.0 * (tbo - 1.0))) * n_extra;
          }
        }
      }
    }
    
    // Is the ring itself aromatic?  Apply Huckel's rule to test whether there are 4n + 2 pi
    // electrons.
    if (((pi_electrons + lp_electrons) & 0x1) == 0 && (pi_electrons + lp_electrons) / 4 >= 1) {
      for (int j = fused_ring_atom_bounds[i]; j < fused_ring_atom_bounds[i + 1]; j++) {
        tmp_aromatic_groups->push_back(fused_ring_atoms[j]);
      }
      tmp_aromatic_group_bounds->push_back(tmp_aromatic_groups->size());
      tmp_aromatic_pi_electrons->push_back(pi_electrons + lp_electrons);
    }
  }
  for (int i = 0; i < ring_count; i++) {

    // Immediately eliminate rings that have less than five atoms, or that have one or more
    // atoms which cannot take on sp2 character.
    const int rb_llim = tmp_ring_atom_bounds[i];
    const int rb_hlim = tmp_ring_atom_bounds[i + 1];
    const int nring_atom = rb_hlim - rb_llim;
    if (nring_atom < 5) {
      continue;
    }
    bool ring_qualifies = true;
    for (int j = rb_llim; j < rb_hlim; j++) {
      ring_qualifies = (ring_qualifies && ring_sp2_character[j]);
    }
    if (ring_qualifies == false) {
      continue;
    }    
    const int pi_electrons = pi_electron_count[i];
    const int lp_electrons = lp_electron_count[i];
    
    // Is the ring itself aromatic?  Apply Huckel's rule to test whether there are 4n + 2 pi
    // electrons.
    if (((pi_electrons + lp_electrons) & 0x1) == 0 && (pi_electrons + lp_electrons) / 4 >= 1) {
      for (int j = rb_llim; j < rb_hlim; j++) {
        tmp_aromatic_groups->push_back(tmp_ring_atoms[j]);
      }
      tmp_aromatic_group_bounds->push_back(tmp_aromatic_groups->size());
    }
  }

  // Record the total number of groups, for convenience
  aromatic_group_count = tmp_aromatic_group_bounds->size() - 1;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::detailChiralCenters(const NonbondedKit<double> &nbk,
                                                       const ValenceKit<double> &vk,
                                                       const ChemicalDetailsKit &cdk) {
  
  // Prepare to construct tree structures for each molecule, similar to what was done above.
  // However, the trees will not trace rings this time.  Rather, they will continue along the path
  // radiating outwards from any given atom which has four unique substituents, trying to determine
  // whether one is superior to another in terms of atomic numbers and bond orders.
  int max_branches = 0;
  for (int i = 0; i < atom_count; i++) {
    max_branches = std::max(nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i], max_branches);
  }
  const int max_molecule = ag_pointer->getLargestMoleculeSize();
  std::vector<std::vector<int>> all_branch_atoms(4, std::vector<int>(max_branches * max_molecule));
  std::vector<std::vector<BondedNode>> links(4, std::vector<BondedNode>(atom_count));
  std::vector<std::vector<int>> tree_positions(4, std::vector<int>(atom_count, -1));
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < max_molecule; j++) {
      links[i][j].setBranchPointer(&all_branch_atoms[i], j, max_branches);
    }
  }

  // Loop over all atoms, proceeding to explore all available bonds, until all atoms have been
  // either used as the start of chain / ring exploration or have been included in the exploration
  // initiated by some other atom.
  std::vector<std::vector<bool>> atom_touched(4, std::vector<bool>(atom_count, false));
  std::vector<int> result;
  std::vector<int4> tmp_chiral_arms;
  std::vector<bool> valid_atom_mask(atom_count);
  for (int i = 0; i < atom_count; i++) {
    valid_atom_mask[i] = (cdk.z_numbers[i] > 0);
  }
  for (int i = 0; i < atom_count; i++) {

    // Test whether this atom is a chiral candidate
    if (nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i] < 4) {
      continue;
    }
    bool candidate = true;
    int n_hydrogen = 0;
    int n_vs = 0;
    for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
      const int ex_znum = cdk.z_numbers[nbk.nb12x[j]];
      n_hydrogen += (ex_znum == 1);
      n_vs += (ex_znum == 0);
      candidate = (candidate && n_hydrogen < 2);
    }
    if (candidate == false || nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i] - n_vs != 4) {
      continue;
    }
    
    // Initiate the chains
    int n_real = 0;
    std::vector<int> layer_llim(4, 0);
    std::vector<int> layer_hlim(4, 1);
    std::vector<int> node_count(4, 1);
    std::vector<int> branch_ranks(4, 0);
    std::vector<int> branch_scores(4, 0);
    std::vector<int> touch_min(4);
    std::vector<int> touch_max(4);
    std::vector<bool> score_branch(4);
    int current_layer = 1;
    for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
      if (cdk.z_numbers[nbk.nb12x[j]] == 0) {
        continue;
      }
      links[n_real][0].addToTree(i, nbk.nb12x[j], 0, nbk, valid_atom_mask);
      links[n_real][0].addBondOrder(vk, bond_orders);
      atom_touched[n_real][i] = true;
      atom_touched[n_real][nbk.nb12x[j]] = true;
      touch_min[n_real] = std::min(i, nbk.nb12x[j]);
      touch_max[n_real] = std::max(i, nbk.nb12x[j]);
      n_real++;
    }

    // Once a branch beats another, that dominance relationship must be maintained.  Keep a matrix
    // of the dominance relationships (i, j) : 0 = tie, +1 = i beats j, -1 = j beats i.  The matrix
    // will be antisymmetric.
    std::vector<int> chiral_dominance(16, 0);
    std::vector<int> parallel_growth(16, 0);
    bool advance_chains = scoreChiralBranches(links, layer_llim, layer_hlim, cdk,
                                              &chiral_dominance, &parallel_growth);
    while (advance_chains) {
      
      // Add more nodes
      for (int j = 0; j < 4; j++) {
        for (int k = layer_llim[j]; k < layer_hlim[j]; k++) {
          const int k_atom = links[j][k].getAtom();
          const int k_branch_count = links[j][k].getBranchCount();
          for (int m = 0; m < k_branch_count; m++) {
            const int m_atom = links[j][k].getBranchAtom(m);
            if (! atom_touched[j][m_atom]) {
              atom_touched[j][m_atom] = true;
              touch_min[j] = std::min(touch_min[j], m_atom);
              touch_max[j] = std::max(touch_max[j], m_atom);
              links[j][node_count[j]].addToTree(k_atom, m_atom, current_layer, nbk,
                                                valid_atom_mask);
              links[j][node_count[j]].addBondOrder(vk, bond_orders);
              tree_positions[j][m_atom] = node_count[j];
              node_count[j] += 1;
            }
          }
        }
        layer_llim[j] = layer_hlim[j];
        layer_hlim[j] = node_count[j];
      }
      current_layer++;
      advance_chains = scoreChiralBranches(links, layer_llim, layer_hlim, cdk, &chiral_dominance,
                                           &parallel_growth);
    }

    // Check the dominance matrix: are the off-diagonal elements all nonzero?
    if (chiral_dominance[ 1] == 0 || chiral_dominance[ 2] == 0 || chiral_dominance[ 3] == 0 ||
        chiral_dominance[ 6] == 0 || chiral_dominance[ 7] == 0 || chiral_dominance[11] == 0) {
      continue;      
    }

    // Give the next chirality check a clean slate.  Use the touch_min and touch_max bounds
    // arrays on the range of atoms that became part of each branch's tree to keep us in O(N)
    // territory when doing things like the chirality of protein CA atoms.
    for (int j = 0; j < 4; j++) {
      for (int k = touch_min[j]; k < touch_max[j]; k++) {
        atom_touched[j][k] = false;
      }
    }

    // Check the dominance matrix for branch priorities.  Higher score = higher priority.
    std::vector<int2> priority(4, {0, 0});
    for (int j = 0; j < 4; j++) {
      priority[j].x += (chiral_dominance[j     ] == 1) + (chiral_dominance[j +  4] == 1) +
                       (chiral_dominance[j +  8] == 1) + (chiral_dominance[j + 12] == 1);
      priority[j].y = j;
    }
    std::sort(priority.begin(), priority.end(), [](int2 a, int2 b) { return a.x < b.x; });

    // Make the coordinate origin the atom center.  Align the atoms with respect to the axis
    // extending between the center and the root of the lowest priority branch.
    int4 ichiral_arms = { links[priority[0].y][0].getAtom(),
                          links[priority[3].y][0].getAtom(),
                          links[priority[2].y][0].getAtom(),
                          links[priority[1].y][0].getAtom() };
    tmp_chiral_arms.push_back(ichiral_arms);
    result.push_back(i + 1);
  }

  // Post the chiral arm information now, in its own ARRAY-kind Hybrid object.
  chiral_arm_atoms.resize(tmp_chiral_arms.size());
  chiral_arm_atoms.putHost(tmp_chiral_arms);
  
  // Return the list of chiral atom indices.  The atom indices have all been inflated by one, as
  // if they are all L-orientation.  The orientations will be corrected in the next call to
  // findChiralOrientations().
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ChemicalFeatures::findChiralInversionMethods(const std::vector<int> &tmp_chiral_centers,
                                             const std::vector<int> &tmp_ring_atoms,
                                             const std::vector<int> &tmp_ring_atom_bounds) {
  std::vector<int> result(chiral_center_count);

  // Count the number of rings, irrespective of size, that each atom participates in.
  // Chiral centers which participate in two rings can be inverted using the standard protocol
  // of C2 symmetry rotation so long as those two rings only share one atom, the chiral center
  // itself.  If the two rings share more atoms than just the chiral center, then the additional
  // bindings between the rings prevent a clean rotation of two branches coming out of the center
  // and force inversion by complete reflection of the molecule across some plane.  Chiral atoms
  // which participate in three rings cannot be inverted by C2 symmtery rotation and must be
  // inverted by complete reflection of the entire molecule across a plane.
  std::vector<int> ring_occupancy_bounds(atom_count + 1, 0);
  for (int i = 0; i < tmp_ring_atom_bounds[ring_count]; i++) {
    ring_occupancy_bounds[tmp_ring_atoms[i]] += 1;
  }
  prefixSumInPlace(&ring_occupancy_bounds, PrefixSumType::EXCLUSIVE, "ChemicalFeatures");
  std::vector<int> ring_occupancy(ring_occupancy_bounds[atom_count]);
  for (int i = 0; i < ring_count; i++) {
    for (int j = tmp_ring_atom_bounds[i]; j < tmp_ring_atom_bounds[i + 1]; j++) {
      const int cpos = ring_occupancy_bounds[tmp_ring_atoms[j]];
      ring_occupancy[cpos] = i;
      ring_occupancy_bounds[tmp_ring_atoms[j]] = cpos + 1;
    }
  }
  for (int i = atom_count; i > 0; i--) {
    ring_occupancy_bounds[i] = ring_occupancy_bounds[i - 1];
  }
  ring_occupancy_bounds[0] = 0;
  bool reflection_taken = false;
  for (int i = 0; i < chiral_center_count; i++) {

    // The chiral center indices are encoded with D- and L- orientations by offsetting all of
    // their atom indices by +1 and then multiplying D- center indices by -1.
    const int chatom = abs(tmp_chiral_centers[i] - 1);
    const int nrings = (ring_occupancy_bounds[chatom + 1] - ring_occupancy_bounds[chatom]);
    if (nrings >= 3) {
      result[i] = (reflection_taken) ? static_cast<int>(ChiralInversionProtocol::DO_NOT_INVERT) :
                                       static_cast<int>(ChiralInversionProtocol::REFLECT);
      reflection_taken = true;
    }
    else if (nrings == 2) {
      int nfusion = 0;
      const int r1_idx = ring_occupancy[ring_occupancy_bounds[chatom]];
      const int r2_idx = ring_occupancy[ring_occupancy_bounds[chatom] + 1];
      for (int j = tmp_ring_atom_bounds[r1_idx]; j < tmp_ring_atom_bounds[r1_idx + 1]; j++) {
        const int jring_atom = tmp_ring_atoms[j];
        for (int k = tmp_ring_atom_bounds[r2_idx]; k < tmp_ring_atom_bounds[r2_idx + 1]; k++) {
          nfusion += (tmp_ring_atoms[k] == jring_atom);
        }
      }
      if (nfusion == 1) {
        result[i] = static_cast<int>(ChiralInversionProtocol::ROTATE);
      }
      else {
        result[i] = (reflection_taken) ? static_cast<int>(ChiralInversionProtocol::DO_NOT_INVERT) :
                                         static_cast<int>(ChiralInversionProtocol::REFLECT);
        reflection_taken = true;
      }
    }
    else {
      result[i] = static_cast<int>(ChiralInversionProtocol::ROTATE);
    }
  }
  
  return result;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findRotatableBonds(const ValenceKit<double> &vk,
                                          const ChemicalDetailsKit &cdk,
                                          const NonbondedKit<double> &nbk,
                                          const std::vector<int> &ring_atoms,
                                          const std::vector<int> &ring_atom_bounds,
                                          std::vector<int> *tmp_rotatable_groups,
                                          std::vector<int> *tmp_rotatable_group_bounds,
                                          std::vector<int> *tmp_cis_trans_groups,
                                          std::vector<int> *tmp_cis_trans_group_bounds) {
  
  // Scan over all bonds and accumulate a result
  Approx near_one(1.0, 0.21);
  std::vector<int2> rotators, cistrans;
  std::vector<std::vector<int>> rot_moving_lists, ctx_moving_lists;
  for (int pos = 0; pos < vk.nbond; pos++) {
    
    // Omit bonds within rings.  Those will be handled separately.  Otherwise, allow a liberal
    // definition of a single bond--a small amount of double-bond character might be permissible.
    if (bond_in_ring[pos] || bond_orders.readHost(pos) > 2.01) {
      continue;
    }

    // If the bond order is near one, this is a freely rotatable bond.  Otherwise, the bond is a
    // cis-trans invertible group.
    const bool is_standard_rotator = (bond_orders.readHost(pos) == near_one);

    // Ensure that both ends have more branching from them, and that the branches are worth
    // rotating (more than just hydrogen atoms branching from them)
    const int atom_i = vk.bond_i_atoms[pos];
    const int atom_j = vk.bond_j_atoms[pos];
    int nbranch_i = 0;
    for (int i = nbk.nb12_bounds[atom_i]; i < nbk.nb12_bounds[atom_i + 1]; i++) {
      nbranch_i += (nbk.nb12x[i] != atom_j && cdk.z_numbers[nbk.nb12x[i]] > 1);
    }
    int nbranch_j = 0;
    for (int i = nbk.nb12_bounds[atom_j]; i < nbk.nb12_bounds[atom_j + 1]; i++) {
      nbranch_j += (nbk.nb12x[i] != atom_i && cdk.z_numbers[nbk.nb12x[i]] > 1);
    }    
    if (nbranch_i < 1 || nbranch_j < 1) {
      continue;
    }

    // Test the number of rotatable bonds with atom_i as the root and atom_j as the pivot (the
    // pivot atom is closest to the atoms that will move), then vice-versa.  List the atoms such
    // that the smallest number of atoms will move as a consequence of rotation about the bond.
    // Store the rotating atoms.
    const std::vector<int> ifirst = selectRotatingAtoms(ag_pointer, atom_i, atom_j);
    const std::vector<int> jfirst = selectRotatingAtoms(ag_pointer, atom_j, atom_i);
    bool i_is_fixed = false;
    bool j_is_fixed = false;
    const int n_iatom = ifirst.size();
    const int n_jatom = jfirst.size();
    for (int pos = 0; pos < n_iatom; pos++) {
      i_is_fixed = (i_is_fixed || ag_pointer->getAtomMobility(ifirst[pos]) == false);
    }
    for (int pos = 0; pos < n_jatom; pos++) {
      j_is_fixed = (j_is_fixed || ag_pointer->getAtomMobility(jfirst[pos]) == false);
    }
    if (j_is_fixed == false && (ifirst.size() >= jfirst.size() || i_is_fixed)) {
      if (is_standard_rotator) {
        rotators.push_back({atom_j, atom_i});
        rot_moving_lists.push_back(jfirst);
      }
      else {
        cistrans.push_back({atom_j, atom_i});
        ctx_moving_lists.push_back(jfirst);
      }
    }
    else if (i_is_fixed == false) {
      if (is_standard_rotator) {
        rotators.push_back({atom_i, atom_j});
        rot_moving_lists.push_back(ifirst);
      }
      else {
        cistrans.push_back({atom_i, atom_j});
        ctx_moving_lists.push_back(ifirst);
      }
    }
  }

  // Log rotatable bond groups
  unpackRotatableGroups(rotators, rot_moving_lists, tmp_rotatable_groups,
                        tmp_rotatable_group_bounds);
  unpackRotatableGroups(cistrans, ctx_moving_lists, tmp_cis_trans_groups,
                        tmp_cis_trans_group_bounds);

  // Signal that the rotatable groups have been mapped
  rotating_groups_mapped = true;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findInvertibleGroups(const std::vector<int> &tmp_chiral_centers,
                                            const std::vector<int> &tmp_inversion_methods,
                                            std::vector<int> *tmp_anchor_a_branches,
                                            std::vector<int> *tmp_anchor_b_branches,
                                            std::vector<int> *tmp_invertible_groups,
                                            std::vector<int> *tmp_invertible_group_bounds) {

  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  
  // Loop over all chiral centers and determine the number of atoms in each branch
  tmp_anchor_a_branches->resize(chiral_center_count);
  tmp_anchor_b_branches->resize(chiral_center_count);
  tmp_invertible_group_bounds->resize(chiral_center_count + 1);
  int* anchor_a_ptr = tmp_anchor_a_branches->data();
  int* anchor_b_ptr = tmp_anchor_b_branches->data();
  int* bounds_ptr   = tmp_invertible_group_bounds->data();
  const size_t mask_len = (nbk.natom + uint_bit_count_int - 1) / uint_bit_count_int;
  std::vector<std::vector<uint>> marked(4);
  std::vector<int3> branch_counts(4);
  for (int i = 0; i < 4; i++) {
    marked[i].resize(mask_len);
  }
  std::vector<int> prev_atoms(16), new_atoms(16), tmp_igroup(16);
  int all_grp_size = 0;
  for (int i = 0; i < chiral_center_count; i++) {
    
    // The chiral orientation is encoded in the index of the chiral center by making D-chiral
    // centers the negative of the original index, but this leaves an ambiguity if the atom with
    // topological index 0 is chiral.  Therefore, 1 is added to the index and D-chiral centers
    // have their index values multiplied by -1.  Unroll this arrangement to recover the correct
    // topological index.
    const int chatom = abs(tmp_chiral_centers[i]) - 1;

    // Switch over the various inversion methods--methods imply group content
    switch (static_cast<ChiralInversionProtocol>(tmp_inversion_methods[i])) {
    case ChiralInversionProtocol::ROTATE:
      {
        for (int j = 0; j < 4; j++) {
          branch_counts[j].x = -1;
          branch_counts[j].y = j;
        }
        int brpos = 0;
        int nring = 0;
        int n_j_vs_encountered = 0;
        for (int j = nbk.nb12_bounds[chatom]; j < nbk.nb12_bounds[chatom + 1]; j++) {

          // Skip "branches" that are actually virtual sites
          if (cdk.z_numbers[nbk.nb12x[j]] == 0) {
            n_j_vs_encountered++;
            continue;
          }

          // Skip branches making up the other end of a loop which has already been determined
          if (branch_counts[brpos].x >= 0) {
            brpos++;
            continue;
          }
          bool ring_completed;
          branch_counts[brpos].x = colorConnectivity(nbk, cdk, chatom, nbk.nb12x[j],
                                                     &marked[brpos], &ring_completed);
          
          // Search for one of the remaining branch roots in the marked atoms if there was a ring
          // completion.  This will accomplish one of the following searches.  Mark the fact that
          // this branch (and its partner, from the other end of the ring) are part of the same
          // ring.  If one of these two branches ends up as one of the ones that rotates, the
          // other must rotate along with it.
          if (ring_completed) {
            nring++;
            branch_counts[brpos].z = nring;
            int n_k_vs_encountered = n_j_vs_encountered;
            for (int k = j + 1; k < nbk.nb12_bounds[chatom + 1]; k++) {

              // Again, skip "branches" that are actually virtual sites
              if (cdk.z_numbers[nbk.nb12x[k]] == 0) {
                n_k_vs_encountered++;
                continue;
              }
              if (readBitFromMask(marked[brpos], nbk.nb12x[k])) {
                const size_t mirror_idx = k - nbk.nb12_bounds[chatom] - n_k_vs_encountered;
                for (size_t m = 0; m < mask_len; m++) {
                  marked[mirror_idx][m] = marked[brpos][m];
                }
                branch_counts[mirror_idx].x = branch_counts[brpos].x;
                branch_counts[mirror_idx].z = nring;
              }
            }
          }
          else {
            branch_counts[brpos].z = 0;
          }
          brpos++;
        }
        
        // Unset the chiral atom itself in the masks just created.  It will not rotate during the
        // inversion.  The number of atoms in each branch is one more than the number of "rotators"
        // counted by the colorConnectivity function, which was originally written for mapping the
        // moving atom dependencies of rotatable bonds.
        for (int j = 0; j < 4; j++) {
          branch_counts[j].x += 1;
          unsetBitInMask(&marked[j], chatom);
        }

        // Sort the branches in descending order of size.  The list is only four elements long.
        std::sort(branch_counts.begin(), branch_counts.end(),
                  [](int3 a, int3 b) { return a.x > b.x; });

        // In absense of rings, the shortest two branches will do.  If one of the branches is
        // part of a ring, however, the second branch in that ring must also be taken.
        int invr_a_idx, invr_b_idx;
        if (nring == 0 || (branch_counts[0].z == 0 && branch_counts[1].z == 0) ||
            branch_counts[0].z == branch_counts[1].z) {
          anchor_a_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[0].y];
          anchor_b_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[1].y];
          invr_a_idx = 2;
          invr_b_idx = 3;
        }
        else if (nring == 1) {

          // If there is one ring and branch zero were part of it, branch one would also be a part
          // of the ring due to the ordering and the anchors would have been determined above.
          // Both branches of a ring will have the same number of dependencies and will therefor
          // appear back-to-back in the ordered list.
          if (branch_counts[0].x + branch_counts[3].x <= branch_counts[1].x + branch_counts[2].x) {
            anchor_a_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[0].y];
            anchor_b_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[3].y];  
            invr_a_idx = 1;
            invr_b_idx = 2;
          }
          else {
            anchor_a_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[1].y];
            anchor_b_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[2].y];  
            invr_a_idx = 0;
            invr_b_idx = 3;
          }
        }
        else if (nring >= 2) {

          // If there are two rings that are fused in more places than just the chiral atom, the
          // center will ultimately be marked for reflection rather than rotation.  Otherwise, take
          // the ring with fewer atoms and mark that as the inversion group (this has already been
          // done in the first case above).  If there are more than two rings, the center will have
          // to be reflected.
          anchor_a_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[0].y];
          anchor_b_ptr[i] = nbk.nb12x[nbk.nb12_bounds[chatom] + branch_counts[1].y];
          invr_a_idx = 2;
          invr_b_idx = 3;
        }
        const int invr_a_branch = branch_counts[invr_a_idx].y;
        const int invr_b_branch = branch_counts[invr_b_idx].y;
        tmp_igroup.resize(branch_counts[invr_a_idx].x + branch_counts[invr_b_idx].x + 2);
    
        // Quickly scan the relevant masks for any marked atoms.  This uses the unsigned integers
        // to account for 32 atoms at a time.  In large topologies with relatively few atoms in the
        // inversion mask, this is a big speed advantage.  Make a list of atoms in either branch,
        // then reduce this list to its unique values, unique atoms in the smallest two branches
        // coming out of the chiral center.
        tmp_igroup[0] = nbk.nb12x[nbk.nb12_bounds[chatom] + invr_a_branch];
        tmp_igroup[1] = nbk.nb12x[nbk.nb12_bounds[chatom] + invr_b_branch];
        int ntg = 2;
        for (int j = 0; j < mask_len; j++) {
          if (marked[invr_a_branch][j] > 0) {
            const int jbase = j * uint_bit_count_int;
            for (int k = 0; k < uint_bit_count_int; k++) {
              if (readBitFromMask(marked[invr_a_branch], jbase + k)) {
                tmp_igroup[ntg] = jbase + k;
                ntg++;
              }
            }
          }
          if (marked[invr_b_branch][j] > 0) {
            const int jbase = j * uint_bit_count_int;
            for (int k = 0; k < uint_bit_count_int; k++) {
              if (readBitFromMask(marked[invr_b_branch], jbase + k)) {
                tmp_igroup[ntg] = jbase + k;
                ntg++;
              }
            }
          }
        }
        reduceUniqueValues(&tmp_igroup);
      }
      break;
    case ChiralInversionProtocol::REFLECT:
      {
        // In the special case of a chiral center than can only be inverted by mirroring the whole
        // molecule, list all atoms.
        const int center_home = cdk.mol_home[chatom];
        tmp_igroup.resize(cdk.mol_limits[center_home + 1] - cdk.mol_limits[center_home] + 2);
        int k = 0;
        for (int j = cdk.mol_limits[center_home]; j < cdk.mol_limits[center_home + 1]; j++) {
          tmp_igroup[k] = cdk.mol_contents[j];
          k++;
        }
      }
      break;
    case ChiralInversionProtocol::DO_NOT_INVERT:
      tmp_igroup.resize(0);
      break;
    }

    // Add the temporary atom group to the growing list of inversion instructions
    const int igrp_size = tmp_igroup.size();
    all_grp_size += igrp_size;
    bounds_ptr[i + 1] = all_grp_size;

    // This resizing operation is potentially inefficient as it happens for each chiral center,
    // but it is very hard to anticipate what a reasonable size for the entire array might be.
    tmp_invertible_groups->resize(all_grp_size);
    int* group_ptr = tmp_invertible_groups->data();
    for (int j = 0; j < igrp_size; j++) {
      group_ptr[bounds_ptr[i] + j] = tmp_igroup[j];
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findHydrogenBondElements(const NonbondedKit<double> &nbk,
                                                const ChemicalDetailsKit &cdk,
                                                std::vector<int> *tmp_polar_h,
                                                std::vector<int> *tmp_hb_don,
                                                std::vector<int> *tmp_hb_acc) {
  std::vector<bool> polarh_covered(cdk.natom, false);
  std::vector<bool> donor_covered(cdk.natom, false);
  for (int i = 0; i < cdk.natom; i++) {
    if (cdk.z_numbers[i] == 1) {
      for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
        const int donor_atom = nbk.nb12x[j];
        const int donor_z = cdk.z_numbers[donor_atom];
        if (donor_z == 7 || donor_z == 8 || donor_z == 15 || donor_z == 16) {
          if (polarh_covered[i] == false) {
            tmp_polar_h->push_back(i);
            polarh_covered[i] = true;
          }
          if (donor_covered[donor_atom] == false) {
            tmp_hb_don->push_back(donor_atom);
            donor_covered[donor_atom] = true;
          }
        }
      }
    }
    else if ((cdk.z_numbers[i] ==  7 || cdk.z_numbers[i] ==  8 || cdk.z_numbers[i] == 15 |
              cdk.z_numbers[i] == 16) && free_electrons.readHost(i) > 1.25) {

      // Nitrogen, oxygen, phosphorus, and sulfur atoms with at least some lone pair occupancy
      // are considered hydrogen bond acceptors.
      tmp_hb_acc->push_back(i);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::allocateMutableData(const std::vector<int> &tmp_rotatable_groups,
                                           const std::vector<int> &tmp_rotatable_group_bounds,
                                           const std::vector<int> &tmp_cis_trans_groups,
                                           const std::vector<int> &tmp_cis_trans_group_bounds,
                                           const std::vector<int> &tmp_invertible_groups,
                                           const std::vector<int> &tmp_invertible_group_bounds,
                                           const std::vector<int> &tmp_anchor_a_branches,
                                           const std::vector<int> &tmp_anchor_b_branches) {
  const size_t nmuta = roundUp(tmp_rotatable_groups.size(), warp_size_zu) +
                       roundUp(tmp_rotatable_group_bounds.size(), warp_size_zu) +
                       roundUp(tmp_cis_trans_groups.size(), warp_size_zu) +
                       roundUp(tmp_cis_trans_group_bounds.size(), warp_size_zu) +
                       roundUp(tmp_invertible_groups.size(), warp_size_zu) +
                       roundUp(tmp_invertible_group_bounds.size(), warp_size_zu) +
                       (2 * roundUp(static_cast<size_t>(chiral_center_count), warp_size_zu));
  mutable_group_data.resize(nmuta);
  size_t ic = rotatable_groups.putHost(&mutable_group_data, tmp_rotatable_groups, 0, warp_size_zu);
  ic = rotatable_group_bounds.putHost(&mutable_group_data, tmp_rotatable_group_bounds, ic,
                                      warp_size_zu);
  ic = cis_trans_groups.putHost(&mutable_group_data, tmp_cis_trans_groups, ic, warp_size_zu);
  ic = cis_trans_group_bounds.putHost(&mutable_group_data, tmp_cis_trans_group_bounds, ic,
                                      warp_size_zu);
  ic = invertible_groups.putHost(&mutable_group_data, tmp_invertible_groups, ic, warp_size_zu);
  ic = invertible_group_bounds.putHost(&mutable_group_data, tmp_invertible_group_bounds, ic,
                                       warp_size_zu);
  ic = anchor_a_branches.putHost(&mutable_group_data, tmp_anchor_a_branches, ic, warp_size_zu);
  ic = anchor_b_branches.putHost(&mutable_group_data, tmp_anchor_b_branches, ic, warp_size_zu);
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::repairPointers() {
  formal_charges.swapTarget(&double_data);
  free_electrons.swapTarget(&double_data);
  bond_orders.swapTarget(&double_data);
  planar_centers.swapTarget(&int_data);
  ring_atom_bounds.swapTarget(&int_data);
  ring_atoms.swapTarget(&int_data);
  aromatic_group_bounds.swapTarget(&int_data);
  aromatic_pi_electrons.swapTarget(&int_data);
  aromatic_groups.swapTarget(&int_data);
  polar_hydrogens.swapTarget(&int_data);
  hydrogen_bond_donors.swapTarget(&int_data);
  hydrogen_bond_acceptors.swapTarget(&int_data);
  chiral_centers.swapTarget(&int_data);
  chiral_inversion_methods.swapTarget(&int_data);
  rotatable_groups.swapTarget(&mutable_group_data);
  rotatable_group_bounds.swapTarget(&mutable_group_data);
  cis_trans_groups.swapTarget(&mutable_group_data);
  cis_trans_group_bounds.swapTarget(&mutable_group_data);
  invertible_groups.swapTarget(&mutable_group_data);
  invertible_group_bounds.swapTarget(&mutable_group_data);
  anchor_a_branches.swapTarget(&mutable_group_data);
  anchor_b_branches.swapTarget(&mutable_group_data);
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::upload() {
  ring_inclusion.upload();
  int_data.upload();
  mutable_group_data.upload();
  double_data.upload();
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::download() {
  ring_inclusion.download();
  int_data.download();
  mutable_group_data.download();
  double_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getPlanarAtomCount() const {
  return planar_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getRingCount() const {
  return ring_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getFusedRingCount() const {
  return fused_ring_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getMutableRingCount() const {
  return twistable_ring_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getAromaticGroupCount() const {
  return aromatic_group_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getPolarHydrogenCount() const {
  return polar_hydrogen_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getHydrogenBondDonorCount() const {
  return hbond_donor_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getHydrogenBondAcceptorCount() const {
  return hbond_acceptor_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getChiralCenterCount() const {
  return chiral_center_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getRotatableBondCount() const {
  return rotatable_bond_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getCisTransBondCount() const {
  return cis_trans_bond_count;
}

//-------------------------------------------------------------------------------------------------
bool ChemicalFeatures::rotatableGroupsMapped() const {
  return rotating_groups_mapped;
}

//-------------------------------------------------------------------------------------------------
bool ChemicalFeatures::chiralitiesComputed() const {
  return chiralities_computed;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getRingMask(const int min_ring_size,
                                                const int max_ring_size) const {
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (atom_count + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);
  
  // Loop over all known rings in the system and select those meeting the size criterion
  const int* rb_ptr = ring_atom_bounds.data();
  const int* ra_ptr = ring_atoms.data();
  for (int i = 0; i < ring_count; i++) {
    const int rsize = rb_ptr[i + 1] - rb_ptr[i];
    if (rsize >= min_ring_size && rsize <= max_ring_size) {
      for (int j = rb_ptr[i]; j < rb_ptr[i + 1]; j++) {
        accumulateBitmask(&result, ra_ptr[j]);
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getAromaticMask(const int min_pi_electrons,
                                                    const int max_pi_electrons) const {
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (atom_count + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);
  
  // Loop over all known rings in the system and select those meeting the size criterion
  const int* amb_ptr = aromatic_group_bounds.data();
  const int* amp_ptr = aromatic_pi_electrons.data();
  const int* ama_ptr = aromatic_groups.data();
  for (int i = 0; i < ring_count; i++) {
    if (amp_ptr[i] >= min_pi_electrons && amp_ptr[i] <= max_pi_electrons) {
      for (int j = amb_ptr[i]; j < amb_ptr[i + 1]; j++) {
        accumulateBitmask(&result, ama_ptr[j]);
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getPolarHydrogenList() const {
  return polar_hydrogens.readHost(0, polar_hydrogen_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getHydrogenBondDonorList() const {
  return hydrogen_bond_donors.readHost(0, hbond_donor_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getHydrogenBondAcceptorList() const {
  return hydrogen_bond_acceptors.readHost(0, hbond_acceptor_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getPolarHydrogenMask() const {
  return numberSeriesToBitMask(polar_hydrogens, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getHydrogenBondDonorMask() const {
  return numberSeriesToBitMask(hydrogen_bond_donors, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getHydrogenBondAcceptorMask() const {
  return numberSeriesToBitMask(hydrogen_bond_acceptors, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getChiralCenters(const ChiralOrientation direction) const {
  std::vector<int> result;
  const int* chi_ptr = chiral_centers.data();
  int ncen;
  switch (direction) {
  case ChiralOrientation::RECTUS:
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      ncen += (chi_ptr[i] < 0);
    }
    result.resize(ncen);
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      if (chi_ptr[i] < 0) {
        result[ncen] = 1 - chi_ptr[i];
        ncen++;
      }
    }
    break;
  case ChiralOrientation::SINISTER:
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      ncen += (chi_ptr[i] > 0);
    }
    result.resize(ncen);
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      if (chi_ptr[i] > 0) {
        result[ncen] = chi_ptr[i] - 1;
        ncen++;
      }
    }
    break;
  case ChiralOrientation::NONE:
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      ncen += (chi_ptr[i] != 0);
    }
    result.resize(ncen);
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      if (chi_ptr[i] != 0) {
        result[ncen] = abs(chi_ptr[i]) - 1;
        ncen++;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int4> ChemicalFeatures::getChiralArmBaseAtoms() const {
  return chiral_arm_atoms.readHost();
}

//-------------------------------------------------------------------------------------------------
ChiralOrientation ChemicalFeatures::getAtomChirality(const int atom_index) const {
  if (atom_index < 0 || atom_index >= atom_count) {
    rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a system containing " +
          std::to_string(atom_count) + " atoms.", "ChemicalFeatures", "getAtomChirality");
  }
  const int irec = atom_index + 1;
  const int* chir_ptr = chiral_centers.data();
  for (int i = 0; i < chiral_center_count; i++) {
    if (chir_ptr[i] == irec) {
      return ChiralOrientation::SINISTER;
    }
    else if (chir_ptr[i] == -irec) {
      return ChiralOrientation::RECTUS;
    }
  }
  return ChiralOrientation::NONE;
}

//-------------------------------------------------------------------------------------------------
std::vector<ChiralOrientation> ChemicalFeatures::getAtomChirality(const int low_index,
                                                                  const int high_index) const {
  if (low_index < 0 || low_index >= high_index || high_index > atom_count) {
    rtErr("Atom index range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " is invalid for a system containing " + std::to_string(atom_count) + " atoms.",
          "ChemicalFeatures", "getAtomChirality");
  }
  std::vector<ChiralOrientation> result(high_index - low_index, ChiralOrientation::NONE);
  const int* chir_ptr = chiral_centers.data();
  for (int i = 0; i < chiral_center_count; i++) {
    const int irec = abs(chir_ptr[i]) - 1;
    if (irec >= low_index && irec < high_index) {
      if (chir_ptr[i] > 0) {
        result[irec - low_index] = ChiralOrientation::SINISTER;
      }
      else {
        result[irec - low_index] = ChiralOrientation::RECTUS;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ChiralOrientation> ChemicalFeatures::getAtomChirality() const {
  return getAtomChirality(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getChiralityMask(const ChiralOrientation direction) const {
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (atom_count + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);
  
  // Loop over all known rings in the system and select those meeting the size criterion
  const int* chi_ptr = chiral_centers.data();
  for (int i = 0; i < chiral_center_count; i++) {
    const int atom_index = abs(chi_ptr[i]) - 1;
    const int mask_idx = atom_index / n_bits;
    const int bit_idx  = atom_index - (mask_idx * n_bits);
    switch (direction) {
    case ChiralOrientation::RECTUS:
      if (chi_ptr[i] < 0) {
        result[mask_idx] |= (0x1 << bit_idx);
      }
      break;
    case ChiralOrientation::SINISTER:
      if (chi_ptr[i] > 0) {
        result[mask_idx] |= (0x1 << bit_idx);
      }
      break;
    case ChiralOrientation::NONE:      

      // Here, this indicates no preference for chirality.  The list still contains a distinct
      // subset of the aotms displaying chirality, not atoms that are achiral.
      result[mask_idx] |= (0x1 << bit_idx);
      break;      
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ChemicalFeatures::getFormalCharges() const {
  return formal_charges.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ChemicalFeatures::getBondOrders() const {
  return bond_orders.readHost(0, ag_pointer->getBondTermCount());
}

//-------------------------------------------------------------------------------------------------
double ChemicalFeatures::getFreeElectrons(const int atom_index) const {
  return free_electrons.readHost(atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ChemicalFeatures::getFreeElectrons(const int low_index,
                                                       const int high_index) const {
  return free_electrons.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ChemicalFeatures::getFreeElectrons() const {
  return free_electrons.readHost();
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& ChemicalFeatures::getZeroKelvinFormalCharges() const {
  return zerok_formal_charges;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& ChemicalFeatures::getZeroKelvinBondOrders() const {
  return zerok_bond_orders;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& ChemicalFeatures::getZeroKelvinFreeElectrons() const {
  return zerok_free_electrons;
}

//-------------------------------------------------------------------------------------------------
ullint ChemicalFeatures::getRingInclusion(const int atom_index) const {
  return ring_inclusion.readHost(atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<ullint> ChemicalFeatures::getRingInclusion(const int low_index,
                                                       const int high_index) const {
  return ring_inclusion.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<ullint> ChemicalFeatures::getRingInclusion() const {
  return ring_inclusion.readHost();
}

//-------------------------------------------------------------------------------------------------
bool ChemicalFeatures::bondIsInRing(const int atom_i, const int atom_j) const {

  // Check that the two atoms are valid and form a bond
  if (atom_i < 0 || atom_i >= atom_count || atom_j < 0 || atom_j >= atom_count) {
    rtErr("Atom indices " + std::to_string(atom_i) + " and " + std::to_string(atom_j) + " are "
          "invalid for a topology with " + std::to_string(atom_count) + " atoms.",
          "ChemicalFeatures", "bondIsInRing");
  }
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const int ai_lim = vk.bond_asgn_bounds[atom_i + 1];
  const int aj_lim = vk.bond_asgn_bounds[atom_j + 1];
  for (int i = vk.bond_asgn_bounds[atom_i]; i < ai_lim; i++) {
    if (vk.bond_asgn_atoms[i] == atom_j) {
      return bond_in_ring[vk.bond_asgn_terms[i]];
    }
  }
  for (int i = vk.bond_asgn_bounds[atom_j]; i < aj_lim; i++) {
    if (vk.bond_asgn_atoms[i] == atom_i) {
      return bond_in_ring[vk.bond_asgn_terms[i]];
    }
  }
  rtErr("No bond between atoms " + char4ToString(ag_pointer->getAtomName(atom_i)) + " (index " +
        std::to_string(atom_i) + ") and " + char4ToString(ag_pointer->getAtomName(atom_j)) +
        " (index " + std::to_string(atom_j) + ") exists in topology " +
        getBaseName(ag_pointer->getFileName()) + ".", "ChemicalFeatures", "bondIsInRing");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ChemicalFeatures::bondIsInRing(const int bond_index) const {
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  if (bond_index < 0 || bond_index >= vk.nbond) {
    rtErr("Bond index " + std::to_string(bond_index) + " is invalid for a topology of " +
          std::to_string(vk.nbond) + " bonds.", "ChemicalFeatures", "bondIsInRing");
  }
  return bond_in_ring[bond_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<IsomerPlan> ChemicalFeatures::getRotatableBondGroups() const {
  std::vector<IsomerPlan> result;
  result.reserve(rotatable_bond_count);
  const int* rg_ptr = rotatable_groups.data();
  for (int i = 0; i < rotatable_bond_count; i++) {
    const int llim = rotatable_group_bounds.readHost(i);
    const int hlim = rotatable_group_bounds.readHost(i + 1);
    std::vector<int> moving_atom_idx(hlim - llim - 2);
    int k = 0;
    for (int j = llim + 2; j < hlim; j++) {
      moving_atom_idx[k] = rg_ptr[j];
      k++;
    }
    result.emplace_back(ConformationEdit::BOND_ROTATION, rg_ptr[llim], rg_ptr[llim + 1],
                        moving_atom_idx, ag_pointer);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<IsomerPlan>
ChemicalFeatures::getRotatableBondGroups(const int cutoff, const int mol_index) const {

  // Collect all rotatable groups on a particular molecule larger than a stated cutoff size.
  // Order the results in descending order of the number of atoms that rotate as a consequence
  // of twisting the bond.
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  if (mol_index < 0 || mol_index >= cdk.nmol) {
    rtErr("Molecule index " + std::to_string(mol_index) + " is invalid for a system with " +
          std::to_string(cdk.nmol) + " molecules.", "ChemicalFeatures", "getRotatableBondGroups");
  }
  int n_rgroup = 0;
  for (int i = 0; i < rotatable_bond_count; i++) {
    const int llim = rotatable_group_bounds.readHost(i);
    const int hlim = rotatable_group_bounds.readHost(i + 1);
    n_rgroup += (cdk.mol_home[llim + 1] == mol_index && hlim - llim - 2 >= cutoff);
  }
  std::vector<IsomerPlan> result;
  result.reserve(n_rgroup);
  const int* rg_ptr = rotatable_groups.data();
  n_rgroup = 0;
  for (int i = 0; i < rotatable_bond_count; i++) {
    const int llim = rotatable_group_bounds.readHost(i);
    const int hlim = rotatable_group_bounds.readHost(i + 1);
    if (cdk.mol_home[llim + 1] == mol_index && hlim - llim - 2 >= cutoff) {
      std::vector<int> moving_atom_idx(hlim - llim - 2);
      int k = 0;
      for (int j = llim + 2; j < hlim; j++) {
        moving_atom_idx[k] = rg_ptr[j];
        k++;
      }
      result.emplace_back(ConformationEdit::BOND_ROTATION, rg_ptr[llim], rg_ptr[llim + 1],
                          moving_atom_idx, ag_pointer);
      n_rgroup++;
    }
  }
  std::sort(result.begin(), result.end(),
            [](IsomerPlan a, IsomerPlan b) {
              return a.getMovingAtomCount() > b.getMovingAtomCount();
            });
  
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<IsomerPlan> ChemicalFeatures::getCisTransIsomerizationGroups() const {
  std::vector<IsomerPlan> result;
  result.reserve(cis_trans_bond_count);
  const int* rg_ptr = cis_trans_groups.data();
  for (int i = 0; i < cis_trans_bond_count; i++) {
    const int llim = cis_trans_group_bounds.readHost(i);
    const int hlim = cis_trans_group_bounds.readHost(i + 1);
    std::vector<int> moving_atom_idx(hlim - llim - 2);
    int k = 0;
    for (int j = llim + 2; j < hlim; j++) {
      moving_atom_idx[k] = rg_ptr[j];
      k++;
    }
    result.emplace_back(ConformationEdit::CIS_TRANS_FLIP, rg_ptr[llim], rg_ptr[llim + 1],
                        moving_atom_idx, ag_pointer);    
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<IsomerPlan> ChemicalFeatures::getChiralInversionGroups() const {
  std::vector<IsomerPlan> result;
  result.reserve(chiral_center_count);
  const int* groups_ptr   = invertible_groups.data();
  const int* bounds_ptr   = invertible_group_bounds.data();
  const int* anchor_a_ptr = anchor_a_branches.data();
  const int* anchor_b_ptr = anchor_b_branches.data();
  for (int i = 0; i < chiral_center_count; i++) {
    std::vector<int> moving_atom_idx(bounds_ptr[i + 1] - bounds_ptr[i]);
    int k = 0;
    for (int j = bounds_ptr[i]; j < bounds_ptr[i + 1]; j++) {
      moving_atom_idx[k] = groups_ptr[j];
      k++;
    }
    result.emplace_back(ConformationEdit::CHIRAL_INVERSION,
                        static_cast<ChiralInversionProtocol>(chiral_inversion_methods.readHost(i)),
                        anchor_a_ptr[i], anchor_b_ptr[i], moving_atom_idx, ag_pointer);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ChiralInversionProtocol> ChemicalFeatures::getChiralInversionMethods() const {
  std::vector<ChiralInversionProtocol> result(chiral_center_count);
  for (int i = 0; i < chiral_center_count; i++) {
    result[i] = static_cast<ChiralInversionProtocol>(chiral_inversion_methods.readHost(i));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ChiralInversionProtocol ChemicalFeatures::getChiralInversionMethods(const int index) const {
  return static_cast<ChiralInversionProtocol>(chiral_inversion_methods.readHost(index));
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ChemicalFeatures::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const ChemicalFeatures* ChemicalFeatures::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
ChemicalFeaturesReader ChemicalFeatures::data(const HybridTargetLevel tier) const {
  return ChemicalFeaturesReader(atom_count, planar_atom_count, ring_count, fused_ring_count,
                                twistable_ring_count, conjugated_group_count, aromatic_group_count,
                                polar_hydrogen_count, hbond_donor_count, hbond_acceptor_count,
                                chiral_center_count, rotatable_bond_count, cis_trans_bond_count,
                                double_bond_count, triple_bond_count, max_ring_size, temperature,
                                rotating_groups_mapped, chiralities_computed,
                                planar_centers.data(tier), ring_atoms.data(tier),
                                ring_atom_bounds.data(tier), aromatic_pi_electrons.data(tier),
                                aromatic_groups.data(tier), aromatic_group_bounds.data(tier),
                                polar_hydrogens.data(tier), hydrogen_bond_donors.data(tier),
                                hydrogen_bond_acceptors.data(tier), chiral_centers.data(tier),
                                chiral_inversion_methods.data(tier), rotatable_groups.data(tier),
                                rotatable_group_bounds.data(tier), cis_trans_groups.data(tier),
                                cis_trans_group_bounds.data(tier), invertible_groups.data(tier),
                                invertible_group_bounds.data(tier), chiral_arm_atoms.data(tier),
                                formal_charges.data(tier), bond_orders.data(tier),
                                free_electrons.data(tier), ring_inclusion.data(tier), ag_pointer);
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findChiralOrientations(const CoordinateFrameReader &cfr) {
  int chiral_center_timing;
  if (timer != nullptr) {
    chiral_center_timing = timer->addCategory("[ChemicalFeatures] Chiral center detection");
  }
  const int4* arm_ptr = chiral_arm_atoms.data();
  for (int i = 0; i < chiral_center_count; i++) {
    const int chatom = abs(chiral_centers.readHost(i)) - 1;
    const int chiral_direction = getChiralOrientation(cfr, chatom, arm_ptr[i].x, arm_ptr[i].y,
                                                      arm_ptr[i].z, arm_ptr[i].w);

    // Combine the chirality determination with the atom index (offset by one to prevent the
    // abiguity of zero and -zero), then contribute the result.
    chiral_centers.putHost(chiral_direction * (chatom + 1), i);
  }
  if (timer != nullptr) timer->assignTime(chiral_center_timing);
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findRotatableBondGroups(StopWatch *timer) {
  std::vector<int> tmp_rotatable_groups, tmp_rotatable_group_bounds;
  std::vector<int> tmp_cis_trans_groups, tmp_cis_trans_group_bounds;
  std::vector<int> tmp_invertible_groups, tmp_invertible_group_bounds;
  std::vector<int> tmp_anchor_a_branches, tmp_anchor_b_branches;
  
  // Obtain abstracts from the topology
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  int rotatable_bond_timing;
  if (timer != nullptr) {
    rotatable_bond_timing = timer->addCategory("[ChemicalFeatures] Rotatable bond mapping");
  }
  findRotatableBonds(vk, cdk, nbk, ring_atoms.readHost(), ring_atom_bounds.readHost(),
                     &tmp_rotatable_groups, &tmp_rotatable_group_bounds, &tmp_cis_trans_groups,
                     &tmp_cis_trans_group_bounds);
  rotatable_bond_count = static_cast<int>(tmp_rotatable_group_bounds.size()) - 1;
  cis_trans_bond_count = static_cast<int>(tmp_cis_trans_group_bounds.size()) - 1;
  findInvertibleGroups(chiral_centers.readHost(), chiral_inversion_methods.readHost(),
                       &tmp_anchor_a_branches, &tmp_anchor_b_branches, &tmp_invertible_groups,
                       &tmp_invertible_group_bounds);
  if (timer != nullptr) timer->assignTime(rotatable_bond_timing);
  allocateMutableData(tmp_rotatable_groups, tmp_rotatable_group_bounds, tmp_cis_trans_groups,
                      tmp_cis_trans_group_bounds, tmp_invertible_groups,
                      tmp_invertible_group_bounds, tmp_anchor_a_branches, tmp_anchor_b_branches);
  if (timer != nullptr) timer->assignTime(0);
}

//-------------------------------------------------------------------------------------------------
void unpackRotatableGroups(const std::vector<int2> &bond_endpoints,
                           const std::vector<std::vector<int>> &moving_lists,
                           std::vector<int> *tmp_groups, std::vector<int> *tmp_group_bounds) {
  const size_t ngroup = bond_endpoints.size();
  tmp_group_bounds->resize(ngroup + 1);
  int* gbounds_ptr = tmp_group_bounds->data();
  int acc_size = 0;
  for (size_t i = 0; i < ngroup; i++) {
    gbounds_ptr[i] = acc_size;
    acc_size += 2 + moving_lists[i].size();
  }
  gbounds_ptr[ngroup] = acc_size;
  tmp_groups->resize(acc_size);
  int* grp_ptr = tmp_groups->data();
  size_t k = 0;
  for (size_t i = 0; i < ngroup; i++) {
    grp_ptr[k] = bond_endpoints[i].x;
    k++;
    grp_ptr[k] = bond_endpoints[i].y;
    k++;
    const size_t gsize = moving_lists[i].size();
    for (size_t j = 0; j < gsize; j++) {
      grp_ptr[k] = moving_lists[i][j];
      k++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
bool scoreChiralBranches(const std::vector<std::vector<BondedNode>> &links,
                         const std::vector<int> &layer_llim, const std::vector<int> &layer_hlim,
                         const ChemicalDetailsKit &cdk, std::vector<int> *chiral_dominance,
                         std::vector<int> *parallel_growth) {
  std::vector<bool> score_branch(4);
  int* cdom_ptr = chiral_dominance->data();
  score_branch[0] = (cdom_ptr[ 1] == 0 || cdom_ptr[ 2] == 0 || cdom_ptr[ 3] == 0);
  score_branch[1] = (cdom_ptr[ 1] == 0 || cdom_ptr[ 6] == 0 || cdom_ptr[ 7] == 0);
  score_branch[2] = (cdom_ptr[ 2] == 0 || cdom_ptr[ 6] == 0 || cdom_ptr[11] == 0);
  score_branch[3] = (cdom_ptr[ 3] == 0 || cdom_ptr[ 7] == 0 || cdom_ptr[11] == 0);

  // Determine if there is any work to do.  If there are no new layers to offer new scores,
  // return false to indicate that branches no longer need advancement.
  int max_leaves = 0;
  for (int i = 0; i < 4; i++) {
    if (score_branch[i]) {
      max_leaves = std::max(layer_hlim[i] - layer_llim[i], max_leaves);
    }
  }
  if (max_leaves == 0) {
    return false;
  }

  // Check for parallel growth.  Two consecutive rounds of branches growing with the same atoms
  // implies that neither branch is superior to the other and thus the center is not chiral.
  int* pgrow_ptr = parallel_growth->data();
  for (int i = 0; i < 3; i++) {
    const int lyl_i = layer_llim[i];
    const int lyh_i = layer_hlim[i];
    for (int j = i + 1; j < 4; j++) {
      if (layer_llim[j] == lyl_i && layer_hlim[j] == lyh_i) {
        bool leaves_identical = true;
        for (int k = lyl_i; k < lyh_i; k++) {
          leaves_identical = (leaves_identical && links[i][k].getAtom() == links[j][k].getAtom());
        }
        if (leaves_identical) {
          pgrow_ptr[(4 * j) + i] += 1;
          pgrow_ptr[(4 * i) + j] += 1;
        }
        else {
          pgrow_ptr[(4 * j) + i] = 0;
          pgrow_ptr[(4 * i) + j] = 0;
        }
      }
    }
  }

  // Set up tables to score the various branches
  int max_znum = 0;
  for (int i = 0; i < 4; i++) {
    if (score_branch[i]) {
      for (int j = layer_llim[i]; j < layer_hlim[i]; j++) {
        max_znum = std::max(cdk.z_numbers[links[i][j].getAtom()], max_znum);
      }
    }
  }
  std::vector<std::vector<double>> element_count(4, std::vector<double>(max_znum + 1, 0.0));
  for (int i = 0; i < 4; i++) {
    if (score_branch[i]) {
      for (int j = layer_llim[i]; j < layer_hlim[i]; j++) {
        element_count[i][cdk.z_numbers[links[i][j].getAtom()]] += links[i][j].getRootBondOrder();
      }
    }
  }
  for (int i = 0; i < 3; i++) {
    if (score_branch[i]) {
      for (int j = i + 1; j < 4; j++) {
        const int j4pi = (j * 4) + i;
        if (score_branch[j] && cdom_ptr[j4pi] == 0 && pgrow_ptr[j4pi] < 2) {
          const int i4pj = (i * 4) + j;
          for (int k = max_znum; k >= 1; k--) {
            if (cdom_ptr[j4pi] == 0) {
              if (element_count[i][k] > element_count[j][k] + constants::tiny) {
                cdom_ptr[j4pi] =  1;
                cdom_ptr[i4pj] = -1;
              }
              else if (element_count[i][k] < element_count[j][k] - constants::tiny) {
                cdom_ptr[j4pi] = -1;
                cdom_ptr[i4pj] =  1;
              }
            }
          }
        }
      }
    }
  }

  // If the dominance matrix still contains undecided elements, the search must keep going.
  return ((cdom_ptr[ 1] == 0 && pgrow_ptr[ 1] < 2) || (cdom_ptr[ 2] == 0 && pgrow_ptr[ 2] < 2) ||
          (cdom_ptr[ 3] == 0 && pgrow_ptr[ 3] < 2) || (cdom_ptr[ 6] == 0 && pgrow_ptr[ 6] < 2) ||
          (cdom_ptr[ 7] == 0 && pgrow_ptr[ 7] < 2) || (cdom_ptr[11] == 0 && pgrow_ptr[11] < 2));
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const CoordinateFrameReader &cfr, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(cfr.xcrd, cfr.ycrd, cfr.zcrd, center_atom, root_atom, branch_a_atom,
                              branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const CoordinateFrame *cf, const int center_atom, const int root_atom,
                         const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(cf->data(), center_atom, root_atom, branch_a_atom, branch_b_atom,
                              branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const CoordinateFrame &cf, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(cf.data(), center_atom, root_atom, branch_a_atom, branch_b_atom,
                              branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const PhaseSpaceReader &psr, const int center_atom, const int root_atom,
                         const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(psr.xcrd, psr.ycrd, psr.zcrd, center_atom, root_atom, branch_a_atom,
                              branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const PhaseSpace *ps, const int center_atom, const int root_atom,
                         const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(ps->data(), center_atom, root_atom, branch_a_atom, branch_b_atom,
                              branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const PhaseSpace &ps, const int center_atom, const int root_atom,
                         const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(ps.data(), center_atom, root_atom, branch_a_atom, branch_b_atom,
                              branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const CondensateReader &cdnsr, const int system_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  const size_t offset = cdnsr.atom_starts[system_index];
  switch (cdnsr.mode) {
  case PrecisionModel::DOUBLE:
    return getChiralOrientation(&cdnsr.xcrd[offset], &cdnsr.ycrd[offset], &cdnsr.zcrd[offset],
                                center_atom, root_atom, branch_a_atom, branch_b_atom,
                                branch_c_atom);
  case PrecisionModel::SINGLE:
    return getChiralOrientation(&cdnsr.xcrd_sp[offset], &cdnsr.ycrd_sp[offset],
                                &cdnsr.zcrd_sp[offset], center_atom, root_atom, branch_a_atom,
                                branch_b_atom, branch_c_atom);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const Condensate *cdns, const int system_index, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(cdns->data(), system_index, center_atom, root_atom, branch_a_atom,
                              branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const Condensate &cdns, const int system_index, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {
  return getChiralOrientation(cdns.data(), system_index, center_atom, root_atom, branch_a_atom,
                              branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const PsSynthesisReader &poly_psr, const int system_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  const size_t offset = poly_psr.atom_starts[system_index];
  if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
    std::vector<double> tmp_xcrd(5), tmp_ycrd(5), tmp_zcrd(5);
    std::vector<int> markers = { center_atom, root_atom, branch_a_atom, branch_b_atom,
      branch_c_atom };
    for (int i = 0; i < 5; i++) {
      const size_t im = markers[i] + offset;
      tmp_xcrd[i] = hostInt95ToDouble(poly_psr.xcrd[im], poly_psr.xcrd_ovrf[im]) *
                    poly_psr.inv_gpos_scale;
      tmp_ycrd[i] = hostInt95ToDouble(poly_psr.ycrd[im], poly_psr.ycrd_ovrf[im]) *
                    poly_psr.inv_gpos_scale;
      tmp_zcrd[i] = hostInt95ToDouble(poly_psr.zcrd[im], poly_psr.zcrd_ovrf[im]) *
                    poly_psr.inv_gpos_scale;
    }
    return getChiralOrientation<double>(tmp_xcrd.data(), tmp_ycrd.data(), tmp_zcrd.data(),
                                        0, 1, 2, 3, 4);
  }
  else {
    return getChiralOrientation<llint>(&poly_psr.xcrd[offset], &poly_psr.ycrd[offset],
                                       &poly_psr.zcrd[offset], nullptr, nullptr, nullptr,
                                       center_atom, root_atom, branch_a_atom, branch_b_atom,
                                       branch_c_atom, poly_psr.inv_gpos_scale);
  }
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const PhaseSpaceSynthesis *poly_ps, const int system_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  return getChiralOrientation(poly_ps->data(), system_index, center_atom, root_atom, branch_a_atom,
                              branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const PhaseSpaceSynthesis &poly_ps, const int system_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  return getChiralOrientation(poly_ps.data(), system_index, center_atom, root_atom, branch_a_atom,
                              branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph &ag, const ChemicalFeatures &chemfe, const int atom_a,
                         const int atom_b) {
  return matchBondingPattern(ag, chemfe.getFormalCharges(), chemfe.getFreeElectrons(),
                             chemfe.getRingInclusion(), chemfe.getAtomChirality(), atom_a, atom_b);
}

//-------------------------------------------------------------------------------------------------
bool permutationsAreLinked(const std::vector<IsomerPlan> &isomerizers, const int permi,
                           const int permj, const NonbondedKit<double> &nbk) {
  const int root_i = isomerizers[permi].getRootAtom();
  const int pivt_i = isomerizers[permi].getPivotAtom();
  const int root_j = isomerizers[permj].getRootAtom();
  const int pivt_j = isomerizers[permj].getPivotAtom();
  switch (isomerizers[permi].getMotion()) {
  case ConformationEdit::BOND_ROTATION:
  case ConformationEdit::CIS_TRANS_FLIP:

    // If two isomerizations share atoms or have their root and pivot atoms bonded to one
    // another, consider them coupled and add the product of the number of possible states for
    // each to the running sum.
    switch (isomerizers[permj].getMotion()) {
    case ConformationEdit::BOND_ROTATION:
    case ConformationEdit::CIS_TRANS_FLIP:
      return (root_i == root_j || root_i == pivt_j || pivt_i == root_j || pivt_i == pivt_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j) ||
              isBonded(nbk, pivt_i, root_j) || isBonded(nbk, pivt_i, pivt_j));
    case ConformationEdit::CHIRAL_INVERSION:
      return (root_i == root_j || pivt_i == root_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j));
    }
    break;
  case ConformationEdit::CHIRAL_INVERSION:
    switch (isomerizers[permj].getMotion()) {
    case ConformationEdit::BOND_ROTATION:
    case ConformationEdit::CIS_TRANS_FLIP:
      return (root_i == root_j || root_i == pivt_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j));
    case ConformationEdit::CHIRAL_INVERSION:
      return (root_i == root_j || isBonded(nbk, root_i, root_j));
    }
    break;
  }
  __builtin_unreachable();
}

} // namespace chemistry
} // namespace stormm
