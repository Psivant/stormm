#include "copyright.h"
#include "Constants/generalized_born.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "UnitTesting/stopwatch.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Topology/atomgraph_refinement.h"
#include "atomgraph_synthesis.h"
#include "nonbonded_workunit.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using diskutil::getBaseName;
using parse::realToString;
using parse::NumberFormat;
using restraints::RestraintStage;
using stmath::buildReductionWorkUnits;
using stmath::incrementingSeries;
using stmath::maxAbsoluteDifference;
using stmath::maxValue;
using stmath::minValue;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::RdwuAbstractMap;
using stmath::rdwu_abstract_length;
using stmath::roundUp;
using stmath::sum;
using topology::CmapAccessories;
using topology::computeCmapDerivatives;
using topology::ConstraintKit;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using testing::Approx;
using testing::ComparisonType;
using testing::StopWatch;
using namespace generalized_born_defaults;
  
//-------------------------------------------------------------------------------------------------
std::vector<int>
AtomGraphSynthesis::createDummyRestraints(const std::vector<int> &restraint_indices_in,
                                          const std::vector<int> &topology_indices_in) {
  const int ntop = topologies.size();
  restraint_dummies.reserve(ntop);
  for (int i = 0; i < ntop; i++) {
    restraint_dummies.emplace_back(topologies[i]);
  }
  const int n_orig_networks = restraint_networks.size();
  for (int i = 0; i < ntop; i++) {
    restraint_networks.push_back(&restraint_dummies[i]);
  }
  const int nsys = restraint_indices_in.size();
  std::vector<int> new_restraint_indices(nsys);
  for (int i = 0; i < nsys; i++) {
    if (restraint_networks[restraint_indices_in[i]] == nullptr) {
      new_restraint_indices[i] = n_orig_networks + topology_indices_in[i];
    }
    else {
      new_restraint_indices[i] = restraint_indices_in[i];
    }
  }
  restraint_network_count = restraint_networks.size();
  return new_restraint_indices;
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
AtomGraphSynthesis::checkTopologyList(const std::vector<int> &topology_indices_in) {

  // Check that no system indexes a topology outside of the given list
  if (maxValue(topology_indices_in) >= topology_count || minValue(topology_indices_in) < 0) {
    rtErr("One or more systems references a topology index outside of the list supplied.",
          "AtomGraphSynthesis", "checkTopologyList");
  }
  
  // Check that each topology in the supplied list is being used.  Roll that result into an
  // analysis of the uniqueness of each topology.
  std::vector<bool> topology_unique(topology_count, false);
  for (int i = 0; i < system_count; i++) {
    topology_unique[topology_indices_in[i]] = true;
  }

  // Check that all topologies are, in fact, unique.  Compact the list if necessary and update
  // the system topology indexing.  Rebase the list of topology indices.
  std::vector<int> topology_index_rebase(topology_count, -1);
  int n_unique_top = 0;
  for (int i = 0; i < topology_count; i++) {
    if (topology_unique[i]) {
      const int topi_natom = topologies[i]->getAtomCount();
      for (int j = i + 1; j < topology_count; j++) {
        if (topologies[j] == topologies[i] ||
            (topologies[j]->getAtomCount() == topi_natom &&
             topologies[i]->getFileName() == topologies[j]->getFileName())) {
          topology_unique[j] = false;
          topology_index_rebase[j] = n_unique_top;
        }
      }
      topology_index_rebase[i] = n_unique_top;
      n_unique_top++;
    }
  }

  // Warn if there are unused topologies
  int n_unused_topology = 0;
  for (int i = 0; i < topology_count; i++) {
    n_unused_topology += (topology_unique[i] == false);
  }
  if (n_unused_topology > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Out of " + std::to_string(topology_count) + " topologies, " +
             std::to_string(n_unused_topology) + " are not referenced by any systems in this "
             "synthesis.", "AtomGraphSynthesis", "checkTopologyList");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Condense the list of unique topologies
  n_unique_top = 0;
  for (int i = 0; i < topology_count; i++) {
    if (topology_unique[i]) {
      topologies[n_unique_top] = topologies[i];
      n_unique_top++;
    }
  }
  topology_count = n_unique_top;
  if (topology_count == 0) {
    rtErr("No topologies were detected to describe " + std::to_string(system_count) + " systems.",
          "AtomGraphSynthesis", "checkTopologyList");
  }
  if (system_count == 0) {
    rtErr("No systems making use of any of the " + std::to_string(topology_count) +
          " topologies were detected.", "AtomGraphSynthesis", "checkTopologyList");
  }
  return topology_index_rebase;
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
AtomGraphSynthesis::checkRestraintList(const std::vector<int> &restraint_indices_in,
                                       const std::vector<int> &topology_indices_in,
                                       const std::vector<int> &topology_index_rebase) {
  
  // For each system in the synthesis, there is a restraint group and a topology.  Each restraint
  // group references a topology, to ensure consistency in future applications.  Check that the
  // restraint apparatus for each system references the same topology as the system itself.
  // A reference to a restraint apparatus index less than zero implies that the system uses no
  // restraints.
  if (static_cast<int>(restraint_indices_in.size()) != system_count) {
    rtErr("A restraint network index must be provided for each system that the synthesis is "
          "designed to describe, even if the restraint apparatus pointer it indicates is a "
          "null pointer.  The number of systems in the synthesis is determined by the length of "
          "the topology index array, so there must be as many indices into a list of restraint "
          "networks as there are indices into distinct system topologies.  " +
          std::to_string(static_cast<int>(restraint_indices_in.size())) + " restraint network "
          "pointer indices were provided for a synthesis that should contain " +
          std::to_string(system_count) + " systems.", "AtomGraphSynthesis", "checkRestraintList");
  }
  if (maxValue(restraint_indices_in) >= restraint_network_count) {
    rtErr("One or more systems references a restraint apparatus index outside of the list "
          "supplied.", "AtomGraphSynthesis", "checkRestraintList");
  }
  for (int i = 0; i < system_count; i++) {
    if (restraint_indices_in[i] < 0 || restraint_networks[restraint_indices_in[i]] == nullptr) {
      continue;
    }
    const AtomGraph *topref = topologies[topology_index_rebase[topology_indices_in[i]]];    
    const AtomGraph *rstref = restraint_networks[restraint_indices_in[i]]->getTopologyPointer();
    if (rstref != topref && rstref->getFileName() != topref->getFileName()) {
      rtErr("Mismatch in topologies referenced by the restraint apparatus for system " +
            std::to_string(i) + " and the system itself.  Atom counts of the topologies are " +
            std::to_string(rstref->getAtomCount()) + " and " +
            std::to_string(topref->getAtomCount()) + ".  Topology names are " +
            rstref->getFileName() + " and " + topref->getFileName() + ".", "AtomGraphSynthesis",
            "checkRestraintList");
    }
  }

  // Create a list of unique restraint apparatuses in use by the AtomGraphSynthesis, analogous
  // to the list created for topologies.
  std::vector<bool> network_unique(restraint_network_count, false);
  for (int i = 0; i < system_count; i++) {
    network_unique[restraint_indices_in[i]] = true;
  }
  for (int i = 0; i < restraint_network_count; i++) {
    if (network_unique[i]) {
      for (int j = i + 1; j < restraint_network_count; j++) {
        if (restraint_networks[j] == restraint_networks[i]) {
          network_unique[j] = false;
          continue;
        }
      }
    }
  }
  
  // Warn if there are unused restraint apparatuses
  int n_unused_network = 0;
  for (int i = 0; i < restraint_network_count; i++) {
    n_unused_network += (network_unique[i] == false && restraint_networks[i] != nullptr &&
                         restraint_networks[i]->getTotalRestraintCount() > 0);
  }
  if (n_unused_network > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Out of " + std::to_string(restraint_network_count) + " networks, " +
             std::to_string(n_unused_network) + " are not referenced by any systems in this "
             "synthesis.", "AtomGraphSynthesis", "checkRestraintList");
      break;



    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Rebase the list of restraint network indices
  std::vector<int> restraint_index_rebase(restraint_network_count, -1);
  int n_unique_network = 0;
  for (int i = 0; i < restraint_network_count; i++) {
    if (network_unique[i]) {
      restraint_index_rebase[i] = n_unique_network;
      for (int j = i + 1; j < restraint_network_count; j++) {
        if (restraint_networks[j] == restraint_networks[i]) {
          restraint_index_rebase[j] = n_unique_network;
        }
      }
      restraint_networks[n_unique_network] = restraint_networks[i];
      n_unique_network++;
    }
  }
  restraint_network_count = n_unique_network;
  if (restraint_network_count == 0) {
    rtErr("No restraint apparatuses were detected to describe " + std::to_string(system_count) +
          " systems.", "AtomGraphSynthesis", "checkRestraintList");
  }
  return restraint_index_rebase;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::checkCommonSettings() {

  // Check that all topologies contain compatible boundary conditions and solvent models
  periodic_box_class = topologies[0]->getUnitCellType();
  gb_style = topologies[0]->getImplicitSolventModel();
  dielectric_constant = topologies[0]->getDielectricConstant();
  salt_concentration = topologies[0]->getSaltConcentration();
  coulomb_constant = topologies[0]->getCoulombConstant();
  const Approx appr_dielcon(dielectric_constant, constants::tiny);
  const Approx appr_saltcon(salt_concentration, constants::tiny);
  const Approx appr_coulomb(coulomb_constant, constants::tiny);
  bool ism_problem  = false;
  bool diel_problem = false;
  bool salt_problem = false;
  bool coul_problem = false;
  bool box_problem = false;
  for (int i = 1; i < topology_count; i++) {
    ism_problem  = (ism_problem  || (topologies[i]->getImplicitSolventModel() != gb_style));
    diel_problem = (diel_problem ||
                    appr_dielcon.test(topologies[i]->getDielectricConstant()) == false);
    salt_problem = (salt_problem ||
                    appr_saltcon.test(topologies[i]->getSaltConcentration()) == false);
    coul_problem = (coul_problem ||
                    appr_coulomb.test(topologies[i]->getCoulombConstant()) == false);
    box_problem  = (box_problem  ||
                    (periodic_box_class == UnitCellType::NONE &&
                     (topologies[i]->getUnitCellType() == UnitCellType::ORTHORHOMBIC ||
                      topologies[i]->getUnitCellType() == UnitCellType::TRICLINIC)) ||
                    ((periodic_box_class == UnitCellType::ORTHORHOMBIC ||
                      periodic_box_class == UnitCellType::TRICLINIC) &&
                     topologies[i]->getUnitCellType() == UnitCellType::NONE));
  }
  if (ism_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All topologies must have a consistent implicit solvent model setting.  The first, "
            "set to " + getEnumerationName(gb_style) +
            ", does not agree with subsequent topologies.", "AtomGraphSynthesis",
            "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All topologies must have a consistent implicit solvent model setting.  The first, "
             "set to " + getEnumerationName(gb_style) +
             ", will be applied to all systems.", "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (diel_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All topologies must have a consistent dielectric constant, but values conflicting "
            "with the primary system's " + realToString(dielectric_constant) + " were found.",
            "AtomGraphSynthesis", "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All topologies must have a consistent dielectric constant.  A value of " +
             realToString(dielectric_constant) + ", taken from the first topology, will be "
             "applied throughout.", "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (salt_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All systems must use the same salt concentration, but values conflicting with the "
            "primary system's " + realToString(salt_concentration) + " were found.",
            "AtomGraphSynthesis", "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All systems must use the same salt concentration.  A value of " +
             realToString(salt_concentration) + ", found in the first topology, will be applied "
             "throughout.", "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (coul_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All topologies must make use of the same definitions of physical constants.  "
            "Coulomb's constant is defined as " + realToString(coulomb_constant, 6) +
            " in the first topology but differently in subsequent topologies.",
            "AtomGraphSynthesis", "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All topologies must make use of the same definitions of physical constants.  "
             "Coulomb's constant is defined as " + realToString(coulomb_constant, 6) +
             " in the first topology and this definition will be applied throughout.",
             "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (box_problem) {
    rtErr("All topologies must have consistent periodic or isolated boundary conditions.  "
          "The first topology is \"" + getEnumerationName(periodic_box_class) + "\" but "
          "subsequent topologies are different.  No remedy is available.", "AtomGraphSynthesis",
          "checkCommonSettings");
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::buildAtomAndTermArrays(const std::vector<int> &topology_indices_in,
                                                const std::vector<int> &topology_index_rebase,
                                                const std::vector<int> &new_restraint_indices,
                                                const std::vector<int> &restraint_index_rebase) {

  // Allocate memory and set POINTER-kind arrays for the small packets of data
  const int padded_system_count = roundUp(system_count, warp_size_int);
  int_system_data.resize(46 * padded_system_count);
  int pivot = 0;
  topology_indices.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  restraint_indices.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  residue_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  molecule_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  largest_residue_sizes.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  last_solute_residues.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  last_solute_atoms.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  first_solvent_molecules.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  ubrd_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cimp_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cmap_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  virtual_site_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  posn_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  lj_type_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  total_exclusion_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  rigid_water_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_constraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  degrees_of_freedom.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cnst_degrees_of_freedom.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  nonrigid_particle_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_bit_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  residue_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  molecule_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  ubrd_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cimp_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cmap_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  virtual_site_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  posn_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  sett_group_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cnst_group_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  nb_exclusion_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  lennard_jones_abc_offsets.setPointer(&int_system_data, pivot, system_count);
  
  // Load the topology indexing first
  for (int i = 0; i < system_count; i++) {
    topology_indices.putHost(topology_index_rebase[topology_indices_in[i]], i);
    restraint_indices.putHost(restraint_index_rebase[new_restraint_indices[i]], i);
  }

  // Loop over all systems, fill in the above details, and compute the sizes of various arrays
  // for topological valence terms, virtual sites, exclusions, and restraints
  int atom_offset = 0;
  int abit_offset = 0;
  int resi_offset = 0;
  int mole_offset = 0;
  int ubrd_offset = 0;
  int cimp_offset = 0;
  int cmap_offset = 0;
  int bond_offset = 0;
  int angl_offset = 0;
  int dihe_offset = 0;
  int vste_offset = 0;
  int excl_offset = 0;
  int sett_offset = 0;
  int cnst_offset = 0;
  int cnst_bounds_offset = 0;
  int cnst_atom_offset = 0;
  int rposn_offset = 0;
  int rbond_offset = 0;
  int rangl_offset = 0;
  int rdihe_offset = 0;
  for (int i = 0; i < system_count; i++) {
    const AtomGraph *ag_ptr = topologies[topology_indices.readHost(i)];
    const RestraintApparatus* ra_ptr = restraint_networks[restraint_indices.readHost(i)];
    const ChemicalDetailsKit cdk     = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    const ConstraintKit<double> cnk  = ag_ptr->getDoublePrecisionConstraintKit();
    const int nvsite = ag_ptr->getVirtualSiteCount();
    atom_counts.putHost(cdk.natom, i);
    total_atoms += cdk.natom;
    virtual_site_counts.putHost(nvsite, i);
    total_virtual_sites += nvsite;
    residue_counts.putHost(cdk.nres, i);
    molecule_counts.putHost(cdk.nmol, i);
    largest_residue_sizes.putHost(ag_ptr->getLargestResidueSize(), i);
    last_solute_residues.putHost(ag_ptr->getLastSoluteResidue(), i);
    last_solute_atoms.putHost(ag_ptr->getLastSoluteAtom(), i);
    first_solvent_molecules.putHost(ag_ptr->getFirstSolventMolecule(), i);
    ubrd_term_counts.putHost(vk.nubrd, i);
    cimp_term_counts.putHost(vk.ncimp, i);
    cmap_term_counts.putHost(vk.ncmap, i);
    bond_term_counts.putHost(vk.nbond, i);
    angl_term_counts.putHost(vk.nangl, i);
    dihe_term_counts.putHost(vk.ndihe, i);
    total_ubrd_terms += vk.nubrd;
    total_cimp_terms += vk.ncimp;
    total_cmap_terms += vk.ncmap;
    total_bond_terms += vk.nbond;
    total_angl_terms += vk.nangl;
    total_dihe_terms += vk.ndihe;
    lj_type_counts.putHost(nbk.n_lj_types, i);
    total_exclusion_counts.putHost(ag_ptr->getTotalExclusions(), i);
    rigid_water_counts.putHost(ag_ptr->getRigidWaterCount(), i);
    bond_constraint_counts.putHost(ag_ptr->getBondConstraintCount(), i);
    degrees_of_freedom.putHost(ag_ptr->getDegreesOfFreedom(), i);
    cnst_degrees_of_freedom.putHost(ag_ptr->getConstrainedDegreesOfFreedom(), i);
    nonrigid_particle_counts.putHost(ag_ptr->getNonrigidParticleCount(), i);

    // Record various offsets, and increment the counters
    atom_offsets.putHost(atom_offset, i);
    atom_offset += roundUp(cdk.natom, warp_size_int);
    const int uint_bits = sizeof(uint) * 8;
    atom_bit_offsets.putHost(abit_offset, i);
    abit_offset += roundUp((cdk.natom + uint_bits - 1) / uint_bits, warp_size_int);
    virtual_site_offsets.putHost(vste_offset, i);
    vste_offset += roundUp(nvsite, warp_size_int);
    residue_offsets.putHost(resi_offset, i);
    resi_offset += roundUp(cdk.nres, warp_size_int);
    molecule_offsets.putHost(mole_offset, i);
    mole_offset += roundUp(cdk.nmol, warp_size_int);
    ubrd_term_offsets.putHost(ubrd_offset, i);
    ubrd_offset += roundUp(vk.nubrd, warp_size_int);
    cimp_term_offsets.putHost(cimp_offset, i);
    cimp_offset += roundUp(vk.ncimp, warp_size_int);
    cmap_term_offsets.putHost(cmap_offset, i);
    cmap_offset += roundUp(vk.ncmap, warp_size_int);
    bond_term_offsets.putHost(bond_offset, i);
    bond_offset += roundUp(vk.nbond, warp_size_int);
    angl_term_offsets.putHost(angl_offset, i);
    angl_offset += roundUp(vk.nangl, warp_size_int);
    dihe_term_offsets.putHost(dihe_offset, i);
    dihe_offset += roundUp(vk.ndihe, warp_size_int);
    nb_exclusion_offsets.putHost(excl_offset, i);
    excl_offset += roundUp(ag_ptr->getTotalExclusions(), warp_size_int);
    if (ra_ptr != nullptr) {
      const RestraintKit<double, double2, double4> rar = ra_ptr->dpData();
      posn_restraint_counts.putHost(rar.nposn, i);
      bond_restraint_counts.putHost(rar.nbond, i);
      angl_restraint_counts.putHost(rar.nangl, i);
      dihe_restraint_counts.putHost(rar.ndihe, i);
      total_position_restraints += rar.nposn;
      total_distance_restraints += rar.nbond;
      total_angle_restraints    += rar.nangl;
      total_dihedral_restraints += rar.ndihe;
      posn_restraint_offsets.putHost(rposn_offset, i);
      rposn_offset += roundUp(rar.nposn, warp_size_int);
      bond_restraint_offsets.putHost(rbond_offset, i);
      rbond_offset += roundUp(rar.nbond, warp_size_int);
      angl_restraint_offsets.putHost(rangl_offset, i);
      rangl_offset += roundUp(rar.nangl, warp_size_int);
      dihe_restraint_offsets.putHost(rdihe_offset, i);
      rdihe_offset += roundUp(rar.ndihe, warp_size_int);
    }
    sett_group_offsets.putHost(sett_offset, i);
    sett_offset += roundUp(cnk.nsettle, warp_size_int);
    cnst_group_offsets.putHost(cnst_offset, i);

    // While all of the other offsets pertain to terms in each topology, the parameter sets
    // themselves are not part of the picture yet as the consensus tables have yet to be
    // determined.  However, all of the terms save for these constraint groups have a set size.
    // Additional counters are needed to track the constraint group storage sizes.
    cnst_offset += roundUp(cnk.ngroup, warp_size_int);
    cnst_bounds_offset += roundUp(cnk.ngroup + 1, warp_size_int);
    if (cnk.ngroup > 0) {
      cnst_atom_offset += roundUp(cnk.group_bounds[cnk.ngroup], warp_size_int);
    }
  }

  // Allocate detailed arrays for each descriptor, then collate all topologies.  This
  // fills the "atom and residue details" arrays listed in atomgraph_synthesis.h.
  chem_int_data.resize(9 * atom_offset);
  chem_int2_data.resize(resi_offset + mole_offset);
  chem_double_data.resize(8 * atom_offset);
  chem_float_data.resize(8 * atom_offset);
  chem_char4_data.resize(resi_offset + (2 * atom_offset));
  pivot = 0;
  residue_limits.setPointer(&chem_int2_data, pivot, resi_offset);
  pivot += resi_offset;
  molecule_limits.setPointer(&chem_int2_data, pivot, mole_offset);
  pivot = 0;
  atom_struc_numbers.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  residue_numbers.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  atomic_numbers.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  mobile_atoms.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  molecule_membership.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  molecule_contents.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  charge_indices.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  lennard_jones_indices.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  neck_gb_indices.setPointer(&chem_int_data, pivot, atom_offset);
  pivot = 0;
  atomic_charges.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  atomic_masses.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  inverse_atomic_masses.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  atomic_pb_radii.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  gb_screening_factors.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  gb_alpha_parameters.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  gb_beta_parameters.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  gb_gamma_parameters.setPointer(&chem_double_data, pivot, atom_offset);
  pivot = 0;
  sp_atomic_charges.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_atomic_masses.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_inverse_atomic_masses.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_atomic_pb_radii.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_gb_screening_factors.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_gb_alpha_parameters.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_gb_beta_parameters.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_gb_gamma_parameters.setPointer(&chem_float_data, pivot, atom_offset);
  pivot = 0;
  atom_names.setPointer(&chem_char4_data, pivot, atom_offset);
  pivot += atom_offset;
  atom_types.setPointer(&chem_char4_data, pivot, atom_offset);
  pivot += atom_offset;
  residue_names.setPointer(&chem_char4_data, pivot, resi_offset);

  // Fill the above atom and residue descriptor arrays
  int2* residue_limits_ptr     = residue_limits.data();
  int2* molecule_limits_ptr    = molecule_limits.data();
  int* atom_struc_numbers_ptr  = atom_struc_numbers.data();
  int* residue_numbers_ptr     = residue_numbers.data();
  int* atomic_numbers_ptr      = atomic_numbers.data();
  int* molecule_membership_ptr = molecule_membership.data();
  int* molecule_contents_ptr   = molecule_contents.data();
  int* lj_idx_ptr              = lennard_jones_indices.data();
  int* neck_gb_indices_ptr     = neck_gb_indices.data();
  double* atomic_charges_ptr   = atomic_charges.data();
  double* pb_radii_ptr         = atomic_pb_radii.data();
  double* gb_screen_ptr        = gb_screening_factors.data();
  double* gb_alpha_ptr         = gb_alpha_parameters.data();
  double* gb_beta_ptr          = gb_beta_parameters.data();
  double* gb_gamma_ptr         = gb_gamma_parameters.data();
  float* sp_atomic_charges_ptr = sp_atomic_charges.data();
  float* sp_pb_radii_ptr       = sp_atomic_pb_radii.data();
  float* sp_gb_screen_ptr      = sp_gb_screening_factors.data();
  float* sp_gb_alpha_ptr       = sp_gb_alpha_parameters.data();
  float* sp_gb_beta_ptr        = sp_gb_beta_parameters.data();
  float* sp_gb_gamma_ptr       = sp_gb_gamma_parameters.data();
  char4* atom_names_ptr        = atom_names.data();
  char4* atom_types_ptr        = atom_types.data();
  char4* residue_names_ptr     = residue_names.data();
  for (int i = 0; i < system_count; i++) {
    const AtomGraph *ag_ptr = topologies[topology_indices.readHost(i)];
    const ChemicalDetailsKit cdk           = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk         = ag_ptr->getDoublePrecisionNonbondedKit();
    const NonbondedKit<float> nbk_sp       = ag_ptr->getSinglePrecisionNonbondedKit();
    const ImplicitSolventKit<double> isk   = ag_ptr->getDoublePrecisionImplicitSolventKit();
    const ImplicitSolventKit<float> isk_sp = ag_ptr->getSinglePrecisionImplicitSolventKit();
    const ValenceKit<double> vk            = ag_ptr->getDoublePrecisionValenceKit();
    const int synth_bit_base = atom_bit_offsets.readHost(i);
    const int synth_atom_base = atom_offsets.readHost(i);
    const int synth_residue_base = residue_offsets.readHost(i);
    const int synth_molecule_base = molecule_offsets.readHost(i);
    for (int j = 0; j < cdk.nres; j++) {
      residue_limits_ptr[synth_residue_base + j].x = cdk.res_limits[j    ] + synth_atom_base;
      residue_limits_ptr[synth_residue_base + j].y = cdk.res_limits[j + 1] + synth_atom_base;
    }
    for (int j = 0; j < cdk.nmol; j++) {
      molecule_limits_ptr[synth_molecule_base + j].x = cdk.mol_limits[j    ] + synth_atom_base;
      molecule_limits_ptr[synth_molecule_base + j].y = cdk.mol_limits[j + 1] + synth_atom_base;
    }
    for (int j = 0; j < cdk.natom; j++) {
      atom_struc_numbers_ptr[synth_atom_base + j] = cdk.atom_numbers[j];
      residue_numbers_ptr[synth_residue_base + j] = cdk.res_numbers[j];
      atomic_numbers_ptr[synth_atom_base + j] = cdk.z_numbers[j];
      molecule_membership_ptr[synth_atom_base + j] = cdk.mol_home[j] + synth_molecule_base;
      molecule_contents_ptr[synth_atom_base + j] = cdk.mol_contents[j] + synth_atom_base;
      atomic_charges_ptr[synth_atom_base + j] = nbk.charge[j];
      sp_atomic_charges_ptr[synth_atom_base + j] = nbk_sp.charge[j];
      neck_gb_indices_ptr[synth_atom_base + j] = isk.neck_gb_idx[j];
      pb_radii_ptr[synth_atom_base + j] = isk.pb_radii[j];
      gb_screen_ptr[synth_atom_base + j] = isk.gb_screen[j];
      gb_alpha_ptr[synth_atom_base + j] = isk.gb_alpha[j];
      gb_beta_ptr[synth_atom_base + j] = isk.gb_beta[j];
      gb_gamma_ptr[synth_atom_base + j] = isk.gb_gamma[j];
      sp_pb_radii_ptr[synth_atom_base + j] = isk_sp.pb_radii[j];
      sp_gb_screen_ptr[synth_atom_base + j] = isk_sp.gb_screen[j];
      sp_gb_alpha_ptr[synth_atom_base + j] = isk_sp.gb_alpha[j];
      sp_gb_beta_ptr[synth_atom_base + j] = isk_sp.gb_beta[j];
      sp_gb_gamma_ptr[synth_atom_base + j] = isk_sp.gb_gamma[j];
      lj_idx_ptr[synth_atom_base + j] = nbk.lj_idx[j];
      atom_names_ptr[synth_atom_base + j].x = cdk.atom_names[j].x;
      atom_names_ptr[synth_atom_base + j].y = cdk.atom_names[j].y;
      atom_names_ptr[synth_atom_base + j].z = cdk.atom_names[j].z;
      atom_names_ptr[synth_atom_base + j].w = cdk.atom_names[j].w;
    }
    for (int j = 0; j < cdk.nres; j++) {
      residue_names_ptr[synth_residue_base + j].x = cdk.res_names[j].x;
      residue_names_ptr[synth_residue_base + j].y = cdk.res_names[j].y;
      residue_names_ptr[synth_residue_base + j].z = cdk.res_names[j].z;
      residue_names_ptr[synth_residue_base + j].w = cdk.res_names[j].w;
    }
    const std::vector<bool> atom_mobility = ag_ptr->getAtomMobility();
    const std::vector<int> iatom_mobility(atom_mobility.begin(), atom_mobility.end());
    mobile_atoms.putHost(iatom_mobility, synth_atom_base, cdk.natom);
    atomic_masses.putHost(ag_ptr->getAtomicMass<double>(MassForm::ORDINARY), synth_atom_base,
                          cdk.natom);
    inverse_atomic_masses.putHost(ag_ptr->getAtomicMass<double>(MassForm::INVERSE),
                                  synth_atom_base, cdk.natom);
    sp_atomic_masses.putHost(ag_ptr->getAtomicMass<float>(MassForm::ORDINARY), synth_atom_base,
                             cdk.natom);
    sp_inverse_atomic_masses.putHost(ag_ptr->getAtomicMass<float>(MassForm::INVERSE),
                                     synth_atom_base, cdk.natom);
    atom_types.putHost(ag_ptr->getAtomType().readHost(), synth_atom_base, cdk.natom);
  }

  // Fill in the valence term atom indexing arrays.
  valence_int_data.resize((3 * ubrd_offset) + (5 * cimp_offset) + (6 * cmap_offset) +
                          (3 * bond_offset) + (4 * angl_offset) + (5 * dihe_offset));
  pivot = 0;
  ubrd_i_atoms.setPointer(&valence_int_data, pivot, ubrd_offset);
  pivot += ubrd_offset;
  ubrd_k_atoms.setPointer(&valence_int_data, pivot, ubrd_offset);
  pivot += ubrd_offset;
  ubrd_param_idx.setPointer(&valence_int_data, pivot, ubrd_offset);
  pivot += ubrd_offset;
  cimp_i_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_j_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_k_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_l_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_param_idx.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cmap_i_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_j_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_k_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_l_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_m_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_param_idx.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  bond_i_atoms.setPointer(&valence_int_data, pivot, bond_offset);
  pivot += bond_offset;
  bond_j_atoms.setPointer(&valence_int_data, pivot, bond_offset);
  pivot += bond_offset;
  bond_param_idx.setPointer(&valence_int_data, pivot, bond_offset);
  pivot += bond_offset;
  angl_i_atoms.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  angl_j_atoms.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  angl_k_atoms.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  angl_param_idx.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  dihe_i_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_j_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_k_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_l_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_param_idx.setPointer(&valence_int_data, pivot, dihe_offset);
  for (int sysid = 0; sysid < system_count; sysid++) {
    const AtomGraph *ag_ptr = topologies[topology_indices.readHost(sysid)];
    const ValenceKit<double> vk = ag_ptr->getDoublePrecisionValenceKit();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_ubrd_offset = ubrd_term_offsets.readHost(sysid);
    const int synth_cimp_offset = cimp_term_offsets.readHost(sysid);
    const int synth_cmap_offset = cmap_term_offsets.readHost(sysid);
    const int synth_bond_offset = bond_term_offsets.readHost(sysid);
    const int synth_angl_offset = angl_term_offsets.readHost(sysid);
    const int synth_dihe_offset = dihe_term_offsets.readHost(sysid);
    for (int pos = 0; pos < vk.nubrd; pos++) {
      ubrd_i_atoms.putHost(vk.ubrd_i_atoms[pos] + synth_atom_base, synth_ubrd_offset + pos);
      ubrd_k_atoms.putHost(vk.ubrd_k_atoms[pos] + synth_atom_base, synth_ubrd_offset + pos);
    }
    for (int pos = 0; pos < vk.ncimp; pos++) {
      cimp_i_atoms.putHost(vk.cimp_i_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
      cimp_j_atoms.putHost(vk.cimp_j_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
      cimp_k_atoms.putHost(vk.cimp_k_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
      cimp_l_atoms.putHost(vk.cimp_l_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
    }
    for (int pos = 0; pos < vk.ncmap; pos++) {
      cmap_i_atoms.putHost(vk.cmap_i_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_j_atoms.putHost(vk.cmap_j_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_k_atoms.putHost(vk.cmap_k_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_l_atoms.putHost(vk.cmap_l_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_m_atoms.putHost(vk.cmap_m_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
    }
    for (int pos = 0; pos < vk.nbond; pos++) {
      bond_i_atoms.putHost(vk.bond_i_atoms[pos] + synth_atom_base, synth_bond_offset + pos);
      bond_j_atoms.putHost(vk.bond_j_atoms[pos] + synth_atom_base, synth_bond_offset + pos);
    }
    for (int pos = 0; pos < vk.nangl; pos++) {
      angl_i_atoms.putHost(vk.angl_i_atoms[pos] + synth_atom_base, synth_angl_offset + pos);
      angl_j_atoms.putHost(vk.angl_j_atoms[pos] + synth_atom_base, synth_angl_offset + pos);
      angl_k_atoms.putHost(vk.angl_k_atoms[pos] + synth_atom_base, synth_angl_offset + pos);
    }
    for (int pos = 0; pos < vk.ndihe; pos++) {
      dihe_i_atoms.putHost(vk.dihe_i_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
      dihe_j_atoms.putHost(vk.dihe_j_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
      dihe_k_atoms.putHost(vk.dihe_k_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
      dihe_l_atoms.putHost(vk.dihe_l_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
    }
  }

  // Fill in the restraint term indexing arrays
  nmr_int_data.resize((5 * rposn_offset) + (4 * rbond_offset) + (5 * rangl_offset) +
                      (6 * rdihe_offset));
  pivot = 0;
  rposn_atoms.setPointer(&nmr_int_data, pivot, rposn_offset);
  pivot += rposn_offset;
  rposn_kr_param_idx.setPointer(&nmr_int_data, pivot, rposn_offset);
  pivot += rposn_offset;
  rposn_xyz_param_idx.setPointer(&nmr_int_data, pivot, rposn_offset);
  pivot += rposn_offset;
  rbond_i_atoms.setPointer(&nmr_int_data, pivot, rbond_offset);
  pivot += rbond_offset;
  rbond_j_atoms.setPointer(&nmr_int_data, pivot, rbond_offset);
  pivot += rbond_offset;
  rbond_param_idx.setPointer(&nmr_int_data, pivot, rbond_offset);
  pivot += rbond_offset;
  rangl_i_atoms.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rangl_j_atoms.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rangl_k_atoms.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rangl_param_idx.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rdihe_i_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_j_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_k_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_l_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_param_idx.setPointer(&nmr_int_data, pivot, rdihe_offset);
  for (int sysid = 0; sysid < system_count; sysid++) {
    const RestraintApparatus* ra_ptr = restraint_networks[restraint_indices.readHost(sysid)];
    if (ra_ptr == nullptr) {
      continue;
    }
    const RestraintKit<double, double2, double4> rar = ra_ptr->dpData();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_rposn_offset = posn_restraint_offsets.readHost(sysid);
    const int synth_rbond_offset = bond_restraint_offsets.readHost(sysid);
    const int synth_rangl_offset = angl_restraint_offsets.readHost(sysid);
    const int synth_rdihe_offset = dihe_restraint_offsets.readHost(sysid);
    for (int pos = 0; pos < rar.nposn; pos++) {
      rposn_atoms.putHost(rar.rposn_atoms[pos] + synth_atom_base, synth_rposn_offset + pos);
    }
    for (int pos = 0; pos < rar.nbond; pos++) {
      rbond_i_atoms.putHost(rar.rbond_i_atoms[pos] + synth_atom_base, synth_rbond_offset + pos);
      rbond_j_atoms.putHost(rar.rbond_j_atoms[pos] + synth_atom_base, synth_rbond_offset + pos);
    }
    for (int pos = 0; pos < rar.nangl; pos++) {
      rangl_i_atoms.putHost(rar.rangl_i_atoms[pos] + synth_atom_base, synth_rangl_offset + pos);
      rangl_j_atoms.putHost(rar.rangl_j_atoms[pos] + synth_atom_base, synth_rangl_offset + pos);
      rangl_k_atoms.putHost(rar.rangl_k_atoms[pos] + synth_atom_base, synth_rangl_offset + pos);
    }
    for (int pos = 0; pos < rar.ndihe; pos++) {
      rdihe_i_atoms.putHost(rar.rdihe_i_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
      rdihe_j_atoms.putHost(rar.rdihe_j_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
      rdihe_k_atoms.putHost(rar.rdihe_k_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
      rdihe_l_atoms.putHost(rar.rdihe_l_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
    }
  }

  // Fill in the virtual site indexing arrays
  int vste_map_size = 0;
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph *ag_ptr = topologies[i];
    const int ni_vs = ag_ptr->getVirtualSiteParameterSetCount();
    vste_map_size += roundUp(ni_vs, warp_size_int);
  }
  vsite_int_data.resize(6 * vste_offset);
  pivot = 0;
  virtual_site_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame1_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame2_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame3_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame4_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_parameter_indices.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  for (int sysid = 0; sysid < system_count; sysid++) {
    const AtomGraph *ag_ptr = topologies[topology_indices.readHost(sysid)];
    const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_vste_offset = virtual_site_offsets.readHost(sysid);
    for (int pos = 0; pos < vsk.nsite; pos++) {
      const int svpos = synth_vste_offset + pos;
      virtual_site_atoms.putHost(vsk.vs_atoms[pos] + synth_atom_base, svpos);
      virtual_site_frame1_atoms.putHost(vsk.frame1_idx[pos] + synth_atom_base, svpos);
      virtual_site_frame2_atoms.putHost(vsk.frame2_idx[pos] + synth_atom_base, svpos);
      virtual_site_frame3_atoms.putHost(vsk.frame3_idx[pos] + synth_atom_base, svpos);
      virtual_site_frame4_atoms.putHost(vsk.frame4_idx[pos] + synth_atom_base, svpos);
    }
  }

  // Fill in the constraint group indexing arrays
  settle_group_indexing.resize(sett_offset);
  constraint_group_bounds.resize(cnst_bounds_offset);
  constraint_group_indices.resize(cnst_atom_offset);
  constraint_group_param_idx.resize(cnst_offset);
  int cnst_atom_count = 0;
  for (int sysid = 0; sysid < system_count; sysid++) {
    const AtomGraph *ag_ptr = topologies[topology_indices.readHost(sysid)];
    const ConstraintKit<double> cnk = ag_ptr->getDoublePrecisionConstraintKit();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_sett_offset = sett_group_offsets.readHost(sysid);
    const int synth_cnst_offset = cnst_group_offsets.readHost(sysid);
    for (int pos = 0; pos < cnk.nsettle; pos++) {

      // The SETTLE parameter group type is not yet known--but it is best to package four integers
      // for the oxygen, two hydrogens, and the parameter set into one tuple, so use a placeholder
      // for the type index in the w member.
      const int4 tmp_i4 = { synth_atom_base + cnk.settle_ox_atoms[pos],
                            synth_atom_base + cnk.settle_h1_atoms[pos],
                            synth_atom_base + cnk.settle_h2_atoms[pos], -1 };
      settle_group_indexing.putHost(tmp_i4, synth_sett_offset + pos);
    }
    for (int pos = 0; pos < cnk.ngroup; pos++) {
      for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
        constraint_group_indices.putHost(synth_atom_base + cnk.group_list[j], cnst_atom_count + j);
      }
      const int grp_len = cnk.group_bounds[pos + 1] - cnk.group_bounds[pos];
      const int2 tmp_i2 = { cnst_atom_count, cnst_atom_count + grp_len };
      constraint_group_bounds.putHost(tmp_i2, synth_cnst_offset + pos);
    }
    cnst_atom_count += roundUp(cnk.group_bounds[cnk.ngroup], warp_size_int);
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::condenseParameterTables() {

  // Compute the numbers of unique parameters.  Take the opportunity to compute offsets (starting
  // bounds) for various sets of terms.  While similar to arrays like bond_term_offsets already
  // stored in the AtomGraphSynthesis, these arrays pertain to term offsets for a list of unique
  // topologies, not a list of all systems (which could include multiple coordinate sets
  // referencing a single topology).
  int bond_offset = 0;
  int angl_offset = 0;
  int dihe_offset = 0;
  int ubrd_offset = 0;
  int cimp_offset = 0;
  int cmap_offset = 0;
  int attn_offset = 0;
  int chrg_offset = 0;
  int vste_offset = 0;
  int sett_offset = 0;
  int cnst_offset = 0;
  std::vector<int2> tmp_ubrd_param_map_bounds(topology_count);
  std::vector<int2> tmp_cimp_param_map_bounds(topology_count);
  std::vector<int2> tmp_cmap_param_map_bounds(topology_count);
  std::vector<int2> tmp_bond_param_map_bounds(topology_count);
  std::vector<int2> tmp_angl_param_map_bounds(topology_count);
  std::vector<int2> tmp_dihe_param_map_bounds(topology_count);
  std::vector<int2> tmp_attn_param_map_bounds(topology_count);
  std::vector<int2> tmp_vste_param_map_bounds(topology_count);
  std::vector<int2> tmp_sett_param_map_bounds(topology_count);
  std::vector<int2> tmp_cnst_param_map_bounds(topology_count);
  std::vector<int> topology_bond_table_offsets(topology_count);
  std::vector<int> topology_angl_table_offsets(topology_count);
  std::vector<int> topology_dihe_table_offsets(topology_count);
  std::vector<int> topology_ubrd_table_offsets(topology_count);
  std::vector<int> topology_cimp_table_offsets(topology_count);
  std::vector<int> topology_cmap_table_offsets(topology_count);
  std::vector<int> topology_attn_table_offsets(topology_count);
  std::vector<int> topology_chrg_table_offsets(topology_count);
  std::vector<int> topology_vste_table_offsets(topology_count);
  std::vector<int> topology_sett_table_offsets(topology_count);
  std::vector<int> topology_cnst_table_offsets(topology_count);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph *ag_ptr = topologies[i];
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
    const ConstraintKit<double> cnk  = ag_ptr->getDoublePrecisionConstraintKit();
    topology_bond_table_offsets[i] = bond_offset;
    topology_angl_table_offsets[i] = angl_offset;
    topology_dihe_table_offsets[i] = dihe_offset;
    topology_ubrd_table_offsets[i] = ubrd_offset;
    topology_cimp_table_offsets[i] = cimp_offset;
    topology_cmap_table_offsets[i] = cmap_offset;
    topology_attn_table_offsets[i] = attn_offset;
    topology_chrg_table_offsets[i] = chrg_offset;
    topology_vste_table_offsets[i] = vste_offset;
    topology_sett_table_offsets[i] = sett_offset;
    topology_cnst_table_offsets[i] = cnst_offset;
    tmp_bond_param_map_bounds[i].x = bond_offset;
    tmp_angl_param_map_bounds[i].x = angl_offset;
    tmp_dihe_param_map_bounds[i].x = dihe_offset;
    tmp_ubrd_param_map_bounds[i].x = ubrd_offset;
    tmp_cimp_param_map_bounds[i].x = cimp_offset;
    tmp_cmap_param_map_bounds[i].x = cmap_offset;
    tmp_attn_param_map_bounds[i].x = attn_offset;
    tmp_vste_param_map_bounds[i].x = vste_offset;
    tmp_sett_param_map_bounds[i].x = sett_offset;
    tmp_cnst_param_map_bounds[i].x = cnst_offset;
    tmp_bond_param_map_bounds[i].y = bond_offset + vk.nbond_param;
    tmp_angl_param_map_bounds[i].y = angl_offset + vk.nangl_param;
    tmp_dihe_param_map_bounds[i].y = dihe_offset + vk.ndihe_param;
    tmp_ubrd_param_map_bounds[i].y = ubrd_offset + vk.nubrd_param;
    tmp_cimp_param_map_bounds[i].y = cimp_offset + vk.ncimp_param;
    tmp_cmap_param_map_bounds[i].y = cmap_offset + vk.ncmap_surf;
    tmp_attn_param_map_bounds[i].y = attn_offset + vk.nattn14_param;
    tmp_vste_param_map_bounds[i].y = vste_offset + vsk.nframe_set;
    tmp_sett_param_map_bounds[i].y = sett_offset + cnk.nsett_param;
    tmp_cnst_param_map_bounds[i].y = cnst_offset + cnk.ncnst_param;
    bond_offset += roundUp(vk.nbond_param, warp_size_int);
    angl_offset += roundUp(vk.nangl_param, warp_size_int);
    dihe_offset += roundUp(vk.ndihe_param, warp_size_int);
    ubrd_offset += roundUp(vk.nubrd_param, warp_size_int);
    cimp_offset += roundUp(vk.ncimp_param, warp_size_int);
    cmap_offset += roundUp(vk.ncmap_surf, warp_size_int);
    attn_offset += roundUp(vk.nattn14_param, warp_size_int);
    chrg_offset += roundUp(nbk.n_q_types, warp_size_int);
    vste_offset += roundUp(vsk.nframe_set, warp_size_int);
    sett_offset += roundUp(cnk.nsett_param, warp_size_int);
    cnst_offset += roundUp(cnk.ncnst_param, warp_size_int);
  }
  
  // Pre-compute some quantities relating to CMAPs that will help distinguish these gargantuan
  // "parameters" in the inner loops that follow.
  std::vector<double> cmap_parameter_sums(cmap_offset);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph *iag_ptr = topologies[i];
    const ValenceKit<double> i_vk = iag_ptr->getDoublePrecisionValenceKit();
    for (int j = 0; j < i_vk.ncmap_surf; j++) {
      cmap_parameter_sums[topology_cmap_table_offsets[i] + j] =
        sum<double>(&i_vk.cmap_surf[i_vk.cmap_surf_bounds[j]],
                    i_vk.cmap_dim[j] * i_vk.cmap_dim[j]);
    }
  }
  
  // Create lists of unique parameters for the valence and non-bonded calculations.
  std::vector<int> bond_synthesis_index(bond_offset, -1);
  std::vector<int> angl_synthesis_index(angl_offset, -1);
  std::vector<int> dihe_synthesis_index(dihe_offset, -1);
  std::vector<int> ubrd_synthesis_index(ubrd_offset, -1);
  std::vector<int> cimp_synthesis_index(cimp_offset, -1);
  std::vector<int> cmap_synthesis_index(cmap_offset, -1);
  std::vector<int> attn_synthesis_index(attn_offset, -1);
  std::vector<int> chrg_synthesis_index(chrg_offset, -1);
  std::vector<int> vste_synthesis_index(vste_offset, -1);
  std::vector<int> sett_synthesis_index(sett_offset, -1);
  std::vector<int> cnst_synthesis_index(cnst_offset, -1);
  std::vector<double> filtered_bond_keq, filtered_bond_leq;
  std::vector<float> sp_filtered_bond_keq, sp_filtered_bond_leq;
  std::vector<double> filtered_angl_keq, filtered_angl_theta;
  std::vector<float> sp_filtered_angl_keq, sp_filtered_angl_theta;
  std::vector<double> filtered_dihe_amp, filtered_dihe_freq, filtered_dihe_phi;
  std::vector<float> sp_filtered_dihe_amp, sp_filtered_dihe_freq, sp_filtered_dihe_phi;
  std::vector<double> filtered_ubrd_keq,filtered_ubrd_leq;
  std::vector<float> sp_filtered_ubrd_keq, sp_filtered_ubrd_leq;
  std::vector<double> filtered_cimp_keq, filtered_cimp_phi;
  std::vector<float> sp_filtered_cimp_keq, sp_filtered_cimp_phi;
  std::vector<int> filtered_cmap_dim;
  std::vector<int> filtered_cmap_surf_bounds(1, 0);
  std::vector<double> filtered_cmap_surf;
  std::vector<float> sp_filtered_cmap_surf;
  std::vector<double> filtered_attn14_elec, filtered_attn14_vdw;
  std::vector<float> sp_filtered_attn14_elec, sp_filtered_attn14_vdw;
  std::vector<double> filtered_chrg;
  std::vector<float> sp_filtered_chrg;
  std::vector<double4> filtered_vste_params, filtered_sett_geom;
  std::vector<float4> sp_filtered_vste_params, sp_filtered_sett_geom;
  std::vector<double2> filtered_sett_mass, filtered_cnst_group_params;
  std::vector<float2> sp_filtered_sett_mass, sp_filtered_cnst_group_params;
  std::vector<int> filtered_cnst_group_bounds(1, 0);
  int n_unique_bond = 0;
  int n_unique_angl = 0;
  int n_unique_dihe = 0;
  int n_unique_ubrd = 0;
  int n_unique_cimp = 0;
  int n_unique_cmap = 0;
  int n_unique_attn = 0;
  int n_unique_chrg = 0;
  int n_unique_vste = 0;
  int n_unique_sett = 0;
  int n_unique_cnst = 0;
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph *iag_ptr = topologies[i];
    const ChemicalDetailsKit i_cdk       = iag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> i_nbk     = iag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> i_vk        = iag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> i_vsk   = iag_ptr->getDoublePrecisionVirtualSiteKit();
    const ConstraintKit<double> i_cnk    = iag_ptr->getDoublePrecisionConstraintKit();
    const NonbondedKit<float> i_nbk_sp   = iag_ptr->getSinglePrecisionNonbondedKit();
    const ValenceKit<float> i_vk_sp      = iag_ptr->getSinglePrecisionValenceKit();
    const VirtualSiteKit<float> i_vsk_sp = iag_ptr->getSinglePrecisionVirtualSiteKit();
    const ConstraintKit<float> i_cnk_sp  = iag_ptr->getSinglePrecisionConstraintKit();

    // Seek out unique bond parameters
    for (int j = 0; j < i_vk.nbond_param; j++) {
      if (bond_synthesis_index[topology_bond_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_bond_keq(i_vk.bond_keq[j], constants::verytiny);
      const Approx ij_bond_leq(i_vk.bond_leq[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nbond_param; m++) {
          if (bond_synthesis_index[topology_bond_table_offsets[k] + m] < 0 &&
              ij_bond_keq.test(k_vk.bond_keq[m]) && ij_bond_leq.test(k_vk.bond_leq[m])) {
            bond_synthesis_index[topology_bond_table_offsets[k] + m] = n_unique_bond;
          }
        }
      }

      // Catalog this unique bond and increment the counter
      filtered_bond_keq.push_back(i_vk.bond_keq[j]);
      filtered_bond_leq.push_back(i_vk.bond_leq[j]);
      sp_filtered_bond_keq.push_back(i_vk_sp.bond_keq[j]);
      sp_filtered_bond_leq.push_back(i_vk_sp.bond_leq[j]);
      n_unique_bond++;
    }

    // Seek out unique angle parameters
    for (int j = 0; j < i_vk.nangl_param; j++) {
      if (angl_synthesis_index[topology_angl_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_angl_keq(i_vk.angl_keq[j], constants::verytiny);
      const Approx ij_angl_theta(i_vk.angl_theta[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nangl_param; m++) {
          if (angl_synthesis_index[topology_angl_table_offsets[k] + m] < 0 &&
              ij_angl_keq.test(k_vk.angl_keq[m]) && ij_angl_theta.test(k_vk.angl_theta[m])) {
            angl_synthesis_index[topology_angl_table_offsets[k] + m] = n_unique_angl;
          }
        }
      }

      // Catalog this unique angle and increment the counter
      filtered_angl_keq.push_back(i_vk.angl_keq[j]);
      filtered_angl_theta.push_back(i_vk.angl_theta[j]);
      sp_filtered_angl_keq.push_back(i_vk_sp.angl_keq[j]);
      sp_filtered_angl_theta.push_back(i_vk_sp.angl_theta[j]);
      n_unique_angl++;
    }

    // Seek out unique dihedral parameters
    for (int j = 0; j < i_vk.ndihe_param; j++) {
      if (dihe_synthesis_index[topology_dihe_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_dihe_amp(i_vk.dihe_amp[j], constants::verytiny);
      const Approx ij_dihe_freq(i_vk.dihe_freq[j], constants::verytiny);
      const Approx ij_dihe_phi(i_vk.dihe_phi[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.ndihe_param; m++) {
          if (dihe_synthesis_index[topology_dihe_table_offsets[k] + m] < 0 &&
              ij_dihe_amp.test(k_vk.dihe_amp[m]) && ij_dihe_freq.test(k_vk.dihe_freq[m]) &&
              ij_dihe_phi.test(k_vk.dihe_phi[m])) {
            dihe_synthesis_index[topology_dihe_table_offsets[k] + m] = n_unique_dihe;
          }
        }
      }

      // Catalog this unique dihedral and increment the counter
      filtered_dihe_amp.push_back(i_vk.dihe_amp[j]);
      filtered_dihe_freq.push_back(i_vk.dihe_freq[j]);
      filtered_dihe_phi.push_back(i_vk.dihe_phi[j]);
      sp_filtered_dihe_amp.push_back(i_vk_sp.dihe_amp[j]);
      sp_filtered_dihe_freq.push_back(i_vk_sp.dihe_freq[j]);
      sp_filtered_dihe_phi.push_back(i_vk_sp.dihe_phi[j]);
      n_unique_dihe++;
    }
    
    // Seek out unique Urey-Bradley parameters
    for (int j = 0; j < i_vk.nubrd_param; j++) {
      if (ubrd_synthesis_index[topology_ubrd_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_ubrd_keq(i_vk.ubrd_keq[j], constants::verytiny);
      const Approx ij_ubrd_leq(i_vk.ubrd_leq[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nubrd_param; m++) {
          if (ubrd_synthesis_index[topology_ubrd_table_offsets[k] + m] < 0 &&
              ij_ubrd_keq.test(k_vk.ubrd_keq[m]) && ij_ubrd_leq.test(k_vk.ubrd_leq[m])) {
            ubrd_synthesis_index[topology_ubrd_table_offsets[k] + m] = n_unique_ubrd;
          }
        }
      }

      // Catalog this unique Urey-Bradley parameter set and increment the counter
      filtered_ubrd_keq.push_back(i_vk.ubrd_keq[j]);
      filtered_ubrd_leq.push_back(i_vk.ubrd_leq[j]);
      sp_filtered_ubrd_keq.push_back(i_vk_sp.ubrd_keq[j]);
      sp_filtered_ubrd_leq.push_back(i_vk_sp.ubrd_leq[j]);
      n_unique_ubrd++;
    }

    // Seek out unique CHARMM improper dihedral parameters
    for (int j = 0; j < i_vk.ncimp_param; j++) {
      if (cimp_synthesis_index[topology_cimp_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_cimp_keq(i_vk.cimp_keq[j], constants::verytiny);
      const Approx ij_cimp_phi(i_vk.cimp_phi[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.ncimp_param; m++) {
          if (cimp_synthesis_index[topology_cimp_table_offsets[k] + m] < 0 &&
              ij_cimp_keq.test(k_vk.cimp_keq[m]) && ij_cimp_phi.test(k_vk.cimp_phi[m])) {
            cimp_synthesis_index[topology_cimp_table_offsets[k] + m] = n_unique_cimp;
          }
        }
      }

      // Catalog this unique CHARMM improper dihedral and increment the counter
      filtered_cimp_keq.push_back(i_vk.cimp_keq[j]);
      filtered_cimp_phi.push_back(i_vk.cimp_phi[j]);
      sp_filtered_cimp_keq.push_back(i_vk_sp.cimp_keq[j]);
      sp_filtered_cimp_phi.push_back(i_vk_sp.cimp_phi[j]);
      n_unique_cimp++;
    }

    // Seek out unique CMAP surfaces
    for (int j = 0; j < i_vk.ncmap_surf; j++) {
      if (cmap_synthesis_index[topology_cmap_table_offsets[i] + j] >= 0) {
        continue;
      }

      // Comparing CMAP surfaces is a much more involved process than comparing other parameters.
      // Use a pre-arranged table of surface sums to expedite the process.
      const Approx ij_cmap_sum(cmap_parameter_sums[topology_cmap_table_offsets[i] + j],
                               constants::verytiny);
      const int ij_cmap_dim = i_vk.cmap_dim[j];
      const double* ij_surf_ptr = &i_vk.cmap_surf[i_vk.cmap_surf_bounds[j]];
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk = kag_ptr->getDoublePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.ncmap_surf; m++) {
          if (cmap_synthesis_index[topology_cmap_table_offsets[k] + m] < 0 &&
              ij_cmap_dim == k_vk.cmap_dim[m] &&
              ij_cmap_sum.test(cmap_parameter_sums[topology_cmap_table_offsets[k] + m]) &&
              maxAbsoluteDifference(ij_surf_ptr, &k_vk.cmap_surf[k_vk.cmap_surf_bounds[m]],
                                    ij_cmap_dim * ij_cmap_dim) < constants::verytiny) {
            cmap_synthesis_index[topology_cmap_table_offsets[k] + m] = n_unique_cmap;
          }
        }
      }

      // Catalog this unique CMAP and increment the counter
      for (int k = 0; k < ij_cmap_dim * ij_cmap_dim; k++) {
        filtered_cmap_surf.push_back(i_vk.cmap_surf[i_vk.cmap_surf_bounds[j] + k]);
        sp_filtered_cmap_surf.push_back(i_vk_sp.cmap_surf[i_vk_sp.cmap_surf_bounds[j] + k]);
      }
      filtered_cmap_dim.push_back(ij_cmap_dim);
      filtered_cmap_surf_bounds.push_back(filtered_cmap_surf.size());
      n_unique_cmap++;
    }

    // Seek out unique attenuated 1:4 scaling factor pairs
    for (int j = 0; j < i_vk.nattn14_param; j++) {
      if (attn_synthesis_index[topology_attn_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_attn_qq(i_vk.attn14_elec[j], constants::verytiny);
      const Approx ij_attn_lj(i_vk.attn14_vdw[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nattn14_param; m++) {
          if (attn_synthesis_index[topology_attn_table_offsets[k] + m] < 0 &&
              ij_attn_qq.test(k_vk.attn14_elec[m]) && ij_attn_lj.test(k_vk.attn14_vdw[m])) {
            attn_synthesis_index[topology_attn_table_offsets[k] + m] = n_unique_attn;            
          }
        }
      }

      // Catalog this unique attenuation scaling factor pair
      filtered_attn14_elec.push_back(i_vk.attn14_elec[j]);
      filtered_attn14_vdw.push_back(i_vk.attn14_vdw[j]);
      sp_filtered_attn14_elec.push_back(i_vk_sp.attn14_elec[j]);
      sp_filtered_attn14_vdw.push_back(i_vk_sp.attn14_vdw[j]);
      n_unique_attn++;
    }

    // Seek out unique charges
    for (int j = 0; j < i_nbk.n_q_types; j++) {
      if (chrg_synthesis_index[topology_chrg_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_chrg(i_nbk.q_parameter[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const NonbondedKit<double> k_nbk   = kag_ptr->getDoublePrecisionNonbondedKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_nbk.n_q_types; m++) {
          if (chrg_synthesis_index[topology_chrg_table_offsets[k] + m] < 0 &&
              ij_chrg.test(k_nbk.q_parameter[m])) {
            chrg_synthesis_index[topology_chrg_table_offsets[k] + m] = n_unique_chrg;
          }
        }
      }

      // Catalog this unique charge and increment the counter
      filtered_chrg.push_back(i_nbk.q_parameter[j]);
      sp_filtered_chrg.push_back(i_nbk_sp.q_parameter[j]);
      n_unique_chrg++;
    }

    // Seek out unique virtual site frames
    for (int j = 0; j < i_vsk.nframe_set; j++) {
      if (vste_synthesis_index[topology_vste_table_offsets[i] + j] >= 0) {
        continue;
      }
      const int ij_frame_type = i_vsk.vs_types[j];
      const Approx ij_dim1(i_vsk.dim1[j], constants::verytiny);
      const Approx ij_dim2(i_vsk.dim2[j], constants::verytiny);
      const Approx ij_dim3(i_vsk.dim3[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const VirtualSiteKit<double> k_vsk = kag_ptr->getDoublePrecisionVirtualSiteKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vsk.nframe_set; m++) {
          if (vste_synthesis_index[topology_vste_table_offsets[k] + m] < 0 &&
              ij_frame_type == k_vsk.vs_types[m] && ij_dim1.test(k_vsk.dim1[m]) &&
              ij_dim2.test(k_vsk.dim2[m]) && ij_dim3.test(k_vsk.dim3[m])) {
            vste_synthesis_index[topology_vste_table_offsets[k] + m] = n_unique_vste;
          }
        }
      }

      // Catalog this unique virtual site frame and increment the counter
      filtered_vste_params.push_back({ i_vsk.dim1[j], i_vsk.dim2[j],
                                       i_vsk.dim3[j], static_cast<double>(ij_frame_type) });
      sp_filtered_vste_params.push_back({ i_vsk_sp.dim1[j], i_vsk_sp.dim2[j],
                                          i_vsk_sp.dim3[j], static_cast<float>(ij_frame_type) });
      n_unique_vste++;
    }
    
    // Seek out unique SETTLE parameter sets
    for (int j = 0; j < i_cnk.nsett_param; j++) {
      if (sett_synthesis_index[topology_sett_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_mormt(i_cnk.settle_mormt[j], constants::verytiny);
      const Approx ij_mhrmt(i_cnk.settle_mhrmt[j], constants::verytiny);
      const Approx ij_ra(i_cnk.settle_ra[j], constants::verytiny);
      const Approx ij_rb(i_cnk.settle_rb[j], constants::verytiny);
      const Approx ij_rc(i_cnk.settle_rc[j], constants::verytiny);
      const Approx ij_invra(i_cnk.settle_invra[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ConstraintKit<double> k_cnk = kag_ptr->getDoublePrecisionConstraintKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_cnk.nsett_param; m++) {
          if (sett_synthesis_index[topology_sett_table_offsets[k] + m] < 0 &&
              ij_mormt.test(k_cnk.settle_mormt[m]) && ij_mhrmt.test(k_cnk.settle_mhrmt[m]) &&
              ij_ra.test(k_cnk.settle_ra[m]) && ij_rb.test(k_cnk.settle_rb[m]) &&
              ij_rc.test(k_cnk.settle_rc[m]) && ij_invra.test(k_cnk.settle_invra[m])) {
            sett_synthesis_index[topology_sett_table_offsets[k] + m] = n_unique_sett;
          }
        }
      }

      // Catalog this unique SETTLE group geometry
      filtered_sett_mass.push_back({ i_cnk.settle_mormt[j], i_cnk.settle_mhrmt[j] });
      filtered_sett_geom.push_back({ i_cnk.settle_ra[j], i_cnk.settle_rb[j], i_cnk.settle_rc[j],
                                     i_cnk.settle_invra[j] });
      sp_filtered_sett_mass.push_back({ i_cnk_sp.settle_mormt[j], i_cnk_sp.settle_mhrmt[j] });
      sp_filtered_sett_geom.push_back({ i_cnk_sp.settle_ra[j], i_cnk_sp.settle_rb[j],
                                        i_cnk_sp.settle_rc[j], i_cnk_sp.settle_invra[j] });
      n_unique_sett++;
    }

    // Seek out unique bond length constraint parameter sets.  Track the largest such constraint
    // group.
    for (int j = 0; j < i_cnk.ncnst_param; j++) {
      if (cnst_synthesis_index[topology_cnst_table_offsets[i] + j] >= 0) {
        continue;
      }
      const int nelem = i_cnk.group_param_bounds[j + 1] - i_cnk.group_param_bounds[j];
      largest_constraint_group = std::max(largest_constraint_group, nelem - 1);
      std::vector<Approx> lengths, inv_masses;
      for (int k = i_cnk.group_param_bounds[j]; k < i_cnk.group_param_bounds[j + 1]; k++) {
        lengths.emplace_back(i_cnk.group_sq_lengths[k], ComparisonType::ABSOLUTE,
                             constants::verytiny);
        inv_masses.emplace_back(i_cnk.group_inv_masses[k], ComparisonType::ABSOLUTE,
                                constants::verytiny);
      }
      for (int k = i; k < topology_count; k++) {
        const AtomGraph *kag_ptr = topologies[k];
        const ConstraintKit<double> k_cnk = kag_ptr->getDoublePrecisionConstraintKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_cnk.ncnst_param; m++) {
          if (k_cnk.group_param_bounds[m + 1] - k_cnk.group_param_bounds[m] != nelem) {
            continue;
          }
          bool same = true;
          for (int n = 0; n < nelem; n++) {
            const double km_length = k_cnk.group_sq_lengths[k_cnk.group_param_bounds[m] + n];
            const double km_inv_mass = k_cnk.group_inv_masses[k_cnk.group_param_bounds[m] + n];
            same = (same && lengths[n].test(km_length) && inv_masses[n].test(km_inv_mass));
          }
          if (same) {
            cnst_synthesis_index[topology_cnst_table_offsets[k] + m] = n_unique_cnst;
          }
        }
      }

      // Catalog this unique hub-and-spoke constraint group
      for (int k = 0; k < nelem; k++) {
        filtered_cnst_group_params.push_back({ lengths[k].getValue(), inv_masses[0].getValue() +
                                                                      inv_masses[k].getValue() });
        float2 tmp_f2;
        tmp_f2.x = lengths[k].getValue();
        tmp_f2.y = (inv_masses[0].getValue() + inv_masses[k].getValue());
        sp_filtered_cnst_group_params.push_back(tmp_f2);
      }
      filtered_cnst_group_bounds.push_back(static_cast<int>(filtered_cnst_group_params.size()));
      n_unique_cnst++;
    }
  }

  // Post-process the condensed array of CMAP surfaces into patches
  const CmapAccessories cmap_digest = computeCmapDerivatives(n_unique_cmap, filtered_cmap_dim,
                                                             filtered_cmap_surf_bounds,
                                                             filtered_cmap_surf);

  // The charge parameters get their own arrays due to the time they are finally condensed and the
  // uncertainty of their size.  Other non-bonded parameter sets have dedicated ARRAY-kind Hybrid
  // objects as well.  The indexing array (charge_indices) is a POINTER-kind Hybrid targeting
  // chem_int_data.
  charge_parameters.resize(filtered_chrg.size());
  charge_parameters.putHost(filtered_chrg);
  sp_charge_parameters.resize(sp_filtered_chrg.size());
  sp_charge_parameters.putHost(sp_filtered_chrg);
  
  // Record synthesis parameters found thus far.  Pad the relevant ARRAY-kind Hybrid objects in
  // case some developer tries to have all lanes of a warp access an element.
  total_bond_params   = n_unique_bond;
  total_angl_params   = n_unique_angl;
  total_dihe_params   = n_unique_dihe;
  total_ubrd_params   = n_unique_ubrd;
  total_cimp_params   = n_unique_cimp;
  total_cmap_surfaces = n_unique_cmap;
  total_attn14_params = n_unique_attn;
  total_charge_types  = n_unique_chrg;
  total_vste_params   = n_unique_vste;
  const int rn_space = (2 * roundUp(total_bond_params, warp_size_int)) +
                       (2 * roundUp(total_angl_params, warp_size_int)) +
                       (3 * roundUp(total_dihe_params, warp_size_int)) +
                       (2 * roundUp(total_ubrd_params, warp_size_int)) +
                       (2 * roundUp(total_cimp_params, warp_size_int)) +
                       (2 * roundUp(total_attn14_params, warp_size_int)) +
                       roundUp(static_cast<int>(filtered_cmap_surf.size()), warp_size_int) +
                       roundUp(static_cast<int>(cmap_digest.patch_matrix_form.size()),
                               warp_size_int);
  valparam_double_data.resize(rn_space);
  valparam_float_data.resize(rn_space);
  valparam_int_data.resize(roundUp(total_cmap_surfaces, warp_size_int) +
                           (2 * roundUp(total_cmap_surfaces + 1, warp_size_int)) + ubrd_offset +
                           cimp_offset + cmap_offset + bond_offset + angl_offset + dihe_offset +
                           attn_offset + vste_offset + sett_offset + cnst_offset);
  size_t ic = 0LLU;
  ic = bond_stiffnesses.putHost(&valparam_double_data, filtered_bond_keq, ic, warp_size_zu);
  ic = bond_equilibria.putHost(&valparam_double_data, filtered_bond_leq, ic, warp_size_zu);
  ic = angl_stiffnesses.putHost(&valparam_double_data, filtered_angl_keq, ic, warp_size_zu);
  ic = angl_equilibria.putHost(&valparam_double_data, filtered_angl_theta, ic, warp_size_zu);
  ic = dihe_amplitudes.putHost(&valparam_double_data, filtered_dihe_amp, ic, warp_size_zu);
  ic = dihe_periodicities.putHost(&valparam_double_data, filtered_dihe_freq, ic, warp_size_zu);
  ic = dihe_phase_angles.putHost(&valparam_double_data, filtered_dihe_phi, ic, warp_size_zu);
  ic = ubrd_stiffnesses.putHost(&valparam_double_data, filtered_ubrd_keq, ic, warp_size_zu);
  ic = ubrd_equilibria.putHost(&valparam_double_data, filtered_ubrd_leq, ic, warp_size_zu);
  ic = cimp_stiffnesses.putHost(&valparam_double_data, filtered_cimp_keq, ic, warp_size_zu);
  ic = cimp_phase_angles.putHost(&valparam_double_data, filtered_cimp_phi, ic, warp_size_zu);
  ic = attn14_elec_factors.putHost(&valparam_double_data, filtered_attn14_elec, ic, warp_size_zu);
  ic = attn14_vdw_factors.putHost(&valparam_double_data, filtered_attn14_vdw, ic, warp_size_zu);
  ic = cmap_surfaces.putHost(&valparam_double_data, filtered_cmap_surf, ic, warp_size_zu);
  ic = cmap_patches.putHost(&valparam_double_data, cmap_digest.patch_matrix_form, ic,
                            warp_size_zu);
  ic = 0LLU;
  ic = sp_bond_stiffnesses.putHost(&valparam_float_data, sp_filtered_bond_keq, ic, warp_size_zu);
  ic = sp_bond_equilibria.putHost(&valparam_float_data, sp_filtered_bond_leq, ic, warp_size_zu);
  ic = sp_angl_stiffnesses.putHost(&valparam_float_data, sp_filtered_angl_keq, ic, warp_size_zu);
  ic = sp_angl_equilibria.putHost(&valparam_float_data, sp_filtered_angl_theta, ic, warp_size_zu);
  ic = sp_dihe_amplitudes.putHost(&valparam_float_data, sp_filtered_dihe_amp, ic, warp_size_zu);
  ic = sp_dihe_periodicities.putHost(&valparam_float_data, sp_filtered_dihe_freq, ic,
                                     warp_size_zu);
  ic = sp_dihe_phase_angles.putHost(&valparam_float_data, sp_filtered_dihe_phi, ic, warp_size_zu);
  ic = sp_ubrd_stiffnesses.putHost(&valparam_float_data, sp_filtered_ubrd_keq, ic, warp_size_zu);
  ic = sp_ubrd_equilibria.putHost(&valparam_float_data, sp_filtered_ubrd_leq, ic, warp_size_zu);
  ic = sp_cimp_stiffnesses.putHost(&valparam_float_data, sp_filtered_cimp_keq, ic, warp_size_zu);
  ic = sp_cimp_phase_angles.putHost(&valparam_float_data, sp_filtered_cimp_phi, ic, warp_size_zu);
  ic = sp_attn14_elec_factors.putHost(&valparam_float_data, sp_filtered_attn14_elec, ic,
                                      warp_size_zu);
  ic = sp_attn14_vdw_factors.putHost(&valparam_float_data, sp_filtered_attn14_vdw, ic,
                                     warp_size_zu);
  ic = sp_cmap_surfaces.putHost(&valparam_float_data, sp_filtered_cmap_surf, ic, warp_size_zu);
  ic = sp_cmap_patches.putHost(&valparam_float_data,
                               std::vector<float>(cmap_digest.patch_matrix_form.begin(),
                                                  cmap_digest.patch_matrix_form.end()), ic,
                               warp_size_zu);
  ic = 0LLU;
  ic = cmap_surface_dimensions.putHost(&valparam_int_data, filtered_cmap_dim, ic, warp_size_zu);
  ic = cmap_surface_bounds.putHost(&valparam_int_data, filtered_cmap_surf_bounds, ic,
                                   warp_size_zu);
  ic = cmap_patch_bounds.putHost(&valparam_int_data, cmap_digest.patch_matrix_bounds, ic,
                                 warp_size_zu);
  ic = ubrd_param_map.putHost(&valparam_int_data, ubrd_synthesis_index, ic, warp_size_zu);
  ic = cimp_param_map.putHost(&valparam_int_data, cimp_synthesis_index, ic, warp_size_zu);
  ic = cmap_param_map.putHost(&valparam_int_data, cmap_synthesis_index, ic, warp_size_zu);
  ic = bond_param_map.putHost(&valparam_int_data, bond_synthesis_index, ic, warp_size_zu);
  ic = angl_param_map.putHost(&valparam_int_data, angl_synthesis_index, ic, warp_size_zu);
  ic = dihe_param_map.putHost(&valparam_int_data, dihe_synthesis_index, ic, warp_size_zu);
  ic = attn14_param_map.putHost(&valparam_int_data, attn_synthesis_index, ic, warp_size_zu);
  ic = vste_param_map.putHost(&valparam_int_data, vste_synthesis_index, ic, warp_size_zu);
  ic = sett_param_map.putHost(&valparam_int_data, sett_synthesis_index, ic, warp_size_zu);
  ic = cnst_param_map.putHost(&valparam_int_data, cnst_synthesis_index, ic, warp_size_zu);
  const int tc_offset = roundUp(topology_count, warp_size_int);
  valparam_int2_data.resize(10 * tc_offset);
  ic = 0LLU;
  ic = ubrd_param_map_bounds.putHost(&valparam_int2_data, tmp_ubrd_param_map_bounds, ic,
                                     warp_size_zu);
  ic = cimp_param_map_bounds.putHost(&valparam_int2_data, tmp_cimp_param_map_bounds, ic,
                                     warp_size_zu);
  ic = cmap_param_map_bounds.putHost(&valparam_int2_data, tmp_cmap_param_map_bounds, ic,
                                     warp_size_zu);
  ic = bond_param_map_bounds.putHost(&valparam_int2_data, tmp_bond_param_map_bounds, ic,
                                     warp_size_zu);
  ic = angl_param_map_bounds.putHost(&valparam_int2_data, tmp_angl_param_map_bounds, ic,
                                     warp_size_zu);
  ic = dihe_param_map_bounds.putHost(&valparam_int2_data, tmp_dihe_param_map_bounds, ic,
                                     warp_size_zu);
  ic = attn14_param_map_bounds.putHost(&valparam_int2_data, tmp_attn_param_map_bounds, ic,
                                       warp_size_zu);
  ic = vste_param_map_bounds.putHost(&valparam_int2_data, tmp_vste_param_map_bounds, ic,
                                     warp_size_zu);
  ic = sett_param_map_bounds.putHost(&valparam_int2_data, tmp_sett_param_map_bounds, ic,
                                     warp_size_zu);
  ic = cnst_param_map_bounds.putHost(&valparam_int2_data, tmp_cnst_param_map_bounds, ic,
                                     warp_size_zu);

  // Some parameter arrays have unusual types and are held in ARRAY-kind Hybrid objects
  virtual_site_parameters.resize(filtered_vste_params.size());
  virtual_site_parameters.putHost(filtered_vste_params);
  sp_virtual_site_parameters.resize(sp_filtered_vste_params.size());
  sp_virtual_site_parameters.putHost(sp_filtered_vste_params);
  settle_group_masses.resize(filtered_sett_mass.size());
  settle_group_masses.putHost(filtered_sett_mass);
  settle_group_geometry.resize(filtered_sett_geom.size());
  settle_group_geometry.putHost(filtered_sett_geom);
  sp_settle_group_masses.resize(sp_filtered_sett_mass.size());
  sp_settle_group_masses.putHost(sp_filtered_sett_mass);
  sp_settle_group_geometry.resize(sp_filtered_sett_geom.size());
  sp_settle_group_geometry.putHost(sp_filtered_sett_geom);
  constraint_param_bounds.resize(filtered_cnst_group_bounds.size());
  constraint_param_bounds.putHost(filtered_cnst_group_bounds);
  constraint_group_params.resize(filtered_cnst_group_params.size());
  constraint_group_params.putHost(filtered_cnst_group_params);
  sp_constraint_group_params.resize(sp_filtered_cnst_group_params.size());
  sp_constraint_group_params.putHost(sp_filtered_cnst_group_params);
  
  // Loop back over all systems and copy the known mapping of individual topologies to the
  // synthesis as a whole.  Fill out the parameter maps.
  for (int i = 0; i < system_count; i++) {
    const int tp_index = topology_indices.readHost(i);
    const AtomGraph *iag_ptr = topologies[tp_index];
    const NonbondedKit<double> i_nbk   = iag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> i_vk      = iag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> i_vsk = iag_ptr->getDoublePrecisionVirtualSiteKit();
    const int ag_bond_table_offset = topology_bond_table_offsets[tp_index];
    const int ag_angl_table_offset = topology_angl_table_offsets[tp_index];
    const int ag_dihe_table_offset = topology_dihe_table_offsets[tp_index];
    const int ag_ubrd_table_offset = topology_ubrd_table_offsets[tp_index];
    const int ag_cimp_table_offset = topology_cimp_table_offsets[tp_index];
    const int ag_cmap_table_offset = topology_cmap_table_offsets[tp_index];
    const int ag_vste_table_offset = topology_vste_table_offsets[tp_index];
    const int ag_chrg_table_offset = topology_chrg_table_offsets[tp_index];
    const int synth_bond_term_offset = bond_term_offsets.readHost(i);
    const int synth_angl_term_offset = angl_term_offsets.readHost(i);
    const int synth_dihe_term_offset = dihe_term_offsets.readHost(i);
    const int synth_ubrd_term_offset = ubrd_term_offsets.readHost(i);
    const int synth_cimp_term_offset = cimp_term_offsets.readHost(i);
    const int synth_cmap_term_offset = cmap_term_offsets.readHost(i);
    const int synth_vste_term_offset = virtual_site_offsets.readHost(i);
    const int synth_atom_offset = atom_offsets.readHost(i);
    for (int j = 0; j < i_vk.nbond; j++) {
      bond_param_idx.putHost(bond_synthesis_index[ag_bond_table_offset + i_vk.bond_param_idx[j]],
                             synth_bond_term_offset + j);
    }
    for (int j = 0; j < i_vk.nangl; j++) {
      angl_param_idx.putHost(angl_synthesis_index[ag_angl_table_offset + i_vk.angl_param_idx[j]],
                             synth_angl_term_offset + j);
    }
    for (int j = 0; j < i_vk.ndihe; j++) {
      dihe_param_idx.putHost(dihe_synthesis_index[ag_dihe_table_offset + i_vk.dihe_param_idx[j]],
                             synth_dihe_term_offset + j);
    }
    for (int j = 0; j < i_vk.nubrd; j++) {
      ubrd_param_idx.putHost(ubrd_synthesis_index[ag_ubrd_table_offset + i_vk.ubrd_param_idx[j]],
                             synth_ubrd_term_offset + j);
    }
    for (int j = 0; j < i_vk.ncimp; j++) {
      cimp_param_idx.putHost(cimp_synthesis_index[ag_cimp_table_offset + i_vk.cimp_param_idx[j]],
                             synth_cimp_term_offset + j);
    }
    for (int j = 0; j < i_vk.ncmap; j++) {
      cmap_param_idx.putHost(cmap_synthesis_index[ag_cmap_table_offset + i_vk.cmap_surf_idx[j]],
                             synth_cmap_term_offset + j);
    }
    for (int j = 0; j < i_nbk.natom; j++) {
      charge_indices.putHost(chrg_synthesis_index[ag_chrg_table_offset + i_nbk.q_idx[j]],
                             synth_atom_offset + j);
    }
    for (int j = 0; j < i_vsk.nsite; j++) {
      virtual_site_parameter_indices.putHost(vste_synthesis_index[ag_vste_table_offset +
                                                                  i_vsk.vs_param_idx[j]],
                                             synth_vste_term_offset + j);
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::extendLJMatrices() {

  // Identify systems with identical Lennard-Jones tables and have their offsets point to the
  // same table in the synthesis.  This may conserve some L1 cache if everything else is done
  // right.  Track the necessary array size to hold all tables.
  int table_acc = 0;
  std::vector<NonbondedKit<double>> nbkvec;
  nbkvec.reserve(topology_count);
  std::vector<int> table_idx(topology_count, -1);
  for (int i = 0; i < topology_count; i++) {
    nbkvec.push_back(topologies[i]->getDoublePrecisionNonbondedKit());
  }
  int n_unique_tables = 0;
  for (int i = 0; i < topology_count; i++) {
    const int ni_lj_types = nbkvec[i].n_lj_types;
    for (int j = i + 1; j < topology_count; j++) {
      if (ni_lj_types != nbkvec[j].n_lj_types) {
        continue;
      }
      bool diags_identical = true;
      for (int k = 0; k < ni_lj_types; k++) {
        const int dgidx = k * (ni_lj_types + 1);
        diags_identical = (diags_identical &&
                           nbkvec[i].lja_coeff[dgidx] == nbkvec[j].lja_coeff[dgidx] &&
                           nbkvec[i].ljb_coeff[dgidx] == nbkvec[j].ljb_coeff[dgidx] &&
                           nbkvec[i].ljc_coeff[dgidx] == nbkvec[j].ljc_coeff[dgidx] &&
                           nbkvec[i].lja_14_coeff[dgidx] == nbkvec[j].lja_14_coeff[dgidx] &&
                           nbkvec[i].ljb_14_coeff[dgidx] == nbkvec[j].ljb_14_coeff[dgidx] &&
                           nbkvec[i].ljc_14_coeff[dgidx] == nbkvec[j].ljc_14_coeff[dgidx]);
      }
      if (diags_identical == false) {
        continue;
      }
      bool tables_identical = true;
      for (int k = 0; k < ni_lj_types * ni_lj_types; k++) {
        tables_identical = (tables_identical &&
                            nbkvec[i].lja_coeff[k] == nbkvec[j].lja_coeff[k] &&
                            nbkvec[i].ljb_coeff[k] == nbkvec[j].ljb_coeff[k] &&
                            nbkvec[i].ljc_coeff[k] == nbkvec[j].ljc_coeff[k] &&
                            nbkvec[i].lja_14_coeff[k] == nbkvec[j].lja_14_coeff[k] &&
                            nbkvec[i].ljb_14_coeff[k] == nbkvec[j].ljb_14_coeff[k] &&
                            nbkvec[i].ljc_14_coeff[k] == nbkvec[j].ljc_14_coeff[k]);
      }
      if (tables_identical) {
        table_idx[j] = n_unique_tables;      
      }
    }
    table_idx[i] = n_unique_tables;
    n_unique_tables++;
  }

  // Loop back over systems to obtain the necessary allocation base on the unique tables
  int seek_idx = 0;
  int alloc_size = 0;
  std::vector<int> tmp_table_offsets(topology_count, 0);
  for (int i = 0; i < topology_count; i++) {
    for (int j = i; j < topology_count; j++) {
      if (table_idx[i] == seek_idx) {
        tmp_table_offsets[j] = alloc_size;
      }
    }
    if (table_idx[i] == seek_idx) {
      alloc_size += roundUp(nbkvec[i].n_lj_types * nbkvec[i].n_lj_types, warp_size_int);
      seek_idx++;
    }
  }

  // Allocate and fill the Lennard-Jones tables
  lennard_jones_ab_coeff.resize(alloc_size);
  lennard_jones_c_coeff.resize(alloc_size);
  lennard_jones_14_a_coeff.resize(alloc_size);
  lennard_jones_14_b_coeff.resize(alloc_size);
  lennard_jones_14_c_coeff.resize(alloc_size);
  lennard_jones_sigma.resize(alloc_size);
  lennard_jones_14_sigma.resize(alloc_size);
  sp_lennard_jones_ab_coeff.resize(alloc_size);
  sp_lennard_jones_c_coeff.resize(alloc_size);
  sp_lennard_jones_14_a_coeff.resize(alloc_size);
  sp_lennard_jones_14_b_coeff.resize(alloc_size);
  sp_lennard_jones_14_c_coeff.resize(alloc_size);
  sp_lennard_jones_sigma.resize(alloc_size);
  sp_lennard_jones_14_sigma.resize(alloc_size);
  double2* ab_ptr = lennard_jones_ab_coeff.data();
  double* c_ptr = lennard_jones_c_coeff.data();
  double* a_14_ptr = lennard_jones_14_a_coeff.data();
  double* b_14_ptr = lennard_jones_14_b_coeff.data();
  double* c_14_ptr = lennard_jones_14_c_coeff.data();
  double* sigma_ptr = lennard_jones_sigma.data();
  double* sigma_14_ptr = lennard_jones_14_sigma.data();
  float2* sp_ab_ptr = sp_lennard_jones_ab_coeff.data();
  float* sp_c_ptr = sp_lennard_jones_c_coeff.data();
  float* sp_a_14_ptr = sp_lennard_jones_14_a_coeff.data();
  float* sp_b_14_ptr = sp_lennard_jones_14_b_coeff.data();
  float* sp_c_14_ptr = sp_lennard_jones_14_c_coeff.data();
  float* sp_sigma_ptr = sp_lennard_jones_sigma.data();
  float* sp_sigma_14_ptr = sp_lennard_jones_14_sigma.data();
  seek_idx = 0;
  for (int i = 0; i < topology_count; i++) {
    if (table_idx[i] == seek_idx) {
      const int offset = tmp_table_offsets[i];
      for (int j = 0; j < nbkvec[i].n_lj_types * nbkvec[i].n_lj_types; j++) {
        const int joffset = offset + j;
        ab_ptr[joffset] = { nbkvec[i].lja_coeff[j], nbkvec[i].ljb_coeff[j] };
        c_ptr[joffset] = nbkvec[i].ljc_coeff[j];
        a_14_ptr[joffset] = nbkvec[i].lja_14_coeff[j];
        b_14_ptr[joffset] = nbkvec[i].ljb_14_coeff[j];
        c_14_ptr[joffset] = nbkvec[i].ljc_14_coeff[j];
        sigma_ptr[joffset] = nbkvec[i].lj_sigma[j];
        sigma_14_ptr[joffset] = nbkvec[i].lj_14_sigma[j];
        sp_ab_ptr[joffset] = { static_cast<float>(nbkvec[i].lja_coeff[j]),
                               static_cast<float>(nbkvec[i].ljb_coeff[j]) };
        sp_c_ptr[joffset] = nbkvec[i].ljc_coeff[j];
        sp_a_14_ptr[joffset] = nbkvec[i].lja_14_coeff[j];
        sp_b_14_ptr[joffset] = nbkvec[i].ljb_14_coeff[j];
        sp_c_14_ptr[joffset] = nbkvec[i].ljc_14_coeff[j];
        sp_sigma_ptr[joffset] = nbkvec[i].lj_sigma[j];
        sp_sigma_14_ptr[joffset] = nbkvec[i].lj_14_sigma[j];
      }
      seek_idx++;
    }
  }

  // Expand the topology-oriented offsets to become system-oriented offsets
  std::vector<int> tmp_lj_system_offsets(system_count, 0);
  for (int i = 0; i < system_count; i++) {
    tmp_lj_system_offsets[i] = tmp_table_offsets[topology_indices.readHost(i)];
  }
  lennard_jones_abc_offsets.putHost(tmp_lj_system_offsets);
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::mapUniqueRestraintKRSeries(const int order,
                                                   const std::vector<int> &network_table_offsets,
                                                   std::vector<int> *synthesis_index,
                                                   std::vector<int2> *filtered_step_bounds,
                                                   std::vector<double2> *filtered_init_keq,
                                                   std::vector<double2> *filtered_finl_keq,
                                                   std::vector<double4> *filtered_init_r,
                                                   std::vector<double4> *filtered_finl_r) {
  int n_unique_term = 0;
  int* synthesis_index_ptr = synthesis_index->data();
  for (int i = 0; i < restraint_network_count; i++) {
    const RestraintApparatus* ira_ptr = restraint_networks[i];
    if (ira_ptr == nullptr) {
      continue;
    }
    const RestraintStage istage = RestraintStage::INITIAL;
    const RestraintStage fstage = RestraintStage::FINAL;
    const int* irstr_init_step    = ira_ptr->getApplicationStepPointer(order, istage);
    const int* irstr_finl_step    = ira_ptr->getApplicationStepPointer(order, fstage);
    const double2* irstr_init_keq = ira_ptr->getHarmonicStiffnessPointer(order, istage);
    const double2* irstr_finl_keq = ira_ptr->getHarmonicStiffnessPointer(order, fstage);
    const double4* irstr_init_r   = ira_ptr->getDisplacementPointer(order, istage);
    const double4* irstr_finl_r   = ira_ptr->getDisplacementPointer(order, fstage);
    int jmax;
    switch (order) {
    case 1:
      jmax = ira_ptr->getPositionalRestraintCount();
      break;
    case 2:
      jmax = ira_ptr->getDistanceRestraintCount();
      break;
    case 3:
      jmax = ira_ptr->getAngleRestraintCount();
      break;
    case 4:
      jmax = ira_ptr->getDihedralRestraintCount();
      break;
    default:
      break;
    }
    
    // Seek out unique positional restraints.  This process differs from filtering unique valence
    // parameters because with the restraints, there is just a list of parameters and atoms to
    // which they apply, whereas with the valence parameters there were lists of atoms to which
    // each term applied, a corresponding index into a parameter table, and then the separate
    // parameter tables.  This process, in contrast, must just loop over all restraints, rather
    // than all entries in the tables, to determine which have unique parameters.  In many cases,
    // all restraints will be unique, but the instructions for doing the restraints will have atoms
    // and a restraint parameter index regardless.
    for (int j = 0; j < jmax; j++) {

      // Time dependence, and the exact steps at which that time dependence starts, can signify
      // that both k / r series as well as x / y / z targets are unique.
      const int ij_init_step = irstr_init_step[j];
      const int ij_finl_step = irstr_finl_step[j];
      const bool ij_time_dep = (ij_finl_step == 0);
      
      // Skip restraints that have already been determined to have a known k / r series
      if (synthesis_index_ptr[network_table_offsets[i] + j] < 0) {
        const Approx ij_posn_init_k2(irstr_init_keq[j].x, constants::verytiny);
        const Approx ij_posn_init_k3(irstr_init_keq[j].y, constants::verytiny);
        const Approx ij_posn_finl_k2(irstr_finl_keq[j].x, constants::verytiny);
        const Approx ij_posn_finl_k3(irstr_finl_keq[j].y, constants::verytiny);
        const Approx ij_posn_init_r1(irstr_init_r[j].x, constants::verytiny);
        const Approx ij_posn_init_r2(irstr_init_r[j].y, constants::verytiny);
        const Approx ij_posn_init_r3(irstr_init_r[j].z, constants::verytiny);
        const Approx ij_posn_init_r4(irstr_init_r[j].w, constants::verytiny);
        const Approx ij_posn_finl_r1(irstr_finl_r[j].x, constants::verytiny);
        const Approx ij_posn_finl_r2(irstr_finl_r[j].y, constants::verytiny);
        const Approx ij_posn_finl_r3(irstr_finl_r[j].z, constants::verytiny);
        const Approx ij_posn_finl_r4(irstr_finl_r[j].w, constants::verytiny);
        for (int k = i; k < restraint_network_count; k++) {
          const RestraintApparatus* kra_ptr = restraint_networks[k];
          if (kra_ptr == nullptr) {
            continue;
          }

          // Another switch to define constants and pointers based on the order of restraint
          // being handled, this time for the abstract of the kth restraint apparatus
          const int* krstr_init_step    = kra_ptr->getApplicationStepPointer(order, istage);
          const int* krstr_finl_step    = kra_ptr->getApplicationStepPointer(order, fstage);
          const double2* krstr_init_keq = kra_ptr->getHarmonicStiffnessPointer(order, istage);
          const double2* krstr_finl_keq = kra_ptr->getHarmonicStiffnessPointer(order, fstage);
          const double4* krstr_init_r   = kra_ptr->getDisplacementPointer(order, istage);
          const double4* krstr_finl_r   = kra_ptr->getDisplacementPointer(order, fstage);          
          int mmax;
          switch (order) {
          case 1:
            mmax = kra_ptr->getPositionalRestraintCount();
            break;
          case 2:
            mmax = kra_ptr->getDistanceRestraintCount();
            break;
          case 3:
            mmax = kra_ptr->getAngleRestraintCount();
            break;
          case 4:
            mmax = kra_ptr->getDihedralRestraintCount();
            break;
          default:
            break;
          }          
          const int mstart = (k == i) ? j : 0;
          for (int m = mstart; m < mmax; m++) {
            if (synthesis_index_ptr[network_table_offsets[k] + m] >= 0) {
              continue;
            }
            const bool km_time_dep = (krstr_finl_step[m] == 0);
            if (ij_time_dep != km_time_dep) {
              continue;
            }
            if (ij_time_dep) {
              if (krstr_init_step[m] == ij_init_step &&
                  krstr_finl_step[m] == ij_finl_step &&
                  ij_posn_init_k2.test(krstr_init_keq[m].x) &&
                  ij_posn_init_k3.test(krstr_init_keq[m].y) &&
                  ij_posn_finl_k2.test(krstr_finl_keq[m].x) &&
                  ij_posn_finl_k3.test(krstr_finl_keq[m].y) &&
                  ij_posn_init_r1.test(krstr_init_r[m].x) &&
                  ij_posn_init_r2.test(krstr_init_r[m].y) &&
                  ij_posn_init_r3.test(krstr_init_r[m].z) &&
                  ij_posn_init_r4.test(krstr_init_r[m].w) &&
                  ij_posn_finl_r1.test(krstr_finl_r[m].x) &&
                  ij_posn_finl_r2.test(krstr_finl_r[m].y) &&
                  ij_posn_finl_r3.test(krstr_finl_r[m].z) &&
                  ij_posn_finl_r4.test(krstr_finl_r[m].w)) {
                synthesis_index_ptr[network_table_offsets[k] + m] = n_unique_term;
              }
            }
            else {
              if (ij_posn_init_k2.test(krstr_init_keq[m].x) &&
                  ij_posn_init_k3.test(krstr_init_keq[m].y) &&
                  ij_posn_init_r1.test(krstr_init_r[m].x) &&
                  ij_posn_init_r2.test(krstr_init_r[m].y) &&
                  ij_posn_init_r3.test(krstr_init_r[m].z) &&
                  ij_posn_init_r4.test(krstr_init_r[m].w)) {
                synthesis_index_ptr[network_table_offsets[k] + m] = n_unique_term;
              }
            }
          }
        }

        // Catalog this unique positional restraint parameter set (the x, y, and z coordinates
        // are a separate parameter set and will be considered next);
        filtered_step_bounds->push_back({ij_init_step, ij_finl_step});
        filtered_init_keq->push_back(irstr_init_keq[j]);
        filtered_finl_keq->push_back(irstr_finl_keq[j]);
        filtered_init_r->push_back(irstr_init_r[j]);
        filtered_finl_r->push_back(irstr_finl_r[j]);
        n_unique_term++;
      }
    }
  }
  return n_unique_term;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::condenseRestraintNetworks() {

  // As with topology parameters, compute the numbers of unique positional, distance, angle, and
  // dihedral angle restraints.  Take the opportunity to compute offsets (starting bounds) for
  // various sets of terms.  One difference in this implementation is that the restraint apparatus
  // pointer for any particular network may be the null pointer, indicating that there is no
  // restraint system in place for that case (or perhaps any cases).
  int rposn_offset = 0;
  int rbond_offset = 0;
  int rangl_offset = 0;
  int rdihe_offset = 0;
  std::vector<int2> network_rposn_param_map_bounds(restraint_network_count);
  std::vector<int2> network_rbond_param_map_bounds(restraint_network_count);
  std::vector<int2> network_rangl_param_map_bounds(restraint_network_count);
  std::vector<int2> network_rdihe_param_map_bounds(restraint_network_count);
  std::vector<int> network_rposn_table_offsets(restraint_network_count, 0);
  std::vector<int> network_rbond_table_offsets(restraint_network_count, 0);
  std::vector<int> network_rangl_table_offsets(restraint_network_count, 0);
  std::vector<int> network_rdihe_table_offsets(restraint_network_count, 0);
  for (int i = 0; i < restraint_network_count; i++) {
    const RestraintApparatus* ra_ptr = restraint_networks[i];
    network_rposn_table_offsets[i] = rposn_offset;
    network_rbond_table_offsets[i] = rbond_offset;
    network_rangl_table_offsets[i] = rangl_offset;
    network_rdihe_table_offsets[i] = rdihe_offset;
    network_rposn_param_map_bounds[i].x = rposn_offset;
    network_rbond_param_map_bounds[i].x = rbond_offset;
    network_rangl_param_map_bounds[i].x = rangl_offset;
    network_rdihe_param_map_bounds[i].x = rdihe_offset;
    if (ra_ptr != nullptr) {
      network_rposn_param_map_bounds[i].y = rposn_offset + ra_ptr->getPositionalRestraintCount();
      network_rbond_param_map_bounds[i].y = rbond_offset + ra_ptr->getDistanceRestraintCount();
      network_rangl_param_map_bounds[i].y = rangl_offset + ra_ptr->getAngleRestraintCount();
      network_rdihe_param_map_bounds[i].y = rdihe_offset + ra_ptr->getDihedralRestraintCount();
      rposn_offset += roundUp(ra_ptr->getPositionalRestraintCount(), warp_size_int);
      rbond_offset += roundUp(ra_ptr->getDistanceRestraintCount(), warp_size_int);
      rangl_offset += roundUp(ra_ptr->getAngleRestraintCount(), warp_size_int);
      rdihe_offset += roundUp(ra_ptr->getDihedralRestraintCount(), warp_size_int);
    }
    else {
      network_rposn_param_map_bounds[i].y = rposn_offset;
      network_rbond_param_map_bounds[i].y = rbond_offset;
      network_rangl_param_map_bounds[i].y = rangl_offset;
      network_rdihe_param_map_bounds[i].y = rdihe_offset;
    }
  }

  // Create lists of unique parameters for the various types of restraints.
  std::vector<int> rposn_synthesis_kr_index(rposn_offset, -1);  
  std::vector<int> rposn_synthesis_xyz_index(rposn_offset, -1);  
  std::vector<int> rbond_synthesis_index(rbond_offset, -1);
  std::vector<int> rangl_synthesis_index(rangl_offset, -1);
  std::vector<int> rdihe_synthesis_index(rdihe_offset, -1);
  std::vector<int2> filtered_rposn_step_bounds, filtered_rbond_step_bounds;
  std::vector<int2> filtered_rangl_step_bounds, filtered_rdihe_step_bounds;
  std::vector<double2> filtered_rposn_init_keq, filtered_rposn_finl_keq;
  std::vector<double4> filtered_rposn_init_r, filtered_rposn_finl_r;
  std::vector<double2> filtered_rposn_init_xy, filtered_rposn_finl_xy;
  std::vector<double> filtered_rposn_init_z, filtered_rposn_finl_z;
  std::vector<double2> filtered_rbond_init_keq, filtered_rbond_finl_keq;
  std::vector<double4> filtered_rbond_init_r, filtered_rbond_finl_r;
  std::vector<double2> filtered_rangl_init_keq, filtered_rangl_finl_keq;
  std::vector<double4> filtered_rangl_init_r, filtered_rangl_finl_r;
  std::vector<double2> filtered_rdihe_init_keq, filtered_rdihe_finl_keq;
  std::vector<double4> filtered_rdihe_init_r, filtered_rdihe_finl_r;
  const int n_unique_posn_kr = mapUniqueRestraintKRSeries(1, network_rposn_table_offsets,
                                                          &rposn_synthesis_kr_index,
                                                          &filtered_rposn_step_bounds,
                                                          &filtered_rposn_init_keq,
                                                          &filtered_rposn_finl_keq,
                                                          &filtered_rposn_init_r,
                                                          &filtered_rposn_finl_r);
  const int n_unique_bond    = mapUniqueRestraintKRSeries(2, network_rbond_table_offsets,
                                                          &rbond_synthesis_index,
                                                          &filtered_rbond_step_bounds,
                                                          &filtered_rbond_init_keq,
                                                          &filtered_rbond_finl_keq,
                                                          &filtered_rbond_init_r,
                                                          &filtered_rbond_finl_r);
  const int n_unique_angl    = mapUniqueRestraintKRSeries(3, network_rangl_table_offsets,
                                                          &rangl_synthesis_index,
                                                          &filtered_rangl_step_bounds,
                                                          &filtered_rangl_init_keq,
                                                          &filtered_rangl_finl_keq,
                                                          &filtered_rangl_init_r,
                                                          &filtered_rangl_finl_r);
  const int n_unique_dihe    = mapUniqueRestraintKRSeries(4, network_rdihe_table_offsets,
                                                          &rdihe_synthesis_index,
                                                          &filtered_rdihe_step_bounds,
                                                          &filtered_rdihe_init_keq,
                                                          &filtered_rdihe_finl_keq,
                                                          &filtered_rdihe_init_r,
                                                          &filtered_rdihe_finl_r);

  // Seek out unique positional restraint x, y, and z targets.  This process differs from
  // filtering unique valence parameters because with the restraints, there is just a list of
  // parameters and atoms to which they apply, whereas with the valence parameters there were
  // lists of atoms to which each term applied, a corresponding index into a parameter table,
  // and then the separate parameter tables.  This process, in contrast, must just loop over all
  // restraints, rather than all entries in the tables, to determine which have unique
  // parameters.  In many cases, all restraints will be unique, but the instructions for doing
  // the restraints will have atoms and a restraint parameter index regardless.
  int n_unique_posn_xyz = 0;
  for (int i = 0; i < restraint_network_count; i++) {
    const RestraintApparatus* ira_ptr = restraint_networks[i];
    if (ira_ptr == nullptr) {
      continue;
    }
    const RestraintKit<double, double2, double4> irar_dp = ira_ptr->dpData();
    for (int j = 0; j < irar_dp.nposn; j++) {

      // Skip restraints that have been determined to have unique x / y / z targets.  This is
      // better done by checking whether the x / y / z target has not yet been matched to some
      // previous target, as the k / r combination can also be unique.
      if (rposn_synthesis_xyz_index[network_rposn_table_offsets[i] + j] < 0) {
        const int ij_init_step = irar_dp.rposn_init_step[j];
        const int ij_finl_step = irar_dp.rposn_finl_step[j];
        const bool ij_time_dep = (ij_finl_step == 0);
        const Approx ij_posn_init_x(irar_dp.rposn_init_xy[j].x, constants::verytiny);
        const Approx ij_posn_init_y(irar_dp.rposn_init_xy[j].y, constants::verytiny);
        const Approx ij_posn_init_z(irar_dp.rposn_init_z[j],    constants::verytiny);
        const Approx ij_posn_finl_x(irar_dp.rposn_finl_xy[j].x, constants::verytiny);
        const Approx ij_posn_finl_y(irar_dp.rposn_finl_xy[j].y, constants::verytiny);
        const Approx ij_posn_finl_z(irar_dp.rposn_finl_z[j],    constants::verytiny);
        for (int k = i; k < restraint_network_count; k++) {
          const RestraintApparatus* kra_ptr = restraint_networks[k];
          if (kra_ptr == nullptr) {
            continue;
          }
          const RestraintKit<double,
                             double2, double4> krar_dp = kra_ptr->dpData();
          const int mstart = (k == i) ? j : 0;
          for (int m = mstart; m < krar_dp.nposn; m++) {
            if (rposn_synthesis_xyz_index[network_rposn_table_offsets[k] + m] >= 0) {
              continue;
            }
            const bool km_time_dep = (krar_dp.rposn_finl_step[m] == 0);
            if (ij_time_dep != km_time_dep) {
              continue;
            }
            if (ij_time_dep) {
              if (ij_posn_init_x.test(krar_dp.rposn_init_xy[m].x) &&
                  ij_posn_init_y.test(krar_dp.rposn_init_xy[m].y) &&
                  ij_posn_init_z.test(krar_dp.rposn_init_z[m]) &&
                  ij_posn_finl_x.test(krar_dp.rposn_finl_xy[m].x) &&
                  ij_posn_finl_y.test(krar_dp.rposn_finl_xy[m].y) &&
                  ij_posn_finl_z.test(krar_dp.rposn_finl_z[m])) {
                rposn_synthesis_xyz_index[network_rposn_table_offsets[k] + m] = n_unique_posn_xyz;
              }
            }
            else {
              if (ij_posn_init_x.test(krar_dp.rposn_init_xy[m].x) &&
                  ij_posn_init_y.test(krar_dp.rposn_init_xy[m].y) &&
                  ij_posn_init_z.test(krar_dp.rposn_init_z[m])) {
                rposn_synthesis_xyz_index[network_rposn_table_offsets[k] + m] = n_unique_posn_xyz;
              }
            }
          }
        }

        // Catalog this unique positional restraint target coordinate set
        filtered_rposn_init_xy.push_back({irar_dp.rposn_init_xy[j].x, irar_dp.rposn_init_xy[j].y});
        filtered_rposn_init_z.push_back(irar_dp.rposn_init_z[j]);
        filtered_rposn_finl_xy.push_back({irar_dp.rposn_finl_xy[j].x, irar_dp.rposn_finl_xy[j].y});
        filtered_rposn_finl_z.push_back(irar_dp.rposn_finl_z[j]);
        n_unique_posn_xyz++;
      }
    }
  }

  // Resize the Hybrid arrays, set pointers, and fill Hybrid objects in the AtomGraphSynthesis
  const int rnc_offset = roundUp(restraint_network_count, warp_size_int);
  const size_t i2c = roundUp(n_unique_posn_kr, warp_size_int) +
                     roundUp(n_unique_bond, warp_size_int) +
                     roundUp(n_unique_angl, warp_size_int) +
                     roundUp(n_unique_dihe, warp_size_int) + (4 * rnc_offset);
  nmr_int2_data.resize(i2c);
  nmr_double_data.resize(2 * roundUp(n_unique_posn_xyz, warp_size_int));
  nmr_double2_data.resize(2 * (i2c + n_unique_posn_xyz));
  nmr_double4_data.resize(2 * i2c);
  nmr_float_data.resize(2 * roundUp(n_unique_posn_xyz, warp_size_int));
  nmr_float2_data.resize(2 * (i2c + n_unique_posn_xyz));
  nmr_float4_data.resize(2 * i2c);
  size_t ic = 0LLU;
  ic = rposn_step_bounds.putHost(&nmr_int2_data, filtered_rposn_step_bounds, ic, warp_size_zu);
  ic = rbond_step_bounds.putHost(&nmr_int2_data, filtered_rbond_step_bounds, ic, warp_size_zu);
  ic = rangl_step_bounds.putHost(&nmr_int2_data, filtered_rangl_step_bounds, ic, warp_size_zu);
  ic = rdihe_step_bounds.putHost(&nmr_int2_data, filtered_rdihe_step_bounds, ic, warp_size_zu);
  ic = rposn_param_map_bounds.putHost(&nmr_int2_data, network_rposn_param_map_bounds, ic,
                                      warp_size_zu);
  ic = rbond_param_map_bounds.putHost(&nmr_int2_data, network_rbond_param_map_bounds, ic,
                                      warp_size_zu);
  ic = rangl_param_map_bounds.putHost(&nmr_int2_data, network_rangl_param_map_bounds, ic,
                                      warp_size_zu);
  ic = rdihe_param_map_bounds.putHost(&nmr_int2_data, network_rdihe_param_map_bounds, ic,
                                      warp_size_zu);
  ic = 0LLU;
  ic = rposn_init_z.putHost(&nmr_double_data, filtered_rposn_init_z, ic, warp_size_zu);
  ic = rposn_final_z.putHost(&nmr_double_data, filtered_rposn_finl_z, ic, warp_size_zu);
  ic = 0LLU;
  ic = rposn_init_k.putHost(&nmr_double2_data, filtered_rposn_init_keq, ic, warp_size_zu);
  ic = rposn_final_k.putHost(&nmr_double2_data, filtered_rposn_finl_keq, ic, warp_size_zu);
  ic = rposn_init_xy.putHost(&nmr_double2_data, filtered_rposn_init_xy, ic, warp_size_zu);
  ic = rposn_final_xy.putHost(&nmr_double2_data, filtered_rposn_finl_xy, ic, warp_size_zu);
  ic = rbond_init_k.putHost(&nmr_double2_data, filtered_rbond_init_keq, ic, warp_size_zu);
  ic = rbond_final_k.putHost(&nmr_double2_data, filtered_rbond_finl_keq, ic, warp_size_zu);
  ic = rangl_init_k.putHost(&nmr_double2_data, filtered_rangl_init_keq, ic, warp_size_zu);
  ic = rangl_final_k.putHost(&nmr_double2_data, filtered_rangl_finl_keq, ic, warp_size_zu);
  ic = rdihe_init_k.putHost(&nmr_double2_data, filtered_rdihe_init_keq, ic, warp_size_zu);
  ic = rdihe_final_k.putHost(&nmr_double2_data, filtered_rdihe_finl_keq, ic, warp_size_zu);
  ic = 0LLU;
  ic = rposn_init_r.putHost(&nmr_double4_data, filtered_rposn_init_r, ic, warp_size_zu);
  ic = rposn_final_r.putHost(&nmr_double4_data, filtered_rposn_finl_r, ic, warp_size_zu);  
  ic = rbond_init_r.putHost(&nmr_double4_data, filtered_rbond_init_r, ic, warp_size_zu);
  ic = rbond_final_r.putHost(&nmr_double4_data, filtered_rbond_finl_r, ic, warp_size_zu);
  ic = rangl_init_r.putHost(&nmr_double4_data, filtered_rangl_init_r, ic, warp_size_zu);
  ic = rangl_final_r.putHost(&nmr_double4_data, filtered_rangl_finl_r, ic, warp_size_zu);
  ic = rdihe_init_r.putHost(&nmr_double4_data, filtered_rdihe_init_r, ic, warp_size_zu);
  ic = rdihe_final_r.putHost(&nmr_double4_data, filtered_rdihe_finl_r, ic, warp_size_zu);

  // The nmr_int_data array has already been resized in condenseParameterTables(), with extra
  // space allocated to hold the parameter maps.  Set those pointers now.
  ic = (3 * rposn_offset) + (3 * rbond_offset) + (4 * rangl_offset) + (5 * rdihe_offset);
  ic = rposn_kr_param_map.putHost(&nmr_int_data, rposn_synthesis_kr_index, ic, warp_size_zu);
  ic = rposn_xyz_param_map.putHost(&nmr_int_data, rposn_synthesis_xyz_index, ic, warp_size_zu);
  ic = rbond_param_map.putHost(&nmr_int_data, rbond_synthesis_index, ic, warp_size_zu);
  ic = rangl_param_map.putHost(&nmr_int_data, rangl_synthesis_index, ic, warp_size_zu);
  ic = rdihe_param_map.putHost(&nmr_int_data, rdihe_synthesis_index, ic, warp_size_zu);
  
  // Create temporary arrays for floating-point data.  This is done more in line with what happens
  // in an AtomGraph, with each POINTER-kind Hybrid object targeting a handful of larger arrays.
  // The likelihood that most simulations will have few to no restraints makes it advantageous to
  // condense these arrays as much as possible.
  const std::vector<float> spfil_rposn_init_z(filtered_rposn_init_z.begin(),
                                              filtered_rposn_init_z.end());
  const std::vector<float> spfil_rposn_finl_z(filtered_rposn_finl_z.begin(),
                                              filtered_rposn_finl_z.end());
  const std::vector<float2> spfil_rposn_init_xy = vtConv2f(filtered_rposn_init_xy);
  const std::vector<float2> spfil_rposn_finl_xy = vtConv2f(filtered_rposn_finl_xy);
  const std::vector<float2> spfil_rposn_init_keq = vtConv2f(filtered_rposn_init_keq);
  const std::vector<float2> spfil_rposn_finl_keq = vtConv2f(filtered_rposn_finl_keq);
  const std::vector<float2> spfil_rbond_init_keq = vtConv2f(filtered_rbond_init_keq);
  const std::vector<float2> spfil_rbond_finl_keq = vtConv2f(filtered_rbond_finl_keq);
  const std::vector<float2> spfil_rangl_init_keq = vtConv2f(filtered_rangl_init_keq);
  const std::vector<float2> spfil_rangl_finl_keq = vtConv2f(filtered_rangl_finl_keq);
  const std::vector<float2> spfil_rdihe_init_keq = vtConv2f(filtered_rdihe_init_keq);
  const std::vector<float2> spfil_rdihe_finl_keq = vtConv2f(filtered_rdihe_finl_keq);
  const std::vector<float4> spfil_rposn_init_r = vtConv4f(filtered_rposn_init_r);
  const std::vector<float4> spfil_rposn_finl_r = vtConv4f(filtered_rposn_finl_r);
  const std::vector<float4> spfil_rbond_init_r = vtConv4f(filtered_rbond_init_r);
  const std::vector<float4> spfil_rbond_finl_r = vtConv4f(filtered_rbond_finl_r);
  const std::vector<float4> spfil_rangl_init_r = vtConv4f(filtered_rangl_init_r);
  const std::vector<float4> spfil_rangl_finl_r = vtConv4f(filtered_rangl_finl_r);
  const std::vector<float4> spfil_rdihe_init_r = vtConv4f(filtered_rdihe_init_r);
  const std::vector<float4> spfil_rdihe_finl_r = vtConv4f(filtered_rdihe_finl_r);
  ic = 0LLU;
  ic = sp_rposn_init_z.putHost(&nmr_float_data, spfil_rposn_init_z, ic, warp_size_zu);
  ic = sp_rposn_final_z.putHost(&nmr_float_data, spfil_rposn_finl_z, ic, warp_size_zu);
  ic = 0LLU;
  ic = sp_rposn_init_k.putHost(&nmr_float2_data, spfil_rposn_init_keq, ic, warp_size_zu);
  ic = sp_rposn_final_k.putHost(&nmr_float2_data, spfil_rposn_finl_keq, ic, warp_size_zu);
  ic = sp_rposn_init_xy.putHost(&nmr_float2_data, spfil_rposn_init_xy, ic, warp_size_zu);
  ic = sp_rposn_final_xy.putHost(&nmr_float2_data, spfil_rposn_finl_xy, ic, warp_size_zu);
  ic = sp_rbond_init_k.putHost(&nmr_float2_data, spfil_rbond_init_keq, ic, warp_size_zu);
  ic = sp_rbond_final_k.putHost(&nmr_float2_data, spfil_rbond_finl_keq, ic, warp_size_zu);
  ic = sp_rangl_init_k.putHost(&nmr_float2_data, spfil_rangl_init_keq, ic, warp_size_zu);
  ic = sp_rangl_final_k.putHost(&nmr_float2_data, spfil_rangl_finl_keq, ic, warp_size_zu);
  ic = sp_rdihe_init_k.putHost(&nmr_float2_data, spfil_rdihe_init_keq, ic, warp_size_zu);
  ic = sp_rdihe_final_k.putHost(&nmr_float2_data, spfil_rdihe_finl_keq, ic, warp_size_zu);
  ic = 0LLU;
  ic = sp_rposn_init_r.putHost(&nmr_float4_data, spfil_rposn_init_r, ic, warp_size_zu);
  ic = sp_rposn_final_r.putHost(&nmr_float4_data, spfil_rposn_finl_r, ic, warp_size_zu);  
  ic = sp_rbond_init_r.putHost(&nmr_float4_data, spfil_rbond_init_r, ic, warp_size_zu);
  ic = sp_rbond_final_r.putHost(&nmr_float4_data, spfil_rbond_finl_r, ic, warp_size_zu);
  ic = sp_rangl_init_r.putHost(&nmr_float4_data, spfil_rangl_init_r, ic, warp_size_zu);
  ic = sp_rangl_final_r.putHost(&nmr_float4_data, spfil_rangl_finl_r, ic, warp_size_zu);
  ic = sp_rdihe_init_r.putHost(&nmr_float4_data, spfil_rdihe_init_r, ic, warp_size_zu);
  ic = sp_rdihe_final_r.putHost(&nmr_float4_data, spfil_rdihe_finl_r, ic, warp_size_zu);
  
  // With the restraint parameter tables assembled, mark the restraints for each system in
  // terms of the synthesis tables.
  for (int sysid = 0; sysid < system_count; sysid++) {
    const int ra_index = restraint_indices.readHost(sysid);
    const RestraintApparatus *ra_ptr = restraint_networks[ra_index];
    if (ra_ptr == nullptr) {
      continue;
    }
    const RestraintKit<double, double2, double4> rar = ra_ptr->dpData();
    const int ra_posn_table_offset = network_rposn_table_offsets[ra_index];
    const int ra_bond_table_offset = network_rbond_table_offsets[ra_index];
    const int ra_angl_table_offset = network_rangl_table_offsets[ra_index];
    const int ra_dihe_table_offset = network_rdihe_table_offsets[ra_index];
    const int synth_posn_table_offset = posn_restraint_offsets.readHost(sysid);
    const int synth_bond_table_offset = bond_restraint_offsets.readHost(sysid);
    const int synth_angl_table_offset = angl_restraint_offsets.readHost(sysid);
    const int synth_dihe_table_offset = dihe_restraint_offsets.readHost(sysid);
    for (int j = 0; j < rar.nposn; j++) {
      rposn_kr_param_idx.putHost(rposn_synthesis_kr_index[ra_posn_table_offset + j],
                                 synth_posn_table_offset + j);
      rposn_xyz_param_idx.putHost(rposn_synthesis_xyz_index[ra_posn_table_offset + j],
                                  synth_posn_table_offset + j);
    }
    for (int j = 0; j < rar.nbond; j++) {
      rbond_param_idx.putHost(rbond_synthesis_index[ra_bond_table_offset + j],
                              synth_bond_table_offset + j);
    }
    for (int j = 0; j < rar.nangl; j++) {
      rangl_param_idx.putHost(rangl_synthesis_index[ra_angl_table_offset + j],
                              synth_angl_table_offset + j);
    }
    for (int j = 0; j < rar.ndihe; j++) {
      rdihe_param_idx.putHost(rdihe_synthesis_index[ra_dihe_table_offset + j],
                              synth_dihe_table_offset + j);
    }
  }
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::setVwuAbstractLimits(const int item_counter, const int vwu_counter,
                                             const VwuAbstractMap slot, const int item_quantity) {
  int2* vwu_insr_ptr = vwu_instruction_sets.data();
  const int vv_offset = vwu_counter * vwu_abstract_length + static_cast<int>(slot);
  vwu_insr_ptr[vv_offset].x = item_counter;
  vwu_insr_ptr[vv_offset].y = item_counter + item_quantity;
  return roundUp(item_counter + item_quantity, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::loadValenceWorkUnits(const int vwu_atom_limit) {

  // Loop over all systems and create lists of valence work units
  std::vector<std::vector<ValenceWorkUnit>> all_vwu;
  all_vwu.reserve(system_count);
  std::vector<int> total_task_counts(static_cast<size_t>(VwuTask::ALL_TASKS), 0);
  int total_vwu = 0;
  for (int i = 0; i < system_count; i++) {
    const AtomGraph *ag = topologies[topology_indices.readHost(i)];

    // Check the restraint network--if the restraint network is blank then construct a 
    // blank network.
    const int rn = restraint_indices.readHost(i);
    const RestraintApparatus ra_blank = RestraintApparatus(ag);
    const RestraintApparatus *ra = (restraint_networks[rn] == nullptr) ? &ra_blank :
                                                                         restraint_networks[rn];
    all_vwu.push_back(buildValenceWorkUnits(ag, ra, vwu_atom_limit));
    total_vwu += all_vwu[i].size();    
  }

  // Allocate the valence work unit instruction set array, an abstract of each ValenceWorkUnit
  // in a stream of integers.
  total_valence_work_units = total_vwu;
  vwu_instruction_sets.resize(total_vwu * vwu_abstract_length);
  int2* vwu_insr_ptr = vwu_instruction_sets.data();
  int atom_import_count = 0;
  int manipulation_mask_count = 0;
  int cbnd_count = 0;
  int angl_count = 0;
  int cdhe_count = 0;
  int cmap_count = 0;
  int infr14_count = 0;
  int rposn_count = 0;
  int rbond_count = 0;
  int rangl_count = 0;
  int rdihe_count = 0;
  int vsite_count = 0;
  int settle_count = 0;
  int cgroup_count = 0;
  int cbnd_nrg_count = 0;
  int angl_nrg_count = 0;
  int cdhe_nrg_count = 0;
  int cmap_nrg_count = 0;
  int infr14_nrg_count = 0;
  int rposn_nrg_count = 0;
  int rbond_nrg_count = 0;
  int rangl_nrg_count = 0;
  int rdihe_nrg_count = 0;
  int vwu_count = 0;
  for (int i = 0; i < system_count; i++) {
    const size_t nvwu = all_vwu[i].size(); 
    for (size_t j = 0LLU; j < nvwu; j++) {
      const std::vector<int> task_counts = all_vwu[i][j].getTaskCounts();
      
      // Record the appropriate limits for the present system
      const int ntatom = all_vwu[i][j].getImportedAtomCount();
      const int manip_uints = (ntatom + uint_bit_count_int - 1) / uint_bit_count_int;
      atom_import_count = setVwuAbstractLimits(atom_import_count, vwu_count,
                                               VwuAbstractMap::IMPORT, ntatom);
      manipulation_mask_count = setVwuAbstractLimits(manipulation_mask_count, vwu_count,
                                                     VwuAbstractMap::MANIPULATE, manip_uints);
      const int ntcbnd  = task_counts[static_cast<size_t>(VwuTask::CBND)];
      const int ntangl  = task_counts[static_cast<size_t>(VwuTask::ANGL)];
      const int ntcdhe  = task_counts[static_cast<size_t>(VwuTask::CDHE)];
      const int ntcmap  = task_counts[static_cast<size_t>(VwuTask::CMAP)];
      const int ntinfr  = task_counts[static_cast<size_t>(VwuTask::INFR14)];
      const int ntrposn = task_counts[static_cast<size_t>(VwuTask::RPOSN)];
      const int ntrbond = task_counts[static_cast<size_t>(VwuTask::RBOND)];
      const int ntrangl = task_counts[static_cast<size_t>(VwuTask::RANGL)];
      const int ntrdihe = task_counts[static_cast<size_t>(VwuTask::RDIHE)];
      const int ntvste  = task_counts[static_cast<size_t>(VwuTask::VSITE)];
      const int ntsett  = task_counts[static_cast<size_t>(VwuTask::SETTLE)];
      const int ntcnst  = task_counts[static_cast<size_t>(VwuTask::CGROUP)];
      cbnd_count   = setVwuAbstractLimits(cbnd_count, vwu_count, VwuAbstractMap::CBND, ntcbnd);
      angl_count   = setVwuAbstractLimits(angl_count, vwu_count, VwuAbstractMap::ANGL, ntangl);
      cdhe_count   = setVwuAbstractLimits(cdhe_count, vwu_count, VwuAbstractMap::CDHE, ntcdhe);
      cmap_count   = setVwuAbstractLimits(cmap_count, vwu_count, VwuAbstractMap::CMAP, ntcmap);
      infr14_count = setVwuAbstractLimits(infr14_count, vwu_count, VwuAbstractMap::INFR14, ntinfr);
      rposn_count  = setVwuAbstractLimits(rposn_count, vwu_count, VwuAbstractMap::RPOSN, ntrposn);
      rbond_count  = setVwuAbstractLimits(rbond_count, vwu_count, VwuAbstractMap::RBOND, ntrbond);
      rangl_count  = setVwuAbstractLimits(rangl_count, vwu_count, VwuAbstractMap::RANGL, ntrangl);
      rdihe_count  = setVwuAbstractLimits(rdihe_count, vwu_count, VwuAbstractMap::RDIHE, ntrdihe);
      vsite_count  = setVwuAbstractLimits(vsite_count, vwu_count, VwuAbstractMap::VSITE, ntvste);
      settle_count = setVwuAbstractLimits(settle_count, vwu_count, VwuAbstractMap::SETTLE, ntsett);
      cgroup_count = setVwuAbstractLimits(cgroup_count, vwu_count, VwuAbstractMap::CGROUP, ntcnst);
      const int cbnd_uints   = (ntcbnd + uint_bit_count_int - 1) / uint_bit_count_int;
      const int angl_uints   = (ntangl + uint_bit_count_int - 1) / uint_bit_count_int;
      const int cdhe_uints   = (ntcdhe + uint_bit_count_int - 1) / uint_bit_count_int;
      const int cmap_uints   = (ntcmap + uint_bit_count_int - 1) / uint_bit_count_int;
      const int infr14_uints = (ntinfr + uint_bit_count_int - 1) / uint_bit_count_int;
      const int rposn_uints  = (ntrposn + uint_bit_count_int - 1) / uint_bit_count_int;
      const int rbond_uints  = (ntrbond + uint_bit_count_int - 1) / uint_bit_count_int;
      const int rangl_uints  = (ntrangl + uint_bit_count_int - 1) / uint_bit_count_int;
      const int rdihe_uints  = (ntrdihe + uint_bit_count_int - 1) / uint_bit_count_int;
      cbnd_nrg_count   = setVwuAbstractLimits(cbnd_nrg_count, vwu_count, VwuAbstractMap::CBND_NRG,
                                              cbnd_uints);
      angl_nrg_count   = setVwuAbstractLimits(angl_nrg_count, vwu_count, VwuAbstractMap::ANGL_NRG,
                                              angl_uints);
      cdhe_nrg_count   = setVwuAbstractLimits(cdhe_nrg_count, vwu_count, VwuAbstractMap::CDHE_NRG,
                                              cdhe_uints);
      cmap_nrg_count   = setVwuAbstractLimits(cmap_nrg_count, vwu_count, VwuAbstractMap::CMAP_NRG,
                                              cmap_uints);
      infr14_nrg_count = setVwuAbstractLimits(infr14_nrg_count, vwu_count,
                                              VwuAbstractMap::INFR14_NRG, infr14_uints);
      rposn_nrg_count  = setVwuAbstractLimits(rposn_nrg_count, vwu_count,
                                              VwuAbstractMap::RPOSN_NRG, rposn_uints);
      rbond_nrg_count  = setVwuAbstractLimits(rbond_nrg_count, vwu_count,
                                              VwuAbstractMap::RBOND_NRG, rbond_uints);
      rangl_nrg_count  = setVwuAbstractLimits(rangl_nrg_count, vwu_count,
                                              VwuAbstractMap::RANGL_NRG, rangl_uints);
      rdihe_nrg_count  = setVwuAbstractLimits(rdihe_nrg_count, vwu_count,
                                              VwuAbstractMap::RDIHE_NRG, rdihe_uints);
      const size_t klim = static_cast<size_t>(VwuTask::ALL_TASKS);
      for (size_t k = 0LLU; k < klim; k++) {
        total_task_counts[k] += roundUp(task_counts[k], warp_size_int);        
      }
      setVwuAbstractLimits(i, vwu_count, VwuAbstractMap::SYSTEM_ID, 0);
      vwu_count++;
    }
  }
  
  // Allocate instruction arrays in the topological synthesis
  vwu_import_lists.resize(atom_import_count);
  insr_uint_data.resize(cdhe_count + infr14_count + cbnd_nrg_count + angl_nrg_count +
                        cdhe_nrg_count + cmap_nrg_count + infr14_nrg_count + rposn_nrg_count +
                        rbond_nrg_count + rangl_nrg_count + rdihe_nrg_count);
  insr_uint2_data.resize(manipulation_mask_count + cbnd_count + angl_count + cdhe_count +
                         cmap_count + rposn_count + rbond_count + rangl_count + rdihe_count +
                         vsite_count + settle_count + cgroup_count);
  size_t pivot = 0LLU;
  cdhe_overtones.setPointer(&insr_uint_data, pivot, cdhe_count);
  pivot += cdhe_count;
  infr14_instructions.setPointer(&insr_uint_data, pivot, infr14_count);
  pivot += infr14_count;
  accumulate_cbnd_energy.setPointer(&insr_uint_data, pivot, cbnd_nrg_count);
  pivot += cbnd_nrg_count;
  accumulate_angl_energy.setPointer(&insr_uint_data, pivot, angl_nrg_count);
  pivot += angl_nrg_count;
  accumulate_cdhe_energy.setPointer(&insr_uint_data, pivot, cdhe_nrg_count);
  pivot += cdhe_nrg_count;
  accumulate_cmap_energy.setPointer(&insr_uint_data, pivot, cmap_nrg_count);
  pivot += cmap_nrg_count;
  accumulate_infr14_energy.setPointer(&insr_uint_data, pivot, infr14_nrg_count);
  pivot += infr14_nrg_count;
  accumulate_rposn_energy.setPointer(&insr_uint_data, pivot, rposn_nrg_count);
  pivot += rposn_nrg_count;
  accumulate_rbond_energy.setPointer(&insr_uint_data, pivot, rbond_nrg_count);
  pivot += rbond_nrg_count;
  accumulate_rangl_energy.setPointer(&insr_uint_data, pivot, rangl_nrg_count);
  pivot += rangl_nrg_count;
  accumulate_rdihe_energy.setPointer(&insr_uint_data, pivot, rdihe_nrg_count);  
  pivot = 0LLU;
  vwu_manipulation_masks.setPointer(&insr_uint2_data, pivot, manipulation_mask_count);
  pivot += manipulation_mask_count;
  cbnd_instructions.setPointer(&insr_uint2_data, pivot, cbnd_count);
  pivot += cbnd_count;
  angl_instructions.setPointer(&insr_uint2_data, pivot, angl_count);
  pivot += angl_count;
  cdhe_instructions.setPointer(&insr_uint2_data, pivot, cdhe_count);
  pivot += cdhe_count;
  cmap_instructions.setPointer(&insr_uint2_data, pivot, cmap_count);
  pivot += cmap_count;
  rposn_instructions.setPointer(&insr_uint2_data, pivot, rposn_count);
  pivot += rposn_count;
  rbond_instructions.setPointer(&insr_uint2_data, pivot, rbond_count);
  pivot += rbond_count;
  rangl_instructions.setPointer(&insr_uint2_data, pivot, rangl_count);
  pivot += rangl_count;
  rdihe_instructions.setPointer(&insr_uint2_data, pivot, rdihe_count);
  pivot += rdihe_count;
  vste_instructions.setPointer(&insr_uint2_data, pivot, vsite_count);
  pivot += vsite_count;
  sett_instructions.setPointer(&insr_uint2_data, pivot, settle_count);
  pivot += settle_count;
  cnst_instructions.setPointer(&insr_uint2_data, pivot, cgroup_count);
  pivot += cgroup_count;

  // Populate the instruction arrays with appropriate mapping and imported atom offsets
  int atom_import_idx = 0;
  int manipulation_mask_idx = 0;
  int cbnd_idx = 0;
  int angl_idx = 0;
  int cdhe_idx = 0;
  int cmap_idx = 0;
  int infr14_idx = 0;
  int rposn_idx = 0;
  int rbond_idx = 0;
  int rangl_idx = 0;
  int rdihe_idx = 0;
  int vsite_idx = 0;
  int settle_idx = 0;
  int cgroup_idx = 0;
  int cbnd_nrg_idx = 0;
  int angl_nrg_idx = 0;
  int cdhe_nrg_idx = 0;
  int cmap_nrg_idx = 0;
  int infr14_nrg_idx = 0;
  int rposn_nrg_idx = 0;
  int rbond_nrg_idx = 0;
  int rangl_nrg_idx = 0;
  int rdihe_nrg_idx = 0;
  int vwu_idx = 0;
  const std::vector<int> synth_cnst_param_bounds = constraint_param_bounds.readHost();
  for (int i = 0; i < system_count; i++) {
    const size_t nvwu = all_vwu[i].size();
    const int tp_idx = topology_indices.readHost(i);
    const int rn_idx = restraint_indices.readHost(i);
    const int2 bond_maplim  = bond_param_map_bounds.readHost(tp_idx);
    const int2 angl_maplim  = angl_param_map_bounds.readHost(tp_idx);
    const int2 dihe_maplim  = dihe_param_map_bounds.readHost(tp_idx);
    const int2 ubrd_maplim  = ubrd_param_map_bounds.readHost(tp_idx);
    const int2 cimp_maplim  = cimp_param_map_bounds.readHost(tp_idx);
    const int2 cmap_maplim  = cmap_param_map_bounds.readHost(tp_idx);
    const int2 attn_maplim  = attn14_param_map_bounds.readHost(tp_idx);
    const int2 rposn_maplim = rposn_param_map_bounds.readHost(rn_idx);
    const int2 rbond_maplim = rbond_param_map_bounds.readHost(rn_idx);
    const int2 rangl_maplim = rangl_param_map_bounds.readHost(rn_idx);
    const int2 rdihe_maplim = rdihe_param_map_bounds.readHost(rn_idx);
    const int2 vste_maplim  = vste_param_map_bounds.readHost(tp_idx);
    const int2 sett_maplim  = sett_param_map_bounds.readHost(tp_idx);
    const int2 cnst_maplim  = cnst_param_map_bounds.readHost(tp_idx);
    const std::vector<int> sysi_bond_map = bond_param_map.readHost(bond_maplim.x,
                                                                   bond_maplim.y - bond_maplim.x);
    const std::vector<int> sysi_angl_map = angl_param_map.readHost(angl_maplim.x,
                                                                   angl_maplim.y - angl_maplim.x);
    const std::vector<int> sysi_dihe_map = dihe_param_map.readHost(dihe_maplim.x,
                                                                   dihe_maplim.y - dihe_maplim.x);
    const std::vector<int> sysi_ubrd_map = ubrd_param_map.readHost(ubrd_maplim.x,
                                                                   ubrd_maplim.y - ubrd_maplim.x);
    const std::vector<int> sysi_cimp_map = cimp_param_map.readHost(cimp_maplim.x,
                                                                   cimp_maplim.y - cimp_maplim.x);
    const std::vector<int> sysi_cmap_map = cmap_param_map.readHost(cmap_maplim.x,
                                                                   cmap_maplim.y - cmap_maplim.x);
    const std::vector<int> sysi_attn_map =
      attn14_param_map.readHost(attn_maplim.x, attn_maplim.y - attn_maplim.x);
    const std::vector<int> sysi_rposn_kr_map =
      rposn_kr_param_map.readHost(rposn_maplim.x, rposn_maplim.y - rposn_maplim.x);
    const std::vector<int> sysi_rposn_xyz_map =
      rposn_xyz_param_map.readHost(rposn_maplim.x, rposn_maplim.y - rposn_maplim.x);
    const std::vector<int> sysi_rbond_map =
      rbond_param_map.readHost(rbond_maplim.x, rbond_maplim.y - rbond_maplim.x);
    const std::vector<int> sysi_rangl_map =
      rangl_param_map.readHost(rangl_maplim.x, rangl_maplim.y - rangl_maplim.x);
    const std::vector<int> sysi_rdihe_map =
      rdihe_param_map.readHost(rdihe_maplim.x, rdihe_maplim.y - rdihe_maplim.x);
    const std::vector<int> sysi_vste_map = vste_param_map.readHost(vste_maplim.x,
                                                                   vste_maplim.y - vste_maplim.x);
    const std::vector<int> sysi_sett_map = sett_param_map.readHost(sett_maplim.x,
                                                                   sett_maplim.y - sett_maplim.x);
    const std::vector<int> sysi_cnst_map = cnst_param_map.readHost(cnst_maplim.x,
                                                                   cnst_maplim.y - cnst_maplim.x);
    for (size_t j = 0LLU; j < nvwu; j++) {

      // Reconfigure the work units to store instructions appropriate for the AtomGraphSynthesis.
      all_vwu[i][j].storeCompositeBondInstructions(sysi_bond_map, sysi_ubrd_map);
      all_vwu[i][j].storeAngleInstructions(sysi_angl_map);
      all_vwu[i][j].storeCompositeDihedralInstructions(sysi_dihe_map, sysi_attn_map,
                                                       sysi_cimp_map);
      all_vwu[i][j].storeCmapInstructions(sysi_cmap_map);
      all_vwu[i][j].storeInferred14Instructions(sysi_attn_map);
      all_vwu[i][j].storePositionalRestraintInstructions(sysi_rposn_kr_map, sysi_rposn_xyz_map);
      all_vwu[i][j].storeDistanceRestraintInstructions(sysi_rbond_map);
      all_vwu[i][j].storeAngleRestraintInstructions(sysi_rangl_map);
      all_vwu[i][j].storeDihedralRestraintInstructions(sysi_rdihe_map);
      all_vwu[i][j].storeVirtualSiteInstructions(sysi_vste_map);
      all_vwu[i][j].storeSettleGroupInstructions(sysi_sett_map);
      all_vwu[i][j].storeConstraintGroupInstructions(sysi_cnst_map, synth_cnst_param_bounds);

      // Get limits for all data components of the valence work unit
      const int2* vwu_abs_ptr      = &vwu_instruction_sets.data()[vwu_idx * vwu_abstract_length];
      const int2 import_limits     = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::IMPORT)];
      const int2 manipulate_limits = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::MANIPULATE)];
      const int2 cbnd_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CBND)];
      const int2 angl_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::ANGL)];
      const int2 cdhe_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CDHE)];
      const int2 cmap_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CMAP)];
      const int2 infr14_limits     = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::INFR14)];
      const int2 rposn_limits      = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RPOSN)];
      const int2 rbond_limits      = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RBOND)];
      const int2 rangl_limits      = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RANGL)];
      const int2 rdihe_limits      = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RDIHE)];
      const int2 vste_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::VSITE)];
      const int2 sett_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::SETTLE)];
      const int2 cnst_limits       = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CGROUP)];
      const int2 cbnd_nrg_limits   = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CBND_NRG)];
      const int2 angl_nrg_limits   = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::ANGL_NRG)];
      const int2 cdhe_nrg_limits   = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CDHE_NRG)];
      const int2 cmap_nrg_limits   = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::CMAP_NRG)];
      const int2 infr14_nrg_limits = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::INFR14_NRG)];
      const int2 rposn_nrg_limits  = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RPOSN_NRG)];
      const int2 rbond_nrg_limits  = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RBOND_NRG)];
      const int2 rangl_nrg_limits  = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RANGL_NRG)];
      const int2 rdihe_nrg_limits  = vwu_abs_ptr[static_cast<int>(VwuAbstractMap::RDIHE_NRG)];

      // Upload the atom import array
      const int atom_start = atom_offsets.readHost(i);
      int* int_ptr = vwu_import_lists.data();
      for (int k = import_limits.x; k < import_limits.y; k++) {
        int_ptr[k] = all_vwu[i][j].getImportedAtomIndex(k - import_limits.x, atom_start);
      }

      // Upload the manipulation and various energy accumulation bitmasks
      vwu_manipulation_masks.putHost(all_vwu[i][j].getAtomManipulationMasks(), manipulate_limits.x,
                                     manipulate_limits.y - manipulate_limits.x);
      accumulate_cbnd_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::CBND),
                                     cbnd_nrg_limits.x, cbnd_nrg_limits.y - cbnd_nrg_limits.x);
      accumulate_angl_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::ANGL),
                                     angl_nrg_limits.x, angl_nrg_limits.y - angl_nrg_limits.x);
      accumulate_cdhe_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::CDHE),
                                     cdhe_nrg_limits.x, cdhe_nrg_limits.y - cdhe_nrg_limits.x);
      accumulate_cmap_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::CMAP),
                                     cmap_nrg_limits.x, cmap_nrg_limits.y - cmap_nrg_limits.x);
      accumulate_infr14_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::INFR14),
                                       infr14_nrg_limits.x,
                                       infr14_nrg_limits.y - infr14_nrg_limits.x);
      accumulate_rposn_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::RPOSN),
                                      rposn_nrg_limits.x, rposn_nrg_limits.y - rposn_nrg_limits.x);
      accumulate_rbond_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::RBOND),
                                      rbond_nrg_limits.x, rbond_nrg_limits.y - rbond_nrg_limits.x);
      accumulate_rangl_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::RANGL),
                                      rangl_nrg_limits.x, rangl_nrg_limits.y - rangl_nrg_limits.x);
      accumulate_rdihe_energy.putHost(all_vwu[i][j].getAccumulationFlags(VwuTask::RDIHE),
                                      rdihe_nrg_limits.x, rdihe_nrg_limits.y - rdihe_nrg_limits.x);

      // Upload the term instructions
      cbnd_instructions.putHost(all_vwu[i][j].getCompositeBondInstructions(), cbnd_limits.x,
                                cbnd_limits.y - cbnd_limits.x);
      angl_instructions.putHost(all_vwu[i][j].getAngleInstructions(), angl_limits.x,
                                angl_limits.y - angl_limits.x);
      const std::vector<uint3> cdhe_insr = all_vwu[i][j].getCompositeDihedralInstructions();
      uint2* uint2_ptr = cdhe_instructions.data();
      uint* uint_ptr = cdhe_overtones.data();
      for (int k = cdhe_limits.x; k < cdhe_limits.y; k++) {
        uint2_ptr[k].x = cdhe_insr[k - cdhe_limits.x].x;
        uint2_ptr[k].y = cdhe_insr[k - cdhe_limits.x].y;
        uint_ptr[k]    = cdhe_insr[k - cdhe_limits.x].z;
      }
      cmap_instructions.putHost(all_vwu[i][j].getCmapInstructions(), cmap_limits.x,
                                cmap_limits.y - cmap_limits.x);
      infr14_instructions.putHost(all_vwu[i][j].getInferred14Instructions(), infr14_limits.x,
                                  infr14_limits.y - infr14_limits.x);
      rposn_instructions.putHost(all_vwu[i][j].getPositionalRestraintInstructions(),
                                 rposn_limits.x, rposn_limits.y - rposn_limits.x);
      rbond_instructions.putHost(all_vwu[i][j].getDistanceRestraintInstructions(), rbond_limits.x,
                                 rbond_limits.y - rbond_limits.x);
      rangl_instructions.putHost(all_vwu[i][j].getAngleRestraintInstructions(), rangl_limits.x,
                                 rangl_limits.y - rangl_limits.x);
      rdihe_instructions.putHost(all_vwu[i][j].getDihedralRestraintInstructions(), rdihe_limits.x,
                                 rdihe_limits.y - rdihe_limits.x);
      vste_instructions.putHost(all_vwu[i][j].getVirtualSiteInstructions(), vste_limits.x,
                                vste_limits.y - vste_limits.x);
      sett_instructions.putHost(all_vwu[i][j].getSettleGroupInstructions(), sett_limits.x,
                                sett_limits.y - sett_limits.x);
      cnst_instructions.putHost(all_vwu[i][j].getConstraintGroupInstructions(), cnst_limits.x,
                                cnst_limits.y - cnst_limits.x);
      vwu_idx++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::loadNonbondedWorkUnits(const StaticExclusionMaskSynthesis &poly_se,
                                                const InitializationTask init_request,
                                                const int random_cache_depth,
                                                const GpuDetails &gpu) {
  
  // Build a list of non-bonded work units from the exclusion mask.  Translate those work units
  // into non-bonded abstracts and instructions.  All work units are assumed to have the same
  // target tile count.
  const std::vector<NonbondedWorkUnit> nbwu_list = buildNonbondedWorkUnits(poly_se, init_request,
                                                                           random_cache_depth,
                                                                           gpu);
  const int nbwu_count = nbwu_list.size();
  total_nonbonded_work_units = nbwu_count;
  int max_tile_count = 0;
  size_t total_tiles = 0LLU;
  size_t padded_tile_count = 0LLU;
  for (int i = 0; i < nbwu_count; i++) {
    const int n_tiles = nbwu_list[i].getTileCount();
    max_tile_count = std::max(n_tiles, max_tile_count);
    total_tiles += n_tiles;
    padded_tile_count += roundUp(n_tiles, warp_size_int);
  }
  nonbonded_work_type = (max_tile_count <= large_nbwu_tiles) ? NbwuKind::TILE_GROUPS :
                                                               NbwuKind::SUPERTILES;
  switch (nonbonded_work_type) {
  case NbwuKind::TILE_GROUPS:
    nonbonded_abstracts.resize(nbwu_count * tile_groups_wu_abstract_length);
    nbwu_instructions.resize(padded_tile_count);
    break;
  case NbwuKind::SUPERTILES:
    nonbonded_abstracts.resize(nbwu_count * supertile_wu_abstract_length);
    break;
  case NbwuKind::HONEYCOMB:
  case NbwuKind::UNKNOWN:
    break;
  }
  int insr_offset = 0;
  int abstract_offset = 0;
  for (int i = 0; i < nbwu_count; i++) {
    switch (nonbonded_work_type) {
    case NbwuKind::TILE_GROUPS:
      {
        const int n_tiles = nbwu_list[i].getTileCount();
        nonbonded_abstracts.putHost(nbwu_list[i].getAbstract(insr_offset),
                                    tile_groups_wu_abstract_length * i,
                                    tile_groups_wu_abstract_length);
        nbwu_instructions.putHost(nbwu_list[i].getTileInstructions(), insr_offset, n_tiles);
        insr_offset += roundUp(n_tiles, warp_size_int);
      }
      break;
    case NbwuKind::SUPERTILES:
      {
        nonbonded_abstracts.putHost(nbwu_list[i].getAbstract(insr_offset),
                                    supertile_wu_abstract_length * i,
                                    supertile_wu_abstract_length);
      }
      break;
    case NbwuKind::HONEYCOMB:
    case NbwuKind::UNKNOWN:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::loadReductionWorkUnits(const GpuDetails &gpu) {
  const std::vector<int> vatom_offsets = atom_offsets.readHost();
  const std::vector<int> vatom_counts  = atom_counts.readHost();
  const std::vector<ReductionWorkUnit> rdwu_vec = buildReductionWorkUnits(vatom_offsets,
                                                                          vatom_counts, gpu, 3);
  total_reduction_work_units = rdwu_vec.size();
  bool multi_wu = false;
  reduction_abstracts.resize(total_reduction_work_units * rdwu_abstract_length);
  const int atom_start_pos = static_cast<int>(RdwuAbstractMap::ATOM_START);
  const int atom_end_pos   = static_cast<int>(RdwuAbstractMap::ATOM_END);
  const int result_pos     = static_cast<int>(RdwuAbstractMap::RESULT_INDEX);
  const int depn_start_pos = static_cast<int>(RdwuAbstractMap::DEPN_START);
  const int depn_end_pos   = static_cast<int>(RdwuAbstractMap::DEPN_END);
  const int system_pos     = static_cast<int>(RdwuAbstractMap::SYSTEM_ID);
  int* abstract_ptr = reduction_abstracts.data();
  for (int i = 0; i < total_reduction_work_units; i++) {
    const int wu_pos = (rdwu_abstract_length * i);
    abstract_ptr[wu_pos + atom_start_pos] = rdwu_vec[i].getAtomStart(); 
    abstract_ptr[wu_pos + atom_end_pos]   = rdwu_vec[i].getAtomEnd(); 
    abstract_ptr[wu_pos + result_pos]     = rdwu_vec[i].getResultIndex(); 
    abstract_ptr[wu_pos + depn_start_pos] = rdwu_vec[i].getDependencyStart(); 
    abstract_ptr[wu_pos + depn_end_pos]   = rdwu_vec[i].getDependencyEnd(); 
    abstract_ptr[wu_pos + system_pos]     = rdwu_vec[i].getSystemIndex(); 
    multi_wu = (multi_wu ||
                (abstract_ptr[wu_pos + depn_end_pos] - abstract_ptr[wu_pos + depn_start_pos] > 1));
  }
  rdwu_per_system = (multi_wu) ? RdwuPerSystem::MULTIPLE : RdwuPerSystem::ONE;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getUniqueTopologyCount() const {
  return topology_count;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* AtomGraphSynthesis::getSystemTopologyPointer(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("A synthesis with " + std::to_string(system_count) + " systems cannot produce a "
          "topology relating to system index " + std::to_string(system_index) + ".",
          "AtomGraphSynthesis", "getSystemTopologyPointer");
  }
  return topologies[topology_indices.readHost(system_index)];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> AtomGraphSynthesis::getSystemTopologyPointer() const {
  std::vector<AtomGraph*> result(system_count);
  const int* top_idx_ptr = topology_indices.data();
  for (int i = 0; i < system_count; i++) {
    result[i] = topologies[top_idx_ptr[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*>& AtomGraphSynthesis::getUniqueTopologies() const {
  return topologies;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraphSynthesis::getTopologyIndices() const {
  return topology_indices.readHost();
}

//-------------------------------------------------------------------------------------------------
RestraintApparatus* AtomGraphSynthesis::getSystemRestraintPointer(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("A synthesis with " + std::to_string(system_count) + " systems cannot produce a "
          "restraint apparatus relating to system index " + std::to_string(system_index) + ".",
          "AtomGraphSynthesis", "getSystemRestraintPointer");
  }
  return restraint_networks[restraint_indices.readHost(system_index)];
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAtomCount() const {
  return total_atoms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAtomCount(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("A synthesis with " + std::to_string(system_count) + " systems cannot produce an atom "
          "count for system index " + std::to_string(system_index) + ".", "AtomGraphSynthesis",
          "getAtomCount");
  }
  return atom_counts.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getPaddedAtomCount() const {
  return atom_offsets.readHost(system_count - 1) +
         roundUp(atom_counts.readHost(system_count - 1), warp_size_int);
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAtomOffset(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("A synthesis with " + std::to_string(system_count) + " systems cannot produce an atom "
          "offset for system index " + std::to_string(system_index) + ".", "AtomGraphSynthesis",
          "getAtomOffset");
  }
  return atom_offsets.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getVirtualSiteCount() const {
  return total_virtual_sites;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getBondTermCount() const {
  return total_bond_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAngleTermCount() const {
  return total_angl_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getDihedralTermCount() const {
  return total_dihe_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getUreyBradleyTermCount() const {
  return total_ubrd_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCharmmImproperTermCount() const {
  return total_cimp_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCmapTermCount() const {
  return total_cmap_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getLJTypeCount() const {
  return total_lj_types;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getChargeTypeCount() const {
  return total_charge_types;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getBondParameterCount() const {
  return total_bond_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAngleParameterCount() const {
  return total_angl_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getDihedralParameterCount() const {
  return total_dihe_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getUreyBradleyParameterCount() const {
  return total_ubrd_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCharmmImproperParameterCount() const {
  return total_cimp_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCmapSurfaceCount() const {
  return total_cmap_surfaces;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int>& AtomGraphSynthesis::getSystemAtomCounts() const {
  return atom_counts;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int>& AtomGraphSynthesis::getSystemAtomOffsets() const {
  return atom_offsets;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int>& AtomGraphSynthesis::getDegreesOfFreedom() const {
  return degrees_of_freedom;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int>& AtomGraphSynthesis::getConstrainedDegreesOfFreedom() const {
  return cnst_degrees_of_freedom;
}

//-------------------------------------------------------------------------------------------------
UnitCellType AtomGraphSynthesis::getUnitCellType() const {
  return periodic_box_class;
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventModel AtomGraphSynthesis::getImplicitSolventModel() const {
  return gb_style;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphSynthesis::getDielectricConstant() const {
  return dielectric_constant;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphSynthesis::getSaltConcentration() const {
  return salt_concentration;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphSynthesis::getCoulombConstant() const {
  return coulomb_constant;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomicRadiusSet> AtomGraphSynthesis::getPBRadiiSet() const {
  return pb_radii_sets;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomicRadiusSet> AtomGraphSynthesis::getPBRadiiSet(const int low_limit,
                                                           const int high_limit) const {
  if (high_limit < low_limit || high_limit < 0 || low_limit < 0 || high_limit > system_count ||
      low_limit > system_count) {
    rtErr("The requested range " + std::to_string(low_limit) + " - " + std::to_string(high_limit) +
          " is invalid for a collection of " + std::to_string(system_count) + " systems.",
          "AtomGraphSynthesis", "getPBRadiiSet");
  }
  std::vector<AtomicRadiusSet> result(high_limit - low_limit);
  for (int i = low_limit; i < high_limit; i++) {
    result[i - low_limit] = pb_radii_sets[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
AtomicRadiusSet AtomGraphSynthesis::getPBRadiiSet(const int index) const {
  if (index < 0 || index > system_count) {
    rtErr("The index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "AtomGraphSynthesis", "getPBRadiiSet");
  }
  return pb_radii_sets[index];
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getValenceWorkUnitCount() const {
  return total_valence_work_units;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getValenceWorkUnitSize() const {
  return valence_work_unit_size;
}

//-------------------------------------------------------------------------------------------------
ValenceKernelSize AtomGraphSynthesis::getValenceThreadBlockSize() const {
  return valence_thread_block_size;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int2>& AtomGraphSynthesis::getValenceWorkUnitAbstracts() const {
  return vwu_instruction_sets;
}

//-------------------------------------------------------------------------------------------------
NbwuKind AtomGraphSynthesis::getNonbondedWorkType() const {
  return nonbonded_work_type;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getNonbondedWorkUnitCount() const {
  return total_nonbonded_work_units;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getRandomCacheDepth() const {
  if (total_nonbonded_work_units == 0) {
    rtWarn("No non-bonded work units have yet been constructed for this topology synthesis.  A "
           "cache depth of zero will be assumed until the proper value can be read from the first "
           "work unit.", "AtomGraphSynthesis", "getRandomCacheDepth");
    return 0;
  }
  else {
    const int refresh_code = nonbonded_abstracts.readHost((2 * small_block_max_imports) + 10);
    return ((refresh_code >> 8) & 0xff);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getReductionWorkUnitCount() const {
  return total_reduction_work_units;
}

//-------------------------------------------------------------------------------------------------
RdwuPerSystem AtomGraphSynthesis::getRdwuPerSystem() const {
  return rdwu_per_system;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<int>& AtomGraphSynthesis::getReductionWorkUnitAbstracts() const {
  return reduction_abstracts;  
}

//-------------------------------------------------------------------------------------------------
SyValenceKit<double>
AtomGraphSynthesis::getDoublePrecisionValenceKit(const HybridTargetLevel tier) const {
  return SyValenceKit<double>(total_valence_work_units, coulomb_constant,
                              bond_stiffnesses.data(tier), bond_equilibria.data(tier),
                              angl_stiffnesses.data(tier), angl_equilibria.data(tier),
                              dihe_amplitudes.data(tier), dihe_periodicities.data(tier),
                              dihe_phase_angles.data(tier), attn14_elec_factors.data(tier),
                              attn14_vdw_factors.data(tier), atomic_charges.data(tier),
                              lennard_jones_14_a_coeff.data(tier),
                              lennard_jones_14_b_coeff.data(tier),
                              lennard_jones_14_c_coeff.data(tier),
                              lennard_jones_14_sigma.data(tier), lennard_jones_indices.data(tier),
                              lj_type_counts.data(tier), lennard_jones_abc_offsets.data(tier),
                              ubrd_stiffnesses.data(tier), ubrd_equilibria.data(tier),
                              cimp_stiffnesses.data(tier), cimp_phase_angles.data(tier),
                              cmap_surface_dimensions.data(tier), cmap_patches.data(tier),
                              cmap_patch_bounds.data(tier), vwu_instruction_sets.data(tier),
                              vwu_import_lists.data(tier), cbnd_instructions.data(tier),
                              angl_instructions.data(tier), cdhe_instructions.data(tier),
                              cdhe_overtones.data(tier), cmap_instructions.data(tier),
                              infr14_instructions.data(tier), accumulate_cbnd_energy.data(tier),
                              accumulate_angl_energy.data(tier), accumulate_cdhe_energy.data(tier),
                              accumulate_cmap_energy.data(tier),
                              accumulate_infr14_energy.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyValenceKit<float>
AtomGraphSynthesis::getSinglePrecisionValenceKit(const HybridTargetLevel tier) const {
  return SyValenceKit<float>(total_valence_work_units, coulomb_constant,
                             sp_bond_stiffnesses.data(tier), sp_bond_equilibria.data(tier),
                             sp_angl_stiffnesses.data(tier), sp_angl_equilibria.data(tier),
                             sp_dihe_amplitudes.data(tier), sp_dihe_periodicities.data(tier),
                             sp_dihe_phase_angles.data(tier), sp_attn14_elec_factors.data(tier),
                             sp_attn14_vdw_factors.data(tier), sp_atomic_charges.data(tier),
                             sp_lennard_jones_14_a_coeff.data(tier),
                             sp_lennard_jones_14_b_coeff.data(tier),
                             sp_lennard_jones_14_c_coeff.data(tier),
                             sp_lennard_jones_14_sigma.data(tier),
                             lennard_jones_indices.data(tier), lj_type_counts.data(tier),
                             lennard_jones_abc_offsets.data(tier), sp_ubrd_stiffnesses.data(tier),
                             sp_ubrd_equilibria.data(tier), sp_cimp_stiffnesses.data(tier),
                             sp_cimp_phase_angles.data(tier), cmap_surface_dimensions.data(tier),
                             sp_cmap_patches.data(tier), cmap_patch_bounds.data(tier),
                             vwu_instruction_sets.data(tier), vwu_import_lists.data(tier),
                             cbnd_instructions.data(tier), angl_instructions.data(tier),
                             cdhe_instructions.data(tier), cdhe_overtones.data(tier),
                             cmap_instructions.data(tier), infr14_instructions.data(tier),
                             accumulate_cbnd_energy.data(tier), accumulate_angl_energy.data(tier),
                             accumulate_cdhe_energy.data(tier), accumulate_cmap_energy.data(tier),
                             accumulate_infr14_energy.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyRestraintKit<double, double2, double4>
AtomGraphSynthesis::getDoublePrecisionRestraintKit(const HybridTargetLevel tier) const {
  return SyRestraintKit<double,
                        double2,
                        double4>(rposn_step_bounds.data(tier), rbond_step_bounds.data(tier),
                                 rangl_step_bounds.data(tier), rdihe_step_bounds.data(tier),
                                 rposn_init_k.data(tier), rposn_final_k.data(tier),
                                 rposn_init_r.data(tier), rposn_final_r.data(tier),
                                 rposn_init_xy.data(tier), rposn_init_z.data(tier),
                                 rposn_final_xy.data(tier), rposn_final_z.data(tier),
                                 rbond_init_k.data(tier), rbond_final_k.data(tier),
                                 rbond_init_r.data(tier), rbond_final_r.data(tier),
                                 rangl_init_k.data(tier), rangl_final_k.data(tier),
                                 rangl_init_r.data(tier), rangl_final_r.data(tier),
                                 rdihe_init_k.data(tier), rdihe_final_k.data(tier),
                                 rdihe_init_r.data(tier), rdihe_final_r.data(tier),
                                 rposn_instructions.data(tier), rbond_instructions.data(tier),
                                 rangl_instructions.data(tier), rdihe_instructions.data(tier),
                                 accumulate_rposn_energy.data(tier),
                                 accumulate_rbond_energy.data(tier),
                                 accumulate_rangl_energy.data(tier),
                                 accumulate_rdihe_energy.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyRestraintKit<float, float2, float4>
AtomGraphSynthesis::getSinglePrecisionRestraintKit(const HybridTargetLevel tier) const {
  return SyRestraintKit<float,
                        float2,
                        float4>(rposn_step_bounds.data(tier), rbond_step_bounds.data(tier),
                                rangl_step_bounds.data(tier), rdihe_step_bounds.data(tier),
                                sp_rposn_init_k.data(tier), sp_rposn_final_k.data(tier),
                                sp_rposn_init_r.data(tier), sp_rposn_final_r.data(tier),
                                sp_rposn_init_xy.data(tier), sp_rposn_init_z.data(tier),
                                sp_rposn_final_xy.data(tier), sp_rposn_final_z.data(tier),
                                sp_rbond_init_k.data(tier), sp_rbond_final_k.data(tier),
                                sp_rbond_init_r.data(tier), sp_rbond_final_r.data(tier),
                                sp_rangl_init_k.data(tier), sp_rangl_final_k.data(tier),
                                sp_rangl_init_r.data(tier), sp_rangl_final_r.data(tier),
                                sp_rdihe_init_k.data(tier), sp_rdihe_final_k.data(tier),
                                sp_rdihe_init_r.data(tier), sp_rdihe_final_r.data(tier),
                                rposn_instructions.data(tier), rbond_instructions.data(tier),
                                rangl_instructions.data(tier), rdihe_instructions.data(tier),
                                accumulate_rposn_energy.data(tier),
                                accumulate_rbond_energy.data(tier),
                                accumulate_rangl_energy.data(tier),
                                accumulate_rdihe_energy.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyNonbondedKit<double, double2>
AtomGraphSynthesis::getDoublePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  return SyNonbondedKit<double,
                        double2>(system_count, periodic_box_class, total_nonbonded_work_units,
                                 nonbonded_abstracts.data(tier), nbwu_instructions.data(tier),
                                 atom_offsets.data(tier), atom_counts.data(tier), coulomb_constant,
                                 gb_style, neck_table_size, dielectric_constant, is_kappa,
                                 salt_concentration, gb_offset, gb_neckscale, gb_neckcut,
                                 atomic_charges.data(tier), charge_indices.data(tier),
                                 charge_parameters.data(tier), lennard_jones_indices.data(tier),
                                 lj_type_counts.data(tier), lennard_jones_abc_offsets.data(tier),
                                 lennard_jones_ab_coeff.data(tier),
                                 lennard_jones_c_coeff.data(tier), lennard_jones_sigma.data(tier),
                                 neck_gb_indices.data(tier), atomic_pb_radii.data(tier),
                                 gb_screening_factors.data(tier), gb_alpha_parameters.data(tier),
                                 gb_beta_parameters.data(tier), gb_gamma_parameters.data(tier),
                                 neck_limit_tables.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyNonbondedKit<float, float2>
AtomGraphSynthesis::getSinglePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  return SyNonbondedKit<float,
                        float2>(system_count, periodic_box_class, total_nonbonded_work_units,
                                nonbonded_abstracts.data(tier), nbwu_instructions.data(tier),
                                atom_offsets.data(tier), atom_counts.data(tier), coulomb_constant,
                                gb_style, neck_table_size, dielectric_constant, is_kappa,
                                salt_concentration, gb_offset, gb_neckscale, gb_neckcut,
                                sp_atomic_charges.data(tier), charge_indices.data(tier),
                                sp_charge_parameters.data(tier), lennard_jones_indices.data(tier),
                                lj_type_counts.data(tier), lennard_jones_abc_offsets.data(tier),
                                sp_lennard_jones_ab_coeff.data(tier),
                                sp_lennard_jones_c_coeff.data(tier),
                                sp_lennard_jones_sigma.data(tier), neck_gb_indices.data(tier),
                                sp_atomic_pb_radii.data(tier), sp_gb_screening_factors.data(tier),
                                sp_gb_alpha_parameters.data(tier),
                                sp_gb_beta_parameters.data(tier),
                                sp_gb_gamma_parameters.data(tier),
                                sp_neck_limit_tables.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyAtomUpdateKit<double, double2, double4>
AtomGraphSynthesis::getDoublePrecisionAtomUpdateKit(const HybridTargetLevel tier) const {
  return SyAtomUpdateKit<double,
                         double2,
                         double4>(atomic_masses.data(tier), inverse_atomic_masses.data(tier),
                                  largest_constraint_group, virtual_site_parameters.data(tier),
                                  settle_group_geometry.data(tier), settle_group_masses.data(tier),
                                  constraint_group_params.data(tier), vste_instructions.data(tier),
                                  sett_instructions.data(tier), cnst_instructions.data(tier),
                                  vwu_manipulation_masks.data(tier));
}

//-------------------------------------------------------------------------------------------------
SyAtomUpdateKit<float, float2, float4>
AtomGraphSynthesis::getSinglePrecisionAtomUpdateKit(const HybridTargetLevel tier) const {
  return SyAtomUpdateKit<float,
                         float2,
                         float4>(sp_atomic_masses.data(tier), sp_inverse_atomic_masses.data(tier),
                                 largest_constraint_group, sp_virtual_site_parameters.data(tier),
                                 sp_settle_group_geometry.data(tier),
                                 sp_settle_group_masses.data(tier),
                                 sp_constraint_group_params.data(tier),
                                 vste_instructions.data(tier), sett_instructions.data(tier),
                                 cnst_instructions.data(tier), vwu_manipulation_masks.data(tier));
}

//-------------------------------------------------------------------------------------------------
const AtomGraphSynthesis* AtomGraphSynthesis::getSelfPointer() const {
  return this;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::upload() {
  int_system_data.upload();
  chem_int_data.upload();
  chem_int2_data.upload();
  chem_double_data.upload();
  chem_float_data.upload();
  chem_char4_data.upload();
  valparam_double_data.upload();
  valparam_float_data.upload();
  valparam_int_data.upload();
  valparam_int2_data.upload();
  valence_int_data.upload();
  charge_parameters.upload();
  lennard_jones_ab_coeff.upload();
  lennard_jones_c_coeff.upload();
  lennard_jones_14_a_coeff.upload();
  lennard_jones_14_b_coeff.upload();
  lennard_jones_14_c_coeff.upload();
  lennard_jones_sigma.upload();
  lennard_jones_14_sigma.upload();
  sp_charge_parameters.upload();
  sp_lennard_jones_ab_coeff.upload();
  sp_lennard_jones_c_coeff.upload();
  sp_lennard_jones_14_a_coeff.upload();
  sp_lennard_jones_14_b_coeff.upload();
  sp_lennard_jones_14_c_coeff.upload();
  sp_lennard_jones_sigma.upload();
  sp_lennard_jones_14_sigma.upload();
  neck_limit_tables.upload();
  sp_neck_limit_tables.upload();
  nmr_int2_data.upload();
  nmr_double_data.upload();
  nmr_double2_data.upload();
  nmr_double4_data.upload();
  nmr_float_data.upload();
  nmr_float2_data.upload();
  nmr_float4_data.upload();
  nmr_int_data.upload();
  nmr_int2_data.upload();
  virtual_site_parameters.upload();
  sp_virtual_site_parameters.upload();
  vsite_int_data.upload();
  settle_group_geometry.upload();
  settle_group_masses.upload();
  sp_settle_group_geometry.upload();
  sp_settle_group_masses.upload();
  settle_group_indexing.upload();
  constraint_group_indices.upload();
  constraint_group_bounds.upload();
  constraint_group_param_idx.upload();
  constraint_param_bounds.upload();
  constraint_group_params.upload();
  sp_constraint_group_params.upload();
  vwu_import_lists.upload();
  vwu_instruction_sets.upload();
  insr_uint_data.upload();
  insr_uint2_data.upload();
  nonbonded_abstracts.upload();
  nbwu_instructions.upload();
  reduction_abstracts.upload();
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::download() {
  int_system_data.download();
  chem_int_data.download();
  chem_int2_data.download();
  chem_double_data.download();
  chem_float_data.download();
  chem_char4_data.download();
  valparam_double_data.download();
  valparam_float_data.download();
  valparam_int_data.download();
  valparam_int2_data.download();
  valence_int_data.download();
  charge_parameters.download();
  lennard_jones_ab_coeff.download();
  lennard_jones_c_coeff.download();
  lennard_jones_14_a_coeff.download();
  lennard_jones_14_b_coeff.download();
  lennard_jones_14_c_coeff.download();
  lennard_jones_sigma.download();
  lennard_jones_14_sigma.download();
  sp_charge_parameters.download();
  sp_lennard_jones_ab_coeff.download();
  sp_lennard_jones_c_coeff.download();
  sp_lennard_jones_14_a_coeff.download();
  sp_lennard_jones_14_b_coeff.download();
  sp_lennard_jones_14_c_coeff.download();
  sp_lennard_jones_sigma.download();
  sp_lennard_jones_14_sigma.download();
  neck_limit_tables.download();
  sp_neck_limit_tables.download();
  nmr_int2_data.download();
  nmr_double_data.download();
  nmr_double2_data.download();
  nmr_double4_data.download();
  nmr_float_data.download();
  nmr_float2_data.download();
  nmr_float4_data.download();
  nmr_int_data.download();
  nmr_int2_data.download();
  virtual_site_parameters.download();
  sp_virtual_site_parameters.download();
  vsite_int_data.download();
  settle_group_geometry.download();
  settle_group_masses.download();
  sp_settle_group_geometry.download();
  sp_settle_group_masses.download();
  settle_group_indexing.download();
  constraint_group_indices.download();
  constraint_group_bounds.download();
  constraint_group_param_idx.download();
  constraint_param_bounds.download();
  constraint_group_params.download();
  sp_constraint_group_params.download();
  vwu_import_lists.download();
  vwu_instruction_sets.download();
  insr_uint_data.download();
  insr_uint2_data.download();
  nonbonded_abstracts.download();
  nbwu_instructions.download();
  reduction_abstracts.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::importImplicitSolventAtomParameters(const int system_index) {
  const int top_idx = topology_indices.readHost(system_index);
  const AtomGraph *ag_ptr = topologies[top_idx];
  const ImplicitSolventKit<double> isk = ag_ptr->getDoublePrecisionImplicitSolventKit();
  const int imin = atom_offsets.readHost(system_index);
  const int imax = imin + atom_counts.readHost(system_index);
  int* neck_gb_indices_ptr = neck_gb_indices.data();
  double* atomic_pb_radii_ptr = atomic_pb_radii.data();
  double* gb_screen_ptr = gb_screening_factors.data();
  double* gb_alpha_ptr = gb_alpha_parameters.data();
  double* gb_beta_ptr = gb_beta_parameters.data();
  double* gb_gamma_ptr = gb_gamma_parameters.data();
  float* sp_atomic_pb_radii_ptr = sp_atomic_pb_radii.data();
  float* sp_gb_screen_ptr = sp_gb_screening_factors.data();
  float* sp_gb_alpha_ptr = sp_gb_alpha_parameters.data();
  float* sp_gb_beta_ptr = sp_gb_beta_parameters.data();
  float* sp_gb_gamma_ptr = sp_gb_gamma_parameters.data();
  for (int i = imin; i < imax; i++) {
    neck_gb_indices_ptr[i] = isk.neck_gb_idx[i - imin];
    atomic_pb_radii_ptr[i] = isk.pb_radii[i - imin];
    gb_screen_ptr[i] = isk.gb_screen[i - imin];
    gb_alpha_ptr[i] = isk.gb_alpha[i - imin];
    gb_beta_ptr[i] = isk.gb_beta[i - imin];
    gb_gamma_ptr[i] = isk.gb_gamma[i - imin];
    sp_atomic_pb_radii_ptr[i] = isk.pb_radii[i - imin];
    sp_gb_screen_ptr[i] = isk.gb_screen[i - imin];
    sp_gb_alpha_ptr[i] = isk.gb_alpha[i - imin];
    sp_gb_beta_ptr[i] = isk.gb_beta[i - imin];
    sp_gb_gamma_ptr[i] = isk.gb_gamma[i - imin];
  }  
}

//-------------------------------------------------------------------------------------------------
void
AtomGraphSynthesis::setImplicitSolventNeckParameters(const NeckGeneralizedBornTable &ngb_tab) {
  const NeckGeneralizedBornKit<double> ngbk = ngb_tab.dpData();
  neck_table_size = ngbk.table_size;
  neck_limit_tables.resize(neck_table_size * neck_table_size);
  sp_neck_limit_tables.resize(neck_table_size * neck_table_size);
  double2* neck_limit_ptr = neck_limit_tables.data();
  float2* sp_neck_limit_ptr = sp_neck_limit_tables.data();
  for (int i = 0; i < neck_table_size * neck_table_size; i++) {
    neck_limit_ptr[i].x = ngbk.max_separation[i];
    neck_limit_ptr[i].y = ngbk.max_value[i];
    sp_neck_limit_ptr[i].x = ngbk.max_separation[i];
    sp_neck_limit_ptr[i].y = ngbk.max_value[i];
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::setImplicitSolventNeckParameters() {
  const NeckGeneralizedBornTable ngb_tab;
  setImplicitSolventNeckParameters(ngb_tab);
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::setImplicitSolventModel() {

  // Set the global GB parameters based on the first topology with a valid implicit solvent model.
  // Check subsequent topologies to ensure that there is agreement where required.
  if (system_count == 0) {
    return;
  }
  const NeckGeneralizedBornTable ngb_tab;
  const NeckGeneralizedBornKit<double> ngbk = ngb_tab.dpData();
  int sys_idx = 0;
  int gb_setting_system;
  gb_style = ImplicitSolventModel::NONE;
  while (sys_idx < system_count && gb_style == ImplicitSolventModel::NONE) {
    const AtomGraph* ag_base = getSystemTopologyPointer(sys_idx);
    if (ag_base->getImplicitSolventModel() != ImplicitSolventModel::NONE) {
      gb_style = ag_base->getImplicitSolventModel();
      dielectric_constant = ag_base->getDielectricConstant();
      salt_concentration = ag_base->getSaltConcentration();
      is_kappa = (salt_concentration > constants::tiny) ?
                 sqrt(default_salt_kappa_dependence * salt_concentration) : 0.0;
      gb_offset = (gb_style == ImplicitSolventModel::NECK_GB_II) ?
                  static_cast<double>(default_neck_ii_gb_radii_offset) :
                  static_cast<double>(default_gb_radii_offset);
      gb_neckscale = (gb_style == ImplicitSolventModel::NECK_GB_II) ?
                     static_cast<double>(default_gb_neck_ii_scale) :
                     static_cast<double>(default_gb_neck_scale);
      gb_neckcut = ngbk.neck_cut;
      gb_setting_system = sys_idx;
    }
    sys_idx++;
  }

  // Check that all other topologies contain the same critical implicit solvent parameters
  while (sys_idx < system_count) {
    const AtomGraph* ag_base = getSystemTopologyPointer(sys_idx);
    const ImplicitSolventModel test_ism = ag_base->getImplicitSolventModel();
    switch (test_ism) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
      switch (gb_style) {
      case ImplicitSolventModel::NONE:
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        rtErr("A 'Neck' GB model of the form " + getEnumerationName(gb_style) + " set for the "
              "synthesis by system index " + std::to_string(gb_setting_system) + " (" +
              getBaseName(getSystemTopologyPointer(gb_setting_system)->getFileName()) + "), is "
              "not compatible with other GB models (" + getEnumerationName(test_ism) +
              " in system index " + std::to_string(sys_idx) + ", " +
              getBaseName(getSystemTopologyPointer(sys_idx)->getFileName()) + ").",
              "AtomGraphSynthesis", "setImplicitSolventModel");
        break;
      }
      break;
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (gb_style) {
      case ImplicitSolventModel::NONE:
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        rtErr("A GB model of the form " + getEnumerationName(gb_style) + " set for the "
              "synthesis by system index " + std::to_string(gb_setting_system) + " (" +
              getBaseName(getSystemTopologyPointer(gb_setting_system)->getFileName()) + "), is "
              "not compatible with a 'Neck' GB model (" + getEnumerationName(test_ism) +
              " in system index " + std::to_string(sys_idx) + ", " +
              getBaseName(getSystemTopologyPointer(sys_idx)->getFileName()) + ").",
              "AtomGraphSynthesis", "setImplicitSolventModel");
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        if (gb_style != test_ism) {
          rtErr("Multiple 'Neck' GB models cannot be placed in the same synthesis (system index " +
                std::to_string(gb_setting_system) + ", of style " + getEnumerationName(gb_style) +
                " and system index " + std::to_string(sys_idx) + " of style " +
                getEnumerationName(test_ism) + ", taken from " +
                getBaseName(getSystemTopologyPointer(gb_setting_system)->getFileName()) + " and " +
                getBaseName(getSystemTopologyPointer(sys_idx)->getFileName()) + ", respectively).",
                "AtomGraphSynthesis", "setImplicitSolventModel");
        }
        break;
      }
    }
    const double test_dielectric = ag_base->getDielectricConstant();
    if (fabs(test_dielectric - dielectric_constant) > 1.0e-4) {
      rtErr("Multiple dielectric constants (" +
            realToString(test_dielectric, 8, 4, NumberFormat::STANDARD_REAL) + " and " +
            realToString(dielectric_constant, 8, 4, NumberFormat::STANDARD_REAL) + " from " +
            getBaseName(getSystemTopologyPointer(sys_idx)->getFileName()) + " and " +
            getBaseName(getSystemTopologyPointer(gb_setting_system)->getFileName()) + " in system "
            "indices " + std::to_string(sys_idx) + " and " + std::to_string(gb_setting_system) +
            ", respectively) cannot be combined in the same synthesis.", "AtomGraphSynthesis",
            "setImplicitSolventModel");
    }
    const double test_salt = ag_base->getSaltConcentration();
    if (fabs(test_salt - salt_concentration) > constants::small) {
      rtErr("Multiple salt concentrations (" +
            realToString(test_salt, 8, 4, NumberFormat::STANDARD_REAL) + " and " +
            realToString(salt_concentration, 8, 4, NumberFormat::STANDARD_REAL) + " from " +
            getBaseName(getSystemTopologyPointer(sys_idx)->getFileName()) + " and " +
            getBaseName(getSystemTopologyPointer(gb_setting_system)->getFileName()) + " in system "
            "indices " + std::to_string(sys_idx) + " and " + std::to_string(gb_setting_system) +
            ", respectively) cannot be combined in the same synthesis.", "AtomGraphSynthesis",
            "setImplicitSolventModel");
    }
    sys_idx++;
  }
  
  // The original GB settings from the input topologies will be imported into the synthesis.
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    for (int i = 0; i < system_count; i++) {
      importImplicitSolventAtomParameters(i);
    }
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < system_count; i++) {
      importImplicitSolventAtomParameters(i);
    }
    setImplicitSolventNeckParameters();
    break;
  }
}
 
//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                                 const NeckGeneralizedBornKit<double> &ngbk,
                                                 const std::vector<AtomicRadiusSet> &radii_sets_in,
                                                 const double dielectric_in,
                                                 const double saltcon_in,
                                                 const ExceptionResponse policy) {

  // Check that there are the appropriate number of atomic radius sets
  if (static_cast<int>(radii_sets_in.size()) != system_count) {
    rtErr("For " + std::to_string(system_count) + " systems, " +
          std::to_string(radii_sets_in.size()) + " is an invalid number of unique atomic radii "
          "sets.", "AtomGraphSynthesis", "setImplicitSolventModel");
  }
  
  // Set features of the synthesis.  Here, the GB offset, neck scaling, and neck cutoff are
  // explicit members of the object, as opposed to the single-system AtomGraph where it is inferred
  // when an ImplicitSolventRecipe is constructed.  In either case, the values of the GB offset are
  // hard-coded.
  gb_style = igb_in;
  dielectric_constant = dielectric_in;
  salt_concentration = saltcon_in;
  is_kappa = (salt_concentration > constants::tiny) ?
             sqrt(default_salt_kappa_dependence * salt_concentration) : 0.0;
  gb_offset = (igb_in == ImplicitSolventModel::NECK_GB_II) ?
              static_cast<double>(default_neck_ii_gb_radii_offset) :
              static_cast<double>(default_gb_radii_offset);
  gb_neckscale = (igb_in == ImplicitSolventModel::NECK_GB_II) ?
                 static_cast<double>(default_gb_neck_ii_scale) :
                 static_cast<double>(default_gb_neck_scale);
  gb_neckcut = ngbk.neck_cut;

  // Set the implicit solvent models in each reference topology one by one, then port the results
  // to the appropriate array positions in the synthesis.
  std::vector<AtomicRadiusSet> top_rad_settings(topology_count);
  std::vector<bool> top_coverage(topology_count, false);
  for (int i = 0; i < system_count; i++) {
    const int top_idx = topology_indices.readHost(i);
    AtomGraph *iag_ptr = topologies[top_idx];
    if (top_coverage[top_idx]) {
      if (top_rad_settings[top_idx] != radii_sets_in[i]) {
        iag_ptr->setImplicitSolventModel(igb_in, dielectric_in, saltcon_in, radii_sets_in[i],
                                         policy);
        top_rad_settings[top_idx] = radii_sets_in[i];
      }
    }
    else {
      top_coverage[top_idx] = true;
      iag_ptr->setImplicitSolventModel(igb_in, dielectric_in, saltcon_in, radii_sets_in[i],
                                       policy);
      top_rad_settings[top_idx] = radii_sets_in[i];
    }
    importImplicitSolventAtomParameters(i);
  }

  // Load the neck Generalized Born tables, if applicable
  switch (igb_in) {
  case ImplicitSolventModel::NONE:
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    setImplicitSolventNeckParameters();
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                                 const NeckGeneralizedBornKit<double> &ngbk,
                                                 const AtomicRadiusSet radii_set_in,
                                                 const double dielectric_in,
                                                 const double saltcon_in,
                                                 const ExceptionResponse policy) {
  const std::vector rsets_in = std::vector<AtomicRadiusSet>(system_count, radii_set_in);
  setImplicitSolventModel(igb_in, ngbk, rsets_in, dielectric_in, saltcon_in, policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                                 const NeckGeneralizedBornTable &ngb_tab,
                                                 const std::vector<AtomicRadiusSet> &radii_sets_in,
                                                 const double dielectric_in,
                                                 const double saltcon_in,
                                                 const ExceptionResponse policy) {
  setImplicitSolventModel(igb_in, ngb_tab.dpData(), radii_sets_in, dielectric_in, saltcon_in,
                          policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                                 const NeckGeneralizedBornTable &ngb_tab,
                                                 const AtomicRadiusSet radii_set_in,
                                                 const double dielectric_in,
                                                 const double saltcon_in,
                                                 const ExceptionResponse policy) {
  const std::vector rsets_in = std::vector<AtomicRadiusSet>(system_count, radii_set_in);
  setImplicitSolventModel(igb_in, ngb_tab.dpData(), rsets_in, dielectric_in, saltcon_in, policy); 
}

} // namespace synthesis
} // namespace stormm
