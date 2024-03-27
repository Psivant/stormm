#include <algorithm>
#include <sys/time.h>
#include "copyright.h"
#include "Chemistry/periodic_table.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Math/series_ops.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "atomgraph_refinement.h"
#include "atomgraph_stage.h"

namespace stormm {
namespace topology {

using constants::ExceptionResponse;
using chemistry::elemental_masses;
using parse::char4ToString;
using parse::NumberFormat;
using parse::realToString;
using stmath::applyAssociatedSort;
using stmath::findBin;
using stmath::incrementingSeries;
using stmath::indexingArray;
using stmath::maxValue;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::roundUp;

//-------------------------------------------------------------------------------------------------
ScreeningFactorTag::ScreeningFactorTag() :
    dihe_idx{-1}, i_at{unknown_atom_type}, j_at{wildcard_atom_type}, k_at{wildcard_atom_type},
    l_at{unknown_atom_type}, elec14{1.0}, vdw14{1.0} 
{}

//-------------------------------------------------------------------------------------------------
ScreeningFactorTag::ScreeningFactorTag(const int dihe_idx_in, const double elec14_in,
                                       const double vdw14_in) :
    dihe_idx{dihe_idx_in}, i_at{unknown_atom_type}, j_at{wildcard_atom_type},
    k_at{wildcard_atom_type}, l_at{unknown_atom_type}, elec14{elec14_in}, vdw14{vdw14_in}
{}

//-------------------------------------------------------------------------------------------------
ScreeningFactorTag::ScreeningFactorTag(const char4 i_at_in, const char4 l_at_in,
                                       const double elec14_in, const double vdw14_in) :
    dihe_idx{-1}, i_at{i_at_in}, j_at{wildcard_atom_type}, k_at{wildcard_atom_type}, l_at{l_at_in},
    elec14{elec14_in}, vdw14{vdw14_in}
{}

//-------------------------------------------------------------------------------------------------
ScreeningFactorTag::ScreeningFactorTag(const char4 i_at_in, const char4 j_at_in,
                                       const char4 k_at_in, const char4 l_at_in,
                                       const double elec14_in, const double vdw14_in) :
    dihe_idx{-1}, i_at{i_at_in}, j_at{j_at_in}, k_at{k_at_in}, l_at{l_at_in}, elec14{elec14_in},
    vdw14{vdw14_in}
{}

//-------------------------------------------------------------------------------------------------
AtomGraphStage::AtomGraphStage(const int atom_count_in, const std::vector<int> &residue_limits_in,
                               const ExceptionResponse policy_in) :
    policy{policy_in},

    // Metadata
    title{}, output_file{}, force_fields{},

    // Sizing constants and arrays of limits
    atom_count{atom_count_in},
    residue_count{(residue_limits_in.size() > 0) ?
                  static_cast<int>(residue_limits_in.size()) - 1 : 0},
    molecule_count{0},
    residue_limits{residue_limits_in},
    molecule_limits{std::vector<int>(molecule_count + 1, 0)},
    molecule_membership{std::vector<int>(atom_count, 0)},
    molecule_contents{std::vector<int>(atom_count, 0)},

    // Atom and residue descriptors
    atomic_numbers{std::vector<int>(atom_count, 1)},
    atom_names{std::vector<char4>(atom_count, unknown_atom_name)},
    atom_types{std::vector<char4>(atom_count, unknown_atom_type)},
    residue_names{std::vector<char4>(residue_count, unknown_residue_name)},

    // CHARMM force field term and parameter set counts
    urey_bradley_term_count{0}, charmm_impr_term_count{0}, cmap_term_count{0},
    urey_bradley_parameter_count{0}, charmm_impr_parameter_count{0}, cmap_surface_count{0},

    // Atom and parameter indexing in CHARMM force field terms
    urey_bradley_i_atoms{}, urey_bradley_k_atoms{}, urey_bradley_parameter_indices{},
    charmm_impr_i_atoms{}, charmm_impr_j_atoms{}, charmm_impr_k_atoms{}, charmm_impr_l_atoms{},
    charmm_impr_parameter_indices{}, cmap_i_atoms{}, cmap_j_atoms{}, cmap_k_atoms{},
    cmap_l_atoms{}, cmap_m_atoms{}, cmap_surface_indices{},

    // Valence parameters for CHARMM force field terms
    cmap_surface_dimensions{}, cmap_surface_bounds{}, urey_bradley_stiffnesses{},
    urey_bradley_equilibria{}, charmm_impr_stiffnesses{}, charmm_impr_phase_angles{},
    cmap_surfaces{},

    // Basic force field term and parameter set counts
    bond_term_count{0}, angl_term_count{0}, dihe_term_count{0}, bond_parameter_count{0},
    angl_parameter_count{0}, dihe_parameter_count{0}, attenuated_14_type_count{0},
    inferred_14_attenuations{0},

    // Atom and parameter indexing for basic force field terms
    bond_i_atoms{}, bond_j_atoms{}, bond_parameter_indices{}, angl_i_atoms{}, angl_j_atoms{},
    angl_k_atoms{}, angl_parameter_indices{}, dihe_i_atoms{}, dihe_j_atoms{}, dihe_k_atoms{},
    dihe_l_atoms{}, dihe_parameter_indices{}, dihe14_parameter_indices{}, infr14_i_atoms{},
    infr14_l_atoms{}, infr14_parameter_indices{}, bond_modifiers{}, angl_modifiers{},
    dihe_modifiers{},

    /// Valence parameters for basic force field terms
    bond_stiffnesses{}, bond_equilibria{}, angl_stiffnesses{}, angl_equilibria{},
    dihe_amplitudes{}, dihe_periodicities{}, dihe_phase_angles{},
    
    // Information relevant to virtual site placement
    virtual_site_count{0}, virtual_site_parameter_set_count{0}, virtual_site_atoms{},
    vs_frame_types{}, virtual_site_frame1_atoms{}, virtual_site_frame2_atoms{},
    virtual_site_frame3_atoms{}, virtual_site_frame4_atoms{}, virtual_site_parameter_indices{},
    virtual_site_frame_dim1{}, virtual_site_frame_dim2{}, virtual_site_frame_dim3{},

    // Non-bonded parameter counts
    charge_type_count{0}, lj_type_count{0}, total_exclusions{0},
    periodic_box_class{UnitCellType::NONE}, coulomb_constant{amber_ancient_bioq},
    elec14_screening_factor{amber_default_elec14_screen},
    vdw14_screening_factor{amber_default_vdw14_screen},

    // Arrays of non-bonded parameters
    charge_indices{}, lennard_jones_indices{},
    nb11_exclusion_bounds{std::vector<int>(atom_count + 1, 0)},
    nb11_exclusion_list{},
    nb12_exclusion_bounds{std::vector<int>(atom_count + 1, 0)},
    nb12_exclusion_list{},
    nb13_exclusion_bounds{std::vector<int>(atom_count + 1, 0)},
    nb13_exclusion_list{},
    nb14_exclusion_bounds{std::vector<int>(atom_count + 1, 0)},
    nb14_exclusion_list{},
    atomic_charges{std::vector<double>(atom_count, 0.0)},
    charge_parameters{}, lj_a_values{}, lj_b_values{}, lj_c_values{}, lj_14_a_values{},
    lj_14_b_values{}, lj_14_c_values{}, lj_sigma_values{}, lj_epsilon_values{},
    lj_14_sigma_values{}, lj_14_epsilon_values{}, attn14_tags{},
    atom_overflow_names{}, atom_overflow_types{}, residue_overflow_names{}
{
  // If the residue limits are nonzero and do not match the number of atoms, raise an exception.
  // Otherwise, set the residue limits to one residue encompassing all atoms.
  if (residue_limits.size() > 0) {
    residue_count = residue_limits.size() - 1;
    if (residue_limits[residue_count] != atom_count) {
      rtErr("The upper bound of the residue limits (" +
            std::to_string(residue_limits[residue_count]) + ") much match the atom count.",
            "AtomGraphStage");
    }
    for (int i = 0; i < residue_count; i++) {
      if (residue_limits[i + 1] - residue_limits[i] < 0) {
        rtErr("Residue index " + std::to_string(i) + " holds a negative number of atoms (" +
              std::to_string(residue_limits[i + 1] - residue_limits[i]) + ".", "AtomGraphStage");
      }
    }
  }
  else {
    residue_count = 1;
    residue_limits.resize(2);
    residue_limits[0] = 0;
    residue_limits[1] = atom_count;
  }
}

//-------------------------------------------------------------------------------------------------
AtomGraphStage::AtomGraphStage(const AtomGraph *ag_in, const ExceptionResponse policy_in) :
    AtomGraphStage(ag_in->getAtomCount(), ag_in->getResidueLimits().readHost(), policy_in)
{
  // Load the content from the topology.
  title = ag_in->title;
  output_file = ag_in->source;
  force_fields = ag_in->force_fields;

  // Load information about atoms.
  atomic_numbers = ag_in->atomic_numbers.readHost();
  atom_names = ag_in->atom_names.readHost();
  atom_types = ag_in->atom_types.readHost();
  charge_indices = ag_in->charge_indices.readHost();
  lennard_jones_indices = ag_in->lennard_jones_indices.readHost();
  atomic_charges = ag_in->atomic_charges.readHost();

  // Load information about residues.
  residue_count = ag_in->residue_count;
  residue_limits = ag_in->residue_limits.readHost();
  residue_names = ag_in->residue_names.readHost();

  // Load information about molecules.
  molecule_limits = ag_in->molecule_limits.readHost();
  molecule_membership = ag_in->molecule_membership.readHost();
  molecule_contents = ag_in->molecule_contents.readHost();

  // Load sizing constants for CHARMM force field terms.
  urey_bradley_term_count = ag_in->urey_bradley_term_count;
  charmm_impr_term_count = ag_in->charmm_impr_term_count;
  cmap_term_count = ag_in->cmap_term_count;
  urey_bradley_parameter_count = ag_in->urey_bradley_parameter_count;
  charmm_impr_parameter_count = ag_in->charmm_impr_parameter_count;
  cmap_surface_count = ag_in->cmap_surface_count;

  // Load information about CHARMM force field term indexing.
  urey_bradley_i_atoms = ag_in->urey_bradley_i_atoms.readHost();
  urey_bradley_k_atoms = ag_in->urey_bradley_k_atoms.readHost();
  urey_bradley_parameter_indices = ag_in->urey_bradley_parameter_indices.readHost();
  charmm_impr_i_atoms = ag_in->charmm_impr_i_atoms.readHost();
  charmm_impr_j_atoms = ag_in->charmm_impr_j_atoms.readHost();
  charmm_impr_k_atoms = ag_in->charmm_impr_k_atoms.readHost();
  charmm_impr_l_atoms = ag_in->charmm_impr_l_atoms.readHost();
  charmm_impr_parameter_indices = ag_in->charmm_impr_parameter_indices.readHost();
  cmap_i_atoms = ag_in->cmap_i_atoms.readHost();
  cmap_j_atoms = ag_in->cmap_j_atoms.readHost();
  cmap_k_atoms = ag_in->cmap_k_atoms.readHost();
  cmap_l_atoms = ag_in->cmap_l_atoms.readHost();
  cmap_m_atoms = ag_in->cmap_m_atoms.readHost();
  cmap_surface_indices = ag_in->cmap_surface_indices.readHost();

  // Load information about CHARMM force field parameters
  cmap_surface_dimensions = ag_in->cmap_surface_dimensions.readHost();  
  cmap_surface_bounds = ag_in->cmap_surface_bounds.readHost();
  urey_bradley_stiffnesses = ag_in->urey_bradley_stiffnesses.readHost();
  urey_bradley_equilibria = ag_in->urey_bradley_equilibria.readHost();
  charmm_impr_stiffnesses = ag_in->charmm_impr_stiffnesses.readHost();
  charmm_impr_phase_angles = ag_in->charmm_impr_phase_angles.readHost();
  cmap_surfaces = ag_in->cmap_surfaces.readHost();

  // Load sizing constants for common force field terms.
  bond_term_count = ag_in->bond_term_count;
  angl_term_count = ag_in->angl_term_count;
  dihe_term_count = ag_in->dihe_term_count;
  bond_parameter_count = ag_in->bond_parameter_count;
  angl_parameter_count = ag_in->angl_parameter_count;
  dihe_parameter_count = ag_in->dihe_parameter_count;
  attenuated_14_type_count = ag_in->attenuated_14_type_count;
  inferred_14_attenuations = ag_in->inferred_14_attenuations;

  // Load information about common force field term indexing.
  bond_i_atoms = ag_in->bond_i_atoms.readHost();
  bond_j_atoms = ag_in->bond_j_atoms.readHost();
  bond_parameter_indices = ag_in->bond_parameter_indices.readHost();
  angl_i_atoms = ag_in->angl_i_atoms.readHost();
  angl_j_atoms = ag_in->angl_j_atoms.readHost();
  angl_k_atoms = ag_in->angl_k_atoms.readHost();
  angl_parameter_indices = ag_in->angl_parameter_indices.readHost();
  dihe_i_atoms = ag_in->dihe_i_atoms.readHost();
  dihe_j_atoms = ag_in->dihe_j_atoms.readHost();
  dihe_k_atoms = ag_in->dihe_k_atoms.readHost();
  dihe_l_atoms = ag_in->dihe_l_atoms.readHost();
  dihe_parameter_indices = ag_in->dihe_parameter_indices.readHost();
  infr14_i_atoms = ag_in->infr14_i_atoms.readHost();
  infr14_l_atoms = ag_in->infr14_l_atoms.readHost();
  infr14_parameter_indices = ag_in->infr14_parameter_indices.readHost();
  bond_modifiers = ag_in->bond_modifiers.readHost();
  angl_modifiers = ag_in->angl_modifiers.readHost();
  dihe_modifiers = ag_in->dihe_modifiers.readHost();

  // Load common force field valence parameters.
  bond_stiffnesses = ag_in->bond_stiffnesses.readHost();
  bond_equilibria = ag_in->bond_equilibria.readHost();
  angl_stiffnesses = ag_in->angl_stiffnesses.readHost();
  angl_equilibria = ag_in->angl_equilibria.readHost();
  dihe_amplitudes = ag_in->dihe_amplitudes.readHost();
  dihe_periodicities = ag_in->dihe_periodicities.readHost();
  dihe_phase_angles = ag_in->dihe_phase_angles.readHost();
  
  // Load information about virtual sites.
  virtual_site_count = ag_in->virtual_site_count;
  virtual_site_parameter_set_count = ag_in->virtual_site_parameter_set_count;
  virtual_site_atoms = ag_in->virtual_site_atoms.readHost();
  vs_frame_types.resize(virtual_site_parameter_set_count);
  for (int i = 0; i < virtual_site_parameter_set_count; i++) {
    vs_frame_types[i] = static_cast<VirtualSiteKind>(ag_in->virtual_site_frame_types.readHost(i));
  }
  virtual_site_frame1_atoms = ag_in->virtual_site_frame1_atoms.readHost();
  virtual_site_frame2_atoms = ag_in->virtual_site_frame2_atoms.readHost();
  virtual_site_frame3_atoms = ag_in->virtual_site_frame3_atoms.readHost();
  virtual_site_frame4_atoms = ag_in->virtual_site_frame4_atoms.readHost();
  virtual_site_parameter_indices = ag_in->virtual_site_parameter_indices.readHost();
  virtual_site_frame_dim1 = ag_in->virtual_site_frame_dim1.readHost();
  virtual_site_frame_dim2 = ag_in->virtual_site_frame_dim2.readHost();
  virtual_site_frame_dim3 = ag_in->virtual_site_frame_dim3.readHost();
  
  // Prepare a non-bonded exclusions table and use that to determine the molecular structure(s) of
  // the system.  The data in this workspace object will be transferred to the AtomGraphStage in
  // short order.
  Map1234 nb_excl;
  nb_excl.nb11_excl_bounds = ag_in->nb11_exclusion_bounds.readHost();
  nb_excl.nb12_excl_bounds = ag_in->nb12_exclusion_bounds.readHost();
  nb_excl.nb13_excl_bounds = ag_in->nb13_exclusion_bounds.readHost();
  nb_excl.nb14_excl_bounds = ag_in->nb14_exclusion_bounds.readHost();
  nb_excl.nb11_excl_list = ag_in->nb11_exclusion_list.readHost();
  nb_excl.nb12_excl_list = ag_in->nb12_exclusion_list.readHost();
  nb_excl.nb13_excl_list = ag_in->nb13_exclusion_list.readHost();
  nb_excl.nb14_excl_list = ag_in->nb14_exclusion_list.readHost();
  mapMolecules(atom_count, &molecule_count, nb_excl, &molecule_membership, &molecule_limits,
               &molecule_contents);
  
  // Load physical constants and sizing constants for non-bonded parameters.
  charge_type_count = ag_in->charge_type_count;
  lj_type_count = ag_in->lj_type_count;
  total_exclusions = ag_in->total_exclusions;
  periodic_box_class = ag_in->periodic_box_class;
  coulomb_constant = ag_in->coulomb_constant;
  elec14_screening_factor = ag_in->elec14_screening_factor;
  vdw14_screening_factor = ag_in->vdw14_screening_factor;
  charge_indices = ag_in->charge_indices.readHost();
  lennard_jones_indices = ag_in->lennard_jones_indices.readHost();
  nb11_exclusion_bounds = std::move(nb_excl.nb11_excl_bounds);
  nb11_exclusion_list   = std::move(nb_excl.nb11_excl_list);
  nb12_exclusion_bounds = std::move(nb_excl.nb12_excl_bounds);
  nb12_exclusion_list   = std::move(nb_excl.nb12_excl_list);
  nb13_exclusion_bounds = std::move(nb_excl.nb13_excl_bounds);
  nb13_exclusion_list   = std::move(nb_excl.nb13_excl_list);
  nb14_exclusion_bounds = std::move(nb_excl.nb14_excl_bounds);
  nb14_exclusion_list   = std::move(nb_excl.nb14_excl_list);
  atomic_charges = ag_in->atomic_charges.readHost();
  charge_parameters = ag_in->charge_parameters.readHost();
  lj_a_values = ag_in->lj_a_values.readHost();
  lj_b_values = ag_in->lj_b_values.readHost();
  lj_c_values = ag_in->lj_c_values.readHost();
  lj_14_a_values = ag_in->lj_14_a_values.readHost();
  lj_14_b_values = ag_in->lj_14_b_values.readHost();
  lj_14_c_values = ag_in->lj_14_c_values.readHost();
  lj_sigma_values = ag_in->lj_sigma_values.readHost();
  lj_14_sigma_values = ag_in->lj_14_sigma_values.readHost();
  lj_epsilon_values.resize(lj_sigma_values.size());
  lj_14_epsilon_values.resize(lj_sigma_values.size());
  for (int i = 0; i < lj_type_count; i++) {
    lj_epsilon_values[i] = 0.25 * lj_b_values[i] / pow(lj_sigma_values[i], 6.0);
    lj_14_epsilon_values[i] = 0.25 * lj_14_b_values[i] / pow(lj_14_sigma_values[i], 6.0);
  }

  // The 1:4 screening factors are not taken at face value from a prior topology.  In the staging
  // class, these arrays represent pairs of atom types with specific 1:4 scaling factors to be
  // applied between them.  In the AtomGraph, these arrays represent condensed sets of screening
  // pairs able to represent all dihedral-specific screened 1:4 interactions.

  // Load atom overflow naming terms.
  atom_overflow_names = ag_in->atom_overflow_names.readHost();
  atom_overflow_types = ag_in->atom_overflow_types.readHost();
  residue_overflow_names = ag_in->residue_overflow_names.readHost();
}

//-------------------------------------------------------------------------------------------------
AtomGraphStage::AtomGraphStage(const AtomGraph &ag_in, const ExceptionResponse policy_in) :
    AtomGraphStage(ag_in.getSelfPointer(), policy_in)
{}

//-------------------------------------------------------------------------------------------------
AtomGraphStage::AtomGraphStage(const MdlMol *cmpd_in, const ExceptionResponse policy_in) :
    AtomGraphStage(cmpd_in->getAtomCount(), {}, policy_in)
{}

//-------------------------------------------------------------------------------------------------
AtomGraphStage::AtomGraphStage(const MdlMol &cmpd_in, const ExceptionResponse policy_in) :
    AtomGraphStage(cmpd_in.getSelfPointer(), policy_in)
{}

//-------------------------------------------------------------------------------------------------
const std::string& AtomGraphStage::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
const std::string& AtomGraphStage::getOutputFile() const {
  return output_file;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphStage::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphStage::getResidueCount() const {
  return residue_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphStage::getMoleculeCount() const {
  return molecule_count;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphStage::getElectrostatic14Screening() const {
  return elec14_screening_factor;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphStage::getVanDerWaals14Screening() const {
  return vdw14_screening_factor;
}

//-------------------------------------------------------------------------------------------------
AtomGraph
AtomGraphStage::exportTopology(const int minimum_solute_size,
                               const std::vector<char4> &solute_included_residue_names,
                               const double hmass_repartition_factor) const {
  AtomGraph result;
  time_t rawtime;
  tm *current_time;
  time(&rawtime);
  current_time = gmtime(&rawtime);

  // Set general information
  snprintf(result.version_stamp, 15, "STORMM");
  result.date = *current_time;
  result.title = title;
  result.source = output_file;
  result.force_fields = force_fields;

  // Set information about atoms, residues, and molecules
  result.atom_count = atom_count;
  result.residue_count = residue_count;
  result.molecule_count = molecule_count;
  std::vector<bool> solute_mask = maskSoluteAtoms(minimum_solute_size,
                                                  solute_included_residue_names);
  int biggest_res = 0;
  int last_solute_res = 0;
  int first_solvent_res = residue_count;
  for (int i = 0; i < residue_count; i++) {
    if (solute_mask[i] == false && i < first_solvent_res) {
      first_solvent_res = i;
    }
    if (solute_mask[i] && i > last_solute_res) {
      last_solute_res = i;
    }
    biggest_res = std::max(residue_limits[i + 1] - residue_limits[i], biggest_res);
  }
  result.largest_residue_size = biggest_res;
  result.last_solute_residue = last_solute_res;
  result.last_solute_atom = residue_limits[last_solute_res + 1] - 1;
  result.first_solvent_molecule = (first_solvent_res == residue_count) ?
                                  0 : molecule_membership[residue_limits[first_solvent_res]];
  result.last_atom_before_cap = atom_count - 1;
  result.implicit_copy_count = 0;
  int biggest_mol = 0;
  for (int i = 0; i < molecule_count; i++) {
    biggest_mol = std::max(molecule_limits[i + 1] - molecule_limits[i], biggest_mol);
  }
  result.largest_molecule_size = biggest_mol;
  const std::vector<int> tmp_atom_struc_numbers = incrementingSeries(1, atom_count + 1);
  std::vector<int> tmp_residue_numbers(atom_count);
  for (int i = 0; i < residue_count; i++) {
    for (int j = residue_limits[i]; j < residue_limits[i + 1]; j++) {
      tmp_residue_numbers[j] = i + 1;
    }
  }
  const int mobile_atom_mask_size =
    roundUp<int>(std::max(1, atom_count / (static_cast<int>(sizeof(int)) * 8)), warp_size_int);
  const std::vector<int> tmp_mobile_atoms(mobile_atom_mask_size, -1);
  const std::vector<int> all_atom_zeros(atom_count, 0);
  const std::vector<double> all_atom_dzeros(atom_count, 0);
  const std::vector<char4> all_atom_unknown(atom_count, unknown_atom_name);
  const std::vector<double> tmp_atomic_masses = computeAtomicMasses(hmass_repartition_factor);

  // Construct various tables needed by the constructor by copying some existing arrays.  This
  // cannot be done with std::move in the interest of making this public member function callable
  // in a const object.
  Map1234 tmp_nb_excl;
  tmp_nb_excl.nb11_excl_bounds = nb11_exclusion_bounds;
  tmp_nb_excl.nb12_excl_bounds = nb12_exclusion_bounds;
  tmp_nb_excl.nb13_excl_bounds = nb13_exclusion_bounds;
  tmp_nb_excl.nb14_excl_bounds = nb14_exclusion_bounds;
  tmp_nb_excl.nb11_excl_list = nb11_exclusion_list;
  tmp_nb_excl.nb12_excl_list = nb12_exclusion_list;
  tmp_nb_excl.nb13_excl_list = nb13_exclusion_list;
  tmp_nb_excl.nb14_excl_list = nb14_exclusion_list;
  const CmapAccessories tmp_cmap_tbl = computeCmapDerivatives(cmap_surface_count,
                                                              cmap_surface_dimensions,
                                                              cmap_surface_bounds, cmap_surfaces);
  const CondensedExclusions tmp_cond_excl = calculatePrmtopExclusions(tmp_nb_excl);
  const int zero_excl = countPrmtopZeroExclusions(tmp_nb_excl);
  const std::vector<int> tmp_desc = createPrmtopDescriptors(biggest_res, zero_excl +
                                                            tmp_cond_excl.total_exclusions);
  const BasicValenceTable tmp_bvt(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                                  bond_i_atoms, bond_j_atoms, bond_parameter_indices, angl_i_atoms,
                                  angl_j_atoms, angl_k_atoms, angl_parameter_indices, dihe_i_atoms,
                                  dihe_j_atoms, dihe_k_atoms, dihe_l_atoms,
                                  dihe_parameter_indices);
  const CharmmValenceTable tmp_cvt(atom_count, urey_bradley_term_count, charmm_impr_term_count,
                                   cmap_surface_count, urey_bradley_i_atoms, urey_bradley_k_atoms,
                                   urey_bradley_parameter_indices, charmm_impr_i_atoms,
                                   charmm_impr_j_atoms, charmm_impr_k_atoms, charmm_impr_l_atoms,
                                   charmm_impr_parameter_indices, cmap_i_atoms, cmap_j_atoms,
                                   cmap_k_atoms, cmap_l_atoms, cmap_m_atoms, cmap_surface_indices);
  std::vector<double> attn14_elec_factors, attn14_vdw_factors;
  buildNonbonded14Screen(&attn14_elec_factors, &attn14_vdw_factors);
  const AttenuationParameterSet tmp_attn_parm = condenseScreeningFactors(tmp_bvt,
                                                                         attn14_elec_factors,
                                                                         attn14_vdw_factors,
                                                                         elec14_screening_factor,
                                                                         vdw14_screening_factor);
  VirtualSiteTable tmp_vs_tbl;
  tmp_vs_tbl.vs_count = virtual_site_count;
  tmp_vs_tbl.vs_numbers = std::vector<int>(atom_count, -1);
  tmp_vs_tbl.frame_types = std::vector<int>(virtual_site_parameter_set_count);
  for (int i = 0; i < virtual_site_parameter_set_count; i++) {
    tmp_vs_tbl.vs_numbers[virtual_site_atoms[i]] = i;
    tmp_vs_tbl.frame_types[i] = static_cast<int>(vs_frame_types[i]);
  }
  tmp_vs_tbl.vs_atoms = virtual_site_atoms;
  tmp_vs_tbl.frame1_atoms = virtual_site_frame1_atoms;
  tmp_vs_tbl.frame2_atoms = virtual_site_frame2_atoms;
  tmp_vs_tbl.frame3_atoms = virtual_site_frame3_atoms;
  tmp_vs_tbl.frame4_atoms = virtual_site_frame4_atoms;
  tmp_vs_tbl.param_idx = virtual_site_parameter_indices;
  tmp_vs_tbl.frame_dim1 = virtual_site_frame_dim1;
  tmp_vs_tbl.frame_dim2 = virtual_site_frame_dim2;
  tmp_vs_tbl.frame_dim3 = virtual_site_frame_dim3;
  const ConstraintTable tmp_cnst_tbl(atomic_numbers, tmp_atomic_masses, molecule_limits,
                                     molecule_contents, molecule_membership, tmp_bvt, tmp_nb_excl,
                                     bond_equilibria, angl_equilibria);

  // Various sizing constants of the resulting AtomGraph must be set explicitly, even though they
  // are part of the various objects that were just constructed.  These are the sizing constants
  // for valence terms.
  result.urey_bradley_term_count = urey_bradley_term_count;
  result.charmm_impr_term_count = charmm_impr_term_count;
  result.cmap_term_count = cmap_term_count;
  result.urey_bradley_parameter_count = urey_bradley_parameter_count;
  result.charmm_impr_parameter_count = charmm_impr_parameter_count;
  result.cmap_surface_count = cmap_surface_count;
  result.bond_term_count = bond_term_count;
  result.angl_term_count = angl_term_count;
  result.dihe_term_count = dihe_term_count;
  result.bond_parameter_count = bond_parameter_count;
  result.angl_parameter_count = angl_parameter_count;
  result.dihe_parameter_count = dihe_parameter_count;
  result.attenuated_14_type_count = tmp_attn_parm.total_14_sets;
  result.inferred_14_attenuations = inferred_14_attenuations;

  // Set sizing constants for virtual sites.
  result.virtual_site_count = virtual_site_count;
  result.virtual_site_parameter_set_count = virtual_site_parameter_set_count;

  // Set sizing constants for non-bonded terms.
  result.charge_type_count = charge_type_count;
  result.lj_type_count = lj_type_count;
  result.total_exclusions = tmp_cond_excl.total_exclusions + zero_excl;
  result.periodic_box_class = periodic_box_class;
  result.coulomb_constant = coulomb_constant;
  result.elec14_screening_factor = elec14_screening_factor;
  result.vdw14_screening_factor = vdw14_screening_factor;

  // Perturbations of the resulting system are not handled.  These sizing constants were set to
  // zero when the result was initialized.  However, the AtomGraph does track the number of bond
  // terms with or without hydrogen outside of the tmp_desc array that was just constructed.
  const size_t bwh_idx  = static_cast<size_t>(TopologyDescriptor::BONDS_WITH_HYDROGEN);
  const size_t awh_idx  = static_cast<size_t>(TopologyDescriptor::ANGLES_WITH_HYDROGEN);
  const size_t hwh_idx  = static_cast<size_t>(TopologyDescriptor::DIHEDRALS_WITH_HYDROGEN);
  const size_t bwoh_idx = static_cast<size_t>(TopologyDescriptor::BONDS_WITHOUT_HYDROGEN);
  const size_t awoh_idx = static_cast<size_t>(TopologyDescriptor::ANGLES_WITHOUT_HYDROGEN);
  const size_t hwoh_idx = static_cast<size_t>(TopologyDescriptor::DIHEDRALS_WITHOUT_HYDROGEN);
  result.bond_term_with_hydrogen = tmp_desc[bwh_idx];
  result.angl_term_with_hydrogen = tmp_desc[awh_idx];
  result.dihe_term_with_hydrogen = tmp_desc[hwh_idx];
  result.bond_term_without_hydrogen = tmp_desc[bwoh_idx];
  result.angl_term_without_hydrogen = tmp_desc[awoh_idx];
  result.dihe_term_without_hydrogen = tmp_desc[hwoh_idx];
  result.loadHybridArrays(tmp_desc, residue_limits, tmp_atom_struc_numbers, tmp_residue_numbers,
                          molecule_limits, atomic_numbers, molecule_membership, tmp_mobile_atoms,
                          molecule_contents,  cmap_surface_dimensions, cmap_surface_bounds,
                          charge_indices, lennard_jones_indices, infr14_i_atoms, infr14_l_atoms,
                          infr14_parameter_indices, all_atom_zeros, all_atom_zeros, all_atom_zeros,
                          atomic_charges, tmp_atomic_masses, urey_bradley_stiffnesses,
                          urey_bradley_equilibria, charmm_impr_stiffnesses,
                          charmm_impr_phase_angles, cmap_surfaces, bond_stiffnesses,
                          bond_equilibria, angl_stiffnesses, angl_equilibria, dihe_amplitudes,
                          dihe_periodicities, dihe_phase_angles, charge_parameters, lj_a_values,
                          lj_b_values, lj_c_values, lj_14_a_values, lj_14_b_values, lj_14_c_values,
                          all_atom_dzeros, all_atom_dzeros, all_atom_dzeros, all_atom_dzeros, {},
                          {}, {}, atom_names, atom_types, residue_names, all_atom_unknown,
                          tmp_cmap_tbl, tmp_cond_excl, tmp_bvt, tmp_cvt, tmp_attn_parm, tmp_vs_tbl,
                          tmp_nb_excl, tmp_cnst_tbl);
  return result;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::setTitle(const std::string &title_in) {
  if (title_in.size() > 80) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The topology title can be at most 80 characters, but " +
            std::to_string(title_in.size()) + " were provided.", "AtomGraphStage", "setTitle");
    case ExceptionResponse::WARN:
      rtWarn("The topology title can be at most 80 characters, but " +
             std::to_string(title_in.size()) + " were provided.  Excess characters will be "
             "trimmed.", "AtomGraphStage", "setTitle");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  title = title_in;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::setOutputFile(const std::string &output_file_in) {
  output_file = output_file_in;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::addUbrdParameters(const double k_eq_in, const double l_eq_in,
                                       const char4 atom_i_type, const char4 atom_k_type) {

  // If the atom types are not unknown or wildcards, check that no other bond parameter with these
  // atom types has already been specified.
  int target_index = urey_bradley_parameter_count;
  if ((atom_i_type != unknown_atom_type && atom_i_type != wildcard_atom_type) ||
      (atom_k_type != unknown_atom_type && atom_k_type != wildcard_atom_type)) {
    for (int pos = 0; pos < bond_parameter_count; pos++) {
      if ((urey_bradley_i_atom_types[pos] == atom_i_type &&
           urey_bradley_k_atom_types[pos] == atom_k_type) ||
          (urey_bradley_i_atom_types[pos] == atom_k_type &&
           urey_bradley_k_atom_types[pos] == atom_i_type)) {
        const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
        const std::string msg("A Urey-Bradley stretching term between atom types " +
                              char4ToString(urey_bradley_i_atom_types[pos]) + " and " +
                              char4ToString(urey_bradley_k_atom_types[pos]) + " has already been "
                              "specified.  Old equilibrium length / stiffness: " +
                              realToString(urey_bradley_equilibria[pos], 9, 4, nfmt) + " A / " +
                              realToString(urey_bradley_stiffnesses[pos], 9, 4, nfmt) +
                              " kcal/mol-A^2, new equilibrium length / stiffness: " +
                              realToString(l_eq_in, 9, 4, nfmt) + " A / " +
                              realToString(k_eq_in, 9, 4, nfmt) + " kcal/mol-A^2.");
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(msg, "AtomGraphStage", "addUbrdParameters");
        case ExceptionResponse::WARN:
          rtWarn(msg + "  The latest definition will take precedence.", "AtomGraphStage",
                 "addUbrdParameters");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        target_index = pos;
      }
    }
  }
  if (target_index == urey_bradley_parameter_count) {
    urey_bradley_stiffnesses.push_back(k_eq_in);
    urey_bradley_equilibria.push_back(l_eq_in);
    urey_bradley_i_atom_types.push_back(atom_i_type);
    urey_bradley_k_atom_types.push_back(atom_k_type);
  }
  else {
    urey_bradley_stiffnesses[target_index]  = k_eq_in;
    urey_bradley_equilibria[target_index]   = l_eq_in;
    urey_bradley_i_atom_types[target_index] = atom_i_type;
    urey_bradley_k_atom_types[target_index] = atom_k_type;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::addCImpParameters(const double k_eq_in, const double phi_in,
                                       const char4 atom_i_type, const char4 atom_j_type,
                                       const char4 atom_k_type, const char4 atom_l_type) {

  // If the atom types are not unknown or wildcards, check that no other bond parameter with these
  // atom types has already been specified.
  int target_index = dihe_parameter_count;
  if ((atom_i_type != unknown_atom_type && atom_i_type != wildcard_atom_type) ||
      (atom_j_type != unknown_atom_type && atom_j_type != wildcard_atom_type) ||
      (atom_k_type != unknown_atom_type && atom_k_type != wildcard_atom_type) ||
      (atom_l_type != unknown_atom_type && atom_l_type != wildcard_atom_type)) {
    for (int pos = 0; pos < charmm_impr_parameter_count; pos++) {
      if ((charmm_impr_i_atom_types[pos] == atom_i_type &&
           charmm_impr_j_atom_types[pos] == atom_j_type &&
           charmm_impr_k_atom_types[pos] == atom_k_type &&
           charmm_impr_l_atom_types[pos] == atom_l_type) ||
          (charmm_impr_i_atom_types[pos] == atom_l_type &&
           charmm_impr_j_atom_types[pos] == atom_k_type &&
           charmm_impr_k_atom_types[pos] == atom_j_type &&
           charmm_impr_l_atom_types[pos] == atom_i_type)) {
        const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
        const std::string err_msg("A CHARMM improper dihedral between atom types " +
                                  char4ToString(charmm_impr_i_atom_types[pos]) + ", " + 
                                  char4ToString(charmm_impr_j_atom_types[pos]) + ", " +
                                  char4ToString(charmm_impr_k_atom_types[pos]) + ", and " +
                                  char4ToString(charmm_impr_l_atom_types[pos]) + ", has already "
                                  "been specified.  Old stiffness / phase angle : " +
                                  realToString(charmm_impr_stiffnesses[pos], 9, 4, nfmt) +
                                  " kcal/mol-rad^2 / " +
                                  realToString(charmm_impr_phase_angles[pos], 9, 4, nfmt) +
                                  " rad, new stiffness / phase angle : " +
                                  realToString(k_eq_in, 9, 4, nfmt) + " kcal/mol-rad^2 / " +
                                  realToString(phi_in, 9, 4, nfmt) + " rad.");
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(err_msg, "AtomGraphStage", "addCImpParameters");
        case ExceptionResponse::WARN:
          rtWarn(err_msg + "  The latest definition will take precedence.", "AtomGraphStage",
                 "addCImpParameters");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        target_index = pos;
      }
    }
  }
  if (target_index == charmm_impr_parameter_count) {
    charmm_impr_stiffnesses.push_back(k_eq_in);
    charmm_impr_phase_angles.push_back(phi_in);
    charmm_impr_i_atom_types.push_back(atom_i_type);
    charmm_impr_j_atom_types.push_back(atom_j_type);
    charmm_impr_k_atom_types.push_back(atom_k_type);
    charmm_impr_l_atom_types.push_back(atom_l_type);
  }
  else {
    charmm_impr_stiffnesses[target_index]  = k_eq_in;
    charmm_impr_phase_angles[target_index] = phi_in;
    charmm_impr_i_atom_types[target_index] = atom_i_type;
    charmm_impr_j_atom_types[target_index] = atom_j_type;
    charmm_impr_k_atom_types[target_index] = atom_k_type;
    charmm_impr_l_atom_types[target_index] = atom_l_type;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::addCmapParameters(const int dimension_in, const std::vector<double> &surf_in,
                                       const char4 atom_i_type, const char4 atom_j_type,
                                       const char4 atom_k_type, const char4 atom_l_type,
                                       const char4 atom_m_type) {

  // If the atom types are not unknown or wildcards, check that no other bond parameter with these
  // atom types has already been specified.
  int target_index = dihe_parameter_count;
  if ((atom_i_type != unknown_atom_type && atom_i_type != wildcard_atom_type) ||
      (atom_j_type != unknown_atom_type && atom_j_type != wildcard_atom_type) ||
      (atom_k_type != unknown_atom_type && atom_k_type != wildcard_atom_type) ||
      (atom_l_type != unknown_atom_type && atom_l_type != wildcard_atom_type) ||
      (atom_m_type != unknown_atom_type && atom_m_type != wildcard_atom_type)) {
    for (int pos = 0; pos < cmap_surface_count; pos++) {
      if ((cmap_i_atom_types[pos] == atom_i_type && cmap_j_atom_types[pos] == atom_j_type &&
           cmap_k_atom_types[pos] == atom_k_type && cmap_l_atom_types[pos] == atom_l_type &&
           cmap_m_atom_types[pos] == atom_m_type) ||
          (cmap_i_atom_types[pos] == atom_m_type && cmap_j_atom_types[pos] == atom_l_type &&
           cmap_k_atom_types[pos] == atom_k_type && cmap_l_atom_types[pos] == atom_j_type &&
           cmap_m_atom_types[pos] == atom_i_type)) {
        const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
        const std::string err_msg("A CHARMM improper dihedral between atom types " +
                                  char4ToString(cmap_i_atom_types[pos]) + ", " + 
                                  char4ToString(cmap_j_atom_types[pos]) + ", " +
                                  char4ToString(cmap_k_atom_types[pos]) + ", " +
                                  char4ToString(cmap_l_atom_types[pos]) + ", and " +
                                  char4ToString(cmap_l_atom_types[pos]) + ", has already "
                                  "been specified.  Old dimension: " +
                                  std::to_string(cmap_surface_dimensions[pos]) +
                                  ", new dimension " + std::to_string(dimension_in) + ".");
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(err_msg, "AtomGraphStage", "addCmapParameters");
        case ExceptionResponse::WARN:
          rtWarn(err_msg + "  The latest definition will take precedence.", "AtomGraphStage",
                 "addCmapParameters");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        target_index = pos;
      }
    }
  }
  if (target_index == cmap_surface_count) {
    cmap_surface_dimensions.push_back(dimension_in);
    if (cmap_surface_count == 0) {
      cmap_surface_bounds.push_back(0);
    }
    cmap_surface_bounds.push_back(dimension_in * dimension_in);
    cmap_surfaces.insert(cmap_surfaces.end(), surf_in.begin(), surf_in.end());
    cmap_i_atom_types.push_back(atom_i_type);
    cmap_j_atom_types.push_back(atom_j_type);
    cmap_k_atom_types.push_back(atom_k_type);
    cmap_l_atom_types.push_back(atom_l_type);
    cmap_m_atom_types.push_back(atom_m_type);
  }
  else {

    // Rebuild the CMAP surfaces and bounds arrays
    std::vector<int> tmp_cmap_surface_dimensions = cmap_surface_dimensions;
    std::vector<int> tmp_cmap_surface_bounds(cmap_surface_count + 1, 0);
    for (int i = 0; i < cmap_surface_count; i++) {
      tmp_cmap_surface_bounds[i] = (cmap_surface_dimensions[i] * cmap_surface_dimensions[i]);
    }
    tmp_cmap_surface_dimensions[target_index] = dimension_in;
    tmp_cmap_surface_bounds[target_index] = dimension_in * dimension_in;
    prefixSumInPlace(&tmp_cmap_surface_bounds, PrefixSumType::EXCLUSIVE);
    std::vector<double> tmp_cmap_surfaces(tmp_cmap_surface_bounds[cmap_surface_count]);
    int ic = 0;
    for (int i = 0; i < cmap_surface_count; i++) {
      if (i == target_index) {
        for (int j = 0; j < dimension_in * dimension_in; j++) {
          tmp_cmap_surfaces[ic] = surf_in[j];
          ic++;
        }
      }
      else {
        const int i_offset = cmap_surface_bounds[i];
        const int idim_sq = cmap_surface_dimensions[i] * cmap_surface_dimensions[i];
        for (int j = 0; j < idim_sq; j++) {
          tmp_cmap_surfaces[ic] = cmap_surfaces[i_offset + j];
          ic++;
        }
      }
    }
    cmap_surface_dimensions = std::move(tmp_cmap_surface_dimensions);
    cmap_surface_bounds = std::move(tmp_cmap_surface_bounds);
    cmap_surfaces = std::move(tmp_cmap_surfaces);
    cmap_i_atom_types[target_index] = atom_i_type;
    cmap_j_atom_types[target_index] = atom_j_type;
    cmap_k_atom_types[target_index] = atom_k_type;
    cmap_l_atom_types[target_index] = atom_l_type;
    cmap_m_atom_types[target_index] = atom_m_type;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::addBondParameters(const double k_eq_in, const double l_eq_in,
                                       const char4 atom_i_type, const char4 atom_j_type) {

  // If the atom types are not unknown or wildcards, check that no other bond parameter with these
  // atom types has already been specified.
  int target_index = bond_parameter_count;
  if ((atom_i_type != unknown_atom_type && atom_i_type != wildcard_atom_type) ||
      (atom_j_type != unknown_atom_type && atom_j_type != wildcard_atom_type)) {
    for (int pos = 0; pos < bond_parameter_count; pos++) {
      if ((bond_i_atom_types[pos] == atom_i_type && bond_j_atom_types[pos] == atom_j_type) ||
          (bond_i_atom_types[pos] == atom_j_type && bond_j_atom_types[pos] == atom_i_type)) {
        const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
        const std::string msg("A bond between atom types " +
                              char4ToString(bond_i_atom_types[pos]) + " and " +
                              char4ToString(bond_j_atom_types[pos]) + " has already been "
                              "specified.  Old equilibrium length / stiffness: " +
                              realToString(bond_equilibria[pos], 9, 4, nfmt) + " A / " +
                              realToString(bond_stiffnesses[pos], 9, 4, nfmt) +
                              " kcal/mol-A^2, new equilibrium length / stiffness: " +
                              realToString(l_eq_in, 9, 4, nfmt) + " A / " +
                              realToString(k_eq_in, 9, 4, nfmt) + " kcal/mol-A^2.");
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(msg, "AtomGraphStage", "addBondParameters");
        case ExceptionResponse::WARN:
          rtWarn(msg + "  The latest definition will take precedence.", "AtomGraphStage",
                 "addBondParameters");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        target_index = pos;
      }
    }
  }
  if (target_index == bond_parameter_count) {
    bond_stiffnesses.push_back(k_eq_in);
    bond_equilibria.push_back(l_eq_in);
    bond_i_atom_types.push_back(atom_i_type);
    bond_j_atom_types.push_back(atom_j_type);
  }
  else {
    bond_stiffnesses[target_index]  = k_eq_in;
    bond_equilibria[target_index]   = l_eq_in;
    bond_i_atom_types[target_index] = atom_i_type;
    bond_j_atom_types[target_index] = atom_j_type;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::addAnglParameters(const double k_eq_in, const double th_eq_in,
                                       const char4 atom_i_type, const char4 atom_j_type,
                                       const char4 atom_k_type) {

  // If the atom types are not unknown or wildcards, check that no other bond parameter with these
  // atom types has already been specified.
  int target_index = angl_parameter_count;
  if ((atom_i_type != unknown_atom_type && atom_i_type != wildcard_atom_type) ||
      (atom_j_type != unknown_atom_type && atom_j_type != wildcard_atom_type) ||
      (atom_k_type != unknown_atom_type && atom_k_type != wildcard_atom_type)) {
    for (int pos = 0; pos < angl_parameter_count; pos++) {
      if ((angl_i_atom_types[pos] == atom_i_type && angl_j_atom_types[pos] == atom_j_type &&
           angl_k_atom_types[pos] == atom_k_type) ||
          (angl_i_atom_types[pos] == atom_k_type && angl_j_atom_types[pos] == atom_j_type &&
           angl_k_atom_types[pos] == atom_i_type)) {
        const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
        const std::string err_msg("An angle between atom types " +
                                  char4ToString(angl_i_atom_types[pos]) + ", " + 
                                  char4ToString(angl_j_atom_types[pos]) + ", and " +
                                  char4ToString(angl_k_atom_types[pos]) + ", has already been "
                                  "specified.  Old equilibrium value / stiffness: " +
                                  realToString(angl_equilibria[pos], 9, 4, nfmt) + " rad / " +
                                  realToString(angl_stiffnesses[pos], 9, 4, nfmt) +
                                  " kcal/mol-rad^2, new equilibrium value / stiffness: " +
                                  realToString(th_eq_in, 9, 4, nfmt) + " rad / " +
                                  realToString(k_eq_in, 9, 4, nfmt) + " kcal/mol-rad^2.");
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(err_msg, "AtomGraphStage", "addAnglParameters");
        case ExceptionResponse::WARN:
          rtWarn(err_msg + "  The latest definition will take precedence.", "AtomGraphStage",
                 "addAnglParameters");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        target_index = pos;
      }
    }
  }
  if (target_index == bond_parameter_count) {
    angl_stiffnesses.push_back(k_eq_in);
    angl_equilibria.push_back(th_eq_in);
    angl_i_atom_types.push_back(atom_i_type);
    angl_j_atom_types.push_back(atom_j_type);
    angl_k_atom_types.push_back(atom_k_type);
  }
  else {
    angl_stiffnesses[target_index]  = k_eq_in;
    angl_equilibria[target_index]   = th_eq_in;
    angl_i_atom_types[target_index] = atom_i_type;
    angl_j_atom_types[target_index] = atom_j_type;
    angl_k_atom_types[target_index] = atom_k_type;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::addDiheParameters(const double ampl_in, const double period_in,
                                       const double phi_in, const char4 atom_i_type,
                                       const char4 atom_j_type, const char4 atom_k_type,
                                       const char4 atom_l_type) {

  // If the atom types are not unknown or wildcards, check that no other bond parameter with these
  // atom types has already been specified.
  int target_index = dihe_parameter_count;
  if ((atom_i_type != unknown_atom_type && atom_i_type != wildcard_atom_type) ||
      (atom_j_type != unknown_atom_type && atom_j_type != wildcard_atom_type) ||
      (atom_k_type != unknown_atom_type && atom_k_type != wildcard_atom_type) ||
      (atom_l_type != unknown_atom_type && atom_l_type != wildcard_atom_type)) {
    for (int pos = 0; pos < dihe_parameter_count; pos++) {
      if ((dihe_i_atom_types[pos] == atom_i_type && dihe_j_atom_types[pos] == atom_j_type &&
           dihe_k_atom_types[pos] == atom_k_type && dihe_l_atom_types[pos] == atom_l_type) ||
          (dihe_i_atom_types[pos] == atom_l_type && dihe_j_atom_types[pos] == atom_k_type &&
           dihe_k_atom_types[pos] == atom_j_type && dihe_l_atom_types[pos] == atom_i_type)) {
        const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
        const std::string err_msg("A dihedral between atom types " +
                                  char4ToString(dihe_i_atom_types[pos]) + ", " + 
                                  char4ToString(dihe_j_atom_types[pos]) + ", " +
                                  char4ToString(dihe_k_atom_types[pos]) + ", and " +
                                  char4ToString(dihe_l_atom_types[pos]) + ", has already been "
                                  "specified.  Old amplitude / periodicity / phase angle : " +
                                  realToString(dihe_amplitudes[pos], 9, 4, nfmt) + " kcal/mol / " +
                                  realToString(dihe_periodicities[pos], 9, 4, nfmt) + " / " +
                                  realToString(dihe_phase_angles[pos], 9, 4, nfmt) + " rad, new "
                                  "amplitude / periodicity / phase angle : " +
                                  realToString(ampl_in, 9, 4, nfmt) + " kcal/mol / " +
                                  realToString(period_in, 9, 4, nfmt) + " / " +
                                  realToString(phi_in, 9, 4, nfmt) + " rad.");
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(err_msg, "AtomGraphStage", "addDiheParameters");
        case ExceptionResponse::WARN:
          rtWarn(err_msg + "  The latest definition will take precedence.", "AtomGraphStage",
                 "addDiheParameters");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        target_index = pos;
      }
    }
  }
  if (target_index == bond_parameter_count) {
    dihe_amplitudes.push_back(ampl_in);
    dihe_periodicities.push_back(period_in);
    dihe_phase_angles.push_back(phi_in);
    dihe_i_atom_types.push_back(atom_i_type);
    dihe_j_atom_types.push_back(atom_j_type);
    dihe_k_atom_types.push_back(atom_k_type);
    dihe_l_atom_types.push_back(atom_l_type);
  }
  else {
    dihe_amplitudes[target_index] = ampl_in;
    dihe_periodicities[target_index] = period_in;
    dihe_phase_angles[target_index] = phi_in;
    dihe_i_atom_types[target_index] = atom_i_type;
    dihe_j_atom_types[target_index] = atom_j_type;
    dihe_k_atom_types[target_index] = atom_k_type;
    dihe_l_atom_types[target_index] = atom_l_type;
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraphStage::addAtoms(const std::vector<int> &z_numbers, const int placement,
                                          const std::vector<int> &added_residue_limits,
                                          const std::vector<char4> &added_atom_names,
                                          const std::vector<char4> &added_atom_types,
                                          const std::vector<char4> &added_residue_names,
                                          const std::vector<int> &added_charge_indices,
                                          const std::vector<int> &added_lennard_jones_indices,
                                          const std::vector<double> &added_atomic_charges) {

  // Return immediately if there is nothing to do.
  if (z_numbers.size() == 0) {
    return std::vector<int>();
  }
  
  // Check that no atoms are virtual sites.
  const int added_count = z_numbers.size();
  for (int i = 0; i < added_count; i++) {
    if (z_numbers[i] == 0) {
      rtErr("Virtual sites can only be added with the addVirtualSites() member function.",
            "AtomGraphStage", "addAtoms");
    }
  }
  
  // When adding atoms within a residue, the added residue limits must be a blank array and all
  // added atoms will go into the existing residue.  Otherwise, adding residues at the back of the
  // list or in between other residues will create new residues with the stated limits, adapted
  // to fit the topology.
  const int added_residue_count = std::max(static_cast<int>(added_residue_limits.size()) - 1, 1);
  const int updated_residue_count = residue_count + added_residue_count;
  int insert_ridx;
  if (placement < 0 || placement >= atom_count) {

    // The atoms are being added to the back of the list.  Add the correct number of residues,
    // each with zero atoms, and fill out the array of residue homes for each added atom.
    insert_ridx = residue_count;
  }
  else {
    insert_ridx = findBin(residue_limits, placement, ExceptionResponse::DIE);
    if (placement != residue_limits[insert_ridx]) {

      // The atoms are being added inside of the one of the existing residues.  This will not
      // permit the construction of new residues.
      if (added_residue_limits.size() > 0) {
        rtErr("Adding atoms within residue " + std::to_string(insert_ridx) + " does not accept "
              "new residue limits.", "AtomGraphStage", "addAtoms");
      }
    }
    else {

      // The atoms being added will go between residues.  Adjust the residue arrays now.  If no
      // residue limits are provided for the new atoms, assume that they are all one residue.
      std::vector<int> tmp_res_lims(updated_residue_count + 1, 0);
      for (int i = 0; i < insert_ridx; i++) {
        tmp_res_lims[i] = residue_limits[i + 1] - residue_limits[i];
      }
      for (int i = 0; i < added_residue_count; i++) {

        // The new residues will obtain their atom counts in the insertParticles() function.  Set
        // the atom counts to zero until then.
        tmp_res_lims[insert_ridx + i] = 0;
      }
      for (int i = insert_ridx + added_residue_count; i < updated_residue_count; i++) {
        const int j = i - added_residue_count;
        tmp_res_lims[i] = residue_limits[j + 1] - residue_limits[j];
      }
      prefixSumInPlace(&tmp_res_lims, PrefixSumType::EXCLUSIVE);
      residue_limits = tmp_res_lims;

      // Update the residue names array.
      residue_names.resize(updated_residue_count);
      for (int i = updated_residue_count - 1; i >= insert_ridx + added_residue_count; i--) {
        residue_names[i] = residue_names[i - added_residue_count];
      }
      if (added_residue_names.size() == added_residue_count) {
        for (int i = insert_ridx; i < insert_ridx + added_residue_count; i++) {
          residue_names[i] = added_residue_names[i - insert_ridx];
        }
      }
      else {
        for (int i = insert_ridx; i < insert_ridx + added_residue_count; i++) {
          residue_names[i] = unknown_residue_name;
        }
      }

      // Update the residue count in the topology.
      residue_count = updated_residue_count;
    }
  }

  // Fill out the residue homes of the new atoms.
  std::vector<int> residue_homes(added_count);
  if (added_residue_count == 1) {
    for (int i = 0; i < added_count; i++) {
      residue_homes[i] = insert_ridx;
    }
  }
  else {
    for (int i = 0; i < added_residue_count; i++) {
      for (int j = added_residue_limits[i]; j < added_residue_limits[i + 1]; j++) {
        residue_homes[j] = insert_ridx + i;
      }
    }
  }
  
  // Insert particles.  This will increase the added residues' atom counts as intended, update any
  // existing connectivity and valence term indexing, and move various atom properties to make room
  // for the new atoms.
  const std::vector<int> atom_positions = insertParticles(z_numbers, residue_homes);

  // Check atom property inputs to ensure they are either blank or account for all added atoms.
  const bool apply_names   = applyProperty(added_atom_names.size(), added_count, "names");
  const bool apply_types   = applyProperty(added_atom_types.size(), added_count, "types");
  const bool apply_qidx    = applyProperty(added_charge_indices.size(), added_count,
                                           "charge indices");
  const bool apply_ljidx   = applyProperty(added_lennard_jones_indices.size(), added_count,
                                           "Lennard-Jones indices");
  const bool apply_charges = applyProperty(added_atomic_charges.size(), added_count,
                                           "partial charges");

  // Add the atom names and properties to the topology.
  for (int i = 0; i < added_count; i++) {
    const int posi = atom_positions[i];
    if (apply_names) {
      atom_names[posi] = added_atom_names[i];
    }
    else {
      atom_names[posi] = unknown_atom_name;
    }
    if (apply_types) {
      atom_types[posi] = added_atom_types[i];
    }
    else {
      atom_types[posi] = unknown_atom_type;
    }
    charge_indices[posi]        = (apply_qidx)    ? added_charge_indices[i] : -1;
    lennard_jones_indices[posi] = (apply_ljidx)   ? added_lennard_jones_indices[i] : -1;
    atomic_charges[posi]        = (apply_charges) ? added_atomic_charges[i] : 0.0;
  }
  return atom_positions;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphStage::addAtom(const int z_number, const int placement, const char4 atom_name,
                            const char4 atom_type, const char4 res_name, const int charge_index,
                            const int lennard_jones_index, const double atomic_charge) {
  const std::vector<int> z_numbers(1, z_number);
  const std::vector<int> tmp_res_lims = { 0, 1 };
  const std::vector<char4> added_atom_names(1, atom_name);
  const std::vector<char4> added_atom_types(1, atom_type);
  const std::vector<char4> added_res_names(1, res_name);
  const std::vector<int> added_charge_indices(1, charge_index);
  const std::vector<int> added_lennard_jones_indices(1, lennard_jones_index);
  const std::vector<double> added_atomic_charges(1, atomic_charge);
  const std::vector<int> result = addAtoms(z_numbers, placement, tmp_res_lims, added_atom_names,
                                           added_atom_types, added_res_names, added_charge_indices,
                                           added_lennard_jones_indices, added_atomic_charges);
  return result[0];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraphStage::addVirtualSites(const std::vector<int> &parameter_indices,
                                                 const std::vector<int> &frame1_atoms,
                                                 const std::vector<int> &frame2_atoms,
                                                 const std::vector<int> &frame3_atoms,
                                                 const std::vector<int> &frame4_atoms) {
  if (parameter_indices.size() != frame1_atoms.size() ||
      parameter_indices.size() != frame2_atoms.size()) {
    rtErr("At least two frame atoms must be provided for every virtual site frame (" +
          std::to_string(parameter_indices.size()) + " frames, " +
          std::to_string(frame1_atoms.size()) + " parent atoms, and " +
          std::to_string(frame2_atoms.size()) + " secondary frame atoms provided.",
          "AtomGraphStage", "addVirtualSites");
  }
  const int nvs = parameter_indices.size();
  for (int i = 0; i < nvs; i++) {
    validateFrameIndex(parameter_indices[i]);
    switch (vs_frame_types[parameter_indices[i]]) {
    case VirtualSiteKind::NONE:
      rtErr("A virtual site cannot be specified with no frame type.");
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      if (frame3_atoms.size() <= nvs) {
        rtErr("A sufficient number of third frame atoms must be provided, even if set to "
              "placeholder values, to provide a third frame atom for frame type " +
              getEnumerationName(vs_frame_types[parameter_indices[i]]) + " in element " +
              std::to_string(i) + " of the requested set of virtual sites.", "AtomGraphStage",
              "addVirtualSites");
      }
      break;
    case VirtualSiteKind::FIXED_4:
      if (frame3_atoms.size() <= nvs || frame4_atoms.size() <= nvs) {
        rtErr("A sufficient number of third and fourth frame atoms must be provided, even if set "
              "to placeholder values, to provide a third frame atom for frame type " +
              getEnumerationName(vs_frame_types[parameter_indices[i]]) + " in element " +
              std::to_string(i) + " of the requested set of virtual sites.", "AtomGraphStage",
              "addVirtualSites");
      }
      break;
    }
  }
  std::vector<int> vs_residue_homes(nvs);
  for (int i = 0; i < nvs; i++) {
    validateAtomIndex(frame1_atoms[i]);
    validateAtomIndex(frame2_atoms[i]);
    validateFrameIndex(parameter_indices[i]);
    vs_residue_homes[i] = findBin(residue_limits, frame1_atoms[i], ExceptionResponse::DIE);
    switch (vs_frame_types[parameter_indices[i]]) {
    case VirtualSiteKind::NONE:
      rtErr("A virtual site cannot be specified with no frame type.");
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      if (frame3_atoms[i] >= 0 || frame4_atoms[i] >= 0) {
        rtErr("Frame type " + getEnumerationName(vs_frame_types[parameter_indices[i]]) +
              " has only two frame atoms, but more were specified.", "AtomGraphStage",
              "addVirtualSite");
      }
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      validateAtomIndex(frame3_atoms[i]);
      if (frame4_atoms[i] >= 0) {
        rtErr("Frame type " + getEnumerationName(vs_frame_types[parameter_indices[i]]) +
              " has three frame atoms, but a fourth was specified.", "AtomGraphStage",
              "addVirtualSite");
      }
      break;
    case VirtualSiteKind::FIXED_4:
      validateAtomIndex(frame3_atoms[i]);
      validateAtomIndex(frame4_atoms[i]);
      break;
    }
    
  }
  const std::vector<int> vs_positions = insertParticles(std::vector<int>(nvs, 0), vs_residue_homes,
                                                        parameter_indices, frame1_atoms,
                                                        frame2_atoms, frame3_atoms, frame4_atoms);
  
  return vs_positions;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphStage::addVirtualSite(const int parameter_index, const int frame1_atom,
                                   const int frame2_atom, const int frame3_atom,
                                   const int frame4_atom) {
  const std::vector<int> result = addVirtualSites(std::vector<int>(1, parameter_index),
                                                  std::vector<int>(1, frame1_atom),
                                                  std::vector<int>(1, frame2_atom),
                                                  std::vector<int>(1, frame3_atom),
                                                  std::vector<int>(1, frame4_atom));
  return result[0];
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::setBonds(const std::vector<int> &atom_i, const std::vector<int> &atom_j,
                              const std::vector<int> &parameter_index) {
  if (parameter_index.size() != atom_i.size()) {
    rtErr("Parameter indices (" + std::to_string(parameter_index.size()) + ") must be provided "
          "for each atom (" + std::to_string(atom_i.size()) + ").", "AtomGraphStage", "setBonds");
  }
  const std::vector<int> cxns = addConnections(atom_i, atom_j);
  const int ilim = atom_i.size();
  for (int i = 0; i < ilim; i++) {
    bond_parameter_indices[cxns[i]] = parameter_index[i];
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::setBond(const int atom_i, const int atom_j, const int parameter_index) {
  const int cxns = addConnection(atom_i, atom_j);
  bond_parameter_indices[cxns] = parameter_index;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::setElectrostatic14Screening(const double screening_factor_in) {
  elec14_screening_factor = screening_factor_in;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::setVanDerWaals14Screening(const double screening_factor_in) {
  vdw14_screening_factor = screening_factor_in;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::validateAtomIndex(const int atom_index) const {
  if (atom_index < 0 || atom_index >= atom_count) {
    rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a collection of " +
          std::to_string(atom_count) + " atoms.", "AtomGraphStage", "validateAtomIndex");
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::validateFrameIndex(const int set_index) const {
  if (set_index < 0 || set_index >= virtual_site_parameter_set_count) {
    rtErr("Atom index " + std::to_string(set_index) + " is invalid for a collection of " +
          std::to_string(virtual_site_parameter_set_count) + " unique virtual site frames.",
          "AtomGraphStage", "validateFrameIndex");
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraphStage::insertParticles(const std::vector<int> &z_numbers,
                                                 const std::vector<int> &residue_homes,
                                                 const std::vector<int> &added_vs_parameter_idx,
                                                 const std::vector<int> &added_vs_frame1_atoms,
                                                 const std::vector<int> &added_vs_frame2_atoms,
                                                 const std::vector<int> &added_vs_frame3_atoms,
                                                 const std::vector<int> &added_vs_frame4_atoms) {
  const int nres_homes = residue_homes.size();
  const int nparticles = z_numbers.size();
  if (nparticles == 0) {
    return std::vector<int>();
  }
  
  // Compute the number of atoms in each residue
  int max_res;
  if (nres_homes == nparticles) {
    max_res = std::max(residue_count, maxValue(residue_homes) + 1);
  }
  else if (nres_homes > 0) {
    max_res = std::max(residue_count, maxValue(residue_homes) + 2);
  }
  else {
    max_res = residue_count + 1;
  }
  std::vector<int> residue_atom_counts(max_res, 0);
  for (int i = 0; i < residue_count; i++) {
    residue_atom_counts[i] = residue_limits[i + 1] - residue_limits[i];
  }
  std::vector<int> orig_residue_atom_counts = residue_atom_counts;

  // Determine how the residue limits will evolve.  If no residue indices are provided for some or
  // all of the atoms, place them in a single "residue" at the back of the list.
  for (int i = 0; i < nres_homes; i++) {
    if (i >= nres_homes) {
      residue_atom_counts[max_res - 1] += 1;
      continue;
    }
    if (residue_homes[i] < 0) {
      rtErr("No atom may be assigned to a residue of index less than zero (" +
            std::to_string(residue_homes[i]) + ").", "AtomGraphStage", "insertParticles");
    }
    residue_atom_counts[residue_homes[i]] += 1;
  }
  residue_atom_counts.push_back(0);
  std::vector<int> atoms_added_by_residue(max_res);
  for (int i = 0; i < max_res; i++) {
    atoms_added_by_residue[i] = residue_atom_counts[i] - orig_residue_atom_counts[i];
  }

  // Expand the residues array
  prefixSumInPlace(&residue_atom_counts, PrefixSumType::EXCLUSIVE);

  // Prepare to keep a map of the new atom indices, in order to transform valence term indexing
  // and other information that invovles the atom numbering.  The new residue limits are now
  // contained in residue_atom_counts, which can be compared to the object's residue_limits array
  // to understand how many places forward a particular atom will jump.  
  std::vector<int> revised_atom_indices(atom_count);
  for (int i = 0; i < residue_count; i++) {
    const int res_adv = residue_atom_counts[i] - residue_limits[i];
    for (int j = residue_limits[i]; j < residue_limits[i + 1]; j++) {
      revised_atom_indices[i] = j + res_adv;
    }
  }
  residue_limits = residue_atom_counts;
  residue_count = max_res;
  const int orig_atom_count = atom_count;
  atom_count += nparticles;
  
  // Adjust the arrays of atom descriptors in step with the additions.  Stepping backwards through
  // the arrays works because revised_atom_indices[i] >= i.
  atomic_numbers.resize(atom_count);
  atom_names.resize(atom_count);
  atom_types.resize(atom_count);
  charge_indices.resize(atom_count);
  lennard_jones_indices.resize(atom_count);
  atomic_charges.resize(atom_count);
  for (int i = orig_atom_count; i >= 0; i--) {
    atomic_numbers[revised_atom_indices[i]] = atomic_numbers[i];
    atom_names[revised_atom_indices[i]] = atom_names[i];
    atom_types[revised_atom_indices[i]] = atom_types[i];
    charge_indices[revised_atom_indices[i]] = charge_indices[i];
    lennard_jones_indices[revised_atom_indices[i]] = lennard_jones_indices[i];
    atomic_charges[revised_atom_indices[i]] = atomic_charges[i];
  }

  // Insert particle with the specified atomic numbers into the spaces created amidst the
  // pre-existing atoms.  This is, in effect, placing the new atoms into the system and giving them
  // valid topological indices.  Record where each of the added atoms resides.
  std::vector<int> added_atom_destinations(nparticles);
  int pos = 0;
  for (int i = 0; i < residue_count; i++) {
    const int ioffset = residue_limits[i] + orig_residue_atom_counts[i];
    for (int j = 0; j < atoms_added_by_residue[i]; j++) {
      atomic_numbers[ioffset + j] = z_numbers[pos];
      added_atom_destinations[pos] = ioffset + j;
      pos++;
    }
  }

  // Adjust valence term indexing.  While the atomic_numbers and other atom descriptor arrays are
  // ordered with the topology and therefore must be adjusted in a way that accommodates new
  // additions, the valence indexing lists are merely numbered atoms in the topology and can be
  // adjusted in any order.  This routine is about inserting particles, not connections or bonded
  // terms.  The parameter indices remain the same as no new parameters are being added.
  for (int pos = 0; pos < urey_bradley_term_count; pos++) {
    urey_bradley_i_atoms[pos] = revised_atom_indices[urey_bradley_i_atoms[pos]];
    urey_bradley_k_atoms[pos] = revised_atom_indices[urey_bradley_k_atoms[pos]];
  }
  for (int pos = 0; pos < charmm_impr_term_count; pos++) {
    charmm_impr_i_atoms[pos] = revised_atom_indices[charmm_impr_i_atoms[pos]];
    charmm_impr_j_atoms[pos] = revised_atom_indices[charmm_impr_j_atoms[pos]];
    charmm_impr_k_atoms[pos] = revised_atom_indices[charmm_impr_k_atoms[pos]];
    charmm_impr_l_atoms[pos] = revised_atom_indices[charmm_impr_l_atoms[pos]];
  }
  for (int pos = 0; pos < cmap_term_count; pos++) {
    cmap_i_atoms[pos] = revised_atom_indices[cmap_i_atoms[pos]];
    cmap_j_atoms[pos] = revised_atom_indices[cmap_j_atoms[pos]];
    cmap_k_atoms[pos] = revised_atom_indices[cmap_k_atoms[pos]];
    cmap_l_atoms[pos] = revised_atom_indices[cmap_l_atoms[pos]];
    cmap_m_atoms[pos] = revised_atom_indices[cmap_m_atoms[pos]];
  }
  for (int pos = 0; pos < bond_term_count; pos++) {
    bond_i_atoms[pos] = revised_atom_indices[bond_i_atoms[pos]];
    bond_j_atoms[pos] = revised_atom_indices[bond_j_atoms[pos]];
  }
  for (int pos = 0; pos < angl_term_count; pos++) {
    angl_i_atoms[pos] = revised_atom_indices[angl_i_atoms[pos]];
    angl_j_atoms[pos] = revised_atom_indices[angl_j_atoms[pos]];
    angl_k_atoms[pos] = revised_atom_indices[angl_k_atoms[pos]];
  }
  for (int pos = 0; pos < dihe_term_count; pos++) {
    dihe_i_atoms[pos] = revised_atom_indices[dihe_i_atoms[pos]];
    dihe_j_atoms[pos] = revised_atom_indices[dihe_j_atoms[pos]];
    dihe_k_atoms[pos] = revised_atom_indices[dihe_k_atoms[pos]];
    dihe_l_atoms[pos] = revised_atom_indices[dihe_l_atoms[pos]];
  }
  for (int pos = 0; pos < inferred_14_attenuations; pos++) {
    infr14_i_atoms[pos] = revised_atom_indices[infr14_i_atoms[pos]];
    infr14_l_atoms[pos] = revised_atom_indices[infr14_l_atoms[pos]];
  }

  // Adjust virtual sites, beginning with the indexing of sites that already exist.
  for (int i = 0; i < virtual_site_count; i++) {
    virtual_site_atoms[i]        = revised_atom_indices[virtual_site_atoms[i]];
    virtual_site_frame1_atoms[i] = revised_atom_indices[virtual_site_frame1_atoms[i]];
    virtual_site_frame2_atoms[i] = revised_atom_indices[virtual_site_frame2_atoms[i]];
    virtual_site_frame3_atoms[i] = revised_atom_indices[virtual_site_frame3_atoms[i]];
    virtual_site_frame4_atoms[i] = revised_atom_indices[virtual_site_frame4_atoms[i]];
  }
  
  // Detect whether any of the added atoms are virtual sites and add them to the arrays of virtual
  // sites.  Push new virtual sites to the back of the list and sort it post-hoc to reorder the
  // virtual sites by increasing topological indices.
  int n_added_vs = 0;
  for (int i = 0; i < nparticles; i++) {
    n_added_vs += (z_numbers[i] == 0);
  }
  virtual_site_count += n_added_vs;
  virtual_site_atoms.reserve(virtual_site_count);
  virtual_site_parameter_indices.reserve(virtual_site_count);
  virtual_site_frame1_atoms.reserve(virtual_site_count);
  virtual_site_frame2_atoms.reserve(virtual_site_count);
  virtual_site_frame3_atoms.reserve(virtual_site_count);
  virtual_site_frame4_atoms.reserve(virtual_site_count);
  n_added_vs = 0;
  for (int i = 0; i < nparticles; i++) {
    if (z_numbers[i] != 0) {
      continue;
    }
    const VirtualSiteKind added_vs_type = vs_frame_types[added_vs_parameter_idx[i]];
    if (added_vs_parameter_idx.size() <= i) {
      rtErr("A frame type must be available for added particle index " + std::to_string(i) +
            ", which has atomic number 0 (virtual site).", "AtomGraphStage", "insertParticles");
    }
    if (added_vs_frame1_atoms.size() <= i || added_vs_frame2_atoms.size() <= i) {
      rtErr("At least two frame atoms must be present for a virtual site in added particle "
            "index " + std::to_string(i) + ".", "AtomGraphStage", "insertParticles");
    }
    virtual_site_atoms.push_back(added_atom_destinations[i]);
    virtual_site_parameter_indices.push_back(added_vs_parameter_idx[i]);
    virtual_site_frame1_atoms.push_back(added_vs_frame1_atoms[i]);
    virtual_site_frame2_atoms.push_back(added_vs_frame2_atoms[i]);
    switch (added_vs_type) {
    case VirtualSiteKind::NONE:

      // Error trapped above
      break;
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:

      // All necessary information has been added, but placeholders are still needed to keep the
      // atoms in the proper order in arrays for the three- and four-atom frame types.
      virtual_site_frame3_atoms.push_back(-1);
      virtual_site_frame4_atoms.push_back(-1);
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      if (added_vs_frame3_atoms.size() <= i) {
        rtErr("A third frame atoms must be present for a virtual site in added particle index " +
              std::to_string(i) + ".", "AtomGraphStage", "insertParticles");
      }
      virtual_site_frame3_atoms.push_back(added_vs_frame3_atoms[i]);
      virtual_site_frame4_atoms.push_back(-1);
      break;
    case VirtualSiteKind::FIXED_4:
      if (added_vs_frame3_atoms.size() <= i || added_vs_frame4_atoms.size() <= i) {
        rtErr("A fourth frame atoms must be present for a virtual site in added particle index " +
              std::to_string(i) + ".", "AtomGraphStage", "insertParticles");
      }
      virtual_site_frame3_atoms.push_back(added_vs_frame3_atoms[i]);
      virtual_site_frame4_atoms.push_back(added_vs_frame4_atoms[i]);
      break;
    }

    // Re-order the virtual sites to arrange them in order of increasing topological indices.
    std::vector<int2> vs_order(virtual_site_count);
    for (int i = 0; i < virtual_site_count; i++) {
      vs_order[i] = { virtual_site_atoms[i], i };
    }
    std::sort(vs_order.begin(), vs_order.end(), [](int2 a, int2 b) { return a.x < b.x; });

    // A follow-the-cycles approach might accomplish this in place, but the extra memory is not
    // significant when a single array can be usedfor out-of-place re-arrangement of all integer
    // components.
    virtual_site_atoms             = applyAssociatedSort(virtual_site_atoms, vs_order);
    virtual_site_parameter_indices = applyAssociatedSort(virtual_site_parameter_indices, vs_order);
    virtual_site_frame1_atoms      = applyAssociatedSort(virtual_site_frame1_atoms, vs_order);
    virtual_site_frame2_atoms      = applyAssociatedSort(virtual_site_frame2_atoms, vs_order);
    virtual_site_frame3_atoms      = applyAssociatedSort(virtual_site_frame3_atoms, vs_order);
    virtual_site_frame4_atoms      = applyAssociatedSort(virtual_site_frame4_atoms, vs_order);
  }
  return added_atom_destinations;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraphStage::addConnections(const std::vector<int> &atom_i,
                                                const std::vector<int> &atom_j) {
  if (atom_i.size() != atom_j.size()) {
    rtErr("The number of atoms in each array (" + std::to_string(atom_i.size()) + ", " +
          std::to_string(atom_j.size()) + ") must match.", "AtomGraphStage", "addConnection");
  }
  const int nrequested = atom_i.size();
  int nconnect = nrequested;
  if (nconnect == 0) {
    return std::vector<int>();
  }
  for (int pos = 0; pos < nconnect; pos++) {
    validateAtomIndex(atom_i[pos]);
    validateAtomIndex(atom_j[pos]);
  }
  
  // Prune the input lists of duplicates, including duplicates of existing bonds.  Maintain the
  // ordering of connections such that the ith atom is of lower topological index that the jth
  // atom and the list of new connections goes in ascending order of the ith atoms.
  std::vector<int4> connectors(bond_term_count + nconnect);
  for (int pos = 0; pos < bond_term_count; pos++) {
    connectors[pos] = { bond_i_atoms[pos], bond_j_atoms[pos], pos, bond_parameter_indices[pos] };
  }
  for (int pos = 0; pos < nconnect; pos++) {
    if (atom_i[pos] < atom_j[pos]) {
      connectors[bond_term_count + pos] = { atom_i[pos], atom_j[pos], bond_term_count + pos, -1 };
    }
    else if (atom_i[pos] > atom_j[pos]) {
      connectors[bond_term_count + pos] = { atom_j[pos], atom_i[pos], bond_term_count + pos, -1 };
    }
    else {
      continue;
    }
  }
  nconnect += bond_term_count;

  // The lambda function in this sort will order connections by ascending order of the bond I
  // atoms, then bond J atoms, and finally put new requests for connections already covered ahead
  // of existing connections to give priority to the last statement of any particular bond.  The
  // parameter indices of existing connections are added to the sort, but no parameter indices are
  // assumed for new connections.  Requesting a connection between atoms that are already connected
  // and parameterized will erase the existing parameterization. 
  std::sort(connectors.begin(), connectors.end(),
            [](int4 a, int4 b) { return (a.x < b.x || (a.x == b.x && a.y < b.y) ||
                                         (a.x == b.x && a.y == b.y && a.z > b.z)); });
  std::vector<int> cxn_positions(nrequested, -1);
  int4 last_unique = connectors[0];
  if (connectors[0].z >= bond_term_count) {
    cxn_positions[connectors[0].z - bond_term_count] = 0;
  }
  int nunique = 1;
  for (int pos = 1; pos < nconnect; pos++) {
    if (connectors[pos].x != last_unique.x || connectors[pos].y != last_unique.y) {
      connectors[nunique] = connectors[pos];
      last_unique = connectors[pos];
      if (connectors[pos].z >= bond_term_count) {
        cxn_positions[connectors[pos].z - bond_term_count] = nunique;
      }
      nunique++;
    }
  }
  connectors.resize(nunique);
  bond_term_count = nunique;
  bond_i_atoms.resize(bond_term_count);
  bond_j_atoms.resize(bond_term_count);
  bond_parameter_indices.resize(bond_term_count);
  for (int pos = 0; pos < bond_term_count; pos++) {
    bond_i_atoms[pos] = connectors[pos].x;
    bond_j_atoms[pos] = connectors[pos].y;
    bond_parameter_indices[pos] = connectors[pos].w;
  }

  // Rebuild the list of exclusions.  The list of virtual sites is expected to be up-to-date.
  nb11_exclusion_bounds.resize(0);
  nb12_exclusion_bounds.resize(0);
  nb13_exclusion_bounds.resize(0);
  nb14_exclusion_bounds.resize(0);
  nb11_exclusion_bounds.shrink_to_fit();
  nb12_exclusion_bounds.shrink_to_fit();
  nb13_exclusion_bounds.shrink_to_fit();
  nb14_exclusion_bounds.shrink_to_fit();
  nb11_exclusion_list.resize(0);
  nb12_exclusion_list.resize(0);
  nb13_exclusion_list.resize(0);
  nb14_exclusion_list.resize(0);
  nb11_exclusion_list.shrink_to_fit();
  nb12_exclusion_list.shrink_to_fit();
  nb13_exclusion_list.shrink_to_fit();
  nb14_exclusion_list.shrink_to_fit();
  const Map1234 excl = mapExclusions(atom_count, virtual_site_atoms, virtual_site_frame1_atoms,
                                     bond_i_atoms, bond_j_atoms);

  // Map the molecular structure with the exclusions at hand.
  mapMolecules(atom_count, &molecule_count, excl, &molecule_membership, &molecule_limits,
               &molecule_contents);
  
  // Return the arrays of exclusions to the object.
  nb11_exclusion_bounds = std::move(excl.nb11_excl_bounds);
  nb11_exclusion_list = std::move(excl.nb11_excl_list);
  nb12_exclusion_bounds = std::move(excl.nb12_excl_bounds);
  nb12_exclusion_list = std::move(excl.nb12_excl_list);
  nb13_exclusion_bounds = std::move(excl.nb13_excl_bounds);
  nb13_exclusion_list = std::move(excl.nb13_excl_list);
  nb14_exclusion_bounds = std::move(excl.nb14_excl_bounds);
  nb14_exclusion_list = std::move(excl.nb14_excl_list);  

  
  // Return the list of indices into the bonds array where the requested connections fall.
  // Applying parameters to these bonds after this function is called will override any parameters
  // already set for specific bonds, if the requested bonds are already present in the bond list,
  // even though none of the overlapping requested connections will result in duplicates.
  return cxn_positions;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphStage::addConnection(const int atom_i, const int atom_j) {
  std::vector<int> result = addConnections(std::vector<int>(1, atom_i),
                                           std::vector<int>(1, atom_j));
  return result[0];
}

//-------------------------------------------------------------------------------------------------
bool AtomGraphStage::applyProperty(const int array_size, const int added_count,
                                   const std::string &desc) {
  if (array_size > 0 && array_size != added_count) {
    rtErr("A total of " + std::to_string(array_size) + " " + desc + " were provided for " +
          std::to_string(added_count) + " atoms.", "AtomGraphStage", "applyProperty");
  }
  return (array_size > 0);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool>
AtomGraphStage::maskSoluteAtoms(const int minimum_solute_size,
                                const std::vector<char4> &solute_included_residue_names) const {
  std::vector<bool> result(atom_count, false);
  
  // Mask off molecules that are too big to be considered solvent
  for (int i = 0; i < molecule_count; i++) {
    if (molecule_limits[i + 1] - molecule_limits[i] > minimum_solute_size) {
      for (int j = molecule_limits[i]; j < molecule_limits[i + 1]; j++) {
        result[molecule_contents[j]] = true;
      }
    }
  }
  
  // Mask off additional residues based on an optional list
  const size_t ninc = solute_included_residue_names.size();
  if (ninc > 0) {
    for (int i = 0; i < residue_count; i++) {

      // Determine whether residues are included in the solute
      const size_t rllim = residue_limits[i];
      const size_t rhlim = residue_limits[i + 1];
      int j = 0;
      while (j < ninc && result[rllim] == false) {
        for (size_t j = 0; j < ninc; j++) {
          if (solute_included_residue_names[j] == residue_names[i]) {
            for (int k = rllim; k < rhlim; k++) {
              result[k] = true;
            }
          }
        }
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraphStage::createPrmtopDescriptors(const int largest_residue_atoms,
                                                         const int total_exclusions) const {
  std::vector<int> result(static_cast<int>(TopologyDescriptor::N_VALUES), 0);
  result[static_cast<int>(TopologyDescriptor::ATOM_COUNT)] = atom_count;
  result[static_cast<int>(TopologyDescriptor::ATOM_TYPE_COUNT)] = lj_type_count;
  int nbonh = 0;
  for (int pos = 0; pos < bond_term_count; pos++) {
    nbonh += (atomic_numbers[bond_i_atoms[pos]] == 1 || atomic_numbers[bond_j_atoms[pos]] == 1);
  }
  const int mbona = bond_term_count - nbonh;
  result[static_cast<int>(TopologyDescriptor::BONDS_WITH_HYDROGEN)] = nbonh;
  result[static_cast<int>(TopologyDescriptor::BONDS_WITHOUT_HYDROGEN)] = mbona;
  int ntheth = 0;
  for (int pos = 0; pos < angl_term_count; pos++) {
    ntheth += (atomic_numbers[angl_i_atoms[pos]] == 1 || atomic_numbers[angl_j_atoms[pos]] == 1 ||
               atomic_numbers[angl_k_atoms[pos]] == 1);
  }
  const int mtheta = angl_term_count - ntheth;
  result[static_cast<int>(TopologyDescriptor::ANGLES_WITH_HYDROGEN)] = ntheth;
  result[static_cast<int>(TopologyDescriptor::ANGLES_WITHOUT_HYDROGEN)] = mtheta;
  int nphih = 0;
  for (int pos = 0; pos < angl_term_count; pos++) {
    nphih += (atomic_numbers[dihe_i_atoms[pos]] == 1 || atomic_numbers[dihe_j_atoms[pos]] == 1 ||
              atomic_numbers[dihe_k_atoms[pos]] == 1 || atomic_numbers[dihe_l_atoms[pos]] == 1);
  }
  const int mphia = dihe_term_count - nphih;
  result[static_cast<int>(TopologyDescriptor::DIHEDRALS_WITH_HYDROGEN)] = nphih;
  result[static_cast<int>(TopologyDescriptor::DIHEDRALS_WITHOUT_HYDROGEN)] = mphia;

  // The topology stage does not accommodate Locally Enhanced Sampling (LES).  The total excluded
  // atom count includes unique 1:1, 1:2, 1:3, and 1:4 exclusions of all atoms j by atom i such
  // that j > i.  Atoms that have no forward exclusions are counted as having a single "exclusion"
  // (index -1) which will be listed as 0 in the Amber prmtop format.
  result[static_cast<int>(TopologyDescriptor::TOTAL_EXCLUDED_ATOMS)] = total_exclusions;
  result[static_cast<int>(TopologyDescriptor::RESIDUE_COUNT)] = residue_count;
  result[static_cast<int>(TopologyDescriptor::BOND_TYPE_COUNT)] = bond_parameter_count;
  result[static_cast<int>(TopologyDescriptor::ANGLE_TYPE_COUNT)] = angl_parameter_count;
  result[static_cast<int>(TopologyDescriptor::DIHEDRAL_TYPE_COUNT)] = dihe_parameter_count;

  // Perturbations of the system are unsupported.
  const int pbc_class = static_cast<int>(periodic_box_class);
  result[static_cast<int>(TopologyDescriptor::BOX_TYPE_INDEX)] = pbc_class;
  result[static_cast<int>(TopologyDescriptor::ATOM_COUNT_LARGEST_RESIDUE)] = largest_residue_atoms;
  result[static_cast<int>(TopologyDescriptor::EXTRA_POINT_COUNT)] = virtual_site_count;
  
  // Any descriptors not explicitly set above were initialized to zero.
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double>
AtomGraphStage::computeAtomicMasses(const double hmass_repartition_factor) const {

  // Seed the masses with their natural abundances.
  std::vector<double> result(atom_count);
  for (int i = 0; i < atom_count; i++) {
    result[i] = elemental_masses[atomic_numbers[i]];
  }

  // Repartition masses as requested.
  if (fabs(hmass_repartition_factor - 1.0) > constants::tiny) {
    const double heavy_hydrogen_mass = hmass_repartition_factor * elemental_masses[1];
    const double borrow = (hmass_repartition_factor - 1.0) * elemental_masses[1];
    for (int i = 0; i < atom_count; i++) {
      if (atomic_numbers[i] == 1) {

        // The atom must bind to one and only one heavy atom.  However, if that heavy atom has
        // attached virtual sites, there could be more than one 1:2 exclusion.
        int heavy_atom_id;
        int nheavy_found = 0;
        for (int j = nb12_exclusion_bounds[i]; j < nb12_exclusion_bounds[i + 1]; j++) {
          if (atomic_numbers[nb12_exclusion_list[j]] > 1) {
            heavy_atom_id = nb12_exclusion_list[j];
            nheavy_found++;
          }
        }
        if (nheavy_found != 1) {
          rtErr("Hydrogen mass repartitioning requires that hydrogens be bound to one and only "
                "one parent heavy atom.  " + std::to_string(nheavy_found) + " were found.",
                "AtomGraphStage", "computeAtomicMasses");
        }
        result[i] = heavy_hydrogen_mass;
        result[heavy_atom_id] -= borrow;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphStage::buildNonbonded14Screen(std::vector<double> *attn14_elec_factors,
                                            std::vector<double> *attn14_vdw_factors) const {
  if (attn14_tags.size() == 0) {
    attn14_elec_factors->resize(0);
    attn14_vdw_factors->resize(0);
    return;
  }
  attn14_elec_factors->resize(dihe_parameter_count);
  attn14_vdw_factors->resize(dihe_parameter_count);
  double* attn14_elec_factor_ptr = attn14_elec_factors->data();
  double* attn14_vdw_factor_ptr = attn14_vdw_factors->data();

  // Initialize the screening factors to the general electrostatic and van-der Waals settings.
  for (int i = 0; i < dihe_parameter_count; i++) {
    attn14_elec_factor_ptr[i] = elec14_screening_factor;
    attn14_vdw_factor_ptr[i] = vdw14_screening_factor;
  }

  // Check for tags based on atom types.  Prepare tables to help the application of these tags if
  // necessary.
  const int ntags = attn14_tags.size();
  bool atom_type_tags = false;
  for (int i = 0; i < ntags; i++) {
    atom_type_tags = (atom_type_tags || (attn14_tags[i].dihe_idx < 0));
  }
  if (atom_type_tags) {

    // Check that each dihedral parameter set applies between consistent atom types.
    std::vector<bool> coverage(dihe_term_count, false);
    std::vector<int> dihe_parm_bounds(dihe_parameter_count + 1);
    std::vector<int> dihe_parm_locations(dihe_term_count);
    indexingArray(dihe_parameter_indices, &dihe_parm_locations, &dihe_parm_bounds);
    bool consistent = true;
    for (int pidx = 0; pidx < dihe_parameter_count; pidx++) {
      const int llim = dihe_parm_bounds[pidx];
      const int hlim = dihe_parm_bounds[pidx + 1];
      const char4 at_i = atom_types[dihe_i_atoms[dihe_parm_locations[llim]]];
      const char4 at_j = atom_types[dihe_j_atoms[dihe_parm_locations[llim]]];
      const char4 at_k = atom_types[dihe_k_atoms[dihe_parm_locations[llim]]];
      const char4 at_l = atom_types[dihe_l_atoms[dihe_parm_locations[llim]]];
      for (int pos = llim + 1; pos < hlim; pos++) {
        const char4 cmp_at_i = atom_types[dihe_i_atoms[dihe_parm_locations[pos]]];
        const char4 cmp_at_j = atom_types[dihe_j_atoms[dihe_parm_locations[pos]]];
        const char4 cmp_at_k = atom_types[dihe_k_atoms[dihe_parm_locations[pos]]];
        const char4 cmp_at_l = atom_types[dihe_l_atoms[dihe_parm_locations[pos]]];
        consistent = (consistent &&
                      ((at_i == cmp_at_i && at_j == cmp_at_j &&
                        at_k == cmp_at_k && at_l == cmp_at_l) ||
                       (at_l == cmp_at_i && at_k == cmp_at_j &&
                        at_j == cmp_at_k && at_i == cmp_at_l)));
      }
    }
    if (consistent == false) {
      rtErr("In order to apply 1:4 non-bonded screening factors based on atom types, all "
            "dihedrals parameter sets must apply to a consistent quartet of atom types.",
            "AtomGraphStage", "buildNonbonded14Screen");
    }

    // Apply each tag
    for (int tidx = 0; tidx < ntags; tidx++) {
      const char4 at_i = attn14_tags[tidx].i_at;
      const char4 at_j = attn14_tags[tidx].j_at;
      const char4 at_k = attn14_tags[tidx].k_at;
      const char4 at_l = attn14_tags[tidx].l_at;
      for (int pos = 0; pos < dihe_parameter_count; pos++) {
        const int term_idx = dihe_parm_locations[dihe_parm_bounds[pos]];
        const char4 cmp_at_i = atom_types[dihe_i_atoms[term_idx]];
        const char4 cmp_at_j = atom_types[dihe_j_atoms[term_idx]];
        const char4 cmp_at_k = atom_types[dihe_k_atoms[term_idx]];
        const char4 cmp_at_l = atom_types[dihe_l_atoms[term_idx]];
        if ((at_i == cmp_at_i && at_j == cmp_at_j && at_k == cmp_at_k && at_l == cmp_at_l) ||
            (at_l == cmp_at_i && at_k == cmp_at_j && at_j == cmp_at_k && at_i == cmp_at_l)) {
          attn14_elec_factor_ptr[pos] = attn14_tags[tidx].elec14;
          attn14_vdw_factor_ptr[pos]  = attn14_tags[tidx].vdw14;
        }
      }
    }
  }
  else {
    for (int i = 0; i < ntags; i++) {
      if (attn14_tags[i].dihe_idx >= 0 && attn14_tags[i].dihe_idx < dihe_parameter_count) {
        attn14_elec_factor_ptr[attn14_tags[i].dihe_idx] = attn14_tags[i].elec14;
        attn14_vdw_factor_ptr[attn14_tags[i].dihe_idx]  = attn14_tags[i].vdw14;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> atomWithVirtualSiteChildren(const int atom_index,
                                             const std::vector<int> &nb11_exclusion_list,
                                             const std::vector<int> &nb11_exclusion_bounds) {
  if (atom_index < 0 || atom_index >= nb11_exclusion_bounds.size() - 1) {
    rtErr("Atom index " + std::to_string(atom_index) + " is out of bounds for a collection of " +
          std::to_string(nb11_exclusion_bounds.size() - 1) + " atoms.",
          "atomWithVirtualSiteChildren");
  }
  const int llim = nb11_exclusion_bounds[atom_index];
  const int hlim = nb11_exclusion_bounds[atom_index + 1];
  std::vector<int> result(1 + hlim - llim);
  result[0] = atom_index;
  for (int k = llim; k < hlim; k++) {
    result[1 + k - llim] = nb11_exclusion_list[k];
  }
  return result;
}
 
} // namespace topology
} // namespace stormm
