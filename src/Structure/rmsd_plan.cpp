#include "copyright.h"
#include "Chemistry/chemical_features.h"
#include "Chemistry/chemistry_enumerators.h"
#include "Constants/hpc_bounds.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/summation.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_abstracts.h"
#include "rmsd_plan.h"

namespace stormm {
namespace structure {

using card::HybridKind;
using chemistry::ChemicalDetailsKit;
using chemistry::ChemicalFeatures;
using diskutil::getBaseName;
using stmath::roundUp;
using synthesis::SynthesisMapReader;
  
//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const RMSDMethod strategy_in, const double rmf_in,
		   const PhaseSpaceSynthesis *poly_ps_in, const SystemCache *sc_in,
		   const SynthesisCacheMap *scmap_in) :
    plan_count{0}, general_strategy{strategy_in},
    required_mass_fraction{rmf_in},
    masses{HybridKind::ARRAY, "rplan_masses"},
    sp_masses{HybridKind::ARRAY, "rplan_masses_sp"},
    atom_counts{HybridKind::ARRAY, "rplan_atom_counts"},
    atom_starts{HybridKind::ARRAY, "rplan_atom_starts"},
    asymmetric_core_atoms{HybridKind::ARRAY, "rplan_core_atoms"},
    asymmetric_core_counts{HybridKind::ARRAY, "rplan_core_counts"},
    asymmetric_core_starts{HybridKind::ARRAY, "rplan_core_starts"},
    symmetry_group_atoms{HybridKind::ARRAY, "rplan_symm_atoms"},
    symmetry_group_bounds{HybridKind::ARRAY, "rplan_symm_bounds"},
    alignment_steps{HybridKind::ARRAY, "rplan_align_steps"},
    symmetry_group_ranges{HybridKind::ARRAY, "rplan_symm_ranges"},
    alignment_prerequisites{HybridKind::ARRAY, "rplan_prereqs"},
    alignment_protocols{HybridKind::ARRAY, "rplan_protocols"},
    symmetry_depth{HybridKind::ARRAY, "rplan_depths"},
    symmetry_group_membership{HybridKind::ARRAY, "rplan_membership"},
    symmetry_group_membership_bounds{HybridKind::ARRAY, "rplan_mmb_bounds"},
    ag_pointers{},
    poly_ps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps_in)},
    sc_ptr{const_cast<SystemCache*>(sc_in)},
    scmap_ptr{const_cast<SynthesisCacheMap*>(scmap_in)}
{}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                   const RMSDMethod strategy_in, const double rmf_in, const GpuDetails &gpu,
                   const int low_mol_idx, const int high_mol_idx) :
    RMSDPlan(strategy_in, rmf_in)
{
  // There is one plan that can be defined by the given information.
  plan_count = 1;
  ag_pointers.resize(1, const_cast<AtomGraph*>(ag_in.getSelfPointer()));

  // Get the formal charges, free electron content, chirality, and ring inclusions for the one
  // system.  Each system can, in fact, hold many molecules, but the constructor is set to 
  ChemicalFeatures chemfe(ag_in, cf_in);
  std::vector<AtomEquivalence> eq_found;
  eq_found.reserve(1);
  eq_found.emplace_back(chemfe, nullptr, low_mol_idx, high_mol_idx);
  chartSymmetryGroups(eq_found);
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const AtomEquivalence &eq_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu) :
    RMSDPlan(strategy_in, rmf_in)
{
  // There is one plan defined by the given information
  plan_count = 1;
  ag_pointers.resize(1, const_cast<AtomGraph*>(eq_in.getTopologyPointer()));
  std::vector<AtomEquivalence> eq_tables(1, eq_in);
  chartSymmetryGroups(eq_tables);
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis *poly_ps_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu, const int low_mol_idx,
                   const int high_mol_idx) :
    RMSDPlan(strategy_in, rmf_in, poly_ps_in)
{
  // There will be as many plans as the PhaseSpaceSynthesis has unique topologies.  It will be
  // assumed that a common range of molecules (again, defaulting to the first molecule) applies to
  // the desired RMSD calculations across all systems.
  ag_pointers = poly_ps_in->getUniqueTopologies();
  plan_count = ag_pointers.size();
  const std::vector<int>& sys_top_indices = poly_ps_in->getUniqueTopologyExampleIndices();
  std::vector<AtomEquivalence> eq_tables;
  eq_tables.reserve(plan_count);
  for (int i = 0; i < plan_count; i++) {
    const CoordinateFrame cf_example = poly_ps_in->exportCoordinates(sys_top_indices[i]);
    ChemicalFeatures chemfe(ag_pointers[i], cf_example);
    eq_tables.emplace_back(chemfe, nullptr, low_mol_idx, high_mol_idx);
  }
  chartSymmetryGroups(eq_tables);
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis *poly_ps_in,
                   const std::vector<AtomEquivalence> &eq_list_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu) :
    RMSDPlan(strategy_in, rmf_in, poly_ps_in)
{
  // Check that there are enough plans to cover the various topologies in the synthesis, and the
  // topologies referenced by each plan should match those in the synthesis, at least in terms of
  // the atom counts.
  ag_pointers = poly_ps_in->getUniqueTopologies();
  if (ag_pointers.size() != eq_list_in.size()) {
    rtErr("One atom equivalence object must be provided for each unique topology in the "
          "synthesis.  " + std::to_string(ag_pointers.size()) + " unique topologies were "
          "detected, with " + std::to_string(eq_list_in.size()) + " atom equivalence table to "
          "support them.", "RMSDPlan");
  }
  plan_count = ag_pointers.size();
  for (int i = 0; i < plan_count; i++) {
    if (ag_pointers[i] != eq_list_in[i].getTopologyPointer() &&
        ag_pointers[i]->getAtomCount() != eq_list_in[i].getTopologyPointer()->getAtomCount()) {
      rtErr("The atom count for the unique topology index " + std::to_string(i) + " within the "
            "synthesis, (" + std::to_string(ag_pointers[i]->getAtomCount()) + " originating in "
            "file " + getBaseName(ag_pointers[i]->getFileName()) + "), does not match "
            "the atom count in the atom equivalence table (" +
            std::to_string(eq_list_in[i].getTopologyPointer()->getAtomCount()) + ", originating "
            "in file " + getBaseName(eq_list_in[i].getTopologyPointer()->getFileName()) + ").",
            "RMSDPlan");
    }
  }
  chartSymmetryGroups(eq_list_in);
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu, const int low_mol_idx,
                   const int high_mol_idx) :
    RMSDPlan(poly_ps_in.getSelfPointer(), strategy_in, rmf_in, gpu, low_mol_idx, high_mol_idx)
{}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in,
                   const std::vector<AtomEquivalence> &eq_list_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu) :
    RMSDPlan(poly_ps_in.getSelfPointer(), eq_list_in, strategy_in, rmf_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis *poly_ps_in, const SystemCache *sc_in,
                   const SynthesisCacheMap *scmap_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu, const int low_mol_idx,
                   const int high_mol_idx) :
    RMSDPlan(strategy_in, rmf_in, poly_ps_in, sc_in, scmap_in)
{
  // Order the topologies against the systems cache, even if there are extra topologies in it.
  ag_pointers = sc_in->getTopologyPointer();
  plan_count = ag_pointers.size();
  std::vector<AtomEquivalence> eq_tables;
  eq_tables.reserve(plan_count);
  for (int i = 0; i < plan_count; i++) {
    const int iexample = sc_in->getSystemExampleIndex(i);
    eq_tables.emplace_back(sc_in->getFeatures(iexample), nullptr, low_mol_idx, high_mol_idx);
  }
  chartSymmetryGroups(eq_tables);
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, const SystemCache &sc_in,
                   const SynthesisCacheMap &scmap_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu, const int low_mol_idx,
                   const int high_mol_idx) :
  RMSDPlan(poly_ps_in.getSelfPointer(), sc_in.getSelfPointer(), scmap_in.getSelfPointer(),
           strategy_in, rmf_in, gpu, low_mol_idx, high_mol_idx)
{}

//-------------------------------------------------------------------------------------------------
int RMSDPlan::getPlanCount() const {
  return plan_count;
}

//-------------------------------------------------------------------------------------------------
int RMSDPlan::getPlanIndex(const int index, const SystemGrouping organization) const {
  if (scmap_ptr != nullptr) {
    const SynthesisMapReader scmapr = scmap_ptr->data();
    int synthesis_index;
    switch (organization) {
    case SystemGrouping::SOURCE:
      synthesis_index = scmapr.csystem_proj[scmapr.csystem_bounds[index]];
      break;
    case SystemGrouping::TOPOLOGY:
      synthesis_index = scmapr.ctopol_proj[scmapr.ctopol_bounds[index]];
      break;
    case SystemGrouping::LABEL:
      synthesis_index = scmapr.clabel_proj[scmapr.clabel_bounds[index]];
      break;
    }
    return scmapr.topology_origins[synthesis_index];
  }
  else if (poly_ps_ptr != nullptr) {
    return poly_ps_ptr->getUniqueTopologyIndex(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
RMSDMethod RMSDPlan::getGeneralStrategy() const {
  return general_strategy;
}

//-------------------------------------------------------------------------------------------------
double RMSDPlan::getRequiredMassFraction() const {
  return required_mass_fraction;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* RMSDPlan::getTopologyPointer(const int plan_index) const {
  return ag_pointers[plan_index];
}

//-------------------------------------------------------------------------------------------------
RMSDAlignmentProtocol RMSDPlan::getAlignmentProtocol(const int plan_index) const {
  return static_cast<RMSDAlignmentProtocol>(alignment_protocols.readHost(plan_index));
}

//-------------------------------------------------------------------------------------------------
const RMSDPlanReader<double> RMSDPlan::dpData(const HybridTargetLevel tier) const {
  return RMSDPlanReader<double>(plan_count, general_strategy, required_mass_fraction,
                                masses.data(tier), atom_counts.data(tier), atom_starts.data(tier),
                                alignment_steps.data(tier), asymmetric_core_atoms.data(tier),
                                asymmetric_core_counts.data(tier),
                                asymmetric_core_starts.data(tier),
                                symmetry_group_atoms.data(tier), symmetry_group_bounds.data(tier),
                                symmetry_group_ranges.data(tier));
}

//-------------------------------------------------------------------------------------------------
const RMSDPlanReader<float> RMSDPlan::spData(const HybridTargetLevel tier) const {
  return RMSDPlanReader<float>(plan_count, general_strategy, required_mass_fraction,
                               sp_masses.data(tier), atom_counts.data(tier),
                               atom_starts.data(tier), alignment_steps.data(tier),
                               asymmetric_core_atoms.data(tier), asymmetric_core_counts.data(tier),
                               asymmetric_core_starts.data(tier), symmetry_group_atoms.data(tier),
                               symmetry_group_bounds.data(tier), symmetry_group_ranges.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* RMSDPlan::getCoordinateSynthesisPointer() const {
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
const SystemCache* RMSDPlan::getCachePointer() const {
  return sc_ptr;
}

//-------------------------------------------------------------------------------------------------
const SynthesisCacheMap* RMSDPlan::getSynthesisMapPointer() const {
  return scmap_ptr;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void RMSDPlan::upload() {
  masses.upload();
  sp_masses.upload();
  atom_counts.upload();
  atom_starts.upload();
  asymmetric_core_atoms.upload();
  asymmetric_core_counts.upload();
  asymmetric_core_starts.upload();
  symmetry_group_atoms.upload();
  symmetry_group_bounds.upload();
  alignment_steps.upload();
  symmetry_group_ranges.upload();
  alignment_prerequisites.upload();
  alignment_protocols.upload();
  symmetry_depth.upload();
  symmetry_group_membership.upload();
  symmetry_group_membership_bounds.upload();
}

//-------------------------------------------------------------------------------------------------
void RMSDPlan::download() {
  masses.download();
  sp_masses.download();
  atom_counts.download();
  atom_starts.download();
  asymmetric_core_atoms.download();
  asymmetric_core_counts.download();
  asymmetric_core_starts.download();
  symmetry_group_atoms.download();
  symmetry_group_bounds.download();
  alignment_steps.download();
  symmetry_group_ranges.download();
  alignment_prerequisites.download();
  alignment_protocols.download();
  symmetry_depth.download();
  symmetry_group_membership.download();
  symmetry_group_membership_bounds.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void RMSDPlan::addSystemCache(const SystemCache *sc, const SynthesisCacheMap *scmap) {
  sc_ptr = const_cast<SystemCache*>(sc);
  scmap_ptr = const_cast<SynthesisCacheMap*>(scmap);
}

//-------------------------------------------------------------------------------------------------
void RMSDPlan::addSystemCache(const SystemCache &sc, const SynthesisCacheMap &scmap) {
  addSystemCache(sc.getSelfPointer(), scmap.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void RMSDPlan::chartSymmetryGroups(const std::vector<AtomEquivalence> &eq_tables) {

  // Populate each Hybrid object related to the overall atom content of each system
  const int ntab  = eq_tables.size();
  atom_counts.resize(ntab);
  atom_starts.resize(ntab);
  int* atom_counts_ptr = atom_counts.data();
  int* atom_starts_ptr = atom_starts.data();
  int total_atoms = 0;
  for (int i = 0; i < ntab; i++) {
    const int natom = eq_tables[i].getTopologyPointer()->getAtomCount();
    atom_counts_ptr[i] = natom;
    atom_starts_ptr[i] = total_atoms;
    total_atoms += roundUp(natom, warp_size_int);
  }
  masses.resize(total_atoms);
  sp_masses.resize(total_atoms);
  double* mass_ptr = masses.data();
  float* sp_mass_ptr = sp_masses.data();
  for (int i = 0; i < ntab; i++) {
    const ChemicalDetailsKit cdk = eq_tables[i].getTopologyPointer()->getChemicalDetailsKit();
    const int i_offset = atom_starts_ptr[i];
    for (int j = 0; j < cdk.natom; j++) {
      mass_ptr[i_offset + j]    = cdk.masses[j];
      sp_mass_ptr[i_offset + j] = cdk.sp_masses[j];
    }
  }
  
  // Populate each Hybrid object related to the asymmetric cores of each system.
  asymmetric_core_counts.resize(ntab);
  asymmetric_core_starts.resize(ntab);
  std::vector<int> tmp_asymmetric_core_starts(ntab);
  std::vector<int> tmp_asymmetric_core_counts(ntab);
  int asym_start_counter = 0;
  for (int i = 0; i < ntab; i++) {
    const int ncore = eq_tables[i].getAsymmetricAtomCount();
    tmp_asymmetric_core_counts[i] = ncore;
    tmp_asymmetric_core_starts[i] = asym_start_counter;
    asym_start_counter += roundUp(ncore, warp_size_int);
  }
  asymmetric_core_counts.putHost(tmp_asymmetric_core_counts);
  asymmetric_core_starts.putHost(tmp_asymmetric_core_starts);
  asymmetric_core_atoms.resize(asym_start_counter);
  for (int i = 0; i < ntab; i++) {
    asymmetric_core_atoms.putHost(eq_tables[i].getAsymmetricAtoms(), tmp_asymmetric_core_starts[i],
                                  tmp_asymmetric_core_counts[i]);
  }

  // Determine the order of symmetry groups to list in each system.  As above, construct temporary
  // vectors in CPU memory (for speed as well as convenience), then push the results to the related
  // Hybrid objects.  The alignment_steps array gets sized and populated in step with
  // symmetry_group_bounds, like a fifth property for each group's tuple (detailing the way in
  // which its various combinations of subunits can be swapped).
  size_t symmetry_group_atom_count = 0;
  size_t symmetry_group_count = 0;
  int membership_bounds_size = 0;
  for (int i = 0; i < ntab; i++) {
    const int ni_groups = eq_tables[i].getGroupCount();
    symmetry_group_count += roundUp(ni_groups, warp_size_int);
    for (int j = 0; j < ni_groups; j++) {
      symmetry_group_atom_count += static_cast<size_t>(eq_tables[i].getGroupSize(j) *
                                                       eq_tables[i].getGroupOrder(j));
    }
    symmetry_group_atom_count = roundUp(symmetry_group_atom_count, warp_size_zu);
    membership_bounds_size += roundUp(eq_tables[i].getTopologyPointer()->getAtomCount(),
                                      warp_size_int);
  }
  std::vector<int> tmp_symmetry_group_atoms, tmp_alignment_steps;
  std::vector<int4> tmp_symmetry_group_bounds;
  std::vector<int2> tmp_symmetry_group_ranges(ntab);
  std::vector<int> tmp_alignment_prerequisites(ntab, 0);
  std::vector<int> tmp_alignment_protocols(ntab), tmp_symmetry_depth(ntab);
  std::vector<int2> tmp_sg_membership(symmetry_group_atom_count);
  std::vector<ullint> tmp_sg_membership_bounds(membership_bounds_size);
  tmp_symmetry_group_atoms.reserve(symmetry_group_atom_count);
  tmp_symmetry_group_bounds.reserve(symmetry_group_count);
  tmp_alignment_steps.reserve(symmetry_group_count);
  symmetry_group_atoms.resize(symmetry_group_atom_count);
  symmetry_group_bounds.resize(symmetry_group_count);
  symmetry_group_ranges.resize(ntab);
  symmetry_depth.resize(ntab);
  symmetry_group_membership.resize(symmetry_group_atom_count, { -1, -1});
  symmetry_group_membership_bounds.resize(membership_bounds_size);
  alignment_steps.resize(symmetry_group_count);
  alignment_prerequisites.resize(ntab);
  alignment_protocols.resize(ntab);
  int membership_bounds_system_offset = 0;
  ullint membership_bounds_counter = 0LLU;
  for (int i = 0; i < ntab; i++) {

    // Track which groups were included as part of the work to establish the core's critical mass
    const int ni_groups = eq_tables[i].getGroupCount();
    std::vector<bool> groups_included(ni_groups, false);
    
    // Mark the lower limit of this system's symmetry group bounds
    tmp_symmetry_group_ranges[i].x = tmp_symmetry_group_bounds.size();

    // Find the largest symmetry groups in the system which are not themselves dependents.
    int ni_head_groups = 0;
    std::vector<bool> is_head_group(ni_groups, true);
    for (int j = 0; j < ni_groups; j++) {
      if (is_head_group[j]) {
        std::vector<int> jdeps = eq_tables[i].getGroupDependencies(j);
        const int n_jdeps = jdeps.size();
        for (int k = 0; k < n_jdeps; k++) {
          is_head_group[jdeps[k]] = false;
        }
        ni_head_groups++;        
      }
    }
    std::vector<int2> head_groups(ni_head_groups);
    ni_head_groups = 0;
    for (int j = 0; j < ni_groups; j++) {
      if (is_head_group[j]) {
        head_groups[ni_head_groups].x = j;
        head_groups[ni_head_groups].y = eq_tables[i].getGroupSize(j);
        ni_head_groups++;
      }
    }

    // Determine the number of head groups that is needed to account for the required mass
    // fraction.
    const ChemicalDetailsKit cdk = eq_tables[i].getTopologyPointer()->getChemicalDetailsKit();
    const std::vector<int>& asym_atoms = eq_tables[i].getAsymmetricAtoms();
    const int n_asym_atoms = eq_tables[i].getAsymmetricAtomCount();
    const double total_mass = eq_tables[i].getTopologyPointer()->getTotalMass();
    double core_mass = 0.0;
    for (int j = 0; j < n_asym_atoms; j++) {
      core_mass += cdk.masses[asym_atoms[j]];
    }
    bool core_sgrps_ok = true;
    while (core_mass < required_mass_fraction * total_mass) {
      for (int j = 0; j < ni_head_groups; j++) {

        // Add the group itself and mark that it has been included
        const int hgj = head_groups[j].x;
        const std::vector<int> group_atoms = eq_tables[i].getGroup(hgj);
        const int n_group_atoms = group_atoms.size();
        int4 bounds_ph;
        bounds_ph.x = tmp_symmetry_group_atoms.size();
        for (int k = 0; k < n_group_atoms; k++) {
          const int katom = group_atoms[k];
          core_mass += cdk.masses[katom];
          tmp_symmetry_group_atoms.push_back(katom);
        }
        bounds_ph.y = tmp_symmetry_group_atoms.size();
        bounds_ph.z = eq_tables[i].getGroupSize(hgj);
        bounds_ph.w = eq_tables[i].getGroupOrder(hgj);
        core_sgrps_ok = (core_sgrps_ok && bounds_ph.w < 8);
        tmp_symmetry_group_bounds.push_back(bounds_ph);
        tmp_alignment_steps.push_back(static_cast<int>(eq_tables[i].getGroupRule(hgj)));
        groups_included[hgj] = true;

        // Add the group's dependencies and mark that they have been included
        const std::vector<int> jdeps = eq_tables[i].getGroupDependencies(hgj);
        const int n_jdeps = jdeps.size();
        for (int k = 0; k < n_jdeps; k++) {
          const std::vector<int> dep_atoms = eq_tables[i].getGroup(jdeps[k]);
          const int n_dep_atoms = dep_atoms.size();
          bounds_ph.x = tmp_symmetry_group_atoms.size();
          for (int m = 0; m < n_dep_atoms; m++) {
            const int matom = dep_atoms[m];
            tmp_symmetry_group_atoms.push_back(matom);
          }
          bounds_ph.y = tmp_symmetry_group_atoms.size();
          bounds_ph.z = eq_tables[i].getGroupSize(jdeps[k]);
          bounds_ph.w = eq_tables[i].getGroupOrder(jdeps[k]);
          core_sgrps_ok = (core_sgrps_ok && bounds_ph.w < 8);
          tmp_symmetry_group_bounds.push_back(bounds_ph);
          tmp_alignment_steps.push_back(static_cast<int>(eq_tables[i].getGroupRule(jdeps[k])));
          groups_included[jdeps[k]] = true;
        }
        tmp_alignment_prerequisites[i] += 1 + n_jdeps;
      }
    }

    // It is essential to remember, for each alignment attempted, the settings of each
    // symmetry-related group comprising the core.  If this information cannot be stored locally
    // by the function or thread computing the alignment, it is extremely difficult to allocate
    // and manage a space to serve all of the different RMSD computation protocols that one might
    // attempt.  Plus, because the core must be aligned with a combinatorial enumeration of all
    // groups involved, it becomes prohibitively expensive to do that alignment if there are a
    // large number of groups with large numbers of symmetries.  Therefore, limit the number of
    // core groups to sixteen and ensure that each of them has at most sixteen-fold symmetry.  If
    // either of these conditions is violated, short-circuit the symmetry-aware RMSD calculation
    // and do a straight "ALIGN_ALL" calculation.  Otherwise, specify "BUILD_CORE" if the core
    // includes symmetry-related groups that must be determined, "ALIGN_CORE" if the asymmetric
    // atoms are themselves of sufficient mass to determine the core, and "ALIGN_ALL" if there are
    // no symmetry-related groups to consider.
    if (tmp_alignment_prerequisites[i] > 0) {
      if (tmp_alignment_prerequisites[i] < 16 && core_sgrps_ok) {
        tmp_alignment_protocols[i] = static_cast<int>(RMSDAlignmentProtocol::BUILD_CORE);
      }
      else {
        tmp_alignment_protocols[i] = static_cast<int>(RMSDAlignmentProtocol::ALIGN_ALL);
      }
    }
    else {
      if (ni_groups > 0) {
        tmp_alignment_protocols[i] = static_cast<int>(RMSDAlignmentProtocol::ALIGN_CORE);
      }
      else {
        tmp_alignment_protocols[i] = static_cast<int>(RMSDAlignmentProtocol::ALIGN_ALL);
      }
    }

    // Record other groups
    for (int j = 0; j < ni_groups; j++) {
      if (groups_included[j]) {
        continue;
      }
      const std::vector<int> group_atoms = eq_tables[i].getGroup(j);
      const int n_group_atoms = group_atoms.size();
      int4 bounds_ph;
      bounds_ph.x = tmp_symmetry_group_atoms.size();
      for (int k = 0; k < n_group_atoms; k++) {
        const int katom = group_atoms[k];
        tmp_symmetry_group_atoms.push_back(katom);
      }
      bounds_ph.y = tmp_symmetry_group_atoms.size();
      bounds_ph.z = eq_tables[i].getGroupSize(j);
      bounds_ph.w = eq_tables[i].getGroupOrder(j);
      tmp_symmetry_group_bounds.push_back(bounds_ph);
      tmp_alignment_steps.push_back(static_cast<int>(eq_tables[i].getGroupRule(j)));
      groups_included[j] = true;
    }
    tmp_symmetry_depth[i] = eq_tables[i].getSymmetryDepth();
    
    // Mark the symmetry group membership of all atoms at all levels.  Begin by assembling the
    // membership bounds array for this system
    std::vector<int> iatom_levels(cdk.natom, 0);
    for (int j = 0; j < ni_groups; j++) {
      const int jgroup_order = eq_tables[i].getGroupOrder(j);
      const int jgroup_size  = eq_tables[i].getGroupSize(j);
      for (int k = 0; k < jgroup_order; k++) {
        for (int m = 0; m < jgroup_size; m++) {
          iatom_levels[eq_tables[i].getSymmetryRelatedAtom(j, k, m)] += 1;
        }
      }
    }
    for (int j = 0; j < cdk.natom; j++) {
      const ullint iatom_nlev = iatom_levels[j];
      tmp_sg_membership_bounds[membership_bounds_system_offset + j] = (membership_bounds_counter |
                                                                       (iatom_nlev << 48));
      membership_bounds_counter += iatom_nlev;
    }
    for (int j = 0; j < ni_groups; j++) {
      const int jgroup_order    = eq_tables[i].getGroupOrder(j);
      const int jgroup_size     = eq_tables[i].getGroupSize(j);
      const size_t jgroup_level = eq_tables[i].getGroupLevel(j);
      for (int k = 0; k < jgroup_order; k++) {
        for (int m = 0; m < jgroup_size; m++) {

          // Determine what atom this is, then parse the appropriate bound, and insert information
          // about this group and domain at the appropriate space in the membership array.
          const int jkm_topological_index = eq_tables[i].getSymmetryRelatedAtom(j, k, m);
          const int jkm_synthesis_index = membership_bounds_system_offset + jkm_topological_index;
          const ullint compound_bound_index = tmp_sg_membership_bounds[jkm_synthesis_index];
          const size_t membership_index =
            static_cast<size_t>(compound_bound_index & 0xffffffffffffULL) + jgroup_level;
          tmp_sg_membership[membership_index].x = tmp_symmetry_group_ranges[i].x + j;
          tmp_sg_membership[membership_index].y = k;
        }
      }
    }
    membership_bounds_system_offset += roundUp(cdk.natom, warp_size_int);

    // Mark the upper limit of this system's symmetry group bounds
    tmp_symmetry_group_ranges[i].y = tmp_symmetry_group_bounds.size();
    
    // Use the push_back method to pad the atoms and bounds arrays to ensure that resizing does
    // not inadvertently shrink them.
    size_t jlim = roundUp(tmp_symmetry_group_atoms.size(), warp_size_zu);
    for (size_t j = tmp_symmetry_group_atoms.size(); j < jlim; j++) {
      tmp_symmetry_group_atoms.push_back(-1);
    }
    jlim = roundUp(tmp_symmetry_group_bounds.size(), warp_size_zu);
    for (size_t j = tmp_symmetry_group_bounds.size(); j < jlim; j++) {
      tmp_symmetry_group_bounds.push_back({ -1, -1, -1, -1 });
      tmp_alignment_steps.push_back(-1);
    }
  }
  symmetry_group_atoms.putHost(tmp_symmetry_group_atoms);
  symmetry_group_bounds.putHost(tmp_symmetry_group_bounds);
  symmetry_group_ranges.putHost(tmp_symmetry_group_ranges);
  alignment_steps.putHost(tmp_alignment_steps);
  alignment_prerequisites.putHost(tmp_alignment_prerequisites);
  alignment_protocols.putHost(tmp_alignment_protocols);
  symmetry_depth.putHost(tmp_symmetry_depth);
  symmetry_group_membership.putHost(tmp_sg_membership);
  symmetry_group_membership_bounds.putHost(tmp_sg_membership_bounds);
}

//-------------------------------------------------------------------------------------------------
void RMSDPlan::writeSymmetryGroupCodes() {

  // Each atom may or may not belong to various symmetry groups, which will imply combinatorial
  // tests of the different symmetry-related partner atoms at various stages of computing the
  // best-fit RMSD.  The key is to make this information as compact as possible, and the way to
  // do that is to create bit-packed strings to encode the information about which symmetry groups
  // the atom is a part of.  Let the input coordinates be given { xyz0, xyz1, xyz2, ..., xyzN }
  // for N atoms and three Cartesian dimension x, y, and z.  View each atom in the molecule as a
  // virtual particle (not a virtual site in the sense of a massless particle, but a conceptual
  // entity) that could take the guise of one of a number of different coordinates from within the
  // given representation.  Atoms of topological indices 0, 1, and 2 might be three methyl
  // hydrogens constituting a symmetry group, with atom 3 being the carbon they are bound to.  This
  // methyl group itself may be symmetric with a second methyl group, with hydrogen atoms at
  // topological indices 7, 8, and 9 and the carbon atom at topological index 6.  First, the
  // decision must be made whether to swap the methyl groups, and then the decision can be made
  // in what order to take the hydrogen atoms.  In the final product, the virtual particles take
  // positions from the given coordinates that are then submitted to a standard RMSD calculation.
  // Virtual particle 3 could take coordinates from atoms 3 or 6.  Virtual particles 0, 1, and 2,
  // as well as 7, 8, and 9, could take coordinates from atoms 0, 1, and 2 in any order, without
  // replacement, or form atoms 7, 8, and 9, in any order without replacement.  The dependencies
  // of symmetry groups imply an order in which virtual particles take the coordinates from
  // particular atoms (processing symmetry groups which are not dependent on any other), then take
  // subsequent positions of other virtual particles (for dependent symmetry groups down the line
  // until there are no dependencies left to process).  In the present example, the groups of atoms
  // { 0, 1, 2, 3 } and { 7, 8, 9, 6 } take positions from atoms { 0, 1, 2, 3 } or { 7, 8, 9, 6 }
  // in the first step.  Next, virtual particles { 0, 1, 2 } choose whether to swap positions
  // among themselves, and virtual particles { 7, 8, 9 } go through the same process.  As a first
  // pass, the procedure is to find the best assignment among the highest level symmetry groups
  // (assigning coordinates from atoms in the given list to virtual particles), then filter down,
  // finding further improvements by attempting different assignments of the dependent groups.

}

} // namespace structure
} // namespace stormm
