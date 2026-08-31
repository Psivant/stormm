#include <algorithm>
#include "copyright.h"
#include "Chemistry/atommask.h"
#include "Chemistry/chemical_features.h"
#include "DataTypes/mixed_types.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/statistics.h"
#include "Math/statistical_enumerators.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/ordered_list.h"
#include "Reporting/report_table.h"
#include "Reporting/reporting_enumerators.h"
#include "Reporting/summary_file.h"
#include "Synthesis/systemcache.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "hydrogen_bond_analysis.h"

namespace stormm {
namespace analysis {

using card::HybridFormat;
using card::HybridKind;
using card::HybridLabel;
using chemistry::AtomMask;
using chemistry::ChemicalFeatures;
using diskutil::extractCommonPaths;
using diskutil::listCommonPaths;
using diskutil::splitPath;
using namelist::default_hbond_da_max_separation;
using namelist::default_hbond_dha_min_angle;
using namelist::default_hbond_occ_threshold;
using namelist::default_hbond_act_threshold;
using parse::addTailingWhiteSpace;
using parse::char4ToString;
using parse::findStringInVector;
using parse::intToString;
using parse::JustifyText;
using parse::NumberFormat;
using parse::realToString;
using review::commentSymbol;
using review::ListEnumeration;
using review::OrderedList;
using review::OutputSyntax;
using review::ReportTable;
using stmath::mean;
using stmath::sum;
using stmath::readBitFromMask;
using stmath::roundUp;
using stmath::runningVariance;
using stmath::variance;
using stmath::VarianceMethod;
using synthesis::SystemCache;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using trajectory::CoordinateFrame;

//-------------------------------------------------------------------------------------------------
HBondWriter::HBondWriter(const int n_partners_in, const int stat_blocks_in, const int init_step_in,
                         const int eval_intv_in, const int block_steps_in,
                         const float max_separation_in, const float min_angle_in,
                         const int4* partners_in, const int* partner_bounds_in, int* tallies_in,
                         double* dist_acc_in, double* dist_sqacc_in, double* all_dist_acc_in,
                         double* all_dist_sqacc_in, double* angl_acc_in, double* angl_sqacc_in) :
    n_partners{n_partners_in}, stat_blocks{stat_blocks_in}, init_step{init_step_in},
    eval_intv{eval_intv_in}, block_steps{block_steps_in}, max_separation{max_separation_in},
    min_angle{min_angle_in}, partners{partners_in}, partner_bounds{partner_bounds_in},
    tallies{tallies_in}, dist_acc{dist_acc_in}, dist_sqacc{dist_sqacc_in},
    all_dist_acc{all_dist_acc_in}, all_dist_sqacc{all_dist_sqacc_in}, angl_acc{angl_acc_in},
    angl_sqacc{angl_sqacc_in}
{}

//-------------------------------------------------------------------------------------------------
HBondReader::HBondReader(const int n_partners_in, const int stat_blocks_in, const int init_step_in,
                         const int eval_intv_in, const int block_steps_in,
                         const float max_separation_in, const float min_angle_in,
                         const int4* partners_in, const int* partner_bounds_in,
                         const int* tallies_in, const double* dist_acc_in,
                         const double* dist_sqacc_in, const double* all_dist_acc_in,
                         const double* all_dist_sqacc_in, const double* angl_acc_in,
                         const double* angl_sqacc_in) :
    n_partners{n_partners_in}, stat_blocks{stat_blocks_in}, init_step{init_step_in},
    eval_intv{eval_intv_in}, block_steps{block_steps_in}, max_separation{max_separation_in},
    min_angle{min_angle_in}, partners{partners_in}, partner_bounds{partner_bounds_in},
    tallies{tallies_in}, dist_acc{dist_acc_in}, dist_sqacc{dist_sqacc_in},
    all_dist_acc{all_dist_acc_in}, all_dist_sqacc{all_dist_sqacc_in}, angl_acc{angl_acc_in},
    angl_sqacc{angl_sqacc_in}
{}

//-------------------------------------------------------------------------------------------------
HBondReader::HBondReader(const HBondWriter &w) :
    n_partners{w.n_partners}, stat_blocks{w.stat_blocks}, init_step{w.init_step},
    eval_intv{w.eval_intv}, block_steps{w.block_steps}, max_separation{w.max_separation},
    min_angle{w.min_angle}, partners{w.partners}, partner_bounds{w.partner_bounds},
    tallies{w.tallies}, dist_acc{w.dist_acc}, dist_sqacc{w.dist_sqacc},
    all_dist_acc{w.all_dist_acc}, all_dist_sqacc{w.all_dist_sqacc}, angl_acc{w.angl_acc},
    angl_sqacc{w.angl_sqacc}
{}

//-------------------------------------------------------------------------------------------------
HBondReader::HBondReader(const HBondWriter *w) :
    n_partners{w->n_partners}, stat_blocks{w->stat_blocks}, init_step{w->init_step},
    eval_intv{w->eval_intv}, block_steps{w->block_steps}, max_separation{w->max_separation},
    min_angle{w->min_angle}, partners{w->partners}, partner_bounds{w->partner_bounds},
    tallies{w->tallies}, dist_acc{w->dist_acc}, dist_sqacc{w->dist_sqacc},
    all_dist_acc{w->all_dist_acc}, all_dist_sqacc{w->all_dist_sqacc}, angl_acc{w->angl_acc},
    angl_sqacc{w->angl_sqacc}
{}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis::HydrogenBondAnalysis(const PhaseSpaceSynthesis &poly_ps,
                                           const AtomGraphSynthesis &poly_ag,
                                           const SynthesisCacheMap &scmap,
                                           const AnalysisControls &trkcon,
                                           const int total_simulation_steps) :
    total_partners{0},
    statistical_block_count{trkcon.getHBondStatBlocks()},
    initiation_step{trkcon.getHBondInitiationStep()},
    evaluation_interval{trkcon.getHBondSamplingFrequency()},
    steps_per_block{calcStepsPerBlock(total_simulation_steps, trkcon)},
    maximum_separation{trkcon.getHBondMaxSeparation()},
    minimum_angle{trkcon.getHBondMinAngle()},
    minimum_occupancy{trkcon.getHBondOccupancyThreshold()},
    minimum_activity{trkcon.getHBondActivityThreshold()},
    base_varname{trkcon.getHBondVariableBase()},
    donor_in_first_group{},
    partner_bounds{HybridKind::POINTER, "hbond_partner_bounds"},
    bond_formations{HybridKind::POINTER, "hbond_track_bf"},
    distance_accumulators{HybridKind::POINTER, "hbond_track_dacc"},
    distance_sq_accumulators{HybridKind::POINTER, "hbond_track_dsqacc"},
    full_distance_acc{HybridKind::POINTER, "hbond_track_dacc"},
    full_distance_sq_acc{HybridKind::POINTER, "hbond_track_dsqacc"},
    angle_accumulators{HybridKind::POINTER, "hbond_track_anacc"},
    angle_sq_accumulators{HybridKind::POINTER, "hbond_track_ansqacc"},
    block_sample_counts{std::vector<double>(statistical_block_count)},
    poly_ps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps.getSelfPointer())},
    poly_ag_ptr{const_cast<AtomGraphSynthesis*>(poly_ag.getSelfPointer())},
    scmap_ptr{const_cast<SynthesisCacheMap*>(scmap.getSelfPointer())},
    int_storage{HybridKind::ARRAY, "hbond_idata"},
    dbl_storage{HybridKind::ARRAY, "hbond_ddata"}
{
  // Extract pairs of mask strings and their associated labels.  Find the relevant atom groups and
  // hydrogen bond donors or acceptors within them.  Enumerate all possible pairs along with the
  // polar hydrogens between them.
  const int hb_entries = trkcon.getHBondMaskCount();
  const SystemCache *sc = scmap.getCachePointer();
  std::vector<int4> possible_bonds;
  const int total_systems = poly_ps_ptr->getSystemCount();
  std::vector<std::vector<bool>> counted_as_donor_mask_i(total_systems);
  std::vector<std::vector<bool>> counted_as_donor_mask_ii(total_systems);
  std::vector<std::vector<bool>> counted_as_acceptor_mask_i(total_systems);
  std::vector<std::vector<bool>> counted_as_acceptor_mask_ii(total_systems);
  for (int i = 0; i < total_systems; i++) {
    const int i_natom = poly_ps_ptr->getAtomCount(i);
    counted_as_donor_mask_i[i] = std::vector<bool>(i_natom, false);
    counted_as_donor_mask_ii[i] = std::vector<bool>(i_natom, false);
    counted_as_acceptor_mask_i[i] = std::vector<bool>(i_natom, false);
    counted_as_acceptor_mask_ii[i] = std::vector<bool>(i_natom, false);
  }
  mask_i_donor_counts.resize(total_systems);
  mask_ii_donor_counts.resize(total_systems);
  mask_i_acceptor_counts.resize(total_systems);
  mask_ii_acceptor_counts.resize(total_systems);
  for (int i = 0; i < hb_entries; i++) {

    // Track the number of times this pair of hydrogen bonding masks finds donors and / or
    // acceptors.
    const std::string& i_mask_a = trkcon.getHBondMask(i, 0);
    const std::string& i_mask_b = trkcon.getHBondMask(i, 1);
    const std::string& i_label  = trkcon.getHBondMaskLabel(i);
    const std::vector<int> matching_systems = scmap.getLabelGroup(i_label);
    const size_t i_nsys = matching_systems.size();
    for (size_t j = 0; j < i_nsys; j++) {

      // Identify the system and get references to its details
      const int ij_cache_idx = scmap.getSystemCacheIndex(matching_systems[j]);
      const AtomGraph *ij_ag = poly_ag.getSystemTopologyPointer(matching_systems[j]);
      const int ij_offset = poly_ag.getAtomOffset(matching_systems[j]);
      const ChemicalFeatures& ij_chemfe = sc->getFeatures(ij_cache_idx);
      const CoordinateFrame cf = poly_ps.exportCoordinates(matching_systems[j],
                                                           HybridFormat::HOST_ONLY);
      const NonbondedKit<double> nbk = ij_ag->getDoublePrecisionNonbondedKit();
      const ChemicalDetailsKit cdk = ij_ag->getChemicalDetailsKit();

      // Apply the atom masks from user input to this system
      const AtomMask selection_a(i_mask_a, ij_ag, ij_chemfe, cf);
      const AtomMask selection_b(i_mask_b, ij_ag, ij_chemfe, cf);
      const std::vector<uint> donor_bitwise    = ij_chemfe.getHydrogenBondDonorMask();
      const std::vector<uint> acceptor_bitwise = ij_chemfe.getHydrogenBondAcceptorMask();
      const size_t n_bitblock = donor_bitwise.size();
      const std::vector<int> sela_ids = selection_a.getMaskedAtomList();
      const std::vector<int> selb_ids = selection_b.getMaskedAtomList();
      const int sela_count = sela_ids.size();
      const int selb_count = selb_ids.size();

      // Find donors and acceptors from within the selections.
      std::vector<int> sela_donors, selb_donors, sela_acceptors, selb_acceptors;
      std::vector<int> sela_polarh, selb_polarh;
      for (int k = 0; k < sela_count; k++) {
        if (readBitFromMask(donor_bitwise, sela_ids[k])) {
          sela_donors.push_back(sela_ids[k]);
          
        }
        if (readBitFromMask(acceptor_bitwise, sela_ids[k])) {
          sela_acceptors.push_back(sela_ids[k]);
        }
      }
      for (int k = 0; k < selb_count; k++) {
        if (readBitFromMask(donor_bitwise, selb_ids[k])) {
          selb_donors.push_back(selb_ids[k]);
        }
        if (readBitFromMask(acceptor_bitwise, selb_ids[k])) {
          selb_acceptors.push_back(selb_ids[k]);
        }
      }
      const int na_donor = sela_donors.size();
      const int nb_donor = selb_donors.size();
      const int na_acceptor = sela_acceptors.size();
      const int nb_acceptor = selb_acceptors.size();

      // Make notes of the numbers of donors and acceptors, per system, in each mask.  This will be
      // used in reporting, although the information is not conveyed in the abstract.
      for (int k = 0; k < na_donor; k++) {
        counted_as_donor_mask_i[j][sela_donors[k]] = true;
      }
      for (int k = 0; k < nb_donor; k++) {
        counted_as_donor_mask_ii[j][selb_donors[k]] = true;
      }
      for (int k = 0; k < na_acceptor; k++) {
        counted_as_acceptor_mask_i[j][sela_acceptors[k]] = true;
      }
      for (int k = 0; k < nb_acceptor; k++) {
        counted_as_acceptor_mask_ii[j][selb_acceptors[k]] = true;
      }
      
      // List all possible combinations of donors and acceptors.  Find the relevant polar hydrogens
      // based on the donor indices.  If there are multiple polar hydrogens for a particular donor,
      // list all possible combinations of the donor and its hydrogens crossed with any applicable
      // acceptors.
      for (int k = 0; k < na_donor; k++) {
        const int katom = sela_donors[k];
        for (int m = nbk.nb12_bounds[katom]; m < nbk.nb12_bounds[katom + 1]; m++) {
          const int matom = nbk.nb12x[m];
          if (cdk.z_numbers[matom] == 1) {
            for (int km_acc = 0; km_acc < nb_acceptor; km_acc++) {
              possible_bonds.push_back({ katom + ij_offset,
                                         matom + ij_offset,
                                         selb_acceptors[km_acc] + ij_offset,
                                         matching_systems[j] });
              donor_in_first_group.push_back(true);
            }
          }
        }
      }
      for (int k = 0; k < nb_donor; k++) {
        const int katom = selb_donors[k];
        for (int m = nbk.nb12_bounds[katom]; m < nbk.nb12_bounds[katom + 1]; m++) {
          const int matom = nbk.nb12x[m];
          if (cdk.z_numbers[matom] == 1) {
            for (int km_acc = 0; km_acc < na_acceptor; km_acc++) {
              possible_bonds.push_back({ katom + ij_offset,
                                         matom + ij_offset,
                                         sela_acceptors[km_acc] + ij_offset,
                                         matching_systems[j] });
              donor_in_first_group.push_back(false);
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < total_systems; i++) {
    mask_i_donor_counts[i]  = sum<int>(counted_as_donor_mask_i[i]);
    mask_ii_donor_counts[i] = sum<int>(counted_as_donor_mask_ii[i]);
    mask_i_acceptor_counts[i]  = sum<int>(counted_as_acceptor_mask_i[i]);
    mask_ii_acceptor_counts[i] = sum<int>(counted_as_acceptor_mask_ii[i]);
  }
  total_partners = possible_bonds.size();
  allocate();
  partners.putHost(possible_bonds);

  // Compute the expected numbers of samples in each statistical block
  for (int i = 0; i < statistical_block_count; i++) {
    const int first_step = initiation_step + (i * steps_per_block);
    const int final_step = initiation_step + ((i + 1) * steps_per_block);
    const int first_record = roundUp<int>(first_step, evaluation_interval);
    block_sample_counts[i] = static_cast<double>(((final_step - first_record) +
                                                  evaluation_interval - 1) /
                                                 evaluation_interval);
  }
}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis::HydrogenBondAnalysis(const PhaseSpaceSynthesis &poly_ps,
                                           const AtomGraphSynthesis &poly_ag,
                                           const SynthesisCacheMap &scmap,
                                           const AnalysisControls &trkcon,
                                           const DynamicsControls &dyncon) :
    HydrogenBondAnalysis(poly_ps, poly_ag, scmap, trkcon, dyncon.getStepCount())
{}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis::HydrogenBondAnalysis(const PhaseSpaceSynthesis &poly_ps,
                                           const AtomGraphSynthesis &poly_ag,
                                           const SynthesisCacheMap &scmap,
                                           const AnalysisControls &trkcon,
                                           const MinimizeControls &mincon) :
    HydrogenBondAnalysis(poly_ps, poly_ag, scmap, trkcon, mincon.getTotalCycles())
{}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis::HydrogenBondAnalysis(const HydrogenBondAnalysis &original) :
    total_partners{original.total_partners},
    statistical_block_count{original.statistical_block_count},
    initiation_step{original.initiation_step},
    evaluation_interval{original.evaluation_interval},
    steps_per_block{original.steps_per_block},
    maximum_separation{original.maximum_separation},
    minimum_angle{original.minimum_angle},
    minimum_occupancy{original.minimum_occupancy},
    minimum_activity{original.minimum_activity},
    base_varname{original.base_varname},
    donor_in_first_group{original.donor_in_first_group},
    mask_i_donor_counts{original.mask_i_donor_counts},
    mask_i_acceptor_counts{original.mask_i_acceptor_counts},
    mask_ii_donor_counts{original.mask_ii_donor_counts},
    mask_ii_acceptor_counts{original.mask_ii_acceptor_counts},
    partners{original.partners},
    partner_bounds{original.partner_bounds},
    bond_formations{original.bond_formations},
    distance_accumulators{original.distance_accumulators},
    distance_sq_accumulators{original.distance_sq_accumulators},
    full_distance_acc{original.full_distance_acc},
    full_distance_sq_acc{original.full_distance_sq_acc},
    angle_accumulators{original.angle_accumulators},
    angle_sq_accumulators{original.angle_sq_accumulators},
    block_sample_counts{original.block_sample_counts},
    poly_ps_ptr{original.poly_ps_ptr},
    poly_ag_ptr{original.poly_ag_ptr},
    scmap_ptr{original.scmap_ptr},
    int_storage{original.int_storage},
    dbl_storage{original.dbl_storage}
{
  // As with some other classes, the allocator also serves to rebase POINTER-kind Hybrid objects.
  // The underlying storage arrays will not be further resized, as the initial copy constructions
  // have already set them to the appropriate sizes, with the necessary data content.
  allocate();
}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis::HydrogenBondAnalysis(HydrogenBondAnalysis &&original) :
    total_partners{original.total_partners},
    statistical_block_count{original.statistical_block_count},
    initiation_step{original.initiation_step},
    evaluation_interval{original.evaluation_interval},
    steps_per_block{original.steps_per_block},
    maximum_separation{original.maximum_separation},
    minimum_angle{original.minimum_angle},
    minimum_occupancy{original.minimum_occupancy},
    minimum_activity{original.minimum_activity},
    base_varname{std::move(original.base_varname)},
    donor_in_first_group{std::move(original.donor_in_first_group)},
    mask_i_donor_counts{std::move(original.mask_i_donor_counts)},
    mask_i_acceptor_counts{std::move(original.mask_i_acceptor_counts)},
    mask_ii_donor_counts{std::move(original.mask_ii_donor_counts)},
    mask_ii_acceptor_counts{std::move(original.mask_ii_acceptor_counts)},
    partners{std::move(original.partners)},
    partner_bounds{std::move(original.partner_bounds)},
    bond_formations{std::move(original.bond_formations)},
    distance_accumulators{std::move(original.distance_accumulators)},
    distance_sq_accumulators{std::move(original.distance_sq_accumulators)},
    full_distance_acc{std::move(original.full_distance_acc)},
    full_distance_sq_acc{std::move(original.full_distance_sq_acc)},
    angle_accumulators{std::move(original.angle_accumulators)},
    angle_sq_accumulators{std::move(original.angle_sq_accumulators)},
    block_sample_counts{std::move(original.block_sample_counts)},
    poly_ps_ptr{original.poly_ps_ptr},
    poly_ag_ptr{original.poly_ag_ptr},
    scmap_ptr{original.scmap_ptr},
    int_storage{std::move(original.int_storage)},
    dbl_storage{std::move(original.dbl_storage)}
{}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis& HydrogenBondAnalysis::operator=(const HydrogenBondAnalysis &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // As with the typical copy assignment statement, set scalar values and the values of pointers
  // based on the other object.  Copy POINTER-kind Hybrids and do a deep copy of each ARRAY-kind
  // Hybrid object.
  total_partners = other.total_partners;
  statistical_block_count = other.statistical_block_count;
  initiation_step = other.initiation_step;
  evaluation_interval = other.evaluation_interval;
  steps_per_block = other.steps_per_block;
  maximum_separation = other.maximum_separation;
  minimum_angle = other.minimum_angle;
  minimum_occupancy = other.minimum_occupancy;
  minimum_activity = other.minimum_activity;
  base_varname = other.base_varname;
  donor_in_first_group = other.donor_in_first_group;
  mask_i_donor_counts = other.mask_i_donor_counts;
  mask_i_acceptor_counts = other.mask_i_acceptor_counts;
  mask_ii_donor_counts = other.mask_ii_donor_counts;
  mask_ii_acceptor_counts = other.mask_ii_acceptor_counts;
  partners = other.partners;
  partner_bounds = other.partner_bounds;
  bond_formations = other.bond_formations;
  distance_accumulators = other.distance_accumulators;
  distance_sq_accumulators = other.distance_sq_accumulators;
  full_distance_acc = other.full_distance_acc;
  full_distance_sq_acc = other.full_distance_sq_acc;
  angle_accumulators = other.angle_accumulators;
  angle_sq_accumulators = other.angle_sq_accumulators;
  block_sample_counts = other.block_sample_counts;
  poly_ps_ptr = other.poly_ps_ptr;
  poly_ag_ptr = other.poly_ag_ptr;
  scmap_ptr = other.scmap_ptr;
  int_storage = other.int_storage;
  dbl_storage = other.dbl_storage;

  // The allocation function again resets POINTER-kind Hybrids as needed
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis& HydrogenBondAnalysis::operator=(HydrogenBondAnalysis &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  total_partners = other.total_partners;
  statistical_block_count = other.statistical_block_count;
  initiation_step = other.initiation_step;
  evaluation_interval = other.evaluation_interval;
  steps_per_block = other.steps_per_block;
  maximum_separation = other.maximum_separation;
  minimum_angle = other.minimum_angle;
  minimum_occupancy = other.minimum_occupancy;
  minimum_activity = other.minimum_activity;
  base_varname = std::move(other.base_varname);
  donor_in_first_group = std::move(other.donor_in_first_group);
  mask_i_donor_counts = std::move(other.mask_i_donor_counts);
  mask_i_acceptor_counts = std::move(other.mask_i_acceptor_counts);
  mask_ii_donor_counts = std::move(other.mask_ii_donor_counts);
  mask_ii_acceptor_counts = std::move(other.mask_ii_acceptor_counts);
  partners = std::move(other.partners);
  partner_bounds = std::move(other.partner_bounds);
  bond_formations = std::move(other.bond_formations);
  distance_accumulators = std::move(other.distance_accumulators);
  distance_sq_accumulators = std::move(other.distance_sq_accumulators);
  full_distance_acc = std::move(other.full_distance_acc);
  full_distance_sq_acc = std::move(other.full_distance_sq_acc);
  angle_accumulators = std::move(other.angle_accumulators);
  angle_sq_accumulators = std::move(other.angle_sq_accumulators);
  block_sample_counts = std::move(other.block_sample_counts);
  poly_ps_ptr = other.poly_ps_ptr;
  poly_ag_ptr = other.poly_ag_ptr;
  scmap_ptr = other.scmap_ptr;
  int_storage = std::move(other.int_storage);
  dbl_storage = std::move(other.dbl_storage);
  return *this;
}

//-------------------------------------------------------------------------------------------------
int HydrogenBondAnalysis::getTotalCandidates(const int system_index) const {
  return partner_bounds.readHost(system_index + 1) - partner_bounds.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
int HydrogenBondAnalysis::getTotalCandidates() const {
  return total_partners;
}

//-------------------------------------------------------------------------------------------------
int2 HydrogenBondAnalysis::getSystemPartnerBounds(const int system_index) const {
  return { partner_bounds.readHost(system_index), partner_bounds.readHost(system_index + 1) };
}

//-------------------------------------------------------------------------------------------------
int4 HydrogenBondAnalysis::getPartnerSynthesisIndices(const int pair_index) const {

  // Intercept bad inputs and print a more descriptive error message than the Hybrid bounds check
  // would produce.
  if (pair_index < 0 || pair_index >= total_partners) {
    rtErr("Partner index " + std::to_string(pair_index) + " is invalid for collection of " +
          std::to_string(total_partners) + " possible hydrogen bonds in a synthesis of " +
          std::to_string(poly_ps_ptr->getSystemCount()) + " systems.", "HydrogenBondAnalysis",
          "getPartnerSynthesisIndices");
  }
  return partners.readHost(pair_index);
}
  
//-------------------------------------------------------------------------------------------------
int HydrogenBondAnalysis::getBlockAveragingCount() const {
  return statistical_block_count;
}

//-------------------------------------------------------------------------------------------------
int HydrogenBondAnalysis::getEvaluationInterval() const {
  return evaluation_interval;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* HydrogenBondAnalysis::getCoordinateSynthesisPointer() const {
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
const AtomGraphSynthesis* HydrogenBondAnalysis::getTopologySynthesisPointer() const {
  return poly_ag_ptr;
}

//-------------------------------------------------------------------------------------------------
HBondWriter HydrogenBondAnalysis::data(const HybridTargetLevel tier) {
  return HBondWriter(total_partners, statistical_block_count, initiation_step, evaluation_interval,
                     steps_per_block, maximum_separation, minimum_angle, partners.data(tier),
                     partner_bounds.data(tier), bond_formations.data(tier),
                     distance_accumulators.data(tier), distance_sq_accumulators.data(tier),
                     full_distance_acc.data(tier), full_distance_sq_acc.data(tier),
                     angle_accumulators.data(tier), angle_sq_accumulators.data(tier));
}
  
//-------------------------------------------------------------------------------------------------
const HBondReader HydrogenBondAnalysis::data(const HybridTargetLevel tier) const {
  return HBondReader(total_partners, statistical_block_count, initiation_step, evaluation_interval,
                     steps_per_block, maximum_separation, minimum_angle, partners.data(tier),
                     partner_bounds.data(tier), bond_formations.data(tier),
                     distance_accumulators.data(tier), distance_sq_accumulators.data(tier),
                     full_distance_acc.data(tier), full_distance_sq_acc.data(tier),
                     angle_accumulators.data(tier), angle_sq_accumulators.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::upload() {
  partners.upload();
  int_storage.upload();
  dbl_storage.upload();
}

//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::download() {
  partners.download();
  int_storage.download();
  dbl_storage.download();
}

//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::uploadDefinitions() {
  partners.upload();
  partner_bounds.upload();
}

//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::downloadDefinitions() {
  partners.download();
  partner_bounds.download();
}

//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::uploadAccumulators() {
  bond_formations.upload();
  dbl_storage.upload();
}

//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::downloadAccumulators() {
  bond_formations.download();
  dbl_storage.download();
}
#endif
  
//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::initialize(const HybridTargetLevel tier) {
  const size_t tp_zu = roundUp<int>(total_partners, warp_size_int);
  const size_t fsbc_zu = static_cast<size_t>(4 * statistical_block_count) * tp_zu;
                         
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      int* bond_formation_ptr = bond_formations.data();
      double* storage_ptr = dbl_storage.data();
      for (size_t i = 0; i < tp_zu; i++) {
        bond_formation_ptr[i] = 0;
      }
      for (size_t i = 0; i < fsbc_zu; i++) {
        storage_ptr[i] = 0.0;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    if (cudaMemset(bond_formations.data(HybridTargetLevel::DEVICE), 0, tp_zu * sizeof(int)) !=
        cudaSuccess) {
      rtErr("Error in initializing " + std::string(bond_formations.getLabel().name) + ".",
            "HydrogenBondAnalysis", "initialize");
    }
    if (cudaMemset(dbl_storage.data(HybridTargetLevel::DEVICE), 0, fsbc_zu * sizeof(double)) !=
        cudaSuccess) {
      rtErr("Error in initializing " + std::string(dbl_storage.getLabel().name) + ".",
            "HydrogenBondAnalysis", "initialize");
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
SectionContents HydrogenBondAnalysis::reportResults(const ReportControls &repcon) const {
  
  // Scan for successful hydrogen bonds and systems containing successful hydrogen bonds.  To
  // begin, record the total number of successful bonds found.
  const HBondReader hbr = this->data();
  std::vector<double> occupancies(hbr.stat_blocks);
  std::vector<int> successes;
  const int nsys = poly_ag_ptr->getSystemCount();
  std::vector<int> systemwide_total_participants(nsys, 0);
  std::vector<int> systemwide_strong_participants(nsys, 0);
  std::vector<double> systemwide_bonding(hbr.stat_blocks * nsys, 0.0);
  for (int i = 0; i < total_partners; i++) {
    bool hbond_exists = false;
    for (int j = 0; j < hbr.stat_blocks; j++) {
      const int item_idx = (j * total_partners) + i;
      occupancies[j] = static_cast<double>(hbr.tallies[item_idx]) / block_sample_counts[j];
      hbond_exists = (hbond_exists || hbr.tallies[item_idx] > 0);
    }
    if (hbond_exists) {

      // Add the occupancy to system-as-a-whole counters
      const int4 hb_atoms = hbr.partners[i];
      for (int j = 0; j < hbr.stat_blocks; j++) {
        systemwide_bonding[(j * nsys) + hb_atoms.w] += occupancies[j];
      }
      systemwide_total_participants[hb_atoms.w] += 1;
      
      // Queue the arrangement to receive its own note in the output
      if (mean(occupancies) >= minimum_occupancy) {
        successes.push_back(i);
        systemwide_strong_participants[hb_atoms.w] += 1;
      }
    }
  }
  const int n_found = successes.size();
  const SystemCache* sc_ptr = scmap_ptr->getCachePointer();

  // Most output formats will not tolerate a "+/-" symbol and remain interpretable by the intended
  // third-party packages.
  std::string pm_symb;
  switch (repcon.getOutputSyntax()) {
  case OutputSyntax::MATPLOTLIB:
  case OutputSyntax::MATRIX_PKG:
  case OutputSyntax::JSON:
    pm_symb = std::string(" ");
    break;
  case OutputSyntax::STANDALONE:
    pm_symb = std::string(" +/- ");
    break;
  }
  
  // Allocate tabular data, then fill the tables as appropriate.  The first table is a summary of
  // overall hydrogen bonding and is ordered in terms of overall hydrogen bond formation throughout
  // each simulation, with the strongest systems first.  It has ten columns:
  //
  // 1.)  System index
  // 2.)  Overall hydrogen bond formation rate (total occupancy summed over all possible partners)
  // 3.)  Total hydrogen bond partners participating at any time in the simulation
  // 4.)  Number of "strong" hydrogen bonds observed--those exceeding the occupancy threshold for
  //      individual reporting
  // 5.)  Number of hydrogen bond donors in the first mask
  // 6.)  Number of hydrogen bond acceptors in the first mask
  // 7.)  Number of hydrogen bond donors in the second mask
  // 8.)  Number of hydrogen bond acceptors in the second mask
  // 9.)  Name of the original topology (protected by a comment character)
  // 10.) Name of the original coordinates file
  const double inv_block_count = 1.0 / static_cast<double>(hbr.stat_blocks);
  std::vector<CombineIDp> hbond_best;
  for (int i = 0; i < nsys; i++) {
    double swide_bond_sum = 0.0;
    for (int j = 0; j < hbr.stat_blocks; j++) {
      swide_bond_sum += systemwide_bonding[(j * nsys) + i];
    }
    swide_bond_sum *= inv_block_count;
    if (swide_bond_sum >= minimum_activity) {
      hbond_best.push_back({ i, swide_bond_sum });
    }
  }
  const int n_noteworthy = hbond_best.size();
  std::sort(hbond_best.begin(), hbond_best.end(),
            [](CombineIDp a, CombineIDp b) { return a.y > b.y; });
  std::vector<std::string> summary_data(n_noteworthy * 10);
  const std::string protector = std::string(1, commentSymbol(repcon.getOutputSyntax())) + " ";
  std::vector<std::string> best_topology_names, best_coordinate_names;
  best_topology_names.reserve(n_noteworthy);
  best_coordinate_names.reserve(n_noteworthy);
  for (int i = 0; i < n_noteworthy; i++) {
    const int i_case = hbond_best[i].x;
    const int i_cache_idx = scmap_ptr->getSystemCacheIndex(i_case);
    best_topology_names.push_back(poly_ps_ptr->getSystemTopologyPointer(i_case)->getFileName());
    best_coordinate_names.push_back(sc_ptr->getCoordinatePointer(i_cache_idx)->getFileName());
  }
  std::string before, after;
  std::vector<std::string> topology_suffixes, coordinate_suffixes;
  for (int i = 0; i < n_noteworthy; i++) {
    splitPath(best_topology_names[i], &before, &after);
    best_topology_names[i] = before;
    if (findStringInVector(topology_suffixes, after) == topology_suffixes.size()) {
      topology_suffixes.push_back(after);
    }
    splitPath(best_coordinate_names[i], &before, &after);
    best_coordinate_names[i] = before;
    if (findStringInVector(coordinate_suffixes, after) == coordinate_suffixes.size()) {
      coordinate_suffixes.push_back(after);
    }
  }
  const std::vector<std::string> topl_prefixes = extractCommonPaths(&best_topology_names);
  const std::vector<std::string> coord_prefixes = extractCommonPaths(&best_coordinate_names);
  std::vector<double> swide_buffer(hbr.stat_blocks);
  for (int i = 0; i < n_noteworthy; i++) {
    const int i_case = hbond_best[i].x;
    summary_data[i] = std::to_string(i_case);
    summary_data[i +       n_noteworthy] = realToString(hbond_best[i].y, 7, 4,
                                                        NumberFormat::STANDARD_REAL);
    if (hbr.stat_blocks > 1) {
      for (int j = 0; j < hbr.stat_blocks; j++) {
        swide_buffer[j] = systemwide_bonding[(j * nsys) + i_case];
      }
      const double swide_stdev = variance(swide_buffer, VarianceMethod::STANDARD_DEVIATION);
      summary_data[i + n_noteworthy] += pm_symb + realToString(swide_stdev, 7, 4,
                                                               NumberFormat::STANDARD_REAL);
    }
    summary_data[i + (2 * n_noteworthy)] = std::to_string(systemwide_total_participants[i_case]);
    summary_data[i + (3 * n_noteworthy)] = std::to_string(systemwide_strong_participants[i_case]);
    summary_data[i + (4 * n_noteworthy)] = std::to_string(mask_i_donor_counts[i_case]);
    summary_data[i + (5 * n_noteworthy)] = std::to_string(mask_i_acceptor_counts[i_case]);
    summary_data[i + (6 * n_noteworthy)] = std::to_string(mask_ii_donor_counts[i_case]);
    summary_data[i + (7 * n_noteworthy)] = std::to_string(mask_ii_acceptor_counts[i_case]);
    summary_data[i + (8 * n_noteworthy)] = protector + best_topology_names[i];
    summary_data[i + (9 * n_noteworthy)] = best_coordinate_names[i];
  }
  std::vector<std::string> summary_headings = {
    "System", "Average Hydrogen Bonding", "Unique Hydrogen Bonds", "Strong Hydrogen Bonds",
    "Mask I Donors", "Mask I Acceptors", "Mask II Donors", "Mask II Acceptors", "Topology",
    "Coordinates"
  };
  std::vector<JustifyText> short_justifications(10, JustifyText::RIGHT);
  for (int i = 8; i < 10; i++) {
    short_justifications[i] = JustifyText::LEFT;
  }
  ReportTable short_list(summary_data, summary_headings, base_varname + "_highlights", 300,
                         short_justifications);
  short_list.unprotectContent();
  
  // The comprehensive table of all hydrogen bonds will have up to sixteen columns:
  //
  // 1.)  System index
  // 2.)  Atom index of the first atom, within the system (not the synthesis)
  // 3.)  Atom index of the polar hydrogen, within the system
  // 4.)  Atom index of the second atom, within the system
  // 5.)  Occupancy, as a percentage of the total opportunities to appear in the data.  If block
  //      averaging is in effect, this occupancy will be written with an error bar based on the
  //      standard deviation among all blocks.
  // 6.)  Mean donor : acceptor distance of the hydrogen bond when it is formed, with standard
  //      deviation reported as an error bar
  // 7.)  Mean donor - proton : acceptor angle, with standard deviation
  // 8.)  Mean donor : acceptor distance out of all the data, with standard deviation
  // 9.)  Standard deviation of the mean donor - acceptor distance as computed from block averaging
  //      (if applicable)
  // 10.) Standard deviation of the mean donor - proton : acceptor angle, as computed from block
  //      averaging (if applicable)
  //
  // The following columns appear behind a comment character and are thus omitted from matrix data
  // when the output is interpreted by some third-party package
  // 11.) Residue name, residue number, and atom name of the donor or acceptor from the first group
  // 12.) A double dash (--) or double colon (::) to indicate that the first group's atom was a
  //      proton donor or acceptor, respectively.
  // 13.) Residue name, residue number, and atom name of the proton from the first or second group
  // 14.) A double colon (::) or double dash (--), reflexive of the contents in the second column
  // 15.) Residue name, residue number, and atom name of the donor or acceptor from the second
  //      group
  // 16.) System label of the system where the hydrogen bond is found.  This may denote a group of
  //      molecular systems.
  std::vector<std::string> comprehensive_data(n_found * (14 + (2 * (hbr.stat_blocks > 1))));
  const int max_res = poly_ag_ptr->getLargestResidueCount();
  int max_res_digits = ceil(log10(max_res));
  const int ten_test = pow(10.0, max_res_digits) + 0.5;
  max_res_digits += (max_res == ten_test);
  const std::string double_dash("--");
  const std::string double_coln("::");
  const double total_samples = sum<double>(block_sample_counts);
  for (int i = 0; i < n_found; i++) {
    const int scc_i = successes[i];
    const int4 hb_atoms = hbr.partners[scc_i];
    comprehensive_data[i          ] = std::to_string(hb_atoms.w);

    // Compute the occupancy across all statistical blocks, even if there is just one.
    int total_tally = 0;
    for (int j = 0; j < hbr.stat_blocks; j++) {
      total_tally += hbr.tallies[(j * hbr.n_partners) + scc_i];
    }
    const double occupancy = static_cast<double>(total_tally) / total_samples;
    comprehensive_data[i + (4 * n_found)] = realToString(occupancy, 7, 4,
                                                         NumberFormat::STANDARD_REAL);
    if (hbr.stat_blocks > 1) {
      std::vector<double> occ_binned(hbr.stat_blocks);
      for (int j = 0; j < hbr.stat_blocks; j++) {
        occ_binned[j] = static_cast<double>(hbr.tallies[(j * hbr.n_partners) + scc_i]) /
                        block_sample_counts[j];
      }
      const double occ_dev = variance(occ_binned, VarianceMethod::STANDARD_DEVIATION);
      comprehensive_data[i + (4 * n_found)] += pm_symb + realToString(occ_dev, 7, 4,
                                                                      NumberFormat::STANDARD_REAL);
    }

    // Make super-sums of the distance and angle accumulators
    double dist_acc_sum = 0.0;
    double dist_sqacc_sum = 0.0;
    double angl_acc_sum = 0.0;
    double angl_sqacc_sum = 0.0;
    double full_acc_sum = 0.0;
    double full_sqacc_sum = 0.0;
    for (int j = 0; j < hbr.stat_blocks; j++) {
      dist_acc_sum   += hbr.dist_acc[(j * hbr.n_partners) + scc_i];
      dist_sqacc_sum += hbr.dist_sqacc[(j * hbr.n_partners) + scc_i];
      full_acc_sum   += hbr.all_dist_acc[(j * hbr.n_partners) + scc_i];
      full_sqacc_sum += hbr.all_dist_sqacc[(j * hbr.n_partners) + scc_i];
      angl_acc_sum   += hbr.angl_acc[(j * hbr.n_partners) + scc_i];
      angl_sqacc_sum += hbr.angl_sqacc[(j * hbr.n_partners) + scc_i];
    }
    const double mean_formed_da = dist_acc_sum / static_cast<double>(total_tally);
    const double std_formed_da = runningVariance(dist_sqacc_sum, dist_acc_sum, total_tally,
                                                 VarianceMethod::STANDARD_DEVIATION);
    comprehensive_data[i + (5 * n_found)] = realToString(mean_formed_da, 7, 4,
                                                         NumberFormat::STANDARD_REAL) + pm_symb +
                                            realToString(std_formed_da, 7, 4,
                                                         NumberFormat::STANDARD_REAL);
    const double mean_formed_dah = angl_acc_sum / static_cast<double>(total_tally);
    const double std_formed_dah = runningVariance(angl_sqacc_sum, angl_acc_sum, total_tally,
                                                  VarianceMethod::STANDARD_DEVIATION);
    comprehensive_data[i + (6 * n_found)] = realToString(mean_formed_dah * 180.0 / pi, 7, 2,
                                                         NumberFormat::STANDARD_REAL) + pm_symb +
                                            realToString(std_formed_dah * 180.0 / pi, 7, 2,
                                                         NumberFormat::STANDARD_REAL);
    const double mean_full_da = full_acc_sum / total_samples;
    const double std_full_da = runningVariance(full_sqacc_sum, full_acc_sum, total_samples,
                                               VarianceMethod::STANDARD_DEVIATION);
    comprehensive_data[i + (7 * n_found)] = realToString(mean_full_da, 7, 4,
                                                         NumberFormat::STANDARD_REAL) + pm_symb +
                                            realToString(std_full_da, 7, 4,
                                                         NumberFormat::STANDARD_REAL);
    if (hbr.stat_blocks > 1) {
      std::vector<double> dist_binned(hbr.stat_blocks), angl_binned(hbr.stat_blocks);
      for (int j = 0; j < hbr.stat_blocks; j++) {
        dist_binned[j] = static_cast<double>(hbr.dist_acc[(j * hbr.n_partners) + scc_i]) /
                         block_sample_counts[j];
        angl_binned[j] = static_cast<double>(hbr.angl_acc[(j * hbr.n_partners) + scc_i]) /
                         block_sample_counts[j];
      }
      const double std_block_da  = variance(dist_binned, VarianceMethod::STANDARD_DEVIATION);
      const double std_block_dah = variance(angl_binned, VarianceMethod::STANDARD_DEVIATION);
      comprehensive_data[i + (8 * n_found)] = realToString(std_block_da, 7, 4,
                                                           NumberFormat::STANDARD_REAL);
      comprehensive_data[i + (9 * n_found)] = realToString(std_block_dah, 7, 4,
                                                           NumberFormat::STANDARD_REAL);
    }
    const int stb_offset = 2 * (hbr.stat_blocks > 1);
    
    // The names of atoms are formatted for human comprehensibility
    const int sys_offset = poly_ag_ptr->getAtomOffset(hb_atoms.w);
    const int donor_idx = hb_atoms.x - sys_offset;
    const int protn_idx = hb_atoms.y - sys_offset;
    const int accpt_idx = hb_atoms.z - sys_offset;
    const AtomGraph *iag = poly_ag_ptr->getSystemTopologyPointer(hb_atoms.w);
    const int donor_ridx = iag->getResidueIndex(donor_idx);
    const int protn_ridx = iag->getResidueIndex(protn_idx);
    const int accpt_ridx = iag->getResidueIndex(accpt_idx);
    std::string atom_i, atom_ii;
    int atom_i_idx, atom_ii_idx;
    if (donor_in_first_group[successes[i]]) {
      atom_i = std::string("% ") +
               addTailingWhiteSpace(char4ToString(iag->getResidueName(donor_ridx)), 4) + " " +
               intToString(donor_ridx, max_res_digits) + " " +
               addTailingWhiteSpace(char4ToString(iag->getAtomName(donor_idx)), 4);
      atom_ii = addTailingWhiteSpace(char4ToString(iag->getResidueName(accpt_ridx)), 4) + " " +
                intToString(accpt_ridx, max_res_digits) + " " +
                addTailingWhiteSpace(char4ToString(iag->getAtomName(accpt_idx)), 4);
      atom_i_idx  = donor_idx;
      atom_ii_idx = accpt_idx;
    }
    else {
      atom_i = std::string("% ") +
               addTailingWhiteSpace(char4ToString(iag->getResidueName(accpt_ridx)), 4) + " " +
               intToString(accpt_ridx, max_res_digits) + " " +
               addTailingWhiteSpace(char4ToString(iag->getAtomName(accpt_idx)), 4);
      atom_ii = addTailingWhiteSpace(char4ToString(iag->getResidueName(donor_ridx)), 4) + " " +
                intToString(donor_ridx, max_res_digits) + " " +
                addTailingWhiteSpace(char4ToString(iag->getAtomName(donor_idx)), 4);
      atom_i_idx  = accpt_idx;
      atom_ii_idx = donor_idx;
    }
    const std::string atom_h = char4ToString(iag->getResidueName(protn_ridx)) + " " +
                               intToString(protn_ridx, max_res_digits) + " " +
                               char4ToString(iag->getAtomName(protn_idx));
    comprehensive_data[i +       n_found] = std::to_string(atom_i_idx);
    comprehensive_data[i + (2 * n_found)] = std::to_string(protn_idx);
    comprehensive_data[i + (3 * n_found)] = std::to_string(atom_ii_idx);
    comprehensive_data[i + ((8  + stb_offset) * n_found)] = atom_i;
    comprehensive_data[i + ((10 + stb_offset) * n_found)] = atom_h;
    comprehensive_data[i + ((12 + stb_offset) * n_found)] = atom_ii;
    if (donor_in_first_group[successes[i]]) {
      comprehensive_data[i + ((9  + stb_offset) * n_found)] = double_dash;
      comprehensive_data[i + ((11 + stb_offset) * n_found)] = double_coln;
    }
    else {
      comprehensive_data[i + ((9  + stb_offset) * n_found)] = double_coln;
      comprehensive_data[i + ((11 + stb_offset) * n_found)] = double_dash;
    }

    // Trace the system in question back to the user's input systems cache
    const int lbl_idx = scmap_ptr->getLabelCacheIndex(hb_atoms.w);
    comprehensive_data[i + ((13 + stb_offset) * n_found)] = sc_ptr->getLabel(lbl_idx);
  }

  // Prepare the column justifications for tabulated output
  std::vector<JustifyText> justifications(8 + (2 * (hbr.stat_blocks > 1)), JustifyText::RIGHT);
  for (int i = 0; i < 6; i++) {
    justifications.push_back(JustifyText::LEFT);
  }
  std::vector<std::string> comprehensive_headings = {
    "System Index", "Mask I Atom Index", "Polar H Index", "Mask II Atom Index",
    "Occupancy, Mean / St.Dev.", "Dn::Acc Distance, When Formed, Mean / St.Dev.",
    "Dn--H::Acc Angle, When Formed, Mean / St.Dev.", "Dn::Acc Distance, Overall, Mean / St.Dev."
  };
  if (hbr.stat_blocks > 1) {
    comprehensive_headings.push_back("Dn::Acc Distance, St.Dev. by Block Averaging");
    comprehensive_headings.push_back("Dn--H::Acc Angle, St.Dev. by Block Averaging");
  }
  comprehensive_headings.push_back("Atom, Mask I");
  comprehensive_headings.push_back(" ");
  comprehensive_headings.push_back("Polar H");
  comprehensive_headings.push_back(" ");
  comprehensive_headings.push_back("Atom, Mask II");
  comprehensive_headings.push_back("System Label");
  ReportTable comprehensive_list(comprehensive_data, comprehensive_headings,
                                 base_varname + "_detected", 300, justifications);
  comprehensive_list.unprotectContent();

  // Compile the results, with annotation
  SectionContents result;
  result.setTitle("Analysis of hydrogen bonding");
  result.addNarration("Systems with significant hydrogen bonding are presented, with the most "
                      "prevalent cases at the top.  In the following table, columns show:");
  bool pref_match = true;
  if (topl_prefixes.size() > 0) {
    if (coord_prefixes.size() == topl_prefixes.size()) {
      for (size_t i = 0; i < topl_prefixes.size(); i++) {
        pref_match = (pref_match && topl_prefixes[i] == coord_prefixes[i]);
      }
      if (pref_match) {
        result.addNarration("Prefixes of topology and coordinate files are abbreviated:\n" +
                            listCommonPaths(topl_prefixes));
      }
      else {
        result.addNarration("Prefixes of topology files are abbreviated:\n" +
                            listCommonPaths(topl_prefixes));
      }
    }
  }
  if (coord_prefixes.size() > 0 && pref_match == false) {
    result.addNarration("Prefixes of input coordinate files are abbreviated:\n" +
                        listCommonPaths(coord_prefixes));
  }
  OrderedList short_desc(ListEnumeration::NUMBERED);
  short_desc.addItem("Index of the system, within the synthesis");
  short_desc.addItem("Average number of hydrogen bonds found in any frame of the simulation");
  short_desc.addItem("Number of unique donor and acceptor pairs found to have made a hydrogen "
                     "bond at any time during the simulation");
  short_desc.addItem("Number of donor and acceptor pairs found to have made a hydrogen bond at a "
                     "rate greater than " +
                     realToString(minimum_occupancy * 100.0, 2, NumberFormat::STANDARD_REAL) +
                     " %, as may be specified in the &analysis namelist (keyword \"hb_min_occ\")");
  short_desc.addItem("Number of proton donors (e.g. alcohol oxygen) found in the system's Mask I "
                     "atom selection(s)");
  short_desc.addItem("Number of proton acceptors (e.g. carbonyl oxygen) found in the system's "
                     "Mask I atom selection(s)");
  short_desc.addItem("Number of proton donors found in the system's Mask II atom selection(s)");
  short_desc.addItem("Number of proton acceptors found in the system's Mask II atom selection(s)");
  short_desc.addItem("Abbreviated name of the original system topology file, less its suffix and "
                     "possibly with an abbreviated path");
  short_desc.addItem("Abbreviated name of the original system input coordinates file, less its "
                     "suffix and possibly with an abbreviated path");
  result.addList(short_desc);
  result.addTable(short_list);
  result.addNarration("A complete list of hydrogen bonds, out of the matrices of possible bonds, "
                      "found to exist at any time during the simulations.");
  OrderedList comprehensive_desc(ListEnumeration::NUMBERED);
  comprehensive_desc.addItem("Index of the system, within the synthesis");
  comprehensive_desc.addItem("Index of the donor or acceptor heavy atom in Mask I, by the "
                             "system's topological atom order, beginning at zero");
  comprehensive_desc.addItem("Index of the polar hydrogen atom");
  comprehensive_desc.addItem("Index of the donor or acceptor heavy atom in Mask II");
  if (hbr.stat_blocks > 1) {
    comprehensive_desc.addItem("Mean occupancy of the hydrogen bond, with standard deviation "
                               "reported for " + std::to_string(hbr.stat_blocks) + " simulation "
                               "segments");
  }
  else {
    comprehensive_desc.addItem("Mean occupancy of the hydrogen bond through the simulation");
  }
  comprehensive_desc.addItem("Mean donor-acceptor distance, with standard deviation over all "
                             "frames for which the hydrogen bond is found to be formed");
  comprehensive_desc.addItem("Mean donor-hydrogen / acceptor angle, with standard deviation over "
                             "all frames for which the hydrogen bond is found to be formed");
  comprehensive_desc.addItem("Mean donor-acceptor distance, with standard deviation, over all "
                             "frames of the simulation");
  if (hbr.stat_blocks > 1) {
    comprehensive_desc.addItem("Standard deviation of the mean donor-acceptor distance taken "
                               "over " + std::to_string(hbr.stat_blocks) + " simulation segments");
    comprehensive_desc.addItem("Standard deviation of the mean donor-hydrogen / acceptor angle "
                               "taken over " + std::to_string(hbr.stat_blocks) + " simulation "
                               "segments");
  }
  comprehensive_desc.addItem("Name (residue name, residue index within the system starting at "
                             "zero, atom name) of the donor or acceptor in Mask I");
  comprehensive_desc.addItem("Name of the polar hydrogen forming the bond");
  comprehensive_desc.addItem("Name of the atom in Mask II");
  comprehensive_desc.addItem("System label group, as found in the &files namelist (unique labels "
                             "on every system can assist in identifying a particular case)");
  result.addList(comprehensive_desc);
  result.addNarration("The symbols \"--\" and \"::\" between the above columns indicate that the "
                      "heavy atom to the left or right is a donor or acceptor, respectively.");
  result.addTable(comprehensive_list);
  return result;
}

//-------------------------------------------------------------------------------------------------
int HydrogenBondAnalysis::calcStepsPerBlock(const int total_simulation_steps,
                                            const AnalysisControls &trkcon) const {
  const int stepi  = trkcon.getHBondInitiationStep();
  const int nblock = trkcon.getHBondStatBlocks();
  return (total_simulation_steps - stepi + nblock - 1) / nblock;
}
  
//-------------------------------------------------------------------------------------------------
void HydrogenBondAnalysis::allocate() {
  const size_t nsysp  = roundUp<int>(poly_ps_ptr->getSystemCount() + 1, warp_size_int);
  const size_t tp_zu  = roundUp<int>(total_partners, warp_size_int);
  const size_t sbc_zu = statistical_block_count;
  partners.resize(tp_zu);
  int_storage.resize(nsysp + (sbc_zu * tp_zu));
  partner_bounds.setPointer(&int_storage,      0, nsysp);
  bond_formations.setPointer(&int_storage, nsysp, sbc_zu * tp_zu);
  const size_t six_sbc_zu = 6 * sbc_zu;
  const size_t two_zu     = 2;
  const size_t three_zu   = 3;
  const size_t four_zu    = 4;
  const size_t five_zu    = 5;
  dbl_storage.resize(tp_zu * six_sbc_zu);
  distance_accumulators.setPointer(&dbl_storage,                            0, tp_zu);
  distance_sq_accumulators.setPointer(&dbl_storage,            sbc_zu * tp_zu, tp_zu);
  full_distance_acc.setPointer(&dbl_storage,          two_zu * sbc_zu * tp_zu, tp_zu);
  full_distance_sq_acc.setPointer(&dbl_storage,     three_zu * sbc_zu * tp_zu, tp_zu);
  angle_accumulators.setPointer(&dbl_storage,        four_zu * sbc_zu * tp_zu, tp_zu);
  angle_sq_accumulators.setPointer(&dbl_storage,     five_zu * sbc_zu * tp_zu, tp_zu);
}

} // namespace analysis
} // namespace stormm
