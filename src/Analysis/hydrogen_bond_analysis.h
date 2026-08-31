// -*-c++-*-
#ifndef STORMM_HYDROGEN_BOND_ANALYSIS_H
#define STORMM_HYDROGEN_BOND_ANALYSIS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_enumerators.h"
#include "Namelists/nml_analysis.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_report.h"
#include "Reporting/section_contents.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace analysis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using namelist::AnalysisControls;
using namelist::DynamicsControls;
using namelist::MinimizeControls;
using namelist::ReportControls;
using review::SectionContents;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::SynthesisCacheMap;
  
/// \brief The mutable abstract of the HydrogenBondAnalysis class can be fed to kernels or used by
///        CPU processes for rapid access to the original object's data arrays.
struct HBondWriter {

  /// \brief As with other abstracts, the constructor takes a list of arguments for all member
  ///        variables.
  HBondWriter(int n_partners_in, int stat_blocks_in, int init_step_in, int eval_intv_in,
              int block_steps_in, float max_separation_in, float min_angle_in,
              const int4* partners_in, const int* partner_bounds_in, int* tallies_in,
              double* dist_acc_in, double* dist_sqacc_in, double* all_dist_acc_in,
              double* all_dist_sqacc_in, double* angl_acc_in, double* angl_sqacc_in);

  /// \brief As with other abstracts, the default copy and move constructors are valid but the
  ///        presence of const members invalidates the copy and move assignment operators.  The
  ///        const members of the object, once constructed, cannot be reassigned to other values.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object found on the right hand side of an assignment statement
  /// \{
  HBondWriter(const HBondWriter &original) = default;
  HBondWriter(HBondWriter &&original) = default;
  HBondWriter& operator=(const HBondWriter &original) = delete;
  HBondWriter& operator=(HBondWriter &&original) = delete;
  /// \}
  
  const int n_partners;         ///< Total number of possible hydrogen bonding partners to inspect
  const int stat_blocks;        ///< Number of statistical blocks into which to divide the analysis
  const int init_step;          ///< The step number at which analysis is to begin
  const int eval_intv;          ///< Interval, in time steps, at which to evaluate the hydrogen
                                ///<   bond candidates
  const int block_steps;        ///< The number of steps in each statistical block
  const float max_separation;   ///< Maximum donor : acceptor separation for verifying a hydrogen
                                ///<   bond
  const float min_angle;        ///< Minimum donor - hydrogen : acceptor angle for verifying a
                                ///<   hydrogen bond
  const int4* partners;         ///< The list of possible hydrogen bonds to search
  const int* partner_bounds;    ///< Bounds array for each system's portion of the items in
                                ///<   partners
  int* tallies;                 ///< Number of times each hydrogen bond in the list is confirmed.
                                ///<   These tallies may be partitioned for separate blocks of the
                                ///<   simulation.
  double* dist_acc;             ///< Running sum of the donor : acceptor distances for verified
                                ///<   hydrogen bonds
  double* dist_sqacc;           ///< Running sum of the squared donor : acceptor distances for
                                ///<   verified hydrogen bonds
  double* all_dist_acc;         ///< Running sum of the distances throughout the simulation
  double* all_dist_sqacc;       ///< Running sum of the squared distances throughout the simulation
  double* angl_acc;             ///< Running sum of the angles for verified hydrogen bonds
  double* angl_sqacc;           ///< Running sum of the squared angles for verified hydrogen bonds
};

/// \brief The read-only abstract of the HydrogenBondAnalysis classis available for functions that
///        may compile its data for reporting or reference the results in other analyses.  It is
///        produced when the underlying HydrogenBondAnalysis object is const-qualified.
struct HBondReader {

  /// \brief As with other abstracts, the constructor takes a list of arguments for all member
  ///        variables.
  ///
  /// Overloaded:
  ///   - Construct the object with a list of values for all member variables
  ///   - Construct the object based on an equivalent HBondWriter, the writeable abstract of the
  ///     same underlying class object
  /// \{
  HBondReader(int n_partners_in, int stat_blocks_in, int init_step_in, int eval_intv_in,
              int block_steps_in, float max_separation_in, float min_angle_in,
              const int4* partners_in, const int* partner_bounds_in, const int* tallies_in,
              const double* dist_acc_in, const double* dist_sqacc_in,
              const double* all_dist_acc_in, const double* all_dist_sqacc_in,
              const double* angl_acc_in, const double* angl_sqacc_in);

  HBondReader(const HBondWriter &w);

  HBondReader(const HBondWriter *w);
  /// \}

  /// \brief As with other abstracts, the default copy and move constructors are valid but the
  ///        presence of const members invalidates the copy and move assignment operators.  The
  ///        const members of the object, once constructed, cannot be reassigned to other values.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object found on the right hand side of an assignment statement
  /// \{
  HBondReader(const HBondReader &original) = default;
  HBondReader(HBondReader &&original) = default;
  HBondReader& operator=(const HBondReader &original) = delete;
  HBondReader& operator=(HBondReader &&original) = delete;
  /// \}
  
  const int n_partners;          ///< Total number of possible hydrogen bonding partners to inspect
  const int stat_blocks;         ///< Number of statistical blocks into which to divide the
                                 ///<   analysis
  const int init_step;           ///< The step number at which analysis is to begin
  const int eval_intv;           ///< Interval, in time steps, at which to evaluate the hydrogen
                                 ///<   bond candidates
  const int block_steps;         ///< The number of steps in each statistical block
  const float max_separation;    ///< Maximum donor : acceptor separation for verifying a hydrogen
                                 ///<   bond
  const float min_angle;         ///< Minimum donor - hydrogen : acceptor angle for verifying a
                                 ///<   hydrogen bond
  const int4* partners;          ///< The list of possible hydrogen bonds to search
  const int* partner_bounds;     ///< Bounds array for each system's portion of the items in
                                 ///<   partners
  const int* tallies;            ///< Number of times each hydrogen bond in the list is confirmed
  const double* dist_acc;        ///< Running sum of the donor : acceptor distances for verified
                                 ///<   hydrogen bonds
  const double* dist_sqacc;      ///< Running sum of the squared donor : acceptor distances for
                                 ///<   verified hydrogen bonds
  const double* all_dist_acc;    ///< Running sum of the distances throughout the simulation
  const double* all_dist_sqacc;  ///< Running sum of squared distances throughout the simulation
  const double* angl_acc;        ///< Running sum of the angles for verified hydrogen bonds
  const double* angl_sqacc;      ///< Running sum of the squared angles for verified hydrogen bonds
};

/// \brief Track matrices of potential hydrogen bonding partners over the course of a simulation.
///        The potential donors and acceptors are located from within atom selections derived from
///        systems based on atom masks applied to topologies and system labels from user input.
///        Arrays are allocated to store the characteristics and status of all possible hydrogen
///        bonds from among the selections.
class HydrogenBondAnalysis {
public:

  /// \brief Objects of the class must be tailored to a particular collection of systems, with
  ///        atom selections.
  ///
  /// Overloaded:
  ///   - Set up the analysis based on a given simulation length
  ///   - Set up the analysis for a given set of dynamics input
  ///   - Set up the analysis for a given set of energy minimization input
  ///
  /// \param trkcon                  User input data detailing analysis directives
  /// \param total_simulation_steps  The total length of the calculation, used to plan statistical
  ///                                binning
  /// \param dyncon                  User input data detailing molecular dynamics inputs, including
  ///                                the simulation length
  /// \param mincon                  User input data detailing energy minimization inputs,
  ///                                including the number of cycles which can be taken as the
  ///                                simulation length
  /// \{
  HydrogenBondAnalysis(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                       const SynthesisCacheMap &scmap, const AnalysisControls &trkcon,
                       int total_simulation_steps);

  HydrogenBondAnalysis(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                       const SynthesisCacheMap &scmap, const AnalysisControls &trkcon,
                       const DynamicsControls &dyncon);

  HydrogenBondAnalysis(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                       const SynthesisCacheMap &scmap, const AnalysisControls &trkcon,
                       const MinimizeControls &mincon);
  /// \}

  /// \brief Copy and move constructors, as well as copy and move assignment operators, are all
  ///        valid but must be defined as there are POINTER-kind Hybrid objects.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object, found on the right hand side of the assignment statement
  /// \{
  HydrogenBondAnalysis(const HydrogenBondAnalysis &original);
  HydrogenBondAnalysis(HydrogenBondAnalysis &&original);
  HydrogenBondAnalysis& operator=(const HydrogenBondAnalysis &original);
  HydrogenBondAnalysis& operator=(HydrogenBondAnalysis &&original);
  /// \}
  
  /// \brief Get the total number of hydrogen bonds tracked by the object.
  ///
  /// Overloaded:
  ///   - Get the number of possible (though not always formed) hydrogen bonds in a single system
  ///   - Get the total number of possible hydrogen bonds across all systems
  ///
  /// \param system_index  Index of the system of interest
  /// \{
  int getTotalCandidates(int system_index) const;
  int getTotalCandidates() const;
  /// \}

  /// \brief Get the bounds of a series of possible hydrogen bonds for a particular system.
  ///
  /// \param system_index  Index of the system of interest
  int2 getSystemPartnerBounds(int system_index) const;  
  
  /// \brief Get the atom indices of a hydrogen bond within the synthesis.  The synthesis-wide
  ///        topologicalindex of the donor is returned in the "x" member of the tuple while that
  ///        of the proton is returned in the "y" member and the acceptor is returned in the "z"
  ///        member.  The synthesis system index to which the atoms belong is returned in the "w"
  ///        member.
  ///
  /// \param pair_index  The index of the pair, within the entire list
  int4 getPartnerSynthesisIndices(int pair_index) const;

  /// \brief Get the number of blocks involved in block averaging.
  int getBlockAveragingCount() const;

  /// \brief Get the evaluation interval for detecting hydrogen bonds.
  int getEvaluationInterval() const;
  
  /// \brief Get the pointer to the applicable coordinate synthesis.
  const PhaseSpaceSynthesis* getCoordinateSynthesisPointer() const;

  /// \brief Get the pointer to the applicable topology synthesis.
  const AtomGraphSynthesis* getTopologySynthesisPointer() const;

  /// \brief Get the abstract.
  ///
  /// Overloaded:
  ///   - Return an abstract with writeable data arrays from a non-const object
  ///   - Return a read-only abstract from a const-qualified object
  ///
  /// \param tier  Indicate whether to orient the abstract towards data on the CPU host or data on
  ///              the GPU device
  /// \{
  HBondWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  const HBondReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Upload all data from the CPU host to the GPU device.
  void upload();

  /// \brief Download all data from the GPU device to the CPU host.
  void download();

  /// \brief Upload hydrogen bond definitions from the CPU host to the GPU device.
  void uploadDefinitions();

  /// \brief Download hydrogen bond definitions from the GPU device to the CPU host.
  void downloadDefinitions();

  /// \brief Upload accumulator data from the CPU host to the GPU device.
  void uploadAccumulators();

  /// \brief Download accumulator data from the GPU device to the CPU host.
  void downloadAccumulators();
#endif

  /// \brief Set or reset all accumulators to zero.
  ///
  /// \param tier  Indicate whether to initialize accumulators on the CPU host or GPU device
  void initialize(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Produce a report segment based on the hydrogen bonding analysis.  This report will be
  ///        written based on data which is expected to reside, or have been moved to, memory on
  ///        the CPU host.
  ///
  /// \param repcon  User inputs regarding the format of the output 
  SectionContents reportResults(const ReportControls &repcon) const;
  
private:

  int total_partners;           ///< The total number of possible hydrogen bonding pairs (donor,
                                ///<   acceptor) found across all systems
  int statistical_block_count;  ///< The number of blocks to be used in statistical averaging
  int initiation_step;          ///< The first step of the simulation at which to begin analysis
  int evaluation_interval;      ///< Interval, in time steps, at which to evaluate the hydrogen
                                ///<   bond candidates.  Evaluations will take place whenever the
                                ///<   step number is greater than or equal to initiation_step and
                                ///<   a multiple of the evaluation interval.
  int steps_per_block;          ///< The number of calculation time steps to include per block.
                                ///<   This is not often the number of samples in each statistical
                                ///<   block, and some blocks may obtain different numbers of total
                                ///<   samples if the number of available steps does not divide
                                ///<   evenly into the number of blocks.
  double maximum_separation;    ///< The maximum separation between donor and acceptor heavy atoms
                                ///<   at which a hydrogen bond may be considered to exist
  double minimum_angle;         ///< The minimum angle made by donor, attached polar hydrogen, and
                                ///<   acceptor atoms at which a hydrogen bond may be considered to
                                ///<   exist
  double minimum_occupancy;     ///< Minimum occupancy needed to trigger reporting of a specific
                                ///<   hydrogen bonding arrangement
  double minimum_activity;      ///< Minimum presence of hydrogen bonding (the sum of occupancies
                                ///<   in hydrogen bonds for any of the possible partners) in any
                                ///<   given system necessary to report the system as significant
  std::string base_varname;     ///< Base name of the variables to write into a report file when
                                ///<   analysis is complete

  /// Indicators of whether each tracked hydrogen bond has its proton donor in the first group.
  /// The vector element will read FALSE if the acceptor is in the first group and the donor in the
  /// second.  This array is only available on the CPU host as it is expected to be used when
  /// reporting results.
  std::vector<bool> donor_in_first_group;

  /// Counts of the donors and acceptors in each mask, per system, in the associated synthesis.
  /// This information is used when reporting results.
  /// \{
  std::vector<int> mask_i_donor_counts;
  std::vector<int> mask_i_acceptor_counts;
  std::vector<int> mask_ii_donor_counts;
  std::vector<int> mask_ii_acceptor_counts;
  /// \}
  
  /// The primary array of tuples defining which potential hydrogen bonds the object tracks.  The
  /// index of the of proton donor is stored in the "x" member of the tuple, that of the polar
  /// hydrogen (proton) in the "y" member of the tuple, that of the proton acceptor in the "z"
  /// member, and the system index in the "w" member.
  Hybrid<int4> partners;

  // POINTER-kind arrays store the major analtyic data.
  Hybrid<int> partner_bounds;               ///< Bounds array on the possible partners found in
                                            ///<   each system
  Hybrid<int> bond_formations;              ///< Counters for the numbers of times each potential
                                            ///<   arrangement in the partners member array is
                                            ///<   observed to meet the criteria of a hydrogen
                                            ///<   bond.  If block averaging is in effect, counts
                                            ///<   for the kth block will be recorded with an
                                            ///<   offset of k * total_partners.
  Hybrid<double> distance_accumulators;     ///< Running sums of the distance between the donor and
                                            ///<   acceptor particles named in corresponding
                                            ///<   elements of the partners member array, when
                                            ///<   hydrogen bonds are verified.  If block averaging
                                            ///<   is in effect, distances calculated in the kth
                                            ///<   block will be recorded with an offset of
                                            ///<   k * total_partners.
  Hybrid<double> distance_sq_accumulators;  ///< Running sums of the squared distance between donor
                                            ///<   and acceptor particles named in corresponding
                                            ///<   elements of the partners member array.  Because
                                            ///<   the variance is computed for the simulation as a
                                            ///<   whole, block averaging does not record separate
                                            ///<   squared distance accumulators.
  Hybrid<double> full_distance_acc;         ///< Running sums of the distance between the donor and
                                            ///<   acceptor particles named in corresponding
                                            ///<   elements of the partners member array,
                                            ///<   throughout the entire simulation
  Hybrid<double> full_distance_sq_acc;      ///< Running sums of the squared distance between the
                                            ///<   donor and acceptor particles, throughout the
                                            ///<   entire simulation
  Hybrid<double> angle_accumulators;        ///< Running sums of the donor-H :: acceptor angles
                                            ///<   made by particles named in corresponding
                                            ///<   elements of the partners member array
  Hybrid<double> angle_sq_accumulators;     ///< Running sums of the squared angles made by donors,
                                            ///<   protons, and acceptors named in corresponding
                                            ///<   elements of the partners member array

  /// The expected numbers of samples in each block, calculated when the object is first created.
  std::vector<double> block_sample_counts;
  
  // Pointers to the objects that this analysis supports
  PhaseSpaceSynthesis *poly_ps_ptr;  ///< Pointer to the coordinate synthesis
  AtomGraphSynthesis *poly_ag_ptr;   ///< Pointer to the topology synthesis
  SynthesisCacheMap *scmap_ptr;      ///< Pointer to the synthesis cache map, tying results back
                                     ///<   to user input structures

  /// Storage arrays targeted by various POINTER_kind Hybrid objects
  Hybrid<int> int_storage;     ///< Integer data (used to accumulate tallies as well as overflow to
                               ///<   fixed-precision representations of real-valued running sums)
  Hybrid<double> dbl_storage;  ///< Long long integer data (used to accumulate real-valued running
                               ///<   sums)

  /// \brief Calculate the maximum number of simulation time steps per statistical block used in
  ///        block averaging.  This number, while it does not give a definitive number of samples
  ///        that may be present in each block, informs which block data from any given segment of
  ///        the simulation will contribute to.
  ///
  /// \param total_simulation_steps  The total number of steps in the simulation, nstlim or maxcyc
  ///                                in AMBER input terminology
  /// \param trkcon                  User input containing analysis directives
  int calcStepsPerBlock(const int total_simulation_steps, const AnalysisControls &trkcon) const;
  
  /// \brief Allocate data for the object, based on a calculated number of possible hydrogen bonds
  ///        stored in the total_partners member variable and the number of statistical blocks
  ///        stored in statistical_block_count.
  void allocate();
};

} // namespace analysis
} // namespace stormm

#endif
