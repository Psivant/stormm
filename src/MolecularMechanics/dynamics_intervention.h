// -*-c++-*-
#ifndef STORMM_DYNAMICS_INTERVENTION_H
#define STORMM_DYNAMICS_INTERVENTION_H

#include <string>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/gpu_enumerators.h"
#include "Analysis/hydrogen_bond_analysis.h"
#include "Chemistry/atommask.h"
#include "Debug/watcher.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_util.h"
#include "Math/rounding.h"
#include "Namelists/nml_analysis.h"
#include "Namelists/nml_debug.h"
#include "Namelists/nml_dynamics.h"
#include "Potential/cellgrid.h"
#include "Potential/energy_enumerators.h"
#include "Potential/local_exclusionmask.h"
#include "Potential/scorecard.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/condensate.h"
#include "Synthesis/implicit_solvent_workspace.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "mm_enumerators.h"

namespace stormm {
namespace mm {

using analysis::HydrogenBondAnalysis;
using analysis::HBondWriter;
using card::GpuDetails;
using card::Hybrid;
using card::HybridFormat;
using card::HybridTargetLevel;
using chemistry::AtomMask;
using data_types::getStormmScalarTypeName;
using data_types::getHpcVectorTypeName;
using diskutil::getBaseName;
using energy::CellGrid;
using energy::CellGridWriter;
using energy::LocalExclusionMask;
using energy::LocalExclusionMaskReader;
using energy::NeighborListKind;
using energy::NonbondedTheme;
using energy::restoreType;
using energy::ScoreCard;
using energy::ScoreCardWriter;
using namelist::AnalysisControls;
using namelist::DebugControls;
using namelist::DynamicsControls;
using debug::Watcher;
using debug::WatcherWriter;
using stmath::roundUp;
using structure::ApplyConstraints;
using synthesis::AtomGraphSynthesis;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::CondensateWriter;
using synthesis::ImplicitSolventWorkspace;
using synthesis::ISWorkspaceKit;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using synthesis::SeMaskSynthesisReader;
using synthesis::StaticExclusionMaskSynthesis;
using synthesis::SyAtomUpdateKit;
using synthesis::SyNonbondedKit;
using synthesis::SyRestraintKit;
using synthesis::SyValenceKit;
using topology::UnitCellType;
using trajectory::CoordinateCycle;
using trajectory::getNextCyclePosition;
using trajectory::IntegrationStage;
using DynIntvFuncPtr = void(*)(int, int, HybridTargetLevel);

/// The maximum number of unique intervals that will be considered when determining whether an
/// active intervention is to take place on a particular simulation time step is put in place to
/// guard against making an excessive number of modulo operations (x % y) twice.  If there are more
/// than this number of possible intervention intervals, they will all be replaced with a single
/// interval, 1, which will trigger checks through each attached function on every step.
///
/// The process is to ask, "on this step, is there any intervention that should occur?" and then,
/// if so, "which interventions out of the list will occur?" Whenever a step is found to have some
/// active intervention, each attached function will be queried to determine whether it should fire
/// off.  Therefore, modulo operations will be conducted at the beginning of every function in
/// addition to whatever modulo operations were conducted over the list of relevant intervals to
/// determine whether any attached function might fire off.
constexpr int max_relevant_intervals = 16;
  
/// \brief One global instance of this class will be created at the outset of a program using it,
///        so that any and all information which might be needed by a method which intervenes in
///        the main dynamics loop can reference the information without a need to pass in formal
///        arguments.  This permits such functions to share a common profile and be passed to the
///        dynamics routines as function pointers.  This class carries copies of many abstracts to
///        be dispensed as needed, but has no abstract of its own.
class DynamicsIntervention {
public:

  /// \brief Only the basic constructor, creating an empty object, will be available.  The object
  ///        is expected to be populated once all systems data and other input parameters have been
  ///        read.
  DynamicsIntervention();

  /// \brief Movement and copying of DynamicsIntervention class objects will be forbidden.  The
  ///        one global instance of the class is what functions that intervene in dynamics are
  ///        expected to reference, although others can be created for specific developmental
  ///        purposes.
  /// \{
  DynamicsIntervention(const DynamicsIntervention &original) = delete;
  DynamicsIntervention(DynamicsIntervention &&original) = delete;
  DynamicsIntervention& operator=(const DynamicsIntervention &original) = delete;
  DynamicsIntervention& operator=(DynamicsIntervention &&original) = delete;
  /// \}

  /// \brief Get a signal as to whether interventions are active, at all.
  bool isActive() const;

  /// \brief Get a signal as to whether interventions are active on a particular MD step.  Most
  ///        debugging interventions, for example, will always be active, but other activities
  ///        like simulation analysis will only be active on selected steps.
  ///
  /// \param step_number  The number of the step in question
  bool isActiveOnStep(int step_number) const;

  /// \brief Return a const reference to the vector of active intervals, for inspection or testing
  ///        purposes.
  const std::vector<int>& getActiveIntervals() const;
  
  /// \brief Get the specifications for the GPU available to the calculation.
  const GpuDetails& getGpuDetails() const;

  /// \brief Get the number of streaming multiprocessors on the GPU.
  int getSMPCount() const;

  /// \brief Indicate whether the official problem set's coordinate synthesis is attached to the
  ///        object.
  bool hasPhaseSpaceSynthesis() const;

  /// \brief Indicate whether there is a Condensate attached to the object, holding reserved
  ///        floating-point memory for the coordinate synthesis.
  bool hasCondensate() const;

  /// \brief Indicate whether there is a topology synthesis attached to the object.
  bool hasTopologySynthesis() const;
  
  /// \brief Indicate whether the object has an ImplicitSolventWorkspace attached to it.
  bool hasImplicitSolventWorkspace() const;

  /// \brief Indicate whether the appropriate set of exclusion masks has been attached to the
  ///        object, if a topology synthesis for the set of problems is also attached.  Otherwise,
  ///        indicate whether any set of exclusion masks has been attached.
  bool hasExclusionMasks() const;

  /// \brief Indicate whether an energy tracking object for the set of problems has been attached
  ///        to the DynamicsIntervention object.
  bool hasEnergyTracking() const;

  /// \brief Indicate whether an anomaly reporting object for tracking large forces, failed
  ///        constraint convergences, and other events has been attached.
  bool hasAnomalyReporting() const;

  // The following accessors indicate that specific analyses or debugging operations are in effect,
  // which may require additional attachments or preparations to be made deeper in the program
  // once these objects are available with their size parameters known.
  bool doNeighborListChecks() const;  ///< Check the simulation neighbor list and trigger creation
                                      ///<   of a neighbor list workspace to support the process
  
  /// \brief Get an abstract to the main coordinate synthesis, the split fixed-precision
  ///        PhaseSpaceSynthesis.
  ///
  /// Overloaded:
  ///   - Produce the writeable abstract of the coordinate synthesis from a non-const
  ///     DynamicsIntervention object
  ///   - Produce the read-only abstract of the coordinate synthesis from a const-qualified
  ///     DynamicsIntervention object
  /// 
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  /// \{
  PsSynthesisWriter getCoordinateData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  const PsSynthesisReader
  getCoordinateData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get an abstract to the main coordinate synthesis, the split fixed-precision
  ///        PhaseSpaceSynthesis.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const PsSynthesisReader
  getReadOnlyCoordinateData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get an abstract to the auxiliary coordinate synthesis, the real-valued Condensate.
  ///        Overloading and descriptions of input parameters follow from getCoordinateData(),
  ///        above.
  /// \{
  CondensateWriter getCondensateData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  const CondensateReader getCondensateData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}
  
  /// \brief Get the double-precision valence parameter kit, spanning all force constants and terms
  ///        for all systems in the synthesis.  The abstract will be returned for data on the CPU
  ///        host or GPU device according to the way it was created.  If no such abstract is
  ///        available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyValenceKit<double>
  getDoublePrecisionValenceKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision valence parameter kit.  The abstract will be returned for
  ///        data on the CPU host or GPU device according to the way it was created.  If no such
  ///        abstract is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyValenceKit<float>
  getSinglePrecisionValenceKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the double-precision non-bonded parameter kit.  The abstract will be returned for
  ///        data on the CPU host or GPU device according to the way it was created.  If no such
  ///        abstract is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyNonbondedKit<double, double2>
  getDoublePrecisionNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision non-bonded parameter kit.  The abstract will be returned for
  ///        data on the CPU host or GPU device according to the way it was created.  If no such
  ///        abstract is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyNonbondedKit<float, float2>
  getSinglePrecisionNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get the double-precision atom update kit.  The abstract will be returned for data on
  ///        the CPU host or GPU device according to the way it was created.  If no such abstract
  ///        is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyRestraintKit<double, double2, double4_16a>
  getDoublePrecisionRestraintKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision atom update kit.  The abstract will be returned for data on
  ///        the CPU host or GPU device according to the way it was created.  If no such abstract
  ///        is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyRestraintKit<float, float2, float4>
  getSinglePrecisionRestraintKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the double-precision atom update kit.  The abstract will be returned for data on
  ///        the CPU host or GPU device according to the way it was created.  If no such abstract
  ///        is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyAtomUpdateKit<double, double2, double4_16a>
  getDoublePrecisionAtomUpdateKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision atom update kit.  The abstract will be returned for data on
  ///        the CPU host or GPU device according to the way it was created.  If no such abstract
  ///        is available, this will raise an exception.
  ///
  /// \param tier  Indicate whether to get an abstract with pointers valid for data on the CPU host
  ///              or GPU device
  const SyAtomUpdateKit<float, float2, float4>
  getSinglePrecisionAtomUpdateKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get an abstract to the attached CellGrid object, with its template parameters
  ///        reconstituted.  The developer must specify the return type, but the function itself
  ///        will check that specification against the details of the attached CellGrid.
  ///
  /// Overloaded:
  ///   - Provide the non-bonded theme, in addition to the orientation and tier
  ///   - Do not provide the non-bonded theme, relying on internal checks to confirm that all
  ///     options have the same theme which will then be automatically selected
  ///
  /// \param theme        Selected theme for which to obtain the abstract
  /// \param orientation  Selected point in the coordinate cycle for which to obtain the abstract
  /// \param tier         Indicate whether to obtain the abstract on the CPU host or GPU device
  /// \{
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4>
  getNeighborListData(NonbondedTheme theme, CoordinateCycle orientation = CoordinateCycle::WHITE,
                      HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4>
  getNeighborListData(CoordinateCycle orientation = CoordinateCycle::WHITE,
                      HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get a double-precision implicit solvent workspace abstract.
  ///
  /// \param orientation  Indicate the point in the WHITE / BLACK time cycle for which to obtain
  ///                     the abstract
  /// \param tier         Indicate whether to obtain the abstract for the CPU host or GPU device
  ISWorkspaceKit<double>
  getDoublePrecisionISWorkspace(CoordinateCycle orientation,
                                HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Get a single-precision implicit solvent workspace abstract.
  ///
  /// \param orientation  Indicate the point in the WHITE / BLACK time cycle for which to obtain
  ///                     the abstract
  /// \param tier         Indicate whether to obtain the abstract for the CPU host or GPU device
  ISWorkspaceKit<float>
  getSinglePrecisionISWorkspace(CoordinateCycle orientation,
                                HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Get a local exclusion mask abstract.
  ///
  /// \param tier  Indicate whether to obtain the abstract for the CPU host or GPU device
  const LocalExclusionMaskReader
  getLocalExclusionMaskData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a static exclusion mask abstract.
  ///
  /// \param tier  Indicate whether to obtain the abstract for the CPU host or GPU device
  const SeMaskSynthesisReader
  getStaticExclusionMaskData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a coded signal for the coordinate data type of the attached neighbor list (and
  ///        that of any attached neighbor list workspaces).
  size_t getNeighborListCoordinateType() const;

  /// \brief Get a coded signal for the accumulation data type of the attached neighbor list (and
  ///        that of any attached neighbor list workspaces).
  size_t getNeighborListAccumulationType() const;

  /// \brief Get a coded signal for the calculation data type of the attached neighbor list (and
  ///        that of any attached neighbor list workspaces).
  size_t getNeighborListCalculationType() const;

  /// \brief Get a coded signal for the coordinate and property tuple data type of the attached
  ///        neighbor list (and that of any attached neighbor list workspaces).
  size_t getNeighborListTupleType() const;

  /// \brief Get the layout of the simulation neighbor list, as noted when the relevant objects
  ///        were attached.
  NeighborListKind getNeighborListLayout() const;
  
  /// \brief Get an energy tracking abstract.
  ///
  /// \param tier  Indicate whether to obtain the abstract on the CPU host or GPU device
  ScoreCardWriter getEnergyTrackingData(HybridTargetLevel tier = HybridTargetLevel::HOST);  

  /// \brief Get an abstract for the attached Watcher class object, for reporting anomalous forces
  ///        acting on particles, failed constraint convergences, and other notable events in the
  ///        simulation.
  ///
  /// \param tier  Indicate whether to obtain the abstract on the CPU host or GPU device
  WatcherWriter getAnomalyReportingData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Get an abstract for the attached HydrogenBondAnalysis object, for tracking the
  ///        occupancy of possible hydrogen bonds among all systems in the calculation.
  ///
  /// \param index  Provide the index of a particular analysis
  /// \param tier   Indicate whether to obtain the abstract on the CPU host or GPU device
  HBondWriter getHydrogenBondAnalysisData(int index = 0,
                                          HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Get a pointer to the coordinate PhaseSpaceSynthesis associated with the object.
  ///
  /// Overloaded:
  ///   - Get a mutable pointer from a non-cnost DynamicsIntervention object
  ///   - Get a const-qualified pointer if the DynamicsIntervention itself is immutable
  /// \{
  const PhaseSpaceSynthesis* getPhaseSpaceSynthesisPointer() const;
  PhaseSpaceSynthesis* getPhaseSpaceSynthesisPointer();
  /// \}

  /// \brief Get a pointer to a coordinate workspace, an additonal PhaseSpaceSynthesis object,
  ///        associated with the DynamicsIntervention object under certain circumstances such as
  ///        Debugging operations.  This will be a non-const pointer so that the workspace can be
  ///        modified by any function using it.
  ///
  /// \param index  The index of the workspace of interest (likely only a single alternate copy of
  ///               the coordinates, velocities, and forces is available, though it will also have
  ///               WHITE and BLACK cycle points within it)
  PhaseSpaceSynthesis* getWorkspacePointer(size_t index = 0);

  /// \brief Return a pointer to the attached neighbor list object, as is in use by the simulation.
  ///        The requested templated types will be checked against the pointer which was actually
  ///        stored.
  ///
  /// \param theme  Indicate the theme of the neighbor list of interest
  /// \{
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>* getNeighborListPointer(NonbondedTheme theme);

  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>* getNeighborListPointer(NonbondedTheme theme) const;
  /// \}
  
  /// \brief Get a pointer to the neighbor list workspace, an additional CellGrid object templated
  ///        in the same manner as the main problem set and associated with the
  ///        DynamicsIntervention object for debugging purposes.  As with getWorkspacePointer, the
  ///        function and result are non-const so that the workspace may be modified by any
  ///        function using it.  Overloads and descriptions of input parameters follow from the
  ///        member function getWorkspacePointer(), above.
  /// \{
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>* getNLWorkspacePointer(NonbondedTheme theme);

  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>* getNLWorkspacePointer(NonbondedTheme theme) const;
  /// \}
  
  /// \brief Get a pointer to the coordinate Condensate associated with the object.  This will be
  ///        a non-const pointer so that thte coordiante object behind it can be manipulated.
  Condensate* getCondensatePointer();

  /// \brief Get a pointer to the topology synthesis associated with the object.  This will be a
  ///        non-const pointer so that the topology itself may be manipulated during dynamics
  ///        (i.e. a polarizable force field).
  AtomGraphSynthesis* getTopologySynthesisPointer();

  /// \brief Get a pointer to the implicit solvent workspace associated with the object.  Like
  ///        workspaces for the coordinates, the attached implicit solvent workspace may be (and,
  ///        is likely to be) a different object that the one used in the main force calculations.
  ImplicitSolventWorkspace* getImplicitSolventWorkspacePointer();
  
  /// \brief Get a pointer to the set of exclusion masks attached to the object.
  const LocalExclusionMask* getLocalExclusionMaskPointer() const;

  /// \brief Get a pointer to the set of exclusion masks attached to the object.
  const StaticExclusionMaskSynthesis* getStaticExclusionMaskPointer() const;

  /// \brief Get a mutable pointer to the energy tracking object used in interventions.
  ///
  /// Overloaded:
  ///   - Produce a mutable pointer to the energy tracking object attached to a non-const
  ///     DynamicsIntervention object
  ///   - Produce a const-qualified pointer to the energy tracking object attached to a
  ///     const-qualified DynamicsIntervention object
  /// \{
  ScoreCard* getEnergyTrackingPointer();
  const ScoreCard* getEnergyTrackingPointer() const;
  /// \}

  /// \brief Get a mutable pointer to the anomaly reporting object used during interventions.
  ///        Overloading in this function is similar to getEnergyTrackingPointer(), above.
  /// \{
  Watcher* getAnomalyReportingPointer();
  const Watcher* getAnomalyReportingPointer() const;
  /// \}

  /// \brief Get the number of attached hydrogen bond analyses.
  int getHydrogenBondAnalysisCount() const;

  /// \brief Get a mutable pointer to an attached hydrogen bonding analysis.
  ///
  /// \param index  Index of a particular analysis of interest, if more than one has been attached
  HydrogenBondAnalysis* getHydrogenBondAnalysisPointer(int index = 0);

#ifdef STORMM_USE_HPC
  /// \brief Upload all debugging-specific data.  This will collect processes for each attached
  ///        object related to debugging.
  void uploadDebugging();

  /// \brief Download all debugging-specific data.  This will collect processes for each attached
  ///        object related to debugging.
  void downloadDebugging();
  
  /// \brief Upload all analysis-specific data.  This subsumes processes found in
  ///        uploadAnalysisSetup() and uploadAnalysisData(), below, either of which could be more
  ///        efficient and precise for a particular step in a workflow.
  void uploadAnalysis();

  /// \brief Download all analysis-specific data.  This subsumes processes found in
  ///        downloadAnalysisSetup() and downloadAnalysisData(), below, either of which could be
  ///        more efficient and precise for a particular step in a workflow.
  void downloadAnalysis();
  
  /// \brief Upload analysis-specific setup parameters in ALL of the analysis-related attachments.
  ///        This can be used a general preparatory step, prior to beginning simulations in a
  ///        STORMM application.
  void uploadAnalysisSetup();
  
  /// \brief Download analysis-specific setup parameters in ALL of the analysis-related
  ///        attachments.
  void downloadAnalysisSetup();

  /// \brief Upload analysis-specific data in ALL of the analysis-related attachments.
  void uploadAnalysisData();
  
  /// \brief Download analysis-specific data in ALL of the analysis-related attachments.
  void downloadAnalysisData();
#endif

  /// \brief Make a record of the GPU assigned to carry out calculations.
  ///
  /// \param gpu_in  Details of the available GPU
  void setGpu(const GpuDetails &gpu_in);
  
  /// \brief Attach a coordinate synthesis which interventions will take place upon.
  ///
  /// Overloaded:
  ///   - Attach a PhaseSpaceSynthesis, the "main" coordinate synthesis
  ///   - Attach a Condensate, the "auxiliary" or "staging" coordinate synthesis
  ///
  /// \param poly_ps_in  The main coordinate synthesis
  /// \param cdns_in     The auxiliary coordinate synthesis,
  /// \{
  void setCoordinateSynthesis(PhaseSpaceSynthesis *poly_ps_in);
  void setCoordinateSynthesis(Condensate *cdns_in);
  /// \}
  
  /// \brief Coordinate data for the entire synthesis is passed to an intervention function as
  ///        one of obligatory arguments.  A call to this function will create new topology
  ///        abstracts, in all available precision models, to afford such information when needed.
  ///        Abstracts will be created and stored for all available tiers of the topology
  ///        synthesis.
  ///
  /// \param poly_ag_in  The topology synthesis from which to take abstracts
  void setTopologySynthesis(AtomGraphSynthesis *poly_ag_in);

  /// \brief While the coordinate source of truth in any simulation is the PhaseSpaceSynthesis,
  ///        containing global positions and the highest level of numerical precision, neighbor
  ///        lists may or may not be employed as part of the simulation.  The neighbor list class
  ///        is templated, and as such a pointer or reference to the original class object would
  ///        create unwanted complexity in DynamicsIntervention class objects.  However, the
  ///        neighbor list's template-free abstracts may be stored without forcing templated
  ///        charadteritics onto DynamicsIntervention.  The neighbor list has its own time cycle,
  ///        which will be checked for alignment to the underlying coordinate synthesis.
  ///
  /// \param cg_a  The first of up to two neighbor list objects to reference during interventions
  /// \param cg_b  Second of two neighbor list objects to reference in subsequent interventions
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  void setNeighborList(CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_a,
                       CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_b = nullptr);

  /// \brief Attach in implicit solvent workspace to the object appropriate for the synthesis of
  ///        systems in the current problem set.  Abstracts for the workspace will be created in
  ///        either precision mode by which it can be utilized, and dispensed like other abstracts
  ///        by the DynamicsIntervention class object.  Because this workspace will be used for
  ///        calculations initiated by intervention functions, it is often preferable to attach a
  ///        workspa
  ///
  /// \param isw  The implicit solvent workspace to attach
  void setImplicitSolventWorkspace(ImplicitSolventWorkspace *isw);
  
  /// \brief Like the underlying topology data, pre-computed exclusion masks will be available as
  ///        an array of abstracts at all relevant memory levels.
  ///
  /// Overloaded:
  ///   - Load exclusion masks for periodic systems
  ///   - Load exclusion masks for non-periodic systems
  ///
  /// \param lem  The pre-established set of exclusion masks for periodic systems
  /// \param se   The pre-established set of exclusion masks for non-periodic systems
  /// \{
  void setExclusionMasks(const LocalExclusionMask &lem);
  void setExclusionMasks(const StaticExclusionMaskSynthesis &se);
  /// \}

  /// \brief Set the energy tracking object to be accessible by interventions.
  ///
  /// \param sc_in  The pre-established energy tracking object
  void setEnergyTracking(ScoreCard *sc_in);

  /// \brief Set the anomaly reporting hub to be accessible by interventions.
  ///
  /// \param anom_in  The pre-established anomaly reporting object
  void setAnomalyReporting(Watcher *anom_in);

  /// \brief Set a workspace or workspaces for calculations involving the neighbor list.  Only one
  ///        neighbor list, or a pair of them, will be kept by the object and all analyses or
  ///        debugging operations involving the workspace will work out of index 0.
  ///
  /// \param cg_a  The first neighbor list to attach
  /// \param cg_b  The second neighbor list to attach.  The default value will indicate that no
  ///              secondary neighbor list is present.
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  void setNLWorkspace(CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_a,
                      CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_b = nullptr);

  /// \brief Attach the recording objects and their associated abstract for various analyses.
  ///        Overloads are available for each unique process.
  ///
  /// \param attachment  A mutable pointer to the prepared object which will govern / record the
  ///                    analysis of interest
  /// \param dyncon      User input for molecular dynamics controls, contains information such as
  ///                    whether geometric constraints are in effect
  /// \param tier        Indicate whether calcualtions will run on CPU host or GPU device resources
  /// \{
  void setAnalysis(HydrogenBondAnalysis *attachment, const DynamicsControls &dyncon,
                   HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Following the setter functions above for assigning simulation resources and storing
  ///        their abstracts in a DynamicsIntervention object, this function will 
  ///        debugging operations.
  ///
  /// Overloaded:
  ///   - Provide an arbitrary state of constraint applications
  ///   - Provide the &dynamics namelist control data from which to extract the state of constraint
  ///     applications during the calculation.
  ///
  /// \param anom_in       Modifiable pointer to the object created to catalog large (anomalous)
  ///                      forces, high velocities, and other irregularities
  /// \param dbgcon        User input focused on debugging operations
  /// \param workspace_in  Designate a (new) coordinate synthesis which can be used to hold scratch
  ///                      work as the CPU carries out calculations.  Coordinates will nbe copied
  ///                      into this space and forces recomputed by the CPU stored in it.  The
  ///                      workspace must have a memory layout with allocations on the CPU host.
  /// \param tier          Indicate whether debugging is to focus on calculations which were done
  ///                      on the CPU host or GPU device.  Most, if not all, debugging checks will
  ///                      be performed with CPU resources.
  /// \{
  void setDebugging(Watcher *anom_in, const DebugControls &dbgcon,
                    ApplyConstraints enforce_constraints, PhaseSpaceSynthesis *workspace_in,
                    HybridTargetLevel tier = HybridTargetLevel::HOST);

  void setDebugging(Watcher *anom_in, const DebugControls &dbgcon, const DynamicsControls &dyncon,
                    PhaseSpaceSynthesis *workspace_in,
                    HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Add an action to be taken by the dynamics intervention.  Actions will be taken by
  ///        the object, at such intervals as they are requested, based on the order in which they
  ///        are added to the object.  
  ///
  /// Overloaded:
  ///   - Provide an execution interval
  ///   - Indicate only the compute resource which will execute the action.  Omitting an execution
  ///     frequency implies that the action will occur every step.
  ///
  /// \param action            Function pointer for the action to take
  /// \param process           Stage of the dynamics cycle at which this action shall be taken
  /// \param interval          Optional interval with which to execute the action.  The default
  ///                          behavior is to carry out the action at every step.
  /// \param attachment_index  Index of the attached resource to draw upon from the implied array
  ///                          of such resources' abstracts
  /// \param tier              Indicate whether the CPU host or the GPU device performs the action
  /// \{
  void addAction(DynIntvFuncPtr action, IntegrationStage process, int interval = 1,
                 int attachment_index = 0, HybridTargetLevel tier = HybridTargetLevel::HOST);

  void addAction(DynIntvFuncPtr action, IntegrationStage process, HybridTargetLevel tier);
  /// \}

  /// \brief Execute the series of interventions queued in the object, based on the attached
  ///        resources.
  ///
  /// \param step     Current step of the calculation (e.g. molecular dynamics or energy
  ///                 minimization)
  /// \param process  The stage of the molecular dynamics time cycle at which the execution is
  ///                 occurring.  The execute function will call one of a series of private
  ///                 functions where the actions appropriate to each stage are stored.
  void execute(int step, IntegrationStage process);
  
#ifdef STORMM_USE_HPC
  /// \brief Upload all data to the GPU device
  void upload();
  
  /// \brief Download all data from the GPU device
  void download();
#endif

private:

  // The DynamicsIntervention object records what it will do, including debugging activities.
  bool active_interventions;    ///< Indicate that any kind of intervention is to be performed

  // Some actions can trigger additional layers of debugging or perhaps analysis.  These flags may
  // be carried over from various debugging and analysis namelists, for the sake of having the
  // information ready to dispense directly from the DynamicsIntervention object, in particular its
  // global instance.
  bool investigate_large_forces;  ///< Investigate large forces with an additional CPU-based
                                  ///<   evaluation, based on coordinates copied from the current
                                  ///<   state of the simulation.
  bool check_neighbor_list;       ///< Check the contents of the neighbor list with additional
                                  ///<   CPU-based processes, running independently of the original
                                  ///<   simulation even if that, too, utilizes CPU resources
  
  // The object will indicate whether interventions are active for a given molecular simulation
  // step number.  This is evaluated based on a series of relevant step intervals.
  int relevant_interval_count;          ///< Trusted length of the relevant_intervals array,
                                        ///<   stored separately for speed of access
  std::vector<int> relevant_intervals;  ///< The object keeps a list of intervals over which its
                                        ///<   activities may need to be performed.  If a large
                                        ///<   number of different frequencies are relevant, the
                                        ///<   object will revert to a state in which it attempts
                                        ///<   to execute interventions at every step and
                                        ///<   relevant_interval_count is set to 1.  On any step
                                        ///<   in which interventions may occur, activity-specific
                                        ///<   intervals are then queried.
  
  /// Record details of the available GPU.
  GpuDetails gpu;
  
  /// Indicate whether topological abstracts are taken for data on the CPU host or GPU device.
  /// Each element of these arrays corresponds to an element of one of the {d,f}poly_??s array
  /// member variable, below.  When requesting an abstract to memory at a given tier, these arrays
  /// are searched for the correct match.
  std::vector<HybridTargetLevel> topology_abstract_tiers;

  /// Arrays holding various topology abstracts.  Pointer member variable above are set to the
  /// first element of the corresponding array, and whenever abstracts are assigned the arrays are
  /// first wiped by resizing them to zero elements so that the new abtracts may be added via the
  /// push_back method.  This overcomes problems associated with adding abstracts to the
  /// DynamicsIntervention object after its initial creation, as most abstracts have const members.
  /// The arrays are otherwise expected to hold one element, each.
  /// \{
  std::vector<SyValenceKit<double>> dpoly_vks;
  std::vector<SyValenceKit<float>> fpoly_vks;
  std::vector<SyNonbondedKit<double, double2>> dpoly_nbks;
  std::vector<SyNonbondedKit<float, float2>> fpoly_nbks;
  std::vector<SyRestraintKit<double, double2, double4_16a>> dpoly_rks;
  std::vector<SyRestraintKit<float, float2, float4>> fpoly_rks;
  std::vector<SyAtomUpdateKit<double, double2, double4_16a>> dpoly_auks;
  std::vector<SyAtomUpdateKit<float, float2, float4>> fpoly_auks;
  /// \}

  /// The coordinate abstracts will target data on the CPU host or on the GPU device.  Their
  /// associations are noted in this array.
  std::vector<HybridTargetLevel> coordinate_abstract_tiers;

  /// Abstracts of a PhaseSpaceSynthesis are also associated with one stage of the coordinate time
  /// cycle.  Those associations are provided in this array.  If a neighbor list is also linked to
  /// the DynamicsIntervention object, the cycle stages of CellGrid abstracts should align with
  /// these coordinate abstracts.  There can be two neighbor lists for every coordinate synthesis,
  /// however, so the nieghbor lists get separate arrays to catalog their associations.
  std::vector<CoordinateCycle> coordinate_cycle_stages;
  
  /// The main coordinate abstracts are held in similar arrays, to be dispensed to any function
  /// that needs them.  Both mutalbe and read-only abstracts are kept on hand.
  /// \{
  std::vector<PsSynthesisWriter> poly_psws;
  std::vector<PsSynthesisReader> poly_psrs;
  /// \}
  
  /// The condensate object's abstracts likewise point to data on the CPU host or the GPU device.
  std::vector<HybridTargetLevel> condensate_abstract_tiers;
  
  /// A condensed format coordinate object's abstract can be likewise stored, ready to perform
  /// interventions on the dynamics.
  std::vector<CondensateWriter> cdnsws;

  // The neighbor list depends on four data types, but by stripping out templated character in
  // the abstracts the nature of these types would be lost.  Record them here so that the
  // abstracts may be restored to their proper templated type.
  size_t ngbr_list_tcoord;   ///< Data type for local coordinates and cell transformation matrices
                             ///<   in the neighbor list
  size_t ngbr_list_tacc;     ///< Data type for force accumulations in the neighbor list
  size_t ngbr_list_tcalc;    ///< Data type for calculations involving the neighbor list
  size_t ngbr_list_tcoord4;  ///< Data type for coordinate / property tuple storage in the
                             ///<   neighbor list images

  /// Take note of the layout of the neighbor list.  It may be two separate neighbor lists for
  /// electrostatic and van-der Waals interactions, which would, in turn, have consequences for the
  /// layout of the member array variables neighbor_list_themes and voided_cgws, below.
  NeighborListKind nl_layout;
  
  /// This array marks whether the neighbor list abstracts in each corresponding element of
  /// voided_cgs (below) are pointed to data on the CPU host or the GPU device.
  std::vector<HybridTargetLevel> neighbor_list_abstract_tiers;

  /// This array records the point in the coordinate time cycle at which neighbor list abstracts
  /// are taken, crucial for understanding which abstract to use at a given time step.  The
  /// current point in the time cycle can be obtained from the pointer to the PhaseSpaceSynthesis,
  /// which will be checked for alignment to the underlying CellGrid object(s) as the neighbor
  /// list abstracts are loaded.
  std::vector<CoordinateCycle> neighbor_list_cycle_stages;

  /// Neighbor list themes can include "electrostatics," "van der Waals", or "both." This array
  /// records the nature of each neighbor list abstract so that the correct one can be retrieved.
  std::vector<NonbondedTheme> neighbor_list_themes;
  
  /// An array of neighbor list abstracts, all templateless.  The associated size_t codes for the
  /// various templated types can be used to restore the abstracts to their functional form.
  std::vector<CellGridWriter<void, void, void, void>> voided_cgws;

  /// This array indicates whether abstracts of an attached implicit solvent workspace are pointed
  /// to data on the CPU host or on the GPU device.  The array applies to both single- and double-
  /// precision abstracts of the scratch work object.
  std::vector<HybridTargetLevel> implicit_solvent_abstract_tiers;

  /// Like the main coordinate synthesis and neighbor list objects, the implicit solvent workspace
  /// has separate memory allocations for the WHITE and BLACK stages of the dynamics time cycle.
  std::vector<CoordinateCycle> implicit_solvent_cycle_stages;

  /// There are two precision modes for the implicit solvent workspace, with abstracts create for
  /// either of them.
  /// \{
  std::vector<ISWorkspaceKit<double>> d_iswks;
  std::vector<ISWorkspaceKit<float>> f_iswks;
  /// \}
  
  /// Indicators of whether exclusion masks are stored at the level of the CPU host or GPU device
  std::vector<HybridTargetLevel> exclusion_mask_abstract_tiers;
  
  /// Abstracts of the exclusions mask associated to the topology synthesis.  Like other abstracts,
  /// variants at each memory tier are stored.
  /// \{
  std::vector<LocalExclusionMaskReader> lemrs;
  std::vector<SeMaskSynthesisReader> semrs;
  /// \}

  /// Energy tracking abstracts are stored at all memory tiers relevant to the energy tracking
  /// object.
  std::vector<HybridTargetLevel> utracking_abstract_tiers;

  /// Energy tracking is a critical component of interventions that will report back trial moves
  /// and other measurements of the dynamics.
  std::vector<ScoreCardWriter> scws;

  /// The watcher will have memory on both CPU host and GPU device.  This array records which
  /// tier each abstract points to.
  std::vector<HybridTargetLevel> anomaly_abstract_tiers;
  
  /// Anomaly reporting (including for debugging purposes) is accessible via abstracts to any
  /// Watcher class object associated with the simulation.
  std::vector<WatcherWriter> bugws;

  /// Hydrogen bond tracking can be carried out on either the CPU hot or GPU device based on
  /// abstracts kept for either memory tier.
  std::vector<HybridTargetLevel> hbond_abstract_tiers;

  /// Hydrogen bond tracking abstracts
  std::vector<HBondWriter> hbws;
  
  /// Pointers to the original coordinate syntheses, for convenience
  /// \{
  PhaseSpaceSynthesis *poly_ps_ptr;
  Condensate *cdns_ptr;
  /// \}

  /// A pointer to the original neighbor list, or a pair of pointers to the pair of original
  /// neighbor lists.  The pointers are cast to an arbitrary templating scheme, with the true
  /// list of templated types recorded in the ngbr_list_{tcoord, tacc, tcalc, tcoord4}, member
  /// variables above.  The member variable nl_layout (above) stores the configuration in the
  /// number of grids.  If there is only one neighbor list, the first element in the array will
  /// store its address.  Otherwise, the first element in the array will take the electrostatic
  /// neighbor list while the seocnd element takes the Lennard-Jones (van-der Waals) neighbor list.
  std::vector<CellGrid<float, int, float, float4>*> cg_ptr;
  
  /// A pointer to the original topology synthesis, for convenience
  AtomGraphSynthesis *poly_ag_ptr;

  /// A pointer to an implicit solvent workspace, for scratch work
  ImplicitSolventWorkspace *isw_ptr;
  
  /// Pointer to the original set of exclusion masks, for convenience.  When requested, each
  /// pointer is returned const-qualified, reflecting the assumption that bonding patterns will
  /// not change over the course of a simulation, even though individual particle parameters
  /// might change per time step as a result of interventions (i.e. in a polarizable force field).
  /// \{
  LocalExclusionMask *lem_ptr;
  StaticExclusionMaskSynthesis *se_ptr;
  /// \}

  /// A pointer to the energy tracking object used to store the energetic consequences or findings
  /// of the interventions.  This need not be the same energy tracking object used in standard
  /// simulation analysis, but it must be prepared to serve the same number of systems as are
  /// present in the synthesis.
  ScoreCard *sc_ptr;

  /// A pointer to the main monitoring object used to report large forces, failed convergence, and
  /// other anomalies found in a simulation.
  Watcher *anom_ptr;

  // Pointers to various analysis-oriented objects.  Their abstracts can be found in arrays with
  // similar names, below.
  std::vector<HydrogenBondAnalysis*> hbond_trk_ptr;
  
  // The actions that will be taken fall into separate arrays for each enumeration of the
  // IntegrationStage enum class.
  std::vector<DynIntvFuncPtr> force_calc_actions;     ///< Actions to be taken after standard
                                                      ///<   force calculations are complete,
                                                      ///<   possibly including additional force
                                                      ///<   calculations
  std::vector<DynIntvFuncPtr> velocity_adv_actions;   ///< Actions taken after velocity update 1
                                                      ///<   (but before any constraints may be
                                                      ///<   applied)
  std::vector<DynIntvFuncPtr> velocity_cnst_actions;  ///< Actions taken after applying velocity
                                                      ///<   constraints
  std::vector<DynIntvFuncPtr> kinetic_calc_actions;   ///< Actions taken after computing kinetic
                                                      ///<   energies.  The kinetic energy
                                                      ///<   calculation does not change the
                                                      ///<   state of the system, but these
                                                      ///<   actions may benefit from the
                                                      ///<   available energy value.
  std::vector<DynIntvFuncPtr> position_adv_actions;   ///< Actions taken after velocity update 2
                                                      ///<   and (unconstrained) position updates
  std::vector<DynIntvFuncPtr> geometry_cnst_actions;  ///< Actions taken after applying geometry
                                                      ///<   constraints

  // The frequencies at which these actions are taken are stored alongside the function pointers
  // themselves.
  std::vector<int> force_calc_action_intv;     ///< Force calculation action intervals
  std::vector<int> velocity_adv_action_intv;   ///< Velocity update I action intervals
  std::vector<int> velocity_cnst_action_intv;  ///< Velocity constraint action intervals
  std::vector<int> kinetic_calc_action_intv;   ///< Kinetic calculation action intervals
  std::vector<int> position_adv_action_intv;   ///< Position advancement action intervals
  std::vector<int> geometry_cnst_action_intv;  ///< Geometric constraint action intervals

  // Some actions may draw upon one of a number of attached resources of the same object types,
  // e.g. one hydrogen bonding analysis out of several to be carried out in a particular run.
  // These arrays store such indices for each action, or otherwise hold zero if there is no array
  // of resources to consider.
  std::vector<int> force_calc_action_index;     ///< Force calculation action resource indices
  std::vector<int> velocity_adv_action_index;   ///< Velocity update I action resource indices
  std::vector<int> velocity_cnst_action_index;  ///< Velocity constraint action resource indices
  std::vector<int> kinetic_calc_action_index;   ///< Kinetic calculation action resource indices
  std::vector<int> position_adv_action_index;   ///< Position advancement action resource indices
  std::vector<int> geometry_cnst_action_index;  ///< Geometric constraint action resource indices

  // The compute resources (PU host or GPU device) which will execute each action are stored
  // alongside the function pointers, as well.
  std::vector<HybridTargetLevel> force_calc_action_tier;     ///< Indicate whether each force
                                                             ///<   calculation action will be
                                                             ///<   performed by the CPU host or
                                                             ///<   GPU device
  std::vector<HybridTargetLevel> velocity_adv_action_tier;   ///< Velocity advance action tiers
  std::vector<HybridTargetLevel> velocity_cnst_action_tier;  ///< Velocity constraint action tiers
  std::vector<HybridTargetLevel> kinetic_calc_action_tier;   ///< Kinetic energy action tiers
  std::vector<HybridTargetLevel> position_adv_action_tier;   ///< Position advance action tiers
  std::vector<HybridTargetLevel> geometry_cnst_action_tier;  ///< Geometry constraint action tiers

  /// Force calculations and other interventions may require independent workspaces spanning the
  /// synthesis of systems.
  std::vector<PhaseSpaceSynthesis*> workspaces;

  /// Neighbor list forces and composition checks may require independent workspaces of their own.
  /// Because the pointer must be cast to void along all templated parameters, codes for each of
  /// the templated types given in ngbr_list_tcoord, ngbr_list_tacc, ngbr_list_tcalc, and
  /// ngbr_list_tcoord4 are trusted to describe any pointers attached to this list.
  std::vector<CellGrid<float, int, float, float4>*> nl_workspaces;
  
  /// \brief General-purpose function for checking some aspect of the latest attached resource
  ///        against that in other attached resources in the intervention object.
  ///
  /// Overloaded:
  ///   - Check the system count
  ///   - Check the system count and a series of system sizes
  ///
  /// \param system_count  The number of systems covered by the latest attached resource
  /// \param system_sizes  A series of atom counts in systems held within each attached resource
  /// \param atch          Enumerated type of the attachment, to streamline checks against
  ///                      large lists of system sizes
  /// \param caller        Name of the calling function, for backtracing purposes
  /// \{
  void validateContents(int system_count, const char* caller) const;
  void validateContents(const int padded_atom_count, const AttachmentKind atch,
                        const char* caller) const;
  void validateContents(const std::vector<int> &system_sizes, AttachmentKind atch,
                        const char* caller) const;
  /// \}

  /// \brief Validate the requested index of some analysis.
  ///
  /// \param index         Index of the analysis of interest
  /// \param analysis_set  The set of all analyses of the requested type
  /// \param caller        Name of the calling function.  This is for error tracing, but also
  ///                      provides value in describing the type of analysis in question.
  template <typename T>
  void validateAnalysisIndex(int index, const std::vector<T*> &assays, const char* caller);
  
  /// \brief Query the list of step intervals at which to intervene, and modify it or add a new
  ///        interval as needed.  The function seeks a collection of largest common factors which
  ///        cover the applicable intervals of all attached analyses and other interventions.
  ///
  /// \param next_intv  The next active interval to account for
  void setActiveInterval(int next_intv);

  /// \brief Validate a neighbor list (CellGrid) pointer request based on the data types in which
  ///        the pointer is to be returned.
  ///
  /// \param caller  Name of the calling function, for backtracing purposes
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  void validateNeighborListTypes(const char* caller) const;

  /// \brief Validate the theme of a requested neighbor list attachment.  This will check against
  ///        the known layout of the attached neighbor list from the simulation.
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  void validateNLAttachment(const NonbondedTheme theme,
                            const std::vector<CellGrid<float, int, float, float4>*> &attachments,
                            const char* caller) const;
  
  /// \brief Obtain one of the neighbor list (CellGrid) attachments after restoring its proper
  ///        types.  This generalized implementation is called to retrieve pointers to either
  ///        workspaces or the main simulation neighbor list.
  ///
  /// Overloaded:
  ///   - A const-qualified variant returns read-only pointers from an immutable
  ///     DynamicsIntervention object
  ///   - Mutable pointers are returned which allow writing data to the attached neighbor list if
  ///     the DynamicsIntervention object is likewise mutable
  /// \{
  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
  getNLAttachment(const NonbondedTheme theme,
                  std::vector<CellGrid<float, int, float, float4>*> *attachments,
                  const char* caller);

  template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
  const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
  getNLAttachment(const NonbondedTheme theme,
                  const std::vector<CellGrid<float, int, float, float4>*> &attachments,
                  const char* caller) const;
  /// \}
};

/// Make a single, global instance of the class which may be loaded at will and thereafter become
/// accessible to any customizable function.  In this way, all functions for dynamics intervention
/// may share the same basic signature, even though the sorts of information they can access will
/// include any members of this class loaded at runtime.
extern DynamicsIntervention dyna_tk;

// Functions that can be queued under the global instance of DynamicsIntervention (or any other
// object of the class) will draw upon resources attached to it as well as general resources in
// the DynamicsIntervention class itself.  In most cases, attached resources and associated
// evaluator functions will go hand in hand--attaching the resource implies the need to call the
// evaluator, and the evaluator references a prerequisite resource in DynamicsIntervention.  To
// avoid circular dependencies or the need for forward declarations, most of these evaluator
// functions are declared here, below the DynamicsIntervention class itself.  The implementations
// can be found in various files other than dynamics_intervention.cpp, all including this central
// header file.
  
/// \brief Re-evaluate force calculations using methods written for CPU resources.  If performed on
///        the GPU, particle positions will be downloaded to the CPU.  Checks on CPU calculations
///        may use alternative and simpler (but less efficient) methods, but no matter the pathway
///        checks on CPU-based calculations will take place on a copy of the original coordinates
///        with accumulation in new force arrays, to provide independent verification of the
///        original calculations.
///        
/// \param step   The current step number of the simulation
/// \param index  Provide the index of a particular analysis, if more than one is present in the
///               DynamicsIntervention object
/// \param tier   Indicate whether to carry out CPU-based debugging for calculations first
///               performed on the CPU host or GPU device
void execForceDebug(int step, int index = 0, HybridTargetLevel tier = HybridTargetLevel::HOST);

/// \brief Re-evaluate neighbor list composition based on the particle locations found in the
///        coordinate synthesis, re-assigning the populations of each neighbor list cell.
///        Descriptions of input parameters follow from execForceDebug(), above.
void execNLCompositionDebug(int step, int index = 0,
                            HybridTargetLevel tier = HybridTargetLevel::HOST);

/// \brief Carry out the neighbor list examination called by execNLCompositionDebug, above.
///
/// Overloaded:
///   - A templated version for CPU-only simulations
///   - Supply the GPU specifications to differentiate the GPU-enabled mode, which will pull data
///     from the original object's GPU-resident memory onto host data, then carry out the same
///     CPU-based checks
///        
/// \param step  The current step number of the simulation
/// \param gpu   Specifications of the GPU that will carry out the data transfer
/// \{
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void debugNLComposition(int step);

#ifdef STORMM_USE_HPC
void debugNLComposition(int step, const GpuDetails &gpu);
#endif
/// \}

/// \brief Re-evaluate neighbor list-based force calculations using methods written for CPU
///        resources.  If performed on the GPU, the neighbor list and its particle positions will
///        be downloaded to the CPU in a workspace set aside to handle such calculations.
///        
/// \param step   The current step number of the simulation
/// \param index  Provide the index of a particular analysis, if more than one is present in the
///               DynamicsIntervention object
/// \param tier   Indicate whether to carry out CPU-based debugging for calculations first
///               performed on the CPU host or GPU device
void execNLForceDebug(const int step, const int index, const HybridTargetLevel tier);

/// \brief Check for anomalous large forces among all particles.  Descriptions of input parameters
///        follow from execForceDebug(), above.
void evalForceAnomalies(int step, int index = 0, HybridTargetLevel tier = HybridTargetLevel::HOST);
  
/// \brief Evaluate hydrogen bonding using an attached HydrogenBondAnalysis class.
///
/// \param step   The current step number of the simulation
/// \param index  Provide the index of a particular analysis, if more than one is present in the
///               DynamicsIntervention object
/// \param tier   Indicate whether to carry out calculations on the CPU host or GPU device
void evalHydrogenBondAnalysis(int step, int index = 0,
                              HybridTargetLevel tier = HybridTargetLevel::HOST);
  
} // namespace mm
} // namespace stormm

#include "dynamics_intervention.tpp"
#include "Debug/debug_neighbor_list.tpp"

#endif
