// -*-c++-*-
#ifndef STORMM_MM_CONTROLS_H
#define STORMM_MM_CONTROLS_H

#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Potential/energy_enumerators.h"
#include "Potential/cellgrid.h"
#include "Math/reduction_enumerators.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_minimize.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Topology/atomgraph_enumerators.h"

namespace stormm {
namespace mm {

using card::GpuDetails;
using card::CoreKlManager;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::PrecisionModel;
using energy::CellGrid;
using energy::ClashResponse;
using energy::EvaluateEnergy;
using energy::EvaluateForce;
using energy::NeighborListKind;
using energy::QMapMethod;
using energy::TinyBoxPresence;
using stmath::ReductionStage;
using namelist::default_dynamics_time_step;
using namelist::default_electrostatic_cutoff;
using namelist::default_minimize_dx0;
using namelist::default_minimize_maxcyc;
using namelist::default_minimize_ncyc;
using namelist::default_rattle_tolerance;
using namelist::default_nt_warp_multiplicity;
using namelist::default_van_der_waals_cutoff;
using namelist::DynamicsControls;
using namelist::MinimizeControls;
using synthesis::AtomGraphSynthesis;
using synthesis::VwuGoal;
using topology::ImplicitSolventModel;

/// \brief The C-style, always writeable abstract for the MolecularMechanicsControls object.  To
///        not be able to modify this object's contents would be nonsensical, as it is intended to
///        to keep counters of the simulation time step as well as force evaluation work units.
template <typename T> struct MMControlKit {

  /// \brief The constructor takes a straight list of values and pointers.  The step number is
  ///        left modifiable so that the object can be re-used over successive time steps.
  MMControlKit(int step_in, int sd_cycles_in, int max_cycles_in, T initial_step_in,
               int nt_warp_mult_in, const T elec_cut_in, const T vdw_cut_in, int* vwu_progress_in,
               int* vupt_progress_in, int* vcns_progress_in, int* pupt_progress_in,
               int* gcns_progress_in, int* nbwu_progress_in, int* pmewu_progress_in,
               int* gbrwu_progress_in, int* gbdwu_progress_in, int* gtwu_progress_in,
               int* scwu_progress_in, int* rdwu_progress_in);

  /// \brief The usual copy and move constructors for an abstract apply here.  
  /// \{
  MMControlKit(const MMControlKit &original) = default;
  MMControlKit(MMControlKit &&original) = default;
  /// \}

  int step;                ///< The current simulation step
  const int sd_cycles;     ///< The number of steepest-descent energy minimization cycles
  const int max_cycles;    ///< The total number of energy minimization cycles or dynamics steps
  const T initial_step;    ///< Initial step size to be taken in energy minimization
  const int nt_warp_mult;  ///< The number of warps to use, per neighbor list decomposition cell,
                           ///<   in neutral territory tile processing
  const T elec_cut;        ///< The cutoff for electrostatic interactions (this will be applied
                           ///<   only if dual neighbor list grids are in effect)
  const T vdw_cut;         ///< The cutoff for van-der Waals interactions (this will be the sole
                           ///<   cutoff if dual neighbor list grids are not in effect)
  const T elec_cut_sq;     ///< Squared cutoff for electrostatic interactions
  const T vdw_cut_sq;      ///< Squared cutoff for van-der Waals interactions
  int* vwu_progress;       ///< Progress counters for valence work units
  int* vupt_progress;      ///< Progress counters for standalone velocity update work units
  int* vcns_progress;      ///< Progress counters for standalone velocity constraint work units
  int* pupt_progress;      ///< Progress counters for standalone coordinate update work units
  int* gcns_progress;      ///< Progress counters for standalone positional constraint work units
  int* nbwu_progress;      ///< Progress counters for non-bonded work units
  int* pmewu_progress;     ///< Progress counters for PME long-ranged work units
  int* gbrwu_progress;     ///< Progress counters for Generalized Born radii computations
  int* gbdwu_progress;     ///< Progress counters for Generalized Born derivative computations
  int* gtwu_progress;      ///< Progress counters for gathering work units
  int* scwu_progress;      ///< Progress counters for scattering work units
  int* rdwu_progress;      ///< Progress counters for reduction work units
};

/// \brief A collection of contol data for molecular mechanics simulations, conveying the current
///        step number, progress counters through various work units, the time step, RATTLE
///        tolerance, and other critical parameters to guide calculations.  This common struct
///        accommodates both molecular dynamics and molecular mechanics energy minimizations, and
///        can be created from control objects derived from &dynamics or &minimize namelists.
class MolecularMechanicsControls {
public:

  /// \brief The constructor can create an empty object with default parameters for the time step
  ///        and rattle tolerance, accept user-specified values for critical constants, or be
  ///        built from a namelist-derived object on the CPU.  The reason that those
  ///        namelist-derived control objects don't just get pushed by value to the GPU is that the
  ///        work unit counters need special arrays, which this object manages.
  ///
  /// \param initial_step_in  The desired initial step (for energy minimization)
  /// \param user_input       Namelist-derived molecular dynamics or energy minimization controls
  /// \{
  MolecularMechanicsControls(double initial_step_in = default_minimize_dx0,
                             int sd_cycles_in = default_minimize_ncyc,
                             int max_cycles_in = default_minimize_maxcyc,
                             int nt_warp_multiplicity_in = default_nt_warp_multiplicity,
                             double electrostaitc_cutoff_in = default_electrostatic_cutoff,
                             double van_der_waals_cutoff_in = default_van_der_waals_cutoff);

  MolecularMechanicsControls(const DynamicsControls &user_input);
                             
  MolecularMechanicsControls(const MinimizeControls &user_input);
  /// \}

  /// \brief The copy constructor handles assignment of internal POINTER-kind Hybrid objects.
  ///
  /// \param original  The object to copy
  MolecularMechanicsControls(const MolecularMechanicsControls &original);

  /// \brief The copy assignment operator likewise handles assignment of internal POINTER-kind
  ///        Hybrid objects.
  ///
  /// \param other  Another way to say original, in a different semantic context
  MolecularMechanicsControls& operator=(const MolecularMechanicsControls &other);

  /// \brief The move constructor prepared the original object for destruction.
  ///
  /// \param original  The object from which to take contents
  MolecularMechanicsControls(MolecularMechanicsControls &&original);

  /// \brief Get the current step number.
  int getStepNumber() const;
  
  /// \brief Get the number of steepest descent minimization cycles.
  int getSteepestDescentCycles() const;
  
  /// \brief Get the total number of minimization steps, or total MD cycles.
  int getTotalCycles() const;

  /// \brief Get the neutral-territory warp multiplicity.
  int getNTWarpMultiplicity() const;

  /// \brief Get the initial step for energy minimization.
  double getInitialMinimizationStep() const;

  /// \brief Get the cutoff for non-bonded electrostatic interactions in periodic simulations.
  double getElectrostaticCutoff() const;
  
  /// \brief Get the cutoff for non-bonded van-der Waals interactions in periodic simulations.
  double getVanDerWaalsCutoff() const;
  
  /// \brief Get the value of one of the valence work unit progress counters on the host or the
  ///        HPC device.
  ///
  /// \param counter_index  Index of the work unit counter to retrieve
  /// \param tier           Get the result from the CPU host or the HPC device
  int getValenceWorkUnitProgress(int counter_index,
                                 HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the value of one of the non-bonded work unit progress counters on the host or the
  ///        HPC device.
  ///
  /// \param counter_index  Index of the work unit counter to retrieve
  /// \param tier           Get the result from the CPU host or the HPC device
  int getNonbondedWorkUnitProgress(int counter_index,
                                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the value of one of the Particle-Mesh Ewald work unit progress counters on the
  ///        host or the HPC device.
  ///
  /// \param counter_index  Index of the work unit counter to retrieve
  /// \param tier           Get the result from the CPU host or the HPC device
  int getPmeWorkUnitProgress(int counter_index,
                             HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the value of one of the reduction work unit progress counters on the host or the
  ///        HPC device.
  ///
  /// \param counter_index  Index of the work unit counter to retrieve
  /// \param tier           Get the result from the CPU host or the HPC device
  int getReductionWorkUnitProgress(int counter_index, ReductionStage process,
                                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief The move assignment operator looks much like the copy assignment operator.
  ///
  /// \param other  The object from which to take contents
  MolecularMechanicsControls& operator=(MolecularMechanicsControls &&other);
  
  /// \brief Obtain a double-precision abstract for this object.
  MMControlKit<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Obtain a single-precision abstract for this object.
  MMControlKit<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Prime the work unit counters based on a particular GPU configuration.
  ///
  /// Overloaded:
  ///   - Provide information about the PME approach
  ///   - Make provisions for clashes to be forgiven (or not) with softcore potential functions
  ///   - Provide only basic information for a vacuum-phase MM calculation
  /// 
  /// \param launcher          Object containing launch parameters for all kernels
  /// \param eval_frc          Indicate whether forces are to be evaluated--some kernels allocate
  ///                          different numbers of blocks in the launch grid if forces are not
  ///                          required.
  /// \param eval_nrg          Indicate whether energy is to be evaluated--some kernels allocate
  ///                          different numbers of blocks in the launch grid if the energy is not
  ///                          required.
  /// \param softcore          Indicate whether the control counters should take clash handling
  ///                          kernel variants into account
  /// \param valence_prec      Precision model for valence calculations
  /// \param nonbond_prec      Precision model for non-bonded calculations
  /// \param general_prec      Precision model for all calculations
  /// \param qspread_approach  The approach taken by GPU kernels to spread charge density onto the
  ///                          particle-mesh interaction grids
  /// \param acc_prec          Precision model for accumulations (needed in case the density
  ///                          mapping approach is accumulation in the __shared__ memory space)
  /// \param image_coord_type  Data type used to represent coordinates in the neighbor list images
  /// \param qspread_order     The order of B-spline interpolation used to spread charge density to
  ///                          the particle-mesh interaction grids
  /// \param poly_ag           Compilation of topologies describing the workload (used here for
  ///                          general descriptors such as the non-bonded work unit type)
  /// \{
  void primeWorkUnitCounters(const CoreKlManager &launcher, EvaluateForce eval_frc,
                             EvaluateEnergy eval_nrg, ClashResponse softcore, VwuGoal purpose,
                             PrecisionModel valence_prec, PrecisionModel nonbond_prec,
                             QMapMethod qspread_approach, PrecisionModel acc_prec,
                             size_t image_coord_type, int qspread_order,
                             NeighborListKind nbgr_config, TinyBoxPresence has_tiny_box,
                             const AtomGraphSynthesis &poly_ag);

  void primeWorkUnitCounters(const CoreKlManager &launcher, EvaluateForce eval_frc,
                             EvaluateEnergy eval_nrg, const ClashResponse softcore,
                             VwuGoal purpose, PrecisionModel valence_prec,
                             PrecisionModel nonbond_prec, const AtomGraphSynthesis &poly_ag); 

  void primeWorkUnitCounters(const CoreKlManager &launcher, EvaluateForce eval_frc,
                             EvaluateEnergy eval_nrg, VwuGoal purpose, PrecisionModel valence_prec,
                             PrecisionModel nonbond_prec, const AtomGraphSynthesis &poly_ag); 

  void primeWorkUnitCounters(const CoreKlManager &launcher, EvaluateForce eval_frc,
                             EvaluateEnergy eval_nrg, VwuGoal purpose, PrecisionModel general_prec,
                             const AtomGraphSynthesis &poly_ag); 
  /// \}

  /// \brief Set the neutral-territory warp multiplicity based on one or two gell grids, for a
  ///        given GPU.
  ///
  /// Overloaded:
  ///   - Accept one neighbor list and specifications the GPU that will perform the calculations
  ///   - Accept two neighbor lists and the specifications of the GPU
  ///
  /// \param cg_a     The first cell grid neighbor list
  /// \param cg_b     The second (optional) cell grid neighbor list
  /// \param gpu      Details of the GPU that will evaluate pairwise interactions
  /// \param mult_in  The multiplicity to explicitly set
  /// \{
  template <typename T, typename Tacc, typename Tcalc, typename T4>
  void setNTWarpMultiplicity(const CellGrid<T, Tacc, Tcalc, T4> *cg_a,
                             const CellGrid<T, Tacc, Tcalc, T4> *cg_b, const GpuDetails &gpu);

  template <typename T, typename Tacc, typename Tcalc, typename T4>
  void setNTWarpMultiplicity(const CellGrid<T, Tacc, Tcalc, T4> *cg_a, const GpuDetails &gpu);

  void setNTWarpMultiplicity(int mult_in);
  /// \}
  
  /// \brief Increment the step counter, moving the controls to a different progress counter.
  void incrementStep();
  
#ifdef STORMM_USE_HPC
  /// \brief Upload the object's contents to the device (needed so that CPU-primed work unit
  ///        counters can go into effect)
  void upload();

  /// \brief Download the object's contents from the device (useful for debugging)
  void download();
#endif
  
private:
  int step_number;              ///< The step counter for the simulation
  int sd_cycles;                ///< The number of steepest-descent energy minimization cycles
  int max_cycles;               ///< Total number of energy minimization cycles or dynamics steps
  double initial_step;          ///< The initial step length (for all systems, if these controls
                                ///<   govern a synthesis) to take in energy minimization
                                ///<   calculations
  int nt_warp_multiplicity;     ///< The number of warps that should cooperate to evaluate all pair
                                ///<   interactions in a cell grid neighbor list
  double electrostatic_cutoff;  ///< The cutoff to be applied to electrostatic interactions in
                                ///<   periodic simulations
  double van_der_waals_cutoff;  ///< The cutoff to be applied to van-der Waals interactions in
                                ///<   periodic simulations

  /// Progress counters through valence work units, an array of two times the number of lanes per
  /// warp so that the valence work units kernel can perform resets of elements in this array with
  /// maximum efficiency.  The kernel will work based on the progress counter with index equal to
  /// the step number modulo twice the warp size.  When the step counter is equal to the a
  /// multiple of the warp size, one of the warps in the first block of this kernel will reset
  /// counters in array indices [ 0, warp size ).  When the step counter is equal to a multiple of
  /// twice the warp size, one of the warps in this kernel will reset the counters in array indices
  /// [ warp size, 2 * warp size ).
  Hybrid<int> vwu_progress;

  /// Progress counters for stand-alone velocity update work units.  These will track the valence
  /// work units in terms of their overall number, and even reference valence work unit data, but
  /// the kernels that launch them may be of different sizes and therefore the progress counters
  /// be set and reset differently by the standalone kernels.  In typical simulations, the valence
  /// work units will drive all processes in trajectory integration, bypassing these counters, but
  /// for specific applications requiring intervention amidst a modular and decoupled process,
  /// these progress counters will be made available.
  Hybrid<int> velocity_update_progress;

  /// Stand-alone velocity constraint work unit progress counters.  See the array member variable
  /// velocity_update_progress for motivation.
  Hybrid<int> velocity_constraint_progress;

  /// Stand-alone particle position update work unit progress counters.  See the array member
  /// variable velocity_update_progress for motivation.
  Hybrid<int> position_update_progress;

  /// Stand-alone geometry constraint work unit progress counters.  See the array member variable
  /// velocity_update_progress for motivation.
  Hybrid<int> geometry_constraint_progress;
  
  /// Progress counters for non-bonded work units, arrayed in a manner like the valence work units
  Hybrid<int> nbwu_progress;

  /// Progress counters for long-ranged, mesh-based PME work units
  Hybrid<int> pmewu_progress;

  /// Progress counters for Generalized Born radii computation work units
  Hybrid<int> gbrwu_progress;

  /// Progress counters for Generalized Born derivative computation work units
  Hybrid<int> gbdwu_progress;
  
  /// Progress counters for gathering information across each system.
  Hybrid<int> gather_wu_progress;

  /// Progress counters for scattering gathered information back across each system.
  Hybrid<int> scatter_wu_progress;
  
  /// Progress through all-reduce work units.  While it may seem tedious, these are necessary in
  /// order to provide a general solution to the problem of accumulating results across an
  /// arbitrary number of systems with a fixed number of thread blocks.
  Hybrid<int> all_reduce_wu_progress;

  /// ARRAY-kind Hybrid object targeted by the above POINTER-kind Hybrid objects
  Hybrid<int> int_data;

  /// \brief Re-assign pointers to reference the object's private integer array data.  This is
  ///        invoked after a copy construction or copy assignment.
  void rebasePointers();
};

} // namespace mm
} // namespace stormm

#include "mm_controls.tpp"

#endif
