// -*-c++-*-
#ifndef STORMM_IMPLICIT_SOLVENT_WORKSPACE_H
#define STORMM_IMPLICIT_SOLVENT_WORKSPACE_H

#include <cmath>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::PrecisionModel;
using trajectory::CoordinateCycle;

/// \brief A simple abstract for the implicit solvent workspace.  There are not readers and writers
///        as the only useful application involves this workspace being writeable.
template <typename T> struct ISWorkspaceKit {

  /// \brief The constructor takes a straight list of arguments.
  ISWorkspaceKit(int fp_bits_in, llint* psi_in, int* psi_overflow_in, llint* sum_deijda_in,
                 int* sum_deijda_overflow_in, llint* alt_psi_in, int* alt_psi_overflow_in,
                 llint* alt_sum_deijda_in, int* alt_sum_deijda_overflow_in);

  /// \brief The copy and move constructors are defaulted, assignment operators implicitly deleted
  ///        due to const members in this object.
  ///
  /// \param original  The object to copy or move
  /// \{
  ISWorkspaceKit(const ISWorkspaceKit &original) = default;
  ISWorkspaceKit(ISWorkspaceKit &&original) = default;
  /// \}

  const int fp_bits;         ///< Fixed-precision bits after the decimal
  const T fp_scale;          ///< The "forward" scaling factor to take real values into the
                             ///<   fixed-precision representation
  const T inv_fp_scale;      ///< The "backward" scaling factor to take fixed-precision values
                             ///<   back into their real number representations and units
  llint* psi;                ///< Accumulators for quantities that will determine effective Born
                             ///<   radii (this will be used exclusively for single-precision
                             ///<   mode calculations, even if a GPU kernel's local thread block
                             ///<   accumulation is handled in the split fixed precision method)
  int* psi_ovrf;             ///< Overflow accumulators for psi (used for double-precision mode
                             ///<   calculations, only)
  llint* sum_deijda;         ///< Accumulators for quantities that will become force
                             ///<   contributions due to derivatives of the effective Born radii
  int* sum_deijda_ovrf;      ///< Overflow accumulators for sum_deijda (used for
                             ///<   double-precision mode calculations, only)

  // The following are alternative arrays made available so that one of the GB kernels can
  // initialize them for use in the following calculation, saving a kernel call in the GB cycle.
  llint* alt_psi;            ///< Alternate array for psi accumulators
  int* alt_psi_ovrf;         ///< Alternate array for psi overflow accumulators
  llint* alt_sum_deijda;     ///< Alternate GB radii derivative accumulators
  int* alt_sum_deijda_ovrf;  ///< Alternate GB radii derivative overflow accumulators
};
  
/// \brief A small collection of arrays to manage temporary accumulators for computing Born radii
///        and their derivatives.
class ImplicitSolventWorkspace {
public:
  
  /// \brief The object contains arrays to store Born radii and derivatives to very high precision.
  ///
  /// Overloaded:
  ///   - Accept an explicit bit count for storing Born radii and derivative information in
  ///     fixed-precision format.
  ///   - Accept a precision model and apply an automatic scaling factor for the Born radii and
  ///     derivatives.
  ///
  /// \param atom_starts  Starting positions of the atoms for each system in the coordinate arrays
  ///                     of the related AtomGraphSynthesis and PhaseSpaceSynthesis
  /// \param atom_counts  Atom counts in each system of the related syntheses
  /// \param bit_count    The number of bits to keep after the decimal (units of Angstroms for the
  ///                     radii, or kcal/mol-Angstroms per change in radius)
  /// \param prec         The precision model to prepare for (automates selection of bit_count)
  /// \{
  ImplicitSolventWorkspace(const Hybrid<int> &atom_starts, const Hybrid<int> &atom_counts,
                           int bit_count);

  ImplicitSolventWorkspace(const Hybrid<int> &atom_starts, const Hybrid<int> &atom_counts,
                           PrecisionModel prec);
  /// \}

  /// \brief Copy and move constructors, as well as the respective assignment operators, take the
  ///        usual form for an object with POINTER-kind Hybrids that require repair.
  ///
  /// \param original  The object to copy or move
  /// \param other     Right-hand side object for assignment operators
  /// \{
  ImplicitSolventWorkspace(const ImplicitSolventWorkspace &original);
  ImplicitSolventWorkspace(ImplicitSolventWorkspace &&original);
  ImplicitSolventWorkspace& operator=(const ImplicitSolventWorkspace &other);
  ImplicitSolventWorkspace& operator=(ImplicitSolventWorkspace &&other);
  /// \}
  
  /// \brief Get the precision in which Born radii and sums of the per-atom energy derivatives
  ///        are stored.
  int getFixedPrecisionBits() const;

  /// \brief Get the object's current coordinate cycle position.
  CoordinateCycle getCyclePosition() const;
  
  /// \brief Get the double-precision abstract, containing pointers to data on the host or device.
  ///
  /// Overloaded:
  ///   - Provide the stage of the coordinate cycle towards which pointers will orient
  ///   - Orient pointers towards the object's current coordinate cycle stage
  ///
  /// \param tier         Level at which to retrieve pointers (the CPU host, or GPU device)
  /// \param orientation  Perform the zeroing on psi / sumdeijda or alt_psi / alt_sumdeijda arrays
  /// \{
  ISWorkspaceKit<double> dpData(CoordinateCycle orientation,
                                HybridTargetLevel tier = HybridTargetLevel::HOST);

  ISWorkspaceKit<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Get the single-precision abstract, containing pointers to data on the host or device.
  ///        Overloading and descriptions of input parameters follow from dpData(), above.
  /// \{
  ISWorkspaceKit<float> spData(CoordinateCycle orientation,
                               HybridTargetLevel tier = HybridTargetLevel::HOST);

  ISWorkspaceKit<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

#ifdef STORMM_USE_HPC
  /// \brief Upload the object's contents from the host to the GPU device.  This may be useful in
  ///        debugging scenarios.
  void upload();

  /// \brief Download the object's contents from the GPU device to the host.  This is also useful
  ///        in debugging scenarios.
  void download();
#endif

  /// \brief Increment the object's cycle counter.
  void updateCyclePosition();
  
  /// \brief Set the Generalized Born radii and radii derivative accumulators to zero.
  ///
  /// \param tier         Do the zeroing on the host or on the GPU device
  /// \param orientation  Perform the zeroing on psi / sumdeijda or alt_psi / alt_sumdeijda arrays
  /// \param gpu          Details of the GPU in use
  void initialize(HybridTargetLevel tier = HybridTargetLevel::HOST,
                  CoordinateCycle orientation = CoordinateCycle::WHITE,
                  const GpuDetails &gpu = null_gpu);

#ifdef STORMM_USE_HPC
  /// \brief Launch a small kernel to manage initialization of this object on the HPC device.
  ///
  /// \param orientation  Perform the zeroing on psi / sumdeijda or alt_psi / alt_sumdeijda arrays
  void launchInitialization(const GpuDetails &gpu,
                            CoordinateCycle orientation = CoordinateCycle::WHITE);
#endif

private:
  int fp_bits;                          ///< Fixed-precision bits after the decimal
  int padded_atom_count;                ///< The number of atoms in the system, including padding
                                        ///<   between separate systems (the starting and ending
                                        ///<   points of each system are not stored in this object,
                                        ///<   as it can only function in the context of various
                                        ///<   synthesis objects which will have this information)
  CoordinateCycle cycle_position;       ///< The current position in the coordinate cycle.  This
                                        ///<   can be used to guide the production of abstracts to
                                        ///<   the relevant data arrays.
  Hybrid<llint> psi;                    ///< The accumulates psi values for ecah atom, which are
                                        ///<   then transformed, in place, to become the effective
                                        ///<   GB radii
  Hybrid<int> psi_overflow;             ///< Overflow accumulators for psi
  Hybrid<llint> sum_deijda;             ///< Sums of the energy derivative for each atom i making
                                        ///<   pair contributions with all other atoms j with
                                        ///<   respect to an order parameter
  Hybrid<int> sum_deijda_overflow;      ///< Overflow accumulators for sum_deijda
  Hybrid<llint> alt_psi;                ///< Alternate array for accumulating psi values--one can
                                        ///<   be in use while the other is being reset
  Hybrid<int> alt_psi_overflow;         ///< Overflow accumulators for the alternate psi array
  Hybrid<llint> alt_sum_deijda;         ///< Alternative array for storing GB radii derivatives
  Hybrid<int> alt_sum_deijda_overflow;  ///< Alternative array for GB radial derivative overflows
  Hybrid<llint> llint_data;             ///< ARRAY-kind Hybrid targeted by the preceding
                                        ///<   POINTER-kind objects to store all long long integer
                                        ///<   data
  Hybrid<int> int_data;                 ///< ARRAY-kind Hybrid targeted by the preceding
                                        ///<   POINTER-kind objects to store all overflow data

  /// \brief Repair POINTER-kind Hybrid objects in a newly copied ImplicitSolventWorkspace by
  ///        making them point to the object's own ARRAY-kind Hybrids.
  void rebasePointers();
};

} // namespace synthesis
} // namespace stormm

#include "implicit_solvent_workspace.tpp"

#endif
