// -*-c++-*-
#ifndef STORMM_LINE_MINIMIZATION_H
#define STORMM_LINE_MINIMIZATION_H

#include "copyright.h"
#include "Accelerator/hybrid.h"

namespace stormm {
namespace mm {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief Abstract for the line minimization object incorporating write access.
struct LinMinWriter {

  /// \brief The constructor works like any other abstract, taking arguments in the order that
  ///        fills the member variables.
  LinMinWriter(const int nsys_in, double* l_move_in, double* s_move_in, double* mfac_a_in,
               double* mfac_b_in, double* mfac_c_in, double* nrg_a_in, double* nrg_b_in,
               double* nrg_c_in, double* nrg_d_in);

  /// \brief Take the standard copy and move constructors for a struct with const elements.
  /// \{
  LinMinWriter(const LinMinWriter &original) = default;
  LinMinWriter(LinMinWriter &&original) = default;
  /// \}
  
  // Member variables
  const int nsys;  ///< The number of systems subject to their own independent line minimizations
  double* l_move;  ///< Lengths of the moves for each system after the initial force and energy
                   ///<   computation.  The product of this with mfac_a determines the total
                   ///<   displacement of all atoms in the system.
  double* s_move;  ///< During the first of three moves along the chosen line, the contents of
                   ///<   l_move are copied into this vector for each system.  This vector is then
                   ///<   used as the reference for the baseline move length for the second and
                   ///<   third moves along the chosen line.  During the third move, the contents
                   ///<   of s_move will be multiplied by the final move scaling factor and used
                   ///<   to update l_move in preparation for the next cycle of line minimization.
                   ///<   This protects against race conditions without forcing work on the CPU,
                   ///<   such as explicitly swapping pointers to alternating l_move arrays
                   ///<   between steps.
  double* mfac_a;  ///< Move multiplying factors of the 1st move along each system's line
  double* mfac_b;  ///< Move multiplying factors of the 2nd move along each system's line
  double* mfac_c;  ///< Move multiplying factors of the 3rd move along each system's line
  double* nrg_a;   ///< Energies obtained for each system at the outset of the cycle
  double* nrg_b;   ///< Energies obtained for each system after the 1st move in the cycle
  double* nrg_c;   ///< Energies obtained for each system after the 2nd move in the cycle
  double* nrg_d;   ///< Energies obtained for each system after the 3rd move in the cycle
};

/// \brief Abstract for the line minimization object incorporating write access.
struct LinMinReader {

  /// \brief The constructor works like any other abstract, taking arguments in the order that
  ///        fills the member variables.
  LinMinReader(const int nsys_in, const double* l_move_in, const double* s_move_in,
               const double* mfac_a_in, const double* mfac_b_in, const double* mfac_c_in,
               const double* nrg_a_in, const double* nrg_b_in, const double* nrg_c_in,
               const double* nrg_d_in);

  /// \brief Take the standard copy and move constructors for a struct with const elements.
  /// \{
  LinMinReader(const LinMinReader &original) = default;
  LinMinReader(LinMinReader &&original) = default;
  /// \}
  
  // Member variables
  const int nsys;  ///< The number of systems subject to their own independent line minimizations
  const double* l_move;  ///< Lengths of the moves for each system in the 1st step
  const double* s_move;  ///< Lengths of the moves for each system in the 2nd and 3rd steps
  const double* mfac_a;  ///< Move multiplying factors of the 1st move along each system's line
  const double* mfac_b;  ///< Move multiplying factors of the 2nd move along each system's line
  const double* mfac_c;  ///< Move multiplying factors of the 3rd move along each system's line
  const double* nrg_a;   ///< Energies obtained for each system at the outset of the cycle
  const double* nrg_b;   ///< Energies obtained for each system after the 1st move in the cycle
  const double* nrg_c;   ///< Energies obtained for each system after the 2nd move in the cycle
  const double* nrg_d;   ///< Energies obtained for each system after the 3rd move in the cycle
};

/// \brief This object serves a synthesis of systems.  Hold the move multipliers for all of them
///        and energies obtained by moving various distances along the conjugate gradient vector.
///        This provides a basis for solving a cubic polynomial to identify the best move along
///        the computed gradient.
class LineMinimization {
public:

  /// \brief The constructor takes only the number of systems as input and allocates for three
  ///        moves along the computed gradient.
  ///
  /// \param system_count_in  The number of systems to prepare for line minimizations upon
  /// \param dx0              Initial step size for total movement of all atoms along the line
  LineMinimization(int system_count_in = 0, double dx0 = 0.0);

  /// \brief The copy constructor and copy assignment operator must perform pointer repair.  The
  ///        default move and move assignment operators are the best choice.
  ///
  /// \param original  The original object, to be replicated or moved
  /// \param other     Another object, to be replicated or moved
  /// \{
  LineMinimization(const LineMinimization &original);
  LineMinimization& operator=(const LineMinimization &other);
  LineMinimization(LineMinimization &&original) = default;
  LineMinimization& operator=(LineMinimization &&other) = default;
  /// \}
  
  /// \brief Get the number of systems
  int getSystemCount() const;

  /// \brief Get the current move lengths for one or more systems.
  ///
  /// Overloaded:
  ///   - Get the move lengths for all systems for a particular move.
  ///   - Get the move lengths of a particular system and a specified move.
  ///
  /// \param system_index  Index of the system of interest (if unspecified, the move lengths for
  ///                      all systems will be returned)
  /// \param tier          Obtain results from the host or the device
  std::vector<double> getMoveLength(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double getMoveLength(int system_index,
                       HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the movement multipliers for one or more systems defining a particular move.
  ///        Collectively, the system's particles move such that the magnitude of the displacment
  ///        vector is the instantaneous move length times the current move scaling factor.  A
  ///        move length of 0.0001A might be multiplied by factors of 1.0, 1.05, 1.025, and 1.020
  ///        as the line minimization proceeds along one computed gradient.
  ///
  /// Overloaded:
  ///   - Get the factors of all systems for a particular move.
  ///   - Get the factors of a particular system and a specified move.
  ///
  /// \param move_index    Index of the move (0 to 3)
  /// \param system_index  Index of the system of interest (if unspecified, the move lengths for
  ///                      all systems will be returned)
  /// \param tier          Obtain results from the host or the device
  /// \{
  std::vector<double> getMoveFactor(int move_index,
                                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double getMoveFactor(int move_index, int system_index,
                       HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}
  

  /// \brief Get the energies for one or more systems after a particular move.
  ///
  /// Overloaded:
  ///   - Get all energies after a particular move.
  ///   - Get the energy for a particular system after a particular move.
  ///
  /// \param move_index    Index of the move (0 to 3)
  /// \param system_index  Index of the system of interest (if unspecified, the energies of all
  ///                      systems will be returned)
  /// \param tier          Obtain results from the host or the device
  /// \{
  std::vector<double> getEnergy(int move_index,
                                HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double getEnergy(int move_index, int system_index,
                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get a read-only or writeable abstract, as appropriate.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const LineMinimization object.
  ///   - Get an abstract with write access for a non-const LineMinimization object.
  ///
  /// \param tier  Leve at which to obtain pointers: on the host, or on the HPC device
  /// \{
  LinMinReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  LinMinWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Prime the initial movement lengths for a series of line minimization cycles using a
  ///        default or user-specified value.  In the HPC compilation, this function will
  ///        automatically upload the present move lengths to the device.
  ///
  /// \param dx0  The initial movelength length to take
  void primeMoveLengths(double dx0);

#ifdef STORMM_USE_HPC
  /// \brief Upload the contents of the object from host to device memory.  This is not necessary
  ///        except for testing purposes.  Movement lengths set of the host by primeMoveLengths()
  ///        will be automatically uploaded.
  void upload();

  /// \brief Download the contents of the object from device memory to the host.  This provides
  ///        debugging capability.
  void download();
#endif

private:
  int system_count;             ///< The number of systems and the length of data controlled by
                                ///<   each of the following POINTER-kind Hybrid objects
  Hybrid<double> move_length;   ///< Lengths of the moves to make along the gradient line for each
                                ///<   system in the first move.  Three such steps are made per
                                ///<   line minimization, with this base length, but to allow
                                ///<   updates without race conditions its values for each system
                                ///<   are stored in the following move_save array during the 1st
                                ///<   step.  The 2nd and 3rd steps pull the values from move_save
                                ///<   and the 3rd step stores the updated value of the baseline
                                ///<   move in back in move_length.  These baseline move lengths
                                ///<   evolve incrementally based on the success or failure of
                                ///<   each step to reduce the energy.
  Hybrid<double> save_length;   ///< Baseline lengths of moves along the gradient line, a holding
                                ///<   array to avoid race conditions
  Hybrid<double> move_factor_a; ///< Move multiplication factor of the 1st move along the gradient
  Hybrid<double> move_factor_b; ///< Move multiplication factor of the 2nd move along the gradient
  Hybrid<double> move_factor_c; ///< Move multiplication factor of the 3rd move along the gradient
  Hybrid<double> energy_a;      ///< Energy obtained before the first move along the gradient (at
                                ///<   the outset of the line minimization cycle)
  Hybrid<double> energy_b;      ///< Energy obtained after the first move along the gradient
  Hybrid<double> energy_c;      ///< Energy obtained after the second move along the gradient
  Hybrid<double> energy_d;      ///< Energy obtained after the third move along the gradient
  Hybrid<double> storage;       ///< ARRAY-kind Hybrid object targeted by all previous objects
};

} // namespace mm
} // namespace stormm

#endif
