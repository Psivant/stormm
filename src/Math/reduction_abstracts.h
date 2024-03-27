// -*-c++-*-
#ifndef STORMM_REDUCTION_ABSTRACTS_H
#define STORMM_REDUCTION_ABSTRACTS_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "reduction_bridge.h"
#include "reduction_enumerators.h"

namespace stormm {
namespace stmath {

using card::HybridTargetLevel;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;

/// \brief Collect the simple components needed to guide reductions across all systems in a
///        topology synthesis (or the corresponding compilation of coordinates): the number of
///        work units, the strategy, and the list of work unit abstracts.
struct ReductionKit {

  /// \brief The constructor takes a straight list of values for all member variables, or a
  ///        combination of familiar objects.
  ///
  /// \param poly_ag  Compilation of topologies with reduction work unit counts and abstracts
  /// \param tier     Level at which to obtain the data arrays (needed for passing to subsequent
  ///                 getter functions in the topology and coordinate syntheses)
  /// \{
  ReductionKit(int nrdwu_in, RdwuPerSystem rps_in, const int* rdwu_abstracts_in,
               const int* atom_counts_in);

  ReductionKit(const AtomGraphSynthesis &poly_ag,
               HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Take the typical copy and move constructors for an abstract with constants.
  /// \{
  ReductionKit(const ReductionKit &original) = default;
  ReductionKit(ReductionKit &&original) = default;
  /// \}

  const int nrdwu;            ///< The number of reduction work units in the synthesis
  const RdwuPerSystem rps;    ///< Are there one or multiple work units serving each system?
  const int* rdwu_abstracts;  ///< Reduction work unit abstracts (strides of rdwu_abstract_length)
  const int* atom_counts;     ///< The number of atoms per system (for normalization purposes)
};

/// \brief Collect pointers to data subject to reduction operations.  Reductions can happen on up
///        to three data sets at once, with one consistent fixed-precision scaling factor (if
///        fixed-precision applies).
template <typename T> struct GenericRdSubstrate {

  /// \brief The constructor can take up to three buffers, each with possible overflows to
  ///        accommodate extended fixed-precision formats.
  ///
  /// Overloaded:
  ///   - Prepare to reduce one array of standard data
  ///   - Prepare to reduce three arrays of standard data
  ///   - Accommodate extended fixed-precision formats
  /// \{
  GenericRdSubstrate(const T* x_read_in, double* x_buffer_in, T* x_write_in,
                     const int scale_bits_in = 0);

  GenericRdSubstrate(const T* x_read_in, const T* y_read_in, const T* z_read_in,
                     double* x_buffer_in, double* y_buffer_in, double* z_buffer_in, T* x_write_in,
                     T* y_write_in, T* z_write_in, const int scale_bits_in = 0);

  GenericRdSubstrate(const T* x_read_in, const int* x_read_ovrf_in, double* x_buffer_in,
                     T* x_write_in, int* x_write_ovrf_in, const int scale_bits_in = 0);

  GenericRdSubstrate(const T* x_read_in, const int* x_read_ovrf_in, const T* y_read_in,
                     const int* y_read_ovrf_in, const T* z_read_in, const int* z_read_ovrf_in,
                     double* x_buffer_in, double* y_buffer_in, double* z_buffer_in, T* x_write_in,
                     int* x_write_ovrf_in, T* y_write_in, int* y_write_ovrf_in, T* z_write_in,
                     int* z_write_ovrf_in, const int scale_bits_in = 0);
  /// \}

  /// \brief Take the typical copy and move constructors for an abstract with const elements.
  /// \{
  GenericRdSubstrate(const GenericRdSubstrate &original) = default;
  GenericRdSubstrate(GenericRdSubstrate &&original) = default;
  /// \}

  const int dim;                ///< The number of dimensions to the data involved in the
                                ///<   reduction, i.e. X/Y/Z coordinate reduction has dimension 3
  const int scale_bits;         ///< The number of bits after the decimal in fixed-precision data
  const double fp_scaling;      ///< Scaling factor for fixed-precision data
  const double inv_fp_scaling;  ///< Inverse of the scaling factor for fixed-precision data
  const T* x_read;              ///< Read-only arrays of reducible data for the 1st dimension
  const T* y_read;              ///< Read-only arrays of reducible data for the 2nd dimension
  const T* z_read;              ///< Read-only arrays of reducible data for the 3rd dimension
  const int* x_read_ovrf;       ///< Overflow arrays for extended fixed-precision format in the
                                ///<   1st dimension of reducible data
  const int* y_read_ovrf;       ///< Overflow arrays for extended fixed-precision format in the
                                ///<   2nd dimension of reducible data
  const int* z_read_ovrf;       ///< Overflow arrays for extended fixed-precision format in the
                                ///<   3rd dimension of reducible data
  double* x_buffer;             ///< Buffer for work unit sums along the 1st dimension
  double* y_buffer;             ///< Buffer for work unit sums along the 2nd dimension
  double* z_buffer;             ///< Buffer for work unit sums along the 3rd dimension
  T* x_write;                   ///< Array accepting results of scattering in the 1st dimension
  T* y_write;                   ///< Array accepting results of scattering in the 2nd dimension
  T* z_write;                   ///< Array accepting results of scattering in the 3rd dimension
  int* x_write_ovrf;            ///< Overflow arrays for scatter results in the 1st dimension,
                                ///<   when extended fixed-precision representations are in effect
  int* y_write_ovrf;            ///< Overflow arrays for scatter results in the 2nd dimension,
                                ///<   when extended fixed-precision representations are in effect
  int* z_write_ovrf;            ///< Overflow arrays for scatter results in the 3rd dimension,
                                ///<   when extended fixed-precision representations are in effect
};

/// \brief Collect pointers needed for conjugate gradient reduction operations, normalizing forces
///        and mixing the results with a separate vector containing some memory of previous
///        iterations.  This object re-purposes data in the PhaseSpaceSynthesis object such as the
///        velocities and prior position arrays.  Absent from this object is a sense of where each
///        system within the compiled synthesis starts and stops.
struct ConjGradSubstrate {

  /// \brief The constructor takes a PhaseSpaceSynthesis object and a series of double-precision
  ///        allocations.
  ///
  /// Overloaded:
  ///   - Accept a pointer to a coordinate synthesis object
  ///   - Accept coordinate synthesis abstract (passed by value)
  ///
  /// \param poly_psw  Writeable abstract to a coordinate synthesis
  /// \param poly_ps   Collection of coordinates in fixed-precision representation
  /// \param rbg       Allocations of double-precision reals to hold transitional sums between
  ///                  gather and scatter kernels
  /// \param tier      Get data pointers on the host (the customary default) or on the HPC device
  /// \{
  ConjGradSubstrate(PsSynthesisWriter poly_psw, ReductionBridge *rbg,
                    HybridTargetLevel tier = HybridTargetLevel::HOST);

  ConjGradSubstrate(PhaseSpaceSynthesis *poly_ps, ReductionBridge *rbg,
                    HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// Inverse force scaling constant.  This value is stored to provide a means for unrolling the
  /// fixed-precision scaling when computing squared gradients and gradient evolution quantities,
  /// although the double-precision accumulators are still very unlikely to overflow.
  const double inv_frc_scale;

  // The data pointers often repurpose arrays in the PhaseSpaceSynthesis object.
  llint* xfrc;          ///< Forces acting on all particles in the Cartesian X direction
  llint* yfrc;          ///< Forces acting on all particles in the Cartesian Y direction
  llint* zfrc;          ///< Forces acting on all particles in the Cartesian Z direction
  int* xfrc_ovrf;       ///< Overflow in X forces, when using extended precision models
  int* yfrc_ovrf;       ///< Overflow in Y forces, when using extended precision models
  int* zfrc_ovrf;       ///< Overflow in Z forces, when using extended precision models
  llint* xprv;          ///< Prior forces acting on all particles in the Cartesian X direction
  llint* yprv;          ///< Prior forces acting on all particles in the Cartesian Y direction
  llint* zprv;          ///< Prior forces acting on all particles in the Cartesian Z direction
  int* xprv_ovrf;       ///< Prior overflow in X forces, when using extended precision models
  int* yprv_ovrf;       ///< Prior overflow in Y forces, when using extended precision models
  int* zprv_ovrf;       ///< Prior overflow in Z forces, when using extended precision models
  llint* x_cg_temp;     ///< Mixing components for the X-direction conjugate gradient vector
  llint* y_cg_temp;     ///< Mixing components for the Y-direction conjugate gradient vector
  llint* z_cg_temp;     ///< Mixing components for the Z-direction conjugate gradient vector
  int* x_cg_temp_ovrf;  ///< Overflow for x_cg_temp when using extended precision models
  int* y_cg_temp_ovrf;  ///< Overflow for y_cg_temp when using extended precision models
  int* z_cg_temp_ovrf;  ///< Overflow for z_cg_temp when using extended precision models
  double* gg_buffer;    ///< Squared gradient partial sums
  double* dgg_buffer;   ///< Gradient evolution partial sums 
  double* msum_buffer;  ///< Squared sums of remixed forces on all atoms
};

} // namespace stmath
} // namespace stormm

#include "reduction_abstracts.tpp"

#endif
