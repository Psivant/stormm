// -*-c++-*-
#ifndef STORMM_PPI_TABLE_H
#define STORMM_PPI_TABLE_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Math/log_scale_spline.h"
#include "Math/math_enumerators.h"
#include "energy_enumerators.h"
#include "pme_util.h"

namespace stormm {
namespace energy {
  
using constants::ExceptionResponse;
using constants::PrecisionModel;
using symbols::amber_ancient_bioq;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using stmath::BasisFunctions;
using stmath::LogScaleSpline;
using stmath::LogSplineForm;
using stmath::LogSplineTable;
using stmath::TableIndexing;

/// \brief Abstract for the particle-particle interaction table, with pointers to energies and
///        forces and sizing constants for navigating the excluded versus non-excluded
///        interactions.  Like the LogScaleSpline abstract, this is read-only.
template <typename T4> class PPIKit {

  /// \brief The constructor takes the customary list of pointers and critical constants.
  PPIKit(NonbondedTheme theme_in, BasisFunctions basis_in, TableIndexing lookup_in,
         int index_bound_in, int excl_offset_in, uint sp_detail_mask_in, ullint dp_detail_mask_in,
         float idx_offset_in, const T4* energy_in, const T4* force_in, const T4* energy_excl_in,
         const T4* force_excl_in);

  /// \brief With all const members, the copy and move assignment operators are implicitly
  ///        deleted, but with no pointers to repair the copy and move constructors are legal.
  ///
  /// \param original  The original PPIKit tocopy or move
  /// \{
  PPIKit(const PPIKit &original) = default;
  PPIKit(PPIKit &&original) = default;
  /// \}

  const NonbondedTheme theme;   ///< The type of non-bonded potential and derivatives approximated
                                ///<   by the splines
  const BasisFunctions basis;   ///< The type of basis functions scaled by spline coefficients to
                                ///<   obtain the result
  const TableIndexing lookup;   ///< The method of transforming the argument to the underlying
                                ///<   function
  const int index_bound;        ///< The upper bound of spline indices in any given table.  This
                                ///<   is half the value of excl_offset.
  const int excl_offset;        ///< The offset to add to any spline index in the energy or force
                                ///<   arrays in order to obtain the corresponding excluded energy
                                ///<   or force spline coefficients
  const uint sp_detail_mask;    ///< A pre-computed bitmask useful for converting an argument of
                                ///<   the underlying benhmark function into a quantity for
                                ///<   single-precision splines to operate on
  const ullint dp_detail_mask;  ///< A pre-computed bitmask useful for converting an argument of
                                ///<   the underlying benhmark function into a quantity for
                                ///<   double-precision splines to operate on
  const float idx_offset;       ///< Offset used in transforming arguments of the underlying
                                ///<   benchmark function into spline table indices
  const T4* energy;             ///< Table of energy coefficients.  Add excl_offset to any index
                                ///<   for the raw energy to get the corresponding index in the
                                ///<   excluded energy table.  Add index_bound to any index for the
                                ///<   raw energy to get the raw force.
  const T4* force;              ///< Table of force coefficients.  Add excl_offset to any index
                                ///<   for the raw force to get the corresponding index in the
                                ///<   excluded force table.
  const T4* energy_excl;        ///< Table of excluded energy coefficients
  const T4* force_excl;         ///< Table of excluded force coefficients
};
  
/// \brief A tabulated non-bonded potential, with or without exclusions, to be used in the context
///        of particle-particle, particle-mesh calculations.  The key is to create two tables, one
///        for the non-excluded form of the interaction and the other for the excluded form.  The
///        entries for each table will then be concatenated, such that all non-excluded
///        interactions are contiguous and then all exclude interactions are contiguous.  The
///        offset for accessing an excluded interaction based on an index calculated from a
///        particle-particle distance is stored alongside the tabulated splines in the abstract.
template <typename T4> class PPITable {
public:

  /// \brief The constuctor can accept all of the arguments that might be useful for making a
  ///        LosScaleSpline, or a LogScaleSpline itself.  Tables for both the otential and the
  ///        derivative will be computed.
  /// \{
  PPITable(NonbondedTheme theme_in = NonbondedTheme::ELECTROSTATIC,
           BasisFunctions basis_set_in = BasisFunctions::POLYNOMIAL,
           TableIndexing indexing_method_in = TableIndexing::SQUARED_ARG,
           double cutoff_in = default_pme_cutoff, double dsum_tol_in = default_dsum_tol,
           int mantissa_bits_in = 5, double coulomb_in = amber_ancient_bioq,
           double min_spl_compute_range_in = 0.015625, double min_offset_in = 0.125);

  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
           const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
           double cutoff_in = default_pme_cutoff);
  
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
           const LogScaleSpline<T4> &spl_c);

  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b);

  PPITable(const LogScaleSpline<T4> &spl_a);
  /// \}
  
  /// \brief The presence of POINTER-kind Hybrid objects implies pointer repairs that require the
  ///        copy and move constructors as well as assignment operators to be spelled out.
  ///
  /// \param original  THe original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  PPITable(const PPITable<T4> &original);
  PPITable(PPITable<T4> &&original);
  PPITable& operator=(const PPITable<T4> &original);
  PPITable& operator=(PPITable<T4> &&original);
  /// \}

  /// \brief Get the non-bonded potential described by an object of this class.
  NonbondedTheme getTheme() const;

  /// \brief Get the cutoff on particle-particle interactions.
  double getCutoff() const;

  /// \brief Get the maximum range of the spline.
  double getMaximumRange() const;

  /// \brief Get the direct sum tolerance.
  double getDirectSumTolerance() const;

  /// \brief Get the number of bits of the mantissa used for table indexing.
  int getBitStride() const;

  /// \brief Get the Ewald coefficient used to perform the switching between short- and
  ///        long-ranged potentials.
  double getEwaldCoefficient() const;

  /// \brief Get the Gaussian RMS sigma parameter (one half the inverse of the Ewald coefficient)
  ///        used to perform the switching between short- and long-ranged potentials.
  double getGaussianWidth() const;

  /// \brief Get the abstract for use of the splines in a manner like C programming.
  ///
  /// \param tier  Indicate whether to obtain pointers for data on the CPU host or GPU device
  const PPIKit<T4> data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get an abstract for the object with its templated character cast away to void.
  const PPIKit<void> temlateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

private:

  NonbondedTheme theme;           ///< The type of non-bonded potential encoded in the tabulated
                                  ///<   coefficients.
  BasisFunctions basis_set;       ///< The type of basis functions to use in spline approximations
  TableIndexing indexing_method;  ///< The indexing method that converts the argument of the
                                  ///<   underlying function into a spline table element
  double cutoff;                  ///< Maximum inter-particle distance for which interactions are
                                  ///<   valid
  double max_range;               ///< The maximum interparticle distance for which the table can
                                  ///<   produce a result
  double dsum_tol;                ///< The direct sum tolerance, a proportion of each interaction
                                  ///<   that is discarded beginning at the particle-particle pair
                                  ///<   cutoff
  double ew_coeff;                ///< The Ewald coefficient that governs the potential splitting
                                  ///<   between short- and long-ranged components
  int mantissa_bits;              ///< Number of mantissa bits used to index into various tables
  double coulomb;                 ///< The object's take on Coulomb's constant
  int exclusion_offset;           ///< The first index of the table of excluded interactions.  The
                                  ///<   series of non-excluded interactions appears first in each
                                  ///<   table, whether of potentials or derivatives.

  // The main arrays, in full (double) precision
  Hybrid<T4> energy;     ///< Double-precision coefficients for the energy, without exclusions
  Hybrid<T4> force;      ///< Double-precision coefficients for the force, without exclusions

  /// \brief Double-precision coefficients for the energy, with excluded 1:4 interactions
  Hybrid<T4> energy_with_exclusions;

  /// \brief Double-precision coefficients for the force, with excluded 1:4 interactions
  Hybrid<T4> force_with_exclusions;

  // Storage for all coefficients (these are ARRAY-kind Hybrids, unlike the POINTER-kind Hybrids
  // above)
  Hybrid<T4> coeffs;     ///< Coefficients of the tables computed in double precision.  This is an
                         ///<   ARRAY-kind Hybrid targeted by the  double-precision POINTER-kind
                         ///<   Hybrid arrays above.

  /// \brief Find a potential function, without exclusions, among the spline tables provided.
  ///
  /// \param spl_a  The first spline table provided
  /// \param spl_b  The second spline table provided
  /// \param spl_c  The third spline table provided
  /// \param spl_d  The fourth spline table provided
  const LogScaleSpline<T4>& findNonExclPotential(const LogScaleSpline<T4> &spl_a,
                                                 const LogScaleSpline<T4> &spl_b,
                                                 const LogScaleSpline<T4> &spl_c,
                                                 const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a potential function, with non-bonded 1:4 exclusions, among the spline tables
  ///        provided.  Descriptions of input parameters follow from findNonExclPotential() above.
  const LogScaleSpline<T4>& findExclPotential(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b,
                                              const LogScaleSpline<T4> &spl_c,
                                              const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a force function, without exclusions, among the spline tables provided.
  ///        Descriptions of input parameters follow from findNonExclPotential() above.
  const LogScaleSpline<T4>& findNonExclForce(const LogScaleSpline<T4> &spl_a,
                                             const LogScaleSpline<T4> &spl_b,
                                             const LogScaleSpline<T4> &spl_c,
                                             const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a force function, with non-bonded 1:4 exclusions, among the spline tables
  ///        provided.  Descriptions of input parameters follow from findNonExclPotential() above.
  const LogScaleSpline<T4>& findExclForce(const LogScaleSpline<T4> &spl_a,
                                          const LogScaleSpline<T4> &spl_b,
                                          const LogScaleSpline<T4> &spl_c,
                                          const LogScaleSpline<T4> &spl_d) const;

  /// \brief Check a bit in a mask to signify the presence of a particular, necessary function.
  ///
  /// \param spl_x  A spline containing one function relevant to the PPITable
  uint checkPriority(const LogScaleSpline<T4> &spl_x) const;
  
  /// \brief Determine the highest-priority missing functional form.
  ///
  /// \param holdings  A bitmask checked 1 for the availability of relevant functions: 0x1 for
  ///                  the potential, 0x2 for the excluded potential, 0x4 for the force, 0x8 for
  ///                  the excluded force
  LogSplineForm findMissingForm(const uint holdings) const;
  
  /// \brief Construct the necessary splines for this aggregate table, based on one to three
  ///        other splines which specify the type.  The priority system is first to find the
  ///        non-excluded energy, then the excluded energy, then the non-excluded force, and
  ///        finally the excluded force (between particles).  Various overloads work given up to
  ///        three available logarithmic spline tables.
  ///
  /// \param spl_a  The first of up to three known spline tables
  /// \param spl_b  The second of up to three known spline tables
  /// \param spl_c  The third and last known spline table
  /// \{
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a,
                                      const LogScaleSpline<T4> &spl_b,
                                      const LogScaleSpline<T4> &spl_c) const;

  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a,
                                       const LogScaleSpline<T4> &spl_b) const;

  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a) const;
  /// \}

  /// \brief Fill the coefficients table with the energy, force, excluded energy, and excluded
  ///        force splines.
  ///
  /// \param u    The energy spline table
  /// \param du   The force spline table
  /// \param ux   The excluded energy spline table
  /// \param dux  The excluded force spline table
  void populateCoefficients(const LogScaleSpline<T4> &u, const LogScaleSpline<T4> &du,
                            const LogScaleSpline<T4> &ux, const LogScaleSpline<T4> &dux)
};

/// \brief Restore the type of a particle-particle interaction kit, the abstract of the PPITable
///        class above.
///
/// Overloaded:
///   - Provide a const pointer to the template-free abstract
///   - Provide a const reference the template-free abstract
///
/// \param rasa  The original abstract with all templated pointers cast to void, a "blank slate"
/// \{
template <typename T4> const PPIKit<T4> restoreType(const PPIKit<void> *rasa);
template <typename T4> const PPIKit<T4> restoreType(const PPIKit<void> &rasa);
/// \}

} // namespace energy
} // namespace stormm

#endif
