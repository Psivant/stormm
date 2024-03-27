// -*-c++-*-
#ifndef STORMM_REDUCTION_WORKUNIT_H
#define STORMM_REDUCTION_WORKUNIT_H

#include <cmath>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace stmath {

using card::GpuDetails;
using card::Hybrid;
  
/// \brief The maximum number of slots into which gathering-related kernels can funnel their
///        results.
constexpr int maximum_gathering_results = 1024;

/// \brief Define the abstract length for a ReductionWorkUnit's abstract.  This must be a multiple
///        of 2, and less than the warp width (as small as 16 on Intel GPUs).
constexpr int rdwu_abstract_length = 8;

/// \brief A work unit to describe the manner in which groups of atoms in each structure of a
///        synthesis come together to contribute to a single result.  Reduction work units serve
///        one and only one system apiece.  The maximum size of these work units can be specified
///        when they are created, in a manner analogous to valence or non-bonded work units, but
///        is not limited by space in the GPU L1 cache so that, for practical optimizations,
///        gather and scatter operations can be combined into all-reduce operations with a single
///        kernel launch.
class ReductionWorkUnit {
public:

  /// \brief The constructor takes arguments for all members.
  ReductionWorkUnit(int atom_start_in, int atom_end_in, int result_index_in,
                    int dependency_start_in, int dependency_end_in, int system_index_in);

  /// \brief Take all copy and move constructors and assignment operators.
  /// \{
  ReductionWorkUnit(const ReductionWorkUnit &original) = default;
  ReductionWorkUnit(ReductionWorkUnit &&original) = default;
  ReductionWorkUnit& operator=(const ReductionWorkUnit &other) = default;
  ReductionWorkUnit& operator=(ReductionWorkUnit &&other) = default;
  /// \}
  
  /// \brief Get the atom starting index.
  int getAtomStart() const;

  /// \brief Get the upper limit of atoms in this work unit.
  int getAtomEnd() const;

  /// \brief Get the index of whatever result array where this work unit will put its result.
  int getResultIndex() const;

  /// \brief Get the start of dependencies in the result array which pertain to the same system as
  ///        this reduction work unit.  All reduction work units serving the same system will
  ///        contribute their gathering results to contiguous elements of whatever result arrays.
  int getDependencyStart() const;

  /// \brief Get the upper limit of dependencies in the result array which pertain to the same
  ///        system as this reduction work unit.
  int getDependencyEnd() const;

  /// \brief Get the system to which this reduction work unit pertains (each reduction work unit
  ///        will serve one and only one system in a synthesis)
  int getSystemIndex() const;

  /// \brief Produce an abstract containing all of the information, wrapped in a series of eight
  ///        integers.  See the enum class RdwuAbstractMap to understand which element represents
  ///        which value.
  std::vector<int> getAbstract() const;
  
private:
  int atom_start;        ///< Lower limit of atoms in the unified synthesis array which contribute
                         ///<   to this work unit's reduction operation
  int atom_end;          ///< Upper limit of atoms in the unified synthesis array which contribute
                         ///<   to this work unit's reduction operation
  int result_index;      ///< Index in the accumulation arrays to which results should be written
                         ///<   (if a reduction across two kernels is needed at all)
  int dependency_start;  ///< Staring index, in the accumulation arrays, of results from reduction
                         ///<   work units pertaining to the same system
  int dependency_end;    ///< Upper bounding index, in the accumulation arrays, of results from
                         ///<   reduction work units pertaining to the same system
  int system_index;      ///< The system to which the result pertains (this may or may not be
                         ///<   necessary, given the availability of bounds in atom_start and
                         ///<   atom_end)
};

/// \brief Find the optimal subdivision of various reduction kernel blocks given a series of
///        system sizes.  The reduction kernels have a maximum thread count of small_block_size
///        (256, see Constants/hpc_bounds.h) and are divisible up to two times.
///
/// Overloaded:
///   - Provide the system sizes as a C-style array with a trusted length.
///   - Provide the system sizes as a Standard Template Library vector.
///   - Provide the system sizes as a Hybrid object (host data will be valid).
///
/// \param atom_counts  The numbers of atoms in each system
/// \param n_systems    Number of systems (if a C-style array of the atom counts is provided)
/// \param gpu          Details of the GPU to use in calculations
/// \{
int optReductionKernelSubdivision(const int* atom_counts, int n_systems, const GpuDetails &gpu);
int optReductionKernelSubdivision(const std::vector<int> &atom_counts, const GpuDetails &gpu);
int optReductionKernelSubdivision(const Hybrid<int> &atom_counts, const GpuDetails &gpu);
/// \}

/// \brief Build the reduction (and their components of gathering and scattering) work units for
///        a series of systems of stated sizes.  Only the starting indices and system sizes are
///        needed, and only those are provided (rather than the AtomGraphSynthesis itself, from
///        which these arrays are probably derived), to avoid circular dependencies.
///
/// \param atom_starts     Starting indices for atoms of each system in the collection (the length
///                        of this array will determine the number of systems--provide a vector
///                        with a single entry of zero to get reduction work units for a single
///                        system)
/// \param atom_counts     Atom counts for all systems in the collection
/// \param gpu             Details of the selected GPU to be used for calculations
/// \param launcher        Object to collect wisdom about optimal kernel launch configurations with
///                        the work unit array chosen to suit the collection of systems
/// \param tasks_per_atom  The number of values to reduce across each system.  A center of geometry
///                        computation would have three values (X, Y, and Z), whereas a total
///                        charge summation with a mask would have only one.
std::vector<ReductionWorkUnit> buildReductionWorkUnits(const std::vector<int> &atom_starts,
                                                       const std::vector<int> &atom_counts,
                                                       const GpuDetails &gpu,
                                                       int tasks_per_atom = 1);

} // namespace stmath
} // namespace stormm

#endif

