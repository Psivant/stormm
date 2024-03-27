// -*-c++-*-
#ifndef STORMM_BRICKWORK_H
#define STORMM_BRICKWORK_H

#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using constants::CartesianDimension;
using constants::UnitCellAxis;
  
/// \brief Define one volume of whole-numbered side lengths within a larger unit cell, which itself
///        is part of a synthesis of many systems.
class VolumePartition {
public:

  /// \brief The constructor takes unitless lengths along the unit cell A, B, and C axes plus
  ///        an origin coordinate in the grid at which the volume partition begins.  In general,
  ///        volume partitions do not wrap but rather fill the primary unit cell volume exactly.
  VolumePartition(int a_dim_in = 0, int b_dim_in = 0, int c_dim_in = 0, int a_orig_in = 0,
                  int b_orig_in = 0, int c_orig_in = 0, int halo_under_in = 0,
                  int halo_over_in = 0, int system_index_in = 0);

  /// \brief With no pointers to repair, no const members, and any array data members taken from
  ///        the Standard Template Library, the default copy and move constructors, as well as
  ///        copy and move assignment operators, are appropriate.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  VolumePartition(const VolumePartition &original) = default;
  VolumePartition(VolumePartition &&original) = default;
  VolumePartition& operator=(const VolumePartition &original) = default;
  VolumePartition& operator=(VolumePartition &&original) = default;
  /// \}

  /// \brief Get the length along a specific dimension, or the lengths along all dimensions.
  ///
  /// Overloaded:
  ///   - Provide the unit cell axis designation or its Cartesian analog
  ///   - Get the length along one dimension or all dimensions
  ///
  /// \param dim  The unit cell axis along which to take the length
  /// \{
  int3 getLength() const;
  int getLength(UnitCellAxis dim) const;
  int getLength(CartesianDimension dim) const;
  /// \}

  /// \brief Get the length along a specific dimension, or the lengths along all dimensions, plus
  ///        the halo to give the true dimensions of the work involved with a particular brick.
  ///
  /// Overloaded:
  ///   - Provide the unit cell axis designation or its Cartesian analog
  ///   - Get the length along one dimension or all dimensions
  ///
  /// \param dim  The unit cell axis along which to take the length
  /// \{
  int3 getLengthPlusHalo() const;
  int getLengthPlusHalo(UnitCellAxis dim) const;
  int getLengthPlusHalo(CartesianDimension dim) const;
  /// \}

  /// \brief Get the origin of the volume partition.
  int3 getOrigin() const;

  /// \brief Get the total volume of the work unit, including the halo which also, presumably,
  ///        involves computing effort.
  int getVolumeWithHalo() const;
  
  /// \brief Get the system index to which the volume partition applies.
  int getSystemIndex() const;

  /// \brief Reduce the dimension of the work unit along the unit cell A axis.  A new partition
  ///        is returned with the volume that was cleaved off.
  ///
  /// \param a_dim_update  The new dimension to give the current cell
  VolumePartition split(int a_dim_update);
  
private:
  int a_dim;         ///< Length of the volume partition along the unit cell A axis
  int b_dim;         ///< Length of the volume partition along the unit cell B axis
  int c_dim;         ///< Length of the volume partition along the unit cell C axis
  int a_orig;        ///< Origin of the volume partition in the primary unit cell, of the system
                     ///<   to which it belongs, along the unit cell A axis
  int b_orig;        ///< Origin of the volume partition along the unit cell B axis
  int c_orig;        ///< Origin of the volume partition along the unit cell C axis
  int halo_under;    ///< Extent of the halo, in the same unitless dimensions as the side lengths,
                     ///<   in the positive direction along all axes
  int halo_over;     ///< Extent of the halo, in the same unitless dimensions as the side lengths,
                     ///<   in the negative direction along all axes
  int system_index;  ///< Index of the system within the overall synthesis to which this volume
                     ///<   partition applies
};

/// \brief This class will subdivide a series of simulatiion boxes, whose lengths are whole numbers
///        on all sides, into a set of partitions which keep to specified limits and cover the
///        entire space in a number of rectangular partitions as close as possible to a specified
///        number.  The partitions can be ordered by size.
class Brickwork {
public:

  /// \brief The minimum necessary information is a series of tuples indicating the dimensions of
  ///        each system.
  ///
  /// \param target_multiple        Work to bring the number of bricks as close to a multiple of
  ///                               this number as possible, without going over
  /// \param preferred_a_lengths    A list of preferred lengths for volume partitions along the
  ///                               unit cell A axis.
  /// \param discouraged_a_lengths  A list of discouraged lengths for volume partitions along the
  ///                               unit cell A axis
  Brickwork(const std::vector<int3> &system_dimensions_in = {}, int a_span_max_in = 7,
            int bc_cross_section_max_in = 16, int halo_under_in = 1, int halo_over_in = 0,
            int max_nonhalo_volume_in = 48, int target_multiple = 1,
            const std::vector<int> &preferred_a_lengths = {},
            const std::vector<int> &discouraged_a_lengths = {});

  /// \brief With no pointers to repair, no const members, and any array data members taken from
  ///        the Standard Template Library, the default copy and move constructors, as well as
  ///        copy and move assignment operators, are appropriate.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  Brickwork(const Brickwork &original) = default;
  Brickwork(Brickwork &&original) = default;
  Brickwork& operator=(const Brickwork &original) = default;
  Brickwork& operator=(Brickwork &&original) = default;
  /// \}

  /// \brief Get the number of bricks, the number of volume partitions spanning one system or all
  ///        systems.
  ///
  /// Overloaded:
  ///   - Get the number of bricks in the synthesis of all systems
  ///   - Get the number of bricks in a particular system
  ///
  /// \param system_index  The system of interest
  /// \{
  int getBrickCount() const;
  int getBrickCount(int system_index) const;
  /// \}

  /// \brief Get the number of individual systems.
  int getSystemCount() const;

  /// \brief Get the average dimension of all bricks along one unit cell A axis, with its
  ///        standard deviation.  The results are returned in the "x" and "y" members of the tuple,
  ///        respectively.
  ///
  /// Overloaded:
  ///   - Provide one of the unit cell axis desginations (A, B, C)
  ///   - Provide one of the Cartesian dimensions (x ~ A, y ~ B, z ~ C)
  ///
  /// \param dim  The unit cell axis along which to take the average dimension
  /// \{
  double2 getAverageBrickLength(UnitCellAxis dim) const;
  double2 getAverageBrickLength(CartesianDimension dim) const;
  /// \}

  /// \brief Get the origin of the volume partition, without consideration of the halo volume.
  ///
  /// \param brick_index  Index of the volume partition of interest
  int3 getBrickOrigin(int brick_index) const;

  /// \brief Get the origin of the volume partition, with the halo to the left along each axis
  ///        subtracted.  If this runs out of the primary unit cell, the origin coordinates will
  ///        be wrapped.
  ///
  /// \param brick_index  Index of the volume partition of interest
  int3 getBrickOriginWithHalo(int brick_index) const;

  /// \brief Get the lengths of a volume partition along all sides.
  ///
  /// \param brick_index  Index of the volume partition of interest
  int3 getBrickLengths(int brick_index) const;

  /// \brief Get the lengths of a volume partition along all sides, with accounting for its halo.
  ///
  /// \param brick_index  Index of the volume partition of interest
  int3 getBrickLengthsWithHalo(int brick_index) const;

  /// \brief Get the system index to which a particular brick pertains.
  ///
  /// \param brick_index  Index of the volume partition of interest
  int getSystemMembership(int brick_index) const;
  
  /// \brief Get the system dimensions
  ///
  /// \param system_index  The system of interest
  int3 getSystemDimensions(int system_index) const;
  
  /// \brief Subdivide all systems into a series of bricks.  This starts by creating the fewest
  ///        possible bricks and then further subdividing until the target number, or something
  ///        near to a multiple thereof, is hit.
  ///
  /// \param target_work_unit_multiple  The ideal number of work units to create
  /// \param preferred_a_lengths        A list of preferred lengths for volume partitions along
  ///                                   the unit cell A axis
  /// \param discouraged_a_lengths      A list of discouraged lengths for volume partitions along
  ///                                   the unit cell A axis
  void subdivide(int target_work_unit_multiple = 1,
                 const std::vector<int> &preferred_a_lengths = {},
                 const std::vector<int> &discouraged_a_lengths = {});

private:
  
  // Store the parameters and metrics for later reference.
  int a_span_max;            ///< Maximum allowed dimension of any one brick along the unit cell
                             ///<   A axis, including its halos over and under.
  int bc_cross_section_max;  ///< Maximum allowed cross section of any brick in the BC plane
                             ///<   (the product of its lengths along the unit cell B and C axes,
                             ///<   including their halos over and under)
  int halo_under;            ///< The unitless number of grid elements associated with any brick,
                             ///<   to the left of the brick's origin along any one axis.  Unlike
                             ///<   the limits of the brick itself, the halo can wrap inside the
                             ///<   unit cell.  A brick starting at { 0, 0, 0 } within a unit cell
                             ///<   of lengths { M, N, P } can have a halo beginning at
                             ///<   { M - 1, N - 1, P - 1 }.
  int halo_over;             ///< The unitless number of grid elements associated with any brick,
                             ///<   to the right of the brick's origin along any one axis.
  int max_nonhalo_volume;    ///< The maximum size of each work unit's volume, excluding that
                             ///<   covered by the halo.  For example, in charge mapping, the halo
                             ///<   region might be one unit underneath the non-halo volume's
                             ///<   footprint along all three unit cell axes.  The maximum non-halo
                             ///<   volume may be essential to ensuring that certain memory limits
                             ///<   are not overrun by the work unit.
  
  // Store the overall dimensions of each system's grid
  std::vector<int3> system_dimensions;

  // An array of volume partitions stores the manner in which everything is broken up
  std::vector<VolumePartition> bricks;

  /// \brief Validate a particular system index.
  ///
  /// \param system_index  Index of the system of interest
  /// \param caller        Name of the calling function (for error back-tracing)
  void validateSystemIndex(int system_index, const char* caller) const;
  
  /// \brief Validate a particular brick index.
  ///
  /// \param index   Index of the volume partition of interest.
  /// \param caller   Name of the calling function (for error back-tracing)
  void validateBrickIndex(int index, const char* caller) const;
  
  /// \brief Count the minimal number of work units results under a certain plan for subdividing
  ///        each unit cell volume with given target lengths along the B and C dimensions.
  ///
  /// \param trimmed_bdim  Target length for work units (less halos) along the B unit cell axis
  /// \param trimmed_cdim  Target length for work units (less halos) along the C unit cell axis
  int countMinimalWorkUnits(int trimmed_bdim, int trimmed_cdim) const;
};
  
} // namespace synthesis
} // namespace stormm

#endif
