// -*-c++-*-
#ifndef STORMM_CACHERESOURCE_H
#define STORMM_CACHERESOURCE_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief Abstract for the CacheResource object, accessible as a C-style struct and suitable for
///        passing to GPU kernels as a kernel argument.
template <typename T> struct CacheResourceKit {

  /// \brief The constructor for this abstract takes all of the usual pointer arguments to fill
  ///        out its member variables, the template parameter corresponding to the precision model
  ///        of the charges.
  CacheResourceKit(int max_blocks_in, int max_atoms_in, llint* xcrd_in, llint* ycrd_in,
                   llint* zcrd_in, llint* xvel_in, llint* yvel_in, llint* zvel_in,
                   int* xcrd_ovrf_in, int* ycrd_ovrf_in, int* zcrd_ovrf_in, int* xvel_ovrf_in,
                   int* yvel_ovrf_in, int* zvel_ovrf_in, int* xfrc_ovrf_in, int* yfrc_ovrf_in,
                   int* zfrc_ovrf_in, T* charges_in, int* lj_idx_in);

  /// \brief The usual copy and move constructors for an abstract with one or more const member
  ///        variables apply.
  /// \{
  CacheResourceKit(const CacheResourceKit &original) = default;
  CacheResourceKit(CacheResourceKit &&original) = default;
  /// \}
  
  const int max_blocks;  ///< Maximum number of blocks that resources are designed to accommodate
  const int max_atoms;   ///< Maximum number of atoms that resources are designed to accommodate.
                         ///<   The product of this and max_blocks (above) is what really matters.
  llint* xcrd;           ///< Cartesian X coordinates of locally cached particles
  llint* ycrd;           ///< Cartesian Y coordinates of locally cached particles
  llint* zcrd;           ///< Cartesian Z coordinates of locally cached particles
  llint* xvel;           ///< Cartesian X velocities of locally cached particles
  llint* yvel;           ///< Cartesian Y velocities of locally cached particles
  llint* zvel;           ///< Cartesian Z velocities of locally cached particles
  int* xcrd_ovrf;        ///< Cartesian X coordinate overflow buffers
  int* ycrd_ovrf;        ///< Cartesian Y coordinate overflow buffers
  int* zcrd_ovrf;        ///< Cartesian Z coordinate overflow buffers
  int* xvel_ovrf;        ///< Cartesian X velocity overflow buffers
  int* yvel_ovrf;        ///< Cartesian Y velocity overflow buffers
  int* zvel_ovrf;        ///< Cartesian Z velocity overflow buffers
  int* xfrc_ovrf;        ///< Cartesian X force overflow buffers
  int* yfrc_ovrf;        ///< Cartesian Y force overflow buffers
  int* zfrc_ovrf;        ///< Cartesian Z force overflow buffers
  T* charges;            ///< Charge parameters for locally cached particles (non-const
                         ///<   anticipating some future support of a polarizable model)
  int* lj_idx;           ///< Lennard-Jones indices of locally cached particles
};

/// \brief An object to hold temporary data for a particular work unit (whether bonded or
///        non-bonded), resident in GMEM but private to a particular thread block.  This object
///        must be allocated in such a way as to be ready to hold the private workspaces of any
///        thread blocks that will make use of it.  The object is arranged to support operations
///        based on PhaseSpaceSynthesis respresentations of systems' coordinates and forces.
class CacheResource {
public:

  /// \brief The constructor does not take a GPU description, but instead maximum numbers of
  ///        blocks and atoms per block that might be required, allocating space as appropriate.
  ///
  /// \param block_limit_in  The maximum number of thread blocks that will need resources
  /// \param atom_limit_in   The maximum number of atoms per block requiring local copies
  CacheResource(int block_limit_in, int atom_limit_in);

  /// \brief Basic copy and move constructors
  ///
  /// \param original  The object to copy
  /// \{
  CacheResource(const CacheResource &original);
  CacheResource(CacheResource &&original);
  /// \}

  /// \brief Basic copy and move assignment operators
  ///
  /// \param other  The object to assign against
  /// \{
  CacheResource& operator=(const CacheResource &other);
  CacheResource& operator=(CacheResource &&other);
  /// \}

  /// \brief Get a set of pointers to this object with double-precision representations for the
  ///        charges.
  CacheResourceKit<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Get a set of pointers to this object with single-precision representations for the
  ///        charges.
  CacheResourceKit<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST);

private:
  int block_limit;                    ///< The number of blocks that the resources are designed
                                      ///<   to accommodate
  int atom_limit;                     ///< The maximum number of atoms each block's resources are
                                      ///<   intended to accommodate (half as many blocks could
                                      ///<   also accommodate twice as many atoms)
  Hybrid<llint> x_coordinates;        ///< Cartesian X coordinates of locally cached particles
  Hybrid<llint> y_coordinates;        ///< Cartesian Y coordinates of locally cached particles
  Hybrid<llint> z_coordinates;        ///< Cartesian Z coordinates of locally cached particles
  Hybrid<llint> x_velocities;         ///< Cartesian X velocities of locally cached particles
  Hybrid<llint> y_velocities;         ///< Cartesian Y velocities of locally cached particles
  Hybrid<llint> z_velocities;         ///< Cartesian Z velocities of locally cached particles
  Hybrid<int> x_coordinate_overflow;  ///< Cartesian X coordinate overflow buffers
  Hybrid<int> y_coordinate_overflow;  ///< Cartesian Y coordintae overflow buffers
  Hybrid<int> z_coordinate_overflow;  ///< Cartesian Z coordinate overflow buffers
  Hybrid<int> x_velocity_overflow;    ///< Cartesian X velocity overflow buffers
  Hybrid<int> y_velocity_overflow;    ///< Cartesian Y velocity overflow buffers
  Hybrid<int> z_velocity_overflow;    ///< Cartesian Z velocity overflow buffers
  Hybrid<int> x_force_overflow;       ///< Cartesian X force overflow buffers
  Hybrid<int> y_force_overflow;       ///< Cartesian Y force overflow buffers
  Hybrid<int> z_force_overflow;       ///< Cartesian Z force overflow buffers
  Hybrid<double> charges;             ///< Charge parameters for locally cached particles
  Hybrid<float> sp_charges;           ///< Charge parameters for locally cached particles (single
                                      ///<   precision)
  Hybrid<int> lennard_jones_indices;  ///< Lennard-Jones indices of locally cached particles  
  Hybrid<int> int_data;               ///< Storage array targeted by all other POINTER-kind int
                                      ///<   Hybrid member variables in this object
  Hybrid<llint> llint_data;           ///< Storage array targeted by all other POINTER-kind llint
                                      ///<   Hybrid member variables in this object

  /// \brief Reset the POINTER-kind Hybrid objects of an object that has just been copied, making
  ///        them target the object's own ARRAY-kind storage rather than that of some original
  ///        object.
  void repairPointers();
};

} // namespace energy
} // namespace stormm

#include "cacheresource.tpp"

#endif
