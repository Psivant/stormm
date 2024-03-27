#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Math/rounding.h"
#include "cacheresource.h"

namespace stormm {
namespace energy {

using card::HybridKind;
using stmath::roundUp;
  
//-------------------------------------------------------------------------------------------------
CacheResource::CacheResource(const int block_limit_in, const int atom_limit_in) :
    block_limit{block_limit_in},
    atom_limit{roundUp(atom_limit_in, warp_size_int)},
    x_coordinates{HybridKind::POINTER, "cache_xcrd"},
    y_coordinates{HybridKind::POINTER, "cache_ycrd"},
    z_coordinates{HybridKind::POINTER, "cache_zcrd"},
    x_velocities{HybridKind::POINTER, "cache_xvel"},
    y_velocities{HybridKind::POINTER, "cache_yvel"},
    z_velocities{HybridKind::POINTER, "cache_zvel"},
    x_coordinate_overflow{HybridKind::POINTER, "cache_xcrd_ovrf"},
    y_coordinate_overflow{HybridKind::POINTER, "cache_ycrd_ovrf"},
    z_coordinate_overflow{HybridKind::POINTER, "cache_zcrd_ovrf"},
    x_velocity_overflow{HybridKind::POINTER, "cache_xvel_ovrf"},
    y_velocity_overflow{HybridKind::POINTER, "cache_yvel_ovrf"},
    z_velocity_overflow{HybridKind::POINTER, "cache_zvel_ovrf"},
    x_force_overflow{HybridKind::POINTER, "cache_xfrc_ovrf"},
    y_force_overflow{HybridKind::POINTER, "cache_yfrc_ovrf"},
    z_force_overflow{HybridKind::POINTER, "cache_zfrc_ovrf"},
    charges{static_cast<size_t>(block_limit * atom_limit), "cache_charges"},
    sp_charges{static_cast<size_t>(block_limit * atom_limit), "sp_cache_charges"},
    lennard_jones_indices{HybridKind::POINTER, "cache_ljidx"},
    int_data{static_cast<size_t>(10 * block_limit * atom_limit), "cache_int_data"},
    llint_data{static_cast<size_t>(6 * block_limit * atom_limit), "cache_llint_data"}
{
  if (block_limit <= 0 || atom_limit <= 0) {
    rtErr("Device resources private to individual blocks cannot be allocated for " +
          std::to_string(block_limit) + " blocks with " + std::to_string(atom_limit) + " atoms.",
          "CacheResource");
  }
  const int per_item = block_limit * atom_limit;
  x_coordinates.setPointer(&llint_data,            0, per_item);
  y_coordinates.setPointer(&llint_data,     per_item, per_item);
  z_coordinates.setPointer(&llint_data, 2 * per_item, per_item);
  x_velocities.setPointer(&llint_data,  3 * per_item, per_item);
  y_velocities.setPointer(&llint_data,  4 * per_item, per_item);
  z_velocities.setPointer(&llint_data,  5 * per_item, per_item);
  x_coordinate_overflow.setPointer(&int_data,            0, per_item);
  y_coordinate_overflow.setPointer(&int_data,     per_item, per_item);
  z_coordinate_overflow.setPointer(&int_data, 2 * per_item, per_item);
  x_velocity_overflow.setPointer(&int_data,   3 * per_item, per_item);
  y_velocity_overflow.setPointer(&int_data,   4 * per_item, per_item);
  z_velocity_overflow.setPointer(&int_data,   5 * per_item, per_item);
  x_force_overflow.setPointer(&int_data,      6 * per_item, per_item);
  y_force_overflow.setPointer(&int_data,      7 * per_item, per_item);
  z_force_overflow.setPointer(&int_data,      8 * per_item, per_item);
  lennard_jones_indices.setPointer(&int_data, 9 * per_item, per_item);
}

//-------------------------------------------------------------------------------------------------
CacheResource::CacheResource(const CacheResource &original) :
    block_limit{original.block_limit},
    atom_limit{original.atom_limit},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    x_velocities{original.x_velocities},
    y_velocities{original.y_velocities},
    z_velocities{original.z_velocities},
    x_coordinate_overflow{original.x_coordinate_overflow},
    y_coordinate_overflow{original.y_coordinate_overflow},
    z_coordinate_overflow{original.z_coordinate_overflow},
    x_velocity_overflow{original.x_velocity_overflow},
    y_velocity_overflow{original.y_velocity_overflow},
    z_velocity_overflow{original.z_velocity_overflow},
    x_force_overflow{original.x_force_overflow},
    y_force_overflow{original.y_force_overflow},
    z_force_overflow{original.z_force_overflow},
    charges{original.charges},
    sp_charges{original.sp_charges},
    lennard_jones_indices{original.lennard_jones_indices},
    int_data{original.int_data},
    llint_data{original.llint_data}
{
  repairPointers();
}

//-------------------------------------------------------------------------------------------------
CacheResource::CacheResource(CacheResource &&original) :
    block_limit{original.block_limit},
    atom_limit{original.atom_limit},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    x_velocities{std::move(original.x_velocities)},
    y_velocities{std::move(original.y_velocities)},
    z_velocities{std::move(original.z_velocities)},
    x_coordinate_overflow{std::move(original.x_coordinate_overflow)},
    y_coordinate_overflow{std::move(original.y_coordinate_overflow)},
    z_coordinate_overflow{std::move(original.z_coordinate_overflow)},
    x_velocity_overflow{std::move(original.x_velocity_overflow)},
    y_velocity_overflow{std::move(original.y_velocity_overflow)},
    z_velocity_overflow{std::move(original.z_velocity_overflow)},
    x_force_overflow{std::move(original.x_force_overflow)},
    y_force_overflow{std::move(original.y_force_overflow)},
    z_force_overflow{std::move(original.z_force_overflow)},
    charges{std::move(original.charges)},
    sp_charges{std::move(original.sp_charges)},
    lennard_jones_indices{std::move(original.lennard_jones_indices)},
    int_data{std::move(original.int_data)},
    llint_data{std::move(original.llint_data)}
{}

//-------------------------------------------------------------------------------------------------
CacheResource& CacheResource::operator=(const CacheResource &other) {

  // Guard against self-assignment.  Otherwise, copy the necessary elements.
  if (this == &other) {
    return *this;
  }
  block_limit = other.block_limit;
  atom_limit = other.atom_limit;
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  x_velocities = other.x_velocities;
  y_velocities = other.y_velocities;
  z_velocities = other.z_velocities;
  x_coordinate_overflow = other.x_coordinate_overflow;
  y_coordinate_overflow = other.y_coordinate_overflow;
  z_coordinate_overflow = other.z_coordinate_overflow;
  x_velocity_overflow = other.x_velocity_overflow;
  y_velocity_overflow = other.y_velocity_overflow;
  z_velocity_overflow = other.z_velocity_overflow;
  x_force_overflow = other.x_force_overflow;
  y_force_overflow = other.y_force_overflow;
  z_force_overflow = other.z_force_overflow;
  charges = other.charges;
  sp_charges = other.sp_charges;
  lennard_jones_indices = other.lennard_jones_indices;
  int_data = other.int_data;
  llint_data = other.llint_data;

  // Repair the pointers and return the object
  repairPointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
CacheResource& CacheResource::operator=(CacheResource &&other) {

  // Guard against self-assignment.  Otherwise, move the necessary elements.
  if (this == &other) {
    return *this;
  }
  block_limit = other.block_limit;
  atom_limit = other.atom_limit;
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  x_velocities = std::move(other.x_velocities);
  y_velocities = std::move(other.y_velocities);
  z_velocities = std::move(other.z_velocities);
  x_coordinate_overflow = std::move(other.x_coordinate_overflow);
  y_coordinate_overflow = std::move(other.y_coordinate_overflow);
  z_coordinate_overflow = std::move(other.z_coordinate_overflow);
  x_velocity_overflow = std::move(other.x_velocity_overflow);
  y_velocity_overflow = std::move(other.y_velocity_overflow);
  z_velocity_overflow = std::move(other.z_velocity_overflow);
  x_force_overflow = std::move(other.x_force_overflow);
  y_force_overflow = std::move(other.y_force_overflow);
  z_force_overflow = std::move(other.z_force_overflow);
  charges = std::move(other.charges);
  sp_charges = std::move(other.sp_charges);
  lennard_jones_indices = std::move(other.lennard_jones_indices);
  int_data = std::move(other.int_data);
  llint_data = std::move(other.llint_data);

  // No further pointer repair is needed
  return *this;
}

//-------------------------------------------------------------------------------------------------
CacheResourceKit<double> CacheResource::dpData(const HybridTargetLevel tier) {
  return CacheResourceKit<double>(block_limit, atom_limit, x_coordinates.data(tier),
                                  y_coordinates.data(tier), z_coordinates.data(tier),
                                  x_velocities.data(tier), y_velocities.data(tier),
                                  z_velocities.data(tier), x_coordinate_overflow.data(tier),
                                  y_coordinate_overflow.data(tier),
                                  z_coordinate_overflow.data(tier), x_velocity_overflow.data(tier),
                                  y_velocity_overflow.data(tier), z_velocity_overflow.data(tier),
                                  x_force_overflow.data(tier), y_force_overflow.data(tier),
                                  z_force_overflow.data(tier), charges.data(tier),
                                  lennard_jones_indices.data(tier));
}

//-------------------------------------------------------------------------------------------------
CacheResourceKit<float> CacheResource::spData(const HybridTargetLevel tier) {
  return CacheResourceKit<float>(block_limit, atom_limit, x_coordinates.data(tier),
                                 y_coordinates.data(tier), z_coordinates.data(tier),
                                 x_velocities.data(tier), y_velocities.data(tier),
                                 z_velocities.data(tier), x_coordinate_overflow.data(tier),
                                 y_coordinate_overflow.data(tier),
                                 z_coordinate_overflow.data(tier), x_velocity_overflow.data(tier),
                                 y_velocity_overflow.data(tier), z_velocity_overflow.data(tier),
                                 x_force_overflow.data(tier), y_force_overflow.data(tier),
                                 z_force_overflow.data(tier), sp_charges.data(tier),
                                 lennard_jones_indices.data(tier));
}

//-------------------------------------------------------------------------------------------------
void CacheResource::repairPointers() {
  x_coordinates.swapTarget(&llint_data);
  y_coordinates.swapTarget(&llint_data);
  z_coordinates.swapTarget(&llint_data);
  x_velocities.swapTarget(&llint_data);
  y_velocities.swapTarget(&llint_data);
  z_velocities.swapTarget(&llint_data);
  x_coordinate_overflow.swapTarget(&int_data);
  y_coordinate_overflow.swapTarget(&int_data);
  z_coordinate_overflow.swapTarget(&int_data);
  x_velocity_overflow.swapTarget(&int_data);
  y_velocity_overflow.swapTarget(&int_data);
  z_velocity_overflow.swapTarget(&int_data);
  x_force_overflow.swapTarget(&int_data);
  y_force_overflow.swapTarget(&int_data);
  z_force_overflow.swapTarget(&int_data);
  lennard_jones_indices.swapTarget(&int_data);
}
  
} // namespace energy
} // namespace stormm
