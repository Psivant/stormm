// -*-c++-*-
#ifndef STORMM_COORDINATE_UTIL_H
#define STORMM_COORDINATE_UTIL_H

#include <string>
#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "copyright.h"
#include "coordinateframe.h"
#include "coordinate_series.h"
#include "phasespace.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace trajectory {

using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
  
/// \brief Type indices for all coordinate objects
/// \{
static const size_t cf_type_index = std::type_index(typeid(CoordinateFrame)).hash_code();
static const size_t ps_type_index = std::type_index(typeid(PhaseSpace)).hash_code();
static const size_t cs_dtype_index = std::type_index(typeid(CoordinateSeries<double>)).hash_code();
static const size_t cs_ftype_index = std::type_index(typeid(CoordinateSeries<float>)).hash_code();
static const size_t cs_stype_index = std::type_index(typeid(CoordinateSeries<short>)).hash_code();
static const size_t cs_itype_index = std::type_index(typeid(CoordinateSeries<int>)).hash_code();
static const size_t cs_ltype_index = std::type_index(typeid(CoordinateSeries<llint>)).hash_code();
static const size_t poly_ps_type_index = std::type_index(typeid(PhaseSpaceSynthesis)).hash_code();
static const size_t cdns_type_index = std::type_index(typeid(Condensate)).hash_code();
/// \}
  
/// \brief Translate the detected type index of a coordinate object into a human-readable string
///        for error reporting or other explanatory purposes.
std::string nameCoordinateType(size_t ct_coords);

/// \brief Check the validity of a given copy operation, based on the memory level from which
///        coordinates will originate and the memory level they will go to.  This provides more
///        explanatory error messages for developers making heavy use of coordCopy() between the
///        CPU and GPU.
///
/// \param destination       The object into which the coordinates will go
/// \param origin            The object from which coordinates will be taken
/// \param destination_tier  The memory level (CPU host or GPU device) that will receive
///                          coordinates
/// \param origin_tier       The memory level from which coordinates originate
template <typename Tdest, typename Torig>
void checkCopyValidity(const Tdest *destination, const Torig &origin,
                       HybridTargetLevel destination_tier, HybridTargetLevel origin_tier);
  
} // namespace trajectory
} // namespace stormm

// This library provides automatic inclusion of the coordinate objects' type indices.
namespace stormm {

using trajectory::cf_type_index;
using trajectory::ps_type_index;
using trajectory::cs_dtype_index;
using trajectory::cs_ftype_index;
using trajectory::cs_stype_index;
using trajectory::cs_itype_index;
using trajectory::cs_ltype_index;
using trajectory::poly_ps_type_index;
using trajectory::cdns_type_index;

} // namespace stormm

#include "coordinate_util.tpp"

#endif
