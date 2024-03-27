#include "copyright.h"
#include "mm_enumerators.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ValenceWorkUnitSpecs input) {
  switch (input) {
  case ValenceWorkUnitSpecs::ATOM_INDEX_START:
    return std::string("ATOM_INDEX_START");
  case ValenceWorkUnitSpecs::ATOM_DIRECT_READ_START:
    return std::string("ATOM_DIRECT_READ_START");
  case ValenceWorkUnitSpecs::BOND_INSR_START:
    return std::string("BOND_INSR_START");
  case ValenceWorkUnitSpecs::ANGL_INSR_START:
    return std::string("ANGL_INSR_START");
  case ValenceWorkUnitSpecs::DIHE_INSR_START:
    return std::string("DIHE_INSR_START");
  case ValenceWorkUnitSpecs::ROGUE_14_INSR_START:
    return std::string("ROGUE_14_INSR_START");
  case ValenceWorkUnitSpecs::UBRD_INSR_START:
    return std::string("UBRD_INSR_START");
  case ValenceWorkUnitSpecs::CIMP_INSR_START:
    return std::string("CIMP_INSR_START");
  case ValenceWorkUnitSpecs::CMAP_INSR_START:
    return std::string("CMAP_INSR_START");
  case ValenceWorkUnitSpecs::RESTRAINT_INSR_START:
    return std::string("RESTRAINT_INSR_START");
  case ValenceWorkUnitSpecs::CONSTRAINT_INSR_START:
    return std::string("CONSTRAINT_INSR_START");
  case ValenceWorkUnitSpecs::SETTLE_INSR_START:
    return std::string("SETTLE_INSR_START");
  case ValenceWorkUnitSpecs::VSITE_INSR_START:
    return std::string("VSITE_INSR_START");
  case ValenceWorkUnitSpecs::MOVE_MASK_START:
    return std::string("MOVE_MASK_START");
  case ValenceWorkUnitSpecs::LJ_TABLE_START:
    return std::string("LJ_TABLE_START");
  case ValenceWorkUnitSpecs::LJ_14_TABLE_START:
    return std::string("LJ_14_TABLE_START");
  case ValenceWorkUnitSpecs::NONBONDED_TILE_START:
    return std::string("NONBONDED_TILE_START");
  case ValenceWorkUnitSpecs::ATOM_COUNT:
    return std::string("ATOM_COUNT");
  case ValenceWorkUnitSpecs::BOND_COUNT:
    return std::string("BOND_COUNT");
  case ValenceWorkUnitSpecs::ANGL_COUNT:
    return std::string("ANGL_COUNT");
  case ValenceWorkUnitSpecs::DIHE_COUNT:
    return std::string("DIHE_COUNT");
  case ValenceWorkUnitSpecs::ROGUE_14_COUNT:
    return std::string("ROGUE_14_COUNT");
  case ValenceWorkUnitSpecs::UBRD_COUNT:
    return std::string("UBRD_COUNT");
  case ValenceWorkUnitSpecs::CIMP_COUNT:
    return std::string("CIMP_COUNT");
  case ValenceWorkUnitSpecs::CMAP_COUNT:
    return std::string("CMAP_COUNT");
  case ValenceWorkUnitSpecs::RESTRAINT_COUNT:
    return std::string("RESTRAINT_COUNT");
  case ValenceWorkUnitSpecs::CONSTRAINT_GROUP_COUNT:
    return std::string("CONSTRAINT_GROUP_COUNT");
  case ValenceWorkUnitSpecs::SETTLE_GROUP_COUNT:
    return std::string("SETTLE_GROUP_COUNT");
  case ValenceWorkUnitSpecs::VSITE_COUNT:
    return std::string("VSITE_COUNT");
  case ValenceWorkUnitSpecs::LJ_TYPE_COUNT:
    return std::string("LJ_TYPE_COUNT");
  case ValenceWorkUnitSpecs::NONBONDED_TILE_COUNT:
    return std::string("NONBONDED_TILE_COUNT");
  case ValenceWorkUnitSpecs::DESCRIPTOR_COUNT:
    return std::string("DESCRIPTOR_COUNT");
  }
  __builtin_unreachable();
}

} // namespace mm
} // namespace stormm
