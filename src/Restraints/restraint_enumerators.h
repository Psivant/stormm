// -*-c++-*-
#include <string>
#include "copyright.h"

namespace stormm {
namespace restraints {

/// \brief An enumerator for the various restraint assortments that one can apply to a system.
enum class RestraintEnsemble {
  SPECIFIC_ATOMS,            ///< The restraint applies only to specific atoms in the system and
                             ///<   guides their arrangement towards an explicit target
  PREVENT_HBONDS,            ///< Penalize hydrogen bonds that may form between identified donor
                             ///<   and acceptor pairs.
  PRESERVE_HEAVY_DIHEDRALS,  ///< Apply dihedral restraints to dihedrals involving purely heavy
                             ///<   atoms.
  PRESERVE_POSITIONS,        ///< Preserve the positions of masked atoms in the molecule.
  PRESERVE_DISTANCES         ///< Select heavy atoms throughout the molecule and apply distance
                             ///<   restraints to maintain the relative displacements.  This is
                             ///<   a distinct method of preserving molecular geometry, but works
                             ///<   towards the same goal as PRESERVE_HEAVY_DIHEDRALS
};

/// \brief Produce a human-readable string based on an enumerated value.  Various overloads of this
///        function are found in other libraries and namespaces.
std::string getEnumerationName(RestraintEnsemble input);
  
/// \brief Translate a user input value or other string code into the appropriate RestraintEnsemble
///        enumeration.
///
/// \param rst_group  Human-parseable name for the restraint ensemble enumeration
RestraintEnsemble translateRestraintEnsemble(const std::string &rst_group);

} // namespace restraints
} // namespace stormm
