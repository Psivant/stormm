#ifndef STORMM_TRIPOS_FORMAT_H
#define STORMM_TRIPOS_FORMAT_H

#include <string>
#include "copyright.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace structure {
    
/// \brief The Tripos mol2 format, as annotated in the mol2 .pdf document found within this
///        source directory.
struct TriposMol {

  /// \brief Constructors include a blank constructor, constructors from file names as strings or
  ///        character arrays, and a constructor based on existing topologies.
  /// \{
  TriposMol();
  TriposMol(const std::string &filename);
  TriposMol(const char* filename);
  TriposMol(const AtomGraph &mol_in);
  /// \}
  
  // Getter functions for each member variable
  void getFileName();
  void getTitle();
  int getAtomCount();
  int getBondCount();
  int getSubstructureCount();
  int getFeatureCount();
  int getSetCount();
  TriposMoleculeKind getMoleculeKind();
  TriposChargeKind getChargeKind();

private:
  int atom_count;                ///< Number of atoms in the structure
  int bond_count;                ///< Number of bonds in the structure
  int substructure_count;        ///< Number of substructures defined for this system
  int feature_count;             ///< Number of features defined for this system
  int set_count;                 ///< Number of sets in the molecule or system
  TriposMoleculeKind mol_kind;   ///< Classification of molecule (i.e. protein, DNA)
  TriposChargeKind charge_kind;  ///< Charge style used to develop the MM properties
};

} // namespace structure
} // namespace stormm

#endif
