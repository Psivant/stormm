#include <vector>
#include "copyright.h"
#include "Math/summation.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace topology {

using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
  
//-------------------------------------------------------------------------------------------------
void atomValidityCheck(const int index, const int max_atoms, const char* class_caller,
                       const char* method_caller) {
  if (index >= max_atoms) {
    rtErr("Atom index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(max_atoms) + " atoms.", class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void atomValidityCheck(const int low_index, const int high_index, const int max_atoms,
                       const char* class_caller, const char* method_caller) {
  if (low_index < 0 || high_index > max_atoms || low_index >= high_index) {
    rtErr("Atom range " + std::to_string(low_index) + " - " + std::to_string(high_index) +
          " is invalid for a topology with " + std::to_string(max_atoms) + " atoms.",
          class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void residueValidityCheck(const int index, const int max_residues, const char* class_caller,
                          const char* method_caller) {
  if (index >= max_residues) {
    rtErr("Residue index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(max_residues) + " residues.", class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void residueValidityCheck(const int low_index, const int high_index, const int max_residues,
                          const char* class_caller, const char* method_caller) {
  if (low_index < 0 || high_index > max_residues || low_index >= high_index) {
    rtErr("Residue range " + std::to_string(low_index) + " - " + std::to_string(high_index) +
          " is invalid for a topology with " + std::to_string(max_residues) + " residues.",
          class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void moleculeValidityCheck(const int index, const int max_molecules, const char* class_caller,
                           const char* method_caller) {
  if (index >= max_molecules) {
    rtErr("Molecule index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(max_molecules) + " molecules.", class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void moleculeValidityCheck(const int low_index, const int high_index, const int max_molecules,
                           const char* class_caller, const char* method_caller) {
  if (low_index < 0 || high_index > max_molecules || low_index >= high_index) {
    rtErr("Molecule range " + std::to_string(low_index) + " - " + std::to_string(high_index) +
          " is invalid for a topology with " + std::to_string(max_molecules) + " molecules.",
          class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void affectorBoundsContribution(const int item_count, const int* atoms_from_the_item,
                                int *bounds_ptr) {
  if (atoms_from_the_item != nullptr) {
    for (int pos = 0; pos < item_count; pos++) {
      bounds_ptr[atoms_from_the_item[pos]] += 1;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void affectorListAssembly(const int item_count, const int* atoms_from_the_item,
                          int *bounds_ptr, int* list_ptr) {
  if (atoms_from_the_item != nullptr) {
    for (int pos = 0; pos < item_count; pos++) {
      const int atmi = atoms_from_the_item[pos];
      const int list_idx = bounds_ptr[atmi];
      list_ptr[list_idx] = pos;
      bounds_ptr[atmi] = list_idx + 1;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void markAffectorAtoms(std::vector<int> *affector_bounds, std::vector<int> *affector_list,
                       const int item_count, const int* i_atoms, const int* j_atoms,
                       const int* k_atoms, const int* l_atoms, const int* m_atoms) {

  // Zero the bounds array in case it is not initialized
  const int natom = affector_bounds->size() - 1;
  int* bounds_ptr = affector_bounds->data();
  for (int i = 0; i <= natom; i++) {
    bounds_ptr[i] = 0;
  }

  // Loop over all I, J, K, L, and M atoms, depending on what is available
  affectorBoundsContribution(item_count, i_atoms, bounds_ptr);
  affectorBoundsContribution(item_count, j_atoms, bounds_ptr);
  affectorBoundsContribution(item_count, k_atoms, bounds_ptr);
  affectorBoundsContribution(item_count, l_atoms, bounds_ptr);
  affectorBoundsContribution(item_count, m_atoms, bounds_ptr);

  // Compute the exclusive prefix sum to get the starting locations of the stretch of items
  // affecting each atom.
  prefixSumInPlace(affector_bounds, PrefixSumType::EXCLUSIVE, "markAffectorAtoms");
  int *list_ptr = affector_list->data();
  affectorListAssembly(item_count, i_atoms, bounds_ptr, list_ptr);
  affectorListAssembly(item_count, j_atoms, bounds_ptr, list_ptr);
  affectorListAssembly(item_count, k_atoms, bounds_ptr, list_ptr);
  affectorListAssembly(item_count, l_atoms, bounds_ptr, list_ptr);
  affectorListAssembly(item_count, m_atoms, bounds_ptr, list_ptr);

  // Roll back the prefix sum to recover the original, exclusive sum.
  for (int i = natom; i > 0; i--) {
    bounds_ptr[i] = bounds_ptr[i - 1];
  }
  bounds_ptr[0] = 0;
}

} // namespace topology
} // namespace stormm
