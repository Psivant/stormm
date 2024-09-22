// -*-c++-*-
#ifndef STORMM_TOPOLOGY_LIMITS_H
#define STORMM_TOPOLOGY_LIMITS_H

#include "copyright.h"

namespace stormm {
namespace topology {

/// \brief Check an atom index against the available content, reporting an error if the bounds are
///        exceeded.
///
/// \param index          Index of the atom of interest
/// \param max_atoms      The maximum number of available atoms
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void atomValidityCheck(int index, int max_atoms, const char* class_caller = nullptr,
                       const char* method_caller = nullptr);

/// \brief Check an atom range against the available content, reporting an error if the bounds are
///        exceeded.
///
/// \param low_index      Atom index at the left end of the range of interest
/// \param high_index     The right end of the range of interest (this can be the maximum number of
///                       atoms in the system--the group is assumed to be [low_index, high_index)
/// \param max_atoms      The maximum number of available atoms
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void atomValidityCheck(int low_index, int high_index, int max_atoms,
                       const char* class_caller = nullptr, const char* method_caller = nullptr);

/// \brief Check a residue index against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param index          Index of the residue of interest
/// \param max_atoms      The maximum number of available residues
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void residueValidityCheck(int index, int max_residues, const char* class_caller = nullptr,
                          const char* method_caller = nullptr);

/// \brief Check an residue range against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param low_index      Residue index at the left end of the range of interest
/// \param high_index     The right end of the range of interest (this can be the systems' maximum
///                       number of residues--the group is assumed to be [low_index, high_index)
/// \param max_residues   The maximum number of available residues
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void residueValidityCheck(int low_index, int high_index, int max_residues,
                          const char* class_caller = nullptr, const char* method_caller = nullptr);


/// \brief Check a molecule index against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param index          Index of the molecule of interest
/// \param max_atoms      The maximum number of available molecules
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void moleculeValidityCheck(int index, int max_molecules, const char* class_caller = nullptr,
                           const char* method_caller = nullptr);

/// \brief Check a molecule range against the available content, reporting an error if the bounds
///        are exceeded.
///
/// \param low_index      Molecule index at the left end of the range of interest
/// \param high_index     The right end of the range of interest (this can be the system's maximum
///                       number of molecules--the group is assumed to be [low_index, high_index)
/// \param max_molecules  The maximum number of available molecules
/// \param class_caller   Name of the calling function or class, to be included in error messages
/// \param method_caller  Name of a class member function calling this routine, if applicable
void moleculeValidityCheck(int low_index, int high_index, int max_molecules,
                           const char* class_caller = nullptr,
                           const char* method_caller = nullptr);

/// \brief Accumulate the numbers of some topological item (bond, angle, CMAP) that affect any
///        given atom of a topology.  This can be applied repeatedly for items affecting more than
///        one atom.
///
/// \param item_count           The number of items (trusted length of atoms_from_the_item)
/// \param atoms_from_the_item  Array of atoms within whatever item of interest
/// \param bounds_ptr           The developing bounds array
void affectorBoundsContribution(int item_count, const int* atoms_from_the_item, int *bounds_ptr);

/// \brief Assemble the list of numbered items affecting any given atom, provided a bounds list
///        and a list of atoms to which the items contribute.  This can be applied repeatedly for
///        items affecting more than one atom.
///
/// \param item_count           The number of items (trusted length of atoms_from_the_item)
/// \param atoms_from_the_item  Array of atoms within whatever item of interest
/// \param bounds_ptr           The developing bounds array
/// \param list_ptr             The developing list of items affecting each atom
void affectorListAssembly(int item_count, const int* atoms_from_the_item, int *bounds_ptr,
                          int* list_ptr);

/// \brief Mark all atoms affected by a particular set of terms, each term affecting up to five
///        atoms.
///
/// \param affector_bounds  Pre-allocated array of bounds for the list of all numbered items found
///                         to affect each atom.  Length must be the total number of atoms in the
///                         system, plus one.  Filled and returned.
/// \param affector_list    The list of all terms affecting each atom.  Its bounds array is
///                         affector_bounds.  Filled and returned.
/// \param item_count       The number of items (trusted length of i_atoms, j_atoms, ...)
/// \param i_atoms          The first atoms of each item, i.e. the I atom of each bond
/// \param j_atoms          Second atoms of each item
/// \param k_atoms          Third atoms of each item
/// \param l_atoms          Fourth atoms of each item
/// \param m_atoms          Fifth atoms of each item, i.e. atom M from each CMAP term
void markAffectorAtoms(std::vector<int> *affector_bounds, std::vector<int> *affector_list,
                       int item_count, const int* i_atoms, const int* j_atoms = nullptr,
                       const int* k_atoms = nullptr, const int* l_atoms = nullptr,
                       const int* m_atoms = nullptr);

} // namespace topology
} // namespace stormm

#endif
