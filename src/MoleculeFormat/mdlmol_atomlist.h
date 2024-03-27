// -*-c++-*-
#ifndef STORMM_MOLOBJ_ATOMLIST_H
#define STORMM_MOLOBJ_ATOMLIST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/textfile.h"

namespace stormm {
namespace structure {

using parse::TextFile;

/// \brief An atom list entry (this object can be assembled either from one of the deprecated
///        V2000 format lines after the bonds block, or from one of the "M  ALS" properties)
class MdlMolAtomList {
public:

  /// \brief The constructor can take all member variables (and all come with default values to
  ///        let this form of the constructor serve as the blank object constructor), or a pointer
  ///        to the line of a text file from which the information shall come.
  ///
  /// \param tf           Text of the original .sdf or .mol file, read into RAM
  /// \param line_number  Number of the line on which to read the data
  /// \param title        The title of the structure, if known, for error tracing purposes
  /// \{
  MdlMolAtomList(const std::vector<int> &atomic_numbers_in = {}, bool exclusions_in = false,
                 int atom_attachment_in = 0);

  MdlMolAtomList(const TextFile &tf, int line_number, const std::string &title = std::string(""));
  /// \}

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable to this object which has no const members or pointers to repair.
  /// \{
  MdlMolAtomList(const MdlMolAtomList &original) = default;
  MdlMolAtomList(MdlMolAtomList &&original) = default;
  MdlMolAtomList& operator=(const MdlMolAtomList &other) = default;
  MdlMolAtomList& operator=(MdlMolAtomList &&other) = default;
  /// \}

  /// \brief Get the number of element (atomic number) entries.
  int getEntryCount() const;

  /// \brief Get one of the elemental entries.
  ///
  /// \param index  The entry of interest
  int getEntry(int index) const;

  /// \brief Get a TRUE or FALSE reading on whether this atom list pertains to exclusions.
  bool applyToExclusions() const;

  /// \brief Get the coded letter for a V2000 atom list block entry detailing the exclusion
  ///        behavior.
  char getExclusionCode() const;

  /// \brief Get the central atom to which all listed elements will attach.
  int getAttachmentPoint() const;
  
private:
  int entry_count;                  ///< The number of atomic elements (identified by Z-numbers)
                                    ///<   in this list
  std::vector<int> atomic_numbers;  ///< Atomic (Z-) numbers of atoms that are to be excluded or
                                    ///<   included by processing this list
  bool exclusions;                  ///< Indicate whether this list covers atomic numbers which
                                    ///<   are to be excluded (TRUE) or included (FALSE)
  int atom_attachment;              ///< Attachment point of the list, an index of an atom in the
                                    ///<   molecule itself
};

} // namespace structure
} // namespace stormm

#endif
