// -*-c++-*-
#ifndef STORMM_NAMELIST_ENUMERATORS_H
#define STORMM_NAMELIST_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace namelist {

/// \brief Enumerator to describe all types of namelist variables
enum class NamelistType {
  INTEGER, REAL, STRING, STRUCT
};

/// \brief Enumerate choices on whether to accept multiple values (this exists to make the code
///        more legible and is otherwise translated into a boolean variable)
enum class InputRepeats {
  NO, YES
};

/// \brief Spell out the choice being made when supplying optional input to the various
///        setDefault(...) methods.  A "YES" value will cause the entry_index counter to
///        increment with each contributed default setting.
enum class DefaultIsObligatory {
  NO, YES
};

/// \brief Enumerate the possible ways that input could have come in (or failed to be obtained)
enum class InputStatus {
  USER_SPECIFIED, DEFAULT, MISSING
};

/// \brief Enumeration to indicate whether a restraint encodes its atoms by atom masks (strings) or
///        by topological indices.
enum class RestraintAnchoring {
  ATOMMASK,  ///< Atom masks encode the connected atoms
  INDICES,   ///< Topological indices specify the restrained atoms
  MIXED,     ///< Both masks and topological indices specify the mask
  UNKNOWN    ///< Who knows?  This should never be possible to achieve.
};

/// \brief Enumeration to indicate the requirement of STRUCT subkeys
enum class KeyRequirement {
  OPTIONAL,  ///< The keyword is not required, but can be specified
  REQUIRED,  ///< The keyword is essential and must be supplied in order to create a valid STRUCT
  BOGUS      ///< The option should actually NOT be present--useful for situations in which a
             ///<   pre-compiled namelist is re-used but not intended to take arguments from
             ///<   certain STRUCT options
};

/// \brief List the ways that a namelist can be presented in a transcript file.
enum class NamelistIntroduction {
  HEADER,          ///< Print the namelist title and category headings
  COMPACT_HEADER,  ///< Print a more compact header that will not stand out as much against section
                   ///<   headings
  BLANK_LINE,      ///< Print a single blank line to space one namelist from another (no headings
                   ///<   or title)
  NONE             ///< Prin no headings, title, or spacing
};
  
/// \brief Return a string corresponding to the namelist data type enumerations.  Various overloads
///        of this function are also found in other libraries.
///
/// \param input  The type of namelist keyword in question
/// \{
std::string getEnumerationName(NamelistType input);
std::string getEnumerationName(InputRepeats input);
std::string getEnumerationName(DefaultIsObligatory input);
std::string getEnumerationName(InputStatus input);
std::string getEnumerationName(RestraintAnchoring input);
std::string getEnumerationName(KeyRequirement input);
std::string getEnumerationName(NamelistIntroduction input);
/// \}

/// \brief Interpret various inputs into a KeyRequirement enumeration.
///
/// \param input  The word or phrase to translate
KeyRequirement translateKeyRequirement(const std::string &input);
  
/// \brief Return a string corresponding to an input status setting.
///
/// \param stt  The status of interest
std::string getInputStatusString(InputStatus stt);
  
} // namespace namelist
} // namespace stormm

#endif
