// -*-c++-*-
#ifndef STORMM_MOLECULE_FILE_IO
#define STORMM_MOLECULE_FILE_IO

#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/textfile.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using constants::ExceptionResponse;
using parse::TextFile;
using parse::TextFileReader;
  
/// \brief Find the end of the MDL MOL entry formatting within a specified range.  If the
///        end of the search range takes its default value of -1, the search will continue until
///        the end of the TextFile object.
///
/// Overloaded:
///   - Accept a TextFile object
///   - Accept the abstract of a TextFile object
///
/// \param tfr          Abstract of the text
/// \param tf           Text of the MDL MOL file (or SD file), read into RAM
/// \param line_start   The point at which to assume the MDL MOL format entry begins
/// \param line_end_in  The limit at which to stop searching for the MDL MOL format end.  Defaults
///                     to -1 to indicate "end of file."
/// \{
int getMdlMolSectionEnd(const TextFileReader &tfr, int line_start, int line_end_in = -1);

int getMdlMolSectionEnd(const TextFile &tf, int line_start, int line_end_in = -1);
/// \}

/// \brief Find the end of the SD file compound entry.  If the end of search range takes its
///        default value of -1, the search will continue until the end of the TextFile object.
///        The overloading and formal argument descriptors of this function follow from
///        getMdlMolSectionEnd().
/// \{
int getCompoundSectionEnd(const TextFileReader &tfr, int line_start, int line_end_in = -1);

int getCompoundSectionEnd(const TextFile &tf, int line_start, int line_end_in = -1);
/// \}
  
/// \brief Find the version stamp in an MDL MOL file (or a part of an SD file, .sdf).  Assume the
///        legacy V2000 if no version is found.
///
/// Overloaded:
///   - Accept a pointer to a single line
///   - Accept the original TextFile object (safest method, as it implies a bounds check on the
///     line of interest)
///
/// \param text   Text to search, probably from the fourth line of the MOL file or SD entry
/// \param nchar  The number of characters on the line
/// \{
MdlMolVersion findMolObjVersion(const char* text, const int nchar);
  
MdlMolVersion findMolObjVersion(const TextFile &tf, const int line_number);
/// \}

/// \brief Find the limits for different MDL MOL entries within an SD file.
///
/// \param tf  The SD file, read into RAM
std::vector<int2> findSdfMolEntryLimits(const TextFile &tf);

/// \brief Extract the vector of molecular coordiantes from an SDF file.
///
/// Overloaded:
///   - Accept the line numbers defining the limits of the frame
///   - Accept the frame number to extract
///
/// \param tf            The text file to parse, read from disk into RAM
/// \param line_start    First line at which to begin parsing an MDL MOL entry
/// \param line_end_in   Last line at which to seek data from an MDL MOL entry
/// \param frame_number  The frame number to parse
/// \param policy        Action to take if errors are encountered
/// \{
std::vector<double3> extractSdfCoordinates(const TextFile &tf, int line_start, int line_end_in,
                                           ExceptionResponse policy = ExceptionResponse::WARN);

std::vector<double3> extractSdfCoordinates(const TextFile &tf, int frame_number = 0,
                                           ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

} // namespace structure
} // namespace stormm

#endif
