// -*-c++-*-
#ifndef STORMM_FILE_SHAPSHOT_H
#define STORMM_FILE_SHAPSHOT_H

#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Parsing/parse.h"

namespace stormm {
namespace testing {

using diskutil::PrintSituation;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::TextFile;

/// \brief Read reference data from a snapshot file
///
/// Overloaded:
///   - Read from a TextFile object
///   - Read from a file, given a file name (cascades into the original function)
///
/// \param tf        TextFile object transcribing an ascii text file with reference data
/// \param filename  Name of the file containing reference data
/// \param label     Expected label of the data to read (a check against file corruption or
///                  misdirection)
/// \{
std::vector<PolyNumeric> readSnapshot(const TextFile &tf,
                                      const std::string &label = std::string(""));

std::vector<PolyNumeric> readSnapshot(const std::string &filename,
                                      const std::string &label = std::string(""));
/// \}

/// \brief Write a vector of scalar numbers to a snapshot file.  This can be useful for taking
///        checkpoints of more intensive calculations, when subtle corruption is hard to root out
///        or drift over time is expected.
///
/// \param filename     Name of the file to write
/// \param content      Numerical (or char4) data to transcribe
/// \param label        Label under which to store the data (this will become the name of a
///                     variable in a Matlab- or Octave-readable script)
/// \param tol          Tolerance for future testing.  Setting this tolerance tighter than the
///                     actual tests require is not a problem, for the purposes of correcting
///                     data that is within bounds but shows some slight drift.
/// \param data_format  Format of the data to write.  The CHAR4 format is currently unavailable
///                     due to the fact that vectors of char4 tuples cannot be written to disk in
///                     this manner.
void writeSnapshot(const std::string &filename, const std::vector<PolyNumeric> &content,
                   const std::string &label, const double tol, const NumberFormat data_format,
                   const PrintSituation expectation);

/// \brief Read text segments from a snapshot file
///
/// Overloaded:
///   - Read from a TextFile object (convert a subset of the information to another TextFile)
///   - Read from a file, given a file name (cascades into the previous function)
///
/// \param tf        TextFile object transcribing an ascii text file with reference data
/// \param filename  Name of the file containing reference data
/// \param label     Expected label of the data to read (a check against file corruption or
///                  misdirection)
/// \{
TextFile readTextSnapshot(const TextFile &tf, const std::string &label = std::string(""));

TextFile readTextSnapshot(const std::string &filename, const std::string &label = std::string(""));
/// \}

/// \brief Write a segment of text to a snapshot file.  This can be useful for taking checkpoints
///        of complex sequences and updating the targets as the underlying code for writing the
///        output evolves in supervised increments.
///
/// Overloaded:
///   - Transcribe a formatted string (must contain carriage returns where line breaks are to
///     appear)
///   - Transcribe a TextFile object
///
/// \param filename    Name of the file to write
/// \param content     The text sequence to transcribe
/// \param label       Label under which to store the data (this will have special demarcations,
///                    distinct from those for numeric snapshots)
/// \para expectation  Indicate in what state the snapshot file is expected to be found
/// \{
void writeTextSnapshot(const std::string &filename, const TextFile &content,
                       const std::string &label, PrintSituation expectation);

void writeTextSnapshot(const std::string &filename, const std::string &content,
                       const std::string &label, PrintSituation expectation);
/// \}

} // namespace testing
} // namespace stormm

#endif
