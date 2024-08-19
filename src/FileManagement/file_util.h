// -*-c++-*-
#ifndef STORMM_FILE_UTIL_H
#define STORMM_FILE_UTIL_H

#include <fstream>
#include <iostream>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "Trajectory/trajectory_enumerators.h"
#include "file_enumerators.h"

namespace stormm {
namespace diskutil {

using constants::ExceptionResponse;
using parse::TextFile;
using trajectory::CoordinateFileKind;
  
const char default_amber_crd_extension[] = "crd";
const char default_amber_inpcrd_extension[] = "inpcrd";
const char default_amber_ascii_rst_extension[] = "rst";
const char default_amber_netcdf_extension[] = "cdf";
const char default_amber_netcdf_rst_extension[] = "rcdf";
const char default_sd_file_extension[] = "sdf";

/// \brief Open a file for output writing.  This encapsulates error messages in the event that
///        the file cannot be opened as expected.
///
/// \param filename     Name of the file to write
/// \param expectation  Conditions to look for in order to successfully open the file for writing
/// \param style        The kind of file to write, ascii or binary
/// \param description  Description of the file, to help identify problems if it cannot be written
std::ofstream openOutputFile(const std::string &filename,
                             PrintSituation expectation = PrintSituation::OPEN_NEW,
                             const std::string &description = std::string(""),
                             DataFormat style = DataFormat::ASCII);

/// \brief Remove a file.  This encapsulates error messages in the event that the file cannot be
///        removed as requested.
///
/// \param filename  Name of the file to remove
/// \param policy    Indicates what to do if the file cannot be removed for some reason
int removeFile(const std::string &filename, ExceptionResponse policy = ExceptionResponse::WARN);

/// \brief Get the nature of a trajectory file based on the stated format (this will return
///        binary or ASCII based on the stated trajectory file kind)
///
/// \param cfkind  The trajectory file kind
DataFormat getTrajectoryFormat(CoordinateFileKind cfkind);

/// \brief Detect various coordinate file types.
///
/// Overloaded:
///   - Work from the file name
///   - Work from a pre-translated TextFile object
///
/// \param file_name  Name of the file to test
/// \param caller     Name of the calling function (optional)
/// \param tf         Text of an asci--format file already translated into RAM
/// \{
CoordinateFileKind detectCoordinateFileKind(const std::string &file_name,
                                            const std::string &caller = std::string(""));
CoordinateFileKind detectCoordinateFileKind(const TextFile &tf);
/// \}

/// \brief Get the default file extension for a file based on its type.
///
/// \param kind  Format of the file--different enumerators indicate different purposes for the file
std::string getDefaultFileExtension(CoordinateFileKind kind);

/// \brief Infer the intended type of a coordinate file based on its extension.  This is not a
///        substitute for detectCoordinateFileKind() above, which will look at the actual contents
///        of the file to determine the type.  Rather, it is useful when the name of a file to be
///        written has been provided and the type remains undecided.  This function is insensitive
///        to the letter case and will return UNKNOWN if the file type cannot be identified.
///
/// \param file_name  Name of the file of interest
CoordinateFileKind inferCoordinateFileKind(const std::string &file_name);

/// \brief Assemble the name of a file path given three typical components.
///
/// \param root_path  The root of the path, including any parent directories.  This can be supplied
///                   as an absolute path or, if the file is to be found in some tree within the
///                   current working directory, a relative path.  If no such input is supplied,
///                   the root path and any subsequent path separator will be omitted from the
///                   result.
/// \param base_path  The base file name, as would be found by using an ls command in the file's
///                   current directory
/// \param extn       The file extension, to be added to the rest of the file name following a
///                   '.' character.  If no such input is supplied, no '.' will be appended,
///                   either.
std::string assembleFilePath(const std::string &root_path, const std::string &base_path,
                             const std::string &extn);
  
} // namespace diskutil
} // namespace stormm

#endif
