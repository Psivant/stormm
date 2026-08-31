#include <cstdio>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "file_util.h"
#include "file_listing.h"

namespace stormm {
namespace diskutil {

using constants::CaseSensitivity;
using parse::digestText;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::separateText;
using parse::strcmpCased;
using parse::strncmpCased;
using parse::TextFileReader;
using parse::TextOrigin;
using parse::verifyNumberFormat;

//-------------------------------------------------------------------------------------------------
std::ofstream openOutputFile(const std::string &filename, const PrintSituation expectation,
                             const std::string &description, const DataFormat style) {

  // Check the file conditions: look before you leap
  std::ofstream foutp;
  bool directory_in_the_way = false;
  switch (expectation) {
  case PrintSituation::UNKNOWN:
  case PrintSituation::OPEN_NEW:
    switch (getDrivePathType(filename)) {
    case DrivePathType::FILE:
      rtErr("Unable to open a new file \'" + filename + "\' because it already exists.  If this "
            "file is to be overwritten, the command line option \"-O\" appended or placed amongst "
            "the program options can often override this safety mechanism.  Activity: " +
            description + ".", "openOutputFile");
      break;
    case DrivePathType::DIRECTORY:
      directory_in_the_way = true;
      break;
    case DrivePathType::REGEXP:
      switch (style) {
      case DataFormat::ASCII:
        foutp.open(filename, std::ofstream::out);
        break;
      case DataFormat::BINARY:

	// The std::ofstream specifiers occupy different bits of a single
	// integer-like data type.  Use the | operator to combine them.
        foutp.open(filename, std::ofstream::out | std::ofstream::binary);
        break;
      }
      break;
    }
    break;
  case PrintSituation::APPEND:
    switch (style) {
    case DataFormat::ASCII:
      foutp.open(filename, std::ofstream::ate | std::ofstream::app);
      break;
    case DataFormat::BINARY:
      foutp.open(filename, std::ofstream::ate | std::ofstream::app | std::ofstream::binary);
      break;
    }
    break;
  case PrintSituation::OVERWRITE:
    switch (getDrivePathType(filename)) {
    case DrivePathType::FILE:
    case DrivePathType::REGEXP:
      switch (style) {
      case DataFormat::ASCII:
        foutp.open(filename, std::ofstream::trunc);
        break;
      case DataFormat::BINARY:
        foutp.open(filename, std::ofstream::trunc | std::ofstream::binary);
        break;
      }
      break;
    case DrivePathType::DIRECTORY:
      directory_in_the_way = true;
      break;
    }
    break;
  }
  if (directory_in_the_way) {
    rtErr("Unable to open a new file \'" + filename + "\' because a directory of the same name "
          "exists.  No remedy is provided in the command line -O option for this case.  "
          "Activity: " + description + ".", "openOutputFile");
  }

  // Check that the file has been opened
  if (foutp.is_open() == false) {

    // Check for the existence of the directory in which the file is to be written
    
    rtErr("Attempt to open file \'" + filename + "\' failed.  Bad writing permissions or "
          "insufficient disk space may be the problem.  Activity: " + description + ".",
          "openOutputFile");
  }
  return foutp;
}

#if STORMM_INCLUDE_NETCDF
//-------------------------------------------------------------------------------------------------
int openOutputNetCDF(const std::string &filename, const PrintSituation expectation,
                     const char* caller) {
  int result, chk_code;
  switch (expectation) {
  case PrintSituation::OPEN_NEW:
    chk_code = nc_create(filename.c_str(), NC_NOCLOBBER, &result);
    break;
  case PrintSituation::APPEND:
    rtErr("In general, existing NetCDF files cannot be appended.  The developer must write a "
          "special-purpose routine to read the existing file into memory, then write its contents "
          "back to disk and return the active NetCDF file identifier to continue writing.",
          "openNetCDF");
  case PrintSituation::OVERWRITE:
    chk_code = nc_create(filename.c_str(), NC_CLOBBER, &result);
    break;
  case PrintSituation::UNKNOWN:
    rtErr("The developer must specify whether the file is expected to not exist (OPEN_NEW) or is "
          "irrelevant (OVERWRITE).", "openNetCDF");
  }
  switch (chk_code) {
  case NC_NOERR:
    break;
  case NC_EEXIST:
    rtErr("The NetCDF file " + filename + " already exists, but overwriting is forbidden.",
          caller);
  case NC_ENFILE:
    rtErr("The number of open NetCDF files has exceeded the allowed maximum.", caller);
  default:
    rtErr("Net CDF error: " + std::string(nc_strerror(chk_code)), caller);
  }
  return result;
}
#endif
  
//-------------------------------------------------------------------------------------------------
int removeFile(const std::string &filename, const ExceptionResponse policy) {
  switch (getDrivePathType(filename)) {
  case DrivePathType::FILE:
    break;
  case DrivePathType::DIRECTORY:
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(filename + " is a directory.", "removeFile");
    case ExceptionResponse::WARN:
      rtWarn(filename + " is a directory.  Use stormmRmdir() to remove it.", "removeFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 1;
  case DrivePathType::REGEXP:
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(filename + " does not exist and therefore cannot be removed.", "removeFile");
    case ExceptionResponse::WARN:
      rtWarn(filename + " does not exist and therefore cannot be removed.", "removeFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 1;
  }
  const int result = remove(filename.c_str());
  if (result != 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(filename + " could not be removed as requested.", "removeFile");
    case ExceptionResponse::WARN:
      rtWarn(filename + " could not be removed as requested.", "removeFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
DataFormat getTrajectoryFormat(CoordinateFileKind cfkind) {
  switch (cfkind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::PDB:
    return DataFormat::ASCII;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return DataFormat::BINARY;
  case CoordinateFileKind::UNKNOWN:
    rtWarn("Unable to determine the nature of an UNKNOWN kind trajectory file.  Reporting BINARY.",
           "getTrajectoryFormat");
    return DataFormat::BINARY;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
TextFile readFileSample(const std::string &file_name, const std::string &caller,
                        const size_t nchar) {

  // Read the first 2kB of the file, if that much exists.  Look for chracters that might indicate
  // it is a binary file.  If no such characters are found, try to parse it as one of the known
  // ASCII formats.
  std::ifstream finp;
  finp.open(file_name.c_str());
  if (finp.is_open() == false) {
    if (caller.size() == 0) {
      rtErr(file_name + " was not found.", "readFileSample");
    }
    else {
      rtErr(file_name + " was not found when called from " + caller + ".", "readFileSample");
    }
  }  
  const int maxchar = 2048;
  std::vector<char> buffer(maxchar);
  int pos = 0;
  char c;
  bool is_binary = false;
  while (pos < maxchar && finp.get(c) && is_binary == false) {
    is_binary = (is_binary || static_cast<int>(c) < 0);
    buffer[pos] = c;
    pos++;
  }
  finp.close();
  const size_t nchar_found = (is_binary) ? 0 : pos;
  std::string first_part(nchar_found, ' ');
  for (size_t i = 0; i < nchar_found; i++) {
    first_part[i] = buffer[i];
  }
  const TextFile result(file_name, TextOrigin::RAM, first_part, "readFileSample");
  return result;
}
  
//-------------------------------------------------------------------------------------------------
CoordinateFileKind detectCoordinateFileKind(const TextFile &tf) {
  const TextFileReader tfr = tf.data();

  // Test for an Amber coordinate file
  if (tfr.line_count < 3) {
    return CoordinateFileKind::UNKNOWN;
  }
  const std::vector<std::string> line_two = separateText(&tfr.text[tfr.line_limits[1]],
                                                         tfr.line_limits[2] - tfr.line_limits[1]);
  if (line_two.size() == 2 && verifyNumberFormat(line_two[0].c_str(), NumberFormat::INTEGER) &&
      stoi(line_two[0]) > 0 &&
      (verifyNumberFormat(line_two[1].c_str(), NumberFormat::SCIENTIFIC) ||
       verifyNumberFormat(line_two[1].c_str(), NumberFormat::STANDARD_REAL))) {

    // Check the third line to confirm
    if (tfr.line_limits[3] - tfr.line_limits[2] >= 36 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2], 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 12, 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 24, 12)) {

      // Some input files will still have a time stamp.  A check for the presence of actual
      // velocities will be needed later, as this routine may have only been fed a stub from the
      // head of the file.
      return CoordinateFileKind::AMBER_ASCII_RST;
    }
  }
  else if (line_two.size() == 1 &&
           verifyNumberFormat(line_two[0].c_str(), NumberFormat::INTEGER) &&
           stoi(line_two[0]) > 0) {

    // Check the third line to confirm
    if (tfr.line_limits[3] - tfr.line_limits[2] >= 36 &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2], 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 12, 12) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[2] + 24, 12)) {
      return CoordinateFileKind::AMBER_INPCRD;
    }
  }
  else if (tfr.line_limits[2] - tfr.line_limits[1] >= 24 &&
           verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1], 8) &&
           verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1] + 8, 8) &&
           verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[1] + 16, 8)) {
    return CoordinateFileKind::AMBER_CRD;
  }

  // If the fourth line looks like an MDL MOL-format counts line, return that result after checking
  // the first atom line.
  if (tfr.line_count > 4 &&
      verifyContents(tfr, 3, 0, 3, NumberFormat::INTEGER) &&
      verifyContents(tfr, 3, 3, 3, NumberFormat::INTEGER) &&
      readIntegerValue(&tfr.text[tfr.line_limits[3]], 0, 3) >= 0) {
    if (verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[4], 10) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[4] + 10, 10) &&
        verifyNumberFormat(tfr.text, NumberFormat::STANDARD_REAL, tfr.line_limits[4] + 20, 10)) {
      return CoordinateFileKind::SDF;
    }
  }

  // If the phrase "V2000" or "V3000" can be detected anywhere in the fourth through tenth lines,
  // return that result as well.
  const int nsearch_lines = std::min(10, tfr.line_count);
  for (int i = 3; i < nsearch_lines; i++) {
    const int jlim = tfr.line_limits[i + 1] - 4;
    for (int j = tfr.line_limits[i]; j < jlim; j++) {
      if (tfr.text[j] == 'V' && (tfr.text[j + 1] == '2' || tfr.text[j + 1] == '3') &&
          tfr.text[j + 2] == '0' && tfr.text[j + 3] == '0' && tfr.text[j + 4] == '0') {
        return CoordinateFileKind::SDF;
      }
    }    
  }

  // If the file type has not yet been ascertained, return UNKNOWN
  return CoordinateFileKind::UNKNOWN;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind detectCoordinateFileKind(const std::string &file_name,
                                            const std::string &caller) {
  const TextFile tf = readFileSample(file_name, caller, 2048);
  if (tf.getTextSize() == 0) {

    // TBD: Determine the binary format.
    return CoordinateFileKind::UNKNOWN;
  }
  return detectCoordinateFileKind(tf);
}

//-------------------------------------------------------------------------------------------------
TopologyKind detectTopologyKind(const TextFile &tf) {
  const TextFileReader tfr = tf.data();
  
  // Detect the AMBER, CHARMM, or GROMACS formats
  bool seek_ambr_version = true;
  bool seek_ambr_pointers = true;
  bool seek_chrm_version = true;
  bool seek_chrm_mass = true;
  bool seek_gmx_dir = true;
  bool seek_gmx_incl = true;
  int line_idx = 0;
  while (line_idx < tfr.line_count && (seek_ambr_version || seek_ambr_pointers) &&
         (seek_chrm_version || seek_chrm_mass) && (seek_gmx_dir || seek_gmx_incl)) {
    const std::string ln_digest = digestText(tfr.text, tfr.line_limits[line_idx],
                                             tfr.line_lengths[line_idx], "[]", 2, "", 0);
    const std::vector<std::string> ln_contents = separateText(ln_digest);
    const size_t word_count = ln_contents.size();

    // Check for AMBER hallmarks
    if (seek_ambr_version && strcmpCased(ln_contents[0], "%version", CaseSensitivity::NO)) {
      seek_ambr_version = false;
    }
    if (seek_ambr_pointers && word_count >= 2 &&
        strcmpCased(ln_contents[0], "%flag", CaseSensitivity::NO) &&
        strcmpCased(ln_contents[1], "pointers", CaseSensitivity::NO)) {
      seek_ambr_pointers = false;
    }

    // Check for CHARMM hallmarks
    if (seek_chrm_version) {
      bool all_integer = true;
      for (size_t i = 0; i < word_count; i++) {
        all_integer = (all_integer || verifyContents(ln_contents[i], NumberFormat::INTEGER));
      }
      if (all_integer) {
        seek_chrm_version = false;
      }
    }
    else if (seek_chrm_mass) {
      if (word_count == 5 && strcmpCased(ln_contents[0], "mass", CaseSensitivity::NO) &&
          verifyContents(ln_contents[1], NumberFormat::INTEGER) &&
          verifyContents(ln_contents[2], NumberFormat::INTEGER) == false &&
          verifyContents(ln_contents[2], NumberFormat::STANDARD_REAL) == false &&
          verifyContents(ln_contents[2], NumberFormat::SCIENTIFIC) == false &&
          (verifyContents(ln_contents[3], NumberFormat::STANDARD_REAL) ||
           verifyContents(ln_contents[3], NumberFormat::SCIENTIFIC))) {
        seek_chrm_mass = false;
      }
    }

    // Check for GROMACS hallmarks
    if (word_count == 3 && ln_contents[0].length() == 1 && ln_contents[2].length() == 1 &&
        ln_contents[0][0] == '[' && ln_contents[2][0] == ']') {
      if (strcmpCased(ln_contents[1], "defaults", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "atoms", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "atomtypes", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "bonds", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "bondtypess", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "angles", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "angletypes", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "dihedrals", CaseSensitivity::NO) ||
          strcmpCased(ln_contents[1], "dihedraltypes", CaseSensitivity::NO)) {
        seek_gmx_dir = false;
      }
    }
    else if (word_count >= 2 && ln_contents[0] == "#include") {
      seek_gmx_incl = false;
    }

    // Increment the line counter
    line_idx++;
  }

  // Determine whether any format can be identified
  if (seek_ambr_version == false && seek_ambr_pointers == false) {
    return TopologyKind::AMBER;
  }
  if (seek_chrm_version == false && seek_chrm_mass == false) {
    return TopologyKind::CHARMM;
  }
  if (seek_gmx_dir == false && seek_gmx_incl == false) {
    return TopologyKind::GROMACS;
  }

  // No other topology formats are recognizable  
  rtErr("No known topology format was identified in \"" + tf.getFileName() + "\".",
        "detectTopologyKind");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
TopologyKind detectTopologyKind(const std::string &file_name, const std::string &caller) {
  const TextFile tf = readFileSample(file_name, caller, 2048);
  return detectTopologyKind(tf);
}
  
//-------------------------------------------------------------------------------------------------
std::string getDefaultFileExtension(const CoordinateFileKind kind) {
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
    return std::string(default_amber_crd_extension);
  case CoordinateFileKind::AMBER_INPCRD:
    return std::string(default_amber_inpcrd_extension);
  case CoordinateFileKind::AMBER_ASCII_RST:
    return std::string(default_amber_ascii_rst_extension);
  case CoordinateFileKind::AMBER_NETCDF:
    return std::string(default_amber_netcdf_extension);
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return std::string(default_amber_netcdf_rst_extension);
  case CoordinateFileKind::SDF:
    return std::string(default_sd_file_extension);
  case CoordinateFileKind::PDB:
    return std::string(default_pdb_file_extension);
  case CoordinateFileKind::UNKNOWN:
    rtErr("No default extension is available for files of " + getEnumerationName(kind) + " type.",
          "getDefaultFileExtension");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind inferCoordinateFileKind(const std::string &file_name) {
  std::string before, after;
  splitPath(file_name, &before, &after);
  if (strcmpCased(after, std::string(default_amber_crd_extension), CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_CRD;
  }
  else if (strcmpCased(after, std::string(default_amber_inpcrd_extension), CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_INPCRD;
  }
  else if (strcmpCased(after, std::string(default_amber_ascii_rst_extension),
                       CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_ASCII_RST;
  }
  else if (strcmpCased(after, std::string(default_amber_netcdf_extension), CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_NETCDF;
  }
  else if (strcmpCased(after, std::string(default_amber_netcdf_rst_extension),
                       CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_NETCDF_RST;
  }
  else if (strcmpCased(after, std::string(default_sd_file_extension), CaseSensitivity::NO)) {
    return CoordinateFileKind::SDF;
  }
  else if (strcmpCased(after, std::string(default_pdb_file_extension), CaseSensitivity::NO)) {
    return CoordinateFileKind::PDB;
  }
  else {
    return CoordinateFileKind::UNKNOWN;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string assembleFilePath(const std::string &root_path, const std::string &base_path,
                             const std::string &extn) {
  std::string result;
  if (root_path.size() > 0) {
    result = root_path + osSeparator();
  }
  result += base_path;
  if (extn.size() > 0) {
    result += "." + extn;
  }
  return result;
}

} // namespace diskutil
} // namespace stormm
