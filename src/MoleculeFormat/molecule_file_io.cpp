#include "copyright.h"
#include "Parsing/ascii_numbers.h"
#include "molecule_file_io.h"

namespace stormm {
namespace structure {

using parse::readIntegerValue;
using parse::readRealValue;

//-------------------------------------------------------------------------------------------------
int getMdlMolSectionEnd(const TextFileReader &tfr, const int line_start, const int line_end_in) {
  const int line_end = (line_end_in < 0) ? tfr.line_count : line_end_in;
  int result = line_start;
  bool found = false;
  while (found == false && result < line_end) {
    const char* lptr = &tfr.text[tfr.line_limits[result]];
    if (tfr.line_limits[result + 1] - tfr.line_limits[result] >= 6 &&
        lptr[0] == 'M' && lptr[1] == ' ' && lptr[2] == ' ' && lptr[3] == 'E' && lptr[4] == 'N' &&
        lptr[5] == 'D') {
      found = true;
    }
    else {
      result++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int getMdlMolSectionEnd(const TextFile &tf, const int line_start, const int line_end_in) {
  return getMdlMolSectionEnd(tf.data(), line_start, line_end_in);
}

//-------------------------------------------------------------------------------------------------
int getCompoundSectionEnd(const TextFileReader &tfr, const int line_start, const int line_end_in) {
  const int line_end = (line_end_in < 0) ? tfr.line_count : line_end_in;
  int result = line_start;
  bool found = false;
  while (found == false && result < line_end) {
    const char* lptr = &tfr.text[tfr.line_limits[result]];
    if (tfr.line_limits[result + 1] - tfr.line_limits[result] >= 6 &&
        lptr[0] == '$' && lptr[1] == '$' && lptr[2] == '$' && lptr[3] == '$') {
      found = true;
    }
    else {
      result++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int getCompoundSectionEnd(const TextFile &tf, const int line_start, const int line_end_in) {
  return getCompoundSectionEnd(tf.data(), line_start, line_end_in);
}

//-------------------------------------------------------------------------------------------------
MdlMolVersion findMolObjVersion(const char* text, const int nchar) {
  for (int i = 0; i < nchar; i++) {
    if (text[i] == 'V' || text[i] == 'v') {
      if (i < nchar - 4 && text[i + 2] == '0' && text[i + 3] == '0' && text[i + 4] == '0') {
        if (text[i + 1] == '2') {
          return MdlMolVersion::V2000;
        }
        else if (text[i + 1] == '3') {
          return MdlMolVersion::V3000;
        }
      }
    }
  }
  return MdlMolVersion::UNKNOWN;
}

//-------------------------------------------------------------------------------------------------
MdlMolVersion findMolObjVersion(const TextFile &tf, const int line_number) {
  return findMolObjVersion(tf.getLinePointer(line_number), tf.getLineLength(line_number));
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> findSdfMolEntryLimits(const TextFile &tf) {
  const TextFileReader tfr = tf.data();
  int nsection = 0;
  int last_delimiter_line = -1;
  std::vector<int2> result;
  for (int i = 0; i < tfr.line_count; i++) {
    if (tfr.line_lengths[i] >= 4 &&
        tfr.text[tfr.line_limits[i]    ] == '$' && tfr.text[tfr.line_limits[i] + 1] == '$' &&
        tfr.text[tfr.line_limits[i] + 2] == '$' && tfr.text[tfr.line_limits[i] + 3] == '$') {

      // A MOL entry must have at least four lines, so any text preceding a $$$$ marker must be at
      // least four lines long to count as a valid entry worth parsing
      if ((last_delimiter_line >=  0 && i - last_delimiter_line >= 4) ||
          (last_delimiter_line == -1 && i >= 4)) {

        // Estimate the number of frames based on the length of the first
        if (nsection == 0) {
          const int nframe_est = (tfr.line_count / (i - last_delimiter_line)) + 1;
          result.reserve(nframe_est);
        }
        result.push_back({ last_delimiter_line + 1, i });
        nsection++;
      }
      last_delimiter_line = i;
    }
  }
  if (tfr.line_count - last_delimiter_line >= 4) {
    result.push_back({ last_delimiter_line, tfr.line_count });
  }
  result.shrink_to_fit();
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double3> extractSdfCoordinates(const TextFile &tf, const int line_start,
                                           const int line_end_in, const ExceptionResponse policy) {
  std::vector<double3> result;

  // Check that there is sufficient information to read the SDF
  const TextFileReader tfr = tf.data();
  const int line_end = getMdlMolSectionEnd(tfr, line_start, line_end_in);
  if (line_end - line_start <= 4) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("There is not enough information in " + tf.getFileName() + " to contain atoms in MDL "
            "MOL format.", "extractSdfCoordinates");
    case ExceptionResponse::WARN:
      rtWarn("There is not enough information in " + tf.getFileName() + " to contain atoms in MDL "
            "MOL format.  An empty array will be returned.", "extractSdfCoordinates");
      break;    
    case ExceptionResponse::SILENT:
      break;
    }
    return result;
  }

  // Get the version
  const MdlMolVersion version_no = findMolObjVersion(tf, line_start + 3);
  
  // Read the number of atoms from the counts line, then read coordinates for each atom
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int atom_count = readIntegerValue(&tfr.text[tfr.line_limits[line_start + 3]], 0, 3);
      int iatm = 0;
      const int ilim = line_start + 4 + atom_count;
      result.reserve(atom_count);
      for (int i = line_start + 4; i < ilim; i++) {
        const char* atom_line_ptr = &tfr.text[tfr.line_limits[i]];
        result[iatm] = { readRealValue(atom_line_ptr,  0, 10),
                         readRealValue(atom_line_ptr, 10, 10),
                         readRealValue(atom_line_ptr, 20, 10) };
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    rtErr("No valid MDL MOL version was detected.  Parsing in " + tf.getFileName() + " cannot "
          "proceed.");
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double3> extractSdfCoordinates(const TextFile &tf, const int frame_number,
                                           const ExceptionResponse policy) {
  const std::vector<int2> limits = findSdfMolEntryLimits(tf);
  if (frame_number < 0 || frame_number >= static_cast<int>(limits.size())) {
    rtErr("The requested frame index (" + std::to_string(frame_number) + ") does not exist in an "
          "SD file with " + std::to_string(limits.size()) + " MDL MOL entries.",
          "extractSdfCoordinates");
  }
  return extractSdfCoordinates(tf, limits[frame_number].x, limits[frame_number].y);
}
  
} // namespace structure
} // namespace stormm
