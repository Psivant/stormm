#include <cstring>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parse.h"
#include "Parsing/ascii_numbers.h"
#include "amber_prmtop_util.h"

namespace stormm {
namespace topology {

using parse::CaseSensitivity;
using parse::extractFormattedNumber;
using parse::separateText;
using parse::strncmpCased;
using parse::TextFileReader;
using parse::uppercase;
using parse::verifyNumberFormat;

//-------------------------------------------------------------------------------------------------
int4 parseLineFormat(const std::string &fmt_in, const std::string &file_name, const int line_idx) {
  const int nchar = fmt_in.size();
  int4 result;
  int i = 0;
  bool problem = false;

  // First, get the number of numbers in the format
  std::string nnum;
  while (i < nchar && fmt_in[i] >= '0' && fmt_in[i] <= '9') {
    nnum += fmt_in[i];
    i++;
  }
  if (nnum.size() == 0LLU) {
    result.x = 1;
  }
  else {
    result.x = stol(nnum);
  }

  // Next, get the format (one character) and check it
  result.y = uppercase(fmt_in[i]);
  i++;
  if (result.y != 'I' && result.y != 'E' && result.y != 'F' && result.y != 'A') {
    problem = true;
  }
  else {

    // Now get the length of each number
    std::string width;
    while (i < nchar && fmt_in[i] >= '0' && fmt_in[i] <= '9') {
      width += fmt_in[i];
      i++;
    }
    if (width.size() == 0LLU) {
      problem = true;
    }
    else {
      result.z = stol(width);

      // Finally, get the number of decimal places (if applicable)
      if (result.y == 'I' || result.y == 'A') {
        if (i != nchar) {
          problem = true;
        }
        else {
          result.w = 0;
        }
      }
      else if (result.y == 'E' || result.y == 'F') {
        if (fmt_in[i] != '.') {
          problem = true;
        }
        else {
          i++;
          std::string decimal;
          while (i < nchar && fmt_in[i] >= '0' && fmt_in[i] <= '9') {
            decimal += fmt_in[i];
            i++;
          }
          if (decimal.size() == 0LLU) {
            problem = true;
          }
          else {
            result.w = stol(decimal);
          }
        }
      }
    }
  }
  if (problem) {
    rtErr("Unrecognized format " + fmt_in + " on line " + std::to_string(line_idx) +
          " of file " + file_name + ".", "parseLineFormat");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int scanToFlag(const TextFile &tf, const char* flag, std::vector<int4> *detected_format,
               const TopologyRequirement needed, const int start_line) {
  TextFileReader tfr = tf.data();
  const int flag_len = strlen(flag);

  // Validate the given limits
  const int i_line = (start_line >= 0 && start_line <= tfr.line_count - 1) ? start_line : 0;

  // Check through the specified range
  for (int i = i_line; i < tfr.line_count - 1; i++) {

    // Advance through any whitespace, find two word starts
    int flstart = tfr.line_limits[i];
    while (flstart < tfr.line_limits[i + 1] && tfr.text[flstart] == ' ') {
      flstart++;
    }
    int scstart = flstart + 1;
    while (scstart < tfr.line_limits[i + 1] && tfr.text[scstart] != ' ') {
      scstart++;
    }
    while (scstart < tfr.line_limits[i + 1] && tfr.text[scstart] == ' ') {
      scstart++;
    }
    if (tfr.line_limits[i + 1] >= flstart + 5 &&
        strncmp(&tfr.text[flstart], "%FLAG", 5) == 0 &&
        tfr.line_limits[i + 1] >= scstart + flag_len &&
        strncmp(&tfr.text[scstart], flag, flag_len) == 0 &&
        (tfr.line_limits[i + 1] == scstart + flag_len || tfr.text[scstart + flag_len] == ' ')) {

      // Check the subsequent format
      for (int j = 1; j < tfr.line_count - 1; j++) {

        // Advance through any whitespace
        int fstart = tfr.line_limits[i + j];
        while (fstart < tfr.line_limits[i + j + 1] && tfr.text[fstart] == ' ') {
          fstart++;
        }
        if (tfr.line_limits[i + j + 1] - tfr.line_limits[i + j] >= 7 &&
            strncmpCased(&tfr.text[tfr.line_limits[i + j]], std::string("%FORMAT"),
                         CaseSensitivity::NO, 7)) {

          // Obtain a string of characters between '(' and ')'.  Strip away any further '(' and
          // ')' within the detected format string, i.e. (8(F9.5)) --> 8F9.5.  This is a hack but
          // there are limits to what an Amber prmtop can contain and multiple tuples of different
          // types of data on the same line are not in it (for now).
          int llim = tfr.line_limits[i + j] + 7;
          int hlim = tfr.line_limits[i + j + 1] - 1;
          while (llim < hlim && tfr.text[llim] != '(') {
            llim++;
          }
          llim++;
          while (hlim > llim && tfr.text[hlim] != ')') {
            hlim--;
          }
          hlim--;
          if (hlim > llim) {
            std::string proto_format;
            for (int k = llim; k <= hlim; k++) {
              if (tfr.text[k] != '(' && tfr.text[k] != ')') {
                proto_format += tfr.text[k];
              }
            }
            std::vector<std::string> format_commands = separateText(proto_format,
                                                                    { std::string(",") });
            detected_format->resize(format_commands.size());
            int4* dftmp = detected_format->data();
            for (size_t k = 0; k < format_commands.size(); k++) {
              dftmp[k] = parseLineFormat(format_commands[k], tf.getFileName(), i + j);
            }
          }
          else {
            detected_format->resize(1);
            detected_format->data()[0] = {0, 0, 0, 0};
          }
          return i + j + 1;
        }
      }
    }
  }

  // At this point, the flag must not have been found.  Recursively check the rest of the file.
  if (i_line > 0) {
    return scanToFlag(tf, flag, detected_format, needed);
  }

  // Any missing, essential flag throws an error
  if (needed == TopologyRequirement::ESSENTIAL) {
    std::string sflag(flag);
    rtErr("Unable to find %FLAG " + sflag + " within topology " + tf.getFileName() + ".",
          "scanToFlag");
  }

  // A missing, optional flag returns -1
  return -1;
}

//-------------------------------------------------------------------------------------------------
std::vector<Citation> readForceFieldReferences(const TextFile &tf, const int lstart) {
  std::vector<Citation> fflds;
  TextFileReader tfr = tf.data();
  int i = lstart;
  while (i < tfr.line_count && verifyNumberFormat(tfr.text, NumberFormat::INTEGER,
                                                  tfr.line_limits[i], 2)) {
    PolyNumeric pn = extractFormattedNumber(tfr.text, NumberFormat::INTEGER,
                                            tfr.line_limits[i], 2);
    std::string ffld_desc = "";
    for (int j = tfr.line_limits[i] + 2; j < tfr.line_limits[i + 1]; j++) {
      ffld_desc += tfr.text[j];
    }
    fflds.push_back(Citation(ffld_desc, pn.i));
    i++;
  }
  return fflds;
}

//-------------------------------------------------------------------------------------------------
std::vector<PolyNumeric> amberPrmtopData(const TextFile &tf, const int start_line,
                                         const NumberFormat cform, const int count_per_line,
                                         const int width, const int required_count,
                                         const int possible_count) {
  TextFileReader tfr = tf.data();
  int current_line = start_line;
  std::vector<char> buffer(width + 1, '\0');
  std::string tmp_line((count_per_line * width) + 1, '\0');

  // Look at the format of one of the numbers, if necessary, to get the number of decimal places
  int decimal;
  switch (cform) {
  case NumberFormat::SCIENTIFIC:
    {
      const int llim = tfr.line_limits[current_line];
      const int dot_index = '.';
      const int e_index = 'e';
      const int ecap_index = 'E';
      int dot_position = 0;
      int e_position = 0;
      for (int i = 0; i < width; i++) {
	const int tmpc = tfr.text[llim + i];
	if (tmpc == dot_index) {
	  dot_position = i;
	}
	else if (tmpc == e_index || tmpc == ecap_index) {
	  e_position = i;
	}
      }
      decimal = e_position - dot_position - 1;

      // If the number is misshapen with an E before the dot, don't let that confuse the output
      decimal -= (decimal < 0) * decimal;
    }
    break;
  case NumberFormat::STANDARD_REAL:
    {
      const int llim = tfr.line_limits[current_line];
      const int dot_index = '.';
      int dot_position = 0;
      for (int i = 0; i < width; i++) {
	const int tmpc = tfr.text[llim + i];
	if (tmpc == dot_index) {
	  dot_position = i;
	}
      }
      decimal = width - dot_position - 1;
    }
    break;
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
  case NumberFormat::CHAR4:
    decimal = 0;
    break;
  }

  // Look at the number of viable lines, which are determine by the point at which the end of the
  // file is reached or one begins with whitespace up to a '%'
  int n_viable_lines = 0;
  int count_so_far = 0;
  while (current_line < tfr.line_count &&
         verifyNumberFormat(tfr.text, cform, tfr.line_limits[current_line], width)) {
    const int est_count_this_line = (cform == NumberFormat::CHAR4) ?
                                    (tfr.line_limits[current_line + 1] + (width - 1) -
                                     tfr.line_limits[current_line]) / width :
                                    (tfr.line_limits[current_line + 1] -
                                     tfr.line_limits[current_line]) / width;
    count_so_far += est_count_this_line;

    // CHAR4 is tricky: anything is a potential char4 entry, but a new section or a blank line
    // following the completion of the expected data should end the reading.
    if (cform == NumberFormat::CHAR4 &&
        (count_so_far >= possible_count &&
         (est_count_this_line == 0 ||
          (tfr.line_limits[current_line + 1] - tfr.line_limits[current_line] > 5 &&
           tfr.text[tfr.line_limits[current_line]] == '%' &&
           tfr.text[tfr.line_limits[current_line] + 1] == 'F' &&
           tfr.text[tfr.line_limits[current_line] + 2] == 'L' &&
           tfr.text[tfr.line_limits[current_line] + 3] == 'A' &&
           tfr.text[tfr.line_limits[current_line] + 4] == 'G')))) {
      break;
    }

    // A short line prior to the completion of the series is an error
    if (est_count_this_line < count_per_line && count_so_far < required_count) {
      rtErr("A request for " + std::to_string(required_count) + " entries of type " +
	    nameNumberFormat(cform) + " in " + tf.getFileName() + " cannot be fullfilled due to "
	    "a short line containing at most " + std::to_string(est_count_this_line) + " entries "
	    "(line " + std::to_string(current_line + 1) + ").", "amberPrmtopData");
    }
    n_viable_lines++;
    current_line++;
  }

  // Check the last line, as it may not contain a full complement of entries
  current_line--;
  int n_on_last_line = 0;
  int check_pos = tfr.line_limits[current_line];
  while (check_pos < tfr.line_limits[current_line + 1]) {
    n_on_last_line += verifyNumberFormat(tfr.text, cform, check_pos, width);
    check_pos += width;
  }
  const int available_count = ((n_viable_lines > 0) * (n_viable_lines - 1) * count_per_line) +
                              n_on_last_line;
  if (available_count < required_count) {
    rtErr("A request for " + std::to_string(required_count) + " entries of type " +
          nameNumberFormat(cform) + " in " + tf.getFileName() + " cannot be fullfilled as the "
	  "file has at most " + std::to_string(available_count) + " entries of the specified "
	  "format (starting on line " + std::to_string(start_line + 1) + ") before reaching a "
          "new section or the end of the file.", "amberPrmtopData");
  }
  else if (possible_count > -1 && available_count > possible_count) {
    rtErr("A request for " + std::to_string(required_count) + " entries of type " +
          nameNumberFormat(cform) + " in " + tf.getFileName() + " finds " +
	  std::to_string(available_count) + " entries.  This is greater than the maximum number "
	  "of allowed entries (" + std::to_string(possible_count) + ").", "amberPrmtopData");
  }
  return readNumberSeries(tf, start_line, available_count, count_per_line, width, decimal, cform);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> iAmberPrmtopData(const TextFile &tf, const int start_line,
                                  const int count_per_line, const int width,
                                  const int required_count, const int possible_count) {
  return intFromPolyNumeric(amberPrmtopData(tf, start_line, NumberFormat::INTEGER, count_per_line,
					    width, required_count, possible_count));
}

//-------------------------------------------------------------------------------------------------
std::vector<double> dAmberPrmtopData(const TextFile &tf, const int start_line,
                                     const int count_per_line, const int width,
                                     const int required_count, const int possible_count) {
  return doubleFromPolyNumeric(amberPrmtopData(tf, start_line, NumberFormat::STANDARD_REAL,
                                               count_per_line, width, required_count,
                                               possible_count));
}

//-------------------------------------------------------------------------------------------------
std::vector<double> eAmberPrmtopData(const TextFile &tf, const int start_line,
                                     const int count_per_line, const int width,
                                     const int required_count, const int possible_count) {
  return doubleFromPolyNumeric(amberPrmtopData(tf, start_line, NumberFormat::SCIENTIFIC,
                                               count_per_line, width, required_count,
                                               possible_count));
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> c4AmberPrmtopData(const TextFile &tf, const int start_line,
                                     const int count_per_line, const int required_count,
                                     const int possible_count) {
  return char4FromPolyNumeric(amberPrmtopData(tf, start_line, NumberFormat::CHAR4, count_per_line,
                                              4, required_count, possible_count));
}

} // namespace topology
} // namespace stormm
