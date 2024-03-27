// -*-c++-*-
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/vector_ops.h"
#include "Parsing/ascii_numbers.h"
#include "file_snapshot.h"

namespace stormm {
namespace testing {

using data_types::getStormmScalarTypeName;
using data_types::getStormmHpcVectorTypeName;
using parse::printNumberSeries;
using parse::readNumberSeries;
using parse::separateText;
using parse::TextFileReader;
using parse::TextOrigin;
using stmath::maxAbsValue;
using stmath::minValue;

//-------------------------------------------------------------------------------------------------
std::vector<PolyNumeric> readSnapshot(const TextFile &tf, const std::string &label) {

  // Seek the formatting comment
  TextFileReader tfr = tf.data();
  bool label_found = false;
  bool format_found = false;
  bool bad_type_name = false;
  int iline = 0;
  int n_values, values_per_line, width, decimal;
  NumberFormat data_format;
  std::string tname;
  while (format_found == false && iline < tfr.line_count - 1) {
    if (tfr.line_limits[iline + 1] - tfr.line_limits[iline] >= 17 &&
        tfr.text[tfr.line_limits[iline]] == '%' &&
        strncmp(&tfr.text[tfr.line_limits[iline] + 2], "Snapshot data, ", 15) == 0) {

      // This looks promising.  Separate the words on the line and seek the name of the data type.
      std::vector<std::string> snapshot_info = separateText(&tfr.text[tfr.line_limits[iline]],
                                                            tfr.line_limits[iline + 1] -
                                                            tfr.line_limits[iline]);
      const int n_info = snapshot_info.size();
      n_values = stoi(snapshot_info[3]);
      bool problem = false;
      bool type_name_matches = false;
      if (n_info >= 10) {
        for (int i = 0; i < n_info - 4; i++) {
          if (snapshot_info[i] == "{") {
            values_per_line = stoi(snapshot_info[i + 2]);
            width = stoi(snapshot_info[i + 3]);
            if (snapshot_info[i + 1] == "GENERAL" || snapshot_info[i + 1] == "SCIENTIFIC") {
              decimal = stoi(snapshot_info[i + 4]);
              problem = (snapshot_info[i + 5] != "}" || problem);
              if (snapshot_info[i + 1] == "GENERAL") {
                data_format = NumberFormat::STANDARD_REAL;
              }
              else {
                data_format = NumberFormat::SCIENTIFIC;
              }
              tname = getStormmScalarTypeName<double>();
            }
            else if (snapshot_info[i + 1] == "INTEGER" ||
                     snapshot_info[i + 1] == "UNSIGNED_INTEGER" ||
                     snapshot_info[i + 1] == "LONG_LONG_INTEGER" ||
                     snapshot_info[i + 1] == "UNSIGNED_LONG_LONG_INTEGER" ||
                     snapshot_info[i + 1] == "CHAR4") {
              decimal = 0;
              problem = (snapshot_info[i + 4] != "}" || problem);
              if (snapshot_info[i + 1] == "INTEGER") {
                data_format = NumberFormat::INTEGER;
                tname = getStormmScalarTypeName<int>();
              }
              else if (snapshot_info[i + 1] == "UNSIGNED_INTEGER") {
                data_format = NumberFormat::UNSIGNED_INTEGER;
                tname = getStormmScalarTypeName<uint>();
              }
              else if (snapshot_info[i + 1] == "LONG_LONG_INTEGER") {
                data_format = NumberFormat::LONG_LONG_INTEGER;
                tname = getStormmScalarTypeName<llint>();
              }
              else if (snapshot_info[i + 1] == "UNSIGNED_LONG_LONG_INTEGER"){
                data_format = NumberFormat::UNSIGNED_LONG_LONG_INTEGER;
                tname = getStormmScalarTypeName<ullint>();
              }
	      else {
                data_format = NumberFormat::CHAR4;
                tname = getStormmHpcVectorTypeName<char4>();
	      }
            }
            else {
              problem = true;
            }
          }
        }
	type_name_matches = (snapshot_info[4] == tname);
        
        // This is looking better.  Check the label, if provided, and if
        // the snapshot does not take the default name of "data."
        if (problem == false) {
          iline++;
          snapshot_info = separateText(&tfr.text[tfr.line_limits[iline]],
                                       tfr.line_limits[iline + 1] - tfr.line_limits[iline]);
          if ((label.size() == 0 && snapshot_info[0] == "data") ||
              (label.size() > 0 && snapshot_info[0] == label)) {
            label_found = true;
            if (type_name_matches) {
              format_found = true;
            }
            else {
              bad_type_name = true;
            }
          }
        }
      }
      else {
        problem = true;
      }
      if (problem) {
        std::string format_line;
        for (int j = 0; j < n_info; j++) {
          format_line += snapshot_info[j];
        }
        rtErr("Corrupted snapshot formatting line \"" + format_line + "\".", "readSnapshot");
      }
    }

    // The line should increment whether the correct format and variable were found or not
    iline++;
  }

  // If the format was not found, try to explain why
  if (format_found == false) {
    if (bad_type_name && label.size() > 0) {
      rtErr("A search for \"" + label + "\" in snapshot file " + tf.getFileName() + " yielded a "
            "match but for the correct type (" + tname + " required).", "readSnapshot");
    }
    else if (label_found == false) {
      rtErr("A search for \"" + label + "\" in snapshot file " + tf.getFileName() + " yielded no "
            "matches.", "readSnapshot");
    }
    else {
      rtErr("The proper formatting could not be identified in snapshot file " + tf.getFileName() +
            ".", "readSnapshot");
    }
  }

  // Take in the data
  return readNumberSeries(tf, iline, n_values, values_per_line, width, decimal, data_format,
                          "readSnapshot", "Read snapshot data");
}

//-------------------------------------------------------------------------------------------------
std::vector<PolyNumeric> readSnapshot(const std::string &filename, const std::string &label) {
  const TextFile tf(filename, TextOrigin::DISK, std::string(""), "readSnapshot");
  return readSnapshot(tf, label);
}

//-------------------------------------------------------------------------------------------------
void writeSnapshot(const std::string &filename, const std::vector<PolyNumeric> &content,
                   const std::string &label, const double tol, const NumberFormat data_format,
                   const PrintSituation expectation) {

  // Open the output file
  std::ofstream foutp = openOutputFile(filename, expectation, "Write a snapshot of data for " +
                                       label);

  // Vectors that can hold a proper interpretation of the content
  std::vector<int> icontent;
  std::vector<uint> uicontent;
  std::vector<llint> llicontent;
  std::vector<ullint> ullicontent;
  std::vector<double> dcontent;
  
  // Determine a format for printing this type
  const int necessary_buffer_size = ((2*label.size() + 127) / 64) * 64;
  std::vector<char> buffer(necessary_buffer_size, '\0');
  int width, digits_after_decimal;
  std::string name_of_the_type;
  switch (data_format) {
  case NumberFormat::SCIENTIFIC:
    {
      dcontent = doubleFromPolyNumeric(content);
      const double maxv = maxAbsValue(dcontent) / tol;
      const double abs_maxv = fabs(maxv);
      const int n_digits = (abs_maxv > constants::tiny && fabs(log10(abs_maxv)) > 2) ?
                           fabs(log10(abs_maxv)) : 2;
      width = n_digits + 7;
      digits_after_decimal = n_digits - 1;
      name_of_the_type = getStormmScalarTypeName<double>();
    }
    break;
  case NumberFormat::STANDARD_REAL:
    {
      digits_after_decimal = ceil(fabs(log10(fabs(tol)))) + 0.01;
      dcontent = doubleFromPolyNumeric(content);
      const double maxv = maxAbsValue(dcontent);
      const int digits_before_decimal = (fabs(maxv) > 1.0) ?
                                        ceil(fabs(log10(fabs(maxv)))) + 0.01 : 1;
      width = digits_before_decimal + digits_after_decimal + 3;
      name_of_the_type = getStormmScalarTypeName<double>();
    }
    break;
  case NumberFormat::INTEGER:
    {
      icontent = intFromPolyNumeric(content);
      const int maxv = maxAbsValue(icontent);
      const int n_digits = (abs(maxv) > 0) ?
                           ceil(fabs(log10(fabs(static_cast<double>(maxv))))) + 0.01 : 2;
      width = n_digits + 2 + (minValue(icontent) < 0);
      digits_after_decimal = 0;
      name_of_the_type = getStormmScalarTypeName<int>();
    }
    break;
  case NumberFormat::LONG_LONG_INTEGER:
    {
      llicontent = llintFromPolyNumeric(content);
      const llint maxv = maxAbsValue(llicontent);
      const int n_digits = (std::llabs(maxv) > 0) ?
                           ceil(fabs(log10(fabs(static_cast<double>(maxv) * 1.1)))) + 0.01 : 2;
      width = n_digits + 2 + (minValue(icontent) < 0);
      digits_after_decimal = 0;
      name_of_the_type = getStormmScalarTypeName<llint>();
    } 
    break;
  case NumberFormat::UNSIGNED_INTEGER:
    {
      uicontent = uintFromPolyNumeric(content);
      const llint maxv = maxAbsValue(uicontent);
      const int n_digits = (std::llabs(maxv) > 0LL) ?
                           ceil(fabs(log10(fabs(static_cast<double>(maxv))))) + 0.01 : 2;
      width = n_digits + 1;
      digits_after_decimal = 0;
      name_of_the_type = getStormmScalarTypeName<uint>();
    }
    break;
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    {
      ullicontent = ullintFromPolyNumeric(content);
      const llint maxv = maxAbsValue(ullicontent);
      const int n_digits = (std::llabs(maxv) > 0LL) ?
                           ceil(fabs(log10(fabs(static_cast<double>(maxv))))) + 0.01 : 2;
      width = n_digits + 1;
      digits_after_decimal = 0;
      name_of_the_type = getStormmScalarTypeName<ullint>();
    }
    break;
  case NumberFormat::CHAR4:
    name_of_the_type = getStormmHpcVectorTypeName<char4>();
    width = 5;
    break;
  }

  // Output the format line
  const int values_per_line = (99 / width < 1) ? 1 : 99 / width;
  switch (data_format) {
  case NumberFormat::SCIENTIFIC:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { SCIENTIFIC %d %d "
             "%d } (%d x %%%d.%de per line).\n", content.size(), name_of_the_type.c_str(),
             values_per_line, width, digits_after_decimal, values_per_line, width,
             digits_after_decimal);
    break;
  case NumberFormat::STANDARD_REAL:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { GENERAL %d %d %d "
             "} (%d x %%%d.%df per line).\n", content.size(), name_of_the_type.c_str(),
             values_per_line, width, digits_after_decimal, values_per_line, width,
             digits_after_decimal);
    break;
  case NumberFormat::INTEGER:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { INTEGER %d %d } "
             "(%d x %%%dd per line).\n", content.size(), name_of_the_type.c_str(), values_per_line,
             width, values_per_line, width);
    break;
  case NumberFormat::LONG_LONG_INTEGER:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { LONG_LONG_INTEGER "
             "%d %d } (%d x %%%dlld per line).\n", content.size(), name_of_the_type.c_str(),
             values_per_line, width, values_per_line, width);
    break;
  case NumberFormat::UNSIGNED_INTEGER:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { UNSIGNED_INTEGER "
             "%d %d } (%d x %%%dlu per line).\n", content.size(), name_of_the_type.c_str(),
             values_per_line, width, values_per_line, width);
    break;
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { "
             "UNSIGNED_LONG_LONG_INTEGER %d %d } (%d x %%%dllu per line).\n", content.size(),
             name_of_the_type.c_str(), values_per_line, width, values_per_line, width);
    break;
  case NumberFormat::CHAR4:
    snprintf(buffer.data(), necessary_buffer_size, "%% Snapshot data, %zu %s { CHAR4 %d } (%d x "
             "%%4.4s per line).\n", content.size(), name_of_the_type.c_str(), values_per_line,
             values_per_line);
    break;
  }

  // Print the data, in sections if the vector contains HPC vector tuple types.  Pad the content
  // with zeros in the event that it would print an inconsistent number of values on the last line,
  // so that Matlab and Octave can read the data as a variable and reshape it to a vector just as
  // it was in the STORMM program's memory.
  foutp.write(buffer.data(), strlen(buffer.data()));
  const std::string var_label = (label.size() > 0) ? label : "data";
  snprintf(buffer.data(), necessary_buffer_size, "%s = [\n", var_label.c_str());
  foutp.write(buffer.data(), strlen(buffer.data()));
  const int nval = content.size();
  const int padded_nval = ((nval + values_per_line - 1) / values_per_line) * values_per_line;
  std::vector<PolyNumeric> content_plus(content.begin(), content.end());
  content_plus.resize(padded_nval);
  for (int i = nval; i < padded_nval; i++) {

    // Either unsigned long long int or double will be the largest type.  Set both to zero.
    content_plus[i].ulli = 0ull;
    content_plus[i].d = 0.0;
  }
  printNumberSeries(&foutp, content_plus, values_per_line, width, digits_after_decimal,
                    data_format, "writeSnapshot", "Snapshot program data");
  snprintf(buffer.data(), necessary_buffer_size, "];\n%s = reshape(transpose(%s), 1, %d);\n",
           var_label.c_str(), var_label.c_str(), padded_nval);
  foutp.write(buffer.data(), strlen(buffer.data()));
  if (nval != padded_nval) {
    snprintf(buffer.data(), necessary_buffer_size, "%s = %s(1:%d);\n", var_label.c_str(),
             var_label.c_str(), nval);
    foutp.write(buffer.data(), strlen(buffer.data()));
  }
  
  // Close the output file
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
TextFile readTextSnapshot(const TextFile &tf, const std::string &label) {

  // Seek the formatting comment
  TextFileReader tfr = tf.data();
  bool label_found = false;
  bool format_found = false;
  int iline = 0;
  while (format_found == false && iline < tfr.line_count - 1) {
    if (tfr.line_limits[iline + 1] - tfr.line_limits[iline] >= 11 &&
        tfr.text[tfr.line_limits[iline]] == '|' &&
        strncmp(&tfr.text[tfr.line_limits[iline]], "|>>> Label", 10) == 0) {
      
      // This looks promising.  Separate the words on the line and seek the name of the data type.
      std::vector<std::string> snapshot_info = separateText(&tfr.text[tfr.line_limits[iline] + 11],
                                                            tfr.line_limits[iline + 1] -
                                                            tfr.line_limits[iline] - 11);
      const int n_info = snapshot_info.size();
      if ((label.size() == 0 && snapshot_info[0] == "data") ||
          (label.size() > 0  && snapshot_info[0] == label)) {

        // The label has ben matched.  Read until the next "|>>> End" marker and assemble a new
        // TextFile to hold the result.
        iline++;
        const int trx_start = iline;
        TextFile result;
        while (iline < tfr.line_count &&
               (tfr.line_limits[iline + 1] - tfr.line_limits[iline] < 8 ||
                strncmp(&tfr.text[tfr.line_limits[iline]], "|>>> End", 8) != 0)) {
          iline++;
        }
        if (iline == tfr.line_count) {
          rtErr("File " +  tf.getFileName() + " has no end mark for data in label [" + label +
                "].  Transcription starts on line " + std::to_string(trx_start) + ".",
                "readTextSnapshot");
        }
        const int trx_end = iline;
        const size_t nchar = tfr.line_limits[trx_end] - tfr.line_limits[trx_start];
        std::string proto_result;
        proto_result.reserve(nchar + trx_end - trx_start);
        for (int i = trx_start; i < trx_end; i++) {
          const size_t jlim = tfr.line_limits[i + 1];
          for (size_t j = tfr.line_limits[i]; j < jlim; j++) {
            proto_result.push_back(tfr.text[j]);
          }
          proto_result.push_back('\n');
        }
        return TextFile(proto_result, TextOrigin::RAM);
      }
    }
    iline++;
  }
  rtErr("A search for \"" + label + "\" in snapshot file " + tf.getFileName() + " yielded no "
        "matches.", "readTextSnapshot");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
TextFile readTextSnapshot(const std::string &filename, const std::string &label) {
  const TextFile tf(filename, TextOrigin::DISK, std::string(""), "readTextSnapshot");
  return readTextSnapshot(tf, label);
}

//-------------------------------------------------------------------------------------------------
void writeTextSnapshot(const std::string &filename, const TextFile &content,
                       const std::string &label, const PrintSituation expectation) {
  std::ofstream foutp = openOutputFile(filename, expectation);
  std::string opening_line("|>>> Label ");
  opening_line += (label.size() == 0) ? "data" : label;
  opening_line += "\n";
  foutp.write(opening_line.data(), opening_line.size());
  content.write(&foutp);
  const std::string closing_line("|>>> End\n\n");
  foutp.write(closing_line.data(), closing_line.size());
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void writeTextSnapshot(const std::string &filename, const std::string &content,
                       const std::string &label, const PrintSituation expectation) {
  const TextFile tf(content, TextOrigin::RAM);
  writeTextSnapshot(filename, tf, label, expectation);
}

} // namespace testing
} // namespace stormm
