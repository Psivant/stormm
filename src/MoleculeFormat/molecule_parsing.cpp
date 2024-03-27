#include "copyright.h"
#include "Parsing/parse.h"
#include "molecule_parsing.h"

namespace stormm {
namespace structure {

using parse::separateText;
using parse::stringToChar4;

//-------------------------------------------------------------------------------------------------
std::string condenseSdfItemDataLines(const std::string &item_name, const MdlMol &molecule) {
  std::vector<std::string> data_lines = molecule.getDataItemContent(item_name);
  std::string result;
  const size_t nlines = data_lines.size();
  if (nlines == 0LLU) {
    return result;
  }
  size_t total_chars = nlines - 1LLU;
  for (size_t i = 0LLU; i < nlines; i++) {
    total_chars += data_lines[i].size();
  }
  if (total_chars == 0LLU) {
    return result;
  }
  result.resize(total_chars, '\0');
  size_t cc = 0LLU;
  for (size_t i = 0LLU; i < nlines; i++) {
    const size_t nchar = data_lines[i].size();
    for (size_t j = 0LLU; j < nchar; j++) {
      result[cc] = (data_lines[i][j] == '\n') ? ' ' : data_lines[i][j];
      cc++;
    }
    if (data_lines[i][nchar - 1LLU] != '\n' && data_lines[i][nchar - 1LLU] != ' ') {
      result[cc] = ' ';
      cc++;
    }
  }
  result.resize(cc);
  result.shrink_to_fit();
  return result;
}
  
//-------------------------------------------------------------------------------------------------
AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag, const ChemicalFeatures &chemfe,
                             const CoordinateFrame &cf, const ExceptionResponse policy) {

  // Concatenate the strings, sans carriage returns if they have not already been removed.
  const std::string all_data = condenseSdfItemDataLines(item_name, molecule);
  
  // If there is no data item content, return an empty mask.
  if (all_data.size() == 0LLU) {
    return AtomMask(ag);
  }

  // Try to construct an AtomMask.
  try {
    AtomMask result(all_data, ag, chemfe, cf.data());
    return result;
  }
  catch (std::runtime_error) {

    // If the data item content could not be processed as an atom mask, try parsing it as as list
    // of integers, which would be assumed to indicate atom indices (starting at 1, and thus
    // decremented for C / C++ array indexing).  The text is assumed not to contain comments or
    // quoted strings.  Commas will be counted as delimiters.
    const std::vector<std::string> words = separateText(all_data,
                                                        std::vector<std::string>(1, ","));
    const int word_count = words.size();
    bool all_integer = true;
    bool all_atom_names = true;
    for (int i = 0; i < word_count; i++) {
      all_integer = (all_integer && verifyNumberFormat(words[i].c_str(), NumberFormat::INTEGER));
      all_atom_names = (all_atom_names && words[i].size() <= 4LLU);
    }
    AtomMask result(ag);
    if (all_integer) {
      std::vector<int> atom_adds(word_count);
      for (int i = 0; i < word_count; i++) {
        atom_adds[i] = stoi(words[i]) - 1;
      }
      result.addAtoms(atom_adds, policy);
    }
    else if (all_atom_names) {
      std::vector<char4> atom_adds(word_count);
      for (int i = 0; i < word_count; i++) {
        atom_adds[i] = stringToChar4(words[i]);
      }
      result.addAtoms(atom_adds);
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The text in data item <" + item_name + "> of SD file entry " + molecule.getTitle() +
              " could not be interpreted as a valid atom mask string, a list of integers, or a "
              "list of atom names.", "maskFromSdfDataItem");
      case ExceptionResponse::WARN:
        rtWarn("The text in data item <" + item_name + "> of SD file entry " +
               molecule.getTitle() + " could not be interpreted as a valid atom mask string, a "
               "list of integers, or a list of atom names.  No atoms will be selected.",
               "maskFromSdfDataItem");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    return result;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag, const ExceptionResponse policy) {
  const CoordinateFrame cf = molecule.exportCoordinateFrame();
  const ChemicalFeatures chemfe(ag, cf.data());
  return maskFromSdfDataItem(item_name, molecule, ag, chemfe, cf, policy);
}

//-------------------------------------------------------------------------------------------------
std::vector<PolyNumeric> valuesFromSdfDataItem(const std::string &item_name,
                                               const MdlMol &molecule, const NumberFormat fmt,
                                               const ExceptionResponse policy) {
  const std::string all_data = condenseSdfItemDataLines(item_name, molecule);
  const std::vector<std::string> words = separateText(all_data, std::vector<std::string>(1, ","));
  const int word_count = words.size();
  std::vector<PolyNumeric> result;
  result.resize(word_count);
  int non_conformity = -1;
  int nvalid = 0;
  for (int i = 0; i < word_count; i++) {
    if (verifyNumberFormat(words[i].c_str(), fmt) == false) {
      if (non_conformity < 0) {
        non_conformity = i;
      }
    }
    else {
      switch (fmt) {
      case NumberFormat::SCIENTIFIC:
      case NumberFormat::STANDARD_REAL:
        result[nvalid].d = stod(words[i]);
        break;
      case NumberFormat::INTEGER:
        result[nvalid].i = stoi(words[i]);
        break;
      case NumberFormat::LONG_LONG_INTEGER:
        result[nvalid].lli = stoll(words[i]);
        break;
      case NumberFormat::UNSIGNED_INTEGER:
        result[nvalid].ui = stoi(words[i]);
        break;
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
        result[nvalid].ulli = stoll(words[i]);
        break;
      case NumberFormat::CHAR4:
        result[nvalid].c4 = stringToChar4(words[i]);
        break;
      }
      nvalid++;
    }
  }
  if (nvalid != word_count) {
    result.resize(nvalid);
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A non-integer expression, \"" + words[non_conformity] + "\", was found in the data "
            "item content.", "intFromSdfDataItem");
    case ExceptionResponse::WARN:
      rtWarn("One or more non-integer expressions, e.g. \"" + words[non_conformity] + "\", were "
             "found in the content of data item \"" + item_name + "\" serving MDL MOL entry \"" +
             molecule.getTitle() + "\".  These elements will not translate into integers.  No "
             "comments are interpreted in SD file data item content.", "integerFromSdfDataItem");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> realFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                                        const ExceptionResponse policy) {
  const std::vector<PolyNumeric> as_scinot = valuesFromSdfDataItem(item_name, molecule,
                                                                   NumberFormat::SCIENTIFIC,
                                                                   ExceptionResponse::SILENT);
  const std::vector<PolyNumeric> as_stdreal = valuesFromSdfDataItem(item_name, molecule,
                                                                    NumberFormat::STANDARD_REAL,
                                                                    ExceptionResponse::SILENT);
  return (as_scinot.size() > as_stdreal.size()) ? doubleFromPolyNumeric(as_scinot) :
                                                  doubleFromPolyNumeric(as_stdreal);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> intFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                                    const ExceptionResponse policy) {
  return intFromPolyNumeric(valuesFromSdfDataItem(item_name, molecule, NumberFormat::INTEGER,
                                                  policy));
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> char4FromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                                        const ExceptionResponse policy) {
  return char4FromPolyNumeric(valuesFromSdfDataItem(item_name, molecule, NumberFormat::CHAR4,
                                                    policy));
}

} // namespace structure
} // namespace stormm
