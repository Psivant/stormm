#include "copyright.h"

#include "Chemistry/znumber.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "mdlmol_atomlist.h"

namespace stormm {
namespace structure {

using chemistry::symbolToZNumber;
using constants::CaseSensitivity;
using constants::ExceptionResponse;
using diskutil::getBaseName;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::strncmpCased;
using parse::verifyContents;

//-------------------------------------------------------------------------------------------------
MdlMolAtomList::MdlMolAtomList(const std::vector<int> &atomic_numbers_in, bool const exclusions_in,
                               const int atom_attachment_in) :
  entry_count{static_cast<int>(atomic_numbers_in.size())},
  atomic_numbers{atomic_numbers_in},
  exclusions{exclusions_in},
  atom_attachment{atom_attachment_in}
{}

//-------------------------------------------------------------------------------------------------
MdlMolAtomList::MdlMolAtomList(const TextFile &tf, const int line_number,
                               const std::string &title) :
    MdlMolAtomList()
{
  // Determine whether this originates in an atom list entry or a property of the MDL section
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  if (lnlength < 6) {
    rtErr("Line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) +
          " cannot contain MDL MOL atom list information due to its length being only " +
          std::to_string(lnlength) + ".");
  }
  if (strncmpCased(line_ptr, std::string("M  ALS"))) {

    // The entry originates as a property.  This would invalidate any preceding atom list entries
    // of the deprecated format, but such a contingency has already been taken care of in the
    // parent MdlMolObj reader.
    if (lnlength < 15) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            getBaseName(tf.getFileName()) + " to read a property-based MDL MOL atom list.",
            "MdlMolAtomList");
    }
    if (verifyContents(line_ptr, 7, 3, NumberFormat::INTEGER)) {
      atom_attachment = readIntegerValue(line_ptr, 7, 3);
    }
    if (verifyContents(line_ptr, 10, 3, NumberFormat::INTEGER)) {
      entry_count = readIntegerValue(line_ptr, 10, 3);
    }
    if (entry_count < 0 || entry_count > 16) {
      rtErr("An invalid number of entries, " + std::to_string(entry_count) + ", was recorded at "
            "line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) + ".",
            "MdlMolAtomList");
    }
    switch (line_ptr[14]) {
    case 'T':
      exclusions = true;
      break;
    case 'F':
      exclusions = false;
      break;
    default:
      rtErr("The exclusions flag at line " + std::to_string(line_number) + " of " +
            getBaseName(tf.getFileName()) + ", " + line_ptr[14] + ", is invalid.",
            "MdlMolAtomList");
    }
    if (lnlength < 16 + (4 * entry_count)) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            getBaseName(tf.getFileName()) + " to read an atom list with " +
            std::to_string(entry_count) + "items.", "MdlMolAtomList");
    }
    std::vector<char4> tmp_symbols;
    for (int i = 0; i < entry_count; i++) {
      tmp_symbols[i].x = line_ptr[16 + (4 * i)];
      tmp_symbols[i].y = line_ptr[17 + (4 * i)];
      tmp_symbols[i].z = line_ptr[18 + (4 * i)];
      tmp_symbols[i].w = line_ptr[19 + (4 * i)];
    }
    atomic_numbers = symbolToZNumber(tmp_symbols, CaseSensitivity::YES, ExceptionResponse::DIE);
  }
  else {
    if (lnlength < 10) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            getBaseName(tf.getFileName()) + " to read one of the (deprecated) MDL MOL atom list "
            "entries.", "MdlMolAtomList");
    }

    // Read the deprecated atom list entry format
    if (verifyContents(line_ptr, 0, 3, NumberFormat::INTEGER)) {
      atom_attachment = readIntegerValue(line_ptr, 0, 3);
    }
    switch (line_ptr[4]) {
    case 'T':
      exclusions = true;
      break;
    case 'F':
      exclusions = false;
      break;
    default:
      rtErr("The exclusions flag at line " + std::to_string(line_number) + " of " +
            getBaseName(tf.getFileName()) + ", " + line_ptr[4] + ", is invalid.",
            "MdlMolAtomList");
    }
    if (verifyContents(line_ptr, 5, 5, NumberFormat::INTEGER)) {
      entry_count = readIntegerValue(line_ptr, 5, 5);
    }
    if (entry_count < 0 || entry_count > 5) {
      rtErr("An invalid number of entries, " + std::to_string(entry_count) + ", was recorded at "
            "line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) + ".",
            "MdlMolAtomList");
    }
    if (lnlength < 10 + (4 * entry_count)) {
      rtErr("There is insufficient information on line " + std::to_string(line_number) + " of " +
            getBaseName(tf.getFileName()) + " to read an atom list with " +
            std::to_string(entry_count) + "items.", "MdlMolAtomList");
    }
    atomic_numbers.resize(entry_count);
    for (int i = 0; i < entry_count; i++) {
      if (verifyContents(line_ptr, 11 + (4 * i), 3, NumberFormat::INTEGER)) {
        atomic_numbers[i] = readIntegerValue(line_ptr, 11 + (4 * i), 3);
      }
      else {
        rtErr("An invalid entry was found on line " + std::to_string(line_number) + " of " +
              getBaseName(tf.getFileName()) + ": " +
              tf.extractString(line_number, 11 + (4 * i), 3) + ".  The entries of one of the "
              "(deprecated) MDL MOL format atom lists must be integers corresponding to atomic "
              "numbers.  Use IUPAC element symbols in the new format, placed on MOL property "
              "lines beginning \"M  ALS\".", "MdlMolAtomList");
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
int MdlMolAtomList::getEntryCount() const {
  return entry_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMolAtomList::getEntry(const int index) const {
  if (index < 0 || index >= entry_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for an atom list containing " +
          std::to_string(entry_count) + " entries.", "MdlMolAtomList", "getEntry");
  }
  return atomic_numbers[index];
}

//-------------------------------------------------------------------------------------------------
bool MdlMolAtomList::applyToExclusions() const {
  return exclusions;
}

//-------------------------------------------------------------------------------------------------
char MdlMolAtomList::getExclusionCode() const {
  return (exclusions) ? 'T' : 'F';
}

//-------------------------------------------------------------------------------------------------
int MdlMolAtomList::getAttachmentPoint() const {
  return atom_attachment;
}

} // namespace structure
} // namespace stormm
