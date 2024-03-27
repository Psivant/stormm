// -*-c++-*-
#ifndef STORMM_MOLECULE_PARSING_H
#define STORMM_MOLECULE_PARSING_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/atommask.h"
#include "Chemistry/chemical_features.h"
#include "MoleculeFormat/mdlmol.h"
#include "Parsing/polynumeric.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"

namespace stormm {
namespace structure {

using chemistry::AtomMask;
using chemistry::ChemicalFeatures;
using parse::NumberFormat;
using parse::PolyNumeric;
using structure::MdlMol;
using topology::AtomGraph;
using trajectory::CoordinateFrame;

/// \brief Concatenate the data lines in an SD file data item, and filter out carriage returns.
///
/// \param item_name  Name of the data item
/// \param molecule   The MDL MOL entry with appended data items
std::string condenseSdfItemDataLines(const std::string &item_name, const MdlMol &molecule);
  
/// \brief Derive an atom mask (a list of atoms) from information in the data lines of an SD file
///        data item.
///
/// Overloaded:
///   - Accept a pre-made CoordinateFrame object (for AtomMask production) along with the MDL MOL
///     object and pre-determined chemical features along with the topology
///   - Extract the coordinates from the MDL MOL object and the chemical features from the topology
///     as temporary resources
///
/// \param item_name  Name of the data item
/// \param molecule   The MDL MOL entry with appended data items
/// \param ag         Topology of the system of interest, provided for testing whether the text is
///                   a valid atom mask or a list of atom names or numbers
/// \param chemfe     Chemical features detected for the topology of interest
/// \{
AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag, const ChemicalFeatures &chemfe,
                             const CoordinateFrame &cf,
                             ExceptionResponse policy = ExceptionResponse::DIE);

AtomMask maskFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                             const AtomGraph *ag,
                             ExceptionResponse policy = ExceptionResponse::DIE);
/// \}

/// \brief Return a set of values from an SD file data item's content.  The values will be assumed
///        to make up the entirety of the content, with no comments (although quoted strings are
///        allowed up to the point that their contents conform to one of the alphanumeric strings
///        that could represent one of the enumerated NumberFormat types).  The result is returned
///        in the ecumenical PolyNumeric data type, but subsequent functions will translate it
///        directly into other numeric or character formats.
///
/// \param item_name  Name of the data item
/// \param molecule   The MDL MOL entry with appended data items
/// \param fmt        Format of the numerical (or char4 tuple) data expected
/// \param policy     What to do if problems are encountered
std::vector<PolyNumeric> valuesFromSdfDataItem(const std::string &item_name,
                                               const MdlMol &molecule, const NumberFormat fmt,
                                               ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Return a series of real numbers found in an SD file data item's content.  This wraps the
///        more general valuesFromSdfDataItem() above, and descriptions of its formal arguments
///        follow.
std::vector<double> realFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                                        ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Return a series of integers found in an SD file data item's content.  This wraps the
///        more general valuesFromSdfDataItem() above, and descriptions of its formal arguments
///        follow.
std::vector<int> intFromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                                    ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Return a series of char4 tuples found in an SD file data item's content.  This wraps the
///        more general valuesFromSdfDataItem() above, and descriptions of its formal arguments
///        follow.
std::vector<char4> char4FromSdfDataItem(const std::string &item_name, const MdlMol &molecule,
                                        ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace structure
} // namespace stormm

#endif
