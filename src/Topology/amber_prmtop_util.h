// -*-c++-*-
#ifndef STORMM_AMBER_PMRTOP_UTIL_H
#define STORMM_AMBER_PMRTOP_UTIL_H

#include "copyright.h"
#include "Parsing/citation.h"
#include "Parsing/textfile.h"
#include "Parsing/polynumeric.h"
#include "atomgraph_enumerators.h"

namespace stormm {
namespace topology {

using parse::Citation;
using parse::TextFile;
using parse::PolyNumeric;
using parse::NumberFormat;

/// \brief The charge scaling factor found in Amber topologies, square root of the infamous
///        332.0525 that Amber uses (rather than the more accepted 332.0636) to obtain energy in
///        kcal/mol-e from Coulomb's constant, atomic units of charge, and distances in Angstroms.
constexpr double amber_charge_scaling = 18.2223;
constexpr double inv_amber_charge_scaling = 1.0 / amber_charge_scaling;

/// \brief Function for parsing a format string out of an Amber-format topology.  The exact format
///        of a given section may differ from that in the published standard (see
///        http://ambermd.org/FileFormats.php), so the details must be parsed at runtime.
///
/// \param fmt_in     The detected format string
/// \param file_name  The name of the file being parsed (for error reporting purposes)
/// \param line_idx   Line on the file where the format string is found (for error reporting)
int4 parseLineFormat(const std::string &fmt_in, const std::string &file_name, const int line_idx);

/// \brief Function for reading to a named delimiter within an Amber topology file.
///
/// \param tf               The text file, converted into RAM
/// \param flag             Flag to seek in the topology
/// \param detected_format  Format of the subsequent section found by parsing (returned)
/// \param essential        Indicates whether the delimiter must be found
/// \param start_line
int scanToFlag(const TextFile &tf, const char* flag, std::vector<int4>* detected_format,
               TopologyRequirement needed = TopologyRequirement::ESSENTIAL, int start_line = 0);

/// \brief Read citations / references to various force fields that this topology indicates were
///        used in its construction.
///
/// \param tf      A text file scanned into RAM containing the entire topology
/// \param lstart  The line at which to start reading model references
std::vector<Citation> readForceFieldReferences(const TextFile &tf, const int lstart);

/// \brief Read column-formatted input of a specific type.  This expects an Amber prmtop-like file
///        format whereing comment lines begin with '%'.
///
/// \param tf              The contents of the original text file translated into RAM
/// \param start_line      Line of the file at which to begin searching
/// \param cform           Format of each entry (i.e. scientific notation)
/// \param count_per_line  The number of entries expected on each line, save for the final line
///                        is expected to have up to this many
/// \param width           The text width of each entry
/// \param required_count  The number of entries which must be found
/// \param possible_count  The number of entries that may be encountered.  If left as zero, this
///                        will default to the required number of entries.
std::vector<PolyNumeric> amberPrmtopData(const TextFile &tf, int start_line, NumberFormat cform,
                                         int count_per_line, int width, int required_count,
                                         int possible_count = -1);

/// \brief Wrapper for the generic amberPrmtopData that returns a converted int vector
///
/// \param tf              The contents of the original text file translated into RAM
/// \param start_line      Line of the file at which to begin searching
/// \param count_per_line  The number of entries expected on each line, save for the final line
///                        is expected to have up to this many
/// \param width           The text width of each entry
/// \param required_count  The number of entries which must be found
/// \param possible_count  The number of entries that may be encountered.  If left as zero, this
///                        will default to the required number of entries.
std::vector<int> iAmberPrmtopData(const TextFile &tf, int start_line, int count_per_line,
                                  int width, int required_count, int possible_count = -1);

/// \brief Wrapper for the generic amberPrmtopData that returns a converted double vector
///
/// \param tf              The contents of the original text file translated into RAM
/// \param start_line      Line of the file at which to begin searching
/// \param count_per_line  The number of entries expected on each line, save for the final line
///                        is expected to have up to this many
/// \param width           The text width of each entry
/// \param required_count  The number of entries which must be found
/// \param possible_count  The number of entries that may be encountered.  If left as zero, this
///                        will default to the required number of entries.
std::vector<double> dAmberPrmtopData(const TextFile &tf, int start_line, int count_per_line,
                                     int width, int required_count, int possible_count = -1);

/// \brief Wrapper for the generic amberPrmtopData that returns a converted double vector
///
/// \param tf              The contents of the original text file translated into RAM
/// \param start_line      Line of the file at which to begin searching
/// \param count_per_line  The number of entries expected on each line, save for the final line
///                        is expected to have up to this many
/// \param width           The text width of each entry
/// \param required_count  The number of entries which must be found
/// \param possible_count  The number of entries that may be encountered.  If left as zero, this
///                        will default to the required number of entries.
std::vector<double> eAmberPrmtopData(const TextFile &tf, int start_line, int count_per_line,
                                     int width, int required_count, int possible_count = -1);

/// \brief Wrapper for the generic amberPrmtopData that returns a converted char4 vector
///
/// \param tf              The contents of the original text file translated into RAM
/// \param start_line      Line of the file at which to begin searching
/// \param count_per_line  The number of entries expected on each line, save for the final line
///                        is expected to have up to this many
/// \param width           The text width of each entry
/// \param required_count  The number of entries which must be found
/// \param possible_count  The number of entries that may be encountered.  If left as zero, this
///                        will default to the required number of entries.
std::vector<char4> c4AmberPrmtopData(const TextFile &tf, int start_line, int count_per_line,
                                     int required_count, int possible_count = -1);

} // namespace topology
} // namespace stormm

#endif
