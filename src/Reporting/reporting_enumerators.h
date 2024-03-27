// -*-c++-*-
#ifndef STORMM_REPORTING_ENUMERATORS_H
#define STORMM_REPORTING_ENUMERATORS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"

namespace stormm {
namespace display {

/// \brief Possible modes in which command line program invocations can call for help messages
enum class HelpSignalKind {
  NO_ARGS,       ///< A command line call with no arguments will trigger a help message display
  NO_ARGS_ONLY,  ///< Only a call with no arguments, not a keyword, can trigger a help message
  KEYWORD,       ///< A keyword, by default any case-insensitive version of "help", "-help", or
                 ///<   or "--help" will trigger a help message display
  KEYWORD_ONLY   ///< Help messages are triggered only by an explicit keyword on the command line
};

/// \brief Produce strings detailing each of the enumerations above.
///
/// \param input  The enumeration to describe
std::string getEnumerationName(HelpSignalKind input);

} // namespace display

namespace review {

using constants::ExceptionResponse;

/// \brief Options for reporting results from one or more systems
enum class OutputScope {
  AVERAGES,          ///< Report the average energy of all systems, with standard deviations to
                     ///<   impart error bars.
  OUTLIERS,          ///< Report the average energy of all systems, with error bars, and individual
                     ///<   results for up to three outlier systems
  CLUSTER_AVERAGES,  ///< Report the average energy and error bars of systems governed by each
                     ///<   unique topology
  CLUSTER_OUTLIERS,  ///< Report the average energy and error bars of systems governed by each
                     ///<   unique topology, plus individual results for at most one outlier under
                     ///<   each topology
  FULL               ///< Report complete energies for all systems
};
  
/// \brief Possible formats for organized numerical output
enum class OutputSyntax {
  MATPLOTLIB,  ///< Results will be written to a script which the Python interpreter can run to
               ///<   display graphical results.  The script contains functions to simplify
               ///<   interactive plotting and rendering.
  MATRIX_PKG,  ///< A versatile format that works well with a number of matrix algebra packages,
               ///<   interspersing the data with commented blocks which serve the same purpose as
               ///<   narration in a typical output file.
  STANDALONE   ///< Results will be written to a unique format more typical of MD program output,
               ///<   with no syntax or other commands associated with plotting programs
};

/// \brief List the critical points in the molecular dynamics time step.
enum class IntegrationStage {
  BEGIN,              ///< Point reached prior to non-bonded force calculations (this is the
                      ///<   beginning of the step, but for practical purposes in a repeating MD
                      ///<   cycle this will conincide with the SHAKE or COORDINATE_UPDATE time
                      ///<   points depending on the use of constraints)
  NONBONDED_FORCE,    ///< Point reached after the non-bonded interactions have been computed
  ADDED_FORCE,        ///< Point reached after additional forces have been computed (the forces
                      ///<   on all atoms are expected to be accumulated in the order
                      ///<   non-bonded >> added >> bonded)
  BONDED_FORCE,       ///< Point reached after the bonded interactions have been computed
  VELOCITY_UPDATE,    ///< Point reached after adjusting velocities by half of the accumulated
                      ///<   forces with a Verlet integrator
  RATTLE,             ///< Point reached after applying Velocity constraints
  COORDINATE_UPDATE,  ///< Point reached after the coordinate and second-half velocity update
  SHAKE               ///< Point reached after geometric constraints are applied to all atoms
                      ///<   of the system
};

/// \brief Possible formats for grid file output.  Descriptions of some enumerations follow from
///        the OutputSyntax enumerator above.
enum class GridFileSyntax {
  MATPLOTLIB,  ///< Python-based MatPlotLib output, with isosurfaces
  MATRIX_PKG,  ///< Matlab and GNU Octave, with isosurfaces
  OPEN_DX,     ///< OpenDX format as used by VMD
  CUBEGEN      ///< Gaussian cubegen format
};
  
/// \brief Possible components of each output report file section
enum class SectionComponent {
  NARRATIVE,  ///< A block of narration suitable for processing by the
  LIST,       ///< A bulleted or ordered list within an OrderedList object (see ordered_list.h)
  TABLE,      ///< A formatted table, referencing the ReportTable object (see report_table.h)
  SCRIPT      ///< Additional script to print to the output file with no comment protection
};

/// \brief Possible methods of numbering or identifying separate items in an ordered list
enum class ListEnumeration {
  BULLET,      ///< Use a bulleted list, with every item taking the same single-character symbol
  NUMBERED,    ///< Number each item and place a ")" after the marker
  ALPHABETIC,  ///< Assign a (lowercase) letter to each item and place a ")" after the marker.  If
               ///<   the number of items rises above 26, the next items will have symbol aa, ab,
               ///<   (...).  This implies that there are 26 letters in the fastest incrementing
               ///<   position and 27 (add blank space) in the second, but the way to calculate the
               ///<   number of available combinations is as 26 + 26^2 + 26^3 + ...
  ROMAN,       ///< Assing a (lowercase) Roman numeral to each list item.
  NONE         ///< No enumeration symbol is used.
};

/// \brief Options for terminating a segment of text, such as a pargraph.  In some cases, a
///        new line is desirable, but in others the formatted text should just stop on its last
///        character.
enum class TextEnds {
  AS_IS,   ///< Stop printing with the last submitted character.
  NEWLINE  ///< Ensure that any line containing characters other than white space and a protective
           ///<   marker ends with a carriage return.
};

/// \brief Differtiate the types of ReportTable content with a specific enumeration.
enum class TableContentKind {
  INTEGER,  ///< The table contains data that can be read as integers (does not distinguish 32-bit
            ///<   from 64-bit integers)
  REAL,     ///< The table contains data that can be read as real numbers
  STRING    ///< The table data can only be interpreted as strings, possibly char4
};

/// \brief List the ways in which a surface can be depicted.
enum class SurfaceRender {
  WIRE,     ///< The surface is a wire mesh network with edges colored but no opacity on its faces
  SOLID,    ///< The surface is a collection of opaque or semi-transparent faces
  SCAFFOLD  ///< The surface has both opaque or (semi-transparent) faces and colored edges
};

/// \brief List a series of colors.
enum class PsivantColor {
  BLACK,         ///< Black,  RGB weights [ 0.120, 0.120, 0.120 ], hexcode #111c24
  RED,           ///< Red,    RGB weights [ 0.878, 0.000, 0.204 ], hexcode #e00034
  YELLOW,        ///< Blue,   RGB weights [ 0.000, 0.478, 0.788 ], hexcode #fdc82f
  BLUE,          ///< Yellow, RGB weights [ 0.992, 0.784, 0.184 ], hexcode #007ac9
  PURPLE,        ///< Purple, RGB weights [ 0.557, 0.145, 0.553 ], hexcode #8e258d
  GREY,          ///< Grey,   RGB weights [ 0.651, 0.651, 0.651 ], hexcode #a6a6a6
  LIGHT_RED,     ///< Light red,    RGB weights [ 1.000, 0.553, 0.655 ], hexcode #ff8da7
  LIGHT_YELLOW,  ///< Light yellow, RGB weights [ 0.996, 0.914, 0.675 ], hexcode #fee9ac
  LIGHT_BLUE,    ///< Light blue,   RGB weights [ 0.514, 0.808, 1.000 ], hexcode #83ceff
  LIGHT_PURPLE,  ///< Light purple, RGB weights [ 0.894, 0.588, 0.890 ], hexcode #e496e3
  WHITE          ///< White,        RGB weights [ 1.000, 1.000, 1.000 ], hexcode #ffffff
};

/// \brief Styles of plotting lines, valid in a range of visualization software.
enum class LinePlotStyle {
  SOLID,  ///< The line is solid with whatever weight
  DASHED  ///< The line is dashed with a dash length defined by the plotting program
};
  
/// \brief An array of Roman numerals, queried when using the ListEnumeration ROMAN setting.
const std::vector<std::string> roman_numerals = { "i", "ii", "iii", "iv", "v", "vi", "vii", "viii",
                                                  "ix", "x", "xi", "xii", "xiii", "xiv", "xv",
                                                  "xvi", "xvii", "xviii", "xix", "xx", "xxi",
                                                  "xxii", "xxiii", "xxiv", "xv", " " };

/// \brief The maximum value that STORMM will convert to Roman numerals.  Any integer greater than
///        or equal to this will produce an error.
constexpr int maximum_roman_numeral = 25;

/// \brief Produce strings detailing each of the enumerations above.
///
/// \param input  The enumeration to describe
/// \{
std::string getEnumerationName(OutputScope input);
std::string getEnumerationName(OutputSyntax input);
std::string getEnumerationName(IntegrationStage input);
std::string getEnumerationName(GridFileSyntax input);
std::string getEnumerationName(SectionComponent input);
std::string getEnumerationName(ListEnumeration input);
std::string getEnumerationName(TextEnds input);
std::string getEnumerationName(TableContentKind input);
std::string getEnumerationName(SurfaceRender input);
std::string getEnumerationName(LinePlotStyle input);
std::string getEnumerationName(PsivantColor input);
/// \}

/// \brief Translate a human-readable string indicating the output scope into the appropriate
///        enumeration.
///
/// \param input  The output scope to translate
OutputScope translateOutputScope(const std::string &input);
  
/// \brief Translate various directives for a surface rendering method into their corresponding
///        enumerations.
///
/// \param input  The human-readable string to translate
SurfaceRender translateSurfaceRender(const std::string &input);

/// \brief Traanslate various directives for the line plotting style into their corresponding
///        enumerations.
///
/// \param input  The human-readable string to translate
LinePlotStyle translateLinePlotStyle(const std::string &input);

/// \brief Perform a bounds check on an integer to see if it can be represented with one of the
///        hard-wired Roman numerals in the roman_numerals array above, then convert it if
///        possible.  Produces an error if the integer cannot be converted, or returns a single
///        whitespace character (the final entry of the array) with any degree of fault
///        tolerance.
///
/// \param x       The integer to validate
/// \param policy  The degree of fault tolerance
const std::string& toRoman(int x, ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Write a Psivant color for plotting in one of the supported output formats.
///
/// \param color   The color to display
/// \param syntax  The plotting format to display the color in
std::string encodePsivantColor(PsivantColor color, OutputSyntax syntax);
  
} // namespace review
} // namespace stormm

#endif
