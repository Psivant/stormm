// -*-c++-*-
#ifndef STORMM_PRESENT_FIELD_H
#define STORMM_PRESENT_FIELD_H

#include <fstream>
#include <iostream>
#include <limits.h>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_enumerators.h"
#include "FileManagement/file_listing.h"
#include "FileManagement/file_util.h"
#include "Math/math_enumerators.h"
#include "Math/rounding.h"
#include "Math/tickcounter.h"
#include "Parsing/polynumeric.h"
#include "Structure/background_mesh.h"
#include "Structure/mesh_mechanics.h"
#include "Structure/mesh_parameters.h"
#include "Structure/structure_enumerators.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "reporting_enumerators.h"
#include "render_molecule.h"
#include "render_options.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using constants::ExceptionResponse;
using diskutil::openOutputFile;
using diskutil::PrintSituation;
using diskutil::splitPath;
using parse::findAlignmentWidth;
using parse::intToString;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::realToString;
using stmath::LimitApproach;
using stmath::nearestFactor;
using stmath::primeFactors;
using stmath::TickCounter;
using structure::BackgroundMesh;
using structure::BackgroundMeshReader;
using structure::BoundaryCondition;
using structure::getEnumerationName;
using structure::interpolate;
using structure::MeshFFKit;
using structure::MeshFoundation;
using structure::MeshParameters;
using structure::MeshParamKit;
using structure::NonbondedPotential;
using structure::OffMeshProtocol;
using symbols::angstrom_to_bohr;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using trajectory::coordCopy;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;

/// \brief Default values for displaying a molecular scene, including a potential field
/// \{
constexpr int default_scene_decimal_places = 5;
constexpr char default_border_path_var[] = "border_path";
/// \}

/// \brief Get the scene's data file name based on the name of the script that establishes the
///        scene in one of the supported visualization packages.
///
/// \param scene_file  Name of the script that generates the scene
std::string getSceneDataFileName(const std::string &scene_file);

/// \brief Draw the borders of the mesh on the scene, tracing all edges as well as lines crossing
///        each face between either pair of opposite corners.  Each line is traced once and only
///        once.  This function presumes that the output file already has the mesh element matrix
///        defined as 'invU' and the mesh dimensions defined as 'na', 'nb', and 'nc'.
///
/// Overloaded:
///   - Output to a string
///   - Output directly to an open file stream
///
/// \param foutp      Output file stream to receive the code for border drawing
/// \param mps        Mesh parameters and dimensions for the relevant field
/// \param ropt       Rendering instructions, including the weight, color, and style of the field
///                   border lines
/// \param syntax     Syntax of the output file to write (some formats will be incompatible with
///                   this broder drawing)
/// \param path_name  Name of the variable that will hold the path (default 'border_path')
/// \{
std::string drawFieldBorders(const MeshParameters &mps, const RenderOptions &ropt,
                             GridFileSyntax syntax,
                             const std::string &path_name = std::string(default_border_path_var)); 

void drawFieldBorders(std::ofstream *foutp, const MeshParameters &mps, const RenderOptions &ropt,
                      GridFileSyntax syntax,
                      const std::string &path_name = std::string(default_border_path_var)); 
/// \}

/// \brief Print a BackgroundMesh object to disk using the specified file format.  Isosurfaces in
///        the potential field will be plotted if supported.
///
/// Overloaded:
///   - Provide a secondary set of mesh parameters to print a subset or even a superset of the mesh
///   - Use the mesh's innate parameters and print its entire contents
///
/// \param bgm             The background mesh object to print
/// \param mps             Auxiliary mesh parameters stipulating dimensions of the output.  Points
///                        in the output mesh lattice will be interpolated based on bgm.  If this
///                        is not provided, the parameters from bgm will be used by default.
/// \param file_name       Name of the file to print
/// \param expectation     The state that the named output file is expected to be found in
/// \param syntax          Syntax to use in printing.  This file type includes the Gaussian CubeGen
///                        and OpenDX formats.
/// \param ropt            Drawing details for the molecules and a list of isosurface values to
///                        plot from the associated mesh
/// \param decimal_places  The number of digits after the decimal to include in ASCII-formatted
///                        files
/// \param policy          Prescribed behavior in the event that the file cannot be written as
///                        requested
/// \{
template <typename T>
void printToFile(const BackgroundMesh<T> &bgm, const MeshParameters &mps,
                 const std::string &file_name, PrintSituation expectation, GridFileSyntax syntax,
                 const RenderOptions &ropt = RenderOptions(), int decimal_places = 5,
                 ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void printToFile(const BackgroundMesh<T> &bgm, const std::string &file_name,
                 PrintSituation expectation, GridFileSyntax syntax,
                 const RenderOptions &ropt = RenderOptions(), int decimal_places = 5,
                 ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

} // namespace review
} // namespace stormm

#include "present_field.tpp"

#endif
