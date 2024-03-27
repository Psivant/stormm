// -*-c++-*-
#ifndef STORMM_RENDER_OPTIONS_H
#define STORMM_RENDER_OPTIONS_H

#include <string>
#include <vector>
#include "copyright.h"
#ifndef STORMM_USE_HPC
#  include "DataTypes/stormm_vector_types.h"
#endif
#include "reporting_enumerators.h"

namespace stormm {
namespace review {

/// \brief Default values to be used in the RenderOptions object.  These are defined here, as
///        opposed to in the associated namelist header file for accepting user input to the
///        object.  The namelist header file will include this one, reversing the approach taken
///        for some other objects and their associated namelists.
/// \{
constexpr double default_rec_bond_thickness  = 1.5;
constexpr double default_lig_bond_thickness  = 2.5;
constexpr double default_border_thickness    = 1.5;
constexpr double default_rec_light_atom_size = 12.0;
constexpr double default_rec_heavy_atom_size = 20.0;
constexpr double default_lig_light_atom_size = 16.0;
constexpr double default_lig_heavy_atom_size = 24.0;
constexpr double default_wat_light_atom_size = 8.0;
constexpr double default_wat_heavy_atom_size = 16.0;
constexpr double default_added_point_size    = 8.0;
constexpr LinePlotStyle default_border_style = LinePlotStyle::SOLID;
constexpr uchar4 default_rec_bond_color      = {  48,  48,  48, 255 };
constexpr uchar4 default_lig_bond_color      = { 128, 224, 224, 255 };
constexpr uchar4 default_border_line_color   = { 128, 128, 128, 255 };
constexpr uchar4 default_added_point_color   = { 168,   0,  39, 255 };
constexpr char default_lighting_mode[] = "gouraud";
constexpr char default_camlight_position[] = "right";
/// \}

/// \brief Limits on various settings for displaying molecules and data in such context through a
///        matrix package such as Matlab or GNU Octave.
/// \{
constexpr int maximum_line_width = 16;
constexpr int maximum_marker_size = 64;
/// \}
  
//-------------------------------------------------------------------------------------------------
class RenderOptions {
public:

  /// \brief The default constructor will create an object with reasonable properties.  Various
  ///        setters can then modify particular elements of the scene that will be rendered.
  RenderOptions();

  /// \brief With no pointers to repair and no const member variables, the default copy and move
  ///        constructors as well as copy and move assignment operators are valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     An existing object on the right hand side of the assignment statement
  /// \{
  RenderOptions(const RenderOptions &original) = default;
  RenderOptions(RenderOptions &&original) = default;
  RenderOptions& operator=(const RenderOptions &original) = default;
  RenderOptions& operator=(RenderOptions &&original) = default;
  /// \}
  
  /// \brief Get the receptor bond thickness.  
  double getReceptorLineWidth() const;

  /// \brief Get the ligand bond thickness.
  double getLigandLineWidth() const;

  /// \brief Get the marker size to be used in depictions of light atoms in the receptor.
  double getReceptorLightAtomSize() const;

  /// \brief Get the marker size to be used in depictions of heavy atoms in the receptor.
  double getReceptorHeavyAtomSize() const;

  /// \brief Get the marker size to be used in depictions of light atoms in the ligand.
  double getLigandLightAtomSize() const;

  /// \brief Get the marker size to be used in depictions of heavy atoms in the ligand.
  double getLigandHeavyAtomSize() const;

  /// \brief Get the water hydrogen and virtual site size to display.
  double getWaterLightAtomSize() const;

  /// \brief Get the water oxygen size to display.
  double getWaterHeavyAtomSize() const;

  /// \brief Get the marker size to be used in depictions of added points.
  double getAddedPointSize() const;

  /// \brief Get the atom mask string for highlighting selected atoms of the receptor.
  const std::string& getReceptorHighlight() const;
  
  /// \brief Get the atom mask string for highlighting selected atoms of the ligand.
  const std::string& getLigandHighlight() const;
  
  /// \brief Get the line color to be used in depictions of the receptor.
  uchar4 getReceptorLineColor() const;

  /// \brief Get the color for depicting a line connecting atoms of the receptor in matrix package
  ///        syntax, as a string stating a three-element row vector of reals on a scale of [0.0,
  ///        1.0].  The string will be preceded by "'color'," (the word being in single quotes) and
  ///        followed by another comma, to present a compact way to express the color code as a
  ///        cohesive keyword / parameter pair for a plotting command.
  ///
  /// \param syntax  Indicate the visualization package which will ultimately process the color
  ///                details.
  std::string getFormattedReceptorLineColor(GridFileSyntax syntax) const;
  
  /// \brief Get the line color to be used in depictions of the ligand.
  uchar4 getLigandLineColor() const;

  /// \brief Get the color for depicting a line connecting atoms of the ligand in matrix package
  ///        syntax.  The methodology follows from getFormattedReceptorColor() above.
  ///
  /// \param syntax  Indicate the visualization package which will ultimately process the color
  ///                details.
  std::string getFormattedLigandLineColor(GridFileSyntax syntax) const;

  /// \brief Get the color to be used in placing arbitrary points in the scene.
  uchar4 getAddedPointColor() const;

  /// \brief Get the color for placing added points on the scene in matrix package syntax.  The
  ///        methodology follows from getFormattedReceptorColor() above.
  ///
  /// \param syntax  Indicate the visualization package which will ultimately process the color
  ///                details.
  std::string getFormattedAddedPointColor(GridFileSyntax syntax) const;

  /// \brief Get an indication of whether to display borders of the potential field in the scene.
  bool displayFieldBorders() const;

  /// \brief Get the line weight for plotting borders of a potential field in the scene.
  double getBorderLineWidth() const;

  /// \brief Get the line style for plotting borders of a potential field in the scene.
  LinePlotStyle getBorderLineStyle() const;

  /// \rbief
  
  /// \brief Get the line color for plotting borders of a potential field in the scene.
  uchar4 getBorderLineColor() const;

  /// \brief Get the color for drawing borders around the region containing the potential field,
  ///        formatted for a specific visualization package.  The methodology follows from
  ///        getFormattedReceptorColor() above.
  ///
  /// \param syntax  Indicate the visualization package which will ultimately process the color
  ///                details.
  std::string getFormattedBorderColor(GridFileSyntax syntax) const;
  
  /// \brief Get the lighting mode.
  const std::string& getLightingMode() const;

  /// \brief Get the camlight position.
  const std::string& getCamlightPosition() const;

  /// \brief Get the number of isosurfaces to render.
  int getIsosurfaceCount() const;
  
  /// \brief Get one of the requested isosurface values.  This will be plotted using the matrix
  ///        package's own isosurface() and patch() commands.
  ///
  /// \param index  Index of the isosurface described in one of the member variable arrays.  This
  ///               will be checked for validity.
  double getIsosurfaceValue(int index) const;

  /// \brief Get the color prescribed for depicting a requested isosurface.
  ///
  /// \param index  Index of the isosurface described in one of the member variable arrays
  uchar4 getIsosurfaceColor(int index) const;

  /// \brief Get the method of drawing one of the requested isosurfaces.
  ///
  /// \param index  Index of the isosurface described in one of the member variable arrays
  SurfaceRender getIsosurfaceDrawingMethod(int index) const;

  /// \brief Get the color for depicting a particular element.  The R, G, and B colors will be
  ///        given in the "x", "y", and "z" members of the tuple.  The alpha channel ("w" member)
  ///        is irrelevant, as markers, points, and lines cannot have transparency.
  ///
  /// \param z_number  The atomic number of the element in question (enter zero for virtual sites).
  ///                  This will be checked for validity.
  uchar4 getElementColor(int z_number) const;

  /// \brief Get the color for depicting a particular element in matrix package syntax.  The
  ///        methodology follows from getFormattedReceptorLineColor() above.
  ///
  /// \param z_number  The atomic number of the element in question (enter zero for virtual sites).
  ///                  This will be checked for validity.
  /// \param syntax  Indicate the visualization package which will ultimately process the color
  ///                details.
  std::string getFormattedElementColor(int z_number, GridFileSyntax syntax) const;

  /// \brief Get a const pointer to the object itself.
  const RenderOptions* getSelfPointer() const;
  
  /// \brief Set the line width to be used in plotting receptor bonds.
  ///
  /// \param width_in  The line weight to use
  void setReceptorLineWidth(double width_in);
  
  /// \brief Set the line width to be used in plotting ligand bonds.
  ///
  /// \param width_in  The line weight to use
  void setLigandLineWidth(double width_in);

  /// \brief Set the marker size to be used in plotting light atoms of the receptor.
  ///
  /// \param marker_size_in  The marker weight to use
  void setReceptorLightAtomSize(double marker_size_in);

  /// \brief Set the marker size to be used in plotting heavy atoms of the receptor.
  ///
  /// \param marker_size_in  The marker weight to use
  void setReceptorHeavyAtomSize(double marker_size_in);

  /// \brief Set the marker size to be used in plotting light atoms of the ligand.
  ///
  /// \param marker_size_in  The marker weight to use
  void setLigandLightAtomSize(double marker_size_in);

  /// \brief Set the marker size to be used in plotting heavy atoms of the ligand.
  ///
  /// \param marker_size_in  The marker weight to use
  void setLigandHeavyAtomSize(double marker_size_in);

  /// \brief Set the marker size to be used in plotting hydrogen atoms and virtual sites of water.
  ///
  /// \param marker_size_in  The marker weight to use
  void setWaterLightAtomSize(double marker_size_in);

  /// \brief Set the marker size to be used in plotting oxygen atoms of water.
  ///
  /// \param marker_size_in  The marker weight to use
  void setWaterHeavyAtomSize(double marker_size_in);

  /// \brief Set the marker size to be used in adding arbitrary points to the scene.
  ///
  /// \param marker_size_in  The marker weight to use
  void setAddedPointSize(double marker_size_in);

  /// \brief Set whether to display borders on the potential field.
  ///
  /// \param display_borders_in  Indicate whether to display the borders
  void setPotentialFieldBorderDisplay(bool display_borders_in);

  /// \brief Set the line width for border lines on the potential field.
  ///
  /// \param width_in  The border line width to set
  void setFieldBorderLineWidth(double width_in);

  /// \brief Set the style for displaying border lines on the potential field in the scene.
  ///
  /// \param style_in  The style (e.g. 'solid', 'dashed') to set
  void setFieldBorderStyle(LinePlotStyle style_in);

  /// \brief Set the color for displaying borders to the potential field in the scene.
  ///
  /// \param color_in  The color to set (the alpha channel is irrelevant for a line-based object)
  void setFieldBorderColor(const uchar4 color_in);
  
  /// \brief Set the receptor atom selection mask string.
  ///
  /// \param highlight_in  String specifying an atom mask that will be applied to the receptor
  void setReceptorHighlight(const std::string &highlight_in);
  
  /// \brief Set the ligand atom selection mask string.
  ///
  /// \param highlight_in  String specifying an atom mask that will be applied to the ligand
  void setLigandHighlight(const std::string &highlight_in);

  /// \brief Set the line color to be used in coloring the receptor's bonds.
  ///
  /// Overloaded:
  ///   - Provide an unsigned four-character tuple (includes an alpha channel, which is ignored)
  ///   - Provide three double precision numbers for red, green, and blue intensities on a scale
  ///     of 0.0 (this will be translated to 0 in the unsigned 8-bit integer representation) to
  ///     1.0 (255 in the unsigned 8-bit integer representation).
  ///
  /// \param color_in  Tuple containing R/G/B values
  /// \param red_in    The intensity of red in the color, in the range [0.0, 1.0].
  /// \param green_in  The intensity of green in the color, in the range [0.0, 1.0].
  /// \param blue_in   The intensity of blue in the color, in the range [0.0, 1.0].
  /// \{
  void setReceptorLineColor(const uchar4 color_in);
  void setReceptorLineColor(double red_in, double green_in, double blue_in);
  /// \}

  /// \brief Set the line color to be used in coloring the ligand's bonds.  Overloading and
  ///        descriptions of input parameters follow from setReceptorLineColor() above.
  /// \{
  void setLigandLineColor(const uchar4 color_in);
  void setLigandLineColor(double red_in, double green_in, double blue_in);
  /// \}
  
  /// \brief Set the line color to be used in coloring poitns added to the scene.  Overloading and
  ///        descriptions of input parameters follow from setReceptorLineColor() above.
  /// \{
  void setAddedPointColor(const uchar4 color_in);
  void setAddedPointColor(double red_in, double green_in, double blue_in);
  /// \}

  /// \brief Set the color of a specific element.  Overloading and descriptions of input
  ///        parameters follow from setReceptorLineColor() above, in addition to:
  ///
  /// \param z_number  Atomic number of the atom of interest
  /// \{
  void setAtomicElementColor(int z_number, const uchar4 color_in);
  void setAtomicElementColor(int z_number, double red_in, double green_in, double blue_in);
  /// \}

  /// \brief Set the lighting mode.  The mode will be checked for validity.  Allowed values
  ///        include "none", "flat", and "gouraud".
  ///
  /// \param mode_in  The lighting mode to set
  void setLightingMode(const std::string &mode_in);

  /// \brief Set the camlight position.  Acceptable values include "left", "right", and
  ///        "headlight".
  ///
  /// \param position_in  The position of the light relative to the camera.
  void setCamlightPosition(const std::string &position_in);

  /// \brief Add an isosurface to the list of requests.
  ///
  /// Overloaded:
  ///   - Provide the color and transparency for the surface as a uint_8 four-tuple (RGB/alpha)
  ///   - Provide the color and transparency for the surface as a double four-tuple (RGB/alpha)
  ///
  /// \param value_in     The isosurface value
  /// \param color_in     Color and transparency of the surface to create
  /// \param material_in  The imaginary material from which to build the surface
  /// \{
  void addIsosurface(double value_in, const uchar4 color_in, SurfaceRender material_in);
  void addIsosurface(double value_in, const double4 color_in, SurfaceRender material_in);
  /// \}
  
private:
  double receptor_bond_thickness;     ///< Bond thickness to use in depicting the receptor (this
                                      ///<   will typically be less than or equal to that used
                                      ///<   to depict a single conformation of the ligand, but
                                      ///<   more than that preferred to depict many
                                      ///<   conformations of the ligand)
  double ligand_bond_thickness;       ///< Bond thickness to use in depicting the ligand (ligand
                                      ///<   bonds can also have a different color, to which
                                      ///<   this feature adds distinction)
  double rec_light_atom_marker_size;  ///< Marker size for depicting light atoms (Z <= 4) of the
                                      ///<   receptor
  double rec_heavy_atom_marker_size;  ///< Marker size for depicting heavy atoms (Z > 4) of the
                                      ///<   receptor
  double lig_light_atom_marker_size;  ///< Marker size for depicting light atoms (Z <= 4) of the
                                      ///<   ligand
  double lig_heavy_atom_marker_size;  ///< Marker size for depicting heavy atoms (Z > 4) of the
                                      ///<   ligand
  double wat_light_atom_marker_size;  ///< Marker size for depicting hydrogen or virtual sites of
                                      ///<   water
  double wat_heavy_atom_marker_size;  ///< Marker size for depicting oxygen atoms of water
  double added_point_marker_size;     ///< Marker size to use for added points be depicted based
                                      ///<   on the provided field

  // The following member variables control definition of the borders for a potential field plotted
  // against the scene.
  bool display_field_borders;       ///< Flag to indicate that borders of the potential field
                                    ///<   should be displayed
  double border_line_thickness;     ///< The weight of the border lines
  LinePlotStyle border_line_style;  ///< The style of the border lines
  uchar4 border_line_color;         ///< The color of the border lines
  
  // The following member variables control highlighting, placing circles around additional atom
  // series.  Each string is translated into an atom mask to produce the relevant list of atoms.
  std::string receptor_highlight;  /// Mask string for highlighting atoms in the receptor
  std::string ligand_highlight;    /// Mask string for highlighting atoms in the ligand

  // The following member variables control colors for atoms, bonds, and points.
  uchar4 receptor_line_color;          ///< Line color to use for the receptor bonds with R, G, and
                                       ///<   B color values (0...255) in the "x", "y", and "z"
                                       ///<   members.  The alpha value of a line is irrelevant
                                       ///<   (semi-transparent lines are not implementable).
  uchar4 ligand_line_color;            ///< Line color to use for the ligand bonds
  uchar4 added_point_color;            ///< Marker color to use in depicting added points
  std::vector<uchar4> element_colors;  ///< Colors to use for each element, with the R, G, and B
                                       ///<   values (0...255) in the "x", "y", and "z" members.
                                       ///<   The alpha that could go in the "w" member is ignored
                                       ///<   (semi-transparent markers are not implementable).
  std::string lighting_mode;           ///< The lighting mode to use in the scene
  std::string camlight_position;       ///< Position of the camera light

  // The following member variables control isosurface depictions.
  int isosurface_count;                             ///< The number of requestedisosurfaces to
  std::vector<double> isosurface_values;            ///< Values for each isosurface, in the units
                                                    ///<   of the underlying mesh
  std::vector<uchar4> isosurface_colors;            ///< Colors for each isosurface, with the R, G,
                                                    ///<   and B values (0...255) in the "x", "y",
                                                    ///<   and "z" members and the alpha value
                                                    ///<   (0...255) in the "w" member
  std::vector<SurfaceRender> isosurface_materials;  ///< Drawing methods for each isosurface

  /// \brief Check the validity of the index for a requested isosurface.
  ///
  /// \param index  The index of interest
  void checkIsosurfaceIndex(int index) const;

  /// \brief Check an element's atomic number for validity.
  ///
  /// \param z_number  The atomic number of interest
  /// \param caller    Name of the calling function
  void checkZNumber(int z_number, const char* caller) const;
  
  /// \brief Check the line weight to be used in molecular plotting.
  ///
  /// \param width_in  The proposed line width
  /// \param caller    Name of the calling function
  void checkLineWidth(double width_in, const char* caller = nullptr) const;

  /// \brief Check the marker size to be used in molecular plotting.
  ///
  /// \param marker_size_in  The proposed marker size
  /// \param caller          Name of the calling function
  void checkMarkerSize(double marker_size_in, const char* caller = nullptr) const;

  /// \brief Check the range of red, green, and blue color codes for validity.
  ///
  /// \param red_in    The intensity of red, valid on a scale of [0.0, 1.0]
  /// \param green_in  The intensity of green, valid on a scale of [0.0, 1.0]
  /// \param blue_in   The intensity of blue, valid on a scale of [0.0, 1.0]
  /// \param caller    Name of the calling function
  void checkColorRanges(double red_in, double green_in, double blue_in,
                        const char* caller = nullptr) const;
};

/// \brief Get default colors for all elements to be used in matrix package scenes of molecules.
std::vector<uchar4> getDefaultElementColors();

/// \brief Prepare a formatted color keyword / value pair for matrix package plotting commands.
///
/// \param color_in  The 32-bit color code to process containing R / G / B readings on a scale
///                  of 0 to 255
/// \param syntax    Syntax of the output file for which to prepare the color specification
std::string formatColor(const uchar4 color_in, GridFileSyntax syntax);

} // namespace review
} // namespace stormm

#endif
