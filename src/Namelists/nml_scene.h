// -*-c++-*-
#ifndef STORMM_NML_SCENE_H
#define STORMM_NML_SCENE_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "Reporting/render_options.h"
#include "Reporting/reporting_enumerators.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using parse::TextFile;
using parse::WrapTextSearch;
using review::RenderOptions;
using review::SurfaceRender;
using review::LinePlotStyle;

/// \brief Collect output directives relating to scenes that STORMM might render through a matrix
///        package, MatPlotLib, or PyMol.  While the exact contents of the scene to be rendered
///        will be determined by the application and perhaps a particular process within it, the
///        stylistic details will be contained in one of these namelists.  The namelist contents
///        will primarily go into a RenderOptions object, but one notable element that the namelist
///        contains is a label that can be used by the calling program to control which of multiple
///        RenderOptions objects the contents should be directed into.
class SceneControls {
public:
  
  /// \brief The constructor can prepare an objct with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &minimize namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &minimize namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  SceneControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                WrapTextSearch wrap = WrapTextSearch::NO);
  SceneControls(const TextFile &tf, int *start_line, bool *found_nml,
                ExceptionResponse policy_in = ExceptionResponse::DIE,
                WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  SceneControls(const SceneControls &original) = default;
  SceneControls(SceneControls &&original) = default;
  SceneControls& operator=(const SceneControls &original) = default;
  SceneControls& operator=(SceneControls &&original) = default;
  /// \}

  /// \brief Get the label indicating how this control block is to be connected to program output.
  const std::string& getLabel() const;

  /// \brief Export the relevant information into a RenderOptions object.  Rather than individual
  ///        getters and setters, this object will dump all of its contents at once, into a new
  ///        RenderOptions object if necessary.  The individual getters and setters (with critical
  ///        validation checks) in that object can then be repurposed to filtering the contents of
  ///        this object.
  RenderOptions exportData() const;

  /// \brief Set the label of this object.  While it is expected that the label will come through a
  ///        user input deck, this feature is provided for testing purposes and pairing with the
  ///        specific getter function.
  void setLabel(const std::string &label_in);
  
private:

  /// The label controls which scene these controls will be directed towards.
  std::string label;

  // The following member variables control the rendering weights and cartoon representations of
  // molecules and other objects in the scene.
  double receptor_line_width;       ///< Width of lines or other cartoon elements to be used in
                                    ///<   drawing the receptor
  double ligand_line_width;         ///< Width of lines or other cartoon elements to be used in
                                    ///<   drawing the ligands
  double receptor_light_atom_size;  ///< Sizes of light atom representations in the receptor
  double receptor_heavy_atom_size;  ///< Sizes of heavy atom representations in the receptor
  double ligand_light_atom_size;    ///< Sizes of light atom representations in the ligand(s)
  double ligand_heavy_atom_size;    ///< Sizes of heavy atom representations in the ligand(s)
  double added_point_weight;        ///< Weight of points added and represented in the scene
  bool display_field_borders;       ///< Indicate whether to display the boundaries of the
                                    ///<   potential field in the scene
  double border_line_width;         ///< Width or lines defining the borders of the potential field
  LinePlotStyle border_line_style;  ///< Style to apply to borders of the potential field
  std::string receptor_highlight;   ///< Atom selection mask for the receptor
  std::string ligand_highlight;     ///< Atom selection mask for one or more representations of the
                                    ///<   ligand

  // The following member variables control the coloration and lighting of the scene.
  uchar4 receptor_line_color;     ///< The color of lines or other cartoon objects used to draw the
                                  ///<   receptor in the scene
  uchar4 ligand_line_color;       ///< The color of lines or other cartoon objects used to draw one
                                  ///<   or more ligands in the scene
  uchar4 added_point_color;       ///< The color of points added to the scene
  uchar4 border_line_color;       ///< Color of the border lines for an added potential field
  std::string lighting_mode;      ///< Method of simulating light reflected from rendered surfaces
  std::string camlight_position;  ///< Position of the light relative to the camera.  There are
                                  ///<   ways in matrix packages to place the light in some
                                  ///<   absolute position, but it is more useful to let the user
                                  ///<   manipulate the scene orientation with a lighting position
                                  ///<   that is constant relative to the perspective.

  // The following member variables control isosurfaces to be plotted from some field associated
  // with the scene.
  int isosurface_count;                        ///< The number of isosurfaces
  std::vector<double> isosurface_values;       ///< Value of the potential for which to plot each
                                               ///<   isosurface
  std::vector<uchar4> isosurface_colors;       ///< The colors to be used in each isosurface, which
                                               ///<   will contain red / green / blue components
                                               ///<   and, depending on the material and as enabled
                                               ///<   by the plotting program, an alpha channel
  std::vector<SurfaceRender> isosurface_mats;  ///< Styles in which isosurfaces will be rendered

  /// Store a deep copy of the original namelist emulator as read from the input file.x
  NamelistEmulator nml_transcript;
};
  
/// \brief Pull color information denoted by common subkeys from various STRUCT-type keywords in
///        the &scene namelist.  Both an RGB / alpha code as well as a verbose color name will be
///        sought, with specific RGB values overriding hard-coded color settings and alpha channels
///        being available only through the codified color entry (although, e.g. "color purple
///        -alpha 0.5 would be valid).
///
/// \param t_nml         The namelist at hand, as read from a user input deck
/// \param target_key    Keyword of STRUCT type containing color subkeys
/// \param target_entry  Specific iteration of the keyword to search for the color subkeys
/// \param init_value    Initial value for the color--this can be a default that is then modified,
///                      one member at a time, if valid subkeys are detected
uchar4 extractColorInformation(const NamelistEmulator &t_nml, const std::string &target_key,
                               int target_entry = 0,
                               const uchar4 init_value = { 255, 255, 255, 255 });

/// \brief Produce a namelist for specifying content of a scene that will be rendered as part of
///        STORMM's reporting for a particular program.  This namelist gives users a single
///        interface for controlling the output through multiple visualization programs.
///
/// \param tf
/// \param start_line
/// \param found
/// \param policy
/// \param wrap        
NamelistEmulator sceneInput(const TextFile &tf, int *start_line, bool *found,
                            ExceptionResponse policy = ExceptionResponse::DIE,
                            WrapTextSearch wrap = WrapTextSearch::NO);
  
} // namespace namelist
} // namespace stormm

#endif
