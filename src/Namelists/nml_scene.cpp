#include <cmath>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "namelist_enumerators.h"
#include "nml_scene.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using parse::NumberFormat;
using parse::strcmpCased;
using parse::verifyContents;
using review::default_rec_bond_thickness;
using review::default_lig_bond_thickness;
using review::default_border_thickness;
using review::default_border_style;
using review::default_rec_light_atom_size;
using review::default_rec_heavy_atom_size;
using review::default_lig_light_atom_size;
using review::default_lig_heavy_atom_size;
using review::default_added_point_size;
using review::default_rec_bond_color;
using review::default_lig_bond_color;
using review::default_added_point_color;
using review::default_border_line_color;
using review::default_lighting_mode;
using review::default_camlight_position;
using review::translateSurfaceRender;
using review::translateLinePlotStyle;
  
//-------------------------------------------------------------------------------------------------
SceneControls::SceneControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    label{""},
    receptor_line_width{default_rec_bond_thickness},
    ligand_line_width{default_lig_bond_thickness},
    receptor_light_atom_size{default_rec_light_atom_size},
    receptor_heavy_atom_size{default_rec_heavy_atom_size},
    ligand_light_atom_size{default_lig_light_atom_size},
    ligand_heavy_atom_size{default_lig_heavy_atom_size},
    added_point_weight{default_added_point_size},
    display_field_borders{false},
    border_line_width{default_border_thickness},
    border_line_style{default_border_style},
    border_line_color{default_border_line_color},
    receptor_highlight{""},
    ligand_highlight{""},
    receptor_line_color{default_rec_bond_color},
    ligand_line_color{default_lig_bond_color},
    added_point_color{default_added_point_color},
    lighting_mode{std::string(default_lighting_mode)},
    camlight_position{std::string(default_camlight_position)},
    isosurface_count{0}, isosurface_values{}, isosurface_colors{}, isosurface_mats{},
    nml_transcript{"scene"}
{}

//-------------------------------------------------------------------------------------------------
SceneControls::SceneControls(const TextFile &tf, int *start_line, bool *found_nml,
                             const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    SceneControls(policy_in, wrap)
{
  NamelistEmulator t_nml = sceneInput(tf, start_line, found_nml, policy_in, wrap);
  nml_transcript = t_nml;

  // Hard-code a default color order for scene isosurfaces.  These values can be altered one color
  // component at a time.
  const std::vector<uchar4> color_order = { {   0, 114, 189, 255 },
                                            { 217,  83,  25, 255 },
                                            { 237, 177,  32, 255 },
                                            { 126,  47, 142, 255 },
                                            { 119, 172,  48, 255 },
                                            {  77, 190, 238, 255 },
                                            { 162,  20,  47, 255 } };
  
  // Read structure representation details
  if (t_nml.getKeywordStatus("label") != InputStatus::MISSING) {
    label = t_nml.getStringValue("label");
  }

  // While there are general keywords (e.g. "line_width", "light_atom_size"), they take effect at
  // the time the namelist is read and can be superceded by the corresponding specific keywords.
  // Read pieces of information from values associated with the specific keyword.
  receptor_line_width = t_nml.getRealValue("receptor_line_width");
  ligand_line_width = t_nml.getRealValue("ligand_line_width");
  receptor_light_atom_size = t_nml.getRealValue("receptor_light_atom_size");
  ligand_light_atom_size = t_nml.getRealValue("ligand_light_atom_size");
  receptor_heavy_atom_size = t_nml.getRealValue("receptor_heavy_atom_size");
  ligand_heavy_atom_size = t_nml.getRealValue("ligand_heavy_atom_size");
  added_point_weight = t_nml.getRealValue("added_point_weight");

  // Read atom selection masks
  if (t_nml.getKeywordStatus("receptor_highlight") != InputStatus::MISSING) {
    receptor_highlight = t_nml.getStringValue("receptor_highlight");
  }
  if (t_nml.getKeywordStatus("ligand_highlight") != InputStatus::MISSING) {
    ligand_highlight = t_nml.getStringValue("ligand_highlight");
  }

  // Read various color commands
  receptor_line_color = extractColorInformation(t_nml, "receptor_color", 0, receptor_line_color);  
  ligand_line_color = extractColorInformation(t_nml, "ligand_color", 0, ligand_line_color);  
  added_point_color = extractColorInformation(t_nml, "added_point_color", 0, added_point_color);  
  if (t_nml.getKeywordStatus("lighting_mode") != InputStatus::MISSING) {
    lighting_mode = t_nml.getStringValue("lighting_mode");
  }
  if (t_nml.getKeywordStatus("camlight") != InputStatus::MISSING) {
    camlight_position = t_nml.getStringValue("camlight");
  }

  // Read details of the potential field borders
  if (t_nml.getKeywordStatus("field_borders") != InputStatus::MISSING) {
    display_field_borders = true;
    border_line_width = t_nml.getRealValue("field_borders", "weight");
    border_line_color = extractColorInformation(t_nml, "field_borders", 0, border_line_color);
    border_line_style = translateLinePlotStyle(t_nml.getStringValue("field_borders", "style"));
  }
  else {
    display_field_borders = false;
  }
  
  // Read isosurface commands
  const int n_isosurf = t_nml.getKeywordEntries("isosurface");
  for (int i = 0; i < n_isosurf; i++) {
    isosurface_values.push_back(t_nml.getRealValue("isosurface", "value", i));
    const uchar4 init_color = (i < color_order.size()) ? color_order[i] :
                                                         uchar4({ 255, 255, 255, 255 });
    isosurface_colors.push_back(extractColorInformation(t_nml, "isosurface", i, init_color));
    const std::string i_mat = t_nml.getStringValue("isosurface", "kind", i);
    isosurface_mats.push_back(translateSurfaceRender(i_mat));
  }

  // Immediately attempt to export this data to a RenderOptions object.  This will engage input
  // validation immediately upon reading the input file, rather than waiting for what could be a
  // considerable amount of program wall time before the downstream effects of bad input become
  // apparent.  Follow it with a test that is guaranteed to pass to suppress compiler warnings.
  const RenderOptions mock_output = exportData();
  if (mock_output.getReceptorLightAtomSize() < 0) {
    rtErr("Atoms of the receptor must have a nonzero size.", "SceneControls");
  }
}

//-------------------------------------------------------------------------------------------------
RenderOptions SceneControls::exportData() const {
  RenderOptions result;
  result.setReceptorLineWidth(receptor_line_width);
  result.setLigandLineWidth(ligand_line_width);
  result.setReceptorLightAtomSize(receptor_light_atom_size);
  result.setReceptorHeavyAtomSize(receptor_heavy_atom_size);
  result.setLigandLightAtomSize(ligand_light_atom_size);
  result.setLigandHeavyAtomSize(ligand_heavy_atom_size);
  result.setAddedPointSize(added_point_weight);
  result.setReceptorHighlight(receptor_highlight);
  result.setLigandHighlight(ligand_highlight);
  result.setReceptorLineColor(receptor_line_color);
  result.setLigandLineColor(ligand_line_color);
  result.setAddedPointColor(added_point_color);
  result.setLightingMode(lighting_mode);
  result.setCamlightPosition(camlight_position);
  result.setPotentialFieldBorderDisplay(display_field_borders);
  result.setFieldBorderLineWidth(border_line_width);
  result.setFieldBorderStyle(border_line_style);
  result.setFieldBorderColor(border_line_color);
  for (int i = 0; i < isosurface_count; i++) {
    result.addIsosurface(isosurface_values[i], isosurface_colors[i], isosurface_mats[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
uchar4 extractColorInformation(const NamelistEmulator &t_nml, const std::string &target_key,
                               const int target_entry, const uchar4 init_value) {
  uchar4 result = init_value;
  if (t_nml.getKeywordStatus("target_key", "color", target_entry) != InputStatus::MISSING) {
    const std::string i_color = t_nml.getStringValue("isosurface", "color", target_entry);
    if (strcmpCased(i_color, "red", CaseSensitivity::NO)) {
      result = { 224,   0,  52, 255 };
    }
    else if (strcmpCased(i_color, "burgundy", CaseSensitivity::NO) ||
             strcmpCased(i_color, "dark_red", CaseSensitivity::NO) ||
             strcmpCased(i_color, "merlot", CaseSensitivity::NO)) {
      result = { 168,   0,  39, 255 };
    }
    else if (strcmpCased(i_color, "pink", CaseSensitivity::NO)) {
      result = { 255,  83, 123, 255 };
    }
    else if (strcmpCased(i_color, "blue", CaseSensitivity::NO) ||
             strcmpCased(i_color, "cerulean", CaseSensitivity::NO)) {
      result = {   0, 122, 201, 255 };
    }
    else if (strcmpCased(i_color, "light blue", CaseSensitivity::NO) ||
             strcmpCased(i_color, "light_blue", CaseSensitivity::NO) ||
             strcmpCased(i_color, "sky blue", CaseSensitivity::NO) ||
             strcmpCased(i_color, "sky_blue", CaseSensitivity::NO)) {
      result = {  70, 182, 255, 255 };
    }
    else if (strcmpCased(i_color, "navy", CaseSensitivity::NO) ||
             strcmpCased(i_color, "navy blue", CaseSensitivity::NO) ||
             strcmpCased(i_color, "navy_blue", CaseSensitivity::NO)) {
      result = {   0,  24, 128, 255 };
    }
    else if (strcmpCased(i_color, "green", CaseSensitivity::NO) ||
             strcmpCased(i_color, "forest", CaseSensitivity::NO)) {
      result = {   0, 192,   0, 255 };
    }
    else if (strcmpCased(i_color, "lime", CaseSensitivity::NO) ||
             strcmpCased(i_color, "lime green", CaseSensitivity::NO) ||
             strcmpCased(i_color, "lime_green", CaseSensitivity::NO) ||
             strcmpCased(i_color, "bright green", CaseSensitivity::NO) ||
             strcmpCased(i_color, "bright_green", CaseSensitivity::NO)) {
      result = {   0, 255,   0, 255 };
    }
    else if (strcmpCased(i_color, "gold", CaseSensitivity::NO) ||
             strcmpCased(i_color, "yellow", CaseSensitivity::NO)) {
      result = { 253, 200,  47, 255 };
    }
    else if (strcmpCased(i_color, "purple", CaseSensitivity::NO) ||
             strcmpCased(i_color, "plum", CaseSensitivity::NO)) {
      result = { 142,  37, 131, 255 };
    }
    else if (strcmpCased(i_color, "magenta", CaseSensitivity::NO)) {
      result = { 214,  97, 213, 255 };
    }
    else if (strcmpCased(i_color, "black", CaseSensitivity::NO)) {
      result = {   0,   0,   0, 255 };
    }
    else if (strcmpCased(i_color, "charcoal", CaseSensitivity::NO)) {
      result = {  48,  48,  48, 255 };
    }
    else if (strcmpCased(i_color, "grey", CaseSensitivity::NO) ||
             strcmpCased(i_color, "gray", CaseSensitivity::NO)) {
      result = { 128, 128, 128, 255 };
    }
    else if (strcmpCased(i_color, "white", CaseSensitivity::NO)) {
      result = { 255, 255, 255, 255 };        
    }
    else if (strcmpCased(i_color, "eggshell", CaseSensitivity::NO)) {
      result = { 240, 234, 214, 255 };
    }
    else if (strcmpCased(i_color, "cream", CaseSensitivity::NO)) {
      result = { 255, 253, 208, 255 };
    }
    else if (strcmpCased(i_color, "vanilla", CaseSensitivity::NO)) {
      result = { 243, 229, 171, 255 };
    }
    else if (strcmpCased(i_color, "brown", CaseSensitivity::NO) ||
             strcmpCased(i_color, "chocolate", CaseSensitivity::NO)) {
      result = { 150,  75,   0, 255 };
    }
    else if (strcmpCased(i_color, "tan", CaseSensitivity::NO) ||
             strcmpCased(i_color, "hazelnut", CaseSensitivity::NO)) {
      result = { 210, 180, 140, 255 };
    }
    else {
      rtErr("Unrecognized color " + i_color + ".", "extractColorInformation");
    }
  }
  std::vector<std::string> rgba = { "-r", "-g", "-b", "-alpha" };
  for (size_t i = 0; i < rgba.size(); i++) {
    try {
      if (t_nml.getKeywordStatus(target_key, rgba[i], target_entry) != InputStatus::MISSING) {
        const std::string clr_string = t_nml.getStringValue(target_key, rgba[i], target_entry);
        bool problem = false;
        int clr_val = 0;
        if (verifyContents(clr_string, NumberFormat::INTEGER)) {
          clr_val = stoi(clr_string);
          if (clr_val < 0 || clr_val > 255) {
            problem = true;
          }
        }
        else if (verifyContents(clr_string, NumberFormat::STANDARD_REAL) ||
                 verifyContents(clr_string, NumberFormat::SCIENTIFIC)) {
          const double dclr_val = stod(clr_string);
          if (dclr_val < 0.0 || dclr_val > 1.0) {
            problem = true;
          }
          clr_val = round(dclr_val * 255.0);
        }
        else {
          problem = true;
        }
        if (problem) {
          rtErr("Input for color element " + rgba[i] + " must be an integer in the range [0, 255] "
                "or a real number in the range [0.0, 1.0].", "extractColorInformation");
        }

        // Assign the value to the proper position in the color code          
        if (i == 0) {
          result.x = clr_val;
        }
        else if (i == 1) {
          result.y = clr_val;
        }
        else if (i == 2) {
          result.z = clr_val;
        }
        else if (i == 3) {
          result.w = clr_val;
        }
      }
    }
    catch (std::runtime_error) {

      // If the subkey is not present in the STRUCT-type keyword, proceed to the next element.
      // Not all color settings are compatible with an alpha channel, so rather than require all
      // inputs using this parser to allow for an an alpha channel, this exception handling
      // provides a way to skip missing input.
      continue;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator sceneInput(const TextFile &tf, int *start_line, bool *found,
                            const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("scene", CaseSensitivity::AUTOMATIC, policy, "Collects information "
                         "needed for rendering a scene.  Generally, a scene involves a receptor, "
                         "one or more ligand poses, and possibly a spatial field associated with "
                         "the receptor.  Each &scene namelist can be tagged with a label to "
                         "connect its details to one of several scenes that a program might "
                         "render at runtime.");

  // Keyword: label
  t_nml.addKeyword("label", NamelistType::STRING);
  t_nml.addHelp("label", "Attach a label to the scene controls to indicate one of a number of "
                "possible output scenes that they should apply to.  This specification is not "
                "necessary in applications which only involve a single scene output.");

  // Keyword: line_width
  t_nml.addKeyword("line_width", NamelistType::REAL);
  t_nml.addHelp("line_width", "Specify the width of lines drawn to connect atoms of both ligand "
                "and receptor molecules in the system.  This input will be overridden by "
                "specifying receptor_line_width or ligand_line_width.");

  // Keywords for ligand and receptor drawing
  const std::vector<std::string> molecules = { "receptor", "ligand", "water" };
  const std::vector<std::string> atom_weights = { "light", "heavy" };
  const std::vector<std::string> color_keys_help = {
    "Specify the color of interest by name.  This does not specify an alpha value.",
    "Specify the red component of the RGB color spectrum, either as an integer in the range [0, "
    "255] or as a real number in the range [0.0, 1.0].",
    "Specify the green component of the RGB color spectrum, either as an integer in the range [0, "
    "255] or as a real number in the range [0.0, 1.0].",
    "Specify the blue component of the RGB color spectrum, either as an integer in the range [0, "
    "255] or as a real number in the range [0.0, 1.0].",
    "Specify the alpha channel of the RGB color spectrum, either as an integer in the range [0, "
    "255] or as a real number in the range [0.0, 1.0]. (Zero is completely transparent.)",
  };
  const std::vector<NamelistType> color_keys_types(5, NamelistType::STRING);
  const std::vector<std::string> color_keys_defaults = { std::string(""), std::string("255"),
                                                         std::string("255"), std::string("255") };
  const std::vector<KeyRequirement> color_keys_reqs(5, KeyRequirement::OPTIONAL);
  t_nml.addKeyword("line_width", NamelistType::REAL);
  t_nml.addHelp("line_width", "Specify the weight of bonds between atoms of both the receptor and "
                "ligand.");
  for (size_t i = 0; i < atom_weights.size(); i++) {
    const std::string kw = atom_weights[i] + "_atom_size";
    t_nml.addKeyword(kw, NamelistType::REAL);
    t_nml.addHelp(kw, "Specify the sizing of all " + atom_weights[i] + " in the system.  This "
                  "will be overridden by specifying " + molecules[0] + ", " + molecules[1] +
                  " or " + molecules[2] + "_atom_size.");
  }
  for (size_t i = 0; i < molecules.size(); i++) {

    // Exclude water from certain keyword generation
    if (i < 2) {
      t_nml.addKeyword(molecules[i] + "_line_width", NamelistType::REAL);
      t_nml.addHelp(molecules[i] + "_line_width", "Specify the width of lines drawn to connect "
                    "atoms of the " + molecules[i] + ".");
      t_nml.addKeyword(molecules[i] + "_highlight", NamelistType::STRING);
      t_nml.addHelp(molecules[i] + "_highlight", "Specify an atom mask selecting parts of the " +
                    molecules[i] + " for emphasis.");
      t_nml.addKeyword(molecules[i] + "_color", { "color", "-r", "-g", "-b", "-alpha" },
                       color_keys_types, color_keys_defaults, DefaultIsObligatory::NO,
                       InputRepeats::NO, "Set the color of the " + molecules[i] + " structure.",
                       color_keys_help, color_keys_reqs);
    }
    for (size_t j = 0; j < atom_weights.size(); j++) {
      const std::string kw = molecules[i] + "_" + atom_weights[j] + "_atom_size";
      t_nml.addKeyword(kw, NamelistType::REAL);
      t_nml.addHelp(kw, "Specify the sizing of " + atom_weights[j] + " atoms in the " +
                    molecules[i] + ".");
    }
  }
  t_nml.addKeyword("added_point_size", NamelistType::REAL);
  t_nml.addHelp("added_point_size", "Specify the sizing of points in a scatter plot overlayed on "
                "the scene.");
  t_nml.addKeyword("added_point_color", { "color", "-r", "-g", "-b" },
                   color_keys_types, color_keys_defaults, DefaultIsObligatory::NO,
                   InputRepeats::NO, "Set the color of points added in a scatter plot overlayed "
                   "on the scene.", color_keys_help, color_keys_reqs);
  const std::vector<NamelistType> border_keys_types = { NamelistType::STRING, NamelistType::STRING,
                                                        NamelistType::STRING, NamelistType::STRING,
                                                        NamelistType::STRING, NamelistType::REAL,
                                                        NamelistType::STRING };
  const std::vector<std::string> border_keys_defaults = { "grey", "128", "128", "128", "1.5",
                                                          "solid" };
  const std::vector<std::string> border_keys_help = {
    color_keys_help[0] + "  Default 'grey'.", color_keys_help[1], color_keys_help[2],
    color_keys_help[3],
    "Specify the weight of border lines for the potential field in the scene.  Default 1.5.",
    "Specify the style of border lines for the potential field in the scene.  Options include "
    "'solid' and 'dashed.' Default 'solid'."
  };
  const std::vector<KeyRequirement> border_keys_reqs(6, KeyRequirement::OPTIONAL);
  t_nml.addKeyword("field_borders", { "color", "-r", "-g", "-b", "weight", "style" },
                   border_keys_types, border_keys_defaults, DefaultIsObligatory::NO,
                   InputRepeats::NO, "Set the color and line style of the borders of a potential "
                   "field rendered within the scene.  By default, this feature is OFF, even if "
                   "isosurfaces of some potential are plotted.  Specifying border color and style "
                   "will add this detail to the scene.", border_keys_help, border_keys_reqs);

  // Lighting controls
  t_nml.addKeyword("lighting_mode", NamelistType::STRING);
  t_nml.addHelp("Set the manner in which the scene will be lit, based on the way adjacent faces "
                "reflect incoming light.");
  t_nml.addKeyword("camlight", NamelistType::STRING);
  t_nml.addHelp("Set the position of the light illuminating the scene, relative to the camera "
                "perspective.");

  // Isosurface controls
  const std::vector<NamelistType> isosurf_keys_types(7, NamelistType::STRING);
  const std::vector<KeyRequirement> isosurf_keys_reqs = {
    KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
    KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
    KeyRequirement::OPTIONAL
  };
  const std::vector<std::string> isosurf_keys_defaults = { "", "white", "", "", "", "", "wire" };
  t_nml.addKeyword("isosurface", { "value", "color", "-r", "-g", "-b", "-alpha", "material" },
                   isosurf_keys_types, isosurf_keys_defaults, DefaultIsObligatory::NO,
                   InputRepeats::YES, "Specify an isosurface of a potential field associated with "
                   "the receptor to display.",
                   { "Value of the isosurface to display", color_keys_help[0], color_keys_help[1],
                     color_keys_help[2], color_keys_help[3], color_keys_help[4], "Material to use "
                     "in creating the isosurface, acceptable values including \"wire\", "
                     "\"solid\", and \"scaffold\"." }, isosurf_keys_reqs);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
