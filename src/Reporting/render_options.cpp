#include <cmath>
#include "copyright.h"
#include "Chemistry/periodic_table.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "error_format.h"
#include "render_options.h"

namespace stormm {
namespace review {

using chemistry::element_maximum_count;
using parse::NumberFormat;
using parse::realToString;
using parse::rgbHexCode;

//-------------------------------------------------------------------------------------------------
RenderOptions::RenderOptions() :
    receptor_bond_thickness{default_rec_bond_thickness},
    ligand_bond_thickness{default_lig_bond_thickness},
    rec_light_atom_marker_size{default_rec_light_atom_size},
    rec_heavy_atom_marker_size{default_rec_heavy_atom_size},
    lig_light_atom_marker_size{default_lig_light_atom_size},
    lig_heavy_atom_marker_size{default_lig_heavy_atom_size},
    wat_light_atom_marker_size{default_wat_light_atom_size},
    wat_heavy_atom_marker_size{default_wat_heavy_atom_size},
    added_point_marker_size{default_added_point_size},
    display_field_borders{false},
    border_line_thickness{default_border_thickness},
    border_line_style{default_border_style},
    border_line_color{default_border_line_color},
    receptor_line_color{default_rec_bond_color},
    ligand_line_color{default_lig_bond_color},
    added_point_color{default_added_point_color},
    element_colors{getDefaultElementColors()},
    lighting_mode{std::string(default_lighting_mode)},
    camlight_position{std::string(default_camlight_position)},
    isosurface_count{0}, isosurface_values{}, isosurface_colors{}, isosurface_materials{}
{}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getReceptorLineWidth() const {
  return receptor_bond_thickness;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getLigandLineWidth() const {
  return ligand_bond_thickness;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getReceptorLightAtomSize() const {
  return rec_light_atom_marker_size;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getReceptorHeavyAtomSize() const {
  return rec_heavy_atom_marker_size;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getLigandLightAtomSize() const {
  return lig_light_atom_marker_size;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getLigandHeavyAtomSize() const {
  return lig_heavy_atom_marker_size;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getWaterLightAtomSize() const {
  return wat_light_atom_marker_size;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getWaterHeavyAtomSize() const {
  return wat_heavy_atom_marker_size;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getAddedPointSize() const {
  return added_point_marker_size;
}

//-------------------------------------------------------------------------------------------------
const std::string& RenderOptions::getReceptorHighlight() const {
  return receptor_highlight;
}

//-------------------------------------------------------------------------------------------------
const std::string& RenderOptions::getLigandHighlight() const {
  return ligand_highlight;
}

//-------------------------------------------------------------------------------------------------
uchar4 RenderOptions::getReceptorLineColor() const {
  return receptor_line_color;
}

//-------------------------------------------------------------------------------------------------
std::string RenderOptions::getFormattedReceptorLineColor(const GridFileSyntax syntax) const {
  return formatColor(receptor_line_color, syntax);
}

//-------------------------------------------------------------------------------------------------
uchar4 RenderOptions::getLigandLineColor() const {
  return ligand_line_color;
}

//-------------------------------------------------------------------------------------------------
std::string RenderOptions::getFormattedLigandLineColor(const GridFileSyntax syntax) const {
  return formatColor(ligand_line_color, syntax);
}

//-------------------------------------------------------------------------------------------------
uchar4 RenderOptions::getAddedPointColor() const {
  return added_point_color;
}

//-------------------------------------------------------------------------------------------------
std::string RenderOptions::getFormattedAddedPointColor(const GridFileSyntax syntax) const {
  return formatColor(added_point_color, syntax);
}

//-------------------------------------------------------------------------------------------------
bool RenderOptions::displayFieldBorders() const {
  return display_field_borders;
}
  
//-------------------------------------------------------------------------------------------------
double RenderOptions::getBorderLineWidth() const {
  return border_line_thickness;
}

//-------------------------------------------------------------------------------------------------
LinePlotStyle RenderOptions::getBorderLineStyle() const {
  return border_line_style;
}

//-------------------------------------------------------------------------------------------------
uchar4 RenderOptions::getBorderLineColor() const {
  return border_line_color;
}

//-------------------------------------------------------------------------------------------------
std::string RenderOptions::getFormattedBorderColor(const GridFileSyntax syntax) const {
  return formatColor(border_line_color, syntax);
}

//-------------------------------------------------------------------------------------------------
int RenderOptions::getIsosurfaceCount() const {
  return isosurface_count;
}

//-------------------------------------------------------------------------------------------------
double RenderOptions::getIsosurfaceValue(const int index) const {
  checkIsosurfaceIndex(index);
  return isosurface_values[index];
}

//-------------------------------------------------------------------------------------------------
uchar4 RenderOptions::getIsosurfaceColor(const int index) const {
  checkIsosurfaceIndex(index);
  return isosurface_colors[index];
}

//-------------------------------------------------------------------------------------------------
SurfaceRender RenderOptions::getIsosurfaceDrawingMethod(const int index) const {
  checkIsosurfaceIndex(index);
  return isosurface_materials[index];
}

//-------------------------------------------------------------------------------------------------
uchar4 RenderOptions::getElementColor(const int z_number) const {
  checkZNumber(z_number, "getElementColor");
  return element_colors[z_number];
}

//-------------------------------------------------------------------------------------------------
std::string RenderOptions::getFormattedElementColor(const int z_number,
                                                    const GridFileSyntax syntax) const {
  checkZNumber(z_number, "getFormattedElementColor");
  return formatColor(element_colors[z_number], syntax);
}
  
//-------------------------------------------------------------------------------------------------
const RenderOptions* RenderOptions::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setReceptorLineWidth(const double width_in) {
  checkLineWidth(width_in, "setReceptorLineWidth");
  receptor_bond_thickness = width_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLigandLineWidth(const double width_in) {
  checkLineWidth(width_in, "setLigandLineWidth");
  ligand_bond_thickness = width_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setReceptorLightAtomSize(const double marker_size_in) {
  checkMarkerSize(marker_size_in, "setReceptorLightAtomSize");
  rec_light_atom_marker_size = marker_size_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setReceptorHeavyAtomSize(const double marker_size_in) {
  checkMarkerSize(marker_size_in, "setReceptorHeavyAtomSize");
  rec_heavy_atom_marker_size = marker_size_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLigandLightAtomSize(const double marker_size_in) {
  checkMarkerSize(marker_size_in, "setLigandLightAtomSize");
  lig_light_atom_marker_size = marker_size_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLigandHeavyAtomSize(const double marker_size_in) {
  checkMarkerSize(marker_size_in, "setLigandHeavyAtomSize");
  lig_heavy_atom_marker_size = marker_size_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setAddedPointSize(const double marker_size_in) {
  checkMarkerSize(marker_size_in, "setAddedPointSize");
  added_point_marker_size = marker_size_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setPotentialFieldBorderDisplay(const bool display_borders_in) {
  display_field_borders = display_borders_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setFieldBorderLineWidth(const double width_in) {
  border_line_thickness = width_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setFieldBorderStyle(const LinePlotStyle style_in) {
  border_line_style = style_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setFieldBorderColor(const uchar4 color_in) {
  border_line_color = color_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setReceptorHighlight(const std::string &highlight_in) {
  receptor_highlight = highlight_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLigandHighlight(const std::string &highlight_in) {
  ligand_highlight = highlight_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setReceptorLineColor(const uchar4 color_in) {
  receptor_line_color = color_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setReceptorLineColor(const double red_in, const double blue_in,
                                         const double green_in) {
  checkColorRanges(red_in, green_in, blue_in, "setReceptorLineColor");
  uchar4 icolor;
  icolor.x = round(red_in * 255.0);
  icolor.y = round(green_in * 255.0);
  icolor.z = round(blue_in * 255.0);
  icolor.w = 255;
  setReceptorLineColor(icolor);
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLigandLineColor(const uchar4 color_in) {
  ligand_line_color = color_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLigandLineColor(const double red_in, const double blue_in,
                                       const double green_in) {
  checkColorRanges(red_in, green_in, blue_in, "setLigandLineColor");
  uchar4 icolor;
  icolor.x = round(red_in * 255.0);
  icolor.y = round(green_in * 255.0);
  icolor.z = round(blue_in * 255.0);
  icolor.w = 255;
  setLigandLineColor(icolor);
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setAddedPointColor(const uchar4 color_in) {
  added_point_color = color_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setAddedPointColor(const double red_in, const double blue_in,
                                       const double green_in) {
  checkColorRanges(red_in, green_in, blue_in, "setAddedPointColor");
  uchar4 icolor;
  icolor.x = round(red_in * 255.0);
  icolor.y = round(green_in * 255.0);
  icolor.z = round(blue_in * 255.0);
  icolor.w = 255;
  setAddedPointColor(icolor);
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setAtomicElementColor(const int z_number, const uchar4 color_in) {
  checkZNumber(z_number, "setAtomicElementColor");  
  element_colors[z_number] = color_in;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setAtomicElementColor(const int z_number, const double red_in,
                                          const double blue_in, const double green_in) {
  checkColorRanges(red_in, green_in, blue_in, "setAtomicElementColor");
  uchar4 icolor;
  icolor.x = round(red_in * 255.0);
  icolor.y = round(green_in * 255.0);
  icolor.z = round(blue_in * 255.0);
  icolor.w = 255;
  setAtomicElementColor(z_number, icolor);
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setLightingMode(const std::string &mode_in) {
  if (mode_in == "gouraud" || mode_in == "flat" || mode_in == "none") {
    lighting_mode = mode_in;
  }
  else {
    const std::string phong_error = (mode_in == "phong") ?
      "Lighting mode \"phong\" is not supported for surfaces and patches." : "";
    rtErr("The lighting mode \"" + mode_in + "\" is invalid." + phong_error, "RenderOption",
          "setLightingMode");
  }
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::setCamlightPosition(const std::string &position_in) {
  if (position_in == "left" || position_in == "right" || position_in == "headlight") {
    camlight_position = position_in;
  }
  else {
    rtErr("Camlight position \"" + position_in + "\" is invalid.", "RenderOption",
          "setCamlightPosition");
  }
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::addIsosurface(const double value_in, const uchar4 color_in,
                                  const SurfaceRender material_in) {
  isosurface_values.push_back(value_in);
  isosurface_colors.push_back(color_in);
  isosurface_materials.push_back(material_in);
  isosurface_count += 1;
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::addIsosurface(const double value_in, const double4 color_in,
                                  const SurfaceRender material_in) {
  checkColorRanges(color_in.x, color_in.y, color_in.z, "addIsosurface");
  if (color_in.w < 0.0 || color_in.w > 1.0) {
    rtErr("And alpha intensity of " + realToString(color_in.w, 7, 4, NumberFormat::STANDARD_REAL) +
          " is invalid.  The acceptable range is [0.0, 1.0].", "RenderOptions", "addIsosurface"); 
  }
  uchar4 icolor;
  icolor.x = round(color_in.x * 255.0);
  icolor.y = round(color_in.y * 255.0);
  icolor.z = round(color_in.z * 255.0);
  icolor.w = round(color_in.w * 255.0);
  addIsosurface(value_in, icolor, material_in);
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::checkIsosurfaceIndex(const int index) const {

  // Isosurfaces enter the object as complete definitions, so there is a correspondence enforced
  // between the various array sizes and comparison to the one counter variable is valid.
  if (index >= isosurface_count) {
    rtErr("Isosurface index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(isosurface_count) + " surfaces.", "RenderOptions",
          "checkIsosurfaceIndex");
  }
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::checkZNumber(const int z_number, const char* caller) const {
  if (z_number < 0 || z_number >= element_maximum_count) {
    rtErr("Atomic number " + std::to_string(z_number) + " is invalid.", "RenderOptions",
          "getElementColor");
  }
}
  
//-------------------------------------------------------------------------------------------------
void RenderOptions::checkLineWidth(const double width_in, const char* caller) const {
  if (width_in < 1 || width_in > maximum_line_width) {
    rtErr("A line width of " + std::to_string(width_in) + " is unacceptable.", "RenderOptions",
          caller);
  }
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::checkMarkerSize(const double marker_size_in, const char* caller) const {
  if (marker_size_in < 0.1 || marker_size_in > maximum_marker_size) {
    rtErr("A marker weight of " + realToString(marker_size_in, 8, 4, NumberFormat::STANDARD_REAL) +
          " is unacceptable.", "RenderOptions", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void RenderOptions::checkColorRanges(const double red_in, const double green_in,
                                     const double blue_in, const char* caller) const {
  if (red_in < 0.0 || red_in > 1.0) {
    rtErr("The intensity of red (" + realToString(red_in, 7, 4, NumberFormat::STANDARD_REAL) +
          " must fall in the range [0.0, 1.0].\n", "RenderOptions", caller);
  }
  if (green_in < 0.0 || green_in > 1.0) {
    rtErr("The intensity of green (" + realToString(green_in, 7, 4, NumberFormat::STANDARD_REAL) +
          " must fall in the range [0.0, 1.0].\n", "RenderOptions", caller);
  }
  if (blue_in < 0.0 || green_in > 1.0) {
    rtErr("The intensity of blue (" + realToString(blue_in, 7, 4, NumberFormat::STANDARD_REAL) +
          " must fall in the range [0.0, 1.0].\n", "RenderOptions", caller);
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<uchar4> getDefaultElementColors() {
  std::vector<uchar4> result(element_maximum_count, { 128, 128, 128, 255 });

  // Use Psivant-inspired colors for various common elements.  Virtual sites get a special size as
  // well as a unique color.
  result[ 0] = { 142,  37, 131, 255 };  // [VS] Psivant Plum - accent #1

  // Hydrogen
  result[ 1] = { 229, 229, 229, 255 };  // [H ] Psivant White - background, darker 10%

  // Helium, Neon, and Argon  
  result[ 2] = { 131, 206, 255, 255 };  // [He] Psivant Blue - accent #3, lighter 60%
  result[10] = {  70, 182, 255, 255 };  // [Ne] Psivant Blue - accent #3, lighter 40%
  result[18] = { 193, 231, 255, 255 };  // [Ar] Psivant Blue - accent #3, lighter 80%

  // Lithium, Sodium, and Potassium
  result[ 3] = {  63, 113, 150, 255 };  // [Li] Psivant Ice Blue - accent #6, darker 60%
  result[11] = {  62, 113, 150, 255 };  // [Na] Psivant Ice Blue - accent #6, darker 50%
  result[19] = { 120, 165, 199, 255 };  // [K ] Psivant Ice Blue - accent #6, darker 25%

  // Boron  
  result[ 5] = { 255,  83, 123, 255 };  // [B ] Psivant Red - accent #1, lighter 40%

  // Carbon
  result[ 6] = {  17,  28,  36, 255 };  // [C ] Psivant Blue/Grey - main text

  // Magnesium and calcium
  result[12] = { 228, 150, 227, 255 };  // [Mg] Psivant Plum - accent #1, lighter 60%
  result[20] = { 241, 202, 241, 255 };  // [Ca] Psivant Plum - accent #1, lighter 80%

  // Nitrogen and phosphorus
  result[ 7] = {   0, 122, 201, 255 };  // [N ] Psivant Blue - accent #3
  result[15] = { 214,  97, 213, 255 };  // [P ] Psivant Plum - accent #1, lighter 40%

  // Oxygen, sulfur, and selenium
  result[ 8] = { 224,   0,  52, 255 };  // [O ] Psivant Red - accent #1
  result[16] = { 253, 200,  47, 255 };  // [S ] Psivant Gold - accent #2
  result[34] = { 223, 166,   2, 255 };  // [Se] Psivant Gold - accent #2, darker 25%
  
  // Fluorine, chlorine, and bromine
  result[ 9] = {   0, 176,  80, 255 };  // [F ] Clover green
  result[17] = {  87, 201,  66, 255 };  // [Cl] Green
  result[35] = { 168,   0,  39, 255 };  // [Br] Psivant Red - accent #1, darker 25%

  return result;
}

//-------------------------------------------------------------------------------------------------
std::string formatColor(const uchar4 color_in, const GridFileSyntax syntax) {
  switch (syntax) {
  case GridFileSyntax::MATPLOTLIB:
  case GridFileSyntax::MATRIX_PKG:
    return "'color', [ " + realToString(static_cast<double>(color_in.x) / 255.0,
                                        5, 3, NumberFormat::STANDARD_REAL) + ", " +
           realToString(static_cast<double>(color_in.y) / 255.0,
                        5, 3, NumberFormat::STANDARD_REAL) + ", " +
           realToString(static_cast<double>(color_in.z) / 255.0,
                        5, 3, NumberFormat::STANDARD_REAL) + " ]";
  case GridFileSyntax::OPEN_DX:

    // Assume that PyMol will be used to display the scene.
    return "color " + rgbHexCode(color_in);
  case GridFileSyntax::CUBEGEN:
    rtErr("Gaussian cube files are not formatted for direct visualization with packages for which "
          "STORMM encodes a color plotting syntax.", "formatColor");
  }
  __builtin_unreachable();
}
  
} // namespace review
} // namespace stormm
