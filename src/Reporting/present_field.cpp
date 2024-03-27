#include "copyright.h"
#include "present_field.h"

namespace stormm {
namespace review {

//-------------------------------------------------------------------------------------------------
std::string getSceneDataFileName(const std::string &scene_file) {
  std::string fn_before, fn_after; 
  splitPath(scene_file, &fn_before, &fn_after);
  return fn_before + "_data" + ((fn_after.size() > 0) ? ".dat" : "");
}
  
//-------------------------------------------------------------------------------------------------
std::string drawFieldBorders(const MeshParameters &mps, const RenderOptions &ropt,
                             const GridFileSyntax syntax, const std::string &path_name) {
  std::string result;
  switch (syntax) {
  case GridFileSyntax::MATRIX_PKG:
  case GridFileSyntax::MATPLOTLIB:
    result = protectText("Trace a pattern to define the boundaries of the region containing the "
                         "potential field.", commentSymbol(syntax));
    break;
  case GridFileSyntax::OPEN_DX:
    break;
  case GridFileSyntax::CUBEGEN:
    rtErr("Gaussian cube files are not compatible with visualization in a way that can script the "
          "border representation.", "drawFieldBorders");
  }
  switch (syntax) {
  case GridFileSyntax::MATRIX_PKG:
    result += path_name + " = [\n";
    result += "    0    0    0\n na-1    0 nc-1\n na-1    0    0\n na-1 nb-1 nc-1\n";
    result += " na-1 nb-1    0\n    0 nb-1 nc-1\n    0 nb-1    0\n    0    0 nc-1\n";
    result += " na-1 nb-1 nc-1\n na-1    0 nc-1\n    0 nb-1 nc-1\n    0    0 nc-1\n";
    result += "    0    0    0\n na-1 nb-1    0\n na-1    0    0\n    0 nb-1    0\n";
    result += "    0    0    0\n    0 nb-1 nc-1\n na-1 nb-1 nc-1\n    0 nb-1    0\n";
    result += " na-1 nb-1    0\n na-1    0 nc-1\n    0    0 nc-1\n na-1    0    0\n";
    result += "    0    0    0\n";
    result += "];\n";
    result += path_name + " = (invU * " + path_name + "')';\n";
    result += path_name + "(:,1) = " + path_name + "(:,1) + mg_origin(1);\n";
    result += path_name + "(:,2) = " + path_name + "(:,2) + mg_origin(2);\n";
    result += path_name + "(:,3) = " + path_name + "(:,3) + mg_origin(3);\n";
    result += "plot3(" + path_name + "(:,1), " + path_name + "(:,2), " + path_name + "(:,3), " +
              "'linewidth', " +
              realToString(ropt.getBorderLineWidth(), 4, 2, NumberFormat::STANDARD_REAL) + ", " +
              ropt.getFormattedBorderColor(syntax);
    switch (ropt.getBorderLineStyle()) {
    case LinePlotStyle::SOLID:
      break;
    case LinePlotStyle::DASHED:
      result += "'linestyle', '--'";
      break;
    }
    result += ");\n";
    break;
  case GridFileSyntax::MATPLOTLIB:
    break;
  case GridFileSyntax::OPEN_DX:
  case GridFileSyntax::CUBEGEN:
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void drawFieldBorders(std::ofstream *foutp, const MeshParameters &mps, const RenderOptions &ropt,
                      const GridFileSyntax syntax, const std::string &path_name) {
  const std::string contrib = drawFieldBorders(mps, ropt, syntax, path_name);
  foutp->write(contrib.data(), contrib.size());
}
                      
} // namespace review
} // namespace stormm
