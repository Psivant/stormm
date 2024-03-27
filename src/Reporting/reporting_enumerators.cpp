#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "reporting_enumerators.h"
#include "error_format.h"

namespace stormm {
namespace display {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HelpSignalKind input) {
  switch (input) {
  case HelpSignalKind::NO_ARGS:
    return std::string("NO_ARGS");
  case HelpSignalKind::NO_ARGS_ONLY:
    return std::string("NO_ARGS_ONLY");
  case HelpSignalKind::KEYWORD:
    return std::string("KEYWORD");
  case HelpSignalKind::KEYWORD_ONLY:
    return std::string("KEYWORD_ONLY");
  }
  __builtin_unreachable();
}

} // namespace display
  
namespace review {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const OutputScope input) {
  switch (input) {
  case OutputScope::AVERAGES:
    return std::string("AVERAGES");
  case OutputScope::OUTLIERS:
    return std::string("OUTLIERS");
  case OutputScope::CLUSTER_AVERAGES:
    return std::string("CLUSTER_AVERAGES");
  case OutputScope::CLUSTER_OUTLIERS:
    return std::string("CLUSTER_OUTLIERS");
  case OutputScope::FULL:
    return std::string("FULL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const OutputSyntax input) {
  switch (input) {
  case OutputSyntax::MATPLOTLIB:
    return std::string("MATPLOTLIB");
  case OutputSyntax::MATRIX_PKG:
    return std::string("MATRIX_PKG");
  case OutputSyntax::STANDALONE:
    return std::string("STANDALONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const IntegrationStage input) {
  switch (input) {
  case IntegrationStage::BEGIN:
    return std::string("BEGIN");
  case IntegrationStage::NONBONDED_FORCE:
    return std::string("NONBONDED_FORCE");
  case IntegrationStage::ADDED_FORCE:
    return std::string("ADDED_FORCE");
  case IntegrationStage::BONDED_FORCE:
    return std::string("BONDED_FORCE");
  case IntegrationStage::VELOCITY_UPDATE:
    return std::string("VELOCITY_UPDATE");
  case IntegrationStage::RATTLE:
    return std::string("RATTLE");
  case IntegrationStage::COORDINATE_UPDATE:
    return std::string("COORDINATE_UPDATE");
  case IntegrationStage::SHAKE:
    return std::string("SHAKE");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const GridFileSyntax input) {
  switch (input) {
  case GridFileSyntax::MATPLOTLIB:
    return std::string("MATPLOTLIB");
  case GridFileSyntax::MATRIX_PKG:
    return std::string("MATRIX_PKG");
  case GridFileSyntax::OPEN_DX:
    return std::string("OPEN_DX");
  case GridFileSyntax::CUBEGEN:
    return std::string("CUBEGEN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SectionComponent input) {
  switch (input) {
  case SectionComponent::NARRATIVE:
    return std::string("NARRATIVE");
  case SectionComponent::LIST:
    return std::string("LIST");
  case SectionComponent::TABLE:
    return std::string("TABLE");
  case SectionComponent::SCRIPT:
    return std::string("SCRIPT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ListEnumeration input) {
  switch (input) {
  case ListEnumeration::BULLET:
    return std::string("BULLET");
  case ListEnumeration::NUMBERED:
    return std::string("NUMBERED");
  case ListEnumeration::ALPHABETIC:
    return std::string("ALPHABETIC");
  case ListEnumeration::ROMAN:
    return std::string("ROMAN");
  case ListEnumeration::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TextEnds input) {
  switch (input) {
  case TextEnds::AS_IS:
    return std::string("AS_IS");
  case TextEnds::NEWLINE:
    return std::string("NEWLINE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TableContentKind input) {
  switch (input) {
  case TableContentKind::INTEGER:
    return std::string("INTEGER");
  case TableContentKind::REAL:
    return std::string("REAL");
  case TableContentKind::STRING:
    return std::string("STRING");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SurfaceRender input) {
  switch (input) {
  case SurfaceRender::WIRE:
    return std::string("WIRE");
  case SurfaceRender::SOLID:
    return std::string("SOLID");
  case SurfaceRender::SCAFFOLD:
    return std::string("SCAFFOLD");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const LinePlotStyle input) {
  switch (input) {
  case LinePlotStyle::SOLID:
    return std::string("SOLID");
  case LinePlotStyle::DASHED:
    return std::string("DASHED");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PsivantColor input) {
  switch (input) {
  case PsivantColor::BLACK:
    return std::string("BLACK");
  case PsivantColor::RED:
    return std::string("RED");
  case PsivantColor::YELLOW:
    return std::string("YELLOW");
  case PsivantColor::BLUE:
    return std::string("BLUE");
  case PsivantColor::PURPLE:
    return std::string("PURPLE");
  case PsivantColor::GREY:
    return std::string("GREY");
  case PsivantColor::LIGHT_RED:
    return std::string("LIGHT_RED");
  case PsivantColor::LIGHT_YELLOW:
    return std::string("LIGHT_YELLOW");
  case PsivantColor::LIGHT_BLUE:
    return std::string("LIGHT_BLUE");
  case PsivantColor::LIGHT_PURPLE:
    return std::string("LIGHT_PURPLE");
  case PsivantColor::WHITE:
    return std::string("WHITE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
OutputScope translateOutputScope(const std::string &input) {
  if (strcmpCased(input, "average", CaseSensitivity::NO) ||
      strcmpCased(input, "averages", CaseSensitivity::NO) ||
      strcmpCased(input, "mean", CaseSensitivity::NO) ||
      strcmpCased(input, "means", CaseSensitivity::NO)) {
    return OutputScope::AVERAGES;
  }
  else if (strcmpCased(input, "outlier", CaseSensitivity::NO) ||
           strcmpCased(input, "outliers", CaseSensitivity::NO)) {
    return OutputScope::OUTLIERS;
  }
  else if (strcmpCased(input, "cluster_average", CaseSensitivity::NO) ||
           strcmpCased(input, "cluster_averages", CaseSensitivity::NO) ||
           strcmpCased(input, "clusteraverage", CaseSensitivity::NO) ||
           strcmpCased(input, "clusteraverages", CaseSensitivity::NO)) {
    return OutputScope::CLUSTER_AVERAGES;
  }
  else if (strcmpCased(input, "cluster_outlier", CaseSensitivity::NO) ||
           strcmpCased(input, "cluster_outliers", CaseSensitivity::NO) ||
           strcmpCased(input, "clusteroutlier", CaseSensitivity::NO) ||
           strcmpCased(input, "clusteroutliers", CaseSensitivity::NO)) {
    return OutputScope::CLUSTER_OUTLIERS;
  }
  else if (strcmpCased(input, "all", CaseSensitivity::NO) ||
           strcmpCased(input, "full", CaseSensitivity::NO) ||
           strcmpCased(input, "entire", CaseSensitivity::NO) ||
           strcmpCased(input, "complete", CaseSensitivity::NO)) {
    return OutputScope::FULL;
  }
  else {
    rtErr("No OutputScope enumeration matches \"" + input + "\".", "translateOutputScope");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
IntegrationStage translateIntegrationStage(const std::string &input) {
  if (strcmpCased(input, "begin", CaseSensitivity::NO) ||
      strcmpCased(input, "init", CaseSensitivity::NO)) {
    return IntegrationStage::BEGIN;
  }
  else if (strcmpCased(input, "nonbonded_force", CaseSensitivity::NO) ||
           strcmpCased(input, "nonbonded", CaseSensitivity::NO)) {
    return IntegrationStage::NONBONDED_FORCE;
  }
  else if (strcmpCased(input, "added_force", CaseSensitivity::NO) ||
           strcmpCased(input, "added", CaseSensitivity::NO)) {
    return IntegrationStage::ADDED_FORCE;
  }
  else if (strcmpCased(input, "bonded_force", CaseSensitivity::NO) ||
           strcmpCased(input, "bonded", CaseSensitivity::NO)) {
    return IntegrationStage::BONDED_FORCE;
  }
  else if (strcmpCased(input, "velocity_update", CaseSensitivity::NO) ||
           strcmpCased(input, "velocity", CaseSensitivity::NO)) {
    return IntegrationStage::VELOCITY_UPDATE;
  }
  else if (strcmpCased(input, "rattle", CaseSensitivity::NO) ||
           strcmpCased(input, "velocity_constraint", CaseSensitivity::NO)) {
    return IntegrationStage::RATTLE;
  }
  else if (strcmpCased(input, "coordinate_update", CaseSensitivity::NO) ||
           strcmpCased(input, "coordinate", CaseSensitivity::NO)) {
    return IntegrationStage::COORDINATE_UPDATE;
  }
  else if (strcmpCased(input, "shake", CaseSensitivity::NO) ||
           strcmpCased(input, "coordinate_constraint", CaseSensitivity::NO)) {
    return IntegrationStage::SHAKE;
  }
  else {
    rtErr("No IntegrationStage enumeration matches \"" + input + "\".",
          "translateIntegrationStage");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
SurfaceRender translateSurfaceRender(const std::string &input) {
  if (strcmpCased(input, "wire", CaseSensitivity::NO) ||
      strcmpCased(input, "wire mesh", CaseSensitivity::NO) ||
      strcmpCased(input, "wire_mesh", CaseSensitivity::NO)) {
    return SurfaceRender::WIRE;
  }
  else if (strcmpCased(input, "solid", CaseSensitivity::NO)) {
    return SurfaceRender::SOLID;
  }
  else if (strcmpCased(input, "scaffold", CaseSensitivity::NO)) {
    return SurfaceRender::SCAFFOLD;    
  }
  else {
    rtErr("No SurfaceRender enumeration matches \"" + input + "\".", "translateSurfaceRender");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
LinePlotStyle translateLinePlotStyle(const std::string &input) {
  if (strcmpCased(input, "solid", CaseSensitivity::NO) ||
      strcmpCased(input, "continuous", CaseSensitivity::NO)) {
    return LinePlotStyle::SOLID;
  }
  else if (strcmpCased(input, "dashed", CaseSensitivity::NO) ||
           strcmpCased(input, "broken", CaseSensitivity::NO)) {
    return LinePlotStyle::DASHED;
  }
  else {
    rtErr("No LinePlotStyle enumeration matches \"" + input + "\".", "translateLinePlotStyle");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::string& toRoman(const int x, const ExceptionResponse policy) {
  if (x <= 0 || x > maximum_roman_numeral) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No Roman numeral is available for " + std::to_string(x) + ".", "toRoman");
      break;
    case ExceptionResponse::WARN:
      rtWarn("No Roman numeral is available for " + std::to_string(x) + ".  White space will be "
             "returned instead.", "toRoman");
      return roman_numerals[maximum_roman_numeral];
    case ExceptionResponse::SILENT:
      return roman_numerals[maximum_roman_numeral];
    }
  }
  else {
    return roman_numerals[x - 1];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string encodePsivantColor(const PsivantColor color, const OutputSyntax syntax) {
  std::string result;
  switch (syntax) {
  case OutputSyntax::MATPLOTLIB:
    result += "color = '";
    switch (color) {
    case PsivantColor::BLACK:
      result += "#111c24";
      break;
    case PsivantColor::RED:
      result += "#e00034";
      break;
    case PsivantColor::YELLOW:
      result += "#fdc82f";
      break;
    case PsivantColor::BLUE:
      result += "#007ac9";
      break;
    case PsivantColor::PURPLE:
      result += "#8e258d";
      break;
    case PsivantColor::GREY:
      result += "#a6a6a6";
      break;
    case PsivantColor::LIGHT_RED:
      result += "#ff8da7";
      break;
    case PsivantColor::LIGHT_YELLOW:
      result += "#fee9ac";
      break;
    case PsivantColor::LIGHT_BLUE:
      result += "#83ceff";
      break;
    case PsivantColor::LIGHT_PURPLE:
      result += "#e496e3";
      break;
    case PsivantColor::WHITE:
      result += "#ffffff";
      break;
    }
    result += "'";
    break;
  case OutputSyntax::MATRIX_PKG:
    result += "'color', [ ";
    switch (color) {
    case PsivantColor::BLACK:
      result += "0.067, 0.110, 0.141";
      break;
    case PsivantColor::RED:
      result += "0.878, 0.000, 0.204";
      break;
    case PsivantColor::YELLOW:
      result += "0.992, 0.784, 0.184";
      break;
    case PsivantColor::BLUE:
      result += "0.000, 0.478, 0.788";
      break;
    case PsivantColor::PURPLE:
      result += "0.557, 0.145, 0.553";
      break;
    case PsivantColor::GREY:
      result += "0.651, 0.651, 0.651";
      break;
    case PsivantColor::LIGHT_RED:
      result += "1.000, 0.553, 0.655";
      break;
    case PsivantColor::LIGHT_YELLOW:
      result += "0.996, 0.914, 0.675";
      break;
    case PsivantColor::LIGHT_BLUE:
      result += "0.514, 0.808, 1.000";
      break;
    case PsivantColor::LIGHT_PURPLE:
      result += "0.894, 0.588, 0.890";
      break;
    case PsivantColor::WHITE:
      result += "1.000, 1.000, 1.000";
      break;
    }
    result += " ]";
    break;
  case OutputSyntax::STANDALONE:
    break;
  }
  return result;
}
  
} // namespace review
} // namespace stormm
