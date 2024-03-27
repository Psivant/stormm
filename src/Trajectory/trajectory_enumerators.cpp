#include <fstream>
#include <iostream>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getCoordinateFileKindDescription(const CoordinateFileKind cfkind) {
  switch (cfkind) {
  case CoordinateFileKind::AMBER_CRD:
    return std::string("Amber ascii trajectory");
  case CoordinateFileKind::AMBER_INPCRD:
    return std::string("Amber input coordinates");
  case CoordinateFileKind::AMBER_ASCII_RST:
    return std::string("Amber ascii restart");
  case CoordinateFileKind::AMBER_NETCDF:
    return std::string("Amber NetCDF binary trajectory");
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return std::string("Amber NetCDF binary restart");
  case CoordinateFileKind::SDF:
    return std::string("MDL MOL / SDF");
  case CoordinateFileKind::UNKNOWN:
    return std::string("Unknown coordinate file format");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateFileKind cfkind) {
  switch (cfkind) {
  case CoordinateFileKind::AMBER_CRD:
    return std::string("AMBER_CRD");
  case CoordinateFileKind::AMBER_INPCRD:
    return std::string("AMBER_INPCRD");
  case CoordinateFileKind::AMBER_ASCII_RST:
    return std::string("AMBER_ASCII_RST");
  case CoordinateFileKind::AMBER_NETCDF:
    return std::string("AMBER_NETCDF");
  case CoordinateFileKind::AMBER_NETCDF_RST:
    return std::string("AMBER_NETCDF_RST");
  case CoordinateFileKind::SDF:
    return std::string("SDF");
  case CoordinateFileKind::UNKNOWN:
    return std::string("UNKNOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const AncdfVariable key) {
  switch (key) {
  case AncdfVariable::NCFRAME:
    return std::string("frame");
  case AncdfVariable::NCSPATIAL:
    return std::string("spatial");
  case AncdfVariable::NCATOM:
    return std::string("atom");
  case AncdfVariable::NCCELL_SPATIAL:
    return std::string("cell_spatial");
  case AncdfVariable::NCCELL_LENGTHS:
    return std::string("cell_lengths");
  case AncdfVariable::NCCELL_ANGULAR:
    return std::string("cell_angular");
  case AncdfVariable::NCCELL_ANGLES:
    return std::string("cell_angles");
  case AncdfVariable::NCCOORDS:
    return std::string("coordinates");
  case AncdfVariable::NCVELO:
    return std::string("velocities");
  case AncdfVariable::NCTEMPERATURE:
    return std::string("temp0");
  case AncdfVariable::NCTIME:
    return std::string("time");
  case AncdfVariable::NCLABEL:
    return std::string("label");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TrajectoryKind input) {
  switch (input) {
  case TrajectoryKind::POSITIONS:
    return std::string("POSITIONS");
  case TrajectoryKind::VELOCITIES:
    return std::string("VELOCITIES");
  case TrajectoryKind::FORCES:
    return std::string("FORCES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateCycle orientation) {
  switch (orientation) {
  case CoordinateCycle::WHITE:
    return std::string("WHITE");
  case CoordinateCycle::BLACK:
    return std::string("BLACK");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateLineFormat input) {
  switch (input) {
  case CoordinateLineFormat::FIXED_COLUMN:
    return std::string("FIXED_COLUMN");
  case CoordinateLineFormat::FREE_FORMAT:
    return std::string("FREE_FORMAT");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CoordinateFileRole cpkind) {
  switch (cpkind) {
  case CoordinateFileRole::INITIATE:
    return std::string("INITIATE");
  case CoordinateFileRole::TRAJECTORY:
    return std::string("TRAJECTORY");
  case CoordinateFileRole::CHECKPOINT:
    return std::string("CHECKPOINT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TrajectoryFusion protocol) {
  switch (protocol) {
  case TrajectoryFusion::ON:
    return std::string("ON");
  case TrajectoryFusion::OFF:
    return std::string("OFF");
  case TrajectoryFusion::AUTO:
    return std::string("AUTO");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ThermostatKind input) {
  switch (input) {
  case ThermostatKind::NONE:
    return std::string("NONE");
  case ThermostatKind::ANDERSEN:
    return std::string("ANDERSEN");
  case ThermostatKind::LANGEVIN:
    return std::string("LANGEVIN");
  case ThermostatKind::BERENDSEN:
    return std::string("BERENDSEN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ThermostatPartition input) {
  switch (input) {
  case ThermostatPartition::COMMON:
    return std::string("COMMON");
  case ThermostatPartition::SYSTEMS:
    return std::string("SYSTEMS");
  case ThermostatPartition::ATOMS:
    return std::string("ATOMS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EnforceExactTemperature input) {
  switch (input) {
  case EnforceExactTemperature::YES:
    return std::string("YES");
  case EnforceExactTemperature::NO:
    return std::string("NO");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const IntegrationStage input) {
  switch (input) {
  case IntegrationStage::VELOCITY_ADVANCE:
    return std::string("VELOCITY_ADVANCE");
  case IntegrationStage::VELOCITY_CONSTRAINT:
    return std::string("VELOCITY_CONSTRAINT");
  case IntegrationStage::POSITION_ADVANCE:
    return std::string("POSITION_ADVANCE");
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    return std::string("GEOMETRY_CONSTRAINT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind translateCoordinateFileKind(const std::string &name_in) {
  if (strcmpCased(name_in, "AMBER_CRD", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_CRD;
  }
  else if (strcmpCased(name_in, "AMBER_INPCRD", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_INPCRD;
  }
  else if (strcmpCased(name_in, "AMBER_ASCII_RST", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_ASCII_RST;
  }
  else if (strcmpCased(name_in, "AMBER_NETCDF", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_NETCDF;
  }
  else if (strcmpCased(name_in, "AMBER_NETCDF_RST", CaseSensitivity::NO)) {
    return CoordinateFileKind::AMBER_NETCDF_RST;
  }
  else if (strcmpCased(name_in, "SDF", CaseSensitivity::NO) ||
           strcmpCased(name_in, "MDL_MOL", CaseSensitivity::NO)) {
    return CoordinateFileKind::SDF;
  }
  else if (strcmpCased(name_in, "UNKNOWN", CaseSensitivity::NO)) {
    return CoordinateFileKind::UNKNOWN;
  }
  else {
    rtErr("Unrecognized coordinate file enumeration " + name_in + ".",
          "translateCoordinateFileKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ThermostatKind translateThermostatKind(const std::string &input) {
  if (strcmpCased(input, "none", CaseSensitivity::NO) ||
      strcmpCased(input, "0", CaseSensitivity::NO) ||
      strcmpCased(input, "isoenergetic", CaseSensitivity::NO)) {
    return ThermostatKind::NONE;
  }
  else if (strcmpCased(input, "berendsen", CaseSensitivity::NO) ||
           strcmpCased(input, "1", CaseSensitivity::NO) ||
           strcmpCased(input, "exp_rescale", CaseSensitivity::NO)) {
    return ThermostatKind::BERENDSEN;
  }
  else if (strcmpCased(input, "mass_andersen", CaseSensitivity::NO) ||
           strcmpCased(input, "m_andersen", CaseSensitivity::NO) ||
           strcmpCased(input, "2", CaseSensitivity::NO) ||
           strcmpCased(input, "mass_reset", CaseSensitivity::NO)) {
    return ThermostatKind::ANDERSEN;
  }
  else if (strcmpCased(input, "langevin", CaseSensitivity::NO) ||
           strcmpCased(input, "3", CaseSensitivity::NO)) {
    return ThermostatKind::LANGEVIN;
  }
  else {
    rtErr("Unrecognized thermostat type enumeration " + input + ".",
          "translateThermostatKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CoordinateCycle getNextCyclePosition(const CoordinateCycle orientation) {
  CoordinateCycle result;
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return CoordinateCycle::WHITE;
  case CoordinateCycle::WHITE:
    return CoordinateCycle::BLACK;
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace stormm
