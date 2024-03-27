#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Topology/atomgraph_enumerators.h"
#include "amber_netcdf.h"
#include "netcdf_util.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void AmberNetcdf::AmberNetcdf(const std::string &filename_in, const CoordinateFileKind outkind,
                              const PrintSituation expectation, const int atom_count_in,
                              const UnitCellType box_type, const std::string &title_in,
                              const std::string &application_in, const double temperature_in) :
  file_name{filename_in},
  current_temperature{temperature_in},
  restartTime{0.0},
  output_kind{outkind},
  netcdf_id{-1},
  frame_dimension_id{-1},
  ncframe{-1},
  current_frame{0},
  atom_dimension_id{-1},
  atom_count{atom_count_in},
  coordinate_variable_id{-1},
  velocity_variable_id{-1},
  cell_angle_variable_id{-1},
  cell_length_variable_id{-1},
  spatial_dimension_id{-1},
  label_dimension_id{-1},
  cell_spatial_dimension_id{-1},
  cell_angular_dimension_id{-1},
  spatial_variable_id{-1},
  time_variable_id{-1},
  cell_spatial_variable_id{-1},
  cell_angular_variable_id{-1},
  temperature_variable_id{-1}  
{
  // Determine the creation mode (this becomes an opening mode for the APPEND case)
  netcdf_id = ncdfCreate(file_name, expectation);

  // Begin filling in some critical variables.  Weave the trajectory and restart formats together,
  // as they often take similar attributes.
  int dims[NC_MAX_VAR_DIMS];
  nc_type iprec = (output_kind == CoordinateFileKind::AMBER_NETCDF) ? NC_FLOAT : NC_DOUBLE;
  int ndim;
  switch (output_kind) {
  case CoordinateFileKind::AMBER_NETCDF:
    frame_dimension_id = ncdfDefineDimension(netcdf_id, AncdfVariable::NCFRAME, NC_UNLIMITED,
                                             "defining frame dimensions");
    dims[0] = frame_dimension_id;
    ndim = 1;
    break;
  case CoordinateFileKind::AMBER_NETCDF_RST:
    ndim = 0;
    break;
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::SDF:
    rtErr("The Amber NetCDF writer should not be called to output an " +
          getEnumerationName(output_kind) + " file " + file_name + ".", "AmberNetcdf");
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + file_name + " is unknown.", "AmberNetcdf");
  }
  time_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCTIME, iprec, ndim, dims,
                                        "defining the time variable");

  // Both formats take the same units in the time variable, spatial dimension setup, and atom
  // count setup.
  ncdfPlaceAttributeText(netcdf_id, time_variable_id, "units", "picosecond",
                         "writing time variable units");
  spatial_dimension_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCSPATIAL, 3,
                                            "defining spatial dimensions");
  dims[0] = spatial_dimension_id;
  spatial_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCSPATIAL, NC_CHAR, 1, dims,
                                           "defining spatial variable");
  atom_dimension_id = ncdfDefineDimension(netcdf_id, AncdfVariable::NCATOM, atom_count,
                                          "defining atom dimension");

  // No need for switches beyond this point--there are only two choices and other cases have
  // been trapped.  The trajectory format works in NC_FLOAT while the restart format works in
  // NC_DOUBLE.
  if (output_kind == CoordinateFileKind::AMBER_NETCDF) {
    dims[0] = frame_dimension_id;
    dims[1] = atom_dimension_id;
    dims[2] = spatial_dimension_id;
    ndim = 3;
  }
  else {
    dims[0] = atom_dimension_id;
    dims[1] = spatial_dimension_id;
    ndim = 2;
  }
  coordinate_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCCOORDS, iprec, ndim,
                                              dims, "defining coordinates variable");
  ncdfPlaceAttributeText(netcdf_id, coordinate_variable_id, "units", "angstrom",
                         "writing coordinates variable units");

  // Velocity writing for the restart file
  if (output_kind == CoordinateFileKind::AMBER_NETCDF_RST) {
    velocity_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCVELO, iprec, ndim, dims,
                                              "defining velocities variable");
    ncdfPlaceAttributeText(netcdf_id, velocity_variable_id, "units", "angstrom/picosecond",
                           "writing velocities variable units");
    const double velocityScale = 20.455;
    checkNetcdfStatus(nc_put_add_double(netcdf_id, velocity_variable_id, "scale_factor", iprec, 1,
                                        &velocity_scale), "writing velocities scale factor");
  }

  // Unit cell spatial dimensions
  cell_spatial_dims = ncdfDefineDimension(netcdf_id, AncdfVariable::NCCELL_SPATIAL, 3,
                                          "defining cell spatial dimension.");
  dims[0] = cell_spatial_dimension_id;
  cell_spatial_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCCELL_SPATIAL, NC_CHAR,
                                                1, dims, "defining cell spatial variable");

  // Unit cell angle values
  label_dims = ncdfDefineDimension(netcdf_id, AncdfVariable::NCLABEL,
                                   getAncdfVariableName(AncdfVariable::NCLABEL).size(),
                                   "defining label dimension");
  cell_angular_dimension_id = ncdfDefineDimension(netcdf_id, AncdfVariable::NCCELL_ANGULAR, 3,
                                                  "defining cell angular dimension");
  dims[0] = cell_angular_dimension_id;
  dims[1] = label_dimension_id;
  cell_angular_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCCELL_ANGULAR, NC_CHAR,
                                                2, dims, "defining cell angular variable");

  // Box information
  switch (box_type) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    dims[0] = frame_dimension_id;
    dims[1] = cell_spatial_dimension_id;
    cell_length_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCCELL_LENGTHS,
                                                 NC_DOUBLE, 2, dims,
                                                 "defining cell length variable");
    ncdfPlaceAttributeText(netcdf_id, cell_length_variable_id, "units", 8, "angstrom"),
                      "writing cell length variable units");
    dims[1] = cell_angular_dimension_id;
    cell_angle_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCCELL_ANGLES, NC_DOUBLE,
                                                2, dims, "defining cell angle variable");
    ncdfPlaceAttributeText(netcdf_id, cell_angle_variable_id, "units", 6, "degree"),
                      "writing cell angle variable units");
    break;
  }

  // Title information and other attributes
  ncdfPlaceAttributeText(netcdf_id, NC_GLOBAL, "title", title_in, "writing title");
  ncdfPlaceAttributeText(netcdf_id, NC_GLOBAL, "application", application_in,
                         "writing application");
  ncdfPlaceAttributeText(netcdf_id, NC_GLOBAL, "library", "STORMM", "writing library");
  ncdfPlaceAttributeText(netcdf_id, NC_GLOBAL, "library_version", stormm_version,
                         "writing library version");
  ncdfPlaceAttributeText(netcdf_id, NC_GLOBAL, "conventions", "AMBER", "writing conventions");
  ncdfPlaceAttributeText(netcdf_id, NC_GLOBAL, "convention_version", "1.0",
                         "writing convention version");

  // Replica temperature is for restart files only
  if (output_kind == CoordinateFileKind::AMBER_NETCDF_RST) {
    temperature_variable_id = ncdfDefineVariable(netcdf_id, AncdfVariable::NCTEMPERATURE,
                                                 NC_DOUBLE, 0, dims, "defining replica "
                                                 "temperature in netcdf restart");
    ncdfPlaceAttributeText(netcdf_id, temperature_variable_id, "units", "kelvin",
                           "Defining replica temperature units in netcdf restart");
  }

  // Set the fill mode to "do not fill at this time."
  ncdfSetFillMode(netcdf_id, NC_NOFILL, dimensions, "turning off NetCDF filling");

  // End the definitions for the NetCDF file
  ncdfEndDefinitions(netcdf_id, "error upon capping NetCDF definitions");
}

//-------------------------------------------------------------------------------------------------
void AmberNetcdf::open() {
  checkNetcdfStatus(nc_open(file_name.c_str(), NC_WRITE, &netcdf_id), "(re)opening file");
}

//-------------------------------------------------------------------------------------------------
void AmberNetcdf::close() {
  checkNetcdfStatus(nc_close(netcdf_id), "closing file");
}

} // namespace trajectory
} // namespace stormm
