// -*-c++-*-
#ifndef STORMM_AMBER_NETCDF_H
#define STORMM_AMBER_NETCDF_H

#include <netcdf.h>
#include "copyright.h"

namespace stormm {
namespace trajectory {

class AmberNetcdf {
public:

  /// \brief The constructor takes a file name and other critical information, which will be
  ///        used to open a NetCDF file and write its preamble
  AmberNetcdf();

  /// \brief Open or re-open the output associated with this object
  void open();

  /// \brief Close the output associated with this object
  void close();

  /// \brief Report errors based on the output of NetCDF library functions called by this object
  void checkNetcdfStatus(const int status, const std::string activity = std::string(""));

private:
  std::string file_name;           ///< Name of the file associated with this object
  double current_temperature;      ///< Temperature of current frame (if temperature_variable_id
                                   ///<   is -1)
  double restartTime;              ///< Simulation time if Amber restart
  CoordinateFileKind output_kind;  ///< 0 if trajectory, 1 if restart
  int netcdf_id;                   ///< Unique NetCDF identifier of the associated file 
  int frame_dimension_id;          ///< ID of frame dimension
  int ncframe;                     ///< Number of frames in the file
  int current_frame;               ///< Current frame number
  int atom_dimension_id;           ///< ID of atom dimension
  int atom_count;                  ///< Number of atoms
  int coordinate_variable_id;      ///< ID of coordinates variable
  int velocity_variable_id;        ///< ID of velocities variable
  int cellAngle_variable_id;       ///< ID of box angle variable
  int cellLength_variable_id;      ///< ID of box length variable
  int spatial_dimension_id;
  int label_dimension_id;
  int cell_spatial_dimension_id;
  int cell_angular_dimension_id;
  int spatial_variable_id;
  int time_variable_id;
  int cell_spatial_variable_id;
  int cell_angular_variable_id;
  int temperature_variable_id;
};

} // namespace trajectory
} // namespace stormm

#endif

