// -*-c++-*-
#ifndef STORMM_TRAJECTORY_ENUMERATORS_H
#define STORMM_TRAJECTORY_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace trajectory {

/// \brief Options for the type of coordinate file to write
enum class CoordinateFileKind {
  AMBER_CRD,         ///< The awful .crd format trajectory file, 10 x %8.3f per line, with or
                     ///<   without box coordinates
  AMBER_INPCRD,      ///< The still-used input coordinates format, as output by tleap, 6 x %12.7f
                     ///<   per line, with or without box coordinates
  AMBER_ASCII_RST,   ///< The deprecated ascii formatted Amber restart file format, 6 x %12.7f
                     ///<   per line, with velocities and with or without box coordinates (a
                     ///<   restart file with no velocities is an INPUT_COORDINATES file)
  AMBER_NETCDF,      ///< The binary Amber NetCDF trajectory format
  AMBER_NETCDF_RST,  ///< The binary Amber NetCDF restart format
  SDF,               ///< The trajectory is an MDL MOL file, or a collection of them
  UNKNOWN            ///< The coordinate file kind is not (yet) understood
};

/// \brief An enumerator to track different Amber NetCDF variable identifiers.  The associated
///        translation function produces the strings of interest.
enum class AncdfVariable {
  NCFRAME,         ///< 
  NCSPATIAL,       ///<
  NCATOM,          ///<
  NCCELL_SPATIAL,  ///<
  NCCELL_LENGTHS,  ///<
  NCCELL_ANGULAR,  ///<
  NCCELL_ANGLES,   ///<
  NCCOORDS,        ///<
  NCVELO,          ///<
  NCTEMPERATURE,   ///<
  NCTIME,          ///<
  NCLABEL          ///<
};

/// \brief Enumerate the kind of information that a trajectory can contain
enum class TrajectoryKind {
  POSITIONS,   // Positions of particles
  VELOCITIES,  // Velocities of particles
  FORCES       // Total forces on all particles (the two force arrays are summed in cases where
               //   the forces are partitioned)
};

/// \brief Enumerate the cycle of time for coordinates in objects that support iterative molecular
///        dynamics.  With each time step, future coordinates are built based on past and present
///        coordinates (and forces).  As the step advances, present recedes to past, future becomes
///        present, and what was past becomes the staging ground for the forthcoming future
///        configuration.
enum class CoordinateCycle {
  BLACK, ///< Set the alternate positions arrays (x_alt_coordinates...) as holding the current
         ///<   coordinates.
  WHITE  ///< Set the typical positions arrays (x_coordinates...) as current
};

/// \brief Differentiate between fixed-column and free-format (but ordered) files.
enum class CoordinateLineFormat {
  FIXED_COLUMN,
  FREE_FORMAT
};

/// \brief Enumerate the purposes that a coordinate file can have.  One could imagine a matrix of
///        combinations of these enumerations and those of the CoordinateFileKind enumerator,
///        showing which formats are suitable for which purposes.
enum class CoordinateFileRole {
  INITIATE,    ///< Contents of the file serve to seed the initial states of one or more systems
  TRAJECTORY,  ///< The file will collect or present a trajectory of a single system
  CHECKPOINT   ///< The contents of the file will record final states of one or more systems
};

/// \brief Enumerate the protocols for fusing trajectories and checkpoint files.
enum class TrajectoryFusion {
  ON,   ///< Fuse trajectory files as well as checkpoint files (if possible) for systems grouped
        ///<   under the same label
  OFF,  ///< Do not fuse trajectory or checkpoint files (even if possible)
  AUTO  ///< Fuse checkpoint files (if the format will support it) but not trajectories for systems
        ///<   grouped under the same label
};

/// \brief Enumerate the various thermostats available for simulations
enum class ThermostatKind {
  NONE,      ///< No thermostating takes place.  The simulation explores a microcanonical "NVE"
             ///<   ensemble.
  ANDERSEN,  ///< Andersen thermostating, with periodic velocity reassignments of all particles,
             ///<   is in effect.
  LANGEVIN,  ///< Langevin thermostating, with constant velocity adjustments of all particles,
             ///<   is in effect.
  BERENDSEN  ///< Berendsen thermostating, with velocity rescaling and all of the problems it
             ///<   might confer, is in effect.
};

/// \brief Enumerate the ways in which different atoms of one or more simulations may have unique
///        temperature baths.
enum class ThermostatPartition {
  COMMON,   ///< All atoms of the one system or multiple systems share a common temperature bath.
            ///<   Its initial and final temperatures may differ.
  SYSTEMS,  ///< Different systems of a synthesis will each have their own unique temperature
            ///<   baths.  Each system may have a unique initial and final temperature.
  ATOMS     ///< Individual atoms may have their own temperature baths.  In most cases this will be
            ///<   used to hold the protein and solvent at different temperatures.
};

/// \brief Due to its random nature, the velocity initialization may not set a system at the exact
///        temperature requested by the applied Andersen thermostat.  Enumerate the options for
///        enforcing a particular temperature.
enum class EnforceExactTemperature {
  YES,  ///< The exact temperature will be enforced by a Berendesen velocity rescaling with
        ///<   infinite coupling strength
  NO    ///< The exact temperature will be enforce
};

/// \brief List the basic stages of the integration time step, to aid in differentiating them
///        across various routines.
enum class IntegrationStage {
  VELOCITY_ADVANCE,     ///< Adjust the velocities one half time step according to the forces at
                        ///<   hand.
  VELOCITY_CONSTRAINT,  ///< Adjust the velocities of particles to ensure that those involved in
                        ///<   constrained bonds move orthogonally to their bond axes.
  POSITION_ADVANCE,     ///< Finalize the velocities with one more advancement based on the current
                        ///<   forces, then advance the particle positions based on the velocities
                        ///<   at hand.
  GEOMETRY_CONSTRAINT   ///< Apply geometric constraints to bond lengths, with commensurate
                        ///<   velocity adjustments if required.
};

/// \brief Produce a description of the coordinate file kind.
///
/// \param cfkind  The enumerator instance of interest
std::string getCoordinateFileKindDescription(const CoordinateFileKind cfkind);
  
/// \brief Translate the CoordinateFileKind enumeration.  Various overloads serve different
///        enumerators.
///
/// \param cfkind    The enumerator instance of interest
/// \param key       The enumerator instance of interest
/// \param cpkind    The enumerator instance of interest
/// \param protocol  The enumerator instance of interest
/// \param input     The enumerator instance of interest
/// \{
std::string getEnumerationName(CoordinateFileKind cfkind);
std::string getEnumerationName(AncdfVariable key);
std::string getEnumerationName(TrajectoryKind input);
std::string getEnumerationName(CoordinateCycle orientation);
std::string getEnumerationName(CoordinateLineFormat input);
std::string getEnumerationName(CoordinateFileRole cpkind);
std::string getEnumerationName(TrajectoryFusion protocol);
std::string getEnumerationName(ThermostatKind input);
std::string getEnumerationName(ThermostatPartition input);
std::string getEnumerationName(EnforceExactTemperature input);
std::string getEnumerationName(IntegrationStage input);
/// \}

/// \brief Translate a string into one of the CoordinateFileKind enumerations.
///
/// \param name_in  The string to translate
CoordinateFileKind translateCoordinateFileKind(const std::string &name_in);

/// \brief Translate a string into one of the ThermostatKind enumerations.
///
/// \param input  The string to translate
ThermostatKind translateThermostatKind(const std::string &input);
  
/// \brief Obtain the next point in the coordinate cycle.
///
/// \param orientation  The current point in the coordinate cycle
CoordinateCycle getNextCyclePosition(CoordinateCycle orientation);
  
/// \brief Translate the AncdfVariable enumeration.  This provides keywords to serve as landmarks
///        in an Amber binary trajectory or restart file.
///
/// \param key  The enumerator instance of interest

} // namespace trajectory
} // namespace stormm

#endif
